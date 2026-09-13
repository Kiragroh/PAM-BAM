# -*- coding: utf-8 -*-
"""PAM / BAM viewer 1.2 -- standalone, read-only RayStation file script.

Run in RayStation CPython with an open plan; select a target ROI and Calculate.
Dependencies: numpy, tkinter (no installation or repository import at runtime).
API basis: RayStation v2025 SP2 / 17.2.0. Local validation is still required.

AM = fraction of the projected target outside the MLC/jaw opening.
BAM = sum(AM * Segment.RelativeWeight); PAM = MU-weighted mean of beam BAMs.
RelativeWeight is documented by RS as a fraction of total beam MU, NOT a
cumulative meterset. Dynamic delivery is sampled at the native segments/CPs;
there is no reconstruction of continuous leaf travel or dose calculation.

An available closed native triangle mesh is projected directly. Otherwise the
ROI is voxelized (>=128/255 occupancy) and its exact block surface is extracted
in memory. Both paths project surface polygons onto the isocenter BEV plane;
no per-ray 3D traversal and no TPS representation change. In the voxel fallback,
both voxels and BEV use the selected spacing; for native meshes only BEV changes.
Voxel reads use bounded slabs with neighbour halos; the requested grid is never
silently coarsened. Per-beam detail reports API/surface/projection/MLC timings.
Compare resolutions and the displayed surface method when comparing results.
Lengths from RayStation are cm throughout; only the GUI spacing is in mm.

Supported: photon SMLC/DMLC/DynamicArc/StaticArc, HFS, conventional IEC geometry,
static couch yaw, zero pitch/roll/gimbal, one or two MLC layers projected at
isocenter. The intersection of both layers and native jaw limits is used.
Actual leaf widths/centres are read from the plan and checked against the
machine. Jawless machines need no JawPhysics, or valid native fixed field
limits; zero placeholder limits are rejected. No machine-name heuristics.
Setup beams and zero-MU beams are excluded. Any failed treatment beam prevents
a plan PAM. All beam sets must use the same planning examination and target.
No patient/model modification, saving, file export, networking or logging.

Metric and adapter reference: https://github.com/Kiragroh/PAM-BAM
Reference revision: 51f5dd22f4da2c6dbe1c3ce355401791cd0a9a6c
Hernandez et al., Medical Physics (2025), see repository references.
IEC frame convention also cross-checked against RadiotherapyTransformsIEC:
https://github.com/EBATINCA/RadiotherapyTransformsIEC
This implementation adds a sampled ROI projector; it is not the repository's
exact rectangle-union engine. No clinical thresholds are applied.

MIT License -- Copyright (c) 2026 Maximilian Grohmann
Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:
The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.
THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
"""

import itertools
import math
import time
from contextlib import contextmanager
import numpy as np

VERSION = "1.2"
MAX_VOXELS = 16_000_000
MAX_BEV_PIXELS = 600_000
MAX_SURFACE_FACES = 1_000_000


class CalculationError(Exception):
    """An actionable, patient-independent message safe to display."""


class Cancelled(Exception):
    pass


class Timings:
    """In-memory phase timings only; never logged or exported."""
    def __init__(self):
        self.seconds = dict.fromkeys(("API", "Surface", "Projection", "MLC"), 0.0)
        self.calls = dict.fromkeys(self.seconds, 0)
        self.notes = []

    @contextmanager
    def measure(self, phase):
        start = time.perf_counter()
        try:
            yield
        finally:
            self.seconds[phase] += time.perf_counter()-start
            self.calls[phase] += 1

    def summary(self):
        return " | ".join("{} {:.3f} s".format(k, v) for k, v in self.seconds.items()) + \
            " | {} projections".format(self.calls["Projection"]) + ("\n"+"; ".join(self.notes) if self.notes else "")


def require(condition, message):
    if not condition:
        raise CalculationError(message)


def finite(value):
    result = float(value)
    require(math.isfinite(result), "Non-finite geometry or weight.")
    return result


def xyz(point):
    return np.array([finite(point[k]) for k in "xyz"], dtype=float)


def point_dict(values, cast=float):
    return {k: cast(v) for k, v in zip("xyz", values)}


def beam_frame(gantry, couch, collimator):
    """Columns u,v,w: IEC BLD axes and isocenter-to-source in DICOM HFS.

    At G=C=T=0: u=patient left, v=superior, w=anterior.
    IEC active rotations: HFS @ Rz(-couch) @ Ry(gantry) @ Rz(collimator).
    The negative couch yaw expresses the fixed room in patient coordinates.
    """
    def rz(degrees):
        a = math.radians(finite(degrees))
        c, s = math.cos(a), math.sin(a)
        return np.array([[c, -s, 0], [s, c, 0], [0, 0, 1]])
    a = math.radians(finite(gantry))
    c, s = math.cos(a), math.sin(a)
    ry = np.array([[c, 0, s], [0, 1, 0], [-s, 0, c]])
    hfs = np.array([[1, 0, 0], [0, 0, -1], [0, 1, 0]])
    return hfs @ rz(-couch) @ ry @ rz(collimator)


def voxelize(geometry, spacing):
    require(geometry.HasContours(), "Selected ROI has no geometry on the planning CT.")
    lower, upper = [xyz(p) for p in geometry.GetBoundingBox()]
    require(np.all(upper > lower), "ROI bounding box is empty or invalid.")
    corner = lower - spacing
    counts = np.ceil((upper - lower) / spacing).astype(int) + 2
    require(int(np.prod(counts)) <= MAX_VOXELS,
            "ROI grid too large: choose a coarser resolution.")
    values = geometry.GetRoiGeometryAsVoxels(
        Corner=point_dict(corner), VoxelSize=point_dict([spacing] * 3),
        NrVoxels=point_dict(counts, int))
    values = np.asarray(values, dtype=np.uint8)
    require(values.size == int(np.prod(counts)), "RayStation returned an incomplete ROI grid.")
    mask = values.reshape(tuple(counts[::-1])) >= 128
    require(bool(mask.any()), "ROI is empty at this resolution: choose a finer grid.")
    return mask, corner, spacing


def ray_hits(mask, corner, spacing, source, directions, pulse=lambda: None):
    """Exact voxel traversal of the thresholded ROI for each BEV sample ray.

    Vectorized Amanatides-Woo traversal. Rays start at the source. All axes
    crossing together advance together, so corner-only contacts do not count.
    """
    shape = np.array(mask.shape[::-1])
    upper = corner + spacing * shape
    n = len(directions)
    parallel = np.abs(directions) < 1e-12
    outside = parallel & ((source < corner) | (source >= upper))
    safe_d = np.where(parallel, 1.0, directions)
    a = (corner - source) / safe_d
    b = (upper - source) / safe_d
    near = np.where(parallel, -np.inf, np.minimum(a, b))
    far = np.where(parallel, np.inf, np.maximum(a, b))
    enter = np.maximum(near.max(axis=1), 0)
    leave = far.min(axis=1)
    ids = np.flatnonzero((leave > enter) & ~outside.any(axis=1))
    result = np.zeros(n, dtype=bool)
    if not len(ids):
        return result
    d = directions[ids]
    pos = source + (enter[ids] + 1e-10)[:, None] * d
    idx = np.floor((pos - corner) / spacing).astype(int)
    idx = np.clip(idx, 0, shape - 1)
    step = np.sign(d).astype(int)
    boundary = corner + (idx + (step > 0)) * spacing
    with np.errstate(divide="ignore", invalid="ignore"):
        next_t = np.where(np.abs(d) < 1e-12, np.inf, (boundary - source) / d)
        delta = np.where(np.abs(d) < 1e-12, np.inf, spacing / np.abs(d))
    for iteration in range(int(shape.sum()) + 3):
        if iteration % 32 == 0:
            pulse()
        hit = mask[idx[:, 2], idx[:, 1], idx[:, 0]]
        result[ids[hit]] = True
        keep = ~hit
        ids, idx, step, next_t, delta = (
            arr[keep] for arr in (ids, idx, step, next_t, delta))
        if not len(ids):
            return result
        crossing_t = next_t.min(axis=1)
        axes = np.abs(next_t - crossing_t[:, None]) < 1e-10
        idx += axes * step
        next_t = np.where(axes, next_t + delta, next_t)
        keep = (crossing_t < leave[ids]) & ((idx >= 0) & (idx < shape)).all(axis=1)
        ids, idx, step, next_t, delta = (
            arr[keep] for arr in (ids, idx, step, next_t, delta))
        if not len(ids):
            return result
    raise CalculationError("ROI projection did not converge.")


def project_target(volume, iso, sad, frame, pulse=lambda: None):
    if isinstance(volume, Surface):
        return project_surface(volume, iso, sad, frame, pulse)
    # Retained reference path for independent numerical comparison.
    mask, corner, spacing = volume
    upper = corner + np.array(mask.shape[::-1]) * spacing
    corners = np.array(list(itertools.product(*zip(corner, upper))))
    local = (corners - iso) @ frame
    distance = sad - local[:, 2]
    require(np.all(distance > 1.0), "ROI is at or beyond the source plane.")
    bev = local[:, :2] * (sad / distance[:, None])
    lo = np.floor(bev.min(axis=0) / spacing).astype(int)
    hi = np.ceil(bev.max(axis=0) / spacing).astype(int)
    require(int(np.prod(hi - lo)) <= MAX_BEV_PIXELS,
            "Projection too large: choose a coarser resolution.")
    xx, yy = np.meshgrid((np.arange(lo[0], hi[0]) + .5) * spacing,
                         (np.arange(lo[1], hi[1]) + .5) * spacing)
    x, y = xx.ravel(), yy.ravel()
    source = iso + sad * frame[:, 2]
    directions = (x[:, None] * frame[:, 0] + y[:, None] * frame[:, 1]
                  - sad * frame[:, 2])
    hits = ray_hits(mask, corner, spacing, source, directions, pulse)
    require(bool(hits.any()), "Empty target projection: choose a finer grid.")
    return x[hits], y[hits]


class Surface:
    def __init__(self, faces, spacing, method):
        self.faces = np.asarray(faces, dtype=float)
        require(self.faces.ndim == 3 and self.faces.shape[1] in (3, 4)
                and self.faces.shape[2] == 3 and 0 < len(self.faces) <= MAX_SURFACE_FACES
                and np.isfinite(self.faces).all(), "Invalid target surface.")
        self.spacing = finite(spacing)
        require(self.spacing > 0, "Positive BEV spacing required.")
        self.method = method


def voxel_faces(volume, pulse=lambda: None, core_z=None):
    """Exact exposed voxel faces, merged only within a coplanar rectangle.

    This preserves the thresholded volume, holes and disconnected components.
    It does not smooth, decimate, dilate, or take a convex hull.
    """
    mask, corner, spacing = volume
    faces = []
    for axis in range(3):
        a = np.moveaxis(mask, axis, 0)
        changes = np.diff(np.pad(a.astype(np.int8), ((1, 1), (0, 0), (0, 0))), axis=0)
        other = [k for k in range(3) if k != axis]
        for plane in np.flatnonzero(np.any(changes != 0, axis=(1, 2))):
            pulse()
            if core_z is not None and axis == 0:
                first, last = core_z
                if plane < first or plane > last or (plane == last and last < mask.shape[0]):
                    continue
            boundary = changes[plane] != 0
            if core_z is not None and axis != 0:
                # Z is the first in-plane axis. Halo neighbours determine true
                # boundaries, but only the core contributes side faces.
                boundary[:core_z[0]] = False
                boundary[core_z[1]:] = False
            runs = np.diff(np.pad(boundary.astype(np.int8), ((0, 0), (1, 1))), axis=1)
            active = {}
            rectangles = []
            for row in range(boundary.shape[0] + 1):
                intervals = set() if row == boundary.shape[0] else set(zip(
                    np.flatnonzero(runs[row] == 1), np.flatnonzero(runs[row] == -1)))
                for interval in list(active):
                    if interval not in intervals:
                        rectangles.append((active.pop(interval), row, *interval))
                for interval in intervals:
                    active.setdefault(interval, row)
            for r0, r1, c0, c1 in rectangles:
                require(len(faces) < MAX_SURFACE_FACES, "Target surface too large: choose a coarser resolution.")
                face = np.zeros((4, 3))
                face[:, axis] = plane
                face[:, other[0]] = [r0, r0, r1, r1]
                face[:, other[1]] = [c0, c1, c1, c0]
                faces.append(corner + spacing * face[:, ::-1])
    return faces


def voxel_surface(volume, pulse=lambda: None):
    return Surface(voxel_faces(volume, pulse), volume[2], "Voxel surface")


def read_voxel_surface(geometry, spacing, pulse=lambda: None, timings=None):
    """Read bounded Z slabs at the requested resolution, with one-voxel halos.

    The former 16M limit now bounds each native read, not the whole ROI.
    Halo faces are excluded, so slab interfaces introduce no artificial caps.
    """
    timings = timings if timings is not None else Timings()
    with timings.measure("API"):
        lower, upper = [xyz(p) for p in geometry.GetBoundingBox()]
    require(np.all(upper > lower) and finite(spacing) > 0, "Invalid ROI bounds or spacing.")
    corner = lower-spacing
    counts = np.ceil((upper-lower)/spacing).astype(int)+2
    xy = int(counts[0])*int(counts[1])
    max_depth = MAX_VOXELS//xy
    require(max_depth >= 3, "ROI cross-section exceeds the voxel block limit at this resolution.")
    depth = int(counts[2]) if int(np.prod(counts)) <= MAX_VOXELS else max_depth-2
    timings.notes.append("ROI {}x{}x{} ({:.1f} M voxels), {:.1f} mm, {} block(s)".format(
        *counts, int(np.prod(counts))/1e6, spacing*10, math.ceil(int(counts[2])/depth)))
    faces = []
    for start in range(0, int(counts[2]), depth):
        pulse()
        stop = min(int(counts[2]), start+depth)
        read_start, read_stop = max(0, start-1), min(int(counts[2]), stop+1)
        block_corner = corner+np.array([0., 0., read_start*spacing])
        block_counts = np.array([counts[0], counts[1], read_stop-read_start])
        with timings.measure("API"):
            values = np.asarray(geometry.GetRoiGeometryAsVoxels(
                Corner=point_dict(block_corner), VoxelSize=point_dict([spacing]*3),
                NrVoxels=point_dict(block_counts, int)), dtype=np.uint8)
        require(values.size == int(np.prod(block_counts)), "RayStation returned an incomplete ROI block.")
        with timings.measure("Surface"):
            mask = values.reshape(tuple(block_counts[::-1])) >= 128
            del values
            faces.extend(voxel_faces((mask, block_corner, spacing), pulse, (start-read_start, stop-read_start)))
        require(len(faces) <= MAX_SURFACE_FACES, "Target surface exceeds the surface size limit.")
    require(bool(faces), "ROI is empty at this resolution: choose a finer grid.")
    with timings.measure("Surface"):
        surface = Surface(faces, spacing, "Voxel surface")
    timings.notes.append("{} surface faces".format(len(surface.faces)))
    return surface


def read_surface(geometry, spacing, pulse=lambda: None, timings=None):
    timings = timings if timings is not None else Timings()
    with timings.measure("API"):
        require(geometry.HasContours(), "Selected target has no geometry.")
    pulse()
    # PrimaryShape may be Contours/BinaryRoi instead. Probe only documented
    # read properties; never call SetRepresentation or create a temporary ROI.
    try:
        with timings.measure("API"):
            shape = geometry.PrimaryShape
            native_vertices, native_indices, closed = shape.Vertices, shape.Indices, shape.IsClosed
    except Exception:
        return read_voxel_surface(geometry, spacing, pulse, timings)
    require(bool(closed), "Native target mesh is not closed.")
    require(len(native_indices) <= 3*MAX_SURFACE_FACES, "Native target mesh exceeds the surface size limit.")
    vertices = []
    with timings.measure("API"):
        for i, vertex in enumerate(native_vertices):
            if i % 2048 == 0:
                pulse()
            vertices.append(xyz(vertex))
    vertices = np.array(vertices)
    indices = np.asarray(native_indices)
    require(vertices.ndim == 2 and vertices.shape[1] == 3
            and indices.ndim == 1 and indices.size >= 12 and indices.size % 3 == 0
            and np.issubdtype(indices.dtype, np.integer)
            and np.all((indices >= 0) & (indices < len(vertices))), "Invalid native mesh indices.")
    pulse()
    with timings.measure("Surface"):
        surface = Surface(vertices[indices.reshape(-1, 3)], spacing, "Native mesh")
    timings.notes.append("{} native triangles".format(len(surface.faces)))
    return surface


def project_surface(surface, iso, sad, frame, pulse=lambda: None):
    """Perspective surface projection with vectorized polygon scanlines.

    Each convex face contributes inclusive BEV pixel-centre intervals. Their
    union is accumulated with scanline differences, not a loop over 3D voxels.
    """
    step = surface.spacing
    local = (surface.faces - iso) @ frame
    distance = sad - local[:, :, 2]
    require(np.all(distance > 1.0), "ROI is at or beyond the source plane.")
    bev = local[:, :, :2] * (sad / distance[:, :, None])
    lo = np.floor(bev.min(axis=(0, 1)) / step).astype(int)
    hi = np.ceil(bev.max(axis=(0, 1)) / step).astype(int)
    width, height = hi - lo
    require(width > 0 and height > 0 and int(width)*int(height) <= MAX_BEV_PIXELS,
            "Projection too large: choose a coarser resolution.")
    polygons = bev / step - lo - .5
    differences = np.zeros((height, width + 1), dtype=np.int32)
    batch = max(1, min(512, 100_000 // int(height)))
    for start in range(0, len(polygons), batch):
        pulse()
        poly = polygons[start:start+batch]
        x, y = poly[:, :, 0], poly[:, :, 1]
        xn, yn = np.roll(x, -1, axis=1), np.roll(y, -1, axis=1)
        area = np.abs(np.sum(x*yn-xn*y, axis=1))
        keep = area > 1e-10
        x, y, xn, yn = x[keep], y[keep], xn[keep], yn[keep]
        if not len(x):
            continue
        low = np.maximum(0, np.ceil(y.min(axis=1)-1e-9).astype(int))
        high = np.minimum(height-1, np.floor(y.max(axis=1)+1e-9).astype(int))
        counts = np.maximum(0, high-low+1)
        ids = np.repeat(np.arange(len(x)), counts)
        rows = np.repeat(low, counts) + np.arange(counts.sum()) - np.repeat(np.cumsum(counts)-counts, counts)
        if not len(ids):
            continue
        yy = rows[:, None]
        dy = yn[ids]-y[ids]
        nonhorizontal = np.abs(dy) > 1e-12
        valid = nonhorizontal & (yy >= np.minimum(y[ids], yn[ids])-1e-9) & (yy <= np.maximum(y[ids], yn[ids])+1e-9)
        crossing = x[ids] + (yy-y[ids]) * (xn[ids]-x[ids]) / np.where(nonhorizontal, dy, 1)
        left = np.min(np.where(valid, crossing, np.inf), axis=1)
        right = np.max(np.where(valid, crossing, -np.inf), axis=1)
        horizontal = ~nonhorizontal & (np.abs(yy-y[ids]) <= 1e-9)
        left = np.minimum(left, np.min(np.where(horizontal, np.minimum(x[ids], xn[ids]), np.inf), axis=1))
        right = np.maximum(right, np.max(np.where(horizontal, np.maximum(x[ids], xn[ids]), -np.inf), axis=1))
        valid = np.isfinite(left) & np.isfinite(right)
        rows, left, right = rows[valid], left[valid], right[valid]
        first = np.maximum(0, np.ceil(left-1e-9).astype(int))
        last = np.minimum(width-1, np.floor(right+1e-9).astype(int))
        valid = first <= last
        np.add.at(differences, (rows[valid], first[valid]), 1)
        np.add.at(differences, (rows[valid], last[valid]+1), -1)
    pulse()
    occupied = np.cumsum(differences[:, :-1], axis=1) > 0
    j, i = np.nonzero(occupied)
    require(len(i) > 0, "Empty target projection: choose a finer grid.")
    return (i+lo[0]+.5)*step, (j+lo[1]+.5)*step


def prepare_aperture(x, y, layers, axis):
    """Map target samples to native leaf strips once per target projection."""
    require(axis in ("X", "Y"), "Unsupported MLC movement direction.")
    across = y if axis == "X" else x
    result = []
    for centres, widths in layers:
        low, high = centres-widths/2, centres+widths/2
        order = np.argsort(low)
        if np.all(high[order][:-1] <= low[order][1:]):
            strip = np.searchsorted(low[order], across, side="right")-1
            index = order[np.clip(strip, 0, len(order)-1)]
            valid = (strip >= 0) & (across < high[index])
            extras = (np.array([], dtype=int), np.array([], dtype=int))
        else:
            # Preserve exact union semantics even for numerical overlaps.
            membership = (across[None, :] >= low[:, None]) & (across[None, :] < high[:, None])
            index = membership.argmax(axis=0)
            valid = membership.any(axis=0)
            membership[index, np.arange(len(across))] = False
            extras = np.nonzero(membership)
        result.append((index, valid, extras))
    return result


def aperture_open(x, y, layers, banks, axis, jaws, lookup=None):
    require(axis in ("X", "Y"), "Unsupported MLC movement direction.")
    require(len(banks) == 2 * len(layers), "MLC layer/bank count mismatch.")
    moving = x if axis == "X" else y
    if lookup is None:
        lookup = prepare_aperture(x, y, layers, axis)
    opened = np.ones(len(x), dtype=bool)
    for i, (centres, widths) in enumerate(layers):
        left, right = [np.asarray(b, dtype=float) for b in banks[2*i:2*i+2]]
        require(left.shape == right.shape == centres.shape,
                "MLC leaf count does not match the machine geometry.")
        require(np.isfinite(left).all() and np.isfinite(right).all(), "Invalid MLC positions.")
        require(np.all(left <= right + 1e-8), "Crossed opposing MLC tips are unsupported.")
        index, valid, (extra_leaf, extra_point) = lookup[i]
        layer_open = valid & (moving >= left[index]) & (moving < right[index])
        if len(extra_point):
            np.logical_or.at(layer_open, extra_point,
                (moving[extra_point] >= left[extra_leaf]) & (moving[extra_point] < right[extra_leaf]))
        opened &= layer_open
    if jaws is not None:
        jaws = np.asarray(jaws, dtype=float)
        require(jaws.shape == (4,) and np.isfinite(jaws).all(), "Missing native jaw limits.")
        x1, x2, y1, y2 = jaws
        require(x1 <= x2 and y1 <= y2, "Invalid signed native jaw limits.")
        require(not np.all(np.abs(jaws) < 1e-8),
                "Zero jaw placeholders: jawless field limits need local verification.")
        opened &= (x >= x1) & (x < x2) & (y >= y1) & (y < y2)
    return opened


def weighted_mean(values, weights):
    v, w = np.array(values, dtype=float), np.array(weights, dtype=float)
    require(v.ndim == w.ndim == 1 and v.size == w.size and v.size > 0,
            "Incomplete beam or plan results.")
    require(np.isfinite(v).all() and np.isfinite(w).all() and np.all(w >= 0)
            and np.all((v >= 0) & (v <= 1)) and w.sum() > 0, "Invalid BAM/PAM weights.")
    return float(np.dot(v, w) / w.sum())


def model_for_beam(beam, machine_db):
    machine = machine_db.GetTreatmentMachine(
        machineName=beam.MachineReference.MachineName, lockMode="Read")
    require(machine is not None, "Treatment machine unavailable.")
    expected = beam.MachineReference.CommissioningTime
    actual = machine.CommissionTime
    require(expected is not None and actual is not None,
            "Machine commissioning version cannot be verified.")
    require(expected == actual, "Plan and loaded machine commissioning versions differ.")
    require(str(machine.PatientSupportType) == "Table"
            and not machine.ReplaceCouchRotationByRingRotation,
            "Unsupported treatment machine coordinate system.")
    require(str(machine.RoomViewModel) in ("SchematicLinac", "RingGantry"),
            "Only conventional linac / ring-gantry geometry is supported.")
    physics = machine.Physics
    sad = finite(physics.SourceAxisDistance)
    require(sad > 0, "Invalid source-axis distance.")
    mlc = physics.MlcPhysics
    require(mlc is not None and not mlc.IsApexAddOnMlc
            and str(mlc.LeafProjectionPlane) == "Isocenter",
            "An MLC model projected at isocenter is required.")
    layers = []
    for name in ("UpperLayer", "LowerLayer"):
        planned, commissioned = getattr(beam, name), getattr(mlc, name)
        if commissioned is None:
            require(planned is None or len(planned.LeafWidths) == 0,
                    "Plan has an unexpected MLC layer.")
            continue
        centres = np.array(commissioned.LeafCenterPositions, dtype=float)
        widths = np.array(commissioned.LeafWidths, dtype=float)
        require(planned is not None and centres.size > 0 and centres.shape == widths.shape,
                "Incomplete MLC layer geometry.")
        pc = np.array(planned.LeafCenterPositions, dtype=float)
        pw = np.array(planned.LeafWidths, dtype=float)
        require(pc.shape == centres.shape and pw.shape == widths.shape
                and np.allclose(pc, centres, rtol=0, atol=1e-6)
                and np.allclose(pw, widths, rtol=0, atol=1e-6),
                "Plan and commissioned MLC geometry differ.")
        require(np.isfinite(centres).all() and np.isfinite(widths).all()
                and np.all(widths > 0), "Invalid MLC strip geometry.")
        order = np.argsort(centres)
        c, w = centres[order], widths[order]
        require(np.all(c[:-1] + w[:-1]/2 <= c[1:] - w[1:]/2 + 1e-6),
                "Overlapping leaf strips within one MLC layer.")
        # Retain native order: LeafPositions uses the same leaf indexing.
        layers.append((centres, widths))
    require(len(layers) in (1, 2), "One or two complete MLC layers are required.")
    return sad, layers, str(mlc.MovementDirection), physics.JawPhysics is not None


def beam_samples(beam):
    require(str(beam.PatientPosition) == "HeadFirstSupine", "Only Head First Supine is supported.")
    require(str(beam.DeliveryTechnique) in ("SMLC", "DMLC", "DynamicArc", "StaticArc"),
            "Unsupported delivery technique (e.g. Tomo, CyberKnife or collapsed arc).")
    for name in ("CouchPitchAngle", "CouchRollAngle", "GimbalPanAngle", "GimbalTiltAngle"):
        require(abs(finite(getattr(beam, name))) < 1e-6,
                "Pitch, roll or gimbal rotation is unsupported.")
    for name in ("Cone", "Compensator", "Wedge"):
        require(getattr(beam, name) is None, "Cones, compensators and wedges are unsupported.")
    require(len(beam.Blocks) == 0, "Custom blocks are unsupported.")
    segments = list(beam.Segments)
    require(bool(segments), "Beam has no segments/control points.")
    weights = np.array([finite(s.RelativeWeight) for s in segments])
    require(np.all(weights >= 0) and abs(float(weights.sum()) - 1.0) <= 1e-6,
            "RelativeWeight must be nonnegative and sum to 1 (not cumulative MU).")
    arc = str(beam.ArcRotationDirection) in ("Clockwise", "CounterClockwise")
    gantry0 = finite(beam.GantryAngle)
    couch = finite(beam.CouchRotationAngle)
    angles = []
    for segment in segments:
        if arc:
            require(not segment.IsVirtual, "Virtual arc control points need local weight validation.")
            require(abs(finite(segment.DeltaCouchAngle)) < 1e-6,
                    "Couch motion during a beam is unsupported.")
        angles.append(gantry0 + (finite(segment.DeltaGantryAngle) if arc else 0))
    if arc:
        require(len(angles) >= 2 and abs(angles[0] - gantry0) < 1e-3,
                "Arc first control point does not match beam start.")
        require(beam.ArcStopGantryAngle is not None, "Missing arc stop angle.")
        error = (angles[-1] - finite(beam.ArcStopGantryAngle) + 180) % 360 - 180
        require(abs(error) < 1e-3, "Arc CP angles do not match the native stop angle.")
        sign = 1 if str(beam.ArcRotationDirection) == "Clockwise" else -1
        require(np.all(sign * np.diff(angles) >= -1e-6),
                "Arc delta-angle convention is inconsistent with rotation direction.")
    return segments, weights, angles, couch


def calculate_beam(beam, machine_db, volume, pulse=lambda: None, progress=lambda i, n: None, timings=None):
    timings = timings if timings is not None else Timings()
    with timings.measure("API"):
        segments, weights, angles, couch = beam_samples(beam)
        sad, layers, axis, has_jaws = model_for_beam(beam, machine_db)
        iso = xyz(beam.Isocenter.Position)
    modulation = []
    previous_frame, target, lookup = None, None, None
    for i, (segment, weight, gantry) in enumerate(zip(segments, weights, angles)):
        progress(i+1, len(segments))
        pulse()
        if weight == 0:
            modulation.append(0.0)
            continue
        with timings.measure("API"):
            collimator = finite(segment.CollimatorAngle)
            banks = [np.array(b, dtype=float) for b in segment.LeafPositions]
            jaws = segment.JawPositions if has_jaws else None
        frame = beam_frame(gantry, couch, collimator)
        if previous_frame is None or not np.allclose(frame, previous_frame, rtol=0, atol=1e-12):
            with timings.measure("Projection"):
                target = project_target(volume, iso, sad, frame, pulse)
            with timings.measure("MLC"):
                lookup = prepare_aperture(*target, layers, axis)
            previous_frame = frame
        with timings.measure("MLC"):
            opened = aperture_open(*target, layers, banks, axis, jaws, lookup)
            modulation.append(1.0 - float(np.count_nonzero(opened)) / len(opened))
    return weighted_mean(modulation, weights), len(segments)


class Viewer:
    def __init__(self, root, case, plan, machine_db, tk, ttk):
        self.root, self.case, self.plan, self.machine_db = root, case, plan, machine_db
        self.running = self.cancelled = False
        self.started = self.beam_started = self.stopped = None
        self.active_row = self.timer_job = None
        self.row_timings = {}
        self.last_pulse = 0.0
        root.title("PAM / BAM  " + VERSION)
        root.geometry("1000x640")
        root.minsize(800, 450)
        root.protocol("WM_DELETE_WINDOW", self.close)
        panel = ttk.Frame(root, padding=16)
        panel.pack(fill="both", expand=True)
        ttk.Label(panel, text="Plan Aperture Modulation", font=("Segoe UI", 18, "bold")).pack(anchor="w")
        ttk.Label(panel, text="Read-only research tool | PAM-BAM metric | numerical approximation").pack(anchor="w", pady=(2, 12))
        controls = ttk.Frame(panel)
        controls.pack(fill="x")
        rois = sorted(list(case.PatientModel.RegionsOfInterest), key=lambda r: str(r.Name).lower())
        names = [str(r.Name) for r in rois]
        require(bool(names), "No ROIs available in the current case.")
        default = next((str(r.Name) for r in rois if str(r.Type) == "Ptv"), names[0])
        self.roi = tk.StringVar(value=default)
        self.resolution = tk.StringVar(value="2.0")
        ttk.Label(controls, text="Target ROI").pack(side="left")
        self.roi_box = ttk.Combobox(controls, textvariable=self.roi, values=names, state="readonly", width=30)
        self.roi_box.pack(side="left", padx=(8, 16))
        ttk.Label(controls, text="Grid (mm)").pack(side="left")
        self.res_box = ttk.Combobox(controls, textvariable=self.resolution, values=["2.0", "1.0", "0.5"], state="readonly", width=5)
        self.res_box.pack(side="left", padx=8)
        self.calculate_button = ttk.Button(controls, text="Calculate", command=self.calculate)
        self.calculate_button.pack(side="left", padx=8)
        self.cancel_button = ttk.Button(controls, text="Cancel", command=self.cancel, state="disabled")
        self.cancel_button.pack(side="left")
        table_frame = ttk.Frame(panel)
        table_frame.pack(fill="both", expand=True, pady=(16, 8))
        columns = ("set", "beam", "mu", "cp", "bam", "seconds", "status")
        self.table = ttk.Treeview(table_frame, columns=columns, show="headings", height=10)
        for col, title, width in zip(columns, ("Beam set", "Beam", "MU", "CPs", "BAM", "Seconds", "Status"),
                                     (125, 100, 80, 55, 85, 75, 400)):
            self.table.heading(col, text=title)
            self.table.column(col, width=width, minwidth=45, stretch=col == "status")
        scroll = ttk.Scrollbar(table_frame, orient="vertical", command=self.table.yview)
        self.table.configure(yscrollcommand=scroll.set)
        scroll.pack(side="right", fill="y")
        self.table.pack(side="left", fill="both", expand=True)
        self.pam = tk.StringVar(value="Plan PAM: —")
        ttk.Label(panel, textvariable=self.pam, font=("Segoe UI", 17, "bold")).pack(anchor="w", pady=6)
        self.elapsed = tk.StringVar(value="Elapsed: 0.0 s")
        ttk.Label(panel, textvariable=self.elapsed).pack(anchor="w")
        self.status = tk.StringVar(value="Select the same target ROI for all treatment beams in the open plan.")
        ttk.Label(panel, textvariable=self.status, wraplength=940).pack(anchor="w", pady=4)
        self.detail = tk.StringVar(value="Select a beam row for API / surface / projection / MLC timings.")
        ttk.Label(panel, textvariable=self.detail, wraplength=940).pack(anchor="w", pady=2)
        self.table.bind("<<TreeviewSelect>>", self.show_detail)
        self.roi_box.bind("<<ComboboxSelected>>", self.invalidate)
        self.res_box.bind("<<ComboboxSelected>>", self.invalidate)
        ttk.Label(panel, text="0 = fully open target projection   •   1 = fully blocked target projection\n"
                  "MU-weighted geometry, not a dose or target-coverage score. HFS photon plans only.",
                  wraplength=940).pack(anchor="w", pady=(6, 0))

    def invalidate(self, event=None):
        self.table.delete(*self.table.get_children())
        self.row_timings.clear()
        self.pam.set("Plan PAM: —")
        self.status.set("Selection changed. Press Calculate to update the results.")
        self.detail.set("")
        self.elapsed.set("Elapsed: 0.0 s")

    def show_detail(self, event=None, row=None):
        rows = (row,) if row is not None else self.table.selection()
        if rows:
            values = self.table.item(rows[0], "values")
            timing = self.row_timings.get(rows[0])
            self.detail.set("{} / {}: {}{}".format(values[0], values[1], values[-1],
                "\n"+timing.summary() if timing is not None else ""))

    def cancel(self):
        self.cancelled = True

    def close(self):
        if self.running:
            self.cancel()
        else:
            self.root.destroy()

    def pulse(self):
        now = time.perf_counter()
        if now-self.last_pulse >= .05:
            self.last_pulse = now
            self.root.update()
        if self.cancelled:
            raise Cancelled()

    def tick(self):
        now = self.stopped if self.stopped is not None else time.perf_counter()
        if self.started is not None:
            text = "Elapsed: {:.1f} s".format(now-self.started)
            if self.beam_started is not None:
                seconds = now-self.beam_started
                text += " | Current beam: {:.1f} s".format(seconds)
                if self.active_row is not None:
                    self.table.set(self.active_row, "seconds", "{:.1f}".format(seconds))
            self.elapsed.set(text)
        if self.running:
            self.timer_job = self.root.after(100, self.tick)

    def calculate(self):
        if self.running:
            return
        self.running, self.cancelled = True, False
        self.started = time.perf_counter()
        self.beam_started = self.stopped = None
        self.active_row = None
        self.last_pulse = 0.0
        self.tick()
        self.calculate_button.configure(state="disabled")
        self.cancel_button.configure(state="normal")
        self.roi_box.configure(state="disabled")
        self.res_box.configure(state="disabled")
        self.table.delete(*self.table.get_children())
        self.row_timings.clear()
        self.detail.set("")
        self.pam.set("Plan PAM: calculating…")
        results, mus, exams = [], [], set()
        failed = total = excluded = 0
        volumes = {}
        roi_name, spacing = self.roi.get(), float(self.resolution.get()) / 10
        try:
            beamsets = list(self.plan.BeamSets)
            require(bool(beamsets), "The current plan has no beam sets.")
            for beamset in beamsets:
                for beam in beamset.Beams:
                    self.pulse()
                    label, name = str(beamset.DicomPlanLabel), str(beam.Name)
                    row = self.table.insert("", "end", values=(label, name, "—", "—", "—", "0.0", "Reading…"))
                    self.active_row, self.beam_started = row, time.perf_counter()
                    timings = Timings()
                    self.row_timings[row] = timings
                    mu = None
                    try:
                        if str(beam.DeliveryTechnique) == "Setup":
                            excluded += 1
                            self.table.item(row, values=(label, name, "—", "—", "—", "—", "Setup beam excluded"))
                            continue
                        mu = finite(beam.BeamMU)
                        require(mu >= 0, "Negative beam MU.")
                        if mu == 0:
                            excluded += 1
                            self.table.item(row, values=(label, name, "0", "—", "—", "—", "Zero-MU beam excluded"))
                            continue
                        total += 1
                        require(str(beamset.Modality) == "Photons", "Only photon beam sets are supported.")
                        require(not beamset.EnableDynamicTracking, "Dynamic tracking is unsupported.")
                        exam = beamset.GetPlanningExamination()
                        exams.add(str(exam.Name))
                        if str(exam.Name) not in volumes:
                            self.status.set("Reading the target and preparing its surface…")
                            self.pulse()
                            geometry = self.case.PatientModel.StructureSets[exam.Name].RoiGeometries[roi_name]
                            volumes[str(exam.Name)] = read_surface(geometry, spacing, self.pulse, timings)
                        else:
                            timings.notes.append("Target surface reused from this calculation")
                        def progress(i, n):
                            self.status.set("Beam {} | CP {}/{} | {} mm grid | {}".format(
                                name, i, n, self.resolution.get(), volumes[str(exam.Name)].method))
                        bam, ncp = calculate_beam(beam, self.machine_db, volumes[str(exam.Name)], self.pulse, progress, timings)
                        results.append(bam)
                        mus.append(mu)
                        self.table.item(row, values=(label, name, "{:.2f}".format(mu), ncp,
                                                    "{:.4f}".format(bam), "—", "OK ({})".format(volumes[str(exam.Name)].method)))
                    except Cancelled:
                        self.table.item(row, values=(label, name, "—", "—", "—", "—", "Cancelled"))
                        raise
                    except Exception as exc:
                        failed += 1
                        # Native exceptions can contain patient identifiers; never echo them.
                        message = str(exc) if isinstance(exc, CalculationError) else "RayStation API read failed ({})".format(type(exc).__name__)
                        self.table.item(row, values=(label, name, "—" if mu is None else "{:.2f}".format(mu), "—", "—", "—", message))
                    finally:
                        self.table.set(row, "seconds", "{:.1f}".format(time.perf_counter()-self.beam_started))
                        self.active_row = self.beam_started = None
                        self.show_detail(row=row)
                    self.table.see(row)
            require(failed == 0 and total > 0 and len(results) == total,
                    "Plan PAM unavailable: {} failed beam(s); {} complete beam(s).".format(failed, len(results)))
            require(len(exams) == 1, "Plan PAM unavailable: beam sets use different planning examinations.")
            pam = weighted_mean(results, mus)
            self.pam.set("Plan PAM: {:.4f}  ({:.2f}%)".format(pam, 100*pam))
            self.status.set("{} | {} treatment beams | {:.2f} MU | {} mm grid | {} | {} excluded beam(s).".format(
                roi_name, total, sum(mus), self.resolution.get(), ", ".join(sorted({v.method for v in volumes.values()})), excluded))
        except Cancelled:
            self.pam.set("Plan PAM: unavailable")
            self.status.set("Cancelled. Completed beam rows remain visible; no partial plan PAM is shown.")
        except Exception as exc:
            self.pam.set("Plan PAM: unavailable")
            self.status.set(str(exc) if isinstance(exc, CalculationError) else "RayStation read failed. Check the active plan and ROI.")
        finally:
            volumes.clear()
            self.running = False
            self.stopped = time.perf_counter()
            if self.timer_job is not None:
                self.root.after_cancel(self.timer_job)
                self.timer_job = None
            self.tick()
            self.calculate_button.configure(state="normal")
            self.cancel_button.configure(state="disabled")
            self.roi_box.configure(state="readonly")
            self.res_box.configure(state="readonly")


def main():
    import tkinter as tk
    from tkinter import ttk, messagebox
    root = tk.Tk()
    root.withdraw()
    try:
        try:
            from raystation import get_current
        except ImportError:
            from connect import get_current
        case, plan, machine_db = [get_current(k) for k in ("Case", "Plan", "MachineDB")]
        require(case is not None and plan is not None and machine_db is not None,
                "Open a case and treatment plan in RayStation first.")
        Viewer(root, case, plan, machine_db, tk, ttk)
        root.deiconify()
        root.mainloop()
    except Exception as exc:
        messagebox.showerror("PAM / BAM", str(exc) if isinstance(exc, CalculationError) else
                             "Run this file in RayStation CPython with an open plan.\n"
                             "Required: numpy and tkinter. API basis: v2025 SP2.", parent=root)
        root.destroy()


if __name__ == "__main__":
    main()
