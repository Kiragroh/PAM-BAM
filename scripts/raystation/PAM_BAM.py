# -*- coding: utf-8 -*-
"""PAM / BAM viewer 1.0 -- standalone, read-only RayStation file script.

Run in RayStation CPython with an open plan; select a target ROI and Calculate.
Dependencies: numpy, tkinter (no installation or repository import at runtime).
API basis: RayStation v2025 SP2 / 17.2.0. Local validation is still required.

AM = fraction of the projected target outside the MLC/jaw opening.
BAM = sum(AM * Segment.RelativeWeight); PAM = MU-weighted mean of beam BAMs.
RelativeWeight is documented by RS as a fraction of total beam MU, NOT a
cumulative meterset. Dynamic delivery is sampled at the native segments/CPs;
there is no reconstruction of continuous leaf travel or dose calculation.

The entire ROI is voxelized on the planning examination (>=128/255 occupancy).
Perspective ray tracing projects that binary volume onto the isocenter plane.
Both ROI voxels and BEV sampling use the selected spacing. Values therefore
depend on resolution; compare 1.0 and 0.5 mm when assessing convergence.
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
import numpy as np

VERSION = "1.0"
MAX_VOXELS = 16_000_000
MAX_BEV_PIXELS = 600_000


class CalculationError(Exception):
    """An actionable, patient-independent message safe to display."""


class Cancelled(Exception):
    pass


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


def aperture_open(x, y, layers, banks, axis, jaws):
    require(axis in ("X", "Y"), "Unsupported MLC movement direction.")
    require(len(banks) == 2 * len(layers), "MLC layer/bank count mismatch.")
    moving, across = (x, y) if axis == "X" else (y, x)
    opened = np.ones(len(x), dtype=bool)
    for i, (centres, widths) in enumerate(layers):
        left, right = [np.asarray(b, dtype=float) for b in banks[2*i:2*i+2]]
        require(left.shape == right.shape == centres.shape,
                "MLC leaf count does not match the machine geometry.")
        require(np.isfinite(left).all() and np.isfinite(right).all(), "Invalid MLC positions.")
        require(np.all(left <= right + 1e-8), "Crossed opposing MLC tips are unsupported.")
        layer_open = np.zeros(len(x), dtype=bool)
        for c, width, l, r in zip(centres, widths, left, right):
            layer_open |= ((across >= c-width/2) & (across < c+width/2)
                           & (moving >= l) & (moving < r))
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


def calculate_beam(beam, machine_db, volume, pulse=lambda: None, progress=lambda i, n: None):
    segments, weights, angles, couch = beam_samples(beam)
    sad, layers, axis, has_jaws = model_for_beam(beam, machine_db)
    iso = xyz(beam.Isocenter.Position)
    modulation = []
    previous_frame, target = None, None
    for i, (segment, weight, gantry) in enumerate(zip(segments, weights, angles)):
        progress(i+1, len(segments))
        pulse()
        if weight == 0:
            modulation.append(0.0)
            continue
        frame = beam_frame(gantry, couch, finite(segment.CollimatorAngle))
        if previous_frame is None or not np.allclose(frame, previous_frame, rtol=0, atol=1e-12):
            target = project_target(volume, iso, sad, frame, pulse)
            previous_frame = frame
        banks = [np.array(b, dtype=float) for b in segment.LeafPositions]
        opened = aperture_open(*target, layers, banks, axis,
                               segment.JawPositions if has_jaws else None)
        modulation.append(1.0 - float(np.count_nonzero(opened)) / len(opened))
    return weighted_mean(modulation, weights), len(segments)


class Viewer:
    def __init__(self, root, case, plan, machine_db, tk, ttk):
        self.root, self.case, self.plan, self.machine_db = root, case, plan, machine_db
        self.running = self.cancelled = False
        root.title("PAM / BAM  " + VERSION)
        root.geometry("1000x560")
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
        self.resolution = tk.StringVar(value="1.0")
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
        columns = ("set", "beam", "mu", "cp", "bam", "status")
        self.table = ttk.Treeview(table_frame, columns=columns, show="headings", height=10)
        for col, title, width in zip(columns, ("Beam set", "Beam", "MU", "CPs", "BAM", "Status"),
                                     (125, 130, 80, 55, 105, 430)):
            self.table.heading(col, text=title)
            self.table.column(col, width=width, minwidth=45, stretch=col == "status")
        scroll = ttk.Scrollbar(table_frame, orient="vertical", command=self.table.yview)
        self.table.configure(yscrollcommand=scroll.set)
        scroll.pack(side="right", fill="y")
        self.table.pack(side="left", fill="both", expand=True)
        self.pam = tk.StringVar(value="Plan PAM: —")
        ttk.Label(panel, textvariable=self.pam, font=("Segoe UI", 17, "bold")).pack(anchor="w", pady=6)
        self.status = tk.StringVar(value="Select the same target ROI for all treatment beams in the open plan.")
        ttk.Label(panel, textvariable=self.status, wraplength=940).pack(anchor="w", pady=4)
        self.detail = tk.StringVar(value="Select a beam row to read its full status.")
        ttk.Label(panel, textvariable=self.detail, wraplength=940).pack(anchor="w", pady=2)
        self.table.bind("<<TreeviewSelect>>", self.show_detail)
        self.roi_box.bind("<<ComboboxSelected>>", self.invalidate)
        self.res_box.bind("<<ComboboxSelected>>", self.invalidate)
        ttk.Label(panel, text="0 = fully open target projection   •   1 = fully blocked target projection\n"
                  "MU-weighted geometry, not a dose or target-coverage score. HFS photon plans only.",
                  wraplength=940).pack(anchor="w", pady=(6, 0))

    def invalidate(self, event=None):
        self.table.delete(*self.table.get_children())
        self.pam.set("Plan PAM: —")
        self.status.set("Selection changed. Press Calculate to update the results.")
        self.detail.set("")

    def show_detail(self, event=None):
        rows = self.table.selection()
        if rows:
            values = self.table.item(rows[0], "values")
            self.detail.set("{} / {}: {}".format(values[0], values[1], values[-1]))

    def cancel(self):
        self.cancelled = True

    def close(self):
        if self.running:
            self.cancel()
        else:
            self.root.destroy()

    def pulse(self):
        self.root.update()
        if self.cancelled:
            raise Cancelled()

    def calculate(self):
        if self.running:
            return
        self.running, self.cancelled = True, False
        self.calculate_button.configure(state="disabled")
        self.cancel_button.configure(state="normal")
        self.roi_box.configure(state="disabled")
        self.res_box.configure(state="disabled")
        self.table.delete(*self.table.get_children())
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
                    row = self.table.insert("", "end", values=(label, name, "—", "—", "—", "Reading…"))
                    mu = None
                    try:
                        if str(beam.DeliveryTechnique) == "Setup":
                            excluded += 1
                            self.table.item(row, values=(label, name, "—", "—", "—", "Setup beam excluded"))
                            continue
                        mu = finite(beam.BeamMU)
                        require(mu >= 0, "Negative beam MU.")
                        if mu == 0:
                            excluded += 1
                            self.table.item(row, values=(label, name, "0", "—", "—", "Zero-MU beam excluded"))
                            continue
                        total += 1
                        require(str(beamset.Modality) == "Photons", "Only photon beam sets are supported.")
                        require(not beamset.EnableDynamicTracking, "Dynamic tracking is unsupported.")
                        exam = beamset.GetPlanningExamination()
                        exams.add(str(exam.Name))
                        if str(exam.Name) not in volumes:
                            self.status.set("Reading and voxelizing the selected target ROI…")
                            self.pulse()
                            geometry = self.case.PatientModel.StructureSets[exam.Name].RoiGeometries[roi_name]
                            volumes[str(exam.Name)] = voxelize(geometry, spacing)
                        def progress(i, n):
                            self.status.set("Beam {} | CP {}/{} | {} mm grid".format(name, i, n, self.resolution.get()))
                        bam, ncp = calculate_beam(beam, self.machine_db, volumes[str(exam.Name)], self.pulse, progress)
                        results.append(bam)
                        mus.append(mu)
                        self.table.item(row, values=(label, name, "{:.2f}".format(mu), ncp,
                                                    "{:.4f}".format(bam), "OK (sampled geometry)"))
                    except Cancelled:
                        self.table.item(row, values=(label, name, "—", "—", "—", "Cancelled"))
                        raise
                    except Exception as exc:
                        failed += 1
                        # Native exceptions can contain patient identifiers; never echo them.
                        message = str(exc) if isinstance(exc, CalculationError) else "RayStation API read failed ({})".format(type(exc).__name__)
                        self.table.item(row, values=(label, name, "—" if mu is None else "{:.2f}".format(mu), "—", "—", message))
                    self.table.see(row)
            require(failed == 0 and total > 0 and len(results) == total,
                    "Plan PAM unavailable: {} failed beam(s); {} complete beam(s).".format(failed, len(results)))
            require(len(exams) == 1, "Plan PAM unavailable: beam sets use different planning examinations.")
            pam = weighted_mean(results, mus)
            self.pam.set("Plan PAM: {:.4f}  ({:.2f}%)".format(pam, 100*pam))
            self.status.set("{} | {} treatment beams | {:.2f} MU | {} mm grid | {} excluded beam(s).".format(
                roi_name, total, sum(mus), self.resolution.get(), excluded))
        except Cancelled:
            self.pam.set("Plan PAM: unavailable")
            self.status.set("Cancelled. Completed beam rows remain visible; no partial plan PAM is shown.")
        except Exception as exc:
            self.pam.set("Plan PAM: unavailable")
            self.status.set(str(exc) if isinstance(exc, CalculationError) else "RayStation read failed. Check the active plan and ROI.")
        finally:
            volumes.clear()
            self.running = False
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
