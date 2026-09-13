"""Frozen v1.2 scanline projector for regression comparisons, not a runtime dependency.

Source: Kiragroh/PAM-BAM revision 0c7bd90 (MIT, Maximilian Grohmann).
Only the projector is preserved; its arithmetic/tolerances are unchanged.
"""
import numpy as np
MAX_BEV_PIXELS = 600_000


def require(condition, message):
    if not condition:
        raise ValueError(message)

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
