"""Independent geometric implementation of Hernandez et al., doi:10.1002/mp.70144.

Targets are unions of BEV rectangles (e.g. a supplied silhouette raster), not
bounding boxes inferred from structure names. Finite geometry is required.
"""
from dataclasses import dataclass
import math


def _finite(value):
    value = float(value)
    if not math.isfinite(value):
        raise ValueError('Finite numeric input required')
    return value


def _rectangles(items):
    if items is None:
        raise ValueError('Missing geometry is not a closed aperture')
    result = []
    for item in items:
        if len(item) != 4:
            raise ValueError('Rectangle requires x1, y1, x2, y2')
        x1, y1, x2, y2 = map(_finite, item)
        if x1 > x2 or y1 > y2:
            raise ValueError('Reversed rectangle bounds')
        if x2 > x1 and y2 > y1:
            result.append((x1, y1, x2, y2))
    return tuple(result)


@dataclass(frozen=True)
class Layer:
    boundaries_mm: tuple
    bank_a_mm: tuple
    bank_b_mm: tuple
    axis: str = 'X'


@dataclass(frozen=True)
class BeamResult:
    bam: float
    total_mu: float
    target_key: str
    am: tuple
    weights_mu: tuple
    sampling: str


def _strips(layer):
    if not isinstance(layer, Layer):
        raise ValueError('Every physical MLC layer must be supplied')
    edges = tuple(map(_finite, layer.boundaries_mm))
    left = tuple(map(_finite, layer.bank_a_mm))
    right = tuple(map(_finite, layer.bank_b_mm))
    if (not left or len(edges) != len(left) + 1 or len(right) != len(left)
            or layer.axis not in ('X', 'Y')
            or any(b <= a for a, b in zip(edges, edges[1:]))
            or any(a > b for a, b in zip(left, right))):
        raise ValueError('Invalid leaf boundaries, banks, or travel axis')
    return _rectangles(
        (left[i], edges[i], right[i], edges[i+1]) if layer.axis == 'X'
        else (edges[i], left[i], edges[i+1], right[i])
        for i in range(len(left)))


def union_area(rectangles):
    """Exact union area of axis-aligned rectangles, without overlap double counting."""
    rects = _rectangles(rectangles)
    xs = sorted({r[i] for r in rects for i in (0, 2)})
    area = 0.0
    for x1, x2 in zip(xs, xs[1:]):
        # Each open slab has constant y-interval membership.
        intervals = sorted((r[1], r[3]) for r in rects if r[0] <= x1 and r[2] >= x2)
        length = 0.0
        end = -math.inf
        for low, high in intervals:
            length += max(0.0, high - max(low, end))
            end = max(end, high)
        area += (x2 - x1) * length
    return _finite(area)


def intersect(first, second):
    a, b = _rectangles(first), _rectangles(second)
    result = []
    for r in a:
        for s in b:
            x1, y1 = max(r[0], s[0]), max(r[1], s[1])
            x2, y2 = min(r[2], s[2]), min(r[3], s[3])
            if x2 > x1 and y2 > y1:
                result.append((x1, y1, x2, y2))
    return tuple(sorted(set(result)))


def _limit(openings, jaws, fixed_limits):
    if jaws is not None:
        openings = intersect(openings, [jaws])
    if fixed_limits is not None:
        openings = intersect(openings, [fixed_limits])
    return openings


def single_layer_aperture(layer, jaws=None, fixed_limits=None):
    """Union of physical leaf strips, clipped only by explicitly supplied limits."""
    return _limit(_strips(layer), jaws, fixed_limits)


def dual_layer_aperture(proximal, distal, jaws=None, fixed_limits=None):
    """Effective opening = proximal intersection distal, retaining the layer offsets.

    For a jawless Halcyon pass jaws=None; do not substitute virtual jaws for a layer.
    Both layers must already be projected into the same isocenter plane.
    """
    return _limit(intersect(_strips(proximal), _strips(distal)), jaws, fixed_limits)


def aperture_modulation(target, aperture):
    """AM = 1 - area(target intersect opening) / area(target)."""
    target, aperture = _rectangles(target), _rectangles(aperture)
    total = union_area(target)
    if total <= 0:
        raise ValueError('Target projection must have positive area')
    return min(1.0, max(0.0, 1.0 - union_area(intersect(target, aperture)) / total))


def calculate_bam(targets, apertures, weights_mu, target_key, sampling='weighted-samples'):
    """BAM of one beam from explicit per-sample MU weights (not cumulative weights)."""
    targets, apertures = tuple(targets), tuple(apertures)
    weights = tuple(map(_finite, weights_mu))
    if (not weights or len(targets) != len(weights) or len(apertures) != len(weights)
            or any(w < 0 for w in weights) or not isinstance(target_key, str)
            or not target_key.strip() or not sampling):
        raise ValueError('Complete matched samples, nonnegative MU, and target key required')
    total = _finite(math.fsum(weights))
    if total <= 0:
        raise ValueError('Beam MU must be positive; omit setup beams explicitly')
    am = tuple(aperture_modulation(t, a) for t, a in zip(targets, apertures))
    # All inputs were validated in [0,1]; only bound final arithmetic roundoff.
    bam = min(1.0, max(0.0, math.fsum(a * (w / total) for a, w in zip(am, weights))))
    return BeamResult(bam, total, target_key, am, weights, sampling)


def bam_from_control_points(targets, apertures, cumulative_weights, total_mu, target_key):
    """Trapezoidal MU integration: half each interval MU to either endpoint.

    No delivered timing is inferred. The final cumulative weight need not equal 1.
    The geometric endpoint approximation is not continuous leaf-motion integration.
    """
    cmw = tuple(map(_finite, cumulative_weights))
    mu = _finite(total_mu)
    if (len(cmw) < 2 or cmw[0] != 0 or cmw[-1] <= 0 or mu <= 0
            or any(b < a for a, b in zip(cmw, cmw[1:]))):
        raise ValueError('Ordered cumulative weights starting at zero and positive MU required')
    weights = [0.0] * len(cmw)
    for i in range(1, len(cmw)):
        interval = mu * ((cmw[i] - cmw[i-1]) / cmw[-1])
        weights[i-1] += interval / 2
        weights[i] += interval / 2
    return calculate_bam(targets, apertures, weights, target_key, 'control-point-endpoint-trapezoid')


def calculate_pam(beams):
    """Original plan definition: sum(BAM_b * MU_b) / sum(MU_b), one common target."""
    beams = tuple(beams)
    if not beams:
        raise ValueError('At least one complete treatment beam is required')
    for beam in beams:
        if (not isinstance(beam, BeamResult) or not 0 <= _finite(beam.bam) <= 1
                or _finite(beam.total_mu) <= 0 or not beam.target_key
                or beam.target_key != beams[0].target_key):
            raise ValueError('Complete BAM values with positive MU and identical target required')
    total = _finite(math.fsum(b.total_mu for b in beams))
    return min(1.0, max(0.0, math.fsum(b.bam * (b.total_mu / total) for b in beams)))
