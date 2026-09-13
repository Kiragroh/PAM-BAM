"""Read-only adapter boundary for RayStation photon segments.

RayStation v2025 SP2 API documents LeafPositions bank order as upper-left,
upper-right, lower-left, lower-right. Native lengths are cm; core lengths are mm.
Local boundary arrays, coordinate alignment and weight semantics need commissioning.
"""
import math
from .core import (Layer, _finite, _rectangles, single_layer_aperture, dual_layer_aperture,
                   bam_from_control_points, calculate_bam)


def aperture_from_segment(segment, layer_boundaries_mm, *, jaw_mode, axis='X'):
    boundaries = tuple(tuple(e) for e in layer_boundaries_mm)
    banks = tuple(tuple(_finite(v)*10 for v in bank) for bank in segment.LeafPositions)
    if len(boundaries) not in (1, 2) or len(banks) != 2*len(boundaries):
        raise ValueError('Complete explicit single- or dual-layer bank mapping required')
    layers = [Layer(edges, banks[2*i], banks[2*i+1], axis) for i, edges in enumerate(boundaries)]
    if jaw_mode not in ('none', 'physical', 'fixed'):
        raise ValueError('Specify none, physical, or fixed jaw mode; no hardware inference')
    jaw = None
    if jaw_mode != 'none':
        p = tuple(_finite(v)*10 for v in segment.JawPositions)
        if len(p) != 4:
            raise ValueError('Complete native jaw limits required')
        # RayStation order X1, X2, Y1, Y2 -> rectangle X1, Y1, X2, Y2.
        jaw = (p[0], p[2], p[1], p[3])
    limits = {'jaws': jaw if jaw_mode == 'physical' else None,
              'fixed_limits': jaw if jaw_mode == 'fixed' else None}
    return single_layer_aperture(layers[0], **limits) if len(layers) == 1 else dual_layer_aperture(*layers, **limits)


def calculate_beam(beam, layer_boundaries_mm, target_at_segment, target_key, *,
                   jaw_mode, axis='X', cumulative_weights=None, sample_weights_mu=None):
    """Read one photon beam without saving/modifying it.

    target_at_segment(segment,index) supplies a validated projected target union
    in the beam-limiting-device isocenter plane, in mm, at EACH segment/CP.
    Pass either a validated cumulative control-point sequence (VMAT endpoint
    method) or explicit sample MU (e.g. static segments). RelativeWeight is NOT
    assumed to mean cumulative meterset or endpoint MU in every RS version.
    """
    if (cumulative_weights is None) == (sample_weights_mu is None):
        raise ValueError('Specify exactly one explicit control-point or sample weight basis')
    segments = tuple(beam.Segments)
    mu = _finite(beam.BeamMU)
    if mu <= 0 or not segments:
        raise ValueError('A positive-MU treatment beam with segments is required')
    apertures = [aperture_from_segment(s, layer_boundaries_mm, jaw_mode=jaw_mode, axis=axis) for s in segments]
    # Snapshot inner coordinate buffers immediately: callbacks may reuse them.
    targets = [_rectangles(target_at_segment(s, i)) for i, s in enumerate(segments)]
    if cumulative_weights is not None:
        return bam_from_control_points(targets, apertures, cumulative_weights, mu, target_key)
    result = calculate_bam(targets, apertures, sample_weights_mu, target_key, 'explicit-raystation-sample-mu')
    if not math.isclose(result.total_mu, mu, rel_tol=1e-8, abs_tol=1e-8):
        raise ValueError('Sample MU must sum to the native beam MU; no silent normalization')
    return result
