"""Geometric aperture modulation; all BEV coordinates are isocenter-plane mm."""
from .core import (
    Layer, BeamResult, union_area, intersect, single_layer_aperture,
    dual_layer_aperture, aperture_modulation, calculate_bam,
    bam_from_control_points, calculate_pam,
)

__version__ = '0.1.0'
