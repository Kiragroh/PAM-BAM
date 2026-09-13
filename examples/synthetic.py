"""Run from the repository root: python -m examples.synthetic. No TPS or patient."""
from pambam import Layer, single_layer_aperture, dual_layer_aperture, bam_from_control_points, calculate_pam

target = [(-10, -10, 10, 10)]
upper = Layer((-10, 0, 10), (-10, -10), (10, 10))
lower = Layer((-15, -5, 5, 15), (0, -5, 0), (5, 5, 5))
single = single_layer_aperture(upper)
dual = dual_layer_aperture(upper, lower)
beam_a = bam_from_control_points([target]*2, [single]*2, [0, 1], 100, 'synthetic-target')
beam_b = bam_from_control_points([target]*2, [dual]*2, [0, 1], 300, 'synthetic-target')
print('Single-layer BAM:', beam_a.bam)
print('Dual-layer BAM:', beam_b.bam)
print('Combined synthetic PAM:', calculate_pam([beam_a, beam_b]))
