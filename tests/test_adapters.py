import importlib.util
from types import SimpleNamespace
import unittest


class RayStationAdapterTests(unittest.TestCase):
    def setUp(self):
        self.assertIsNotNone(importlib.util.find_spec('pambam.raystation'), 'RayStation adapter required')
        from pambam import raystation
        self.rs = raystation

    def test_single_cm_to_mm_and_jaw_order(self):
        from pambam import union_area
        seg = SimpleNamespace(LeafPositions=[[-2, -1], [2, 1]], JawPositions=[-1.5, 1.5, -.5, .5])
        self.assertAlmostEqual(250, union_area(self.rs.aperture_from_segment(seg, [(-10, 0, 10)], jaw_mode='physical')))

    def test_four_bank_dual_layer_and_no_jaws(self):
        from pambam import union_area
        # No JawPositions attribute: jawless mode must not even read it.
        seg = SimpleNamespace(LeafPositions=[[-1, -1], [1, 1], [0, -.5, 0], [.5, .5, .5]])
        result = self.rs.aperture_from_segment(seg, [(-10, 0, 10), (-15, -5, 5, 15)], jaw_mode='none')
        self.assertAlmostEqual(150, union_area(result))
        with self.assertRaises(ValueError):
            self.rs.aperture_from_segment(seg, [(-10, 0, 10)], jaw_mode='none')

    def test_explicit_cumulative_weights_not_guessed(self):
        segs = [SimpleNamespace(LeafPositions=[[-1], [v]]) for v in (1, 0)]
        beam = SimpleNamespace(BeamMU=100, Segments=segs)
        t = lambda segment, index: [(-10, -10, 10, 10)]
        r = self.rs.calculate_beam(beam, [(-10,10)], t, 'target',
            jaw_mode='none', cumulative_weights=[0, 1])
        self.assertAlmostEqual(.25, r.bam)
        with self.assertRaises(ValueError):
            self.rs.calculate_beam(beam, [(-10,10)], t, 'target', jaw_mode='none')

    def test_explicit_sample_weights_must_conserve_mu(self):
        beam = SimpleNamespace(BeamMU=100, Segments=[SimpleNamespace(LeafPositions=[[-1], [0]])])
        args = (beam, [(-10,10)], lambda s,i: [(-10,-10,10,10)], 'target')
        r = self.rs.calculate_beam(*args, jaw_mode='none', sample_weights_mu=[100])
        self.assertAlmostEqual(.5, r.bam)
        with self.assertRaises(ValueError):
            self.rs.calculate_beam(*args, jaw_mode='none', sample_weights_mu=[80])

    def test_unknown_jaw_mode_is_rejected(self):
        seg = SimpleNamespace(LeafPositions=[[-1], [1]])
        with self.assertRaises(ValueError):
            self.rs.aperture_from_segment(seg, [(-10,10)], jaw_mode='auto')

    def test_target_provider_buffer_is_snapshotted_immediately(self):
        beam = SimpleNamespace(BeamMU=100, Segments=[SimpleNamespace(LeafPositions=[[-1], [1]])]*2)
        buffer = [[-10,-10,10,10]]
        def provider(segment, index):
            buffer[0][:] = [-10,-10,10,10] if index == 0 else [30,-10,50,10]
            return buffer
        result = self.rs.calculate_beam(beam, [(-10,10)], provider, 'target',
            jaw_mode='none', cumulative_weights=[0,1])
        self.assertAlmostEqual(.5, result.bam)
