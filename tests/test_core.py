import importlib.util
import math
import unittest


class CoreTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.available = importlib.util.find_spec('pambam') is not None
        if cls.available:
            import pambam
            cls.m = pambam

    def setUp(self):
        self.assertTrue(self.available, 'The standalone PAM/BAM library must exist')
        self.t = [(-10, -10, 10, 10)]
        self.half = [(-10, -10, 0, 10)]

    def test_single_layer_and_jaws(self):
        layer = self.m.Layer((-10, 0, 10), (-20, -10), (20, 10))
        a = self.m.single_layer_aperture(layer, jaws=(-15, -5, 15, 5))
        self.assertAlmostEqual(250, self.m.union_area(a))

    def test_dual_layer_staggered(self):
        a = self.m.Layer((-10, 0, 10), (-10, -10), (10, 10))
        b = self.m.Layer((-15, -5, 5, 15), (0, -5, 0), (5, 5, 5))
        aperture = self.m.dual_layer_aperture(a, b)
        self.assertAlmostEqual(150, self.m.union_area(aperture))
        self.assertAlmostEqual(.625, self.m.aperture_modulation(self.t, aperture))
        self.assertEqual(aperture, self.m.dual_layer_aperture(b, a))

    def test_union_not_sum_for_overlapping_targets(self):
        target = [(0, 0, 2, 2), (1, 0, 3, 2)]
        self.assertAlmostEqual(6, self.m.union_area(target))
        self.assertAlmostEqual(1 / 3, self.m.aperture_modulation(target, [(0, 0, 2, 2)]))

    def test_closed_and_open(self):
        self.assertEqual(1, self.m.aperture_modulation(self.t, []))
        self.assertEqual(0, self.m.aperture_modulation(self.t, self.t))

    def test_axis_y_and_fixed_limits(self):
        layer = self.m.Layer((-10, 10), (-20,), (20,), axis='Y')
        self.assertAlmostEqual(200, self.m.union_area(self.m.single_layer_aperture(layer, fixed_limits=(-5, -10, 5, 10))))

    def test_bam_endpoint_weights(self):
        r = self.m.bam_from_control_points([self.t]*3, [self.t, self.half, []], [0, .25, 1], 100, 'target')
        self.assertAlmostEqual(.625, r.bam)
        self.assertEqual((12.5, 50, 37.5), r.weights_mu)
        self.assertEqual(100, r.total_mu)

    def test_nonunit_final_weight(self):
        a = self.m.bam_from_control_points([self.t]*3, [self.t, self.half, []], [0, 25, 100], 100, 'target')
        self.assertAlmostEqual(.625, a.bam)

    def test_zero_mu_intervals(self):
        r = self.m.bam_from_control_points([self.t]*3, [[], self.t, self.t], [0, 0, 1], 100, 'target')
        self.assertEqual(0, r.bam)

    def test_pam_uses_mu_not_equal_beam_weights(self):
        a = self.m.calculate_bam([self.t], [self.half], [100], 'target')
        b = self.m.calculate_bam([self.t], [self.t], [300], 'target')
        self.assertAlmostEqual(.125, self.m.calculate_pam([a, b]))

    def test_invalid_and_mismatched_inputs(self):
        for w in ([0, -1], [0, math.nan], [1, 2], [0, 0]):
            with self.subTest(w=w), self.assertRaises(ValueError):
                self.m.bam_from_control_points([self.t]*2, [self.t]*2, w, 100, 'target')
        for mu in (0, -1, math.nan, math.inf):
            with self.subTest(mu=mu), self.assertRaises(ValueError):
                self.m.bam_from_control_points([self.t]*2, [self.t]*2, [0, 1], mu, 'target')
        with self.assertRaises(ValueError):
            self.m.calculate_bam([self.t], [], [1], 'target')
        with self.assertRaises(ValueError):
            self.m.aperture_modulation([], self.t)
        with self.assertRaises(ValueError):
            self.m.union_area([(0, 0, math.inf, 1)])
        with self.assertRaises(ValueError):
            self.m.single_layer_aperture(self.m.Layer((0, 0), (0,), (1,)))
        with self.assertRaises(ValueError):
            self.m.single_layer_aperture(self.m.Layer((0, 1), (2,), (1,)))
        with self.assertRaises(ValueError):
            self.m.dual_layer_aperture(self.m.Layer((0, 1), (0,), (1,)), None)
        with self.assertRaises(ValueError):
            self.m.calculate_pam([])

    def test_target_identity_must_match(self):
        a = self.m.calculate_bam([self.t], [self.half], [100], 'target-a')
        b = self.m.calculate_bam([self.t], [self.t], [100], 'target-b')
        with self.assertRaises(ValueError):
            self.m.calculate_pam([a, b])

    def test_inputs_not_mutated(self):
        import copy
        target = [self.t, self.t]
        aperture = [self.half, self.t]
        before = copy.deepcopy((target, aperture))
        self.m.bam_from_control_points(target, aperture, [0, 1], 100, 'target')
        self.assertEqual(before, (target, aperture))

    def test_shared_cross_language_fixtures(self):
        import json
        from pathlib import Path
        for case in json.loads(Path(__file__).with_name('geometry_cases.json').read_text()):
            with self.subTest(case=case['name']):
                layers = [self.m.Layer(x['edges'], x['a'], x['b'], x.get('axis','X')) for x in case['layers']]
                fn = self.m.single_layer_aperture if len(layers) == 1 else self.m.dual_layer_aperture
                aperture = fn(*layers, jaws=case.get('jaws'))
                self.assertAlmostEqual(case['area'], self.m.union_area(aperture))
                self.assertAlmostEqual(case['am'], self.m.aperture_modulation(case['target'], aperture))

    def test_closed_beam_roundoff_stays_in_metric_range(self):
        weights = [0.6430194504295635, 0.7875440123144711, 0.003711872322088694,
                   0.09476879162876783, 0.40262496843030005, 0.4084132074884518,
                   0.49099191979852774]
        result = self.m.calculate_bam([self.t]*7, [[]]*7, weights, 'target')
        self.assertEqual(1, self.m.calculate_pam([result]))


if __name__ == '__main__':
    unittest.main()
