"""Exact-mask regression against v1.2, using synthetic geometry only."""
import unittest
from unittest.mock import patch
import numpy as np
from test_pam_bam import p, beam, machine, volume, NS
from reference_projector_v12 import project_surface as reference


class ProjectionParityTests(unittest.TestCase):
    def test_masks_and_modulated_weights_match_v12(self):
        rng = np.random.default_rng(20240913)
        cube = np.ones((8, 9, 10), dtype=bool)
        hole = cube.copy(); hole[:, 3:6, 3:7] = False
        split = cube.copy(); split[:, :, 4:6] = False
        masks = [cube, hole, split] + [rng.random((7, 8, 9)) < .45 for _ in range(7)]
        angles = [(0, 0, 0), (90, 0, 90), (180, 0, 0), (270, 0, 0),
                  (0, 0, 1e-10), (37, 13, 27), (183, 21, 17), (273, 0, 3),
                  (65, 0, 90), (301, 55, 15), (350, 90, 3), (358, 0, 27)]
        layers = [(np.arange(-13.5, 14), np.ones(28)),
                  (np.arange(-14., 14.1), np.ones(29))]
        bams, expected_bams, mus = [], [], []
        comparisons = 0
        for index, mask in enumerate(masks):
            offset = rng.uniform(-150, 150, 3)
            quad = p.voxel_surface((mask, offset-1.5, .3)).faces
            for faces in (quad, np.concatenate((quad[:, [0, 1, 2]], quad[:, [0, 2, 3]]))):
                for spacing in (.2, .1, .05):
                    surface = p.Surface(faces, spacing, 'Synthetic')
                    iso = offset + rng.uniform(-.3, .3, 3)
                    ams, expected_ams = [], []
                    weights = rng.uniform(.01, 1, len(angles)); weights[3] = 0
                    for angle in angles:
                        frame = p.beam_frame(*angle)
                        actual = p.project_surface(surface, iso, 100, frame)
                        expected = reference(surface, iso, 100, frame)
                        np.testing.assert_array_equal(actual, expected)
                        banks = [a for c, _ in layers for a in
                                 np.sort(rng.uniform(-3, 3, (2, len(c))), axis=0)]
                        axis = 'X' if index % 2 else 'Y'
                        jaws = None if index % 3 else [-1.3, 2., -2., 1.8]
                        ams.append(1.-p.aperture_open(*actual, layers, banks, axis, jaws).mean())
                        expected_ams.append(1.-p.aperture_open(*expected, layers, banks, axis, jaws).mean())
                        comparisons += 1
                    bams.append(p.weighted_mean(ams, weights))
                    expected_bams.append(p.weighted_mean(expected_ams, weights))
                    mus.append(rng.uniform(10, 500))
        self.assertEqual(comparisons, 720)
        np.testing.assert_array_equal(bams, expected_bams)
        self.assertEqual(p.weighted_mean(bams, mus), p.weighted_mean(expected_bams, mus))

    def test_pixel_boundaries_and_nearly_horizontal_edges(self):
        for shift in (0., -1e-11, 1e-11, -1e-9, 1e-9):
            face = np.array([[[-.95, -.95, 0], [.95, -.95, 0],
                              [.95, .95, 0], [-.95, .95, 0]]]) + shift
            for tilt in (0., 1e-12, -1e-12, .3):
                face[0, 1, 1] += tilt
                surface = p.Surface(face, .1, 'Synthetic boundary')
                np.testing.assert_array_equal(
                    p.project_surface(surface, np.zeros(3), 100, np.eye(3)),
                    reference(surface, np.zeros(3), 100, np.eye(3)))

    def test_scanline_budget_splits_without_changing_masks(self):
        # Includes faces taller than the budget, zero-area faces and many
        # overlapping rectangles. Every sample still needs the exact union.
        faces = np.array([[[-1, -2, 0], [1, -2, 0], [1, 2, 0], [-1, 2, 0]]]*20)
        faces[1::3, :, 1] = 0
        surface = p.Surface(faces, .1, 'Synthetic tall faces')
        with patch.object(p, 'MAX_SCANLINE_ROWS', 3):
            actual = p.project_surface(surface, np.zeros(3), 100, np.eye(3))
        np.testing.assert_array_equal(actual, reference(surface, np.zeros(3), 100, np.eye(3)))


class CacheTests(unittest.TestCase):
    def test_same_view_reuses_target_but_recalculates_changed_leaves(self):
        surface = p.voxel_surface(volume())
        db = NS(GetTreatmentMachine=lambda **kw: machine())
        opened, closed = beam(), beam()
        closed.Segments[0].LeafPositions = [[0]*4, [0]*4]
        t1, t2 = p.Timings(), p.Timings()
        self.assertEqual(p.calculate_beam(opened, db, surface, timings=t1)[0], 0.)
        self.assertEqual(p.calculate_beam(closed, db, surface, timings=t2)[0], 1.)
        self.assertEqual(t1.calls['Projection'], 1)
        self.assertEqual(t2.calls['Projection'], 0)
        self.assertIn('1 cached projection(s)', t2.summary())

    def test_changed_geometry_cannot_reuse_projection(self):
        surface = p.voxel_surface(volume())
        iso, frame = np.zeros(3), p.beam_frame(0, 0, 0)
        timing = p.Timings()
        p.cached_target(surface, iso, 100, frame, lambda: None, timing)
        _, hit = p.cached_target(surface, iso.copy(), 100, frame.copy(), lambda: None, timing)
        self.assertTrue(hit)
        for pos, sad, orientation in [(iso+.03, 100, frame), (iso, 90, frame),
                                      (iso, 100, p.beam_frame(37, 0, 0))]:
            result, hit = p.cached_target(surface, pos, sad, orientation, lambda: None, timing)
            self.assertFalse(hit)
            np.testing.assert_array_equal(result, reference(surface, pos, sad, orientation))
        surface.spacing = .05
        self.assertFalse(p.cached_target(surface, iso, 100, frame, lambda: None, timing)[1])
        fresh = p.voxel_surface(volume())
        self.assertFalse(p.cached_target(fresh, iso, 100, frame, lambda: None, timing)[1])

    def test_cache_eviction_and_zero_budget(self):
        surface = p.voxel_surface(volume())
        iso, timing = np.zeros(3), p.Timings()
        with patch.object(p, 'MAX_PROJECTION_CACHE_BYTES', 12000):
            for angle in (0, 31, 90):
                p.cached_target(surface, iso, 100, p.beam_frame(angle, 0, 0), lambda: None, timing)
                self.assertLessEqual(surface.projection_bytes, 12000)
            self.assertLess(len(surface.projections), 3)
            self.assertFalse(p.cached_target(surface, iso, 100, p.beam_frame(0, 0, 0), lambda: None, timing)[1])
        fresh = p.voxel_surface(volume())
        with patch.object(p, 'MAX_PROJECTION_CACHE_BYTES', 0):
            for _ in range(2):
                self.assertFalse(p.cached_target(fresh, iso, 100, np.eye(3), lambda: None, timing)[1])
        self.assertEqual(fresh.projection_bytes, 0)
        self.assertFalse(fresh.projections)


if __name__ == '__main__':
    unittest.main(verbosity=2)
