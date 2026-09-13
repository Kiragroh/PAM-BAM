"""Synthetic verification only; no RayStation session or patient data."""
import importlib.util
from pathlib import Path
from types import SimpleNamespace as NS
import unittest
import time
from unittest.mock import patch
import numpy as np

spec = importlib.util.spec_from_file_location("pam_bam", Path(__file__).resolve().parents[1] / "PAM_BAM.py")
p = importlib.util.module_from_spec(spec)
spec.loader.exec_module(p)


def layer():
    return NS(LeafCenterPositions=np.array([-1.5, -.5, .5, 1.5]), LeafWidths=np.ones(4))


def machine():
    return NS(CommissionTime=1, PatientSupportType="Table", ReplaceCouchRotationByRingRotation=False,
              RoomViewModel="SchematicLinac", Physics=NS(SourceAxisDistance=100,
              MlcPhysics=NS(IsApexAddOnMlc=False, LeafProjectionPlane="Isocenter",
                            UpperLayer=layer(), LowerLayer=None, MovementDirection="X"), JawPhysics=NS()))


def beam():
    b = NS(PatientPosition="HeadFirstSupine", DeliveryTechnique="SMLC", Name="Synthetic",
           CouchPitchAngle=0, CouchRollAngle=0, GimbalPanAngle=0, GimbalTiltAngle=0,
           Cone=None, Compensator=None, Wedge=None, Blocks=[], GantryAngle=0,
           CouchRotationAngle=0, ArcRotationDirection="None", BeamMU=100,
           Isocenter=NS(Position=dict(x=0,y=0,z=0)), UpperLayer=layer(), LowerLayer=None,
           MachineReference=NS(MachineName="Synthetic", CommissioningTime=1))
    b.Segments = [NS(RelativeWeight=1, CollimatorAngle=0, JawPositions=[-2,2,-2,2],
                     LeafPositions=[[-2]*4,[2]*4])]
    return b


def volume():
    # Full cube: DICOM [-1,1]^3 cm, sampled at 1 mm.
    return np.ones((20,20,20), dtype=bool), np.array([-1.,-1.,-1.]), .1


class ProjectionTests(unittest.TestCase):
    def test_cardinal_dicom_frames(self):
        expected = [(0,[0,-1,0]), (90,[1,0,0]), (180,[0,1,0]), (270,[-1,0,0])]
        for angle, source_axis in expected:
            f = p.beam_frame(angle, 0, 0)
            np.testing.assert_allclose(f[:,2], source_axis, atol=1e-12)
            np.testing.assert_allclose(f.T@f, np.eye(3), atol=1e-12)
        np.testing.assert_allclose(p.beam_frame(0,0,90)[:,0], [0,0,1], atol=1e-12)
        np.testing.assert_allclose(p.beam_frame(90,90,0)[:,2], [0,0,-1], atol=1e-12)

    def test_perspective_cube_matches_analytical_frustum(self):
        mask, corner, spacing = volume()
        frame = p.beam_frame(0,0,0)
        x,y = p.project_target((mask,corner,spacing), np.zeros(3), 10, frame)
        # Near surface is y=-1: projected half-width = 10/9 cm.
        expected = np.arange(-11,11)*.1 + .05
        np.testing.assert_allclose(np.unique(x), expected, atol=1e-12)
        np.testing.assert_allclose(np.unique(y), expected, atol=1e-12)
        self.assertEqual(len(x), 22*22)

    def test_hole_and_disconnected_components_are_not_filled(self):
        mask,corner,sp = volume()
        mask[:,:,8:12] = False
        x,y = p.project_target((mask,corner,sp),np.zeros(3),100,p.beam_frame(0,0,0))
        self.assertFalse(np.any(abs(x) < .15))
        self.assertTrue(np.any(x < -.5) and np.any(x > .5))

    def test_voxel_rays_match_brute_force_boxes(self):
        rng = np.random.default_rng(824)
        mask = rng.random((5,6,7)) < .15
        corner, sp = np.array([-.7,-.6,-.5]), .2
        occupied = np.argwhere(mask)[:,::-1]
        low = corner + occupied*sp
        high = low + sp
        for g,c,t in [(0,0,0),(73,41,29),(270,0,90)]:
            frame = p.beam_frame(g,c,t)
            source = frame[:,2]*100
            uv = rng.uniform(-1.4,1.4,(150,2))
            d = uv[:,:1]*frame[:,0] + uv[:,1:]*frame[:,1] - source
            actual = p.ray_hits(mask,corner,sp,source,d)
            expected=[]
            for direction in d:
                with np.errstate(divide='ignore',invalid='ignore'):
                    a,b=(low-source)/direction,(high-source)/direction
                enter=np.maximum(np.minimum(a,b).max(axis=1),0)
                leave=np.maximum(a,b).min(axis=1)
                expected.append(np.any(leave > enter))
            np.testing.assert_array_equal(actual,expected)

    def test_off_axis_rotation(self):
        mask,corner,sp=volume()
        corner=corner+np.array([3.,0.,2.])
        x,y=p.project_target((mask,corner,sp),np.zeros(3),100,p.beam_frame(0,0,0))
        self.assertGreater(x.mean(),2.9)
        self.assertGreater(y.mean(),1.9)
        x,y=p.project_target((mask,corner,sp),np.zeros(3),100,p.beam_frame(0,0,90))
        self.assertGreater(x.mean(),1.9)
        self.assertLess(y.mean(),-2.9)

    def test_voxelize_calls_native_full_grid_and_threshold(self):
        class Geometry:
            def HasContours(self): return True
            def GetBoundingBox(self): return [dict(x=0,y=0,z=0),dict(x=1,y=1,z=1)]
            def GetRoiGeometryAsVoxels(self, **kwargs):
                self.kwargs=kwargs
                values=np.zeros((12,12,12),dtype=np.uint8)
                values[5,5,5]=128
                values[5,5,6]=127
                return values.ravel()
        g=Geometry()
        mask,corner,sp=p.voxelize(g,.1)
        self.assertEqual(mask.sum(),1)
        self.assertEqual(g.kwargs['NrVoxels'],dict(x=12,y=12,z=12))
        np.testing.assert_allclose(corner,[-.1]*3)


class SurfaceTests(unittest.TestCase):
    def test_large_fine_grid_is_read_in_bounded_blocks_without_coarsening(self):
        class Geometry:
            def __init__(self): self.requests = []
            def HasContours(self): return True
            def GetBoundingBox(self): return [dict.fromkeys('xyz', 0.), dict.fromkeys('xyz', 13.)]
            def GetRoiGeometryAsVoxels(self, **kw):
                self.requests.append(kw)
                counts = [kw['NrVoxels'][k] for k in 'zyx']
                axes = [(np.arange(kw['NrVoxels'][k])+.5)*kw['VoxelSize'][k]+kw['Corner'][k] for k in 'zyx']
                z, y, x = [(a >= 0) & (a < 13.) for a in axes]
                return (255*(z[:, None, None] & y[None, :, None] & x[None, None, :])).astype(np.uint8).ravel()
        g = Geometry()
        surface = p.read_surface(g, .05)
        self.assertGreater(len(g.requests), 1)
        self.assertTrue(all(np.prod(list(r['NrVoxels'].values())) <= p.MAX_VOXELS for r in g.requests))
        self.assertTrue(all(r['VoxelSize'] == dict.fromkeys('xyz', .05) for r in g.requests))
        np.testing.assert_allclose(surface.faces.min(axis=(0, 1)), [0]*3, atol=1e-10)
        np.testing.assert_allclose(surface.faces.max(axis=(0, 1)), [13]*3, atol=1e-10)

    def test_slab_interfaces_do_not_create_caps_or_fill_holes(self):
        mask = np.zeros((18, 18, 18), dtype=np.uint8)
        mask[1:-1, 1:-1, 1:-1] = 255
        mask[:, 7:11, 7:11] = 0
        class Geometry:
            def HasContours(self): return True
            def GetBoundingBox(self): return [dict.fromkeys('xyz', -2.), dict.fromkeys('xyz', 2.)]
            def GetRoiGeometryAsVoxels(self, **kw):
                start = np.rint((np.array([kw['Corner'][k] for k in 'zyx'])+2.25)/.25).astype(int)
                count = [kw['NrVoxels'][k] for k in 'zyx']
                self.requests.append(kw)
                return mask[tuple(slice(a, a+b) for a, b in zip(start, count))].ravel()
        g = Geometry(); g.requests = []
        with patch.object(p, 'MAX_VOXELS', 1000):
            actual = p.read_surface(g, .25)
        reference = mask >= 128, np.full(3, -2.25), .25
        for angles in ((0, 0, 0), (90, 90, 0), (33, 28, 17)):
            frame = p.beam_frame(*angles)
            np.testing.assert_allclose(p.project_surface(actual, np.zeros(3), 100, frame),
                                       p.project_target(reference, np.zeros(3), 100, frame), atol=1e-10, rtol=0)
        self.assertTrue(all(np.prod(list(r['NrVoxels'].values())) <= 1000 for r in g.requests))

    def test_block_surface_matches_reference_rays_at_oblique_views(self):
        rng = np.random.default_rng(2928)
        mask = rng.random((9, 10, 11)) < .35
        v = mask, np.array([-1.1, -1., -.9]), .2
        surface = p.voxel_surface(v)
        for angle in (0, 13, 91, 153, 271):
            frame = p.beam_frame(angle, 21, 17)
            actual = p.project_surface(surface, np.array([.1, .3, -.4]), 100, frame)
            expected = p.project_target(v, np.array([.1, .3, -.4]), 100, frame)
            np.testing.assert_allclose(actual, expected, atol=1e-12, rtol=0)

    def test_surface_retains_hole_and_disconnected_components(self):
        mask, corner, spacing = volume()
        mask[:, :, 8:12] = False
        surface = p.voxel_surface((mask, corner, spacing))
        x, y = p.project_surface(surface, np.zeros(3), 100, p.beam_frame(0, 0, 0))
        self.assertFalse(np.any(abs(x) < .15))
        self.assertTrue(np.any(x < -.5) and np.any(x > .5))
        mask[:, :, 8:12] = True
        mask[8:12, :, 8:12] = False
        x, y = p.project_surface(p.voxel_surface((mask, corner, spacing)), np.zeros(3), 100, p.beam_frame(0, 0, 0))
        self.assertFalse(np.any((abs(x) < .15) & (abs(y) < .15)))

    def test_native_mesh_read_is_read_only_and_matches_perspective_cube(self):
        cube = p.voxel_surface(volume()).faces
        self.assertEqual(len(cube), 6)
        triangles = np.concatenate((cube[:, [0, 1, 2]], cube[:, [0, 2, 3]]))
        vertices = [dict(zip('xyz', point)) for point in triangles.reshape(-1, 3)]
        geometry = NS(HasContours=lambda: True, PrimaryShape=NS(
            Vertices=vertices, Indices=np.arange(len(vertices), dtype=np.int32), IsClosed=True))
        # No voxel conversion or mutation methods exist on this fixture.
        native = p.read_surface(geometry, .05)
        self.assertEqual(native.method, 'Native mesh')
        self.assertEqual(native.spacing, .05)
        for angles in ((0, 0, 0), (39, 18, 73), (270, 0, 0)):
            frame = p.beam_frame(*angles)
            reference = np.ones((40, 40, 40), dtype=bool), np.full(3, -1.), .05
            np.testing.assert_allclose(p.project_surface(native, np.zeros(3), 10, frame),
                                       p.project_target(reference, np.zeros(3), 10, frame), atol=1e-12, rtol=0)
        geometry.PrimaryShape.IsClosed = False
        with self.assertRaises(p.CalculationError):
            p.read_surface(geometry, .05)

    def test_surface_cancellation_and_invalid_geometry(self):
        def cancel(): raise p.Cancelled()
        with self.assertRaises(p.Cancelled): p.voxel_surface(volume(), cancel)
        surface = p.voxel_surface(volume())
        with self.assertRaises(p.Cancelled): p.project_surface(surface, np.zeros(3), 100, p.beam_frame(0, 0, 0), cancel)
        with self.assertRaises(p.CalculationError): p.Surface([], .1, 'invalid')
        with self.assertRaises(p.CalculationError): p.project_surface(surface, np.zeros(3), 1, p.beam_frame(0, 0, 0))

    def test_surface_arc_uses_new_projection_at_each_angle(self):
        b = beam(); b.DeliveryTechnique = 'DynamicArc'
        b.ArcRotationDirection = 'Clockwise'; b.ArcStopGantryAngle = 90
        b.Segments = [NS(RelativeWeight=w, CollimatorAngle=0, JawPositions=[-2, 2, -2, 2],
                         LeafPositions=[[-2]*4, [2]*4], IsVirtual=False,
                         DeltaGantryAngle=g, DeltaCouchAngle=0) for w, g in ((.25, 0), (.75, 90))]
        mask, corner, spacing = volume()
        surface = p.voxel_surface((mask, corner+np.array([4., 0., 0.]), spacing))
        bam, n = p.calculate_beam(b, NS(GetTreatmentMachine=lambda **kw: machine()), surface)
        self.assertAlmostEqual(bam, .25)
        self.assertEqual(n, 2)


class ApertureTests(unittest.TestCase):
    def test_cached_leaf_lookup_matches_strip_union_including_gaps_and_overlap(self):
        rng = np.random.default_rng(28129)
        x = np.r_[rng.uniform(-5, 5, 4000), [0, 0, 0, 0]]
        y = np.r_[rng.uniform(-3, 3, 4000), [-2, -.5, .5, 2]]
        layers = [(np.array([1.5, -1.5, 0.]), np.array([1., 1., 1.])),
                  (np.array([-.4, .4]), np.array([1., 1.]))]
        for axis in ('X', 'Y'):
            moving, across = (x, y) if axis == 'X' else (y, x)
            lookup = p.prepare_aperture(x, y, layers, axis)
            for _ in range(6):
                banks = []
                expected = np.ones(len(x), dtype=bool)
                for centres, widths in layers:
                    pair = np.sort(rng.uniform(-4, 4, (2, len(centres))), axis=0)
                    banks.extend(pair)
                    opened = np.zeros(len(x), dtype=bool)
                    for centre, width, left, right in zip(centres, widths, *pair):
                        opened |= ((across >= centre-width/2) & (across < centre+width/2)
                                   & (moving >= left) & (moving < right))
                    expected &= opened
                np.testing.assert_array_equal(p.aperture_open(x, y, layers, banks, axis, None, lookup), expected)

    def test_open_closed_half_and_jaws(self):
        x,y=p.project_target(volume(),np.zeros(3),100,p.beam_frame(0,0,0))
        layers=[(np.array([0.]),np.array([4.]))]
        for tips, expected in [([[-2],[2]],1), ([[0],[0]],0), ([[0],[2]],.5)]:
            self.assertAlmostEqual(p.aperture_open(x,y,layers,tips,'X',[-2,2,-2,2]).mean(),expected)
        self.assertAlmostEqual(p.aperture_open(x,y,layers,[[-2],[2]],'X',[-2,0,-2,2]).mean(),.5)

    def test_staggered_dual_layers_intersection_and_y_axis(self):
        x=np.array([-.75,-.75,.75,.75]); y=np.array([-.75,.75,-.75,.75])
        layers=[(np.array([0.]),np.array([2.])),(np.array([.5]),np.array([1.]))]
        banks=[[-2],[2],[0],[2]]
        np.testing.assert_array_equal(p.aperture_open(x,y,layers,banks,'X',None),[False,False,False,True])
        # Swap target coordinates AND motion axis; the aperture is unchanged.
        np.testing.assert_array_equal(p.aperture_open(y,x,layers,banks,'Y',None),[False,False,False,True])

    def test_native_jaw_order_asymmetric(self):
        layers=[(np.array([0.]),np.array([20.]))]
        actual=p.aperture_open(np.array([-1,1,1]),np.array([1,1,3]),layers,[[-5],[5]],'X',[0,2,0,2])
        np.testing.assert_array_equal(actual,[False,True,False])

    def test_invalid_banks_and_jaw_placeholders(self):
        layers=[(np.array([0.]),np.array([2.]))]
        with self.assertRaises(p.CalculationError):
            p.aperture_open(np.array([0]),np.array([0]),layers,[[1],[-1]],'X',None)
        with self.assertRaises(p.CalculationError):
            p.aperture_open(np.array([0]),np.array([0]),layers,[[-1],[1]],'X',[0,0,0,0])


class AdapterTests(unittest.TestCase):
    def test_phase_timings_include_api_delay_and_reuse_static_leaf_lookup(self):
        class SlowSegment(NS):
            def __getattribute__(self, name):
                if name == 'LeafPositions':
                    time.sleep(.005)
                return super().__getattribute__(name)
        b = beam()
        values = vars(b.Segments[0]).copy()
        values['RelativeWeight'] = .5
        b.Segments = [SlowSegment(**values), SlowSegment(**values)]
        timings = p.Timings()
        with patch.object(p, 'prepare_aperture', wraps=p.prepare_aperture) as lookup:
            result, _ = p.calculate_beam(b, NS(GetTreatmentMachine=lambda **kw: machine()),
                                         p.voxel_surface(volume()), timings=timings)
        self.assertEqual(result, 0)
        self.assertEqual(lookup.call_count, 1)
        self.assertEqual(timings.calls['Projection'], 1)
        self.assertGreaterEqual(timings.seconds['API'], .01)
        self.assertIn('MLC', timings.summary())

    def test_weighted_reference_example(self):
        self.assertAlmostEqual(p.weighted_mean([.10,.55,.30,.15],[10,30,40,20]),.325)
        self.assertAlmostEqual(p.weighted_mean([.1,.5],[100,300]),.4)

    def test_native_mu_and_multiple_segments(self):
        b=beam()
        b.Segments[0].RelativeWeight=.25
        b.Segments.append(NS(RelativeWeight=.75,CollimatorAngle=0,JawPositions=[-2,2,-2,2],LeafPositions=[[0]*4,[0]*4]))
        db=NS(GetTreatmentMachine=lambda **kw:machine())
        bam,n=p.calculate_beam(b,db,volume())
        self.assertAlmostEqual(bam,.75)
        self.assertEqual(n,2)

    def test_zero_weight_sample_is_not_dose(self):
        b=beam()
        b.Segments.append(NS(RelativeWeight=0,CollimatorAngle=0,JawPositions=[0]*4,LeafPositions=[]))
        bam,n=p.calculate_beam(b,NS(GetTreatmentMachine=lambda **kw:machine()),volume())
        self.assertEqual(bam,0)
        self.assertEqual(n,2)

    def test_arc_reprojects_target_at_each_control_point(self):
        b=beam(); b.DeliveryTechnique='DynamicArc'
        b.ArcRotationDirection='Clockwise'; b.ArcStopGantryAngle=90
        b.Segments=[NS(RelativeWeight=w,CollimatorAngle=0,JawPositions=[-2,2,-2,2],
                       LeafPositions=[[-2]*4,[2]*4],IsVirtual=False,
                       DeltaGantryAngle=g,DeltaCouchAngle=0) for w,g in [(.25,0),(.75,90)]]
        mask,corner,sp=volume()
        v=(mask,corner+np.array([4.,0.,0.]),sp)
        bam,n=p.calculate_beam(b,NS(GetTreatmentMachine=lambda **kw:machine()),v)
        # At gantry 0 the target is outside the aperture; at 90 it is central.
        self.assertAlmostEqual(bam,.25)
        self.assertEqual(n,2)

    def test_non_normalized_weights_rejected(self):
        b=beam(); b.Segments[0].RelativeWeight=.5
        with self.assertRaises(p.CalculationError): p.beam_samples(b)

    def test_machine_version_and_leaf_mismatch_rejected(self):
        b=beam(); b.MachineReference.CommissioningTime=2
        with self.assertRaises(p.CalculationError): p.model_for_beam(b,NS(GetTreatmentMachine=lambda **kw:machine()))
        b=beam(); b.UpperLayer.LeafWidths[0]=.9
        with self.assertRaises(p.CalculationError): p.model_for_beam(b,NS(GetTreatmentMachine=lambda **kw:machine()))

    def test_arc_direction_wraparound_and_virtual_guard(self):
        for direction,start,stop,deltas in [('Clockwise',350,10,[0,10,20]),('CounterClockwise',10,350,[0,-10,-20])]:
            b=beam(); b.ArcRotationDirection=direction; b.GantryAngle=start; b.ArcStopGantryAngle=stop
            b.DeliveryTechnique='DynamicArc'
            b.Segments=[NS(RelativeWeight=w,IsVirtual=False,DeltaGantryAngle=d,DeltaCouchAngle=0) for w,d in zip([.25,.5,.25],deltas)]
            _,_,angles,_=p.beam_samples(b)
            np.testing.assert_allclose(angles,np.array(deltas)+start)
            b.Segments[1].IsVirtual=True
            with self.assertRaises(p.CalculationError):p.beam_samples(b)

    def test_unsupported_patient_orientation(self):
        b=beam(); b.PatientPosition='FeetFirstSupine'
        with self.assertRaises(p.CalculationError):p.beam_samples(b)


class GuiTests(unittest.TestCase):
    def setUp(self):
        import tkinter as tk
        from tkinter import ttk
        self.root=tk.Tk(); self.root.withdraw()
        self.b=beam()
        self.bs=NS(DicomPlanLabel='Synthetic',Beams=[self.b],Modality='Photons',EnableDynamicTracking=False,
                   GetPlanningExamination=lambda:NS(Name='Synthetic_CT'))
        case=NS(PatientModel=NS(RegionsOfInterest=[NS(Name='Synthetic_PTV',Type='Ptv')]))
        self.viewer=p.Viewer(self.root,case,NS(BeamSets=[self.bs]),NS(GetTreatmentMachine=lambda **kw:machine()),tk,ttk)
        class Geometry:
            def HasContours(self):return True
            def GetBoundingBox(self):return [dict(x=-1,y=-1,z=-1),dict(x=1,y=1,z=1)]
            def GetRoiGeometryAsVoxels(self,**kw):
                size=[kw['NrVoxels'][k] for k in 'zyx']
                vals=np.zeros(size,dtype=np.uint8); vals[1:-1,1:-1,1:-1]=255
                return vals.ravel()
        case.PatientModel.StructureSets={'Synthetic_CT':NS(RoiGeometries={'Synthetic_PTV':Geometry()})}

    def tearDown(self):self.root.destroy()

    def test_full_gui_plan_calculation_and_failed_beam(self):
        self.viewer.calculate()
        self.assertIn('0.0000',self.viewer.pam.get())
        self.assertEqual(len(self.viewer.table.get_children()),1)
        self.bs.Beams.append(beam()); self.bs.Beams[-1].PatientPosition='HeadFirstProne'
        self.viewer.calculate()
        self.assertEqual(self.viewer.pam.get(),'Plan PAM: unavailable')
        self.assertIn('1 failed',self.viewer.status.get())
        self.assertFalse(self.viewer.running)
        self.assertEqual(str(self.viewer.calculate_button['state']),'normal')
        for phase in ('API', 'Surface', 'Projection', 'MLC'):
            self.assertIn(phase, self.viewer.detail.get())
        self.assertIsNone(self.viewer.timer_job)
        self.assertIn('Elapsed:', self.viewer.elapsed.get())
        for row in self.viewer.table.get_children():
            self.assertGreaterEqual(float(self.viewer.table.set(row, 'seconds')), 0)

    def test_no_positive_mu_beams_no_pam(self):
        self.b.BeamMU=0
        self.viewer.calculate()
        self.assertEqual(self.viewer.pam.get(),'Plan PAM: unavailable')

    def test_cancel_no_partial_plan_pam(self):
        self.root.after(0,self.viewer.cancel)
        self.viewer.calculate()
        self.assertEqual(self.viewer.pam.get(),'Plan PAM: unavailable')
        self.assertIn('Cancelled',self.viewer.status.get())
        self.assertIsNone(self.viewer.timer_job)
        self.assertFalse(self.viewer.running)

    def test_live_timer_and_stopped_time(self):
        self.viewer.running = True
        self.viewer.started = 100
        self.viewer.beam_started = 101
        row = self.viewer.table.insert('', 'end', values=('Synthetic', 'Synthetic', 100, 1, 0, 0, 'Running'))
        self.viewer.active_row = row
        with patch.object(p.time, 'perf_counter', return_value=102.5):
            self.viewer.tick()
        self.assertIn('Elapsed: 2.5 s', self.viewer.elapsed.get())
        self.assertIn('Current beam: 1.5 s', self.viewer.elapsed.get())
        self.assertEqual(self.viewer.table.set(row, 'seconds'), '1.5')
        self.root.after_cancel(self.viewer.timer_job)
        self.viewer.timer_job = None
        self.viewer.running = False
        self.viewer.stopped = 103
        with patch.object(p.time, 'perf_counter', return_value=900):
            self.viewer.tick()
        self.assertIn('Elapsed: 3.0 s', self.viewer.elapsed.get())

    def test_mixed_planning_examinations_no_plan_pam(self):
        bs2=NS(**vars(self.bs)); bs2.GetPlanningExamination=lambda:NS(Name='Synthetic_CT_2')
        self.viewer.plan.BeamSets.append(bs2)
        structures=self.viewer.case.PatientModel.StructureSets
        structures['Synthetic_CT_2']=structures['Synthetic_CT']
        self.viewer.calculate()
        self.assertEqual(self.viewer.pam.get(),'Plan PAM: unavailable')
        self.assertIn('different planning examinations',self.viewer.status.get())


if __name__=='__main__':
    unittest.main(verbosity=2)
