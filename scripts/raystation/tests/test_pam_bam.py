"""Synthetic verification only; no RayStation session or patient data."""
import importlib.util
from pathlib import Path
from types import SimpleNamespace as NS
import unittest
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


class ApertureTests(unittest.TestCase):
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

    def test_no_positive_mu_beams_no_pam(self):
        self.b.BeamMU=0
        self.viewer.calculate()
        self.assertEqual(self.viewer.pam.get(),'Plan PAM: unavailable')

    def test_cancel_no_partial_plan_pam(self):
        self.root.after(0,self.viewer.cancel)
        self.viewer.calculate()
        self.assertEqual(self.viewer.pam.get(),'Plan PAM: unavailable')
        self.assertIn('Cancelled',self.viewer.status.get())

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
