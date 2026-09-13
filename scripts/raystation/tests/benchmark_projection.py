"""Synthetic projection benchmark. No RayStation API, patient data or I/O export.

Times surface preparation once, then five views of an 8 cm box containing a
sphere with a through-hole. Compares every occupied BEV sample with the old
voxel-ray algorithm. Timings are local CPU timings, not TPS runtime claims.
"""
import importlib.util
from pathlib import Path
import time
import numpy as np

spec = importlib.util.spec_from_file_location("pam_bam", Path(__file__).resolve().parents[1] / "PAM_BAM.py")
p = importlib.util.module_from_spec(spec)
spec.loader.exec_module(p)


def run():
    for mm in (2.0, 1.0, .5):
        spacing = mm/10
        n = round(8/spacing)
        z, y, x = np.ogrid[:n, :n, :n]
        x, y, z = [(a+.5)*spacing-4 for a in (x, y, z)]
        mask = (x*x+y*y+z*z < 3.5**2) & (x*x+z*z > .7**2)
        volume = mask, np.full(3, -4.), spacing
        start = time.perf_counter()
        surface = p.voxel_surface(volume)
        prep = time.perf_counter()-start
        old_seconds = new_seconds = 0.0
        for gantry, couch, collimator in ((0, 0, 0), (37, 13, 27), (90, 0, 90), (183, 21, 17), (273, 0, 3)):
            frame = p.beam_frame(gantry, couch, collimator)
            start = time.perf_counter()
            old = p.project_target(volume, np.zeros(3), 100, frame)
            old_seconds += time.perf_counter()-start
            start = time.perf_counter()
            new = p.project_surface(surface, np.zeros(3), 100, frame)
            new_seconds += time.perf_counter()-start
            np.testing.assert_allclose(np.column_stack(new), np.column_stack(old), atol=1e-10, rtol=0)
        print(f"{mm:.1f} mm | {len(surface.faces)} faces | prep {prep:.3f} s | "
              f"5 views: rays {old_seconds:.3f} s, surface {new_seconds:.3f} s | "
              f"projection {old_seconds/new_seconds:.1f}x | "
              f"including prep {old_seconds/(new_seconds+prep):.1f}x | identical BEV samples", flush=True)


if __name__ == "__main__":
    run()
