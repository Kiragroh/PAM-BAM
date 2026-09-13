"""Warm, alternating v1.2/v1.3 CPU benchmark; synthetic data only.

No RayStation API, GUI, MLC evaluation or projection-cache hits are timed.
Surface preparation is excluded. Every occupied BEV sample is compared.
"""
import platform
import time
import numpy as np
from test_pam_bam import p
from reference_projector_v12 import project_surface as reference


def run():
    print('Python {} | NumPy {} | {}'.format(platform.python_version(), np.__version__, platform.machine()))
    frames = [p.beam_frame(i*7.2, 0, 27) for i in range(50)]
    for mm in (2., 1., .5):
        step = mm/10
        n = round(8/step)
        z, y, x = np.ogrid[:n, :n, :n]
        x, y, z = [(a+.5)*step-4 for a in (x, y, z)]
        mask = (x*x+y*y+z*z < 3.5**2) & (x*x+z*z > .7**2)
        surface = p.voxel_surface((mask, np.full(3, -4.), step))
        elapsed = {'v1.2': [], 'v1.3': []}
        for repeat in range(4):
            results = {}
            methods = [('v1.2', reference), ('v1.3', p.project_surface)]
            if repeat % 2:
                methods.reverse()
            for name, projector in methods:
                start = time.perf_counter()
                results[name] = [projector(surface, np.zeros(3), 100, frame) for frame in frames]
                seconds = time.perf_counter()-start
                if repeat:
                    elapsed[name].append(seconds)
            for actual, expected in zip(results['v1.3'], results['v1.2']):
                np.testing.assert_array_equal(actual, expected)
        old, new = [float(np.median(elapsed[k])) for k in ('v1.2', 'v1.3')]
        print('{:.1f} mm | {} faces | 50 views: v1.2 {:.3f} s, v1.3 {:.3f} s | {:.2f}x | exact masks'.format(
            mm, len(surface.faces), old, new, old/new), flush=True)


if __name__ == '__main__':
    run()
