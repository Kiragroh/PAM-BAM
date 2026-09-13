# RayStation simple GUI 1.1

Download **[PAM_BAM.py](PAM_BAM.py)** and run it as a file script in RayStation
CPython with a plan open. Select the target ROI, select 2.0 / 1.0 / 0.5 mm and
press **Calculate**. No package installation or sibling files are needed if
NumPy and Tkinter are present in the RayStation Python environment.

- BAM for each treatment beam, PAM for the entire open plan using beam MU.
- Live total/beam timer and a **Seconds** column; the final times remain visible.
- Surface projection with a **2.0 mm default** grid, also selectable at 1.0/0.5 mm.
- Read-only: no patient save, contour change, dose calculation or export.
- Divergent target projection at every CP; actual MLC strip widths and centres;
  effective opening is the intersection of the layers and native jaw limits.
- HFS photon SMLC, DMLC, DynamicArc and StaticArc with supported IEC geometry.
- Static couch yaw; pitch/roll/gimbal, virtual CPs and unsupported accessories
  produce an explanation. Unknown geometry never becomes a valid zero result.
- Zero-MU and setup beams are excluded. A failed treatment beam or different
  planning examinations prevents calculation of a partial plan PAM.

## Numerical method

Native `RelativeWeight` is used as a segment share of total beam MU, as
documented in the RayStation v2025 SP2 API. It is not interpreted as cumulative
meterset. Like the ESAPI viewer, version 1.1 projects a surface instead of
tracing each BEV ray through the volume. The source is shown in each beam row:

- **Native mesh**: read `PrimaryShape.Vertices`, `Indices` and `IsClosed` when
  available. No resampling; the selected grid affects only the BEV plane.
- **Voxel surface**: when no native mesh is exposed, read the entire ROI at the
  selected spacing, threshold at >=128/255, and extract its exact exposed block
  faces in memory. Only coplanar rectangles are merged. Both voxel and BEV grids
  still use the selected spacing, as in 1.0. There is no smoothing or convex hull.

Both paths include beam divergence and retain holes/disconnected targets.
Pixel centres exactly on silhouette edges are included; such boundary samples
can differ from the old ray routine's strict entry/exit convention.
No `SetRepresentation`, temporary ROI or TPS model modification is used. Native
meshes and thresholded voxel surfaces can give different sampled results;
cross-TPS agreement is not implied by matching grid numbers. This also differs
from the library's exact area calculation for supplied rectangle unions.
Compare 1.0 and 0.5 mm results to assess convergence for a particular geometry.

Surface preparation is cached across beams on the same examination; consecutive
identical CP frames reuse the projection. Changing frames get a new projection.
All API access remains on the script thread. The timer updates during the
calculation through Tk event processing (about every 100 ms). A synchronous
native API read can temporarily block repainting; elapsed time includes that
interval and catches up afterward. The first beam's time includes surface setup.

Jawless machines need absent `JawPhysics` or valid native fixed field limits.
Zero placeholder limits are rejected and require local geometry verification.
The script contains the full scope and coordinate assumptions in its header.

## Test status

The user reports a successful RayStation run of **version 1.0**, with about
**4 seconds per field at a 2 mm grid**. That original source is preserved in
[commit dac51a6](https://github.com/Kiragroh/PAM-BAM/blob/dac51a6/scripts/raystation/PAM_BAM.py)
(SHA-256
`D8060D8E376046431D0328DB31258C9C1291FB05846D9C5C4ECAAEBE0C13B0A9`).
Version **1.1 has not yet been run inside RayStation**. Prior feedback for 1.0
does not validate the new surface adapter. Neither establishes dosimetric
accuracy, all-machine compatibility or clinical commissioning.

The **28 synthetic tests** cover independent ray/box intersections, analytical
perspective projection, holes/disconnected targets, MLC/jaw intersections,
MU weights, changing CP angles, API guards and Tk GUI calculation/cancellation,
native mesh reading, surface-vs-ray comparison at oblique angles, and the timer.
They contain only synthetic geometry and names:

```text
python scripts/raystation/tests/test_pam_bam.py
```

## Synthetic speed comparison

One local CPU run of [benchmark_projection.py](tests/benchmark_projection.py),
five views of a sphere with a through-hole in an 8 cm box. Every occupied BEV
sample matched the old ray algorithm. These are numerical-kernel timings;
RayStation API calls, MLC evaluation and GUI overhead are excluded.

| Grid | Surface setup, once | Old rays, 5 views | New surface, 5 views | Projection speedup |
|---|---:|---:|---:|---:|
| 2.0 mm | 0.024 s | 0.027 s | 0.012 s | 2.3× |
| 1.0 mm | 0.096 s | 0.177 s | 0.049 s | 3.6× |
| 0.5 mm | 0.400 s | 1.802 s | 0.182 s | 9.9× |

Including setup, this five-view example is 0.8×/1.2×/3.1× as fast at 2/1/0.5 mm.
Setup can outweigh the gain for a small number of views. Use the new live timer
to measure the complete calculation in RayStation; these timings do not promise
equal wall-clock performance to Eclipse.

```text
python scripts/raystation/tests/benchmark_projection.py
```

Tk tests require a desktop session. API documentation basis: v2025 SP2 / 17.2.0.
The screenshot and clinical plan/ROI names are deliberately not included.
