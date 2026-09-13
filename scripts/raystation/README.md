# RayStation simple GUI

Download **[PAM_BAM.py](PAM_BAM.py)** and run it as a file script in RayStation
CPython with a plan open. Select the target ROI, select 2.0 / 1.0 / 0.5 mm and
press **Calculate**. No package installation or sibling files are needed if
NumPy and Tkinter are present in the RayStation Python environment.

- BAM for each treatment beam, PAM for the entire open plan using beam MU.
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
meterset. Target voxels use >=128/255 occupancy; their perspective projection
is sampled on a BEV grid. Both grids use the selected spacing. This differs
from the library's exact area calculation for supplied rectangle unions.
Compare 1.0 and 0.5 mm results to assess convergence for a particular geometry.

Jawless machines need absent `JawPhysics` or valid native fixed field limits.
Zero placeholder limits are rejected and require local geometry verification.
The script contains the full scope and coordinate assumptions in its header.

## Test status

The user reports a successful RayStation run of **version 1.0**, with about
**4 seconds per field at a 2 mm grid**. The included source is byte-identical to
the tested file (SHA-256
`D8060D8E376046431D0328DB31258C9C1291FB05846D9C5C4ECAAEBE0C13B0A9`).
This establishes runtime feedback for that case, not validated dosimetric
accuracy, all-machine compatibility or clinical commissioning.

The **22 synthetic tests** cover independent ray/box intersections, analytical
perspective projection, holes/disconnected targets, MLC/jaw intersections,
MU weights, changing CP angles, API guards and Tk GUI calculation/cancellation.
They contain only synthetic geometry and names:

```text
python scripts/raystation/tests/test_pam_bam.py
```

Tk tests require a desktop session. API documentation basis: v2025 SP2 / 17.2.0.
The screenshot and clinical plan/ROI names are deliberately not included.
