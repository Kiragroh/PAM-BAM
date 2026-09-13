# PAM and BAM

Small, read-only research functions for **plan aperture modulation (PAM)** and
**beam aperture modulation (BAM)**, with explicit **single-layer** and
**dual-layer MLC** geometry. Python for RayStation scripting and standalone work;
C# source for ESAPI projects. No cloud service, patient export, or TPS mutation.

This is an implementation of published metrics, **not a new metric** and not
clinically commissioned software.

## Simple GUI scripts

| TPS | Standalone script | Current test status |
|---|---|---|
| RayStation | [PAM_BAM.py](scripts/raystation/PAM_BAM.py) | GUI 1.2: bounded voxel reads, cached MLC lookup and phase timings; 32 synthetic tests; new native test pending |
| Eclipse / ESAPI, SD / HD / Halcyon SX | [PAM_BAM.cs](scripts/esapi/PAM_BAM.cs) | GUI 1.2: ESAPI 18 compile; 65 synthetic assertions; fixed user-confirmed 28+29 Halcyon order |

The ESAPI GUI offers Auto, TrueBeam SD (Millennium 120), TrueBeam HD (HD120) and
Halcyon SX dual-layer profiles, using the native target mesh and cumulative CP
weights. The default BEV grid is 2 mm. Halcyon uses the user-confirmed native
order: 28 pairs followed by 29. It reports BAM, plan PAM and elapsed
time per field. See [ESAPI usage and scope](scripts/esapi/README.md).

[RayStation: PAM_BAM.py](scripts/raystation/PAM_BAM.py) is a standalone, read-only
file script: choose the target ROI and grid spacing, then calculate BAM per beam
and MU-weighted PAM for the open plan. It includes a perspective ROI projector
at every native control point. NumPy and Tkinter are required; the repository
package does not have to be installed.

**Runtime feedback:** the user successfully ran GUI version 1.0 in RayStation
and reported approximately **4 seconds per field at a 2 mm grid**. This is one
reported application run, not a universal benchmark or commissioning result.
The tested 1.0 script is preserved in Git history. Version 1.1 added surface
projection and a live timer, but user feedback identified a 0.5 mm grid-size
failure and a remaining speed gap versus Eclipse. Version 1.2 fixes the total
voxel limit using bounded reads and adds phase timings. No clinical screenshot, plan
identifiers or target names are published. See [usage and validation](scripts/raystation/README.md).

## Definitions

At sample j, let T be the target's projected silhouette and O the effective
opening, in the same beam-limiting-device isocenter plane:

`AM_j = 1 - area(T_j intersect O_j) / area(T_j)`

`BAM_b = sum(AM_j * MU_j) / sum(MU_j)` for one beam.

`PAM = sum(BAM_b * MU_b) / sum(MU_b)` for all treatment beams and the **same target**.

Both range from 0 (target projection open throughout) to 1 (fully blocked).
They do not measure dose coverage, delivery accuracy, or MLC transmission.
No clinical accept/reject threshold is supplied. Unlike ClearPlan's separately
declared native PlanCheck-compatible variant, this library does **not** use
`Beam.WeightFactor` for plan aggregation.

For control-point data, `bam_from_control_points` / `BamFromControlPoints`
uses half of each interval's incremental MU at either endpoint. This is a
declared numerical approximation, not an exact integration of leaf motion.
The final cumulative weight is normalized explicitly and need not be 1.
Zero-MU intervals add no artificial weight. Missing geometry raises an error;
an explicitly empty opening means closed. A missing target never becomes PAM=0.

## Single and dual layers

`single_layer_aperture(layer, ...)` builds physical leaf strips.
`dual_layer_aperture(proximal, distal, ...)` intersects the two strip unions.
The C# equivalents are `Metrics.SingleLayerAperture` and `DualLayerAperture`.

Supply each layer's actual leaf boundaries, bank positions and X/Y travel axis.
Staggered layers can have different leaf counts and different boundaries.
For a jawless Halcyon, omit physical jaws. Fixed field limits are a separate
argument. **Do not use the union of both openings, average their areas, assume
equal leaf widths, or apply an extra projection scale to values already at iso.**

All core lengths are **mm at isocenter**, and areas are **mm²**. Targets are unions
of axis-aligned rectangles, typically non-overlapping strips from a validated
BEV silhouette raster. Overlaps are unioned rather than double-counted. Holes and
disconnected components can be represented by the occupied strips. The area
calculation is exact for the supplied rectangles; it is only as accurate as the
target projection/raster supplied by the caller. A target bounding box is NOT a
valid substitute for its silhouette.

## Python quick start

Python 3.8 or later, standard library only at runtime. Put this repository on
the script's Python path, or install locally with `python -m pip install .`.
This project is not published to PyPI; do not assume `pip install pambam` fetches it.

```python
from pambam import Layer, dual_layer_aperture, bam_from_control_points, calculate_pam

# Entirely synthetic geometry, not a commissioned Halcyon profile.
upper = Layer((-10, 0, 10), (-10, -10), (10, 10))
lower = Layer((-15, -5, 5, 15), (0, -5, 0), (5, 5, 5))
opening = dual_layer_aperture(upper, lower)
target = [(-10, -10, 10, 10)]
beam = bam_from_control_points([target, target], [opening, opening],
                               [0, 1], 100, 'synthetic-target')
assert abs(beam.bam - 0.625) < 1e-12
assert abs(calculate_pam([beam]) - 0.625) < 1e-12
```

Run `python -m examples.synthetic` for a two-beam example.

## TPS adapters and what still needs local integration

The following are **library adapters**, distinct from the standalone GUIs above. Neither adapter
creates a patient context, writes a structure, launches a process, saves a plan,
nor calculates the 3D target silhouette automatically. A site-validated target
projection provider at every control point is required. This boundary is
deliberate: a beam-start outline must not be reused throughout a rotating arc.

| Adapter | Native input and explicit requirements |
|---|---|
| [ESAPI C#](csharp/EsapiAdapter.cs) | Reads `ControlPoint.LeafPositions`, `JawPositions`, `MetersetWeight` and beam meterset in MU on the owning STA. An exact MLC profile maps all native leaf indices once to one or two layers. Caller supplies a validated target-BEV provider. |
| [RayStation Python](pambam/raystation.py) | Reads `BeamMU`, `Segments`, `LeafPositions` and, only when configured, `JawPositions`. Converts native cm to mm; dual banks 0/1 and 2/3. Caller supplies leaf boundaries, target provider and an explicit validated cumulative-CP or sample-MU sequence. |

For ESAPI, include `csharp/PamBam.cs` and `csharp/EsapiAdapter.cs` in your existing
read-only ESAPI project and reference your installed vendor assemblies. The core
uses framework-compatible C# syntax; the **test runner** uses .NET 10. The adapter
was compiled against ESAPI 16.1 locally, but live TPS execution is not validated.

For RayStation, see [integration notes](docs/tps-integration.md). The adapter's
bank/length contract was checked against the v2025 SP2 scripting documentation
and tested with synthetic objects. The library adapter itself was **not run in RayStation**; the standalone GUI
has the user-reported run described above. API availability,
projection frames, machine profiles and segment/CP weights require local validation.

## Tests

```text
python -m unittest discover -s tests -v
dotnet run --project csharp/tests/PamBam.Tests.csproj
```

Tests cover analytical apertures, staggered layers, cm/mm conversion, target unions,
endpoint integration, zero-MU intervals, original plan weighting, and rejected
incomplete/nonfinite inputs. No patient data or proprietary vendor DLLs are included.
See [verification](docs/verification.md) for the tested scope and limitations.

## References and related tools

**Metric source:** Hernandez V, Lara-Aristimuño I, Abella R, Saez J.
*Quantification of the aperture modulation in radiotherapy treatment plans.*
Medical Physics. 2025;52(12):e70144.
[doi:10.1002/mp.70144](https://doi.org/10.1002/mp.70144).
Equations 1–2 define AM and PAM; Section 2.2 also describes beam-level BAM.
Please cite the original paper when using these metrics. The authors have not
endorsed or validated this implementation.

**Related independent resource:** [RTplan Complexity Lens metrics](https://rt-complexity-lens.lovable.app/metrics)
and its [source repository](https://github.com/matteomaspero/rt-complexity-lens).
Linked for comparison/discovery, not used as a numerical oracle or runtime service.
No code, screenshots, or data from that project are bundled here; numerical
equivalence has not been established. Never upload clinical data without approval.

## License and safety

MIT for this repository's original source; see [LICENSE](LICENSE).
The publication and linked projects retain their own licenses.
No treatment decisions, QA replacements, auto-approval or clinical safety claims.
Review target selection, all treatment beams, units, layer assignment, projection
convergence and sampling against independent local reference cases before use.
