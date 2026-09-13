# ESAPI simple GUI — TrueBeam / HD120

Download **[PAM_BAM.cs](PAM_BAM.cs)** and run the single-file plug-in from Eclipse
with one external photon plan open. Select the target ROI and press
**Calculate**. Default BEV grid: **2.0 mm**, also selectable at 1.0 and 0.5 mm.
No other repository source, NuGet package or external process is needed.

The GUI shows beam ID, native MLC model, MU, CP count, BAM, elapsed seconds per
field and status, plus the MU-weighted PAM of the complete plan. It provides
progress, cancellation and invalidation when the ROI or grid selection changes.
Setup fields and zero-MU fields are excluded. A failed treatment field prevents
a partial plan PAM; completed field rows remain visible.

## Geometry and scope

- **TrueBeam / HD120, HFS**, static couch yaw, conventional IEC frame. The script
  checks native source positions against that frame, including two independent
  gantry directions. Unknown source geometry is rejected.
- HD120 at isocenter: **60 pairs**, central **32 pairs × 2.5 mm**, outer **14
  pairs × 5 mm on each side**; leaf strip boundaries run from −110 to +110 mm.
  This profile is explicit because ESAPI's public MLC object does not expose
  leaf widths. All native positions and mesh coordinates are in millimetres.
- Exact semantic MLC aliases: `HD120`, `Varian HD120`, `High Definition 120`,
  `Varian High Definition 120` (spaces/punctuation/case do not matter).
  Unknown model names show a status message. Do not add an alias without
  checking commissioned geometry and native leaf ordering.
- Millennium 120 and Halcyon are **not** mapped by leaf count. This GUI targets
  the requested HD120. The library's separate dual-layer adapter remains
  available for sites with validated profiles and index mappings.
- Applicators, custom blocks, wedges, compensators, motion compensation,
  changing table positions and unsupported orientation/geometry are rejected.
- The native ROI surface mesh is copied in memory and projected with beam
  divergence **at every CP**. Projected triangles form a union on the BEV grid;
  holes and disconnected components are retained. No convex hull, bounding-box
  target substitute, start-angle outline reuse, temporary ROI or DICOM export.
- Grid spacing affects BEV sampling only; the native mesh is not revoxelized.
  This differs from the RayStation GUI's voxel-plus-BEV approximation, so
  bit-identical cross-TPS results are not implied. Compare resolutions locally.

## MU weighting

ESAPI `ControlPoint.MetersetWeight` is **cumulative**, unlike the native segment
weights used in the RayStation GUI. Each interval's incremental MU contributes
half to either endpoint. We normalize by the final cumulative weight; repeated
weights add no MU. This is the repository's trapezoidal CP approximation,
not continuous integration of moving leaves. Plan PAM uses native beam MU;
`Beam.WeightFactor` is not multiplied in.

AM/BAM/PAM describe geometry, not dose coverage or clinical pass/fail. The
plug-in is read-only: no writeable attribute, BeginModifications, patient save,
contour changes, dose calculation, network traffic, patient logging or exports.
All ESAPI reads remain on the owning STA thread; UI event pumping does not
transfer API objects to workers.

## Verification

- Compiled as an x64 library using C# 5 syntax and .NET Framework 4.8 references
  against locally available **ESAPI API/Types, API file version 18.0.1.261**.
  This compile check includes the actual ESAPI adapter and WPF GUI.
- **41 vendor-free synthetic assertions passed**: perspective mesh projections
  vs independent ray/box intersections, disconnected targets and a through-hole,
  source-plane rejection, cancellation, IEC cardinal/couch/collimator frames,
  actual HD120 inner/outer widths, jaw clipping, closed/half-open apertures,
  model/bank guards, zero/decreasing cumulative MU and plan aggregation.
- **Not yet executed inside Eclipse.** Compilation and numerical tests do not
  validate native runtime behaviour, every machine configuration or clinical
  accuracy. Test an asymmetric ROI/field at several gantry, collimator and couch
  angles, check MU and native MLC identity/index ordering, and compare 2/1/0.5 mm
  results before relying on the reported numbers.

Run the numerical tests without vendor assemblies:

```text
dotnet run --project scripts/esapi/tests/CoreTests.csproj
```

The test project uses .NET 10; this is **not** the Eclipse runtime requirement.
Alternatively, a .NET Framework C# compiler can build the same tests:

```text
csc /target:exe /platform:x64 /langversion:5 /define:PAMBAM_CORE_TEST /out:CoreTests.exe scripts/esapi/PAM_BAM.cs scripts/esapi/tests/CoreTests.cs
CoreTests.exe
```

For a local native compile check, build `PAM_BAM.cs` as an x64 class library with
the installed `VMS.TPS.Common.Model.API.dll` and `VMS.TPS.Common.Model.Types.dll`,
plus .NET Framework 4.8 System, System.Core, System.Xml, System.Xaml, WindowsBase,
PresentationCore and PresentationFramework references. Vendor DLLs and compiled
binaries are not published here. A binary plug-in, if built locally, must have
the `.esapi.dll` extension; the `.cs` file can be run directly from Eclipse.

## API references

- [HD120 leaf geometry: Boudet et al., JACMP 2022, Section 2.1](https://pmc.ncbi.nlm.nih.gov/articles/PMC9278673/)
- [ControlPoint: coordinates, leaf banks, cumulative meterset](https://docs.developer.varian.com/api/16.1/VMS.TPS.Common.Model.API.ControlPoint.html)
- [MLC: model identity](https://docs.developer.varian.com/api/16.1/VMS.TPS.Common.Model.API.MLC.html)
- [Structure: native MeshGeometry](https://docs.developer.varian.com/api/16.1/VMS.TPS.Common.Model.API.Structure.html)
- [Beam: GetSourceLocation and GetStructureOutlines](https://docs.developer.varian.com/api/16.1/VMS.TPS.Common.Model.API.Beam.html)
- Metric: Hernandez et al., *Medical Physics* 2025;52(12):e70144,
  [doi:10.1002/mp.70144](https://doi.org/10.1002/mp.70144).
