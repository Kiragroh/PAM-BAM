# Verification scope

## Standalone GUI additions

- RayStation GUI 1.0: **22 synthetic tests passed**, including the Tk UI.
  The user subsequently reported successful native execution and approximately
  **4 seconds per field with a 2 mm grid**. The tested 1.0 source is preserved
  in Git history. GUI 1.1 adds surface projection and a live total/beam timer;
  **28 synthetic tests passed**, including oblique comparisons against the old
  ray traversal and native mesh reading. A five-view 0.5 mm synthetic benchmark
  gave 9.9x faster projection (3.1x including surface setup), with identical BEV
  samples. The subsequent user run confirmed the voxel-surface path but found
  a remaining speed gap versus Eclipse and `ROI grid too large` at 0.5 mm.
  GUI 1.2 changes the 16M guard to a per-read bound with halo slabs, caches MLC
  sample-to-strip lookup, and displays per-beam API/surface/projection/MLC times.
  **32 synthetic tests passed**, including an above-16M 0.5 mm ROI, unchanged
  sampling, holes at slab interfaces and phase timing with simulated API delay.
  New execution of 1.2 inside RayStation is pending; no speed parity is claimed.
  See [RayStation details](../scripts/raystation/README.md) for its source hash
  and the distinction between user runtime feedback and commissioning.
- ESAPI GUI 1.2 (Auto / TrueBeam SD / HD / Halcyon SX): native x64/.NET Framework
  4.8 compile against ESAPI API file version **18.0.1.261** passed;
  **65 synthetic numerical assertions passed**. These include SD/HD widths,
  dual-layer intersection and stagger, fixed 28+29 native blocks with asymmetric
  regression inputs independent of the implementation's mapping,
  fixed Halcyon limits and rejection of missing/duplicate/conflicting mappings.
  The user confirms 28+29 in Eclipse and rejects the other choices from 1.1.
  Those unverified assumptions and the order selector have been removed.
  The earlier 71 synthetic assertions did not establish the native leaf order.
  Execution of the new 1.2 version inside Eclipse is pending.
  See [ESAPI details](../scripts/esapi/README.md).
- No patient screenshot, identifiers, target names, DICOM or vendor assemblies
  are included. The earlier library verification below describes the original
  release, not the later reported RayStation GUI run.

## Original library release

Date: 2026-09-13. Version: 0.1.0.

Executed locally with Python 3.12.14 (bundled runtime) and .NET SDK 10.0.204:

- `python -B -m unittest discover -s tests -v`: **20 tests passed**.
- `dotnet run --project csharp/tests/PamBam.Tests.csproj`: **42 analytical assertions passed**.
- Five shared `tests/geometry_cases.json` cases were independently evaluated by
  both implementations: jaw clipping, staggered dual layers, closed leaves,
  overlapping target projections, and Y-axis leaf travel.
- `python -B -m examples.synthetic`: BAM 0.0, BAM 0.625, combined PAM 0.46875.
- Both C# source files compiled as an x64 library with the locally installed
  ESAPI API/Types assemblies, API file version 16.1.3.17. Vendor binaries were
  used only as compiler references and are not distributed.

An independent code review also checked 1,000 synthetic integer-grid cases
against a separate geometry oracle. It found two defects: reuse of a target
callback's mutable buffer and floating-point endpoint overshoot. Both were
reproduced by failing regression tests and then fixed. The reviewer reran the
20 Python tests and 42 C# assertions and reported no remaining important findings
within that review scope.

Limits: no live ESAPI patient or RayStation session was opened for this release.
RayStation objects in tests are synthetic. The project does not implement a
universal native target projector, commissioned machine profiles, or universal
RayStation relative-weight extraction. Those are explicit integration requirements,
not silently substituted values. Tests verify software behavior, not clinical
accuracy, delivery feasibility or safe deployment on every machine/version.

The metric definition and attribution were checked against Hernandez et al.,
Medical Physics 2025;52(12):e70144, DOI 10.1002/mp.70144 (Section 2.2).
The requested RTplan Complexity Lens metrics page is linked as a related tool;
its client-rendered body was not available to the text browser, and no numerical
equivalence or endorsement is claimed.
