# Verification scope

## Standalone GUI additions

- RayStation GUI 1.0: **22 synthetic tests passed**, including the Tk UI.
  The user subsequently reported successful native execution and approximately
  **4 seconds per field with a 2 mm grid**. The exact tested source is included.
  See [RayStation details](../scripts/raystation/README.md) for its source hash
  and the distinction between user runtime feedback and commissioning.
- ESAPI GUI 1.0 (TrueBeam / HD120): native x64/.NET Framework 4.8 compile against
  ESAPI API file version **18.0.1.261** passed; **41 synthetic numerical
  assertions passed**. Execution inside Eclipse is pending.
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
