# PAM BAM implementation plan

Goal: a patient-free, read-only research library for single- and dual-layer aperture modulation, usable from ESAPI and RayStation scripts.

Architecture: independent Python and C# implementations of one explicit BEV geometry contract, verified against shared synthetic inputs. Native wrappers read weights but require site-validated aperture and target-projection providers. No TPS writes, DICOM export, patient lookup, GUI launcher, or machine geometry guesses.

Decisions: BAM is the MU-weighted blocked target fraction of one beam. PAM is the beam-MU-weighted mean of BAM, following Hernandez et al. (2025), not the ClearPlan native Beam.WeightFactor variant. Both layers are intersected in the same isocenter-plane coordinates; rectangle unions prevent double counting. No new metric or clinical threshold is claimed.

- [x] Write Python analytical tests and observe the missing-library failure.
- [x] Implement rectangle union/intersection, explicit leaf strips, single/dual functions, BAM and PAM; run `python -m unittest discover -s tests -v`.
- [x] Add the C# source counterpart and a dependency-free console test harness, then compare deterministic cross-language vectors with Python.
- [x] Provide read-only ESAPI and RayStation adapter examples with an explicit target-projection boundary; validate with synthetic objects and compile the ESAPI adapter against locally installed assemblies, if available.
- [x] Document units, interval semantics, installation, limitations, original paper, related RTplan Complexity Lens link, and license. No vendor DLLs or clinical data enter Git.
- [x] Inspect the exact Git index, create `Kiragroh/PAM-BAM`, push only this new repository and verify remote HEAD.

Separate file organization: retain active code paths; move historical root folders into a dated archive, list every move, check relative file names/sizes/hashes before and after. Keep references distinct. Clinical Downloads require confirmation before entering the synchronized Seafile project; do not infer de-identification from name masking.
