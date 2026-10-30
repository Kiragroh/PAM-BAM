# Figure and screenshot sources

## PAM overview illustration

`pam-explained.png` and `pam-explained.svg` are the author's original PAM
explanation figure, reused with the author's permission for this repository.
The PNG is the 2400 × 2400 export; the SVG is its vector counterpart. Both
were copied without modifying their bytes from the original figure package.

The figure was visually matched to the image in
[Maximilian Grohmann's LinkedIn post](https://www.linkedin.com/posts/maximilian-grohmann-b70588b1_medicalphysics-radiotherapy-opensource-activity-7504920092158259200-Lwjo).
The repository bundles the original exports instead of depending on LinkedIn's
resized image or expiring media URL. The introductory explanation is adapted
from that post, with explicit MU-weighted interpretation.

The illustration uses synthetic geometry, not patient data. It shows one arc,
a single-layer MLC and jaws; the four blocked fractions and MU weights yield
PAM = BAM ≈ 0.325. See the [worked example](../method.md#worked-example-one-arc-four-control-points)
for the calculation and endpoint approximation. The original metric is from
[Hernandez et al. (2025)](https://doi.org/10.1002/mp.70144); the figure is the
repository author's illustration, not a figure reproduced from that paper.

SHA-256 of the bundled original exports:

- `pam-explained.png`: `565bc82958c296ad8b0ba0da301061b665b2c4f7a87ce9db106fdbecebf00cb9`
- `pam-explained.svg`: `fd65491b1bba4dc0e682987ebd31f18f7a4cf0328e0441c5d741e22d4f828152`

## TPS screenshots

- `raystation-pam-bam-1.3.png`: author-provided RayStation 2024 / GUI 1.3 example.
- `eclipse-esapi-pam-bam-1.2.png`: author-provided Eclipse / ESAPI GUI 1.2 example.

Both existing screenshots are retained unchanged. Their result values and
timings describe those example runs; see the [verification record](../verification.md).
