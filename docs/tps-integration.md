# TPS integration contract

## Geometry shared by both platforms

Provide the target silhouette at every sampled gantry/collimator/couch geometry,
already projected into the beam-limiting-device plane at isocenter. Core rectangles
use `(x1, y1, x2, y2)` in mm, with ascending bounds. Use a validated patient-to-beam
transform including source divergence and the actual isocenter. Do not interpret
a structure being found in a DVH as evidence that its projected geometry is known.
Check projection orientation with asymmetric synthetic targets and aperture motion.
Quantify raster convergence; do not fill holes or connect separated targets.

Geometry must be copied before passing detached values to a worker. Neither
adapter starts threads. An STA check alone does not prove it is the ESAPI owning
thread: the host must call the adapter from that same thread and active context.

## ESAPI

Call `EsapiAdapter.CalculateBeam(beam, profile, projectionProvider, targetKey)`.
The supplied beam must be a treatment photon MLC beam without blocking accessories.
`EsapiProfile.ExactMlcModel` must match the native model exactly. Each layer's
`SourceLeafIndices` references columns of the two-bank native LeafPositions array;
every native pair must appear exactly once across the layers. Supply physical
`BoundariesMm`, not equal-width guesses. Explicit `JawMode` is `None`, `Physical`
or `Fixed`. Choose `None` for a jawless device unless a validated fixed boundary
must also clip the opening. A fixed boundary is never labelled as physical jaws.

The callback signature is `Func<Beam, ControlPoint, IEnumerable<PamBam.Rect>>`.
It must produce the target silhouette at that control point. A vendor
`GetStructureOutlines` beam-start outline alone is insufficient for a whole arc.
This release does not supply that site's projector or commissioned machine catalog.

Pass all complete treatment beam results for a single target to
`Metrics.CalculatePam(results)`. Exclude setup/imaging beams explicitly. Do not
silently drop a failed beam or average independently selected targets.

## RayStation

Use `pambam.raystation.calculate_beam` inside the native RayStation scripting
session, with an already selected photon beam. Importing this module does not
import `connect`, query the database or modify the current context.

```python
from pambam.raystation import calculate_beam

def review_selected_beam(beam, commissioned_leaf_boundaries_mm,
                         projected_target_at_segment, validated_cumulative_weights):
    return calculate_beam(
        beam, commissioned_leaf_boundaries_mm, projected_target_at_segment,
        'same-target-for-all-beams', jaw_mode='none',
        cumulative_weights=validated_cumulative_weights)
```

`projected_target_at_segment(segment,index)` returns the rectangle union in mm.
RayStation native leaf banks are cm and are converted once by the adapter. The
boundary arrays passed to it are already mm. For v2025 SP2 dual-layer input,
banks 0/1 are upper left/right and 2/3 lower left/right. Physical jaw input order
is X1, X2, Y1, Y2; the adapter rearranges this into rectangle order.

Weight semantics are deliberately explicit. For verified VMAT endpoint data,
provide the full cumulative sequence (one value for each native sample). For
static segments or a site-computed endpoint quadrature, pass `sample_weights_mu`
instead; their sum must match native `BeamMU`. Do not pass a relative segment
weight as a cumulative weight simply because values happen to be between 0 and 1.
Zero-MU control points are allowed; zero-MU beams are not scored. The sample/CP
weight extraction remains a site-specific integration step in this first release.

The adapter is not suitable for ion, CyberKnife, tomotherapy, applicator or block
geometries. Validate the photon treatment technique and accessories in the host.
The native model's advertised properties do not replace a live acceptance test.

## Reproducibility checklist

Record software version, selected target, included beam count, source units,
machine-profile version/hash, BEV frame convention, projection raster/convergence,
and MU sampling convention with each local analysis. Keep identifiers and geometry
on approved clinical storage; only use synthetic fixtures in public issue reports.
