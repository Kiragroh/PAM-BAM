# PAM/BAM: definitions and geometry

[Back to the overview](../README.md)

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

## Worked example: one arc, four control points

The [overview figure](images/pam-explained.png) uses synthetic geometry: one
target, one arc, a generic single-layer MLC and jaws. The displayed blocked
fractions are rounded. This is a teaching example, not a patient measurement.

| Control point | Blocked target projection (AM) | Numerical MU weight | Contribution to BAM |
|---|---:|---:|---:|
| CP 1 | 0.10 | 10% | 0.010 |
| CP 2 | 0.55 | 30% | 0.165 |
| CP 3 | 0.30 | 40% | 0.120 |
| CP 4 | 0.15 | 20% | 0.030 |
| **Total** | | **100%** | **0.325** |

The cumulative MU sequence is `0, 20, 60, 100`. The interval increments are
`20, 40, 40` MU. Assigning half of each interval to either endpoint gives
`10, 30, 40, 20` MU: these are the numerical weights in the figure, not MU
delivered instantaneously at a control point. The weighted result is
`0.10*0.10 + 0.55*0.30 + 0.30*0.40 + 0.15*0.20 = 0.325`.
For this one-beam plan, **PAM = BAM = 0.325**: on average, 32.5% of the target
projection is blocked. For multiple beams, combine their BAM values using
each beam's total MU and the same target structure.

The small source/target drawing is a projection schematic. The four BEV panels
use the same synthetic 3D target at different angles. Hatching marks the blocked
target projection; it does not depict transmission through leaves. The figure
shows a single layer for readability. A dual-layer calculation uses the
intersection of the two openings, as described above.

## Weighting in the different entry points

| Entry point | Weight input | How it is used |
|---|---|---|
| Python/C# control-point library functions | Cumulative meterset weights | Half of each incremental MU interval at either endpoint |
| Eclipse standalone GUI | Native cumulative `ControlPoint.MetersetWeight` and beam MU | Same endpoint approximation, normalized by the final cumulative weight |
| RayStation standalone GUI | Native segment `RelativeWeight` and `BeamMU` | Segment shares of total beam MU; not interpreted as cumulative weights |

See the [Eclipse guide](../scripts/esapi/README.md#mu-weighting),
[RayStation guide](../scripts/raystation/README.md#numerical-method) and
[integration notes](tps-integration.md) for API-specific requirements.
All plan aggregation uses beam MU; the same grid spacing alone does not
establish equivalent geometry or results between TPS implementations.

## Source

Hernandez V, Lara-Aristimuño I, Abella R, Saez J. *Quantification of the aperture
modulation in radiotherapy treatment plans.* Medical Physics.
2025;52(12):e70144. [doi:10.1002/mp.70144](https://doi.org/10.1002/mp.70144).
Section 2.2 defines AM, PAM and beam-level BAM. The endpoint sampling rule
documented here is this implementation's numerical approximation.
