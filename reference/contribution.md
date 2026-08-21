# Add a conservative attribution map up over a territory

A conservative frame partitions one fixed global budget: each node's
value is its share, not its density, so the shares of a territory add up
(`design/conservative-geometry-contract.md` sections 1.2 and 2).
`contribution()` performs exactly that addition, and nothing else. It is
the ledger reading of an attribution map, and it refuses a detection
map, whose overlapping nodes double-count and whose sum estimates
nothing.

## Usage

``` r
contribution(x, ...)

# Default S3 method
contribution(x, by, using = NULL, ...)
```

## Arguments

- x:

  A view whose per-measurement values are additive contributions: an
  `effect_contrast_view` from
  [`contrast_energy()`](https://bbuchsbaum.github.io/crossform/reference/contrast_energy.md),
  or a query-only `effect_view` from
  [`evaluate_geometry()`](https://bbuchsbaum.github.io/crossform/reference/evaluate_geometry.md)
  or
  [`query_geometry()`](https://bbuchsbaum.github.io/crossform/reference/query_geometry.md).
  Every other result kind is refused with the reason; see *Refusals*.

- ...:

  Passed to the method. The generic dispatches on `x`: a per-measurement
  view is aggregated by the default method documented here, and an
  `effect_population_result` or `effect_population_view` by the
  group-node method in
  [population_views](https://bbuchsbaum.github.io/crossform/reference/population_views.md).

- by:

  The grouping. Either one column name of the per-measurement metadata
  (searched in `using` first, then in `x$index` when that is a table),
  or the grouping itself as a vector or factor with one entry per
  measurement – a region or network label, say. A length-one character
  value is always read as a column name. Every measurement must carry a
  group; an `NA` label is refused rather than dropped.

- using:

  Optional per-measurement metadata table naming the groups a view
  cannot carry itself. A
  [`frame_family()`](https://bbuchsbaum.github.io/crossform/reference/frame_family.md)'s
  `$index` is the intended one: it is joined on its `measurement`
  column, so row order does not matter. A table without that column must
  have one row per measurement, in `$index` order.

## Value

The same class as `x`, with one row per group instead of one per
measurement.

`$index` has one row per group: `measurement`, holding the group key,
and `n_rows`, the number of measurements that went into it. A key that
[`as.character()`](https://rdrr.io/r/base/character.html) would flatten
– a numeric `scale`, say – also gets a typed column under the grouping's
own name. Groups appear in sorted key order, or in level order when `by`
is a factor; an empty group is dropped. `$metadata$aggregation` records
`aggregated_by`, the group count, `frame_relative`, and which components
are budget-exact. `$receipt` is the parent receipt under a derived
`scientific_plan_id`.

## Details

Grouping is **by row**: each row of `x` belongs to exactly one group, so
the group totals partition the whole-domain total exactly. Splitting an
overlapping node's mass across several territories is deliberately not
offered – a row is one node, a node is not divisible, and any split
needs a second partition of the frame weights that can double-count.

## What conserves and what does not

`total` is budget-exact: the group totals add back to the whole-domain
total, which is the value the same contrast takes under
`whole_brain("none")` (contract section 2, claim 2).

`signed` is **not**, and the aggregate reports `NA` for it. A contrast
view's signed marginal is the local weighted *mean* contrast – already
divided by the node's own frame mass – so it is a density, and adding
densities over a territory is the error a conservative frame exists to
avoid (contract section 1.1). The defensible aggregate is a
mass-weighted mean, which needs the node masses a view does not carry,
so the field is masked rather than filled with a number nobody can
interpret.

`coherent` and `configuration` add up as arithmetic, but their sums are
**frame-relative** (contract section 4, claim 4): \\\sum_x G_x^{coh}\\
is not a global quantity, so a coherent budget is a share of *this
frame's* coherent mass and two frames give two incomparable
denominators. The result carries `frame_relative = TRUE` and the print
says so. Never compare a coherent budget across frames.

`coherence_fraction` is recomputed from the aggregated components, not
averaged: a fraction of sums is not a sum of fractions. It is masked by
the same rule a node's fraction is masked by – reported only where the
aggregated total is finite and positive and both components are
nonnegative – so a group whose components do not form a nonnegative
partition reports `NA` rather than a clamped number.

## Refusals

A locally normalized frame signals an `effect_capability_refusal` with
capability `"conservative_frame"` in namespace `"geometry_views"`; so
does a view that does not record its frame normalization. A view whose
values are not additive contributions – an `effect_spectrum_view`,
`effect_rdm_view`, `effect_rsa_view`, `effect_crossnobis_view` – a
stored `effect_geometry`, and an unevaluated `effect_geometry_plan` each
signal capability `"additive_contribution"`. Branch on them with
[`catch_refusal()`](https://bbuchsbaum.github.io/crossform/reference/catch_refusal.md).

## See also

[`contrast_energy()`](https://bbuchsbaum.github.io/crossform/reference/contrast_energy.md)
and
[`evaluate_geometry()`](https://bbuchsbaum.github.io/crossform/reference/evaluate_geometry.md)
for the views this aggregates,
[`frame_family()`](https://bbuchsbaum.github.io/crossform/reference/frame_family.md)
for the `$index` that names scales and members, and
[`frame_conservation()`](https://bbuchsbaum.github.io/crossform/reference/frame_conservation.md)
for the certificate that the frame conserves the budget being divided.

Other geometry plans and views:
[`aggregate_first()`](https://bbuchsbaum.github.io/crossform/reference/aggregate_first.md),
[`bilinear_query()`](https://bbuchsbaum.github.io/crossform/reference/bilinear_query.md),
[`coherence_spectrum()`](https://bbuchsbaum.github.io/crossform/reference/coherence_spectrum.md),
[`compute_policy()`](https://bbuchsbaum.github.io/crossform/reference/compute_policy.md),
[`contrast_energy()`](https://bbuchsbaum.github.io/crossform/reference/contrast_energy.md),
[`crossnobis()`](https://bbuchsbaum.github.io/crossform/reference/crossnobis.md),
[`evaluate_geometry()`](https://bbuchsbaum.github.io/crossform/reference/evaluate_geometry.md),
[`example_fmri_effects()`](https://bbuchsbaum.github.io/crossform/reference/example_fmri_effects.md),
[`geometry_component()`](https://bbuchsbaum.github.io/crossform/reference/geometry_component.md),
[`geometry_spectrum()`](https://bbuchsbaum.github.io/crossform/reference/geometry_spectrum.md),
[`latent_geometry()`](https://bbuchsbaum.github.io/crossform/reference/latent_geometry.md),
[`materialize_geometry()`](https://bbuchsbaum.github.io/crossform/reference/materialize_geometry.md),
[`plan_crossnobis()`](https://bbuchsbaum.github.io/crossform/reference/plan_crossnobis.md),
[`plan_geometry()`](https://bbuchsbaum.github.io/crossform/reference/plan_geometry.md),
[`plot_views`](https://bbuchsbaum.github.io/crossform/reference/plot_views.md),
[`query_geometry()`](https://bbuchsbaum.github.io/crossform/reference/query_geometry.md),
[`rdm()`](https://bbuchsbaum.github.io/crossform/reference/rdm.md),
[`rsa()`](https://bbuchsbaum.github.io/crossform/reference/rsa.md),
[`variation_query()`](https://bbuchsbaum.github.io/crossform/reference/variation_query.md)

## Examples

``` r
# Four regions over eight voxels, read as an attribution map.
domain <- abstract_domain(8, id = "contribution-example")
run1 <- rbind(
  face = c(1, 0.2, 0, 0, 0.5, 0.1, 0, 0),
  house = c(0, 1, 0.1, 0, 0, 0.6, 0.2, 0)
)
run2 <- rbind(
  face = c(0.9, 0.3, 0, 0, 0.4, 0.2, 0, 0),
  house = c(0.1, 0.9, 0, 0, 0.1, 0.5, 0.3, 0)
)
relation <- relation(list(run1 = run1, run2 = run2), domain = domain)
labels <- rep(c("r1", "r2", "r3", "r4"), each = 2)
plan <- plan_geometry(
  relation, compile_frame(regions(labels, "conservative"), domain),
  cross_partitions(relation, independence = "independent")
)
effect <- contrast_energy(plan, c(face = 1, house = -1))

# Two networks, each two regions. The ledger adds up exactly.
network <- c("early", "early", "late", "late")
ledger <- contribution(effect, by = network)
ledger
#> <effect_contrast_view>
#>   measurements: 2
#>   contrast:     face 1, house -1
#>  measurement n_rows signed coherent configuration total coherence_fraction
#>        early      2     NA     0.02          1.26  1.28            0.01562
#>         late      2     NA     0.03          0.33  0.36            0.08333
#>   coherence_fraction: 2 of 2 valid; NA where coherent and configuration are
#>     not a nonnegative partition
#>   aggregated_by: `network`; 2 groups over 4 measurements, summed by row
#>   frame_relative: TRUE -- `coherent`, `configuration` are a share of this
#>     frame's own mass, not of a global one, and are not comparable across
#>     frames
#>   masked: `signed` NA over groups; a local weighted mean is a density, and
#>     densities do not add over a territory
sum(ledger$total) - sum(effect$total)
#> [1] 0

# A detection map is refused: overlapping local nodes double-count, so
# their sum is not a share of anything.
detection <- contrast_energy(
  plan_geometry(
    relation, compile_frame(regions(labels), domain),
    cross_partitions(relation, independence = "independent")
  ),
  c(face = 1, house = -1)
)
refusal <- catch_refusal(contribution(detection, by = network))
refusal$capability
#> [1] "conservative_frame"
```
