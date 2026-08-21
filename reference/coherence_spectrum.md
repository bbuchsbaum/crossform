# Read the coherent share of a conservative frame family against scale

A multiscale conservative frame family divides one fixed global budget
between its scales, and it divides it by the weights the analyst chose:
the `total` component summed over the rows of scale \\s\\ is exactly
\\\alpha_s G\_\Omega\\ whatever the data say
(`design/conservative-geometry-contract.md` section 3.1). What the data
do decide is how each scale's fixed budget splits into a coherent and a
configuration part. `coherence_spectrum()` reports that split, per scale
and – with `by_location = TRUE` – per (scale, location).

## Usage

``` r
coherence_spectrum(
  x,
  weights = NULL,
  by_location = FALSE,
  by = NULL,
  using = NULL
)
```

## Arguments

- x:

  A compiled `effect_geometry_plan` over a frame family, or the
  `effect_contrast_view` that
  [`contrast_energy()`](https://bbuchsbaum.github.io/crossform/reference/contrast_energy.md)
  already returned for one. Every other result kind is refused with the
  reason; see *Refusals*.

- weights:

  One finite contrast weight per experimental effect, as in
  [`contrast_energy()`](https://bbuchsbaum.github.io/crossform/reference/contrast_energy.md).
  Required with a plan, and refused with an evaluated view, which was
  already read at one contrast.

- by_location:

  `FALSE` (the default) returns one row per scale; `TRUE` returns one
  row per (scale, center), the object gap G3 of
  `design/conservative-geometry-contract.md` section 11.4 requires. It
  is the shorthand for `by = c("scale", "center")`.

- by:

  The grouping columns, when they are not the two above: one or more
  column names of the per-measurement metadata, for example
  `by = "family"` for a family whose members have no numeric scale. Pass
  `by` or `by_location`, not both.

- using:

  Optional per-measurement metadata table, joined on its `measurement`
  column so row order does not matter. A
  [`frame_family()`](https://bbuchsbaum.github.io/crossform/reference/frame_family.md)'s
  `$index` is the intended one, and it is required on the view route: a
  result's `$index` is the measurement identifier vector alone, so an
  evaluated view does not carry the scales it was read at. With a plan
  the frame's own `$index` is used unless `using` overrides it.

## Value

An `effect_contrast_view` with one row per group instead of one per
measurement – the same record
[`contribution()`](https://bbuchsbaum.github.io/crossform/reference/contribution.md)
returns, because the spectrum is that aggregation under a scale-resolved
grouping.

`$coherent` is \\E^{\mathrm{coh}}\_s\\, `$configuration` is
\\E^{\mathrm{cfg}}\_s\\, `$total` is their sum, and
`$coherence_fraction` is \\\phi_s\\, masked by
`$coherence_fraction_valid`. `$signed` is `NA`: a signed marginal is a
local weighted mean, and means do not add over a territory.

`$index` has one row per group: `measurement` holds the composite key
(`"1.01::v8"` for a (scale, center) row), the grouping columns keep
their own types, `n_rows` counts the measurements behind the row, and
any other metadata column that takes one value inside every group is
carried through – which is how `family` and `alpha` arrive next to the
energy they fix. `$metadata$aggregation` records
`reduction = "coherence_spectrum"`, `resolved_by`,
`frame_relative = TRUE`, `alpha_invariant = "coherence_fraction"`,
`alpha_fixed = "total"`, `location_collapse = "none"`, and the applied
`alpha` per group. `$receipt` is the parent receipt under a derived
`scientific_plan_id`.

## Details

The reported share is \$\$\phi_s = \frac{\sum\_{x \in s} \langle H,
G^{\mathrm{coh}}\_x\rangle} {\sum\_{x \in s} \langle H, G_x\rangle},\$\$
computed from the **aggregated** components. A fraction of sums is not a
sum of fractions, so averaging per-node fractions would answer a
different question.

## Read the share, not the energy

Every member of a conservative family is column normalized on its own
and the family weights sum to one, so each scale's rows carry exactly
their weight of the whole-domain total. **The `total` column is
therefore the `weights` vector times a constant, and a plot of energy
against scale is a plot of that vector.** The contract makes this
normative: no such panel may be presented as evidence about spatial
scale.

The share is the opposite. Both components are homogeneous of degree one
under a rescaling \\w_x \mapsto \alpha w_x\\ of a row – the total is
linear in \\w\\, and the coherent part \\K\_{\mathrm{coh}} = a a^\top /
(a^\top K_x^{-1} a)\\ has \\a = w_x / \sum w_x\\ invariant to \\\alpha\\
while \\K_x^{-1}\\ scales as \\1/\alpha\\ – so the ratio cancels
\\\alpha\\ exactly. Two families differing only in their weights give
identical shares to machine precision and different energies. The
spectrum is a property of the family's geometry, not of the weighting,
and may be reported without disclosing it (contract sections 3.1 and
3.2).

## What is frame-relative, and what is masked

`total` is budget-exact: the group totals add back to the whole-domain
total (contract section 2). `coherent` and `configuration` add up as
arithmetic, but \\\sum_x G^{\mathrm{coh}}\_x\\ is not a global quantity
(section 4, claim 4), so a coherent energy is a share of *this family's*
coherent mass and two families give two incomparable denominators. The
result carries `frame_relative = TRUE` and the print says so.

`coherence_fraction` is masked, never clamped: it is reported only where
the aggregated total is finite and positive and both aggregated
components are nonnegative, and is `NA` elsewhere. Cross-generalized
components are signed, so a scale whose nodes' common modes
anticorrelate across partitions can carry a negative coherent energy;
that scale reports `NA` rather than a number outside \\\[0,1\]\\.

One consequence is worth expecting. A singleton scale – a point member,
or a radius covering one feature – has exactly zero configuration
(`effect-form-v1` section 7), so its aggregated configuration sits
within one unit in the last place of zero and can land on either side of
it. Where it lands below, the mask fires and the share reads `NA`
instead of `1`. That is the mask working on an exactly-degenerate
partition, and it is the same behaviour a singleton node's own
`coherence_fraction` already has.

## No location-wise collapse

The coherent share is a function of (location, scale) and not a number
(gap G3). `by_location = TRUE` returns that table and stops: the scale
profile at one location is a filter on it, which needs no named
operation, while an alpha-weighted mean over scales or an argmax-scale
map is a declared reduction that would need its own certificate.
`$metadata$aggregation$location_collapse` records `"none"`.

## Refusals

A locally normalized frame signals an `effect_capability_refusal` with
capability `"conservative_frame"` in namespace `"geometry_views"`; so
does a view that does not record its frame normalization. A query-only
`effect_view`, a comparative readout (`effect_spectrum_view`,
`effect_rdm_view`, `effect_rsa_view`, `effect_crossnobis_view`), and
packed geometry each signal capability `"coherence_decomposition"`,
because none of them carries the two components a share is taken of.
Branch on them with
[`catch_refusal()`](https://bbuchsbaum.github.io/crossform/reference/catch_refusal.md).

## See also

[`searchlights()`](https://bbuchsbaum.github.io/crossform/reference/searchlights.md)
and
[`frame_family()`](https://bbuchsbaum.github.io/crossform/reference/frame_family.md)
for the families this reads,
[`contrast_energy()`](https://bbuchsbaum.github.io/crossform/reference/contrast_energy.md)
for the view it aggregates,
[`contribution()`](https://bbuchsbaum.github.io/crossform/reference/contribution.md)
for the same arithmetic over a spatial territory, and
[`frame_conservation()`](https://bbuchsbaum.github.io/crossform/reference/frame_conservation.md)
for the per-block certificate that each scale carries its own weight.

Other geometry plans and views:
[`aggregate_first()`](https://bbuchsbaum.github.io/crossform/reference/aggregate_first.md),
[`bilinear_query()`](https://bbuchsbaum.github.io/crossform/reference/bilinear_query.md),
[`compute_policy()`](https://bbuchsbaum.github.io/crossform/reference/compute_policy.md),
[`contrast_energy()`](https://bbuchsbaum.github.io/crossform/reference/contrast_energy.md),
[`contribution()`](https://bbuchsbaum.github.io/crossform/reference/contribution.md),
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
# A point effect on a line: one voxel carries the whole contrast.
n <- 15L
domain <- abstract_domain(
  n, coordinates = cbind(seq_len(n) - 1, 0),
  feature_ids = paste0("v", seq_len(n)), id = "spectrum-example"
)
pattern <- function(signal) {
  rbind(face = signal, house = rep(0, n))
}
signal <- replace(rep(0, n), 8L, 1)
relation <- relation(
  list(run1 = pattern(signal), run2 = pattern(signal)), domain = domain
)
family <- compile_frame(
  searchlights(c(0.5, 1.01, 2.01), "conservative"), domain
)
plan <- plan_geometry(
  relation, family,
  cross_partitions(relation, independence = "independent")
)

# The share falls as the neighborhood grows past the signal: one voxel of
# evidence looks entirely coherent at a scale that sees only it.
spectrum <- coherence_spectrum(plan, c(face = 1, house = -1))
as.data.frame(spectrum)[, c("scale", "alpha", "total", "coherence_fraction")]
#>   scale     alpha     total coherence_fraction
#> 1  0.50 0.3333333 0.3333333          1.0000000
#> 2  1.01 0.3333333 0.3333333          0.3333333
#> 3  2.01 0.3333333 0.3333333          0.2000000

# The energy column is the family weighting and nothing else, so it is the
# same at every scale here; the share is what varies.
spectrum$total
#> [1] 0.3333333 0.3333333 0.3333333

# One row per (scale, center) instead. A location's scale profile is a
# filter on that table, not a separate reduction.
located <- coherence_spectrum(plan, c(face = 1, house = -1),
  by_location = TRUE)
profile <- as.data.frame(located)
profile[profile$center == "v8", c("scale", "center", "coherence_fraction")]
#>    scale center coherence_fraction
#> 14  0.50     v8          1.0000000
#> 29  1.01     v8          0.3333333
#> 44  2.01     v8          0.2000000
```
