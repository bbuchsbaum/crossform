# Declare a repeated-variation experimental query

A variation query binds a positive semidefinite operator to one
identified experimental space and records the sample axis that gives
coupling its repeated variation. A rank-one operator remains a valid
effect query, but normalized connectivity views reject it as degenerate.

## Usage

``` r
variation_query(
  operator,
  effects,
  sampling_axis,
  construction = c("psd_variation", "joint_covariance"),
  provenance = list()
)
```

## Arguments

- operator:

  A finite symmetric positive-semidefinite matrix.

- effects:

  The exact `effect_space` on both query axes.

- sampling_axis:

  Identifier such as `"time"`, `"trial"`, or `"subject"`.

- construction:

  `"psd_variation"` for a positive variation form, or
  `"joint_covariance"` when the later partition construction is a
  coherent joint covariance.

- provenance:

  Named fixed-construction metadata.

## Value

An `effect_pair_query` carrying the `$operator`, its bound `$left_space`
and `$right_space`, and `$metadata$evidence_capability` recording the
`sampling_axis` and `construction` that later views check.

## See also

[`measurement_form()`](https://bbuchsbaum.github.io/crossform/reference/measurement_form.md)
and
[`coupling()`](https://bbuchsbaum.github.io/crossform/reference/coupling.md),
which take this as `by`;
[`connectivity()`](https://bbuchsbaum.github.io/crossform/reference/connectivity.md),
which refuses a rank-one variation query;
[`pair_query()`](https://bbuchsbaum.github.io/crossform/reference/pair_query.md)
for a fixed query making no variation claim.

Other geometry plans and views:
[`aggregate_first()`](https://bbuchsbaum.github.io/crossform/reference/aggregate_first.md),
[`bilinear_query()`](https://bbuchsbaum.github.io/crossform/reference/bilinear_query.md),
[`coherence_spectrum()`](https://bbuchsbaum.github.io/crossform/reference/coherence_spectrum.md),
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
[`rsa()`](https://bbuchsbaum.github.io/crossform/reference/rsa.md)

## Examples

``` r
# Centering eight repeated time points and dividing by n - 1 is the
# within-session covariance operation, declared on the time axis.
times <- effect_space(paste0("time", 1:8), basis_id = "demo:time:v1")
center <- diag(8) - matrix(1 / 8, 8, 8)
sample_covariance <- variation_query(
  center / 7, times,
  sampling_axis = "time", construction = "joint_covariance",
  provenance = list(estimator = "centered within session")
)
sample_covariance$metadata$evidence_capability$sampling_axis
#> [1] "time"
sample_covariance$metadata$evidence_capability$construction
#> [1] "joint_covariance"

# Centering removes one direction, so eight time points leave rank seven:
# enough repeated variation for a normalized connectivity view.
qr(as.matrix(sample_covariance$operator))$rank
#> [1] 7

# A contrast gives the rank-one query `c c'`. It is a valid effect query,
# but it retains no repeated variation for connectivity to normalize.
direction <- rep(c(-1, 1), each = 4)
rank_one <- variation_query(
  tcrossprod(direction), times, "time", "joint_covariance"
)
qr(as.matrix(rank_one$operator))$rank
#> [1] 1

# The operator must be positive semidefinite to be called variation.
refused <- try(
  variation_query(diag(c(1, -1, rep(1, 6))), times, "time"), silent = TRUE
)
conditionMessage(attr(refused, "condition"))
#> [1] "A variation query construction must be positive semidefinite."
```
