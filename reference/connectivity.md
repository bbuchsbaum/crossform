# Request a validated connectivity view

`connectivity()` is the one entry point for the normalized views, and it
checks their preconditions before reporting a number: repeated variation
of effective rank above one, valid self-blocks, explicit regularization
for the canonical and Gaussian views, and an explicit model declaration
for Gaussian information.

## Usage

``` r
connectivity(
  x,
  view = c("correlation", "canonical", "geometry_alignment", "gaussian_information"),
  ridge = NULL,
  ridge_right = ridge,
  model = NULL,
  units = c("nats", "bits"),
  tolerance = 1e-10
)
```

## Arguments

- x:

  A complete `effect_measurement_form`.

- view:

  One of signed scalar correlation, a canonical spectrum, static
  geometry alignment, or Gaussian mutual information.

- ridge, ridge_right:

  Explicit ridge values for canonical or Gaussian views.

- model:

  A
  [`gaussian_covariance_model()`](https://bbuchsbaum.github.io/crossform/reference/gaussian_covariance_model.md)
  for Gaussian information.

- units:

  Information units, when applicable.

- tolerance:

  Positive numerical tolerance.

## Value

An `effect_coupling_result` whose `$values` data frame carries one row
per edge with the view's column (`correlation`, `canonical_correlation`
per `mode`, `geometry_alignment`, or `information` plus `units`),
alongside `$kind`, `$normalization_axis`, `$regularization`, and
`$terminology`.

## Refusals

Each precondition signals an `effect_capability_refusal` in namespace
`"coupling_views"`, so
[`catch_refusal()`](https://bbuchsbaum.github.io/crossform/reference/catch_refusal.md)
can branch on the cause instead of the prose: capability
`"certified_repeated_variation"` when the form has not established
repeated variation along a named sampling axis,
`"nondegenerate_variation"` when the variation query has effective rank
one, `"declared_regularization"` for a canonical view without `ridge`,
and `"declared_gaussian_model"` for Gaussian information without `model`
or `ridge`.

## See also

[`effect_coupling()`](https://bbuchsbaum.github.io/crossform/reference/coupling_views.md)
for the uninterpreted block,
[`gaussian_covariance_model()`](https://bbuchsbaum.github.io/crossform/reference/gaussian_covariance_model.md)
for the declaration Gaussian information requires, and
[`measurement_form()`](https://bbuchsbaum.github.io/crossform/reference/measurement_form.md)
to build `x`.

Other coupling and connectivity views:
[`control_coupling()`](https://bbuchsbaum.github.io/crossform/reference/control_coupling.md),
[`coupling()`](https://bbuchsbaum.github.io/crossform/reference/coupling.md),
[`coupling_contrast()`](https://bbuchsbaum.github.io/crossform/reference/coupling_contrast.md),
[`coupling_views`](https://bbuchsbaum.github.io/crossform/reference/coupling_views.md),
[`effect_coupling_result`](https://bbuchsbaum.github.io/crossform/reference/effect_coupling_result.md),
[`match_control()`](https://bbuchsbaum.github.io/crossform/reference/match_control.md),
[`match_coupling()`](https://bbuchsbaum.github.io/crossform/reference/match_coupling.md),
[`measurement_components()`](https://bbuchsbaum.github.io/crossform/reference/measurement_components.md),
[`measurement_form()`](https://bbuchsbaum.github.io/crossform/reference/measurement_form.md),
[`pair_lm_query()`](https://bbuchsbaum.github.io/crossform/reference/pair_lm_query.md),
[`pair_query()`](https://bbuchsbaum.github.io/crossform/reference/pair_query.md)

## Examples

``` r
# Two scalar, oriented regional measurements over six repeated time
# points in two sessions.
native <- abstract_domain(4, id = "demo:native:v1")
times <- effect_space(paste0("time", 1:6), basis_id = "demo:time:v1")
trend <- seq(-1.5, 1.5, length.out = 6)
session <- function(shift) cbind(
  trend + shift,
  -0.8 * trend + c(0.2, -0.1, 0.1, 0, -0.1, -0.1),
  sin(seq(shift, pi, length.out = 6)),
  cos(seq(shift, pi, length.out = 6))
)
signals <- relation(
  list(session1 = session(0), session2 = session(0.1)),
  effects = times, domain = native
)
nodes <- measurement_frame(
  list(anterior = matrix(c(1, 0, 0, 0), 1),
       posterior = matrix(c(0, 1, 0, 0), 1)),
  domain = native, id = "demo:regional-means:v1"
)
pairs <- expand.grid(
  from = nodes$node_ids, to = nodes$node_ids, stringsAsFactors = FALSE
)
covariance_over_time <- variation_query(
  (diag(6) - matrix(1 / 6, 6, 6)) / 5, times, "time", "joint_covariance"
)
within_session <- pairing(
  signals$partitions, signals$partitions, directed = TRUE,
  self_pairs = "allow_biased", independence = "not_independent"
)
form <- measurement_form(
  signals, edge_frame(pairs$from, pairs$to, nodes),
  covariance_over_time, within_session
)

# Both nodes are scalar and oriented, so the correlation view returns
# ordinary signed Pearson correlations. The generated data are negatively
# related, and the sign is recovered.
connectivity(form, view = "correlation")$values
#>       edge_id correlation
#> 1 edge_000001   1.0000000
#> 2 edge_000002  -0.9955404
#> 3 edge_000003  -0.9955404
#> 4 edge_000004   1.0000000

# Gaussian mutual information additionally needs an explicit model and
# information units, because it is a modeling claim, not a rescaling.
connectivity(
  form, view = "gaussian_information", ridge = 0.05,
  model = gaussian_covariance_model(
    list(assumption = "joint Gaussian time observations")
  ),
  units = "bits"
)$values
#>       edge_id information units
#> 1 edge_000001    1.869647  bits
#> 2 edge_000002    1.711892  bits
#> 3 edge_000003    1.711892  bits
#> 4 edge_000004    1.697274  bits

# Omitting the model declaration is refused rather than defaulted, and the
# refusal names the missing capability.
refusal <- catch_refusal(
  connectivity(form, view = "gaussian_information", ridge = 0.05)
)
refusal$capability
#> [1] "declared_gaussian_model"
refusal$reasons
#> [1] "gaussian_model_not_declared"
```
