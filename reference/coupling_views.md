# Interpret completed measurement forms

These functions are views over one completed measurement form, not
separate fitting engines. `effect_coupling()` makes no covariance claim.
`covariance_coupling()`, `canonical_coupling()`, and
[`connectivity()`](https://bbuchsbaum.github.io/crossform/reference/connectivity.md)
require certified repeated variation; normalized views additionally
require coherent positive self-covariances. `geometry_alignment()` is
static linear CKA/RV-like alignment and does not redefine established
dynamic informational-connectivity analyses.

## Usage

``` r
effect_coupling(x)

covariance_coupling(x, tolerance = 1e-10)

canonical_coupling(x, ridge, ridge_right = ridge, tolerance = 1e-10)

geometry_alignment(x, tolerance = 1e-10)
```

## Arguments

- x:

  A complete `effect_measurement_form`.

- tolerance:

  Positive numerical tolerance.

- ridge, ridge_right:

  Positive ridge values for the left and right covariance blocks.

## Value

An `effect_coupling_result`. `$values` is one entry per edge (a matrix
block for `effect_coupling()` and `covariance_coupling()`, a data frame
of `canonical_correlation` per mode for `canonical_coupling()`, and a
data frame of `geometry_alignment` for `geometry_alignment()`), with
`$edge_index` naming the edges and `$kind`, `$terminology`,
`$normalization_axis`, and `$regularization` recording what is claimed.

## Refusals

`effect_coupling()` makes no covariance claim and so refuses nothing.
Every other view signals an `effect_capability_refusal` in namespace
`"coupling_views"`, so
[`catch_refusal()`](https://bbuchsbaum.github.io/crossform/reference/catch_refusal.md)
can branch on the cause rather than the prose:

- `"certified_repeated_variation"` (reason
  `"repeated_variation_not_certified"`) when the form has not
  established that its measurement axis repeats;

- `"nondegenerate_variation"` (reason `"rank_one_variation_axis"`) when
  the variation query has effective rank one;

- `"coherent_joint_covariance"`, carrying
  `"joint_covariance_not_certified"` and `"self_blocks_not_validated"`
  as applicable, when the form does not certify a joint covariance
  across both measurement spaces;

- `"nondegenerate_self_variance"` (reason
  `"self_variance_not_strictly_positive"`) when a scalar edge has no
  measured variation to divide by.

`canonical_coupling()` additionally refuses with
`"declared_regularization"` when `ridge` is absent.

## See also

[`measurement_form()`](https://bbuchsbaum.github.io/crossform/reference/measurement_form.md)
and
[`coupling()`](https://bbuchsbaum.github.io/crossform/reference/coupling.md)
to build `x`;
[`connectivity()`](https://bbuchsbaum.github.io/crossform/reference/connectivity.md)
for the same views behind one capability-checked entry point;
[`measurement_components()`](https://bbuchsbaum.github.io/crossform/reference/measurement_components.md)
for decomposed nodes.

Other coupling and connectivity views:
[`connectivity()`](https://bbuchsbaum.github.io/crossform/reference/connectivity.md),
[`control_coupling()`](https://bbuchsbaum.github.io/crossform/reference/control_coupling.md),
[`coupling()`](https://bbuchsbaum.github.io/crossform/reference/coupling.md),
[`coupling_contrast()`](https://bbuchsbaum.github.io/crossform/reference/coupling_contrast.md),
[`effect_coupling_result`](https://bbuchsbaum.github.io/crossform/reference/effect_coupling_result.md),
[`match_control()`](https://bbuchsbaum.github.io/crossform/reference/match_control.md),
[`match_coupling()`](https://bbuchsbaum.github.io/crossform/reference/match_coupling.md),
[`measurement_components()`](https://bbuchsbaum.github.io/crossform/reference/measurement_components.md),
[`measurement_form()`](https://bbuchsbaum.github.io/crossform/reference/measurement_form.md),
[`pair_lm_query()`](https://bbuchsbaum.github.io/crossform/reference/pair_lm_query.md),
[`pair_query()`](https://bbuchsbaum.github.io/crossform/reference/pair_query.md)

## Examples

``` r
# Two sessions of six repeated time points over four native features,
# measured at two multivariate populations.
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
populations <- measurement_frame(
  list(anterior = diag(4)[1:2, , drop = FALSE],
       posterior = diag(4)[3:4, , drop = FALSE]),
  domain = native, id = "demo:populations:v1"
)
pairs <- expand.grid(
  from = populations$node_ids, to = populations$node_ids,
  stringsAsFactors = FALSE
)
form <- measurement_form(
  signals, edge_frame(pairs$from, pairs$to, populations),
  variation_query(
    (diag(6) - matrix(1 / 6, 6, 6)) / 5, times, "time", "joint_covariance"
  ),
  pairing(
    signals$partitions, signals$partitions, directed = TRUE,
    self_pairs = "allow_biased", independence = "not_independent"
  )
)
cross <- form$block_index$edge_id[
  form$block_index$left == "anterior" &
    form$block_index$right == "posterior"
]

# The raw block, with no covariance claim attached to it.
effect_coupling(form)$values[[cross]]
#>             [,1]       [,2]
#> [1,] -0.02298934 -0.9235677
#> [2,]  0.01560129  0.8049703

# The same block, now certified as repeated-sample covariance.
covariance_coupling(form)$kind
#> [1] "covariance_coupling"

# Canonical correlations describe the shared modes in descending order.
# `ridge` is recorded because changing it changes the reported values.
canonical <- canonical_coupling(form, ridge = 0.05)
canonical$values[canonical$values$edge_id == cross, ]
#>       edge_id mode canonical_correlation
#> 5 edge_000003    1            0.94217816
#> 6 edge_000003    2            0.03001232

# Geometry alignment asks a different question: do the two populations
# induce similar geometry over the repeated observations?
alignment <- geometry_alignment(form)
alignment$values$geometry_alignment[alignment$values$edge_id == cross]
#> [1] 0.9423914
```
