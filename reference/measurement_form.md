# Materialize requested measurement-space forms

`measurement_form()` closes the experimental boundary with `by` and
leaves only the explicitly requested neural edges in `between` open.
`left` and `right` always name experimental-neural relations; spatial
endpoints are contained inside `between`.

## Usage

``` r
measurement_form(
  left,
  between,
  by,
  over,
  right = left,
  reducer = aggregate_first(),
  storage = c("memory", "block"),
  storage_path = NULL,
  compute = compute_policy(),
  route = c("auto", "forward_k", "pull_h", "factorized_h", "scalar_stack",
    "multivariate_blocks")
)
```

## Arguments

- left:

  The left `effect_relation`.

- between:

  An explicit
  [`edge_frame()`](https://bbuchsbaum.github.io/crossform/reference/edge_frame.md).

- by:

  An axis-bound
  [`pair_query()`](https://bbuchsbaum.github.io/crossform/reference/pair_query.md)
  or
  [`variation_query()`](https://bbuchsbaum.github.io/crossform/reference/variation_query.md).

- over:

  An explicit partition
  [`pairing()`](https://bbuchsbaum.github.io/crossform/reference/pairing.md).

- right:

  The right `effect_relation`; defaults to `left`.

- reducer:

  A partition-order declaration. Raw forms are bilinear, but the choice
  remains part of scientific plan identity.

- storage:

  Either in-memory or block-backed storage.

- storage_path:

  Required for block-backed storage.

- compute:

  A
  [`compute_policy()`](https://bbuchsbaum.github.io/crossform/reference/compute_policy.md).

- route:

  Internal contraction route. `"auto"` chooses an equivalent bounded
  route without changing the scientific plan.

## Value

An `effect_measurement_form` over the requested edge set: `$block_index`
names one row per edge (`edge_id`, `left`, `right`, block dimensions),
`$diagnostics` reports the `experimental_effective_rank` that gates the
normalized views, and `$capabilities`, `$plan`, and `$receipt` record
what may be claimed. Read the blocks with
[`effect_coupling()`](https://bbuchsbaum.github.io/crossform/reference/coupling_views.md)
and the other views.

## Details

Raw forms are bilinear. Covariance, correlation, CCA, geometry
alignment, and Gaussian information are separate capability-gated views.
The default reducer records that raw partition forms are aggregated
before any later nonlinear connectivity normalization.

## See also

[`effect_coupling()`](https://bbuchsbaum.github.io/crossform/reference/coupling_views.md),
[`connectivity()`](https://bbuchsbaum.github.io/crossform/reference/connectivity.md),
[`measurement_components()`](https://bbuchsbaum.github.io/crossform/reference/measurement_components.md),
and
[`reconstruct_evidence()`](https://bbuchsbaum.github.io/crossform/reference/reconstruct_evidence.md)
for the views;
[`coupling()`](https://bbuchsbaum.github.io/crossform/reference/coupling.md)
to derive the same form from an existing geometry plan.

Other coupling and connectivity views:
[`connectivity()`](https://bbuchsbaum.github.io/crossform/reference/connectivity.md),
[`control_coupling()`](https://bbuchsbaum.github.io/crossform/reference/control_coupling.md),
[`coupling()`](https://bbuchsbaum.github.io/crossform/reference/coupling.md),
[`coupling_contrast()`](https://bbuchsbaum.github.io/crossform/reference/coupling_contrast.md),
[`coupling_views`](https://bbuchsbaum.github.io/crossform/reference/coupling_views.md),
[`effect_coupling_result`](https://bbuchsbaum.github.io/crossform/reference/effect_coupling_result.md),
[`match_control()`](https://bbuchsbaum.github.io/crossform/reference/match_control.md),
[`match_coupling()`](https://bbuchsbaum.github.io/crossform/reference/match_coupling.md),
[`measurement_components()`](https://bbuchsbaum.github.io/crossform/reference/measurement_components.md),
[`pair_lm_query()`](https://bbuchsbaum.github.io/crossform/reference/pair_lm_query.md),
[`pair_query()`](https://bbuchsbaum.github.io/crossform/reference/pair_query.md)

## Examples

``` r
# Two sessions of six repeated time points over four native features.
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

# Two scalar node measurements and all four directed node pairs.
nodes <- measurement_frame(
  list(anterior = matrix(c(1, 0, 0, 0), 1),
       posterior = matrix(c(0, 1, 0, 0), 1)),
  domain = native, id = "demo:regional-means:v1"
)
pairs <- expand.grid(
  from = nodes$node_ids, to = nodes$node_ids, stringsAsFactors = FALSE
)

# `by` closes the experimental axis (within-session covariance over time)
# and `over` says each session multiplies by itself, which is biased and
# must be declared as such.
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
form$block_index[, c("edge_id", "left", "right")]
#> <effect_measurement_block_index>
#>   blocks:   4
#>   edges:    edge_000001, edge_000002, edge_000003 (+1 more)
#>   left:     anterior, posterior
#>   right:    anterior, posterior
#>   widths:   none x none
#>   elements: 0
#>   columns:  edge_id, left, right

# Centering leaves five directions of repeated variation, which is what
# the normalized views require.
form$diagnostics$experimental_effective_rank
#> [1] 5
effect_coupling(form)$values[["edge_000002"]]
#>        [,1]
#> [1,] -1.104
```
