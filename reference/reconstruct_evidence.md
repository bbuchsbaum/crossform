# Reconstruct or project the global neural evidence operator

Lossless reconstruction requires the complete node-edge block form and
full-column-rank identified frames. Rank-deficient frames return only an
explicitly projected operator. A workspace budget is checked before any
measurement block is read.

## Usage

``` r
reconstruct_evidence(
  x,
  between,
  method = c("auto", "parseval", "canonical_dual", "projected_pseudoinverse"),
  tolerance = 1e-10,
  max_condition = 1e+08,
  allow_projection = FALSE,
  workspace_bytes = 512 * 1024^2,
  reference_operator = NULL
)
```

## Arguments

- x:

  A frame-complete `effect_measurement_form`.

- between:

  The exact
  [`edge_frame()`](https://bbuchsbaum.github.io/crossform/reference/edge_frame.md)
  used to construct `x`.

- method:

  Reconstruction method; `"auto"` selects Parseval, canonical dual, or
  projected pseudoinverse according to frame diagnostics.

- tolerance:

  Positive relative singular-value and certification tolerance.

- max_condition:

  Maximum accepted retained condition number.

- allow_projection:

  Whether rank-deficient projected (lossy) reconstruction is admitted.
  It defaults to `FALSE` so a projection is an explicit choice, never a
  silent fallback: rank-deficient frames refuse until projection is
  explicitly admitted.

- workspace_bytes:

  Positive dense-workspace budget. The default is a hard 512 MiB ceiling
  for this explicitly small-node reconstruction path.

- reference_operator:

  Optional finite reference used only to certify the numerical
  reconstruction residual.

## Value

An `effect_tomography_result` with the reconstructed `$operator`, the
`$method` actually used (`"parseval"`, `"canonical_dual"`, or
`"projected_pseudoinverse"`), a `$status` distinguishing exact,
certified, and projected reconstruction, the `$lossless` and
`$certified` flags, `$left_projection`/`$right_projection`, and frame
`$diagnostics`.

## Refusal

A form that is not frame complete — a diagonal-only or requested-edge
map — signals an `effect_capability_refusal` carrying capability
`"complete_edge_set"` in namespace `"tomography"`, with reason
`"edge_set_is_not_frame_complete"`. Inspect it with
[`catch_refusal()`](https://bbuchsbaum.github.io/crossform/reference/catch_refusal.md).

## See also

[`measurement_form()`](https://bbuchsbaum.github.io/crossform/reference/measurement_form.md)
and
[`edge_frame()`](https://bbuchsbaum.github.io/crossform/reference/edge_frame.md),
which must supply every directed node pair for a lossless claim.

Other sampling uncertainty:
[`rdm_sampling_covariance()`](https://bbuchsbaum.github.io/crossform/reference/rdm_sampling_covariance.md),
[`sampling_capabilities()`](https://bbuchsbaum.github.io/crossform/reference/sampling_capabilities.md),
[`sampling_covariance()`](https://bbuchsbaum.github.io/crossform/reference/sampling_covariance.md)

## Examples

``` r
# A Parseval frame: stacking the two node operators gives the identity, so
# the local blocks add back with no correction matrix.
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
halves <- measurement_frame(
  list(first_half = diag(4)[1:2, , drop = FALSE],
       second_half = diag(4)[3:4, , drop = FALSE]),
  native, id = "demo:parseval:v1"
)
pairs <- expand.grid(
  from = halves$node_ids, to = halves$node_ids, stringsAsFactors = FALSE
)
between <- edge_frame(pairs$from, pairs$to, halves)
form <- measurement_form(
  signals, between,
  variation_query(
    (diag(6) - matrix(1 / 6, 6, 6)) / 5, times, "time", "joint_covariance"
  ),
  pairing(
    signals$partitions, signals$partitions, directed = TRUE,
    self_pairs = "allow_biased", independence = "not_independent"
  )
)

# Every directed edge is present and the frame has full column rank, so
# the global 4-by-4 neural operator is recovered exactly.
reconstructed <- reconstruct_evidence(form, between)
c(method = reconstructed$method, status = reconstructed$status,
  lossless = reconstructed$lossless)
#>                           method                           status 
#>                       "parseval" "exact_algebraic_reconstruction" 
#>                         lossless 
#>                           "TRUE" 
round(reconstructed$operator, 3)
#>        [,1]   [,2]   [,3]   [,4]
#> [1,]  1.260 -1.104 -0.023 -0.924
#> [2,] -1.104  0.976  0.016  0.805
#> [3,] -0.023  0.016  0.176  0.013
#> [4,] -0.924  0.805  0.013  0.692

# Diagonal node blocks alone are not enough: dropping the cross-edges
# discards the between-node directions, and the lossless claim is refused.
self_only <- edge_frame(
  halves$node_ids, halves$node_ids, halves
)
partial <- measurement_form(
  signals, self_only,
  variation_query(
    (diag(6) - matrix(1 / 6, 6, 6)) / 5, times, "time", "joint_covariance"
  ),
  pairing(
    signals$partitions, signals$partitions, directed = TRUE,
    self_pairs = "allow_biased", independence = "not_independent"
  )
)
refusal <- catch_refusal(reconstruct_evidence(partial, self_only))
refusal$capability
#> [1] "complete_edge_set"
refusal$remedies
#> [1] "Build the form over every directed node pair, for example with `edge_frame()` on `expand.grid(from = nodes$node_ids, to = nodes$node_ids)`."
```
