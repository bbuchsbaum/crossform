# Summarize a crossed node decomposition on one measurement edge

The returned strengths are invariant to rotations of subspace-only
components. Raw configuration entries are not exposed as basis-free
scientific quantities.

## Usage

``` r
measurement_components(x, edge)
```

## Arguments

- x:

  A complete `effect_measurement_form` whose node measurements carry
  decompositions.

- edge:

  An edge number or edge identifier.

## Value

A data frame with one row per crossed component, carrying
`left_component`/`right_component`, the block dimensions `d_left` and
`d_right`, `left_orientation`/`right_orientation`,
`raw_entries_meaningful` (true only when both sides are oriented),
`frobenius_strength`, and `strongest_singular_value`.

## See also

[`measurement_frame()`](https://bbuchsbaum.github.io/crossform/reference/measurement_frame.md)
with `mode = "coherent_configuration"`, which creates the decomposition,
and
[`effect_coupling()`](https://bbuchsbaum.github.io/crossform/reference/coupling_views.md)
for the undecomposed block.

Other coupling and connectivity views:
[`connectivity()`](https://bbuchsbaum.github.io/crossform/reference/connectivity.md),
[`control_coupling()`](https://bbuchsbaum.github.io/crossform/reference/control_coupling.md),
[`coupling()`](https://bbuchsbaum.github.io/crossform/reference/coupling.md),
[`coupling_contrast()`](https://bbuchsbaum.github.io/crossform/reference/coupling_contrast.md),
[`coupling_views`](https://bbuchsbaum.github.io/crossform/reference/coupling_views.md),
[`effect_coupling_result`](https://bbuchsbaum.github.io/crossform/reference/effect_coupling_result.md),
[`match_control()`](https://bbuchsbaum.github.io/crossform/reference/match_control.md),
[`match_coupling()`](https://bbuchsbaum.github.io/crossform/reference/match_coupling.md),
[`measurement_form()`](https://bbuchsbaum.github.io/crossform/reference/measurement_form.md),
[`pair_lm_query()`](https://bbuchsbaum.github.io/crossform/reference/pair_lm_query.md),
[`pair_query()`](https://bbuchsbaum.github.io/crossform/reference/pair_query.md)

## Examples

``` r
# Two overlapping additive measurements, each split into its
# weighted-mean (coherent) direction and the orthogonal remainder.
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
decomposed <- measurement_frame(
  additive_frame(
    matrix(c(1, 2, 1, 0, 0, 1, 2, 1), 2, 4, byrow = TRUE), domain = native
  ),
  mode = "coherent_configuration"
)
form <- measurement_form(
  signals,
  edge_frame(
    decomposed$node_ids[1], decomposed$node_ids[2], decomposed
  ),
  variation_query(
    (diag(6) - matrix(1 / 6, 6, 6)) / 5, times, "time", "joint_covariance"
  ),
  pairing(
    signals$partitions, signals$partitions, directed = TRUE,
    self_pairs = "allow_biased", independence = "not_independent"
  )
)

# Four crossed components on one edge. Only coherent-to-coherent has a
# fixed orientation on both sides, so only its raw entry is meaningful;
# the others are reported as rotation-invariant strengths.
components <- measurement_components(form, edge = 1)
components[, c(
  "left_component", "right_component", "raw_entries_meaningful",
  "frobenius_strength"
)]
#>   left_component right_component raw_entries_meaningful frobenius_strength
#> 1       coherent        coherent                   TRUE          0.4828991
#> 2  configuration        coherent                  FALSE          1.5874332
#> 3       coherent   configuration                  FALSE          0.3041557
#> 4  configuration   configuration                  FALSE          1.5785237
```
