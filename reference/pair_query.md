# Describe an axis-bound pair query

A pair query is the universal rectangular readout. Its rows are bound to
one ordered left effect space and its columns to one ordered right
effect space; equal dimensions or labels do not substitute for identity.

## Usage

``` r
pair_query(H, left_space, right_space, metadata = list())
```

## Arguments

- H:

  A finite nonempty left-by-right numeric base or `Matrix` matrix.

- left_space, right_space:

  The ordered
  [`effect_space()`](https://bbuchsbaum.github.io/crossform/reference/effect_space.md)
  identities bound to the rows and columns of `H`. Unique character
  coordinates are accepted as shorthand for an unspecified-basis space.

- metadata:

  Optional compact semantic metadata, used by higher-level pair-design
  constructors for balance and design diagnostics.

## Value

An `effect_pair_query` carrying the `$operator`, its bound `$left_space`
and `$right_space`, and any `$metadata` a higher-level constructor
attached.

## See also

[`pair_lm_query()`](https://bbuchsbaum.github.io/crossform/reference/pair_lm_query.md)
and
[`coupling_contrast()`](https://bbuchsbaum.github.io/crossform/reference/coupling_contrast.md),
which compile designed pair operators;
[`bilinear_query()`](https://bbuchsbaum.github.io/crossform/reference/bilinear_query.md)
for the symmetric single-space case;
[`evaluate_geometry()`](https://bbuchsbaum.github.io/crossform/reference/evaluate_geometry.md)
for reading a rectangular plan with one.

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
[`measurement_form()`](https://bbuchsbaum.github.io/crossform/reference/measurement_form.md),
[`pair_lm_query()`](https://bbuchsbaum.github.io/crossform/reference/pair_lm_query.md)

## Examples

``` r
# An encoding-retrieval readout: 2 encoding conditions by 3 retrieval
# conditions, so the query is rectangular and cannot be a square RDM.
encoding <- effect_space(c("encode_a", "encode_b"), basis_id = "demo:enc")
retrieval <- effect_space(
  c("retrieve_a", "retrieve_b", "lure"), basis_id = "demo:ret"
)
query <- pair_query(
  matrix(c(1, -0.5, 0, 0, 0.5, -1), 2, 3, byrow = TRUE),
  encoding, retrieval
)
dim(query$operator)
#> [1] 2 3
query$left_space$coordinates
#> [1] "encode_a" "encode_b"

# Axis identity is checked, not inferred: equal dimensions do not make two
# different experimental spaces interchangeable.
swapped <- try(
  pair_query(matrix(1, 3, 2), encoding, retrieval), silent = TRUE
)
conditionMessage(attr(swapped, "condition"))
#> [1] "The effect-space dimension does not match the declared data."
```
