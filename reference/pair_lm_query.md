# Compile a weighted pair-space linear-model coefficient to an effect query

The returned matrix is the exact linear map from eligible pair values to
a requested weighted least-squares coefficient. Encoding and retrieval
nuisance effects are ordinary design columns, not a new analysis class.

## Usage

``` r
pair_lm_query(
  design,
  coefficient,
  left_space,
  right_space,
  weights = NULL,
  encoding_nuisance = FALSE,
  retrieval_nuisance = FALSE,
  sparse = FALSE
)
```

## Arguments

- design:

  A data frame with `left`, `right`, and one or more finite numeric
  predictor columns. Duplicate pair rows are allowed.

- coefficient:

  One coefficient name, or a named numeric contrast over the compiled
  design columns.

- left_space, right_space:

  Ordered effect-space identities.

- weights:

  Optional positive row weights, or the name of a design weight column.
  A `weight` column is used automatically when present.

- encoding_nuisance, retrieval_nuisance:

  Include fixed-effect nuisance columns for the respective item axis.

- sparse:

  Return `H` as a sparse `Matrix` object.

## Value

An axis-bound
[`pair_query()`](https://bbuchsbaum.github.io/crossform/reference/pair_query.md)
whose `$operator` maps eligible pair values to the requested
coefficient, with `$metadata$coefficient` (the contrast over compiled
columns) and `$metadata$diagnostics` reporting `rank`, `columns`,
`observations`, `unique_pairs`, and operator `balance`.

## See also

[`match_control()`](https://bbuchsbaum.github.io/crossform/reference/match_control.md)
for the matched-versus-control special case, and
[`pair_query()`](https://bbuchsbaum.github.io/crossform/reference/pair_query.md)
for a hand-written operator.

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
[`pair_query()`](https://bbuchsbaum.github.io/crossform/reference/pair_query.md)

## Examples

``` r
# A pair-space regression: how does encoding-retrieval similarity change
# with study-test lag, adjusting for a match indicator?
encoding <- effect_space(
  c("item1", "item2", "item3"), basis_id = "demo:encoding:v1"
)
retrieval <- effect_space(
  c("probe1", "probe2", "probe3"), basis_id = "demo:retrieval:v1"
)
design <- expand.grid(
  left = encoding$coordinates, right = retrieval$coordinates,
  stringsAsFactors = FALSE
)
design$lag <- abs(
  match(design$left, encoding$coordinates) -
    match(design$right, retrieval$coordinates)
)
design$match <- as.numeric(design$lag == 0)

# The result is the exact linear map from pair values to the `lag`
# coefficient, compiled once and reusable as a fixed query.
query <- pair_lm_query(design, "lag", encoding, retrieval)
round(as.matrix(query$operator), 3)
#>       probe1 probe2 probe3
#> item1   0.00  -0.25   0.50
#> item2  -0.25   0.00  -0.25
#> item3   0.50  -0.25   0.00
query$metadata$diagnostics[c("rank", "columns", "observations")]
#> $rank
#> [1] 3
#> 
#> $columns
#> [1] "(Intercept)" "lag"         "match"      
#> 
#> $observations
#> [1] 9
#> 

# A collinear predictor makes the coefficient undefined, and the design is
# rejected before any geometry is read.
design$lag_copy <- design$lag
refused <- try(
  pair_lm_query(design, "lag", encoding, retrieval), silent = TRUE
)
conditionMessage(attr(refused, "condition"))
#> [1] "Pair-space design is rank deficient; revise predictors or nuisances."
```
