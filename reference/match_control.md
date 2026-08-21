# Compile a matched-versus-control pair-space coefficient

This convenience uses the eligible pair set carried by `matches`,
includes both item nuisance families by default, and reports the balance
of the final operator. It does not claim invariance after any later
external weighting.

## Usage

``` r
match_control(
  matches,
  weights = NULL,
  encoding_nuisance = TRUE,
  retrieval_nuisance = TRUE
)
```

## Arguments

- matches:

  A
  [`match_coupling()`](https://bbuchsbaum.github.io/crossform/reference/match_coupling.md)
  value.

- weights:

  Optional positive eligible-pair weights.

- encoding_nuisance, retrieval_nuisance:

  Nuisance-effect flags.

## Value

An axis-bound
[`pair_query()`](https://bbuchsbaum.github.io/crossform/reference/pair_query.md)
for the matched-versus-control coefficient, with `$metadata$constructor`
set to `"match_control"`, `$metadata$diagnostics` (design `rank`,
`columns`, operator `balance`), and a `$metadata$claim` naming exactly
what the reported balance covers.

## See also

[`match_coupling()`](https://bbuchsbaum.github.io/crossform/reference/match_coupling.md)
for the input,
[`coupling_contrast()`](https://bbuchsbaum.github.io/crossform/reference/coupling_contrast.md)
for the plain difference without nuisance effects, and
[`pair_lm_query()`](https://bbuchsbaum.github.io/crossform/reference/pair_lm_query.md)
for the general designed coefficient.

Other coupling and connectivity views:
[`connectivity()`](https://bbuchsbaum.github.io/crossform/reference/connectivity.md),
[`control_coupling()`](https://bbuchsbaum.github.io/crossform/reference/control_coupling.md),
[`coupling()`](https://bbuchsbaum.github.io/crossform/reference/coupling.md),
[`coupling_contrast()`](https://bbuchsbaum.github.io/crossform/reference/coupling_contrast.md),
[`coupling_views`](https://bbuchsbaum.github.io/crossform/reference/coupling_views.md),
[`effect_coupling_result`](https://bbuchsbaum.github.io/crossform/reference/effect_coupling_result.md),
[`match_coupling()`](https://bbuchsbaum.github.io/crossform/reference/match_coupling.md),
[`measurement_components()`](https://bbuchsbaum.github.io/crossform/reference/measurement_components.md),
[`measurement_form()`](https://bbuchsbaum.github.io/crossform/reference/measurement_form.md),
[`pair_lm_query()`](https://bbuchsbaum.github.io/crossform/reference/pair_lm_query.md),
[`pair_query()`](https://bbuchsbaum.github.io/crossform/reference/pair_query.md)

## Examples

``` r
encoding <- effect_space(
  c("item1", "item2", "item3"), basis_id = "demo:encoding:v1"
)
retrieval <- effect_space(
  c("probe1", "probe2", "probe3"), basis_id = "demo:retrieval:v1"
)
matches <- match_coupling(
  c("item1", "item2", "item3"), c("probe1", "probe2", "probe3"),
  encoding, retrieval
)

# Both item nuisance families are included by default, so item-specific
# offsets cannot masquerade as a matching effect.
query <- match_control(matches)
query$metadata$diagnostics$columns
#> [1] "(Intercept)"                    "match"                         
#> [3] "encoding:factor(left_index)2"   "encoding:factor(left_index)3"  
#> [5] "retrieval:factor(right_index)2" "retrieval:factor(right_index)3"
round(as.matrix(query$operator), 3)
#>       probe1 probe2 probe3
#> item1  0.333 -0.167 -0.167
#> item2 -0.167  0.333 -0.167
#> item3 -0.167 -0.167  0.333

# The result records what its balance diagnostics actually claim.
query$metadata$claim
#> [1] "observed_operator_balance_not_downstream_invariance"

# Dropping both nuisance families gives the simpler unadjusted contrast.
unadjusted <- match_control(
  matches, encoding_nuisance = FALSE, retrieval_nuisance = FALSE
)
unadjusted$metadata$diagnostics$columns
#> [1] "(Intercept)" "match"      
```
