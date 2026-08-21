# Mark matched encoding-retrieval effect pairs

Duplicate rows are retained as multiplicity, so repeated retrievals are
explicit rather than silently deduplicated.

## Usage

``` r
match_coupling(left, right, left_space, right_space, eligible = NULL)
```

## Arguments

- left, right:

  Matched coordinate identifiers of equal positive length.

- left_space, right_space:

  Ordered effect-space identities.

- eligible:

  Optional restricted eligible-pair set.

## Value

An `effect_pair_coupling` with `kind = "match"`, a nonnegative `$value`
matrix counting the multiplicity of each matched cell, the logical
`$eligible` mask, and the bound `$left_space`/`$right_space`.

## See also

[`control_coupling()`](https://bbuchsbaum.github.io/crossform/reference/control_coupling.md)
for the complementary cells,
[`coupling_contrast()`](https://bbuchsbaum.github.io/crossform/reference/coupling_contrast.md)
and
[`match_control()`](https://bbuchsbaum.github.io/crossform/reference/match_control.md)
to turn the pair into a query.

Other coupling and connectivity views:
[`connectivity()`](https://bbuchsbaum.github.io/crossform/reference/connectivity.md),
[`control_coupling()`](https://bbuchsbaum.github.io/crossform/reference/control_coupling.md),
[`coupling()`](https://bbuchsbaum.github.io/crossform/reference/coupling.md),
[`coupling_contrast()`](https://bbuchsbaum.github.io/crossform/reference/coupling_contrast.md),
[`coupling_views`](https://bbuchsbaum.github.io/crossform/reference/coupling_views.md),
[`effect_coupling_result`](https://bbuchsbaum.github.io/crossform/reference/effect_coupling_result.md),
[`match_control()`](https://bbuchsbaum.github.io/crossform/reference/match_control.md),
[`measurement_components()`](https://bbuchsbaum.github.io/crossform/reference/measurement_components.md),
[`measurement_form()`](https://bbuchsbaum.github.io/crossform/reference/measurement_form.md),
[`pair_lm_query()`](https://bbuchsbaum.github.io/crossform/reference/pair_lm_query.md),
[`pair_query()`](https://bbuchsbaum.github.io/crossform/reference/pair_query.md)

## Examples

``` r
# Three studied items and their matched retrieval probes.
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
matches$value
#>       probe1 probe2 probe3
#> item1      1      0      0
#> item2      0      1      0
#> item3      0      0      1

# Repeated retrievals of one item are multiplicity, not duplicates, so the
# cell count rises rather than being silently collapsed.
match_coupling(
  c("item1", "item1", "item2"), c("probe1", "probe1", "probe2"),
  encoding, retrieval
)$value
#>       probe1 probe2 probe3
#> item1      2      0      0
#> item2      0      1      0
#> item3      0      0      0

# Identifiers must belong to their declared axis.
refused <- try(
  match_coupling("item1", "item2", encoding, retrieval), silent = TRUE
)
conditionMessage(attr(refused, "condition"))
#> [1] "Every `right` identifier must belong to its effect space."
```
