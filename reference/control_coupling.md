# Mark eligible control pairs

`control_coupling()` names the comparison set for a matched-pair
analysis: the eligible cells that are not matches. Pair it with the
matches through
[`coupling_contrast()`](https://bbuchsbaum.github.io/crossform/reference/coupling_contrast.md)
to obtain the matched-versus-control query.

## Usage

``` r
control_coupling(matches, include_matches = FALSE)
```

## Arguments

- matches:

  A
  [`match_coupling()`](https://bbuchsbaum.github.io/crossform/reference/match_coupling.md)
  value.

- include_matches:

  Whether matched cells are also controls. The default excludes them.

## Value

An `effect_pair_coupling` with `kind = "control"`, a 0/1 `$value`
indicator over the same axes and `$eligible` set as `matches`.

## See also

[`match_coupling()`](https://bbuchsbaum.github.io/crossform/reference/match_coupling.md)
and
[`coupling_contrast()`](https://bbuchsbaum.github.io/crossform/reference/coupling_contrast.md).

Other coupling and connectivity views:
[`connectivity()`](https://bbuchsbaum.github.io/crossform/reference/connectivity.md),
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

# The default controls are the six mismatched cells.
controls <- control_coupling(matches)
controls$value
#>       probe1 probe2 probe3
#> item1      0      1      1
#> item2      1      0      1
#> item3      1      1      0

# Including the matches gives every eligible cell, which is the marginal
# baseline rather than a contrast partner.
control_coupling(matches, include_matches = TRUE)$value
#>       probe1 probe2 probe3
#> item1      1      1      1
#> item2      1      1      1
#> item3      1      1      1

# With no eligible non-matched cell left, there is nothing to contrast.
saturated <- match_coupling(
  rep(c("item1", "item2", "item3"), each = 3),
  rep(c("probe1", "probe2", "probe3"), 3), encoding, retrieval
)
refused <- try(control_coupling(saturated), silent = TRUE)
conditionMessage(attr(refused, "condition"))
#> [1] "No eligible control pairs remain."
```
