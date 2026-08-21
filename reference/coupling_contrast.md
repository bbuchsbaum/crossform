# Contrast matched and control pair couplings

`coupling_contrast()` subtracts the control coupling from the matched
one to give the matched-versus-control readout as a single fixed pair
query. Normalizing first makes the two sides comparable when they
contain different numbers of cells.

## Usage

``` r
coupling_contrast(matches, controls, normalize = TRUE)
```

## Arguments

- matches, controls:

  Compatible match and control couplings.

- normalize:

  Normalize each coupling to unit total mass before subtraction.

## Value

An axis-bound
[`pair_query()`](https://bbuchsbaum.github.io/crossform/reference/pair_query.md)
whose `$operator` is the contrast, and whose `$metadata$balance` reports
the row and column marginals and whether the operator is
`additive_baseline_invariant`.

## See also

[`match_coupling()`](https://bbuchsbaum.github.io/crossform/reference/match_coupling.md),
[`control_coupling()`](https://bbuchsbaum.github.io/crossform/reference/control_coupling.md),
and
[`match_control()`](https://bbuchsbaum.github.io/crossform/reference/match_control.md),
which instead compiles the same comparison as a regression coefficient
with item nuisance effects.

Other coupling and connectivity views:
[`connectivity()`](https://bbuchsbaum.github.io/crossform/reference/connectivity.md),
[`control_coupling()`](https://bbuchsbaum.github.io/crossform/reference/control_coupling.md),
[`coupling()`](https://bbuchsbaum.github.io/crossform/reference/coupling.md),
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
contrast <- coupling_contrast(matches, control_coupling(matches))
round(as.matrix(contrast$operator), 3)
#>       probe1 probe2 probe3
#> item1  0.333 -0.167 -0.167
#> item2 -0.167  0.333 -0.167
#> item3 -0.167 -0.167  0.333

# Both marginals are zero here, so the readout is unchanged by adding a
# constant to any item or probe effect.
contrast$metadata$balance$additive_baseline_invariant
#> [1] TRUE

# That balance is an observed property of this operator, not a promise
# about anything applied downstream.
contrast$metadata$claim
#> [1] "observed_operator_balance_only"
```
