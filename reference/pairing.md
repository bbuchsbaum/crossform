# Construct a partition pairing

`pairing()` declares exactly which partition products a plan may form
and what may be claimed about them. Use it when you need an explicit or
directed edge set; use
[`cross_partitions()`](https://bbuchsbaum.github.io/crossform/reference/cross_partitions.md)
for the usual all-distinct-pairs case.

## Usage

``` r
pairing(
  left,
  right,
  weight = NULL,
  directed = FALSE,
  self_pairs = c("forbid", "allow_biased"),
  independence = NULL,
  generalizes_over = NULL
)
```

## Arguments

- left, right:

  Equal-length vectors naming partition endpoints.

- weight:

  Optional finite nonnegative edge weights. They are normalized to sum
  to one.

- directed:

  Whether endpoint roles have scientific meaning.

- self_pairs:

  Whether diagonal self-products are forbidden or explicitly admitted as
  noise-biased estimates.

- independence:

  Whether distinct endpoint estimates are declared independent. It must
  be stated explicitly for an independence-based interpretation. `NULL`
  records `"undeclared"`; self-products must be marked
  `"not_independent"`.

- generalizes_over:

  Optional name of the sampling axis this pairing generalizes across,
  such as `"run"`, `"session"`, or `"task"`. The axis is part of the
  scientific estimand: cross-run and cross-session reproduction are
  different quantities even at identical fold counts, so the declared
  axis is bound into every plan identity built from this pairing.
  Leaving it `NULL` records the axis as undeclared.

## Value

An `effect_pairing` data frame with `left`, `right`, and unit-mass
`weight` columns, plus the `directed`, `self_pairs`, `independence`,
`generalizes_over`, and derived `estimate` attributes that name the
estimand.

## See also

[`cross_partitions()`](https://bbuchsbaum.github.io/crossform/reference/cross_partitions.md)
for all distinct pairs, and
[`plan_geometry()`](https://bbuchsbaum.github.io/crossform/reference/plan_geometry.md),
which binds the pairing into plan identity.

Other generalization pairings:
[`cross_partitions()`](https://bbuchsbaum.github.io/crossform/reference/cross_partitions.md),
[`reduce_partitions()`](https://bbuchsbaum.github.io/crossform/reference/reduce_partitions.md)

## Examples

``` r
# Directed cross-session edges, declared independent and named as
# generalizing across sessions.
over <- pairing(
  c("ses1", "ses1"), c("ses2", "ses3"),
  directed = TRUE, independence = "independent",
  generalizes_over = "session"
)
over
#> <effect_pairing>
#>   pairs:        2
#>   left:         ses1
#>   right:        ses2, ses3
#>   weights:      equal (0.5)
#>   directed:     yes
#>   self pairs:   forbid
#>   independence: independent
#>   generalizes:  session
#>   estimate:     cross_generalized
attr(over, "estimate")
#> [1] "cross_generalized"

# Weights are normalized to unit mass, so they are estimator weights, not
# counts of independent replicates.
sum(over$weight)
#> [1] 1

# Self-products are biased and must say so; the default forbids them.
refused <- try(pairing("run1", "run1"), silent = TRUE)
conditionMessage(attr(refused, "condition"))
#> [1] "Edge 1 pairs `run1` with itself, which is a noise-biased estimate rather than a cross-generalized one. Declare it with `self_pairs = \"allow_biased\"` if that is what you mean."
attr(
  pairing("run1", "run1", self_pairs = "allow_biased",
    independence = "not_independent"),
  "estimate"
)
#> [1] "self_product_biased"
```
