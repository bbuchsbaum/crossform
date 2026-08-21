# Pair every distinct partition once

Distinct partition estimates are declared independent endpoints for
unbiased cross-products. The resulting pair rows are not independent
sampling replicates: pairs that share a partition also share its
estimation error. In particular,
`sd(pair_values) / sqrt(number_of_pairs)` is not a valid standard error
for their all-pairs mean. Under the fixed-metric, equal-partition
separable model, use
[`rdm_sampling_covariance()`](https://bbuchsbaum.github.io/crossform/reference/rdm_sampling_covariance.md)
for the admitted analytic RDM covariance law.

## Usage

``` r
cross_partitions(partitions, independence = NULL, generalizes_over = NULL)
```

## Arguments

- partitions:

  Partition identifiers or an `effect_relation`. Pass `fit$relation`
  when starting from
  [`lm_relation_fit()`](https://bbuchsbaum.github.io/crossform/reference/lm_relation_fit.md).

- independence:

  Explicit endpoint-independence declaration. Leaving it `NULL`
  preserves a point estimand but does not earn cross-generalized or
  analytic sampling-law capabilities.

- generalizes_over:

  Optional name of the sampling axis the partitions represent, such as
  `"run"` or `"session"`; see
  [`pairing()`](https://bbuchsbaum.github.io/crossform/reference/pairing.md).
  The declared axis is bound into every plan identity built from this
  pairing.

## Value

An undirected pairing containing one row per unordered pair, with
normalized estimator weights. Its rows are computational contributions,
not a declaration of edge-level sampling independence.

## References

Diedrichsen J, Provost S, Zareamoghaddam H (2016), "On the distribution
of cross-validated Mahalanobis distances", especially Eqs. 10, 13, and
35.
[doi:10.48550/arXiv.1607.01371](https://doi.org/10.48550/arXiv.1607.01371)

## See also

[`pairing()`](https://bbuchsbaum.github.io/crossform/reference/pairing.md)
for explicit or directed edge sets,
[`plan_geometry()`](https://bbuchsbaum.github.io/crossform/reference/plan_geometry.md)
which consumes this pairing, and
[`rdm_sampling_covariance()`](https://bbuchsbaum.github.io/crossform/reference/rdm_sampling_covariance.md)
for the admitted uncertainty law.

Other generalization pairings:
[`pairing()`](https://bbuchsbaum.github.io/crossform/reference/pairing.md),
[`reduce_partitions()`](https://bbuchsbaum.github.io/crossform/reference/reduce_partitions.md)

## Examples

``` r
# Four runs give six unordered pairs, each weighted to unit total mass.
example <- example_fmri_effects()
over <- cross_partitions(
  example$fit$relation,
  independence = "independent", generalizes_over = "run"
)
nrow(over)
#> [1] 6
sum(over$weight)
#> [1] 1

# The declared axis and independence become part of the estimand, so they
# travel into every plan built from this pairing.
c(estimate = attr(over, "estimate"),
  axis = attr(over, "generalizes_over"))
#>            estimate                axis 
#> "cross_generalized"               "run" 

# Leaving independence undeclared still yields a point estimand, but it
# does not earn cross-generalized capabilities.
attr(cross_partitions(example$fit$relation), "estimate")
#> [1] "independence_undeclared"
```
