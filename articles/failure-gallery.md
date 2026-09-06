# Failure gallery

Use this guide when a requested result is unavailable, or when a
familiar transformation would change the scientific question. Each case
explains the problem and a supported next step.

A capability refusal is an R condition of class
`effect_capability_refusal`.
[`catch_refusal()`](https://bbuchsbaum.github.io/crossform/reference/catch_refusal.md)
captures it so you can inspect `$capability`, `$reasons`, and
`$remedies`. Ordinary errors still stop execution; a successful
expression returns `NULL` from
[`catch_refusal()`](https://bbuchsbaum.github.io/crossform/reference/catch_refusal.md).

``` r

library(crossform)
example <- example_fmri_effects()
plan <- plan_geometry(
  example$fit$relation, example$frame,
  cross_partitions(example$fit$relation, independence = "independent")
)
```

The examples below execute against this same fit and plan.

## 1. Correlation-style normalization of a signed cross-generalized form

Dividing crossvalidated similarities by crossvalidated “variances” looks
like `1 - Pearson` correlation distance. It is not: cross-partition
diagonal estimates can be zero or negative, so the normalization is
undefined exactly where crossvalidation is doing its job.

``` r

refusal <- catch_refusal({
  rdm(plan, normalize = "correlation")
})
refusal$capability
#> [1] "guaranteed_psd"
cat(paste(strwrap(refusal$message, width = 76), collapse = "\n"))
#> `rdm()` reports signed squared distances and will not apply
#> correlation-style normalization: crossvalidated diagonal estimates can be
#> zero or negative, so dividing by them is not conventional `1 - Pearson`
#> distance and would silently change the estimand. Conventional correlation
#> distance requires a guaranteed positive-semidefinite self form and its own
#> named view; see the correlation-distance policy vignette.
```

Refusal capability: `guaranteed_psd`. The full boundary is [the
correlation-distance
policy](https://bbuchsbaum.github.io/crossform/articles/correlation-distance-policy.md)
([`vignette("correlation-distance-policy", package = "crossform")`](https://bbuchsbaum.github.io/crossform/articles/correlation-distance-policy.md)
offline).

## 2. “Remove the univariate signal”

Demeaning patterns before a multivariate analysis destroys information
while leaving voxelwise mean effects inside the residual subspace, and
then invites the claim that what remains is “purely multivariate.”

``` r

refusal <- catch_refusal({
  contrast_energy(plan, example$contrast, remove_univariate = TRUE)
})
refusal$capability
#> [1] "nondestructive_decomposition"
cat(paste(strwrap(refusal$message, width = 76), collapse = "\n"))
#> `contrast_energy()` does not remove univariate signal: destructive
#> demeaning would change the estimand while leaving voxelwise mean effects in
#> the residual subspace. The returned view already reports `coherent` (the
#> weighted common spatial mode), `configuration` (its orthogonal remainder),
#> and `total` as an exact additive partition; select the component you mean
#> instead of deleting one.
```

Refusal capability: `nondestructive_decomposition`. The accounting it
points to is exact: `total = coherent + configuration` at every
measurement, for contrasts, RDMs, and RSA coefficients alike.

## 3. Clipping negative crossvalidated distances

Crossvalidated squared distances are unbiased around zero under the
null, so negative estimates are informative, not errors. Truncating them
at zero introduces positive bias.

There is nothing to refuse here because no clipping surface exists:
[`rdm()`](https://bbuchsbaum.github.io/crossform/reference/rdm.md) takes
no truncation argument, its values are documented as signed
(“Cross-generalized distances may be negative”), and
[`geometry_spectrum()`](https://bbuchsbaum.github.io/crossform/reference/geometry_spectrum.md)
reports eigenvalues that “are never truncated at zero.” Any
nonnegativity transformation would be a new named view with its own
estimand, not an option on this one.

## 4. Changing the generalization question

What reproduces across runs within a session is a different scientific
quantity from what reproduces across sessions — even at identical fold
counts. The generalization axis is therefore part of estimand identity,
not an estimator knob:

``` r

run_plan <- plan_geometry(
  example$fit$relation, example$frame,
  cross_partitions(
    example$fit$relation,
    independence = "independent",
    generalizes_over = "run"
  )
)
session_plan <- plan_geometry(
  example$fit$relation, example$frame,
  cross_partitions(
    example$fit$relation,
    independence = "independent",
    generalizes_over = "session"
  )
)
identical(run_plan$scientific_plan_id, session_plan$scientific_plan_id)
#> [1] FALSE
```

These declarations give the plans distinct identities even though they
pair the same inputs and therefore compute the same numbers. Relabeling
this fixture does not create session-level evidence: a real session
analysis also needs partitions and pairings that represent the intended
sessions.

## 5. Independence that was never declared

Different run labels and off-diagonal pairing do not prove statistical
independence. Point geometry remains available, but an analytic
covariance law that consumes endpoint independence refuses if the
declaration is absent:

``` r

refusal <- catch_refusal({
  undeclared_plan <- plan_geometry(
    example$fit$relation, example$frame,
    cross_partitions(example$fit$relation)
  )
  rdm_sampling_covariance(
    undeclared_plan, example$fit, target = "null", at = 1L
  )
})
refusal$capability
#> [1] "sampling_covariance"
cat(paste(strwrap(refusal$message, width = 76), collapse = "\n"))
#> Sampling covariance is unavailable because the pairing does not declare its
#> partition endpoints statistically independent.
```

Refusal capability: `sampling_covariance`. The remedy is to declare
`independence = "independent"` only when the acquisition and observation
model justify it; changing that declaration changes the estimand
identity.

## 6. Uncertainty that was never earned

Analytic error bars need both an admitted metric and a fitted error
channel.

**A learned metric wants an analytic law.** The metric was estimated
from data, so the fixed-metric sampling law does not cover the result:

``` r

refusal <- catch_refusal({
  learned_metric_plan <- plan_geometry(
    example$fit$relation, example$frame,
    cross_partitions(example$fit$relation, independence = "independent"),
    metric = neural_metric(
      diag(example$domain$n_features), example$domain,
      estimation = "learned_frozen",
      provenance = list(
        frozen = TRUE,
        training_signature = paste0("sha256:", strrep("ab", 32))
      )
    )
  )
  rdm_sampling_covariance(
    learned_metric_plan, example$fit, target = "null", at = 1L
  )
})
refusal$capability
#> [1] "fixed_metric_sampling_law"
cat(paste(strwrap(refusal$message, width = 76), collapse = "\n"))
#> Analytic RDM sampling covariance is unavailable because this plan's neural
#> metric was learned. Version 0.1 does not propagate metric-estimation
#> uncertainty; use a geometry plan with one common fixed metric or retain
#> this result as a signed point estimate.
```

Refusal capability: `fixed_metric_sampling_law`.

**Precomputed betas want residual uncertainty back.** Beta matrices
alone cannot recover the error channel that was discarded when they were
exported — and the refusal reports every unmet requirement, not only the
first:

``` r

refusal <- catch_refusal({
  rdm_sampling_covariance(plan, example$fit$relation, target = "null", at = 1L)
})
refusal$capability
#> [1] "sampling_covariance"
cat(paste(strwrap(refusal$message, width = 76), collapse = "\n"))
#> Sampling covariance is unavailable; 2 requirements of the admitted analytic
#> law are unmet: * this evidence plan has only a precomputed relation and no
#> error channel. Refit raw observations with `lm_relation_fit()` or supply a
#> validated, identity-bound external error channel; beta matrices alone
#> cannot recover residual uncertainty. * no single sampling axis is declared,
#> or the declared axis conflicts with the error channel's recorded sampling
#> unit.
```

Both reasons are separately actionable. Requirements that only restate
the missing channel — such as whether the partitions share one error
structure — cannot be evaluated without a channel to inspect, so they
are not listed as extra failures.

Refusal capability: `sampling_covariance`. Check availability before
requesting a covariance with `sampling_capabilities(plan, example$fit)`.
It reports the reasons and remedies without computing the covariance.

## Choose a supported next step

Retain signed distances for cross-run geometry. Use `coherent`,
`configuration`, or `total` to select the contrast component you need.
For uncertainty, inspect
[`sampling_capabilities()`](https://bbuchsbaum.github.io/crossform/reference/sampling_capabilities.md)
and address the reported requirements; an independence declaration must
be justified by the acquisition and model. [Fit condition effects from
observations](https://bbuchsbaum.github.io/crossform/articles/from-observations.md)
shows how to retain residuals, and [Reading
results](https://bbuchsbaum.github.io/crossform/articles/interpreting-results.md)
explains the available uncertainty targets.

## Further checks

Every case is covered by executable tests:
[`tests/testthat/test-capability-refusals.R`](https://github.com/bbuchsbaum/crossform/blob/main/tests/testthat/test-capability-refusals.R),
[`tests/testthat/test-integrity-guards.R`](https://github.com/bbuchsbaum/crossform/blob/main/tests/testthat/test-integrity-guards.R),
and
[`tests/testthat/test-generalization-axis.R`](https://github.com/bbuchsbaum/crossform/blob/main/tests/testthat/test-generalization-axis.R).
The Haxby 2001 exemplar records three of these refusals firing against
real data in
[`exemplars/haxby2001/results/smoke-report.md`](https://github.com/bbuchsbaum/crossform/blob/main/exemplars/haxby2001/results/smoke-report.md).
