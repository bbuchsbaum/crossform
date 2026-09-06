# Null-imposed subject-level wild bootstrap for population coefficients

`population_wild_bootstrap()` tests one explicit linear contrast of the
population-model coefficients against one explicit null at every
admitted node and query. It uses a null-imposed, HC3-studentized wild
bootstrap. One random weight is drawn per planned subject and replicate
and reused across **all** nodes and queries, preserving the
within-subject dependence of the jointly analysed result.

## Usage

``` r
population_wild_bootstrap(
  x,
  contrast,
  null = 0,
  replicates = 999L,
  seed = NULL,
  weights = c("rademacher", "mammen"),
  level = 0.95,
  leverage_tolerance = 1e-08
)
```

## Arguments

- x:

  An `effect_population_result` from
  [`estimate_population()`](https://bbuchsbaum.github.io/crossform/reference/estimate_population.md).

- contrast:

  One population coefficient name, or a numeric vector with one weight
  per model term. Named weights are aligned by term name.

- null:

  Finite scalar value of the tested coefficient contrast.

- replicates:

  Number of bootstrap replicates, at least 99.

- seed:

  Required nonnegative integer random seed.

- weights:

  Wild-weight distribution: `"rademacher"` or `"mammen"`.

- level:

  Nominal two-sided test level used for the absolute bootstrap-t
  critical value and `$reject` diagnostic.

- leverage_tolerance:

  Positive tolerance below which \\1-h_i\\ is treated as numerically
  zero.

## Value

An `effect_population_wild_bootstrap` carrying the explicit null,
observed and replicate HC3-t statistics, shared subject weights,
replicate-level failures, plus-one p-values, critical values, rejection
flags, Monte Carlo SEs, exact coverage, conditioning and provenance.

## Algorithm and null

For coefficient contrast \\c\\ and null \\c'\beta = b_0\\, each cell is
first fit under the linear restriction. Bootstrap response \\y^\* =
X\hat\beta_R + w_i e\_{Ri}/(1-h_i)\\ uses the restricted residual and
its HC3 leverage adjustment. Each replicate is refit without the
restriction and studentized by its own HC3 covariance. The two-sided
p-value uses the plus-one rule, so it is never zero; `$monte_carlo_se`
reports its finite-replication uncertainty.

`"rademacher"` weights are symmetric \\-1,+1\\ draws with equal
probability. `"mammen"` uses the two-point, mean-zero, variance-one
Mammen distribution, whose third moment also equals one; it is provided
for asymmetric-error sensitivity analyses, not as an automatic
improvement.

## Coverage, failures, and conditioning

The exact subject set from `$coverage` is used separately at every cell,
while the weights come from one shared planned-subject matrix. Cells
with fewer than six available subjects, rank deficiency, no residual df,
leverage within `leverage_tolerance` of one, or a nonpositive observed
HC3 standard error are refused with a status and reason. Replicate
failures are retained in `$replicate_failure_reason`; a cell is reported
only when at least 90 percent and at least 20 replicates succeed. All
inference is conditional on the realized transport and coverage. It does
not propagate uncertainty from estimating either. Cross-fitting can
limit circular reuse of responses, but does not make this bootstrap
marginal over transport estimation or fold assignment.

The returned `$weights`, `$seed`, `$rng`, `$weight_signature`, replicate
statistics, success flags, failure reasons, contrast, null and parent
receipt provide full resampling provenance. The caller's random-number
state is restored on exit.

## See also

[`population_uncertainty()`](https://bbuchsbaum.github.io/crossform/reference/population_uncertainty.md)
for analytic classical and HC3 covariance without resampling.

Other population transports:
[`anatomical_transport()`](https://bbuchsbaum.github.io/crossform/reference/anatomical_transport.md),
[`estimate_population()`](https://bbuchsbaum.github.io/crossform/reference/estimate_population.md),
[`external_transport()`](https://bbuchsbaum.github.io/crossform/reference/external_transport.md),
[`heterogeneity()`](https://bbuchsbaum.github.io/crossform/reference/heterogeneity.md),
[`location_transport()`](https://bbuchsbaum.github.io/crossform/reference/location_transport.md),
[`materialize_population()`](https://bbuchsbaum.github.io/crossform/reference/materialize_population.md),
[`plan_population()`](https://bbuchsbaum.github.io/crossform/reference/plan_population.md),
[`population_prevalence()`](https://bbuchsbaum.github.io/crossform/reference/population_prevalence.md),
[`population_uncertainty()`](https://bbuchsbaum.github.io/crossform/reference/population_uncertainty.md),
[`population_views`](https://bbuchsbaum.github.io/crossform/reference/population_views.md),
[`transport_values()`](https://bbuchsbaum.github.io/crossform/reference/transport_values.md)
