# Read the group uncertainty layers of an estimated population form

`population_uncertainty()` reports two error bars on an
[`estimate_population()`](https://bbuchsbaum.github.io/crossform/reference/estimate_population.md)
result and keeps them apart. The **between-subject** layer is the
scatter of the participants about the group fit: a standard error per
group node, query and model term from the group OLS residuals, with that
fit's residual degrees of freedom. The **within-subject** layer is the
sampling variance of one participant's transported value, and it exists
only where it is exact.

## Usage

``` r
population_uncertainty(
  x,
  term = NULL,
  level = 0.95,
  estimator = c("classical", "HC3"),
  leverage_tolerance = 1e-08
)
```

## Arguments

- x:

  An `effect_population_result` from
  [`estimate_population()`](https://bbuchsbaum.github.io/crossform/reference/estimate_population.md).

- term:

  Which group model columns to report, by name or position. `NULL` (the
  default) reports every column.

- level:

  Nominal two-sided interval level for `$between$lower` and
  `$between$upper`. The default `0.95` is the level the recorded
  coverage simulation measured; the interval is labelled uncalibrated
  whatever level is asked for.

- estimator:

  Between-subject covariance estimator. `"classical"` is ordinary
  homoskedastic OLS; `"HC3"` is the leverage-adjusted
  heteroskedasticity-robust sandwich estimator.

- leverage_tolerance:

  Positive tolerance below which HC3 treats \\1-h_i\\ as numerically
  zero and refuses that cell.

## Value

An `effect_population_uncertainty`: `$between` holding `$estimate`,
`$covariance`, `$se`, `$t`, `$lower`, `$upper`, leverage, adjusted
residuals and local design diagnostics, estimator identity, assumptions,
status and refusal reasons; `$within` holding the transported layer's
`$admitted`, `$coefficient`, `$source_node`, `$variance` and `$refusal`,
or `NULL` when
[`estimate_population()`](https://bbuchsbaum.github.io/crossform/reference/estimate_population.md)
was not given `uncertainty`; `$separation` stating that the two are
never pooled; and the `$index`, `$queries`, `$ledger`, `$semantics`,
`$normalization` and `$receipt` of the result it read.
`as.data.frame(x, layer = )` returns one layer in long form.

## Classical OLS and HC3

`estimator = "classical"` reproduces the ordinary OLS covariance
\\s^2(X'X)^{-1}\\. `estimator = "HC3"` returns the
heteroskedasticity-robust sandwich covariance
\\(X'X)^{-1}X'\mathrm{diag}\\e_i^2/(1-h_i)^2\\X(X'X)^{-1}\\. The latter
is the recommended sensitivity analysis when subject-level variance
differs with covariates or transport quality. It is not a distributional
calibration theorem, so both routes retain the `"uncalibrated"` label.

Every node-query cell is computed from its exact contributing subject
set. `$between` therefore carries the full coefficient covariance, SEs,
local `n`, design rank and residual df, subject leverage, the residual
adjustment actually used in the sandwich, maximum leverage, conditioning
statement, and a cellwise status/reason. Rank-deficient and saturated
cells are refused. HC3 also refuses a cell when any \\1-h_i\\ is no
larger than `leverage_tolerance`; dividing by an almost-zero leverage
complement would turn a numerical singularity into a confident-looking
result.

## The two layers are never pooled

They answer different questions and are not summands. The
between-subject residual already contains whatever measurement error
survived into each participant's transported value, so adding the
within-subject variance to it would double-count the shared part and
still miss the covariance a variance-components model would need. The
record carries `$between` and `$within` as separate blocks,
[`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html) emits one
`layer` at a time, and there is no field holding their sum.

## The `t` is uncalibrated, and stays uncalibrated

`$between$t` is the estimate over its standard error. Its null
distribution for `estimator = "classical"` has been *measured* against
\\t\_{df}\\ by `benchmarks/run-population-null-coverage.R`, a
2,000-replication null simulation of the group layer. Under a
**correctly specified** group model the nominal 95% interval covered the
null term in **0.9485** of replications at `N = 6`, **0.9500** at
`N = 8`, **0.9520** at `N = 12` and **0.9480** at `N = 24` (Monte Carlo
standard error 0.005), and the Kolmogorov-Smirnov distance between the
null `t` and \\t\_{df}\\ was at most **0.0209** (`p >= 0.34`). The
arithmetic is right, which was never the part in doubt.

The same simulation's second arm is why the label does not move. When
each participant's transported value carries a variance that depends on
the group covariates — the transport-heterogeneity analogue, since a
participant whose warp is poor is noisier and warp quality is not
independent of age, motion or head size — coverage of the same interval
falls to **0.9230** at `N = 6` and **0.8850** at `N = 24`. It gets
**worse** with more participants, because the bias is in the standard
error and not in the sample size: at `N = 24` the nominal 5% test
rejects a true null **11.5%** of the time, and the null `t` is
distinguishable from \\t\_{df}\\ at `p = 1.2e-6`.

Those coverage figures do not transfer to HC3. Its randomized court
proves agreement with the defining sandwich covariance over full-rank
designs and adversarial leverage fixtures; it does not supply a
finite-sample reference distribution. HC3 therefore retains the same
`"uncalibrated"` label. The expanded paired benchmark in
`inst/extdata/certification/population-calibration-results.csv` records
classical, HC3, and wild-bootstrap behavior across eight synthetic
regimes. It supports only the regimes and conditional targets named
there: informative coverage remains ineligible for a marginal claim, and
no arm is marginal over transport estimation.

A real population fit carries misspecification of unknown degree and a
transport whose displacement, entropy and subject coverage vary across
participants (`population-form-v1` section 7.5), which is exactly the
second arm's regime. Report the statistic; do not report a p-value
derived from it without an argument that those diagnostics are benign.

## What the within layer is admitted for

Participant `i`'s transported value at group column `u` is a fixed
linear functional \\w'z_i\\ of that participant's native query values,
so its variance needs the covariance *between* native nodes. D8 refuses
that object (capability `cross_node_sampling_covariance`), so the layer
is admitted exactly where `w` has one nonzero entry: there the cross
terms carry weight zero and \\\mathrm{Var} =
w_x^2\\\mathrm{Var}(z\_{ix})\\ is exact. **No independence assumption is
made anywhere**, and none would be defensible: overlapping searchlight
supports under spatially correlated noise are positively correlated, so
a diagonal sum would be an under-estimate in a known direction.

A future joint covariance may support transported precision or a
separately calibrated spatial procedure, but it is not required by the
cellwise between-subject covariance computed here. Dense/sparse
representation, PSD validation and scaling requirements are specified in
`design/cross-node-covariance-contract.md`. maxT and multiple-comparison
methods are optional later consumers, not part of this gate.

Admission is per participant and per group column, and
`$within$admitted` is the matrix that records it. A hard anatomical
parcellation usually admits nothing, because every group node collects
several native nodes; that is the refusal being visible rather than the
layer being broken.

Two gates come before the per-column one, both recorded in
`$within$refusal$reasons` rather than raised:

- `same_data_ratio_normalization` — a `"unit_budget"` population divides
  each ledger by a total read from the same data, and section 4.3
  records that the standard error of that divisor does not exist.

- `native_node_labels_unaligned` — the covariance batch is named by the
  nodes `rdm_sampling_covariance(at = )` read and the transport by its
  own `native_index`. Binding two differently named node sets by
  position would attach one node's error bar to another's coefficient.

## Refusals

Each is an `effect_capability_refusal` in namespace
`"population_uncertainty"` (see
[`catch_refusal()`](https://bbuchsbaum.github.io/crossform/reference/catch_refusal.md)).

- `population_between_subject_residuals` — a `"complete_form"` result.
  The streamed route does not retain the participant residuals, by
  design.

- `between_subject_residual_df` — a saturated group model. With no
  residual degrees of freedom the participants' scatter is not
  estimable.

## Transport conditioning

Classical and HC3 intervals use the exact subject-level residuals and
may account for homoskedastic or heteroskedastic between-subject
variation. They remain conditional on each realized transport. In
particular, HC3 does not include uncertainty from learning a functional
transport or from assigning its cross-fitting folds. The complete
conditioning record is retained in `$between$conditioning`.

## References

`design/population-form-contract.md` (`population-form-v1`), section 7;
`benchmarks/POPULATION-NULL-COVERAGE.md` for the recorded null
simulation these numbers come from.

## See also

[`estimate_population()`](https://bbuchsbaum.github.io/crossform/reference/estimate_population.md)
for the run this reads, and
[population_views](https://bbuchsbaum.github.io/crossform/reference/population_views.md)
for the group point estimates the standard errors belong to.

Other population transports:
[`anatomical_transport()`](https://bbuchsbaum.github.io/crossform/reference/anatomical_transport.md),
[`estimate_population()`](https://bbuchsbaum.github.io/crossform/reference/estimate_population.md),
[`external_transport()`](https://bbuchsbaum.github.io/crossform/reference/external_transport.md),
[`heterogeneity()`](https://bbuchsbaum.github.io/crossform/reference/heterogeneity.md),
[`location_transport()`](https://bbuchsbaum.github.io/crossform/reference/location_transport.md),
[`materialize_population()`](https://bbuchsbaum.github.io/crossform/reference/materialize_population.md),
[`plan_population()`](https://bbuchsbaum.github.io/crossform/reference/plan_population.md),
[`population_prevalence()`](https://bbuchsbaum.github.io/crossform/reference/population_prevalence.md),
[`population_views`](https://bbuchsbaum.github.io/crossform/reference/population_views.md),
[`population_wild_bootstrap()`](https://bbuchsbaum.github.io/crossform/reference/population_wild_bootstrap.md),
[`transport_values()`](https://bbuchsbaum.github.io/crossform/reference/transport_values.md)

## Examples

``` r
# Six participants on different native frames, one covariate at the group
# level, and a bank of two contrasts.
effects <- effect_space(c("face", "house"), basis_id = "popunc:v1")
subject <- function(id, n, gain) {
  domain <- abstract_domain(n, coordinates = cbind(x = seq_len(n) - 1),
    feature_ids = paste0("f", seq_len(n)), id = id)
  values <- function(divisor) matrix(
    gain * seq_len(2 * n) / (n * divisor), 2, n,
    dimnames = list(c("face", "house"), NULL)
  )
  rel <- relation(list(run1 = values(1), run2 = values(1.7)),
    effects = effects, domain = domain)
  plan_geometry(rel, compile_frame(voxelwise(), domain),
    cross_partitions(rel))
}
carrier <- function(n) anatomical_transport(
  native_coords = cbind(seq_len(n) - 1), group_coords = cbind(c(0, 4)),
  semantics = "budget"
)
sizes <- c(s01 = 5L, s02 = 6L, s03 = 7L, s04 = 8L, s05 = 9L, s06 = 10L)
gains <- c(s01 = 1, s02 = 1.4, s03 = 0.7, s04 = 1.1, s05 = 0.9, s06 = 1.3)
subjects <- stats::setNames(lapply(names(sizes), function(id)
  subject(id, sizes[[id]], gains[[id]])), names(sizes))
covariates <- data.frame(age = c(21, 34, 27, 45, 31, 38))
rownames(covariates) <- names(sizes)
plan <- plan_population(subjects, lapply(sizes, carrier),
  model = ~ age, data = covariates)
fit <- estimate_population(plan, rbind(`face-house` = c(1, -1)))

error_bars <- population_uncertainty(fit)
error_bars
#> <effect_population_uncertainty>
#>   ledger:          transported_total (component "total")
#>   group nodes:     2 + sink
#>   queries:         1 (face-house)
#>   terms:           (Intercept), age
#>   between-subject: classical_ols SE 0 to 0.06273, cell df 4 to 4, |t| up ...
#>   leverage:        max h 0.5972; estimated, degenerate_residual_scale
#>   assumptions:     independent_subjects, fixed_group_design, conditional_...
#>   interval:        95% nominal, uncalibrated
#>   within-subject:  absent (estimate_population() was given no `uncertainty`)
#>   inference:       conditional_on_realized_transport; uncertainty not pro...
#>   normalization:   none
#>   estimand:        population-sha256:e1ed8d16b2bc...
#>   The two layers are reported separately and are never pooled: a
#>     within-subject sampling variance and a between-subject residual
#>     variance answer different questions, and their sum answers neither.
#>   These intervals condition on the realized transport. Cross-fitting limits
#>     circularity but does not propagate uncertainty from estimating that
#>     transport.
#>   The t is UNCALIBRATED for real data. Measured against t_df, the nominal
#>     95% interval covers 0.948 to 0.952 of the time under a correctly
#>     specified group model, and 0.885 to 0.923 when the participants' noise
#>     is linked to the group covariates -- a nominal 5% test rejecting a true
#>     null 11.5% of the time at N = 24. See
#>     benchmarks/POPULATION-NULL-COVERAGE.md.
#>   next:          as.data.frame(x, layer = "between"), x$between$se[, , term]

# The standard error of the group mean at each node, and its uncalibrated t.
error_bars$between$se[, , "(Intercept)"]
#>     group1     group2     <sink> 
#> 0.06272822 0.04709086 0.00000000 
error_bars$between$t[, , "age"]
#>     group1     group2     <sink> 
#> -0.3849895  0.8951872         NA 

# One layer at a time. There is no combined table, because the two layers
# are not summands.
head(as.data.frame(error_bars), 4)
#>     node coord1  sink  units           layer            ledger      query
#> 1 group1      0 FALSE budget between_subject transported_total face-house
#> 2 group2      4 FALSE budget between_subject transported_total face-house
#> 3 <sink>     NA  TRUE budget between_subject transported_total face-house
#> 4 group1      0 FALSE budget between_subject transported_total face-house
#>          term     estimator
#> 1 (Intercept) classical_ols
#> 2 (Intercept) classical_ols
#> 3 (Intercept) classical_ols
#> 4         age classical_ols
#>                                                                                                                                                                    assumptions
#> 1 independent_subjects;fixed_group_design;conditional_on_realized_transport;conditional_on_realized_coverage;homoskedastic_subject_errors;correctly_specified_conditional_mean
#> 2 independent_subjects;fixed_group_design;conditional_on_realized_transport;conditional_on_realized_coverage;homoskedastic_subject_errors;correctly_specified_conditional_mean
#> 3 independent_subjects;fixed_group_design;conditional_on_realized_transport;conditional_on_realized_coverage;homoskedastic_subject_errors;correctly_specified_conditional_mean
#> 4 independent_subjects;fixed_group_design;conditional_on_realized_transport;conditional_on_realized_coverage;homoskedastic_subject_errors;correctly_specified_conditional_mean
#>        estimate          se n design_rank residual_df max_leverage
#> 1  0.0677011006 0.062728222 6           2           4    0.5971698
#> 2  0.0134879486 0.047090859 6           2           4    0.5971698
#> 3  0.0000000000 0.000000000 6           2           4    0.5971698
#> 4 -0.0007196853 0.001869364 6           2           4    0.5971698
#>   max_abs_adjusted_residual        uncertainty_status
#> 1                0.05284663                 estimated
#> 2                0.03987747                 estimated
#> 3                0.00000000 degenerate_residual_scale
#> 4                0.05284663                 estimated
#>           uncertainty_reason            transport_conditioning
#> 1                       <NA> conditional_on_realized_transport
#> 2                       <NA> conditional_on_realized_transport
#> 3 nonpositive_standard_error conditional_on_realized_transport
#> 4                       <NA> conditional_on_realized_transport
#>              coverage_conditioning          t level        lower     upper
#> 1 conditional_on_realized_coverage  1.0792766  0.95 -0.106460365 0.2418626
#> 2 conditional_on_realized_coverage  0.2864239  0.95 -0.117257236 0.1442331
#> 3 conditional_on_realized_coverage         NA  0.95           NA        NA
#> 4 conditional_on_realized_coverage -0.3849895  0.95 -0.005909871 0.0044705
#>    calibration
#> 1 uncalibrated
#> 2 uncalibrated
#> 3 uncalibrated
#> 4 uncalibrated
```
