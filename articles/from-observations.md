# Fit condition effects from observations

Use this guide when you have one scan-by-feature response matrix per
run, event and confound tables, and a design matrix or design compiler.
The goal is to fit named condition effects in every run and pass them to
[`contrast_energy()`](https://bbuchsbaum.github.io/crossform/reference/contrast_energy.md),
[`rdm()`](https://bbuchsbaum.github.io/crossform/reference/rdm.md), or
[`rsa()`](https://bbuchsbaum.github.io/crossform/reference/rsa.md). The
returned `effect_relation_fit` can also retain residuals for
within-participant uncertainty calculations.

The workflow has five steps: bind observations to their clocks and
metadata, declare the condition effects, declare the observation model,
call
[`plan_relation()`](https://bbuchsbaum.github.io/crossform/reference/plan_relation.md),
and fit it with
[`estimate_relation()`](https://bbuchsbaum.github.io/crossform/reference/estimate_relation.md).
A separate
[`plan_geometry()`](https://bbuchsbaum.github.io/crossform/reference/plan_geometry.md)
then states what those fitted effects must reproduce across.

If you already hold a compiled design matrix, or only per-run beta
matrices, you do not need all five steps. Those two shorter routes are
worked at the end of this guide; the only thing that changes along the
way is whether the residual channel — and with it the analytic standard
error — survives.

## What data does the fit require?

Three objects, and nothing else:

1.  **one scan-by-feature response matrix per run**, in a named list —
    scans in rows, neural features in columns;
2.  **an events table** — one row per event, with `partition`, `onset`,
    `duration`, and the condition label, on the same clock as the scan
    index;
3.  **a confounds table** — one row per scan, joined to the responses by
    `partition` and `observation_id`, carrying nuisance regressors and
    an optional Boolean censor column.

The example below generates all three so the vignette can check its own
answers. It has three runs, three named conditions, twelve neural
features, and two censored scans per run. Read the code as a description
of the shape your own data must have; substitute your matrices and
tables for the generated ones.

The design and the planted signal come first. Each run’s responses are
the condition means plus a linear drift, a motion covariate, and noise.

``` r

set.seed(20260815)
partitions <- paste0("run-", 1:3)
n_scan <- 18L
conditions <- c("face", "body", "tool")
domain <- abstract_domain(12L, id = "from-observations:v1")

condition <- factor(
  rep(conditions, length.out = n_scan),
  levels = conditions
)
indicator <- stats::model.matrix(~ condition - 1)
colnames(indicator) <- conditions
motion <- sin(seq(0, 2 * pi, length.out = n_scan))
drift <- seq(-1, 1, length.out = n_scan)
design <- cbind(indicator, drift = drift, motion = motion)
rownames(design) <- as.character(seq_len(n_scan))

truth <- matrix(
  rnorm(length(conditions) * domain$n_features, sd = 0.35),
  nrow = length(conditions),
  dimnames = list(conditions, domain$feature_ids)
)
nuisance_truth <- matrix(
  rnorm(2L * domain$n_features, sd = 0.08),
  nrow = 2L
)

# (1) one scan-by-feature response matrix per run
response_runs <- stats::setNames(lapply(partitions, function(partition) {
  design[, conditions, drop = FALSE] %*% truth +
    design[, c("drift", "motion"), drop = FALSE] %*% nuisance_truth +
    matrix(rnorm(n_scan * domain$n_features, sd = 0.4),
      nrow = n_scan)
}), partitions)
```

Each response matrix needs a clock.
[`observation_index()`](https://bbuchsbaum.github.io/crossform/reference/observation_index.md)
records the scan identifiers, their partition, and their acquisition
times; the events table then speaks the same units.

``` r

scan_indexes <- stats::setNames(lapply(partitions, function(partition) {
  observation_index(
    seq_len(n_scan), partition,
    time = seq(0, by = 2, length.out = n_scan), units = "seconds"
  )
}), partitions)

# (2) one events table for the whole study
event_table <- do.call(rbind, lapply(partitions, function(partition) {
  data.frame(
    partition = partition,
    event_id = paste0(partition, "-", seq_len(n_scan)),
    onset = seq(0, by = 2, length.out = n_scan),
    duration = 0,
    condition = as.character(condition)
  )
}))

# (3) one confounds table, one row per scan, with a censor column
confound_table <- do.call(rbind, lapply(partitions, function(partition) {
  data.frame(
    partition = partition,
    observation_id = seq_len(n_scan),
    motion = motion,
    retained = seq_len(n_scan) <= 16L
  )
}))
```

That is the whole input. Its shape:

``` r

data.frame(
  partition = names(response_runs),
  scans = vapply(response_runs, nrow, integer(1)),
  features = vapply(response_runs, ncol, integer(1))
)
#>       partition scans features
#> run-1     run-1    18       12
#> run-2     run-2    18       12
#> run-3     run-3    18       12
utils::head(event_table, 3)
#>   partition event_id onset duration condition
#> 1     run-1  run-1-1     0        0      face
#> 2     run-1  run-1-2     2        0      body
#> 3     run-1  run-1-3     4        0      tool
utils::head(confound_table, 3)
#>   partition observation_id    motion retained
#> 1     run-1              1 0.0000000     TRUE
#> 2     run-1              2 0.3612417     TRUE
#> 3     run-1              3 0.6736956     TRUE
```

## Bind observations to events and confounds

The constructors below turn the matrices and tables into records with
named axes.
[`study()`](https://bbuchsbaum.github.io/crossform/reference/study.md)
then verifies their correspondence. It rejects missing or extra confound
rows, incompatible identifier types, mismatched clock units, and events
outside the observation window.

``` r

neural_record <- observations(response_runs, scan_indexes, domain)
event_record <- observation_events(event_table)
confound_record <- observation_confounds(
  confound_table, censor = "retained"
)
hierarchy <- partition_hierarchy(data.frame(
  partition = partitions, run = partitions,
  session = "session-1", subject = "subject-1"
))

study_facts <- study(
  neural_record, event_record, confound_record, hierarchy
)
study_capabilities(study_facts)
#>   aligned_observations timing_resolved partition_hierarchy
#> 1                 TRUE            TRUE                TRUE
#>   stable_source_revision
#> 1                   TRUE
```

The returned row reports aligned observations, resolved timing, a bound
partition hierarchy, and stable source revisions. The study’s row
lineage also records which scans censoring retained.
[`study()`](https://bbuchsbaum.github.io/crossform/reference/study.md)
does not decide whether `condition` or `motion` is a target, nuisance
term, or unused column. Those roles belong to the design model.

## Declare effects before choosing coefficient coding

The requested effects are the three condition amplitudes, in arbitrary
BOLD units.
[`condition_space()`](https://bbuchsbaum.github.io/crossform/reference/condition_space.md)
names those coordinates.
[`effect_map()`](https://bbuchsbaum.github.io/crossform/reference/effect_map.md)
requests all three, so later calls can form contrasts without refitting
the run models.

``` r

condition_coordinates <- condition_space(
  conditions,
  basis_id = "scan-level-condition-mean:v1",
  units = "arbitrary-BOLD"
)
effect_weights <- diag(length(conditions))
dimnames(effect_weights) <- list(conditions, conditions)
named_effects <- effect_map(effect_weights, condition_coordinates)
```

The compiled design has one column per condition plus linear drift and
motion.
[`coefficient_parameterization()`](https://bbuchsbaum.github.io/crossform/reference/coefficient_parameterization.md)
declares which linear combinations of those coefficients equal the named
condition amplitudes. This declaration makes the effect request
independent of whether the compiler uses cell-means or treatment coding.

``` r

condition_from_coefficients <- cbind(
  diag(length(conditions)), drift = 0, motion = 0
)
dimnames(condition_from_coefficients) <- list(
  conditions, colnames(design)
)
parameterization <- coefficient_parameterization(
  condition_from_coefficients,
  condition_coordinates,
  coding_id = "cell-means-plus-nuisance"
)

designs <- stats::setNames(
  rep(list(design), length(partitions)), partitions
)
parameterizations <- stats::setNames(
  rep(list(parameterization), length(partitions)), partitions
)
declared_design <- design_model(
  specification = list(
    target = "scan-level condition means",
    nuisance = c("linear drift", "declared motion covariate"),
    censor_policy = "use the bound retained column"
  ),
  conditions = condition_coordinates,
  designs = designs,
  parameterizations = parameterizations,
  row_ids = lapply(scan_indexes, `[[`, "observation_id"),
  solver = "qr",
  protocol = "from-observations-example",
  protocol_version = "1",
  package_version = "1"
)
```

[`design_model()`](https://bbuchsbaum.github.io/crossform/reference/design_model.md)
returns the declared design plus the metadata needed to check and
reproduce its compilation. This example supplies matrices directly. An
external compiler may instead build hemodynamic response functions and
event grammars, but it must provide the same named condition map, row
lineage, ranks, aliasing information, and solver diagnostics.

Nothing declared so far names a coefficient, and that is what makes the
request portable. A compiler that reports treatment coding — an
intercept and two differences from `face` — describes the same three
amplitudes on a different axis.
[`lower_effect_map()`](https://bbuchsbaum.github.io/crossform/reference/lower_effect_map.md)
is the step that moves the request onto whichever axis the compiler
chose;
[`plan_relation()`](https://bbuchsbaum.github.io/crossform/reference/plan_relation.md)
runs it for every partition, and calling it directly is how you inspect
what it produced.

``` r

# A second compiler's coefficients. `recode` is both the change of basis on
# the design and, in its condition rows, the declaration that the intercept
# and the two differences still name `face`, `body`, and `tool`.
recode <- cbind(
  intercept = c(1, 1, 1, 0, 0),
  body_minus_face = c(0, 1, 0, 0, 0),
  tool_minus_face = c(0, 0, 1, 0, 0),
  drift = c(0, 0, 0, 1, 0),
  motion = c(0, 0, 0, 0, 1)
)
rownames(recode) <- colnames(design)
treatment_design <- design %*% recode
treatment_parameterization <- coefficient_parameterization(
  recode[conditions, , drop = FALSE],
  condition_coordinates,
  coding_id = "treatment-plus-nuisance"
)

lowered_cell <- lower_effect_map(named_effects, parameterization)
lowered_treatment <- lower_effect_map(
  named_effects, treatment_parameterization
)

# The effect request keeps its identity; the lowering does not, because the
# two codings put the same request on different coefficient axes.
c(
  same_effect_request = identical(
    lowered_cell$effect_map_id, lowered_treatment$effect_map_id
  ),
  same_lowering = identical(
    lowered_cell$lowering_id, lowered_treatment$lowering_id
  ),
  coding_invariant = lowered_cell$capabilities$coding_invariant
)
#> same_effect_request       same_lowering    coding_invariant 
#>                TRUE               FALSE                TRUE

lowered_cell$target
#>      face body tool drift motion
#> face    1    0    0     0      0
#> body    0    1    0     0      0
#> tool    0    0    1     0      0
lowered_treatment$target
#>      intercept body_minus_face tool_minus_face drift motion
#> face         1               0               0     0      0
#> body         1               1               0     0      0
#> tool         1               0               1     0      0

# Fit one run under each coding and read the amplitudes back out.
run <- response_runs[["run-1"]]
all.equal(
  lowered_cell$target %*% qr.solve(design, run),
  lowered_treatment$target %*% qr.solve(treatment_design, run)
)
#> [1] TRUE
```

The two lowered targets are different matrices over different columns,
yet they extract the same amplitudes from their own design’s
coefficients. That agreement is the coding-invariance claim in the only
form that can be checked, and `$coding_invariant` is the capability flag
that reports it. The escape hatch below,
[`raw_effect_map()`](https://bbuchsbaum.github.io/crossform/reference/raw_effect_map.md),
sets that flag to `FALSE` precisely because a supplied target matrix has
no condition space to be invariant over.

## Declare the fit and its assumptions

The observation model states assumptions that the data tables cannot
prove. Here the estimator is ordinary least squares (OLS), scans are the
sampling unit, and runs are declared independently acquired conditional
on the model. Changing any of these declarations changes the requested
relation.

``` r

noise_model <- observation_model(
  "ols",
  sampling_unit = "scan",
  independence = "runs independently acquired conditional on the model"
)
relation_request <- plan_relation(
  study_facts,
  declared_design,
  named_effects,
  noise_model
)
relation_request
#> relation_plan<3 partitions; semantic; uncertainty analytic>
compiler_conformance(relation_request)
#>   partition semantic_identity regressor_axis condition_lowering row_lineage
#> 1     run-1              TRUE           TRUE               TRUE        TRUE
#> 2     run-2              TRUE           TRUE               TRUE        TRUE
#> 3     run-3              TRUE           TRUE               TRUE        TRUE
#>   rank_and_aliases censor_accounting solver_diagnostics whitening_provenance
#> 1             TRUE              TRUE               TRUE                 TRUE
#> 2             TRUE              TRUE               TRUE                 TRUE
#> 3             TRUE              TRUE               TRUE                 TRUE
#>   source_revision portable_receipt
#> 1            TRUE             TRUE
#> 2            TRUE             TRUE
#> 3            TRUE             TRUE
```

[`plan_relation()`](https://bbuchsbaum.github.io/crossform/reference/plan_relation.md)
compiles and validates the design metadata without reading neural
values. It returns an `effect_relation_plan`. If a requested effect is
not estimable in a run, the call stops and names the run, rank, aliased
regressors, cause, and possible remedies.

The plan, design receipt, and fitted relation have separate identities:

``` r

receipts <- relation_plan_receipts(relation_request)
fit <- estimate_relation(relation_request)
data.frame(
  object = c("relation plan", "run-1 design receipt", "fitted relation"),
  identity = substr(c(
    relation_request$relation_plan_id,
    receipts$`run-1`$design_receipt_id,
    fit$signature
  ), 1, 32)
)
#>                 object                         identity
#> 1        relation plan relation-plan-sha256:deeed8a8d2b
#> 2 run-1 design receipt design-receipt-sha256:9afd43489d
#> 3      fitted relation sha256:71ac64e07fd682d71af377936
```

Changing the solver from QR to SVD changes the design receipt, not the
named condition-space request. A supported change from cell-means to
treatment coding does the same. Changing the requested effects, units,
observation law, or independence declaration changes the relation plan.
`fit$signature` binds the plan, the realized receipts, and the exact
observation revisions used by
[`estimate_relation()`](https://bbuchsbaum.github.io/crossform/reference/estimate_relation.md).

## Ask what reproduces across runs

[`estimate_relation()`](https://bbuchsbaum.github.io/crossform/reference/estimate_relation.md)
emits the same `effect_relation_fit` consumed everywhere else in the
package.
[`cross_partitions()`](https://bbuchsbaum.github.io/crossform/reference/cross_partitions.md)
requests products between different runs; `generalizes_over = "run"`
records that run is the axis across which the effect must reproduce.

The frame decides where results are reported. This example uses three
anatomical regions, so every result below has three rows. Searchlights,
single voxels, and the whole brain are the other choices; the package
calls any one of them a *measurement*, and nothing downstream changes
when you swap one for another.

``` r

frame <- compile_frame(
  regions(rep(paste0("roi-", 1:3), each = 4L), normalization = "local"),
  domain
)
generalization <- cross_partitions(
  fit$relation,
  independence = "independent",
  generalizes_over = "run"
)
geometry_plan <- plan_geometry(fit$relation, frame, generalization)

face_minus_body <- c(face = 1, body = -1, tool = 0)
effect <- contrast_energy(geometry_plan, face_minus_body)
distances <- rdm(geometry_plan)

tool_category <- outer(
  c(face = 0, body = 0, tool = 1),
  c(face = 0, body = 0, tool = 1),
  function(left, right) as.numeric(left != right)
)
model_fit <- rsa(geometry_plan, models = list(tool = tool_category))
```

[`contrast_energy()`](https://bbuchsbaum.github.io/crossform/reference/contrast_energy.md)
returns the signed regional contrast and its coherent, configuration,
and total cross-run energies.
[`rdm()`](https://bbuchsbaum.github.io/crossform/reference/rdm.md)
returns the three pairwise condition distances for each region.
[`rsa()`](https://bbuchsbaum.github.io/crossform/reference/rsa.md)
returns an intercept and the coefficient for the named tool-category
model.

``` r

round(cbind(
  signed = effect$signed,
  coherent = effect$coherent,
  configuration = effect$configuration,
  total = effect$total
), 3)
#>      signed coherent configuration total
#> [1,]  0.278    0.076         0.178 0.255
#> [2,]  0.023   -0.001         0.113 0.112
#> [3,] -0.590    0.339         0.089 0.428
round(distances$values, 3)
#>      face - body face - tool body - tool
#> [1,]       0.255       0.409       0.064
#> [2,]       0.112       0.288       0.138
#> [3,]       0.428       0.047       0.521
round(model_fit$coefficients, 3)
#>      (Intercept)   tool
#> [1,]       0.255 -0.018
#> [2,]       0.112  0.101
#> [3,]       0.428 -0.144
```

For every region, total energy must equal coherent plus configuration
energy. The hidden check below also verifies the expected result
dimensions and every compiler-conformance guarantee.

Because the fit retained residuals and declared a fixed observation
model, the package can calculate analytic standard errors for the three
RDM distances at a selected region:

``` r

distance_covariance <- rdm_sampling_covariance(
  geometry_plan, fit, target = "plugin", at = 1L
)
round(sqrt(sampling_covariance(distance_covariance)), 3)
#> face - body face - tool body - tool 
#>       0.079       0.097       0.054
```

These are within-participant standard errors under the declared fixed,
separable error model, which represents sampling and neural-feature
covariance as separate factors. They are not population-level
uncertainty. A learned whitener would still permit the point estimates
but would not receive this analytic uncertainty capability.

### The assumption those standard errors rest on

The analytic law assumes the observations *within a run* are independent
given the design. fMRI residuals are not — they are temporally
autocorrelated — and the consequence is not a rounding error. Under
AR(1) errors with $`\rho = 0.75`$ (32 trials, four conditions, six runs,
50 features), the true spread of the distance estimator divided by the
standard error crossform reports is:

``` text
randomly interleaved trial order   0.50   (standard error twice too large)
blocked trial order                5.10   (standard error five times too small)
blocked, with a correct whitener   1.03
```

A blocked design without a whitener is the dangerous case: the reported
standard error is anticonservative by a factor of five, because the
estimates in a block share the same slow drift. The direction flips for
rapid event-related designs, where neighboring trials belong to
different conditions and the same autocorrelation *removes* variance
from a contrast.

The remedy is a whitener $`L`$ with $`L^\top L = \Sigma_t^{-1}`$, so the
fit is carried out on the whitened design $`L X`$ and the whitened
responses $`L y`$. For a known AR(1) coefficient that is one line, and
it is declared the same way as every other assumption — as part of the
observation model, which changes the relation’s identity:

``` r

ar1_whitener <- function(observations, rho) {
  covariance <- rho^abs(outer(
    seq_len(observations), seq_len(observations), "-"
  ))
  t(backsolve(chol(covariance), diag(observations)))
}

gls_model <- observation_model(
  "fixed_gls",
  sampling_unit = "scan",
  whitener = ar1_whitener(n_scan, 0.3),
  independence = "runs independently acquired conditional on the model"
)
gls_fit <- estimate_relation(
  plan_relation(study_facts, declared_design, named_effects, gls_model)
)

# A different observation model is a different relation, not a variant of the
# same one, and the fit says so.
c(ols = substr(fit$signature, 1, 24),
  gls = substr(gls_fit$signature, 1, 24))
#>                        ols                        gls 
#> "sha256:71ac64e07fd682d71" "sha256:80844531a3eee55ad"
```

[`lm_relation_fit()`](https://bbuchsbaum.github.io/crossform/reference/lm_relation_fit.md)
takes the same matrix as `observation_whitener = L` on the short route
below. crossform applies the $`L`$ you supply and records its revision
in the observation model’s identity; it cannot check that $`L`$ matches
the autocorrelation actually present in your data, and estimating
$`\rho`$ from the same residuals you are about to calibrate against is a
different procedure that this package does not validate. See
[`?lm_relation_fit`](https://bbuchsbaum.github.io/crossform/reference/lm_relation_fit.md)
and [the evidence-sampling
contract](https://github.com/bbuchsbaum/crossform/blob/main/design/evidence-sampling-contract.md)
for the full statement.

## Shorter routes: a compiled design, or beta matrices

The ladder above declares the first-level model so that the request
survives a change of coding, a change of solver, or a change of
compiler. Not every project starts there. Two shorter routes reach the
same
[`plan_geometry()`](https://bbuchsbaum.github.io/crossform/reference/plan_geometry.md)/[`contrast_energy()`](https://bbuchsbaum.github.io/crossform/reference/contrast_energy.md)
step from less input, and they differ in exactly one consequence:
whether the residual channel survives, and with it the standard errors
calculated above.
[`vignette("introduction")`](https://bbuchsbaum.github.io/crossform/articles/introduction.md)
builds its worked fixture on the first of these and its own-data
walkthrough on the second.

### From raw responses and a compiled design

If you kept the scan-by-feature responses and already have a design
matrix,
[`lm_relation_fit()`](https://bbuchsbaum.github.io/crossform/reference/lm_relation_fit.md)
is the compact route that keeps the error channel. It skips
[`study()`](https://bbuchsbaum.github.io/crossform/reference/study.md),
[`design_model()`](https://bbuchsbaum.github.io/crossform/reference/design_model.md),
and
[`plan_relation()`](https://bbuchsbaum.github.io/crossform/reference/plan_relation.md),
so nothing checks the design’s semantics against events, confounds, and
clocks for you; in exchange it is a single call. The `effects` matrix
says which linear combination of design columns each named condition is.

``` r

set.seed(20260815)
condition_names <- c("face", "body", "house", "tool")
native_domain <- abstract_domain(60L, id = "short-routes:v1")

raw_n_scan <- 40L
scan_condition <- factor(
  rep(condition_names, length.out = raw_n_scan), levels = condition_names
)
design_matrix <- cbind(
  stats::model.matrix(~ 0 + scan_condition),
  drift = seq(-1, 1, length.out = raw_n_scan)
)
colnames(design_matrix) <- c(condition_names, "drift")

target_matrix <- cbind(diag(length(condition_names)), drift = 0)
dimnames(target_matrix) <- list(condition_names, colnames(design_matrix))

# Substitute your own: one scan-by-feature response matrix per run.
signal <- rbind(
  matrix(rnorm(length(condition_names) * 60L, sd = 0.5), length(condition_names)),
  rnorm(60L, sd = 0.1)
)
raw_response_runs <- lapply(seq_len(3), function(run) {
  design_matrix %*% signal + matrix(rnorm(raw_n_scan * 60L), raw_n_scan, 60L)
})
names(raw_response_runs) <- paste0("run", seq_len(3))

raw_fit <- lm_relation_fit(
  raw_response_runs,
  design = design_matrix,
  effects = target_matrix,
  effect_names = effect_space(condition_names, units = "percent-signal"),
  sampling_unit = "trial",
  domain = native_domain
)
raw_plan <- plan_geometry(
  raw_fit$relation,
  compile_frame(whole_brain(), native_domain),
  cross_partitions(
    raw_fit$relation,
    independence = "independent", generalizes_over = "run"
  )
)
sampling_capabilities(raw_plan, raw_fit)
#> <effect_sampling_capabilities>
#>   analytic sampling law: available 
#>   metric: fixed | partitions: equal | error channel: relation_fit
```

[`whole_brain()`](https://bbuchsbaum.github.io/crossform/reference/whole_brain.md)
reports one measurement over every feature in the domain, which is the
simplest frame to check a pipeline against. Because this fit retained
residuals, the analytic law is admitted and the standard errors from the
section above are available here too:

``` r

round(sqrt(sampling_covariance(
  rdm_sampling_covariance(raw_plan, raw_fit, target = "null", at = 1L)
)), 4)
#>  face - body face - house  face - tool body - house  body - tool house - tool 
#>       0.0150       0.0150       0.0151       0.0150       0.0150       0.0150
```

### From beta matrices alone

If the responses are gone and only condition-by-feature estimates per
run remain, build a point relation directly. It needs a named list of
matrices with one row per condition, a named effect space, and a neural
domain.

``` r

# Substitute your own: one condition-by-feature matrix per run, rows named.
beta_runs <- lapply(seq_len(3), function(run) matrix(
  rnorm(length(condition_names) * native_domain$n_features),
  nrow = length(condition_names),
  dimnames = list(condition_names, native_domain$feature_ids)
))
names(beta_runs) <- paste0("run", seq_len(3))

rel <- relation(
  beta_runs,
  effects = effect_space(condition_names, units = "percent-signal"),
  domain = native_domain
)
own_plan <- plan_geometry(
  rel,
  compile_frame(whole_brain(), native_domain),
  cross_partitions(rel, independence = "independent", generalizes_over = "run")
)
own_plan
#> <effect_geometry_plan>
#>   effects:      4 x 4
#>   measurements: 1
#>   features:     60
#>   metric:       implicit identity
#>   generalizes:  3 partition pairs over run, endpoints independent
#>   execution:    query-first, in memory
#>   state:        nothing computed yet
#>   next:         contrast_energy(plan, weights), rdm(plan), rsa(plan)
```

`own_plan` supports
[`contrast_energy()`](https://bbuchsbaum.github.io/crossform/reference/contrast_energy.md),
[`rdm()`](https://bbuchsbaum.github.io/crossform/reference/rdm.md), and
[`rsa()`](https://bbuchsbaum.github.io/crossform/reference/rsa.md)
exactly as the fitted relations above do. What it cannot support is the
analytic standard error: beta matrices carry no residuals and no
residual degrees of freedom, so
[`sampling_capabilities()`](https://bbuchsbaum.github.io/crossform/reference/sampling_capabilities.md)
reports the law unavailable and names `missing_error_channel` as the
unmet requirement. Give the domain a mask and voxel spacing instead of a
bare feature count and searchlights become available;
[`vignette("introduction")`](https://bbuchsbaum.github.io/crossform/articles/introduction.md)
works that variant end to end.

## How do BIDS and external compilers fit?

[`bids_study()`](https://bbuchsbaum.github.io/crossform/reference/bids_study.md)
adapts one task `events.tsv` and one row-aligned confound table per
partition to the same
[`study()`](https://bbuchsbaum.github.io/crossform/reference/study.md)
object used above. It preserves extra event and confound columns. The
design compiler, not the adapter, decides which columns become target or
nuisance terms. If censoring is intended, the call must name the Boolean
censor column. This keeps BIDS at the file boundary; the core study
object is not BIDS-shaped.

The signature mirrors the
[`study()`](https://bbuchsbaum.github.io/crossform/reference/study.md)
call above: an
[`observations()`](https://bbuchsbaum.github.io/crossform/reference/observations.md)
record, then one events file and one confounds file per partition, in
partition order.

``` r

runs <- paste0("run-", 1:3)

bids_facts <- bids_study(
  observations(response_runs, scan_indexes, domain),
  event_files = file.path(
    "sub-01", "func",
    sprintf("sub-01_task-objects_%s_events.tsv", runs)
  ),
  confound_files = file.path(
    "sub-01", "func",
    sprintf("sub-01_task-objects_%s_desc-confounds_timeseries.tsv", runs)
  ),
  partitions = runs,
  censor = "retained",
  units = "seconds"
)
```

[`fmridesign_design_model()`](https://bbuchsbaum.github.io/crossform/reference/fmridesign_design_model.md)
adapts the supported `fmridesign` version into a design that names
condition effects independently of coefficient columns.
[`fmrireg_relation()`](https://bbuchsbaum.github.io/crossform/reference/fmrireg_relation.md)
can execute an OLS relation plan for point-estimate parity. That adapter
does not return an analytic error channel because crossform cannot
verify the required residual covariance contract from its result.

If you have compiled design and target matrices but not their semantic
condition mapping, and you still want the declared ladder rather than
the short route above, use
[`raw_design_model()`](https://bbuchsbaum.github.io/crossform/reference/raw_design_model.md)
with
[`raw_effect_map()`](https://bbuchsbaum.github.io/crossform/reference/raw_effect_map.md).
Their numeric values then become part of the plan identity, and the fit
cannot claim that a different coefficient coding represents the same
request.

The next call depends on the result you need. Use
[`contrast_energy()`](https://bbuchsbaum.github.io/crossform/reference/contrast_energy.md)
for a named effect,
[`rdm()`](https://bbuchsbaum.github.io/crossform/reference/rdm.md) for
condition distances,
[`rsa()`](https://bbuchsbaum.github.io/crossform/reference/rsa.md) for a
fixed linear model, and
[`rdm_sampling_covariance()`](https://bbuchsbaum.github.io/crossform/reference/rdm_sampling_covariance.md)
only after
[`sampling_capabilities()`](https://bbuchsbaum.github.io/crossform/reference/sampling_capabilities.md)
confirms that the fit and plan meet its assumptions.
