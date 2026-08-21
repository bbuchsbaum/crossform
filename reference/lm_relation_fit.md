# Fit a partitioned linear-model relation with a residual error channel

The adapter compiles each supplied design once. Its pure effect map
forms the underlying relation; its orthogonal residual projection,
unscaled effect-coordinate covariance, residual degrees of freedom, and
source revision form a separate separable-GLM error capability. Residual
blocks are produced lazily without constructing a dense observation
residualizer. This error channel is required by the admitted
fixed-metric analytic RDM covariance in
[`rdm_sampling_covariance()`](https://bbuchsbaum.github.io/crossform/reference/rdm_sampling_covariance.md).
The separable GLM is a declared capability: it assumes the supplied
observation model yields a common effect-coordinate covariance structure
across neural features and, for the equal-partition specialization,
across partitions.

## Usage

``` r
lm_relation_fit(
  sources,
  design,
  effects,
  observation_whitener = NULL,
  effect_names = if (is.matrix(effects)) {
     rownames(effects)
 } else {
     NULL
 },
  tolerance = sqrt(.Machine$double.eps),
  source_dims = NULL,
  partitions = NULL,
  domain = NULL,
  domain_id = "abstract",
  capabilities = NULL,
  sampling_unit = "observation",
  provenance = list(),
  whiten = NULL,
  solver = "auto"
)
```

## Arguments

- sources:

  Raw observation-by-feature sources accepted by
  [`relation()`](https://bbuchsbaum.github.io/crossform/reference/relation.md).

- design:

  One design matrix, or one per partition.

- effects:

  One effect target matrix, or one per partition.

- observation_whitener:

  NULL, one observation whitener \\L\\, or one per partition. A finite
  square observation-by-observation matrix; the fit is carried out on
  \\LX\\ and \\Ly\\, so \\L\\ should satisfy \\L^\top L=\Sigma_t^{-1}\\
  for the within-partition error covariance \\\Sigma_t\\. Leaving it
  `NULL` asserts that the observations are already independent given the
  design; see the section above.

- effect_names:

  One common effect space/name vector, or one per partition.

- tolerance:

  Positive rank and estimability tolerance.

- source_dims, partitions, domain, domain_id, capabilities:

  Passed to
  [`relation()`](https://bbuchsbaum.github.io/crossform/reference/relation.md).

- sampling_unit:

  One nonempty sampling-unit label, or one per partition.

- provenance:

  Compact estimator provenance.

- whiten:

  Deprecated alias for `observation_whitener`.

- solver:

  One numerical route, or one per partition: automatic, QR, or SVD. The
  route is recorded in estimator provenance.

## Value

An `effect_relation_fit` whose `$relation` carries the fitted effect map
and whose `$error_models` carry, per partition, the unscaled
`effect_covariance`, `residual_df`, the lazy `residual_source`, the
`observation_whitener` descriptor, and `estimator_provenance`. Every
partition reports `within_participant_calibration`.

## Independent observations, or a whitener

The error channel this function records — the effect covariance \\\Xi=A
A^\top\\ and the residual degrees of freedom \\\nu=n-\mathrm{rank}(X)\\
— describes observations that are independent given the design. fMRI
residuals are not: they are temporally autocorrelated. Fitting without
`observation_whitener` therefore leaves \\\Xi\\ as the plain OLS factor,
which is not the covariance of the estimates under correlated errors,
and leaves \\\nu\\ counting observations rather than independent ones,
so \\\nu\\ overstates the residual information available.

The consequence for
[`rdm_sampling_covariance()`](https://bbuchsbaum.github.io/crossform/reference/rdm_sampling_covariance.md)
runs in *either* direction and is not small. Under AR(1) errors with
\\\rho=0.75\\, 32 trials, four conditions, six runs and 50 features, the
ratio of the true spread of the distance estimator to the standard error
crossform reports is 0.50 for a randomly interleaved trial order (the
standard error is twice too large) and 5.10 for a blocked order (five
times too small). Supplying `observation_whitener = L` with \\L^\top
L=\Sigma_t^{-1}\\ brings the same blocked design to 1.03. crossform
applies the \\L\\ you give it and records its identity; it cannot check
that \\L\\ matches the autocorrelation actually present in your data.

## See also

[`relation()`](https://bbuchsbaum.github.io/crossform/reference/relation.md)
for precomputed effects with no error channel,
[`estimate_relation()`](https://bbuchsbaum.github.io/crossform/reference/estimate_relation.md)
to reach this object from a validated
[`plan_relation()`](https://bbuchsbaum.github.io/crossform/reference/plan_relation.md),
and
[`residual_block()`](https://bbuchsbaum.github.io/crossform/reference/residual_block.md),
[`effect_covariance()`](https://bbuchsbaum.github.io/crossform/reference/effect_covariance.md),
[`residual_df()`](https://bbuchsbaum.github.io/crossform/reference/residual_df.md)
to read the error channel.

Other relation planning and fitting:
[`adapter_version_certificate()`](https://bbuchsbaum.github.io/crossform/reference/adapter_version_certificate.md),
[`compiler_conformance()`](https://bbuchsbaum.github.io/crossform/reference/compiler_conformance.md),
[`effect_extractor()`](https://bbuchsbaum.github.io/crossform/reference/effect_extractor.md),
[`estimate_relation()`](https://bbuchsbaum.github.io/crossform/reference/estimate_relation.md),
[`file_matrix_source()`](https://bbuchsbaum.github.io/crossform/reference/file_matrix_source.md),
[`fmridesign_design_model()`](https://bbuchsbaum.github.io/crossform/reference/fmridesign_design_model.md),
[`fmrireg_relation()`](https://bbuchsbaum.github.io/crossform/reference/fmrireg_relation.md),
[`lm_extractor()`](https://bbuchsbaum.github.io/crossform/reference/lm_extractor.md),
[`plan_relation()`](https://bbuchsbaum.github.io/crossform/reference/plan_relation.md),
[`relation()`](https://bbuchsbaum.github.io/crossform/reference/relation.md),
[`relation_block()`](https://bbuchsbaum.github.io/crossform/reference/relation_block.md),
[`relation_fit()`](https://bbuchsbaum.github.io/crossform/reference/relation_fit.md),
[`relation_fit_capabilities()`](https://bbuchsbaum.github.io/crossform/reference/relation_fit_capabilities.md),
[`relation_plan_receipts()`](https://bbuchsbaum.github.io/crossform/reference/relation_plan_receipts.md),
[`residual_block()`](https://bbuchsbaum.github.io/crossform/reference/residual_block.md),
[`residual_df()`](https://bbuchsbaum.github.io/crossform/reference/residual_df.md),
[`source_capabilities()`](https://bbuchsbaum.github.io/crossform/reference/source_capabilities.md)

## Examples

``` r
set.seed(20260815)
domain <- abstract_domain(6L, id = "lm-relation-fit-example")
condition <- factor(rep(c("face", "body", "tool"), each = 4L))
design <- stats::model.matrix(~ 0 + condition)
colnames(design) <- c("face", "body", "tool")

# Ask for the three condition means themselves.
targets <- diag(3)
dimnames(targets) <- list(colnames(design), colnames(design))
runs <- lapply(c("run-1", "run-2"), function(partition) {
  design %*% matrix(rnorm(18), 3L, 6L) + matrix(rnorm(72), 12L, 6L)
})
names(runs) <- c("run-1", "run-2")

fit <- lm_relation_fit(
  runs, design, targets, domain = domain, sampling_unit = "trial"
)
relation_fit_capabilities(fit)$within_participant_calibration
#> [1] TRUE TRUE

# Twelve trials minus three estimated means leaves nine residual df, which
# is what later analytic uncertainty is calibrated against.
residual_df(fit, "run-1")
#> [1] 9
```
