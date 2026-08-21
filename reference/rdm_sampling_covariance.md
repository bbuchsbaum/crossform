# Construct exact analytic sampling covariance for crossvalidated distances

Builds the local sampling-covariance form from Diedrichsen, Provost, and
Zareamoghaddam (2016), Eq. 13, for an equal-weight all-partition-pairs
RDM. The admitted specialization requires a geometry plan built from an
[`lm_relation_fit()`](https://bbuchsbaum.github.io/crossform/reference/lm_relation_fit.md)
relation, equal partition error structures, independent partition
endpoints, and one common fixed neural metric. The result remains
factorized and queryable; it does not allocate the complete distance-by-
distance covariance.

## Usage

``` r
rdm_sampling_covariance(
  x,
  fit,
  target,
  at,
  residual_strategy = c("node_local", "shared_pair_statistics"),
  residual_workspace_bytes = 512 * 1024^2
)
```

## Arguments

- x:

  An `effect_geometry_plan` using
  [`cross_partitions()`](https://bbuchsbaum.github.io/crossform/reference/cross_partitions.md)
  and a fixed neural metric.

- fit:

  The identity-bound `effect_relation_fit` that supplied
  `x$task$left_relation`.

- target:

  Whether the signal-dependent term is evaluated at the partition-mean
  plug-in estimate (`"plugin"`) or at the fixed zero null (`"null"`).
  These are different calibration policies and are never chosen
  implicitly. `"plugin"` substitutes the partition mean of the
  *estimates* for the unknown signal, and because
  \\E\[\hat\mu_r\Sigma_w\hat\mu_s^\top\] = \mu_r\Sigma_w\mu_s^\top +
  \Xi\_{rs}\\\mathrm{tr}(\Sigma_w\Sigma_w)/M\\, its signal term is
  biased upward by \\4\Xi\_{rs}^2\mathrm{tr}(\Sigma_w\Sigma_w)/M^2\\, an
  inflation that shrinks like \\1/M^2\\ and is largest when the noise
  dominates the true distances. Prefer `"null"` for calibrating a test
  of no effect, where it is exact; use `"plugin"` when reporting
  uncertainty around an estimated nonzero distance and read it as mildly
  conservative.

- at:

  One measurement position, or a vector of positions, in the compiled
  frame. Required, with no default: the analytic law is local, so a
  covariance without a named measurement would be a covariance of
  nothing in particular. A length-1 `at` keeps the historical
  single-node object. A longer `at` compiles the plan, contrast
  transport, and any eligible shared residual statistics once, then
  returns one covariance object per requested node.

- residual_strategy:

  `"node_local"` reads residual blocks only for the requested
  measurement. `"shared_pair_statistics"` explicitly compiles reusable
  residual pair sufficient statistics for a batch of overlapping
  supports.

- residual_workspace_bytes:

  Positive crossform-owned workspace budget for shared residual pair
  statistics.

## Value

An `effect_sampling_covariance` with `basis = "rdm"` for one
measurement, or an `effect_sampling_covariance_batch` list of those
objects when `at` names several measurements. Each object is queryable
by
[`sampling_covariance()`](https://bbuchsbaum.github.io/crossform/reference/sampling_covariance.md).
It contains within-measurement uncertainty only; it does not imply
covariance between spatial locations. Its `$source` records
`residual_df` (\\\nu\\), `residual_effective_dimension`
(\\P\_{\mathrm{eff}}\\), and `noise_trace_estimator`. A batch also
records the shared-resource execution route on each element's
`$source$execution`.

## Details

Writing \\\mu_r\\ for the whitened contrast pattern of distance \\r\\,
\\\Sigma_w\\ for the residual covariance in the same whitened
coordinates, \\\Xi\\ for the effect-coordinate contrast cross-products,
and \\M\\ for the partition count, the law evaluated here is
\$\$\mathrm{Cov}(d_r, d_s) = \frac{4}{M}\\\Xi\_{rs}\\
\mu_r\Sigma_w\mu_s^\top +
\frac{2}{M(M-1)}\\\Xi\_{rs}^2\\\mathrm{tr}(\Sigma_w\Sigma_w).\$\$ The
residual covariance enters the signal term as a metric on the whitened
patterns, so the law is correct for a general anisotropic \\\Sigma_w\\
and not only for the spherical case.

## What is exact and what is estimated

Under `target = "null"` the law is exact *on the variance scale*: the
signal term vanishes and the remaining expression is the true variance
of the estimator, not an approximation to it.

What is estimated is \\\Sigma_w\\ itself. crossform substitutes the
partition-pooled sample residual covariance \\S_w\\, a plug-in with
\\\nu=\sum_m \mathrm{df}\_m\\ degrees of freedom. That matters because
the signal-independent term is *quadratic* in \\\Sigma_w\\. For
\\S_w\sim W_P(\nu,\Sigma_w)/\nu\\, \$\$E\\\mathrm{tr}(S_w^2) =
\frac{\nu+1}{\nu}\mathrm{tr}(\Sigma_w^2) +
\frac{\mathrm{tr}(\Sigma_w)^2}{\nu},\$\$ so taking
\\\mathrm{tr}(S_w^2)\\ at face value overstates the reported standard
error by \\\sqrt{1+(1+P\_{\mathrm{eff}})/\nu}\\, where
\\P\_{\mathrm{eff}}=\mathrm{tr}(\Sigma_w)^2/\mathrm{tr}(\Sigma_w^2)\\ is
the number of residual directions the support actually spends its
variance on. Since 2026-08-16 the quadratic term therefore uses the
Wishart-unbiased estimator \$\$\widehat{\mathrm{tr}}(\Sigma_w^2) =
\frac{\nu^2}{(\nu-1)(\nu+2)}
\left(\mathrm{tr}(S_w^2)-\frac{\mathrm{tr}(S_w)^2}{\nu}\right).\$\$ The
signal term is linear in \\\Sigma_w\\ and needs no correction. Voxelwise
frames, where \\P\_{\mathrm{eff}}=1\\, were never materially affected; a
50-voxel searchlight at \\\nu=168\\ was reporting standard errors 14%
too large and an 800-voxel one more than twice too large.

Both numbers are reported: \\\nu\\ is `$source$residual_df` and
\\P\_{\mathrm{eff}}\\ is `$source$residual_effective_dimension`, and the
[`print()`](https://rdrr.io/r/base/print.html) method shows them. When
\\\nu\<P\_{\mathrm{eff}}\\ there is no usable estimate of the quadratic
term at all, and the call refuses with capability
`"sufficient_residual_df"` rather than returning a confidently small
number.

## Independence within a partition

The sampling law assumes the observations within a partition are
independent given the design. fMRI residuals are not: they are
temporally autocorrelated. Fitting without
`lm_relation_fit(observation_whitener = )` leaves \\\Xi\\ as the plain
OLS factor \\(X^\top X)^{-1}\\, which does not describe the covariance
of the estimates under correlated errors, and leaves \\\nu\\ counting
observations rather than independent ones, so \\\nu\\ overstates the
residual information. The reported standard error can then err in
*either* direction, and the size is not small. Under AR(1) errors with
\\\rho=0.75\\, 32 trials, four conditions, six runs and 50 features, the
ratio of the true spread to the reported standard error is

Passing a whitener \\L\\ with \\L^\top L=\Sigma_t^{-1}\\ makes the
whitened problem satisfy the assumption and restores calibration;
crossform cannot check that the \\L\\ you supply matches the
autocorrelation actually present in your data. See
[`lm_relation_fit()`](https://bbuchsbaum.github.io/crossform/reference/lm_relation_fit.md),
[`observation_model()`](https://bbuchsbaum.github.io/crossform/reference/observation_model.md),
and
[`vignette("from-observations")`](https://bbuchsbaum.github.io/crossform/articles/from-observations.md).

## References

Diedrichsen J, Provost S, Zareamoghaddam H (2016), "On the distribution
of cross-validated Mahalanobis distances", especially Eqs. 10, 13, and
35 and Section 5.1.
[doi:10.48550/arXiv.1607.01371](https://doi.org/10.48550/arXiv.1607.01371)

Srivastava MS (2005), "Some tests concerning the covariance matrix in
high dimensional data", *Journal of the Japan Statistical Society*
35(2), 251–272, for the unbiased estimator of \\\mathrm{tr}(\Sigma^2)\\
used here.

## See also

[`sampling_covariance()`](https://bbuchsbaum.github.io/crossform/reference/sampling_covariance.md)
to query the result, including `queries = ` for the covariance of a bank
of
[`contrast_energy()`](https://bbuchsbaum.github.io/crossform/reference/contrast_energy.md)
queries at one measurement,
[`sampling_capabilities()`](https://bbuchsbaum.github.io/crossform/reference/sampling_capabilities.md)
to ask whether the law is available before provoking a refusal, and
[`rdm()`](https://bbuchsbaum.github.io/crossform/reference/rdm.md) for
the point estimates it describes.

Other sampling uncertainty:
[`reconstruct_evidence()`](https://bbuchsbaum.github.io/crossform/reference/reconstruct_evidence.md),
[`sampling_capabilities()`](https://bbuchsbaum.github.io/crossform/reference/sampling_capabilities.md),
[`sampling_covariance()`](https://bbuchsbaum.github.io/crossform/reference/sampling_covariance.md)

## Examples

``` r
set.seed(23)
design <- model.matrix(~ 0 + factor(rep(c("a", "b"), each = 4)))
colnames(design) <- c("a", "b")
effects <- diag(2)
rownames(effects) <- colnames(design)
truth <- rbind(a = c(0.4, 0.1, 0), b = c(0, 0.2, 0.5))
responses <- setNames(lapply(seq_len(3), function(run) {
  design %*% truth + matrix(rnorm(8 * 3, sd = 0.2), 8, 3)
}), paste0("run", seq_len(3)))
domain <- abstract_domain(3, id = "sampling-example")
fit <- lm_relation_fit(
  responses, design, effects, effect_names = rownames(truth),
  sampling_unit = "trial", domain = domain
)
metric <- noise_precision(
  diag(3), domain, covariance = diag(3),
  provenance = list(source = "fixed-example-metric")
)
plan <- plan_geometry(
  fit$relation, compile_frame(whole_brain(), domain),
  cross_partitions(fit$relation, independence = "independent"),
  metric = metric
)
uncertainty <- rdm_sampling_covariance(plan, fit, target = "null", at = 1L)
sqrt(sampling_covariance(uncertainty))
#>       a - b 
#> 0.005913026 

# The residual channel behind the quadratic noise term is reported, not
# assumed: nu degrees of freedom against P_eff effective directions.
uncertainty$source$residual_df
#> [1] 18
round(uncertainty$source$residual_effective_dimension, 3)
#> [1] 2.971
```
