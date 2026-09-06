# Fit a regularized representational form

Learn a positive semidefinite prediction from signed cross-partition
geometry in a declared model family. Model eigenvalues determine the
inverse-kernel penalty; overlapping models share one union space.
Fitting is biased and adaptive. Independent predictive evidence is a
separate operation on a frozen fit and disjoint observation origins.

## Usage

``` r
fit_geometry(
  plan,
  basis,
  weights = NULL,
  rank,
  penalty,
  component = c("total", "coherent", "configuration"),
  storage = c("memory", "block"),
  storage_path = NULL,
  row_block = 128L,
  retain_source = FALSE,
  loss = "frobenius"
)
```

## Arguments

- plan:

  A geometry plan over the training pairing. Only positive-weight
  pairing endpoints are read. The current complete-form executor admits
  undirected self forms with implicit identity or a fixed SPD neural
  metric.

- basis:

  A data-independent
  [`model_basis()`](https://bbuchsbaum.github.io/crossform/reference/model_basis.md),
  preferably declaring `conditions = relation$effect_space` and
  `normalize = "trace"`.

- weights:

  Named nonnegative model weights summing to one. Required for multiple
  models; a single model defaults to weight one. Zero weights remove
  that model's exclusive directions from admissible support.

- rank:

  Required nonnegative response rank budget. Zero gives the zero
  predictor; an oversized budget clamps to positive pooled support. A
  cut through a positive tied eigenspace refuses rather than choosing an
  arbitrary prediction. This rank is separate from each model's input
  rank.

- penalty:

  Required finite nonnegative inverse-kernel penalty. Declare it before
  evaluation or choose it using training-only validation.

- component:

  The signed target: `"total"`, `"coherent"`, or `"configuration"`. It
  remains fixed when the prediction is evaluated.

- storage:

  `"memory"` or `"block"` for signed computation and retained prediction
  factors/spectra.

- storage_path:

  New durable directory for block storage. A completed record is saved
  as `fit.rds`; existing directories are never overwritten.

- row_block:

  Positive number of measurements processed together.

- retain_source:

  Retain the signed compressed training form as `$source`. Otherwise
  only its signed spectrum and fit diagnostics are retained.

- loss:

  The implemented loss is `"frobenius"` on effect geometry. RDM or
  covariance-weighted loss requires a different numerical solver.

## Value

A frozen `effect_geometry_fit` with model/support specifications,
parameters, exact target identities, actual training dependencies,
reduced prediction factors and spectra, diagnostics, and execution
receipt. It has no raw-geometry unbiasedness or inference capability.
[`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html) reads its
per-measurement diagnostics without reading neural sources. Block
records reopen with `readRDS("<path>/fit.rds")`.

## Details

For pooled kernel `K = Q U D U' Q'` and signed reduced form
`S = Q' G Q`, the estimator is `F = Q A Q'`, where
`A = U [U' S U - penalty * D^-1]_{+,rank} U'`. The bracket keeps only
the largest positive roots within the rank budget. The source is never
clipped before compression or penalty shifting. Diagnostics distinguish
negative source mass, negative shifted mass and truncated shifted mass.

The reported residual and objective are in the model union coordinates.
Compression does not identify the whole-geometry residual or an
explained fraction. A zero-penalty full-span model is labelled a generic
baseline. Positive anisotropic penalties retain model-specific
preferences. With a trace-normalized model kernel the penalty has
geometry units; changing the neural frame normalization changes its
appropriate scale. The source and any retained source spectrum stay
signed. See
[`vignette("predictive-geometry")`](https://bbuchsbaum.github.io/crossform/articles/predictive-geometry.md)
for independent scoring, model-rank choices and the distinction from
new-condition prediction.

## See also

[`model_basis()`](https://bbuchsbaum.github.io/crossform/reference/model_basis.md),
[`score_geometry()`](https://bbuchsbaum.github.io/crossform/reference/score_geometry.md),
[`model_geometry()`](https://bbuchsbaum.github.io/crossform/reference/model_geometry.md)

## Examples

``` r
d <- abstract_domain(3, id = "fit-example")
effect_names <- c("face", "body", "house", "tool")
B <- matrix(seq_len(12)/10, 4, dimnames = list(effect_names, NULL))
rel <- relation(list(run1 = B, run2 = B + 0.05), domain = d)
model <- model_basis(features = matrix(c(1, 1, -1, -1), 4),
  conditions = rel$effect_space, normalize = "trace")
plan <- plan_geometry(rel, compile_frame(whole_brain(), d),
  cross_partitions(rel, independence = "independent", generalizes_over = "run"))
fit <- fit_geometry(plan, model, rank = 1, penalty = 0.01)
fit
#> geometry_fit<1 measurements; model support 1; rank <= 1; penalty 0.01; isotropic_shrinkage>
#> Learned PSD prediction; fitting is biased. Source spectra remain signed.
#> Training: 2 used partitions; component: total.
#> Training ancestry is undeclared; independent scoring is unavailable.
#>  measurement rank prediction_norm_sq compressed_residual_sq
#>  whole_brain    1              9e-04                  1e-04
as.data.frame(fit)
#>   measurement rank prediction_norm_sq penalty_trace compressed_residual_sq
#> 1 whole_brain    1              9e-04          0.03                  1e-04
#>   objective source_negative_mass shifted_negative_mass
#> 1   0.00035                    0                     0
#>   shifted_truncated_positive_mass
#> 1                               0
```
