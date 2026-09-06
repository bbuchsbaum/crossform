# Score a frozen representational prediction on independent geometry

Measure gain over the zero-geometry predictor: twice the signed held-out
inner product minus the fitted prediction's squared Frobenius norm.
Positive gain improves squared geometry prediction error in expectation
under the declared independence assumptions. Gain and mode evidence may
be negative; no fitting or selection occurs during scoring.

## Usage

``` r
score_geometry(
  fit,
  plan,
  storage = c("memory", "block"),
  storage_path = NULL,
  row_block = 128L,
  modes = FALSE,
  components = FALSE
)
```

## Arguments

- fit:

  A frozen
  [`fit_geometry()`](https://bbuchsbaum.github.io/crossform/reference/fit_geometry.md)
  result.

- plan:

  Geometry plan on disjoint observation origins from the same declared
  parent manifest, with the same effects, fixed neural metric,
  generalization axis and spatial measurements. Measurement IDs are
  matched explicitly, so row order may differ. Both neural cross-product
  independence and fitting/evaluation independence are required.

- storage:

  `"memory"` or `"block"` for evaluation forms and mode evidence.

- storage_path:

  A new durable directory for block storage. The complete score is saved
  as `score.rds`; existing directories are never overwritten.

- row_block:

  Maximum number of measurement rows processed together.

- modes:

  Retain signed per-mode evidence and gain. Tied modes have arbitrary
  individual orientations; use the grouped view for invariant evidence
  over a fully retained tied eigenspace.

- components:

  For a total-geometry fit, retain coherent and configuration inner
  products read with the same frozen prediction. They add before the
  prediction cost is subtracted once. Separately fitted latent
  components do not in general conserve the signed decomposition.

## Value

An `effect_geometry_score` with per-measurement gain, rank, inner
product and prediction norm, actual training/evaluation dependencies,
validation assumptions and an execution receipt.
[`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html) returns
the summary; `view = "modes"` or `"groups"` reads retained details. Gain
has squared geometry units. There is no explained fraction, generic
standard error or raw crossvalidated geometry inheritance.

## Details

Conditional on all training and selection data, the mean gain is
`||G*||^2 - ||G* - F||^2` when the independent test form is unbiased for
`G*`. In particular, under a zero signal the expected gain is
`-||F||^2`, not zero. Rank, penalty, model weights and preprocessing
must be fixed without consulting these evaluation outcomes.
Same-condition partition prediction does not establish generalization to
new conditions.

The summary view retains measurement ID, gain, effective rank, held-out
inner product and prediction cost. The mode view separates fitted
amplitude, signed test evidence and gain. The group view sums these over
a fully retained tied eigenspace. `rows` selects measurement rows by
position, preserving order and repetitions; empty selections retain the
corresponding data-frame schema. Mode/group views require `modes = TRUE`
at scoring unless the rank budget is zero.

See
[`vignette("predictive-geometry")`](https://bbuchsbaum.github.io/crossform/articles/predictive-geometry.md)
for a complete fixed-split workflow.

## See also

[`fit_geometry()`](https://bbuchsbaum.github.io/crossform/reference/fit_geometry.md),
[`model_basis()`](https://bbuchsbaum.github.io/crossform/reference/model_basis.md)

## Examples

``` r
d <- abstract_domain(2, id = "score-example")
B <- matrix(c(1, -1, 0, 0, 0, 0, 1, -1), 4,
  dimnames = list(c("face", "body", "house", "tool"), NULL))
origins <- list(id = "four-runs", partitions = list(a = "run1", b = "run2",
  c = "run3", d = "run4"), independence = "independent",
  assumption = "Independent runs with externally fixed preprocessing.")
rel <- relation(list(a = B, b = B, c = B, d = B), domain = d,
  provenance = list(observation_origins = origins))
at <- compile_frame(whole_brain(normalization = "none"), d)
train <- plan_geometry(rel, at, pairing("a", "b",
  independence = "independent", generalizes_over = "run"))
test <- plan_geometry(rel, at, pairing("c", "d",
  independence = "independent", generalizes_over = "run"))
basis <- model_basis(features = B, conditions = rel$effect_space, normalize = "trace")
fit <- fit_geometry(train, basis, rank = 2, penalty = 0.1)
evidence <- score_geometry(fit, test, modes = TRUE)
as.data.frame(evidence)
#>   measurement gain rank inner_product prediction_norm_sq
#> 1 whole_brain 7.92    2           7.2               6.48
as.data.frame(evidence, view = "groups")
#>   measurement group size amplitude_sum evidence_sum prediction_cost gain
#> 1 whole_brain     1    2           3.6            4            6.48 7.92
```
