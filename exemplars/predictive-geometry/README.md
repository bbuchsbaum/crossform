# Predictive geometry on independent runs

Fit model-supported neural amplitudes on two runs, freeze the prediction,
and score it on two independent runs. Run the example with an installed
development version of crossform:

```sh
Rscript exemplars/predictive-geometry/fixed-split.R
```

The example needs no downloaded data. Its deterministic partitions make the
arithmetic visible. The origin manifest records a sampling assumption; equal
numeric values neither prove nor disprove independence.

For a step-by-step version, read `vignette("predictive-geometry")` or the
[predictive geometry guide](https://bbuchsbaum.github.io/crossform/articles/predictive-geometry.html).

```r
basis <- model_basis(kernels = list(model = K),
  conditions = rel$effect_space, normalize = "trace")
fit <- fit_geometry(train_plan, basis, rank = 2, penalty = 0.1)
evidence <- score_geometry(fit, test_plan, modes = TRUE)
as.data.frame(evidence)
as.data.frame(evidence, view = "modes")
```

Model features and explicitly squared Euclidean distances give the same
answer as the kernel. The kernel retains its directional preferences through
the inverse-kernel penalty, after retained-trace normalization.

| Mode | Fitted amplitude | Signed test evidence | Prediction cost | Gain |
|---|---:|---:|---:|---:|
| 1 | 2 | 2 | 4 | 4 |
| 2 | 1 | 0 | 1 | -1 |

Each gain is `2 * amplitude * evidence - amplitude^2`. The rank-two
prediction's total gain is **3 squared geometry units**. The zero predictor
has gain zero. An independently declared rank-one predictor has gain four
on this fixture; the script reports it without selecting rank on these test
outcomes. In an analysis, choose rank, penalty and model weights before final
evaluation or through a separate inner-validation schedule.

Positive gain estimates improvement in squared geometry prediction error
over zero, conditional on training and under the declared sampling
assumptions. It is not a fraction explained. The source geometry and test
evidence remain signed; the learned PSD prediction is biased. Fully retained
tied eigenspaces have invariant grouped evidence available with
`as.data.frame(evidence, view = "groups")`.

The fit records actual training origins. Scoring on the training runs refuses,
including when those observations have been copied or renamed. This example
tests independent runs on the **same conditions**, not prediction for new
conditions. The Frobenius spectral solver also refuses an RDM/GLS loss.

For many measurements, `storage = "block"` writes factors or optional mode
evidence to a new directory and saves a completed `fit.rds` or `score.rds`.
Use `readRDS()` to reopen that record while its directory remains in place.
`row_block` and `compute_policy(workspace_bytes = ...)` bound owned workspace.
The existing complete fixed-neural-metric route admits memory output; block
storage currently uses the implicit identity metric.
