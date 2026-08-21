# Compile an on-demand learned-metric crossnobis plan

The plan freezes the residual-statistics identity, training policy,
regularization recipe, evaluation pairing, and spatial support graph. It
does not retain one covariance or precision matrix per node and fold;
each local solve handle is derived on demand from canonical pair
statistics. The learned metric is itself random. Version 0.1 therefore
returns signed crossnobis point estimates but does not apply the
fixed-metric analytic covariance law in
[`rdm_sampling_covariance()`](https://bbuchsbaum.github.io/crossform/reference/rdm_sampling_covariance.md)
or claim calibrated intervals for this plan; doing so requires
propagation of metric-estimation uncertainty (for example, an admitted
LD-t specialization).

## Usage

``` r
plan_crossnobis(
  x,
  at,
  over,
  metric = shrinkage_precision(),
  training = metric_training_policy("exclude_evaluation"),
  compute = compute_policy(),
  residual_workspace_bytes = NULL
)
```

## Arguments

- x:

  An `effect_relation_fit` with residual-block capability.

- at:

  A support-index-backed compiled spatial frame.

- over:

  Independent cross-partition evaluation edges.

- metric:

  An on-demand metric recipe such as
  [`shrinkage_precision()`](https://bbuchsbaum.github.io/crossform/reference/shrinkage_precision.md).

- training:

  A
  [`metric_training_policy()`](https://bbuchsbaum.github.io/crossform/reference/metric_training_policy.md).

- compute:

  A sequential
  [`compute_policy()`](https://bbuchsbaum.github.io/crossform/reference/compute_policy.md).

- residual_workspace_bytes:

  Positive budget used while accumulating canonical residual pair
  sufficient statistics. It changes cache capacity, never the canonical
  numerical tile shape. Defaults to `compute$workspace_bytes`, or 512
  MiB when the policy declares none.

## Value

An `effect_geometry_plan` reusable across fixed contrasts. It carries
the learned `$metric_schedule` wrapping the frozen schedule, the
`$frame`, `$pairing`, and `$execution_hints`, and a
`$scientific_plan_id` that identifies the estimand independently of
execution choices.

## Details

`plan_crossnobis()` is a convenience wrapper around
[`plan_geometry()`](https://bbuchsbaum.github.io/crossform/reference/plan_geometry.md),
not a deprecated alias: it names the intent, and compiles exactly
`plan_geometry(x, at, over, metric = <recipe>, training = )`. The plan
it returns is an ordinary `effect_geometry_plan` whose metric schedule
has kind `learned_local_before_frame`, lowered by the geometry compiler
like any other. The pairing contract crossnobis requires – independent,
cross-partition, self-product-free – is enforced by
[`crossnobis()`](https://bbuchsbaum.github.io/crossform/reference/crossnobis.md)
when the plan is read, which is where the fixed-metric route enforces
it.

## See also

[`plan_geometry()`](https://bbuchsbaum.github.io/crossform/reference/plan_geometry.md),
which this wraps;
[`crossnobis()`](https://bbuchsbaum.github.io/crossform/reference/crossnobis.md)
to read a contrast from this plan,
[`shrinkage_precision()`](https://bbuchsbaum.github.io/crossform/reference/shrinkage_precision.md)
and
[`metric_training_policy()`](https://bbuchsbaum.github.io/crossform/reference/metric_training_policy.md)
for the metric declarations it freezes, and
[`residual_pair_statistics()`](https://bbuchsbaum.github.io/crossform/reference/residual_pair_statistics.md)
for the sufficient statistics it compiles.

Other geometry plans and views:
[`aggregate_first()`](https://bbuchsbaum.github.io/crossform/reference/aggregate_first.md),
[`bilinear_query()`](https://bbuchsbaum.github.io/crossform/reference/bilinear_query.md),
[`coherence_spectrum()`](https://bbuchsbaum.github.io/crossform/reference/coherence_spectrum.md),
[`compute_policy()`](https://bbuchsbaum.github.io/crossform/reference/compute_policy.md),
[`contrast_energy()`](https://bbuchsbaum.github.io/crossform/reference/contrast_energy.md),
[`contribution()`](https://bbuchsbaum.github.io/crossform/reference/contribution.md),
[`crossnobis()`](https://bbuchsbaum.github.io/crossform/reference/crossnobis.md),
[`evaluate_geometry()`](https://bbuchsbaum.github.io/crossform/reference/evaluate_geometry.md),
[`example_fmri_effects()`](https://bbuchsbaum.github.io/crossform/reference/example_fmri_effects.md),
[`geometry_component()`](https://bbuchsbaum.github.io/crossform/reference/geometry_component.md),
[`geometry_spectrum()`](https://bbuchsbaum.github.io/crossform/reference/geometry_spectrum.md),
[`latent_geometry()`](https://bbuchsbaum.github.io/crossform/reference/latent_geometry.md),
[`materialize_geometry()`](https://bbuchsbaum.github.io/crossform/reference/materialize_geometry.md),
[`plan_geometry()`](https://bbuchsbaum.github.io/crossform/reference/plan_geometry.md),
[`plot_views`](https://bbuchsbaum.github.io/crossform/reference/plot_views.md),
[`query_geometry()`](https://bbuchsbaum.github.io/crossform/reference/query_geometry.md),
[`rdm()`](https://bbuchsbaum.github.io/crossform/reference/rdm.md),
[`rsa()`](https://bbuchsbaum.github.io/crossform/reference/rsa.md),
[`variation_query()`](https://bbuchsbaum.github.io/crossform/reference/variation_query.md)

## Examples

``` r
set.seed(7)
domain <- abstract_domain(
  3, coordinates = cbind(x = 0:2, y = 0),
  id = "learned-crossnobis-example"
)
design <- cbind(intercept = 1, condition = rep(c(-0.5, 0.5), 4))
coefficients <- rbind(
  intercept = c(0, 0, 0),
  condition = c(0.6, -0.4, 0.2)
)
responses <- setNames(lapply(seq_len(3), function(run) {
  design %*% coefficients + matrix(rnorm(8 * 3, sd = 0.25), 8, 3)
}), paste0("run", seq_len(3)))
fit <- lm_relation_fit(
  responses, design, rbind(condition = c(0, 1)), domain = domain
)
plan <- plan_crossnobis(
  fit,
  compile_frame(searchlights(1.01), domain),
  pairing("run1", "run2", independence = "independent"),
  metric = shrinkage_precision(0.2)
)
result <- crossnobis(plan, c(condition = 1))
result
#> <effect_crossnobis_view>
#>   measurements: 3
#>   contrast:     condition 1
#>  measurement crossnobis
#>            1     10.255
#>            2      6.824
#>            3      1.976
as.data.frame(result)
#>   measurement crossnobis
#> 1           1  10.255026
#> 2           2   6.823935
#> 3           3   1.976304
```
