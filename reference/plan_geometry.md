# Compile a reusable query-first geometry plan

`plan_geometry()` validates the relation, spatial frame, pairing, source
capabilities, and compute policy. The plan can answer fixed queries with
[`evaluate_geometry()`](https://bbuchsbaum.github.io/crossform/reference/evaluate_geometry.md)
or be explicitly materialized with
[`materialize_geometry()`](https://bbuchsbaum.github.io/crossform/reference/materialize_geometry.md).
Complete packed geometry is therefore an optional materialization, not
the object that every analysis must allocate.

## Usage

``` r
plan_geometry(
  x,
  at,
  over,
  compute = compute_policy(),
  metric = NULL,
  composition = c("native", "whitened"),
  right = NULL,
  training = metric_training_policy("exclude_evaluation"),
  residual_workspace_bytes = NULL
)
```

## Arguments

- x:

  An `effect_relation` supplying the left experimental axis, or an
  `effect_relation_fit` when `metric` is a learned recipe, which needs
  the fit's residual channel.

- at:

  A compiled additive `effect_frame`.

- over:

  An `effect_pairing`. Rectangular plans require directed pairings whose
  left endpoints identify partitions of `x` and right endpoints identify
  partitions of `right`.

- compute:

  A sequential
  [`compute_policy()`](https://bbuchsbaum.github.io/crossform/reference/compute_policy.md).

- metric:

  Optional fixed
  [`neural_metric()`](https://bbuchsbaum.github.io/crossform/reference/neural_metric.md),
  or an on-demand metric recipe such as
  [`shrinkage_precision()`](https://bbuchsbaum.github.io/crossform/reference/shrinkage_precision.md).
  A domain-wide fixed metric is restricted to each frame support on
  demand; a support-local fixed metric is accepted only for a matching
  one-node frame. A recipe compiles a learned local schedule instead:
  nothing is materialized before the frame, and each support derives its
  own local operator from frozen residual sufficient statistics. Fixed
  metrics and recipes are not yet admitted on rectangular plans.

- composition:

  How the metric composes with the frame, and therefore which estimand
  the plan names. `"native"` (the default, and the only behaviour before
  this argument existed) weights features and then measures them in the
  metric geometry, \\K_x = D(\sqrt{w_x})\\Q\\D(\sqrt{w_x})\\.
  `"whitened"` measures the frame in whitened coordinates instead, \\K_x
  = Q^{1/2}D(w_x)Q^{1/2}\\, which conserves exactly under a conservative
  frame for any positive-definite `Q` where the native composition of a
  dense `Q` does not. The two are **different estimands, not two
  implementations of one**: under whitening a node weights whitened
  coordinates, which are spatially delocalized whenever `Q` is dense, so
  a node's support is no longer its support. Never switch to
  `"whitened"` to repair a failed conservation check. It is admitted for
  a fixed positive-definite domain-wide metric only, uses the symmetric
  positive-definite root \\Q^{1/2}\\, and carries both the composition
  and that root convention into `$scientific_plan_id`. See
  `design/conservative-geometry-contract.md` section 5.

- right:

  Optional second `effect_relation` supplying a distinct right
  experimental axis. Supplying it compiles a rectangular cross-axis
  plan: the resulting form has one row axis per left effect and one
  column axis per right effect, is read with axis-bound
  [`pair_query()`](https://bbuchsbaum.github.io/crossform/reference/pair_query.md)s
  through
  [`evaluate_geometry()`](https://bbuchsbaum.github.io/crossform/reference/evaluate_geometry.md),
  and materializes to a rectangular effect form. Encoding-retrieval
  similarity is the canonical use.

- training:

  A
  [`metric_training_policy()`](https://bbuchsbaum.github.io/crossform/reference/metric_training_policy.md)
  declaring which residual partitions may train a learned `metric`.
  Ignored for fixed metrics.

- residual_workspace_bytes:

  Positive cache budget used while accumulating canonical residual pair
  sufficient statistics for a learned `metric`. It changes cache
  capacity, never the canonical numerical tile shape, and so enters
  `$signature` but not `$scientific_plan_id`. Defaults to
  `compute$workspace_bytes`, or 512 MiB when the policy declares none.
  Admitted only with a recipe.

## Value

An `effect_geometry_plan` recording the compiled `$task`, `$frame`,
`$pairing`, `$metric_schedule`, and `$compute` policy, the
`$logical_shape` and `$measurements` it will produce, and a
`$scientific_plan_id` naming the estimand. Changing block size or
storage changes the execution receipt, not this identity.

## Structure

A plan is a declaration, so every element describes what will be
computed rather than a result.

- `$frame`: the compiled `effect_frame` the geometry is measured at.

- `$pairing`: the `effect_pairing` whose rows are the partition products
  the plan may form.

- `$measurements`: how many spatial measurements every view will return,
  one per row of `$frame$weights`.

- `$logical_shape`: the effect-axis extent as `c(left, right)`. The two
  entries are equal on a self-form plan.

- `$task$left_relation`: the relation supplying the left experimental
  axis. Its `$effects` is the order unnamed contrast weights are read
  in.

- `$compute`: the
  [`compute_policy()`](https://bbuchsbaum.github.io/crossform/reference/compute_policy.md)
  the plan was compiled under.

- `$scientific_plan_id`: the estimand identity. Keep it with the
  analysis record; block size, storage, and machine do not change it.

The rest of `$task`, the `$metric_schedule`, and every other element not
listed here are the lowered form the executor consumes: internal, and
free to change.

## Reads at plan time

With a fixed metric or the implicit identity metric under the native
composition, compilation reads no relation block: everything it checks
is a declaration.

There are two deliberate exceptions. `composition = "whitened"` is the
second: whitening is a global congruence, so the transform \\\tilde B =
BQ^{1/2}\\ is performed once here and every view of the plan reads the
whitened coordinates it produced. Deferring it would repeat one
domain-wide congruence per contrast. The relation is read one feature
block at a time and accumulated, so the source is never fully resident,
but the result is: `$execution_hints` records the resident bytes, and a
declared `compute_policy(workspace_bytes = )` is enforced against them
before the first read. Every refusal that does not need the data — the
metric's domain, its support, its definiteness, the budget — fires
first.

A learned metric recipe is the other, and it is deliberate too. The
schedule freezes canonical residual sufficient statistics, so
`plan_geometry()` accumulates them in one streamed pass over the fit's
residual channel before it returns. Deferring that pass to execution
would re-accumulate the same statistics for every contrast and lose the
plan reuse the frozen schedule exists to provide. The refusals that do
not need the data still fire first: a training-partition shortage and a
workspace budget overflow are both diagnosed before the first residual
read. The accumulation is recorded in `$execution_hints`, so it is
visible in the plan signature rather than only in the execution receipt.

## See also

[`contrast_energy()`](https://bbuchsbaum.github.io/crossform/reference/contrast_energy.md),
[`rdm()`](https://bbuchsbaum.github.io/crossform/reference/rdm.md),
[`rsa()`](https://bbuchsbaum.github.io/crossform/reference/rsa.md), and
[`crossnobis()`](https://bbuchsbaum.github.io/crossform/reference/crossnobis.md)
for the named views of a plan;
[`evaluate_geometry()`](https://bbuchsbaum.github.io/crossform/reference/evaluate_geometry.md)
for an arbitrary fixed query and
[`materialize_geometry()`](https://bbuchsbaum.github.io/crossform/reference/materialize_geometry.md)
for complete packed geometry;
[`coupling()`](https://bbuchsbaum.github.io/crossform/reference/coupling.md)
for the adjoint neural-side closure.

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
[`plan_crossnobis()`](https://bbuchsbaum.github.io/crossform/reference/plan_crossnobis.md),
[`plot_views`](https://bbuchsbaum.github.io/crossform/reference/plot_views.md),
[`query_geometry()`](https://bbuchsbaum.github.io/crossform/reference/query_geometry.md),
[`rdm()`](https://bbuchsbaum.github.io/crossform/reference/rdm.md),
[`rsa()`](https://bbuchsbaum.github.io/crossform/reference/rsa.md),
[`variation_query()`](https://bbuchsbaum.github.io/crossform/reference/variation_query.md)

## Examples

``` r
domain <- abstract_domain(3, id = "plan-example")
run1 <- rbind(a = c(1, 0, 2), b = c(0, 1, 1))
run2 <- rbind(a = c(1.1, 0.1, 1.9), b = c(0.1, 0.9, 1.2))
relation <- relation(list(run1 = run1, run2 = run2), domain = domain)

# The plan validates relation, frame, pairing, and metric agreement
# without reading any neural value.
plan <- plan_geometry(
  relation,
  compile_frame(voxelwise(), domain),
  cross_partitions(relation, independence = "independent")
)
plan
#> <effect_geometry_plan>
#>   effects:      2 x 2
#>   measurements: 3
#>   features:     3
#>   metric:       implicit identity
#>   generalizes:  1 partition pairs (axis undeclared), endpoints independent
#>   execution:    query-first, in memory
#>   state:        nothing computed yet
#>   next:         contrast_energy(plan, weights), rdm(plan), rsa(plan)

# One plan answers many fixed queries. Save its identity with the
# analysis record: it names the estimand, not the execution.
plan$scientific_plan_id
#> [1] "geometry-sha256:db28503f013a522a76efba2c7aeb52c8c4062e7d677e5b60e6ab5f95a3e4fef4"
result <- evaluate_geometry(plan, query = bilinear_query(diag(2)))
as.data.frame(result)
#>   measurement view1
#> 1           1   1.1
#> 2           2   0.9
#> 3           3   5.0

# The usual next step is a named view of the same plan.
contrast_energy(plan, c(a = 1, b = -1))$total
#> [1] 1.0 0.8 0.7
```
