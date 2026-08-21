# Construct the version 0.1 compute policy

Version 0.1 deliberately owns no process pool and accepts exactly one
worker. Participant-level parallelism belongs outside the core geometry
call until a later executor passes its memory and determinism gates.

## Usage

``` r
compute_policy(workers = 1L, block_features = NULL, workspace_bytes = NULL)
```

## Arguments

- workers:

  Number of R workers. Must be exactly one in version 0.1.

- block_features:

  Optional positive feature-block size.

- workspace_bytes:

  Optional positive budget for crossform-owned live workspace. Baseline
  and total process RSS are not charged to this budget.

## Value

An `effect_compute_policy` recording `$workers`, `$block_features`,
`$workspace_bytes`, and the fixed `$process_backend`. It is a
declaration consumed by
[`plan_geometry()`](https://bbuchsbaum.github.io/crossform/reference/plan_geometry.md)
and friends; it never starts a worker itself.

## Refusal

Requesting more than one worker signals an `effect_capability_refusal`
carrying capability `"parallel_execution"` in namespace
`"compute_policy"`, with reason `"worker_pool_not_implemented"`. Branch
on it with
[`catch_refusal()`](https://bbuchsbaum.github.io/crossform/reference/catch_refusal.md)
rather than on the message text.

## See also

[`plan_geometry()`](https://bbuchsbaum.github.io/crossform/reference/plan_geometry.md)
and
[`measurement_form()`](https://bbuchsbaum.github.io/crossform/reference/measurement_form.md),
which accept this policy;
[`catch_refusal()`](https://bbuchsbaum.github.io/crossform/reference/catch_refusal.md)
for inspecting the refusal below.

Other geometry plans and views:
[`aggregate_first()`](https://bbuchsbaum.github.io/crossform/reference/aggregate_first.md),
[`bilinear_query()`](https://bbuchsbaum.github.io/crossform/reference/bilinear_query.md),
[`coherence_spectrum()`](https://bbuchsbaum.github.io/crossform/reference/coherence_spectrum.md),
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
[`plan_geometry()`](https://bbuchsbaum.github.io/crossform/reference/plan_geometry.md),
[`plot_views`](https://bbuchsbaum.github.io/crossform/reference/plot_views.md),
[`query_geometry()`](https://bbuchsbaum.github.io/crossform/reference/query_geometry.md),
[`rdm()`](https://bbuchsbaum.github.io/crossform/reference/rdm.md),
[`rsa()`](https://bbuchsbaum.github.io/crossform/reference/rsa.md),
[`variation_query()`](https://bbuchsbaum.github.io/crossform/reference/variation_query.md)

## Examples

``` r
# The default: one worker, no declared block or workspace ceiling.
default <- compute_policy()
c(workers = default$workers, backend = default$process_backend)
#>      workers      backend 
#>          "1" "sequential" 

# Bound the feature block and the crossform-owned workspace. Neither
# choice changes the estimand, only the execution receipt.
policy <- compute_policy(block_features = 64L, workspace_bytes = 64 * 1024^2)
policy$block_features
#> [1] 64

# Version 0.1 owns no process pool, so more than one worker is refused —
# as a classed capability refusal, not a message to match on.
refusal <- catch_refusal(compute_policy(workers = 4L))
refusal$capability
#> [1] "parallel_execution"
refusal$reasons
#> [1] "worker_pool_not_implemented"
refusal$remedies
#> [1] "Pass `workers = 1` and parallelize across participants outside the geometry call, or bound the work with `compute_policy(block_features = , workspace_bytes = )`."
```
