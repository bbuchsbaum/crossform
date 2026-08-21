# Take the adjoint-side coupling closure of a geometry plan

A geometry plan closes the neural boundary and leaves the experimental
axes open; `coupling()` asks the adjoint question from the same plan
vocabulary: close the experimental boundary with a declared query and
leave two neural measurements open. The plan's spatial frame supplies
the measurement legs and the plan's pairing supplies the partition
products, so no second frame or pairing declaration is required.

## Usage

``` r
coupling(
  x,
  between,
  by,
  over = NULL,
  mode = c("total", "coherent", "coherent_configuration"),
  reducer = aggregate_first(),
  compute = compute_policy(),
  route = "auto"
)
```

## Arguments

- x:

  A self-form `effect_geometry_plan`.

- between:

  A two-column matrix or data frame of frame node pairs, given as node
  names from the plan's frame index or as node indices. Each row becomes
  one ordered measurement edge.

- by:

  The experimental closure: a
  [`variation_query()`](https://bbuchsbaum.github.io/crossform/reference/variation_query.md).

- over:

  Optional `effect_pairing`; defaults to the plan's own pairing.

- mode:

  Measurement mode passed to
  [`measurement_frame()`](https://bbuchsbaum.github.io/crossform/reference/measurement_frame.md).

- reducer:

  Partition reducer for the measurement pipeline.

- compute:

  A
  [`compute_policy()`](https://bbuchsbaum.github.io/crossform/reference/compute_policy.md).

- route:

  Measurement execution route; `"auto"` selects one.

## Value

An `effect_measurement_form` over the requested node pairs, with one
block per edge in `$values` once read. Read it with
[`effect_coupling()`](https://bbuchsbaum.github.io/crossform/reference/coupling_views.md),
[`connectivity()`](https://bbuchsbaum.github.io/crossform/reference/connectivity.md),
[`measurement_components()`](https://bbuchsbaum.github.io/crossform/reference/measurement_components.md),
or the other coupling views.

## Details

The measurement pipeline this compiles into is deliberately small-node:
frames whose dense payload exceeds the admitted ceiling refuse before
reading data, so searchlight-resolved coupling remains explicitly
unadmitted rather than silently truncated.

## Refusal

A rectangular cross-axis plan signals an `effect_capability_refusal`
with capability `"self_form_coupling"` in namespace `"coupling_plans"`;
branch on it with
[`catch_refusal()`](https://bbuchsbaum.github.io/crossform/reference/catch_refusal.md).

## See also

[`plan_geometry()`](https://bbuchsbaum.github.io/crossform/reference/plan_geometry.md)
for the plan this closes the other way,
[`measurement_form()`](https://bbuchsbaum.github.io/crossform/reference/measurement_form.md)
for the general constructor, and
[`variation_query()`](https://bbuchsbaum.github.io/crossform/reference/variation_query.md)
for the experimental closure.

Other coupling and connectivity views:
[`connectivity()`](https://bbuchsbaum.github.io/crossform/reference/connectivity.md),
[`control_coupling()`](https://bbuchsbaum.github.io/crossform/reference/control_coupling.md),
[`coupling_contrast()`](https://bbuchsbaum.github.io/crossform/reference/coupling_contrast.md),
[`coupling_views`](https://bbuchsbaum.github.io/crossform/reference/coupling_views.md),
[`effect_coupling_result`](https://bbuchsbaum.github.io/crossform/reference/effect_coupling_result.md),
[`match_control()`](https://bbuchsbaum.github.io/crossform/reference/match_control.md),
[`match_coupling()`](https://bbuchsbaum.github.io/crossform/reference/match_coupling.md),
[`measurement_components()`](https://bbuchsbaum.github.io/crossform/reference/measurement_components.md),
[`measurement_form()`](https://bbuchsbaum.github.io/crossform/reference/measurement_form.md),
[`pair_lm_query()`](https://bbuchsbaum.github.io/crossform/reference/pair_lm_query.md),
[`pair_query()`](https://bbuchsbaum.github.io/crossform/reference/pair_query.md)

## Examples

``` r
# A plan already fixes the relation, the spatial frame, and the pairing.
# `coupling()` reuses all three and asks the adjoint question: keep two
# neural measurements open and close the experimental axis instead.
set.seed(4)
domain <- abstract_domain(4, id = "coupling-plan-example")
relation <- relation(
  list(
    run1 = matrix(rnorm(8), 2, 4, dimnames = list(c("a", "b"), NULL)),
    run2 = matrix(rnorm(8), 2, 4, dimnames = list(c("a", "b"), NULL))
  ),
  effects = effect_space(c("a", "b")), domain = domain
)
plan <- plan_geometry(
  relation, compile_frame(regions(c("r1", "r1", "r2", "r2")), domain),
  cross_partitions(relation, independence = "independent")
)

# `between` names frame nodes, so no second frame or pairing is declared.
form <- coupling(
  plan, cbind("r1", "r2"),
  by = variation_query(
    diag(2) - 0.5, relation$effect_space, "trial", "psd_variation"
  )
)
effect_coupling(form)$values[[1]]
#>            [,1]        [,2]
#> [1,] 0.05479508 -0.02875980
#> [2,] 0.08096420 -0.07851908

# Node names are resolved against the plan's own frame index.
unknown <- try(coupling(plan, cbind("r1", "r9"), by = variation_query(
  diag(2) - 0.5, relation$effect_space, "trial", "psd_variation"
)), silent = TRUE)
conditionMessage(attr(unknown, "condition"))
#> [1] "`between` to names `r9`, which is not a measurement of the plan's frame. The frame has 2 measurements: `r1`, `r2`."
```
