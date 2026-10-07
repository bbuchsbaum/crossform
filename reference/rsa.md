# Fit multiple-regression RSA as one compiled geometry query

`rsa()` compiles the RDM transform and the least-squares coefficient map
into a single fixed query, so the regression is executed as one pass
over the plan rather than as a second analysis of a stored RDM. Model
rows and columns are aligned to the relation's effect names before any
geometry is read.

## Usage

``` r
rsa(x, ...)

# Default S3 method
rsa(
  x,
  models,
  nuisance = NULL,
  intercept = TRUE,
  component = c("total", "coherent", "configuration"),
  ...
)
```

## Arguments

- x:

  An `effect_geometry_plan` or a complete effect form carrying the
  symmetric self-form capability.

- ...:

  Passed to the method. The generic dispatches on `x`: a geometry plan
  or a complete effect form is read by the default method documented
  here, and an `effect_population_result` from
  [`estimate_population()`](https://bbuchsbaum.github.io/crossform/reference/estimate_population.md)
  by the group-level method in
  [population_views](https://bbuchsbaum.github.io/crossform/reference/population_views.md).

- models:

  Named model RDMs, each a finite symmetric zero-diagonal matrix over
  the experimental effects. Row and column names are optional; when
  supplied they must identify every effect exactly once and are
  reordered to the relation's effect order.

- nuisance:

  Optional named nuisance RDMs, in the same form.

- intercept:

  Whether to include an intercept in RDM space. It is `TRUE` by default,
  so a model RDM that is constant off the diagonal is collinear with it;
  the rank-deficiency message names the columns involved and
  `intercept = FALSE` fits the same models without the constant column.

- component:

  Geometry component to read.

## Value

An `effect_rsa_view`. `$coefficients` has one row per measurement and
one named column per model, nuisance model, and the optional intercept;
`$component`, `$index`, and `$receipt` record what was read.

## Structure

The fit is one measurement-by-term coefficient matrix; the other
elements name its axes and record what was read.

- `$coefficients`: one row per spatial measurement, one named column per
  design term, in `$terms` row order.

- `$terms`: a data frame naming each column of `$coefficients` in `term`
  and labeling it `intercept`, `model`, or `nuisance` in `role`.

- `$component`: the geometry component the models were regressed on.

- `$index`: the measurement identifiers, one per row of `$coefficients`,
  carried from the frame's `$index$measurement`.

- `$receipt`: the execution receipt for the run that produced the fit.

The compiled `$query` and any other element not listed here are internal
and may change.

## Refusal

A rectangular cross-axis plan or a non-symmetric form signals an
`effect_capability_refusal` with capability `"symmetric_self_form"` in
namespace `"geometry_views"`; see
[`catch_refusal()`](https://bbuchsbaum.github.io/crossform/reference/catch_refusal.md).

## See also

[`rdm()`](https://bbuchsbaum.github.io/crossform/reference/rdm.md) for
the distances the regression is fitted to, and
[`plan_geometry()`](https://bbuchsbaum.github.io/crossform/reference/plan_geometry.md)
for the plan.

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
[`plan_geometry()`](https://bbuchsbaum.github.io/crossform/reference/plan_geometry.md),
[`plot_views`](https://bbuchsbaum.github.io/crossform/reference/plot_views.md),
[`query_geometry()`](https://bbuchsbaum.github.io/crossform/reference/query_geometry.md),
[`rdm()`](https://bbuchsbaum.github.io/crossform/reference/rdm.md),
[`variation_query()`](https://bbuchsbaum.github.io/crossform/reference/variation_query.md)

## Examples

``` r
domain <- abstract_domain(4, id = "rsa-example")
run1 <- rbind(
  face = c(1, 0.2, 0, 0), house = c(0, 1, 0.1, 0), tool = c(0, 0, 1, 0.3)
)
run2 <- rbind(
  face = c(0.9, 0.3, 0, 0), house = c(0.1, 0.9, 0, 0), tool = c(0, 0.1, 1.1, 0.2)
)
relation <- relation(list(run1 = run1, run2 = run2), domain = domain)
plan <- plan_geometry(
  relation, compile_frame(regions(c("v1", "v1", "it", "it")), domain),
  cross_partitions(relation, independence = "independent")
)

# An animacy model: the two animate-inanimate pairs are far, the rest near.
conditions <- rownames(run1)
animacy <- matrix(
  c(0, 0, 1, 0, 0, 1, 1, 1, 0), 3, 3,
  dimnames = list(conditions, conditions)
)
fit <- rsa(plan, models = list(animacy = animacy))
round(fit$coefficients, 3)
#>      (Intercept) animacy
#> [1,]        0.64  -0.205
#> [2,]        0.00   0.553
as.data.frame(fit)
#>   measurement  (Intercept) animacy
#> 1          v1 6.400000e-01 -0.2050
#> 2          it 3.525432e-18  0.5525

# A model must be a dissimilarity matrix. Passing a similarity matrix,
# whose diagonal is nonzero, is rejected before geometry is read.
similarity <- 1 - animacy
wrong <- try(rsa(plan, models = list(animacy = similarity)), silent = TRUE)
conditionMessage(attr(wrong, "condition"))
#> [1] "`models` RDM `animacy` has a nonzero diagonal (largest |m[i, i]| is 1). Pass a dissimilarity matrix, not a similarity matrix: an effect is at distance zero from itself."
```
