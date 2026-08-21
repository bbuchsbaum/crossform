# Declare the first-moment observation model

States the estimation assumptions the data tables cannot prove: the
error law, the sampling unit, any whitener, and whether independence was
asserted. Nothing here is inferred from events or confounds, and the
declaration is what determines whether analytic uncertainty is later
available.

## Usage

``` r
observation_model(
  kind = c("ols", "fixed_gls", "learned_frozen_gls"),
  sampling_unit,
  whitener = NULL,
  independence = NULL,
  training_revision = NULL,
  training_provenance = list(),
  assumptions = list(),
  provenance = list()
)
```

## Arguments

- kind:

  `"ols"`, `"fixed_gls"`, or `"learned_frozen_gls"`.

- sampling_unit:

  One sampling-unit axis name, or a named value per partition. It is
  required and is never inferred from events.

- whitener:

  For GLS, one finite square matrix or a named list per partition.

- independence:

  Optional explicit independence declaration. `NULL` records that no
  independence guarantee was asserted.

- training_revision:

  Strong training revision required for learned GLS.

- training_provenance:

  Portable learned-model provenance.

- assumptions:

  Additional portable declared assumptions.

- provenance:

  Portable model provenance.

## Value

An `effect_observation_model`: a list with `$kind`, `$sampling_unit`,
the `$whitener` and its `$whitener_revisions`, `$independence`,
learned-model `$training_revision` and `$training_provenance`,
`$assumptions`, `$provenance`, an `$observation_model_id`, and
`$capabilities` with `fixed_observation_model`,
`learned_observation_model`, `separable_glm_law`,
`analytic_effect_covariance`, and `declared_independence`.

## See also

[`plan_relation()`](https://bbuchsbaum.github.io/crossform/reference/plan_relation.md),
which consumes this declaration, and
[`rdm_sampling_covariance()`](https://bbuchsbaum.github.io/crossform/reference/rdm_sampling_covariance.md),
whose analytic route requires `analytic_effect_covariance`.

Other studies and effect maps:
[`coefficient_parameterization()`](https://bbuchsbaum.github.io/crossform/reference/coefficient_parameterization.md),
[`condition_space()`](https://bbuchsbaum.github.io/crossform/reference/condition_space.md),
[`design_model()`](https://bbuchsbaum.github.io/crossform/reference/design_model.md),
[`effect_map()`](https://bbuchsbaum.github.io/crossform/reference/effect_map.md),
[`effect_space()`](https://bbuchsbaum.github.io/crossform/reference/effect_space.md),
[`lower_effect_map()`](https://bbuchsbaum.github.io/crossform/reference/lower_effect_map.md),
[`raw_design_model()`](https://bbuchsbaum.github.io/crossform/reference/raw_design_model.md),
[`raw_effect_map()`](https://bbuchsbaum.github.io/crossform/reference/raw_effect_map.md),
[`study()`](https://bbuchsbaum.github.io/crossform/reference/study.md),
[`study_axis()`](https://bbuchsbaum.github.io/crossform/reference/study_axis.md),
[`study_capabilities()`](https://bbuchsbaum.github.io/crossform/reference/study_capabilities.md)

## Examples

``` r
# OLS with scans as the sampling unit. Independence must be declared; it is
# never read off the event table.
ols <- observation_model(
  "ols", sampling_unit = "scan",
  independence = "runs independently acquired conditional on the model"
)
ols$capabilities$analytic_effect_covariance
#> [1] TRUE
ols$capabilities$declared_independence
#> [1] TRUE

# A frozen learned whitener still yields point estimates, but the analytic
# error channel is withheld because its training is not accounted for.
learned <- observation_model(
  "learned_frozen_gls", sampling_unit = "scan", whitener = diag(6),
  training_revision = paste0("sha256:", strrep("a", 64)),
  training_provenance = list(method = "AR1", fitted_on = "held-out runs")
)
learned$capabilities$analytic_effect_covariance
#> [1] FALSE
```
