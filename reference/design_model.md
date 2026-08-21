# Declare a semantic design model with compiled routes

`design_model_id` covers the semantic mean-model request. Concrete
design matrices, coding maps, row order, compiler build, and solver
route are kept in `compilation_route_id` and later design receipts.

## Usage

``` r
design_model(
  specification,
  conditions,
  designs,
  parameterizations,
  row_ids = NULL,
  solver = "auto",
  protocol = "crossform-semantic-design",
  protocol_version = "1",
  package = "crossform",
  package_version = "0.0.0.9000",
  provenance = list()
)
```

## Arguments

- specification:

  A portable semantic model declaration.

- conditions:

  The bound
  [`condition_space()`](https://bbuchsbaum.github.io/crossform/reference/condition_space.md).

- designs:

  Named observation-by-coefficient matrices, one per partition.

- parameterizations:

  Named
  [`coefficient_parameterization()`](https://bbuchsbaum.github.io/crossform/reference/coefficient_parameterization.md)
  values for the designs.

- row_ids:

  Ordered observation identifiers per design. Defaults to design row
  names; positional designs are refused.

- solver:

  Numerical route per partition.

- protocol, protocol_version:

  Semantic compiler protocol identity.

- package, package_version:

  Compiler implementation receipt fields.

- provenance:

  Portable semantic model provenance.

## Value

An `effect_design_model`: a list with the bound `$condition_space`, the
`$specification`, `$partitions`, the compiled `$designs`, their
`$parameterizations`, `$row_ids`, the `$compiler` record, the `$solver`
route, `$provenance`, a `$design_model_id` covering the semantic request
only, a `$compilation_route_id` covering the concrete compilation, and
`$capabilities` with `symbolic_model`, `coding_invariant`,
`row_lineage`.

## See also

[`coefficient_parameterization()`](https://bbuchsbaum.github.io/crossform/reference/coefficient_parameterization.md)
for the coding it carries,
[`raw_design_model()`](https://bbuchsbaum.github.io/crossform/reference/raw_design_model.md)
for the route without semantic coding, and
[`plan_relation()`](https://bbuchsbaum.github.io/crossform/reference/plan_relation.md),
which binds this model to a study.

Other studies and effect maps:
[`coefficient_parameterization()`](https://bbuchsbaum.github.io/crossform/reference/coefficient_parameterization.md),
[`condition_space()`](https://bbuchsbaum.github.io/crossform/reference/condition_space.md),
[`effect_map()`](https://bbuchsbaum.github.io/crossform/reference/effect_map.md),
[`effect_space()`](https://bbuchsbaum.github.io/crossform/reference/effect_space.md),
[`lower_effect_map()`](https://bbuchsbaum.github.io/crossform/reference/lower_effect_map.md),
[`observation_model()`](https://bbuchsbaum.github.io/crossform/reference/observation_model.md),
[`raw_design_model()`](https://bbuchsbaum.github.io/crossform/reference/raw_design_model.md),
[`raw_effect_map()`](https://bbuchsbaum.github.io/crossform/reference/raw_effect_map.md),
[`study()`](https://bbuchsbaum.github.io/crossform/reference/study.md),
[`study_axis()`](https://bbuchsbaum.github.io/crossform/reference/study_axis.md),
[`study_capabilities()`](https://bbuchsbaum.github.io/crossform/reference/study_capabilities.md)

## Examples

``` r
conditions <- condition_space(c("face", "body"), basis_id = "cond-mean:v1")
design <- cbind(
  face = c(1, 0, 1, 0, 1, 0), body = c(0, 1, 0, 1, 0, 1),
  drift = seq(-1, 1, length.out = 6L)
)
rownames(design) <- paste0("scan-", 1:6)
map <- cbind(diag(2), drift = 0)
dimnames(map) <- list(conditions$coordinates, colnames(design))
coding <- coefficient_parameterization(
  map, conditions, coding_id = "cell-means-plus-drift"
)
specification <- list(target = "condition means", nuisance = "linear drift")

model <- design_model(
  specification, conditions,
  designs = list(`run-1` = design), parameterizations = list(`run-1` = coding)
)
model$capabilities$coding_invariant
#> [1] TRUE

# Semantic identity ignores the numerical route: switching solvers changes
# the compilation receipt, not what was scientifically requested.
route <- design_model(
  specification, conditions,
  designs = list(`run-1` = design), parameterizations = list(`run-1` = coding),
  solver = "svd"
)
c(same_request = identical(model$design_model_id, route$design_model_id),
  same_route = identical(
    model$compilation_route_id, route$compilation_route_id
  ))
#> same_request   same_route 
#>         TRUE        FALSE 
```
