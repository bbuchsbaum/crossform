# Query an exact factorized sampling-covariance form

Applies exact linear covariance operations without materializing the
full covariance unless `operation = "materialize"` is requested
explicitly. A returned diagonal is a vector of variances, not standard
errors or an automatic confidence interval.

## Usage

``` r
sampling_covariance(
  x,
  operation = c("diagonal", "selected_entries", "apply", "quadratic_form", "transport",
    "materialize"),
  query = NULL,
  max_bytes = 512 * 1024^2,
  queries = NULL,
  scope = c("measurement", "cross_measurement")
)
```

## Arguments

- x:

  An object from
  [`rdm_sampling_covariance()`](https://bbuchsbaum.github.io/crossform/reference/rdm_sampling_covariance.md),
  or a batch of those objects. A batch returns one result per
  measurement.

- operation:

  One exact covariance operation. With `queries` supplied and
  `operation` omitted, the query-bank covariance form itself is returned
  instead of an operation on it.

- query:

  Operation-specific argument: a two-column integer matrix for
  `selected_entries`, an evidence-coordinate vector for `apply` or
  `quadratic_form`, or an output-by-evidence matrix for `transport`.

- max_bytes:

  Positive payload/workspace budget used only for explicit dense
  materialization.

- queries:

  Optional bank of contrast-energy queries: a `K`-by-`q` matrix of
  contrast rows, a list of contrast vectors, or one contrast vector.
  Each contrast must sum to zero. See the section below.

- scope:

  Whether a query bank is read at one measurement (`"measurement"`, the
  only admitted scope) or jointly across measurements
  (`"cross_measurement"`, which refuses).

## Value

A named variance vector (`"diagonal"`), selected covariance vector,
covariance action, scalar quadratic form, transported covariance, or
explicitly materialized covariance matrix, according to `operation`.
Names come from the distance labels carried by `x`, or from the query
labels once `queries` names a bank. With `queries` and no `operation`,
an `effect_sampling_covariance` with `basis = "query_bank"`, or a batch
of them named by measurement.

## Details

The values returned inherit the calibration target chosen when `x` was
built. Under `target = "plugin"` they carry the documented upward bias
of the partition-mean plug-in policy; see the `target` argument of
[`rdm_sampling_covariance()`](https://bbuchsbaum.github.io/crossform/reference/rdm_sampling_covariance.md).
Under `target = "null"` the law is exact on the variance scale.

In both cases the residual covariance behind the law is a plug-in with
\\\nu\\ degrees of freedom, so its *quadratic* contribution carries the
Wishart finite-sample correction described under
[`rdm_sampling_covariance()`](https://bbuchsbaum.github.io/crossform/reference/rdm_sampling_covariance.md).
Read `x$source$residual_df` and `x$source$residual_effective_dimension`
to see \\\nu\\ and \\P\_{\mathrm{eff}}\\ for the measurement you
queried; `print(x)` shows both.

## Banks of contrast-energy queries

`queries` names a bank of \\K\\ contrasts and moves the covariance into
the coordinates of their signed contrast energies – the same estimand
[`contrast_energy()`](https://bbuchsbaum.github.io/crossform/reference/contrast_energy.md)
reports as `$total`. Supply a `K`-by-`q` matrix of contrast rows, a list
of contrast vectors, or one contrast vector; named weights are aligned
to the relation's effect order exactly as
[`contrast_energy()`](https://bbuchsbaum.github.io/crossform/reference/contrast_energy.md)
aligns them, and row or list names become the query labels.

`sampling_covariance(x, queries = )` returns an
`effect_sampling_covariance` with `basis = "query_bank"`: the same
factorized class, over `K` coordinates instead of the pairwise
distances, carrying the aligned bank on `$source$queries` and the
distance lowering it used on `$source$query_lowering`. Adding an
`operation` applies that operation in the new basis, so
`sampling_covariance(x, "materialize", queries = bank)` is the
`K`-by-`K` covariance of the bank at one measurement, and the same call
on a batch returns one `K`-by-`K` matrix per measurement, named by the
measurement it belongs to (`values[[measurement]]`).

This is exactly the `"transport"` of the distance-basis covariance
through the operator that lowers each query onto the distances, and it
is tested against that route rather than derived in parallel: a centred
contrast satisfies
\\cc^\top=-\sum\_{i\<j}c_ic_j(e_i-e_j)(e_i-e_j)^\top\\, so `c'Gc` is a
linear functional of the pairwise distances. A contrast whose weights do
not sum to zero is refused rather than approximated: the crossvalidated
distances are invariant to a shift \\G\to
G+a\mathbf{1}^\top+\mathbf{1}a^\top\\ of the geometry and an uncentred
energy is not, so the distance basis carries no information about it.
Every refusal of the underlying law is inherited unchanged: a learned
metric or a whitened plan is refused by
[`rdm_sampling_covariance()`](https://bbuchsbaum.github.io/crossform/reference/rdm_sampling_covariance.md)
before a bank can be attached to it.

## What a query bank does not give you

Each covariance is local to one measurement.
`scope = "cross_measurement"` refuses with capability
`"cross_node_sampling_covariance"` and names what is missing: an
explicit spatial model, being the cross-support residual second moment
and the frame overlap it induces. A conservative frame does not supply
one either. Conservation is a law about estimates, not about their
uncertainty, so overlapping rows are strongly positively correlated,
per-row variances do not add to the variance of a conserved budget, and
the per-measurement blocks returned here must not be assembled into a
joint covariance.

## See also

[`rdm_sampling_covariance()`](https://bbuchsbaum.github.io/crossform/reference/rdm_sampling_covariance.md)
to build `x`,
[`sampling_capabilities()`](https://bbuchsbaum.github.io/crossform/reference/sampling_capabilities.md)
to check the law is available first, and
[`contrast_energy()`](https://bbuchsbaum.github.io/crossform/reference/contrast_energy.md)
for the point estimates a query bank describes.

Other sampling uncertainty:
[`rdm_sampling_covariance()`](https://bbuchsbaum.github.io/crossform/reference/rdm_sampling_covariance.md),
[`reconstruct_evidence()`](https://bbuchsbaum.github.io/crossform/reference/reconstruct_evidence.md),
[`sampling_capabilities()`](https://bbuchsbaum.github.io/crossform/reference/sampling_capabilities.md)

## Examples

``` r
# Three conditions in three runs, with a fixed identity noise metric.
set.seed(11)
conditions <- c("face", "house", "tool")
design <- model.matrix(~ 0 + factor(rep(conditions, each = 3)))
colnames(design) <- conditions
effects <- diag(3)
rownames(effects) <- conditions
truth <- rbind(
  face = c(0.6, 0.1, 0), house = c(0, 0.5, 0.2), tool = c(0.1, 0, 0.6)
)
responses <- setNames(lapply(1:3, function(run) {
  design %*% truth + matrix(rnorm(9 * 3, sd = 0.3), 9, 3)
}), paste0("run", 1:3))
domain <- abstract_domain(3, id = "sampling-covariance-example")
fit <- lm_relation_fit(
  responses, design, effects, effect_names = conditions,
  sampling_unit = "trial", domain = domain
)
plan <- plan_geometry(
  fit$relation, compile_frame(whole_brain(), domain),
  cross_partitions(fit$relation, independence = "independent"),
  metric = noise_precision(diag(3), domain, covariance = diag(3))
)
uncertainty <- rdm_sampling_covariance(plan, fit, target = "null", at = 1L)

# The diagonal is a vector of variances, so take a square root yourself if
# you want standard errors. No interval is implied.
round(sqrt(sampling_covariance(uncertainty)), 4)
#> face - house  face - tool house - tool 
#>       0.0146       0.0146       0.0146 

# Distances that share a condition covary; read one entry without
# building the full matrix.
uncertainty$labels
#> [1] "face - house" "face - tool"  "house - tool"
round(sampling_covariance(
  uncertainty, "selected_entries", query = cbind(1, 2)
), 6)
#> [1] 5.3e-05

# The variance of a fixed contrast of distances, again without the matrix.
round(sampling_covariance(
  uncertainty, "quadratic_form", query = c(1, -1, 0)
), 6)
#> [1] 0.00032

# Materialization is explicit, never a silent fallback.
round(sampling_covariance(uncertainty, "materialize"), 6)
#>              face - house face - tool house - tool
#> face - house     0.000213    0.000053     0.000053
#> face - tool      0.000053    0.000213     0.000053
#> house - tool     0.000053    0.000053     0.000213

# A bank of contrast-energy queries moves the same covariance into the
# coordinates of the energies `contrast_energy()` reports.
bank <- rbind(
  animacy = c(face = 1, house = -1, tool = 0),
  tools = c(face = -0.5, house = -0.5, tool = 1)
)
energies <- sampling_covariance(uncertainty, queries = bank)
energies$basis
#> [1] "query_bank"
round(sampling_covariance(energies, "materialize"), 6)
#>          animacy   tools
#> animacy 0.000213 0.00000
#> tools   0.000000 0.00012

# It is the transport of the distance covariance through the operator that
# lowers each query onto the distances, and reports itself as such.
round(energies$source$query_lowering, 3)
#>         face - house face - tool house - tool
#> animacy         1.00         0.0          0.0
#> tools          -0.25         0.5          0.5

# Cross-measurement covariance is refused, not approximated.
catch_refusal(
  sampling_covariance(uncertainty, queries = bank,
    scope = "cross_measurement")
)$capability
#> [1] "cross_node_sampling_covariance"
```
