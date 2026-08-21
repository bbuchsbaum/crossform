# Plan a population form over transported conservative geometry

`plan_population()` names a group-level estimand: each participant's
compiled conservative geometry, one declared
[`location_transport()`](https://bbuchsbaum.github.io/crossform/reference/location_transport.md)
each onto one shared set of group nodes, and one group-level model over
the subject axis. It reads no data and fits nothing — it is a sealed
declaration of what a group number will mean, and a set of refusals for
the combinations that would make that meaning false.

## Usage

``` r
plan_population(
  subjects,
  transport,
  model = ~1,
  data = NULL,
  normalization = c("none", "unit_budget", "precision_weighted"),
  compute = compute_policy(),
  allow_nonconservative = FALSE
)
```

## Arguments

- subjects:

  A named list of compiled `effect_geometry_plan` values from
  [`plan_geometry()`](https://bbuchsbaum.github.io/crossform/reference/plan_geometry.md),
  one per participant, over conservative frames on a common experimental
  space. Names are participant identifiers and bind `transport` and
  `data`. The list is stored in a canonical name order, so the plan's
  identity does not depend on the order it was given in.

- transport:

  A named list of
  [`location_transport()`](https://bbuchsbaum.github.io/crossform/reference/location_transport.md)
  values, one per participant with matching names, all landing on one
  shared group node set under one `semantics`.

- model:

  A one-sided group-level formula evaluated against one row per
  participant. The default `~ 1` is the group mean.

- data:

  A data frame with one row per participant, identified by a `subject`
  column or by row names holding the participant identifiers. `NULL`
  (the default) supplies the empty table an intercept-only model needs.

- normalization:

  How incommensurable per-participant budgets are made commensurable:
  `"none"` (the default), `"unit_budget"`, or `"precision_weighted"`
  (gated, see above).

- compute:

  A
  [`compute_policy()`](https://bbuchsbaum.github.io/crossform/reference/compute_policy.md).
  It enters the plan's physical signature, never its scientific
  identity.

- allow_nonconservative:

  Plan on non-conservative subject frames anyway. Recorded on the plan
  and part of its scientific identity, because the result is a different
  estimand rather than a relaxed setting on the same one.

## Value

An `effect_population_plan`: a sealed list carrying `$subjects`,
`$transport`, the shared `$group_index`, the transport `$semantics`, the
`$model` record (formula, canonical text, term labels, model matrix, its
pivoted `$qr`, `$rank` and `$pivot`, and per-covariate content digests),
`$data`, `$normalization`, the `$fit` marker, `$allow_nonconservative`,
a per-participant `$subject_index` audit table, `$compute`, the
`$scientific_plan_id` naming the estimand, and the `$signature` covering
how it will be executed.

## What the plan fixes

Three choices decide what the group number *is*, and all three enter the
plan's scientific identity rather than being applied silently.

- **The transport.** The group form is not "the population geometry"; it
  is the population geometry *as resolved by* `P`. Two transports give
  two different, equally valid estimands.

- **The normalization.** Participants have incommensurable native
  budgets, so combining them requires a declared choice: `"none"` is the
  mean subject ledger in native evidence units (participants with more
  evidence weigh more), and `"unit_budget"` is the mean attribution
  *share*, in which every participant contributes budget one regardless
  of how much evidence they have. On the contract's own fixture the two
  disagree by up to 69 % and pick different group nodes as the maximum.

- **The fit and its order.** The group fit is OLS. Subject-constant
  weights are the commuting class: query, transport and fit are then
  linear maps on distinct axes and every evaluation order gives the same
  answer. With heterogeneous native frames there is no common node axis
  at all, so transport must precede the fit; the order is recorded
  rather than chosen.

## Refusals

Each is an `effect_capability_refusal` carrying the missing capability,
all unmet reasons, and a remedy (see
[`catch_refusal()`](https://bbuchsbaum.github.io/crossform/reference/catch_refusal.md)).

- `common_experimental_space` — participants' geometry indexed by
  different effects, orders, bases, units or scales.

- `symmetric_packed_population` — a rectangular cross-axis subject plan.

- `fixed_metric_population` — a subject plan on a learned neural metric.

- `conservative_subject_geometry` — a subject frame that is not
  conservative, either by declaration or after a diagonal metric fold,
  unless `allow_nonconservative = TRUE`.

- `aligned_subject_transports` — transport and subject names that do not
  match, or a transport whose native rows do not count the participant's
  frame measurements.

- `shared_group_nodes` — transports landing on different group node
  sets, or declaring different `semantics`.

- `aligned_subject_rows` — a `data` table whose rows do not name their
  participants.

- `identified_group_model` — a two-sided formula, an absent or
  incomplete covariate, a non-finite design, or a rank-deficient design
  (the refusal names the aliased columns).

- `precision_weighted_normalization` — see below.

## The `precision_weighted` gate

`"precision_weighted"` weighs each participant by \\\pi_i =
1/\mathrm{Var}(\hat T_i)\\, the precision of a *conserved budget*. That
variance needs the full cross-node sampling covariance rather than a sum
of per-node margins, and crossform has no route that produces one:
[`sampling_covariance()`](https://bbuchsbaum.github.io/crossform/reference/sampling_covariance.md)
reads a query bank at one measurement and refuses
`scope = "cross_measurement"` for want of a spatial model. The mode is
therefore refused rather than approximated — overlapping conservative
rows are strongly correlated, so a precision assembled from per-node
margins would reweight participants by the wrong quantity, and the
contract measures the three normalizations disagreeing by up to 94 %.
The mode stays in the closed set, and in this argument's default,
because the set *is* the plan identity; what is gated is admitting it.

## References

`design/population-form-contract.md` (`population-form-v1`), sections
1–4 and 9.

## See also

[`location_transport()`](https://bbuchsbaum.github.io/crossform/reference/location_transport.md)
and
[`anatomical_transport()`](https://bbuchsbaum.github.io/crossform/reference/anatomical_transport.md)
for the operator,
[`transport_values()`](https://bbuchsbaum.github.io/crossform/reference/transport_values.md)
for the contraction it performs, and
[`plan_geometry()`](https://bbuchsbaum.github.io/crossform/reference/plan_geometry.md)
for the per-participant input.

Other population transports:
[`anatomical_transport()`](https://bbuchsbaum.github.io/crossform/reference/anatomical_transport.md),
[`estimate_population()`](https://bbuchsbaum.github.io/crossform/reference/estimate_population.md),
[`external_transport()`](https://bbuchsbaum.github.io/crossform/reference/external_transport.md),
[`heterogeneity()`](https://bbuchsbaum.github.io/crossform/reference/heterogeneity.md),
[`location_transport()`](https://bbuchsbaum.github.io/crossform/reference/location_transport.md),
[`materialize_population()`](https://bbuchsbaum.github.io/crossform/reference/materialize_population.md),
[`population_prevalence()`](https://bbuchsbaum.github.io/crossform/reference/population_prevalence.md),
[`population_uncertainty()`](https://bbuchsbaum.github.io/crossform/reference/population_uncertainty.md),
[`population_views`](https://bbuchsbaum.github.io/crossform/reference/population_views.md),
[`transport_values()`](https://bbuchsbaum.github.io/crossform/reference/transport_values.md)

## Examples

``` r
# Two participants, different native frames, one shared pair of group nodes.
effects <- effect_space(c("face", "house"), basis_id = "conditions:v1")
subject_plan <- function(id, n) {
  domain <- abstract_domain(n, coordinates = cbind(x = seq_len(n) - 1),
    feature_ids = paste0("f", seq_len(n)), id = id)
  values <- function() matrix(seq_len(2 * n) / n, 2, n,
    dimnames = list(c("face", "house"), NULL))
  rel <- relation(list(run1 = values(), run2 = values() / 2),
    effects = effects, domain = domain)
  plan_geometry(rel, compile_frame(voxelwise(), domain),
    cross_partitions(rel))
}
carrier <- function(n) anatomical_transport(
  native_coords = cbind(seq_len(n) - 1), group_coords = cbind(c(0, 3)),
  semantics = "budget"
)

plan <- plan_population(
  subjects = list(s01 = subject_plan("s01", 5L), s02 = subject_plan("s02", 6L)),
  transport = list(s01 = carrier(5L), s02 = carrier(6L))
)
plan
#> <effect_population_plan>
#>   subjects:      2 (s01, s02)
#>   group nodes:   2 + sink
#>   sink:          empty in every subject (full native coverage)
#>   transport:     budget, anatomical
#>   model:         ~1 -> 1 column, rank 1
#>   normalization: none
#>   fit:           OLS (subject-constant weights), transport then fit
#>   estimand:      population-sha256:7cd6972df1d6...
#>   signature:     sha256:43de83dc28a5...
plan$subject_index[, c("subject", "declared_normalization", "sink_territory")]
#>   subject declared_normalization sink_territory
#> 1     s01           conservative              0
#> 2     s02           conservative              0

# A group covariate: rows name the participant they belong to.
covariates <- data.frame(age = c(24, 31), row.names = c("s01", "s02"))
with_age <- plan_population(
  subjects = list(s01 = subject_plan("s01", 5L), s02 = subject_plan("s02", 6L)),
  transport = list(s01 = carrier(5L), s02 = carrier(6L)),
  model = ~ age, data = covariates
)
with_age$model$columns
#> [1] "(Intercept)" "age"        

# The transport is part of the estimand, so changing it moves the identity.
identical(plan$scientific_plan_id, with_age$scientific_plan_id)
#> [1] FALSE
```
