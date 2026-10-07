# Reader verbs on an estimated population form

[`contrast_energy()`](https://bbuchsbaum.github.io/crossform/reference/contrast_energy.md),
[`rdm()`](https://bbuchsbaum.github.io/crossform/reference/rdm.md),
[`rsa()`](https://bbuchsbaum.github.io/crossform/reference/rsa.md) and
[`contribution()`](https://bbuchsbaum.github.io/crossform/reference/contribution.md)
are S3 generics. Applied to an
[`estimate_population()`](https://bbuchsbaum.github.io/crossform/reference/estimate_population.md)
result they read the *group* coefficients that run already produced: a
population view is a fixed linear combination of estimated query
columns, not a second execution, and no participant's geometry is opened
again.

## Usage

``` r
# S3 method for class 'effect_population_result'
contrast_energy(x, weights, remove_univariate = FALSE, term = NULL, ...)

# S3 method for class 'effect_population_result'
rdm(x, component = NULL, pairs = NULL, normalize = NULL, term = NULL, ...)

# S3 method for class 'effect_population_result'
rsa(
  x,
  models,
  nuisance = NULL,
  intercept = TRUE,
  component = NULL,
  term = NULL,
  ...
)

# S3 method for class 'effect_population_result'
contribution(x, by, using = NULL, term = NULL, ...)

# S3 method for class 'effect_population_view'
contribution(x, by, using = NULL, ...)
```

## Arguments

- x:

  An `effect_population_result` from
  [`estimate_population()`](https://bbuchsbaum.github.io/crossform/reference/estimate_population.md)
  or
  [`materialize_population()`](https://bbuchsbaum.github.io/crossform/reference/materialize_population.md),
  or – for
  [`contribution()`](https://bbuchsbaum.github.io/crossform/reference/contribution.md)
  – an `effect_population_view` from one of the other three verbs.

- weights:

  One finite contrast weight per experimental effect, named or in the
  effect order of `$queries`.

- remove_univariate:

  Must be omitted or `FALSE`; refused exactly as
  [`contrast_energy()`](https://bbuchsbaum.github.io/crossform/reference/contrast_energy.md)
  refuses it.

- term:

  Which column of the group model to read, by name or position. Required
  when the group model fits more than one; the default `~ 1` fits only
  `(Intercept)`, the group mean, and needs no argument.

- ...:

  Unused; present so the methods match their generics.

- component:

  Must be omitted. The geometry component is fixed when the population
  is estimated – it is `component =` on
  [`estimate_population()`](https://bbuchsbaum.github.io/crossform/reference/estimate_population.md)
  – because a transported ledger of one component holds no information
  about another.

- pairs:

  Optional two-column matrix of effect pairs, by name or index. The
  default reports every unordered pair.

- normalize:

  Must be omitted; refused exactly as
  [`rdm()`](https://bbuchsbaum.github.io/crossform/reference/rdm.md)
  refuses it.

- models, nuisance:

  Named model and nuisance dissimilarity matrices over the experimental
  effects, as
  [`rsa()`](https://bbuchsbaum.github.io/crossform/reference/rsa.md)
  takes them.

- intercept:

  Whether the RSA design carries an intercept column.

- by:

  The grouping of group nodes, as one label per group node (the sink
  excluded, since it is always its own group) or one column name of
  `using` or of the view's `$index`.

- using:

  Optional per-group-node metadata table naming the territories a view
  cannot carry itself. It is **joined on its `node` column**, which must
  hold the group node identifiers of `$index`, so row order does not
  matter and a table without that column is refused rather than bound by
  position. Meaningful only when `by` names one of its columns.

## Value

An `effect_population_view`: `$values`, one row per group node plus the
sink and one column per view column; `$columns` naming those columns;
`$index`, the group node table with the sink marked and its units; the
`$ledger` name, `$term`, `$semantics`, `$normalization`; a `$receipt`
carrying the transport, the native frame family, the budget certificate,
the normalization and the basis coefficients that reached this view; and
a `$scientific_plan_id` derived from the population plan and the view's
own parameters.
[`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html) returns
the long table.

## Why a combination is exact

The query, the transport and the group fit act on the experimental,
spatial and participant axes respectively, so under the OLS default they
commute (`population-form-v1` section 3). The fitted coefficient of a
combined query is therefore the same combination of the fitted
coefficients of the bank's queries, exactly rather than approximately.

## What the estimated basis reaches

A result from
[`estimate_population()`](https://bbuchsbaum.github.io/crossform/reference/estimate_population.md)
(`$basis` `"query_bank"`) holds the group coefficient of `K` packed
operators \\\mathrm{svec}(c_k c_k^\top)\\. A view is derivable exactly
when its own packed operator lies in their span, which is a wider set
than the bank's contrast vectors: `2c` is derivable when `c` is in the
bank, and the bank of all \\q(q-1)/2\\ pairwise difference contrasts
spans **exactly** the quadratic forms of zero-sum contrasts. So a full
pairwise bank reaches every
[`rdm()`](https://bbuchsbaum.github.io/crossform/reference/rdm.md) and
every [`rsa()`](https://bbuchsbaum.github.io/crossform/reference/rsa.md)
regression over those effects, and reaches no uncentred contrast at all.
(Other banks reach other things: three uncentred contrasts can reach a
fourth. The span is the criterion; the pairwise case is the one worth
stating because it is the one
[`rdm()`](https://bbuchsbaum.github.io/crossform/reference/rdm.md) and
[`rsa()`](https://bbuchsbaum.github.io/crossform/reference/rsa.md) rest
on.) Outside the span the verb refuses and names the re-estimation
remedy rather than reporting the nearest thing it can reach under the
query's name.

A result from
[`materialize_population()`](https://bbuchsbaum.github.io/crossform/reference/materialize_population.md)
(`$basis` `"complete_form"`) holds the assembled packed coefficient
form, whose basis is the whole coordinate system. Every symmetric query
is in its span by construction, so these verbs query the form directly
and no span refusal can fire.

## `unit_budget` admits one estimated query at a time

Under `normalization = "unit_budget"` each query column is divided by
*that query's* own native total (`population-form-v1` section 4.1), so
the divisor varies along the query axis and the response is not linear
in the query operator. A view that reads **one** estimated query is
exact and is admitted, including a multiple of one: `s H_k` has ledger
`s c_k` over total `s T_{ik}`, so its share is column `k` for every
nonzero `s`, and the weight is set to one rather than to the \\s\\ a
span solve returns. A view that **mixes** estimated queries is refused:
the mixture carries no denominator, and a share of nothing is not a
share.

## The sink

Every view carries the sink as a row of `$index`, labelled `<sink>`,
marked by `$index$sink` and reported in budget units. It is never a
value at a location; it is there because
[`contribution()`](https://bbuchsbaum.github.io/crossform/reference/contribution.md)'s
ledger adds up only when unmapped territory is a number in the table.

## Labelling

A transported component is a `native_coherent_ledger` or a
`native_configuration_ledger`, never a bare `coherent` or
`configuration` (`population-form-v1` section 8.1). `$ledger` carries
the name, `$receipt` carries it beside the native frame family the
ledger belongs to and the transport that carried it, and no output
column or printed label of a population view uses the bare name.
`native_coherent_ledger + native_configuration_ledger = transported_total`
holds exactly; what fails is reading the coherent ledger as a group-node
common mode.

## Refusals

Each is an `effect_capability_refusal` in namespace `"population_views"`
(see
[`catch_refusal()`](https://bbuchsbaum.github.io/crossform/reference/catch_refusal.md)).

- `population_query_in_bank` — the view's packed operator is outside the
  span of the estimated bank. Remedy: re-estimate with the contrast in
  the bank.

- `population_view_query_linearity` — a view mixing estimated queries on
  a `unit_budget` result, or any normalization no view knows how to
  combine under.

- `population_component_fixed_at_estimation` — `component`, which is an
  argument of
  [`estimate_population()`](https://bbuchsbaum.github.io/crossform/reference/estimate_population.md)
  and not of a view.

- `budget_semantics` —
  [`contribution()`](https://bbuchsbaum.github.io/crossform/reference/contribution.md)
  on a `"density"` population. Density satisfies no conservation law, so
  a territory sum of densities is a share of nothing.

- `conservative_frame` —
  [`contribution()`](https://bbuchsbaum.github.io/crossform/reference/contribution.md)
  on a population planned with `allow_nonconservative = TRUE` over a
  frame that is not column normalized.

- `sink_is_not_a_territory` — a `by` that labels a group node `<sink>`.
  The sink is appended as its own row automatically.

- `nondestructive_decomposition` and `guaranteed_psd` —
  `remove_univariate` and `normalize`, refused for the same reasons the
  per-participant views refuse them.

## References

`design/population-form-contract.md` (`population-form-v1`), sections 2,
3, 4 and 8.

## See also

[`estimate_population()`](https://bbuchsbaum.github.io/crossform/reference/estimate_population.md)
and
[`materialize_population()`](https://bbuchsbaum.github.io/crossform/reference/materialize_population.md)
for the runs these read, and
[`contrast_energy()`](https://bbuchsbaum.github.io/crossform/reference/contrast_energy.md),
[`rdm()`](https://bbuchsbaum.github.io/crossform/reference/rdm.md),
[`rsa()`](https://bbuchsbaum.github.io/crossform/reference/rsa.md) and
[`contribution()`](https://bbuchsbaum.github.io/crossform/reference/contribution.md)
for the per-participant verbs of the same names.

Other population transports:
[`anatomical_transport()`](https://bbuchsbaum.github.io/crossform/reference/anatomical_transport.md),
[`estimate_population()`](https://bbuchsbaum.github.io/crossform/reference/estimate_population.md),
[`external_transport()`](https://bbuchsbaum.github.io/crossform/reference/external_transport.md),
[`heterogeneity()`](https://bbuchsbaum.github.io/crossform/reference/heterogeneity.md),
[`location_transport()`](https://bbuchsbaum.github.io/crossform/reference/location_transport.md),
[`materialize_population()`](https://bbuchsbaum.github.io/crossform/reference/materialize_population.md),
[`plan_population()`](https://bbuchsbaum.github.io/crossform/reference/plan_population.md),
[`population_prevalence()`](https://bbuchsbaum.github.io/crossform/reference/population_prevalence.md),
[`population_uncertainty()`](https://bbuchsbaum.github.io/crossform/reference/population_uncertainty.md),
[`transport_values()`](https://bbuchsbaum.github.io/crossform/reference/transport_values.md)

## Examples

``` r
# Three participants on different native frames, two group nodes, and a
# bank spanning all three pairwise contrasts.
effects <- effect_space(c("face", "house", "tool"), basis_id = "popview:v1")
subject <- function(id, n) {
  domain <- abstract_domain(n, coordinates = cbind(x = seq_len(n) - 1),
    feature_ids = paste0("f", seq_len(n)), id = id)
  values <- function(divisor) matrix(
    seq_len(3 * n) / (n * divisor), 3, n,
    dimnames = list(c("face", "house", "tool"), NULL)
  )
  rel <- relation(list(run1 = values(1), run2 = values(1.4)),
    effects = effects, domain = domain)
  plan_geometry(rel, compile_frame(voxelwise(), domain),
    cross_partitions(rel))
}
carrier <- function(n) anatomical_transport(
  native_coords = cbind(seq_len(n) - 1), group_coords = cbind(c(0, 4)),
  semantics = "budget"
)
sizes <- c(s01 = 5L, s02 = 6L, s03 = 7L)
subjects <- stats::setNames(
  lapply(names(sizes), function(id) subject(id, sizes[[id]])), names(sizes)
)
fit <- estimate_population(
  plan_population(subjects, lapply(sizes, carrier)),
  rbind(`face-house` = c(1, -1, 0), `face-tool` = c(1, 0, -1),
    `house-tool` = c(0, 1, -1))
)

# One contrast, selected out of the estimated bank.
contrast_energy(fit, c(face = 1, house = -1, tool = 0))
#> <effect_population_view>
#>   view:          contrast energy
#>   ledger:        transported_total
#>   term:          (Intercept)
#>   group nodes:   2 + sink
#>   columns:       1 (energy)
#>   frame:         undeclared, conservative
#>   transport:     budget, anatomical, cross-fit not declared
#>   normalization: none (mean subject ledger, native evidence units)
#>   basis:         selection over 3 estimated queries (face-house, face-too...
#>   budget:        preserved, worst relative deviation 1.94e-16 against 1e-12
#>   estimand:      population-sha256:4bb36d3c7b5f...
#>   next:          as.data.frame(x), contribution(x, by = ...)

# The three distances, and an RSA regression on them: both are inside the
# span of a full pairwise bank, so neither needs a second run.
rdm(fit)$values
#>        face - house face - tool house - tool
#> group1   0.06298996   0.2519598   0.06298996
#> group2   0.05832523   0.2333009   0.05832523
#> <sink>   0.00000000   0.0000000   0.00000000
animacy <- matrix(c(0, 0, 1, 0, 0, 1, 1, 1, 0), 3, 3,
  dimnames = list(c("face", "house", "tool"), c("face", "house", "tool")))
as.data.frame(rsa(fit, models = list(animacy = animacy)))
#>     node coord1  sink  units view            ledger        term view_column
#> 1 group1      0 FALSE budget  rsa transported_total (Intercept) (Intercept)
#> 2 group2      4 FALSE budget  rsa transported_total (Intercept) (Intercept)
#> 3 <sink>     NA  TRUE budget  rsa transported_total (Intercept) (Intercept)
#> 4 group1      0 FALSE budget  rsa transported_total (Intercept)     animacy
#> 5 group2      4 FALSE budget  rsa transported_total (Intercept)     animacy
#> 6 <sink>     NA  TRUE budget  rsa transported_total (Intercept)     animacy
#>        role   estimate
#> 1 intercept 0.06298996
#> 2 intercept 0.05832523
#> 3 intercept 0.00000000
#> 4     model 0.09448494
#> 5     model 0.08748785
#> 6     model 0.00000000

# The ledger adds up over the group nodes and the sink.
ledger <- contribution(fit, by = c("anterior", "posterior"))
ledger
#> <effect_population_view>
#>   view:          aggregated ledger
#>   ledger:        transported_total
#>   term:          (Intercept)
#>   group nodes:   2 territories + sink, from 2 group nodes
#>   columns:       3 (face-house, face-tool, house-tool)
#>   frame:         undeclared, conservative
#>   transport:     budget, anatomical, cross-fit not declared
#>   normalization: none (mean subject ledger, native evidence units)
#>   basis:         whole bank over 3 estimated queries (face-house, face-to...
#>   budget:        preserved, worst relative deviation 1.94e-16 against 1e-12
#>   aggregation:   by group, budget-exact; the sink is its own row
#>   estimand:      population-sha256:54230de1614b...
#>   next:          as.data.frame(x), x$values

# An uncentred contrast is outside the span of a zero-sum bank, and is
# refused with the re-estimation remedy rather than projected onto it.
refusal <- catch_refusal(contrast_energy(fit, c(1, 1, 0)))
refusal$capability
#> [1] "population_query_in_bank"
```
