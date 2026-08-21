# Estimate a planned population form through a bank of queries

`estimate_population()` is the authorized group-level execution verb: it
takes the estimand
[`plan_population()`](https://bbuchsbaum.github.io/crossform/reference/plan_population.md)
sealed, reads each participant's geometry through one bank of `K`
contrasts, carries the queried ledgers onto the shared group nodes with
that participant's declared transport, and fits the group model once at
every group node and query.

## Usage

``` r
estimate_population(
  plan,
  queries,
  component = c("total", "coherent", "configuration"),
  uncertainty = NULL,
  budget_floor = 0
)
```

## Arguments

- plan:

  An `effect_population_plan` from
  [`plan_population()`](https://bbuchsbaum.github.io/crossform/reference/plan_population.md).

- queries:

  The bank of contrasts to read: a numeric contrast vector, a `K`-by-`q`
  matrix with one contrast per row, or a list of contrast vectors, over
  the participants' shared experimental effects. Named rows name the
  queries; unnamed ones are numbered. This is the same bank
  [`sampling_covariance()`](https://bbuchsbaum.github.io/crossform/reference/sampling_covariance.md)
  takes.

- component:

  Which component of each participant's conservative geometry to carry:
  `"total"` (the default), `"coherent"`, or `"configuration"`. See the
  table above for the name the transported object takes.

- uncertainty:

  Optional named list of
  [`rdm_sampling_covariance()`](https://bbuchsbaum.github.io/crossform/reference/rdm_sampling_covariance.md)
  results, one per participant, lowered onto the query bank and carried
  **untransported**. The covariance of a transported value is refused
  except where a group column is fed by exactly one native row, and both
  the refusal and the admitted carve-out travel with the result. The
  between-subject layer needs none of this and is always available.

- budget_floor:

  Nonnegative relative floor for the `"unit_budget"` divisor: a
  participant contributes a share for query `k` only when \\\|T_i\| \>
  \mathrm{floor} \cdot \sum_x \|c\_{ix}\|\\. Ignored under `"none"`.

## Value

An `effect_population_result`: a sealed record carrying `$coefficients`
(a `node`-by-`query`-by-`term` array), `$values`, `$fitted` and
`$residuals` (`node`-by-`query`-by-`subject`), the `$index` of group
nodes plus the sink, `$queries`, `$ledger`, `$uncertainty` (see the
section above; read it through
[`population_uncertainty()`](https://bbuchsbaum.github.io/crossform/reference/population_uncertainty.md)),
and a `$receipt` recording every participant's read, its transport
signature, the budget certificate, and the normalization.
[`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html) returns
the coefficient table in long form.

## The order, and why it is not a choice

The evaluation order is **query, then transport, then fit**. Query and
transport act on different axes of the same array — the query mixes
experimental coordinates, the transport mixes spatial nodes — so they
commute exactly, and under a subject-constant weight operator the fit
commutes with both. What is *not* free is which may run first:
participants have different native frames, so there is no common node
axis to fit at before transporting, and fit-then-transport is not
definable at all. The order is recorded on the plan rather than selected
here.

Reading the query first is what makes the route cheap. Complete geometry
at a native node is \\q(q+1)/2\\ packed coordinates; a bank of `K`
queries is `K` numbers, and complete geometry is never allocated.

## What the result is a ledger of

`component` selects which part of each participant's conservative
geometry is carried, and the transported object takes its own name
(`population-form-v1` section 8.1) because it is not the group-node
quantity the bare name would suggest.

|  |  |  |
|----|----|----|
| `component` | `$ledger` | what it is |
| `"total"` | `transported_total` | a genuine group-node total for that participant |
| `"coherent"` | `native_coherent_ledger` | native-node coherent evidence carried to a group node |
| `"configuration"` | `native_configuration_ledger` | native-node configuration evidence carried to a group node |

`native_coherent_ledger + native_configuration_ledger = transported_total`
holds exactly; what fails is reading the coherent ledger as a group-node
common mode. Transport carries *nodes*, so there is no group frame and
no group-node weight vector for such a mode to be defined against.

The sink is a row of the result, labelled `<sink>`, not a discarded
remainder: it is fitted like any other column because covariate-linked
sink mass is how a differentially failing transport becomes visible. It
is never a value at a location, and it is always reported in budget
units.

## Budget preservation

Under `"budget"` semantics each participant's transported total, sink
included, equals its native total exactly. The certificate is asserted
per participant and per query at fit time, with a tolerance of `1e-12`
times the ledger's \\L^1\\ norm — not its total, which is a signed sum
that may sit near zero while its summands are large. A participant that
fails it has an ill-formed transport, and the run refuses rather than
renormalizing.

## Normalization

`"none"` fits the ledgers as they are, so participants with more
evidence weigh more. `"unit_budget"` divides each participant's ledger
by its own native total first, so each contributes budget one; the
divisor is a *signed* estimate, so a participant whose total is not
bounded away from zero by `budget_floor` is marked `NA` for that query
rather than emitting a divergent share. The default `budget_floor = 0`
admits every nonzero total, which is the weakest honest criterion and is
recorded as such in the receipt: the constant is an open maintainer
decision, and it needs a per-participant standard error that does not
exist yet. Both normalizations are computed against a total read from
the same data as the ledger, which the receipt declares
(`budget_estimate = "same_data_ratio"`).

## Uncertainty

`$uncertainty` always carries `$between`, the ingredients of the
between-subject layer: the group design's unscaled covariance
\\(X'X)^{-1}\\ and the terms it is indexed by. That is what
[`population_uncertainty()`](https://bbuchsbaum.github.io/crossform/reference/population_uncertainty.md)
needs and the result does not otherwise hold; the standard errors
themselves are computed when they are read rather than stored as three
more `node`-by-`query`-by-`term` arrays.

Supplying `uncertainty` adds two more blocks. `$native` is D8's own
per-native-node covariance per participant, carried **untransported**,
with `$transported` recording the refusal: a group node's value mixes
native nodes, so its covariance needs the covariance *between* them, and
no route produces one. `$within` is the part of that refusal E8 lifts —
a group column fed by exactly one native row, where the cross-node terms
carry weight zero — admitted per participant and per column and absent
elsewhere. No independence assumption is made anywhere.

The two layers are reported separately and are never pooled; there is no
field holding their sum. See
[`population_uncertainty()`](https://bbuchsbaum.github.io/crossform/reference/population_uncertainty.md).

## Refusals

Each is an `effect_capability_refusal` (see
[`catch_refusal()`](https://bbuchsbaum.github.io/crossform/reference/catch_refusal.md)).

- `budget_normalization_semantics` — `"unit_budget"` on a `"density"`
  plan. Density conserves nothing, so the share would be a share of
  nothing.

- `aligned_subject_uncertainty` — `uncertainty` names that do not match
  the participants.

- `distance_basis_query_bank` — an uncentred contrast in the bank when
  `uncertainty` is supplied. Point estimates admit uncentred contrasts;
  their sampling law does not.

## References

`design/population-form-contract.md` (`population-form-v1`), sections 2,
3, 4 and 8.

## See also

[`plan_population()`](https://bbuchsbaum.github.io/crossform/reference/plan_population.md)
for the estimand this executes,
[`transport_values()`](https://bbuchsbaum.github.io/crossform/reference/transport_values.md)
for the carrying step,
[`population_uncertainty()`](https://bbuchsbaum.github.io/crossform/reference/population_uncertainty.md)
for the error bars on what it produced, and
[`contrast_energy()`](https://bbuchsbaum.github.io/crossform/reference/contrast_energy.md)
for the single-participant reading of the same query.

Other population transports:
[`anatomical_transport()`](https://bbuchsbaum.github.io/crossform/reference/anatomical_transport.md),
[`external_transport()`](https://bbuchsbaum.github.io/crossform/reference/external_transport.md),
[`heterogeneity()`](https://bbuchsbaum.github.io/crossform/reference/heterogeneity.md),
[`location_transport()`](https://bbuchsbaum.github.io/crossform/reference/location_transport.md),
[`materialize_population()`](https://bbuchsbaum.github.io/crossform/reference/materialize_population.md),
[`plan_population()`](https://bbuchsbaum.github.io/crossform/reference/plan_population.md),
[`population_prevalence()`](https://bbuchsbaum.github.io/crossform/reference/population_prevalence.md),
[`population_uncertainty()`](https://bbuchsbaum.github.io/crossform/reference/population_uncertainty.md),
[`population_views`](https://bbuchsbaum.github.io/crossform/reference/population_views.md),
[`transport_values()`](https://bbuchsbaum.github.io/crossform/reference/transport_values.md)

## Examples

``` r
# Four participants on different native frames, two shared group nodes.
effects <- effect_space(c("face", "house"), basis_id = "population:v1")
subject <- function(id, n, gain) {
  domain <- abstract_domain(n, coordinates = cbind(x = seq_len(n) - 1),
    feature_ids = paste0("f", seq_len(n)), id = id)
  values <- function(divisor) matrix(
    gain * seq_len(2 * n) / (n * divisor), 2, n,
    dimnames = list(c("face", "house"), NULL)
  )
  rel <- relation(list(run1 = values(1), run2 = values(2)),
    effects = effects, domain = domain)
  plan_geometry(rel, compile_frame(voxelwise(), domain),
    cross_partitions(rel))
}
carrier <- function(n) anatomical_transport(
  native_coords = cbind(seq_len(n) - 1), group_coords = cbind(c(0, 3)),
  semantics = "budget"
)
sizes <- c(s01 = 5L, s02 = 6L, s03 = 7L, s04 = 8L)
gains <- c(s01 = 1, s02 = 1.4, s03 = 0.7, s04 = 1.1)
subjects <- stats::setNames(lapply(names(sizes), function(id)
  subject(id, sizes[[id]], gains[[id]])), names(sizes))
plan <- plan_population(subjects, lapply(sizes, carrier))

# One bank of two contrasts, read at every group node in one pass.
fit <- estimate_population(plan, rbind(
  `face-house` = c(1, -1), `face+house` = c(1, 1)
))
fit
#> <effect_population_result>
#>   ledger:        transported_total (component "total")
#>   subjects:      4 (s01, s02, s03, s04), residual df 3
#>   group nodes:   2 + sink
#>   queries:       2 (face-house, face+house)
#>   frame:         undeclared, conservative
#>   transport:     budget, anatomical, cross-fit not declared
#>   normalization: none (mean subject ledger, native evidence units)
#>   fit:           OLS (subject-constant weights), transport then fit
#>   budget:        preserved, worst relative deviation 2.02e-16 against 1e-12
#>   uncertainty:   between-subject SE, df 3 (uncalibrated)
#>   estimand:      population-sha256:fc0f369ba976...
#>   next:          as.data.frame(x), x$coefficients[, , term], population_uncertainty(x)

# The group mean at each node and query, with the sink as its own row.
fit$coefficients[, , "(Intercept)"]
#>         query
#> node     face-house face+house
#>   group1 0.03083767  0.8942925
#>   group2 0.06265191 21.3178776
#>   <sink> 0.00000000  0.0000000

# Budget preservation: each participant's transported total, sink
# included, is its native total.
fit$receipt$budget$max_relative_deviation < 1e-12
#> [1] TRUE

# The transported object is named for what it is.
fit$ledger
#> [1] "transported_total"
```
