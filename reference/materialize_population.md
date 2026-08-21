# Materialize a planned population form at every group node

`materialize_population()` is the complete-geometry counterpart of
[`estimate_population()`](https://bbuchsbaum.github.io/crossform/reference/estimate_population.md).
Where that verb reads each participant's geometry through a declared
bank of `K` contrasts, this one carries the **whole** \\q \times q\\
form to the group nodes and fits the group model at every packed
coordinate, so the result is a coefficient *form* per group node and
model term rather than a coefficient per query.

## Usage

``` r
materialize_population(
  plan,
  component = c("total", "coherent", "configuration"),
  coordinate_tile = NULL
)
```

## Arguments

- plan:

  An `effect_population_plan` from
  [`plan_population()`](https://bbuchsbaum.github.io/crossform/reference/plan_population.md).

- component:

  Which component of each participant's conservative geometry to carry:
  `"total"` (the default), `"coherent"`, or `"configuration"`. The
  transported object takes the `population-form-v1` section 8.1 name for
  it, as in
  [`estimate_population()`](https://bbuchsbaum.github.io/crossform/reference/estimate_population.md).

- coordinate_tile:

  Optional positive number of packed coordinates held in flight at once.
  `NULL` takes the tile from the plan's
  [`compute_policy()`](https://bbuchsbaum.github.io/crossform/reference/compute_policy.md).
  The tile bounds memory and is not part of the estimand: it does not
  enter `$scientific_plan_id`, and two tile sizes agree to
  [`numerical_contract()`](https://bbuchsbaum.github.io/crossform/reference/numerical_contract.md)`$atol`.
  Not bit-identically — the tile is a blocking choice, and the
  contract's `bitwise_across_blocking` is `FALSE` because it changes the
  column count of one fused contraction.

## Value

An `effect_population_result` with `$basis` `"complete_form"`, carrying
`$coefficient_forms` (a `node`-by-`term`-by-`coordinate` array whose
last axis is the packed `svec` form, so every `[node, term, ]` slice is
one packed \\q \times q\\ coefficient form), `$effects`, the
`$coordinates` table describing the packed codec, the `$index` of group
nodes plus the sink, `$ledger`, and a `$receipt` recording every
participant's read, its transport signature, the budget certificate and
the streaming bound.
[`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html) returns
the coefficient table in long form.

## Details

It stands to
[`estimate_population()`](https://bbuchsbaum.github.io/crossform/reference/estimate_population.md)
exactly as
[`materialize_geometry()`](https://bbuchsbaum.github.io/crossform/reference/materialize_geometry.md)
stands to
[`evaluate_geometry()`](https://bbuchsbaum.github.io/crossform/reference/evaluate_geometry.md):
two verbs, because complete materialization is a different memory regime
and a different estimand, and because a caller who forgets a query bank
should get
[`estimate_population()`](https://bbuchsbaum.github.io/crossform/reference/estimate_population.md)'s
refusal rather than a silently enormous complete form.

## How it streams, and what that bounds

Reading a complete packed form is reading it through the identity, and a
coordinate tile is a column block of the identity. For each tile the
route opens each participant in turn, contracts that participant's
native nodes against the tile (never allocating the full \\n_i \times
p\\ packed block), transports the tile's rows with that participant's
operator, stacks the transported tile across participants, and solves
the plan's prefactored QR once for the tile. Coefficient tiles come out
and are written into the packed coefficient forms.

Nothing of size \\N \times (m+1) \times p\\ is ever allocated. The
arrays in flight are the group stack, \\N \times (m+1) \times
\mathrm{tile}\\, and one native block, \\\max_i n_i \times
\mathrm{tile}\\; the durable output is \\(m+1) \times k \times p\\ for
`k` model terms, which is smaller than the refused array whenever the
design is not saturated. The receipt records all four numbers under
`$streaming`. The cost is \\\lceil p/\mathrm{tile}\rceil\\ passes over
each participant's source instead of one.

## What the result does not carry

There are no `$values`, `$fitted` or `$residuals` arrays. Those are
indexed by participant, group node and packed coordinate — they *are*
the \\N \times (m+1) \times p\\ object this route exists in order not to
build. Read the participant-level transported values through
[`estimate_population()`](https://bbuchsbaum.github.io/crossform/reference/estimate_population.md)
with the contrasts you care about, where the subject axis costs `K` and
not `p`.

## Normalization

Only `normalization = "none"` is admitted; see the refusal below. Budget
and density semantics are both admitted, and the budget certificate of
`population-form-v1` section 2 is asserted per participant and per
packed coordinate, exactly as
[`estimate_population()`](https://bbuchsbaum.github.io/crossform/reference/estimate_population.md)
asserts it per query.

## Refusals

Each is an `effect_capability_refusal` (see
[`catch_refusal()`](https://bbuchsbaum.github.io/crossform/reference/catch_refusal.md)).

- `complete_form_normalization` — a plan declaring `"unit_budget"` or
  `"precision_weighted"`. Their divisor is the native total of a query,
  and a form is not read through one; a per-coordinate divisor would
  break the query-fit commutation the form's queryability rests on.

- `reserved_group_index_columns` — a group index already carrying `sink`
  or `units`, as in
  [`estimate_population()`](https://bbuchsbaum.github.io/crossform/reference/estimate_population.md).

## References

`design/population-form-contract.md` (`population-form-v1`), sections 2,
3, 4, 5 and 8.

## See also

[`estimate_population()`](https://bbuchsbaum.github.io/crossform/reference/estimate_population.md)
for the query-first path over the same plan,
[`plan_population()`](https://bbuchsbaum.github.io/crossform/reference/plan_population.md)
for the estimand both execute, and
[`materialize_geometry()`](https://bbuchsbaum.github.io/crossform/reference/materialize_geometry.md)
for the single-participant analogue.

Other population transports:
[`anatomical_transport()`](https://bbuchsbaum.github.io/crossform/reference/anatomical_transport.md),
[`estimate_population()`](https://bbuchsbaum.github.io/crossform/reference/estimate_population.md),
[`external_transport()`](https://bbuchsbaum.github.io/crossform/reference/external_transport.md),
[`heterogeneity()`](https://bbuchsbaum.github.io/crossform/reference/heterogeneity.md),
[`location_transport()`](https://bbuchsbaum.github.io/crossform/reference/location_transport.md),
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

form <- materialize_population(plan)
form
#> <effect_population_result>
#>   ledger:        transported_total (component "total")
#>   subjects:      4 (s01, s02, s03, s04), residual df 3
#>   group nodes:   2 + sink
#>   readout:       complete 2x2 form, 3 packed coordinates (face:face, hous...
#>   streaming:     coordinate tile 3 of 3, 1 pass/participant; peak group s...
#>   frame:         undeclared, conservative
#>   transport:     budget, anatomical, cross-fit not declared
#>   normalization: none (mean subject ledger, native evidence units)
#>   fit:           OLS (subject-constant weights), transport then fit
#>   budget:        preserved, worst relative deviation 1.65e-16 against 1e-12
#>   estimand:      population-sha256:df8915ebdc51...
#>   No participant-level arrays: indexed by participant, group node and
#>     packed coordinate they are the dense stack this route streams in order
#>     not to build. Read them through estimate_population() with the
#>     contrasts you want.
#>   next:          as.data.frame(x), x$coefficient_forms[node, term, ]

# The group-mean coefficient form at the first group node: three packed
# coordinates, the off-diagonal one carrying the codec's root two.
form$coordinates
#>    coordinate row column    scale
#> 1   face:face   1      1 1.000000
#> 2  house:face   2      1 1.414214
#> 3 house:house   2      2 1.000000
form$coefficient_forms["group1", "(Intercept)", ]
#>   face:face  house:face house:house 
#>   0.1541884   0.3052774   0.3083767 

# Querying the assembled form with a contrast is the same number the
# query-first executor returns for that contrast.
contrast <- c(1, -1)
packed <- crossform:::.svec_symmetric(tcrossprod(contrast))
assembled <- form$coefficient_forms[, "(Intercept)", ] %*% packed
queried <- estimate_population(plan, rbind(`face-house` = contrast))
max(abs(assembled - queried$coefficients[, "face-house", "(Intercept)"]))
#> [1] 1.290634e-15

# The tile bounds memory and nothing else: one coordinate at a time gives
# the same forms, and the same estimand identity.
tiled <- materialize_population(plan, coordinate_tile = 1L)
max(abs(tiled$coefficient_forms - form$coefficient_forms))
#> [1] 0
identical(tiled$scientific_plan_id, form$scientific_plan_id)
#> [1] TRUE
tiled$receipt$streaming
#> $coordinate_tile
#> [1] 1
#> 
#> $packed_width
#> [1] 3
#> 
#> $passes_per_subject
#> [1] 3
#> 
#> $group_stack_doubles
#> [1] 12
#> 
#> $native_block_doubles
#> [1] 8
#> 
#> $output_doubles
#> [1] 9
#> 
#> $refused_dense_doubles
#> [1] 36
#> 
```
