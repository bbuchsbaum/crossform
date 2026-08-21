# Declare a location transport onto a shared group node set

A location transport carries each participant's native conservative
frame nodes onto a common set of group nodes. It is a **declaration**,
not an estimate: `crossform` accepts a typed, sparse, provenance-bearing
transport and evaluates it, and refuses to learn one (registration and
functional-transport learning stay outside the package).

## Usage

``` r
location_transport(
  matrix,
  native_index,
  group_index,
  semantics,
  provenance,
  row_mass = NULL,
  tolerance = 1e-09
)
```

## Arguments

- matrix:

  The `n_native` by `m` nonnegative group block, as a base matrix or a
  `Matrix`. Row sums may be at most one; the shortfall becomes sink
  mass.

- native_index:

  One row per native node: a character vector of node identifiers, or a
  data frame with a `node` column plus any per-row frame metadata
  (`family`, `scale`, `center`, `alpha` are type-checked when present
  and carried through unchanged).

- group_index:

  One row per group node: a character vector of node identifiers, or a
  data frame with a `node` column plus any group-node metadata, such as
  the `coord1`, `coord2`, ... centres the displacement diagnostics need.
  It describes the `m` group nodes only; the sink is column `m + 1` of
  the operator and has no index row.

- semantics:

  `"budget"` or `"density"`. Required; there is no default.

- provenance:

  A list carrying at least `method` and `details`, plus `cross_fit` when
  `method` is `"functional"`.

- row_mass:

  Optional declared positive territory measure, one entry per native
  node. Defaults to the unit vector.

- tolerance:

  Positive row-sum tolerance. Group mass above `1 + tolerance` on any
  row is refused rather than renormalized.

## Value

An `effect_location_transport` carrying `$matrix` (the assembled sparse
operator, sink included), `$native_index`, `$group_index`, `$semantics`,
`$row_mass`, `$provenance`, and a content-addressed `$signature` that
binds all of them — including the operator's own contents — so the
transport can enter a scientific plan identity.

## The operator

`matrix` supplies only the `n_native` by `m` **group** block. The
constructor appends the required sink column itself, with row `x`
receiving sink mass `1 - sum_j matrix[x, j]`, so that the assembled
`n_native x (m + 1)` operator is row-stochastic by construction. A row
whose group mass exceeds one would need negative sink mass and is
refused. The sink column always exists, even when it is identically
zero: unmapped territory has to be a number a reader can see, not an
absence.

## Budget and density

`semantics` has no default, because the two choices estimate different
things and the difference is not a display option.

- `"budget"` transports the signed sum: group value `j` is
  `sum_x P[x, j] * values[x]`. Each participant's transported total over
  the group nodes equals its native total minus the sink mass, exactly.

- `"density"` divides that by the transported row mass,
  `sum_x P[x, j] * row_mass[x]`, giving transported budget per unit
  transported territory. A group node reached by no native mass is `NA`,
  never `0`. The sink is reported in budget units under either
  semantics, because a "density" of unmapped territory has no referent.

`row_mass` is the declared positive territory measure of each native
row. Its default is the unit vector — density is then budget per native
node count — and it is recorded either way, because it is the
denominator of the estimand.

## Provenance

`provenance` must be a list carrying `method`, one of `"anatomical"`
(built from a registration or parcellation that never saw the response
data), `"functional"` (built from data), or `"external"` (built outside
`crossform` by a procedure the caller names), and `details`, one string
saying how the operator was built. A `"functional"` transport
additionally requires `cross_fit`: the runs, sessions or tasks whose
data built it. That field is not advisory. A circular transport — one
fitted on a partition that is also used to evaluate it — reports a
transport gain roughly three times the honest one and is otherwise
indistinguishable from it, so the field that would let a plan exclude
those partitions is required rather than encouraged. Any further keys
are kept as declared.

## Refusal

A `"functional"` transport declared without `cross_fit` signals an
`effect_capability_refusal` carrying capability `"cross_fit_provenance"`
in namespace `"location_transport"`. Branch on it with
[`catch_refusal()`](https://bbuchsbaum.github.io/crossform/reference/catch_refusal.md)
rather than on the message text.

## See also

[`anatomical_transport()`](https://bbuchsbaum.github.io/crossform/reference/anatomical_transport.md)
and
[`external_transport()`](https://bbuchsbaum.github.io/crossform/reference/external_transport.md)
for the two convenience constructors, and
[`transport_values()`](https://bbuchsbaum.github.io/crossform/reference/transport_values.md)
for applying one.

Other population transports:
[`anatomical_transport()`](https://bbuchsbaum.github.io/crossform/reference/anatomical_transport.md),
[`estimate_population()`](https://bbuchsbaum.github.io/crossform/reference/estimate_population.md),
[`external_transport()`](https://bbuchsbaum.github.io/crossform/reference/external_transport.md),
[`heterogeneity()`](https://bbuchsbaum.github.io/crossform/reference/heterogeneity.md),
[`materialize_population()`](https://bbuchsbaum.github.io/crossform/reference/materialize_population.md),
[`plan_population()`](https://bbuchsbaum.github.io/crossform/reference/plan_population.md),
[`population_prevalence()`](https://bbuchsbaum.github.io/crossform/reference/population_prevalence.md),
[`population_uncertainty()`](https://bbuchsbaum.github.io/crossform/reference/population_uncertainty.md),
[`population_views`](https://bbuchsbaum.github.io/crossform/reference/population_views.md),
[`transport_values()`](https://bbuchsbaum.github.io/crossform/reference/transport_values.md)

## Examples

``` r
# Four native nodes onto two group nodes. Node 3 splits its mass; node 4
# keeps 30% of its territory unmapped, so 0.3 lands in the sink.
block <- rbind(c(1, 0), c(1, 0), c(0.6, 0.4), c(0, 0.7))
transport <- location_transport(
  block,
  native_index = c("v1", "v2", "v3", "v4"),
  group_index = c("left", "right"),
  semantics = "budget",
  provenance = list(method = "anatomical", details = "hand-declared")
)
transport
#> <effect_location_transport>
#>   nodes:      4 native -> 2 group + sink
#>   semantics:  budget
#>   sink:       mass 0.3 of 4 rows, 7.5% of territory
#>   provenance: anatomical (cross-fit: none)
#>   built:      hand-declared
#>   signature:  sha256:93779bd18b8b...

# The sink column is materialized and closes every row sum to one.
as.matrix(transport$matrix)
#>      [,1] [,2] [,3]
#> [1,]  1.0  0.0  0.0
#> [2,]  1.0  0.0  0.0
#> [3,]  0.6  0.4  0.0
#> [4,]  0.0  0.7  0.3

# Budget preservation: the transported total, sink included, is the native
# total, for signed ledgers as well as nonnegative ones.
ledger <- c(1.5, -0.5, 2, -1)
carried <- transport_values(transport, ledger)
carried
#>   left  right <sink> 
#>    2.2    0.1   -0.3 
abs(sum(carried) - sum(ledger)) < 1e-12
#> [1] TRUE

# A row that claims more than unit mass would need a negative sink, and is
# refused rather than renormalized.
bad <- try(location_transport(
  rbind(c(0.8, 0.4)), "v1", c("left", "right"),
  semantics = "budget",
  provenance = list(method = "anatomical", details = "hand-declared")
), silent = TRUE)
conditionMessage(attr(bad, "condition"))
#> [1] "Native row 1 carries group mass 1.2, which is more than one, so its sink mass would be negative. Rows are stochastic including the sink, and crossform refuses to renormalize a row rather than silently changing what the transport declares."

# A data-derived transport must name the data that built it.
refusal <- catch_refusal(location_transport(
  block, c("v1", "v2", "v3", "v4"), c("left", "right"),
  semantics = "budget",
  provenance = list(method = "functional", details = "clustered responses")
))
refusal$capability
#> [1] "cross_fit_provenance"
refusal$reasons
#> [1] "cross_fit_partitions_not_declared"
```
