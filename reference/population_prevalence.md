# Count how many participants a transported group value stands on

`population_prevalence()` reports two fractions over the participants of
an
[`estimate_population()`](https://bbuchsbaum.github.io/crossform/reference/estimate_population.md)
result. **Sign prevalence** is, at each group node and query, the
fraction of participants whose transported ledger value exceeds
`threshold`. **Alignment prevalence** is, at each group node, the
fraction of participants whose whole profile across the query bank
points the same way as the rest of the group's.

## Usage

``` r
population_prevalence(
  x,
  query = NULL,
  threshold = 0,
  coverage_floor = NULL,
  term = NULL
)
```

## Arguments

- x:

  An `effect_population_result` from
  [`estimate_population()`](https://bbuchsbaum.github.io/crossform/reference/estimate_population.md).

- query:

  Which queries to count over, by name or position. `NULL` (the default)
  counts over every query in the bank. The alignment inner product is
  taken over exactly the selected coordinates.

- threshold:

  The value a participant's transported ledger must strictly exceed to
  be counted, in the ledger's own units. `0` is the default and counts
  participants with a positive contribution.

- coverage_floor:

  Optional whole number. Group nodes whose smallest
  contributing-participant count falls below it are listed in
  `$coverage$below_floor`. `NULL` (the default) reports the counts and
  marks nothing.

- term:

  Refused. Present so that a caller arriving from
  [`population_uncertainty()`](https://bbuchsbaum.github.io/crossform/reference/population_uncertainty.md)
  gets a reason rather than silence.

## Value

An `effect_population_prevalence`: `$sign` holding `$fraction`,
`$count`, `$resolved` and `$contributing` (each a `node`-by-`query`
matrix); `$alignment` holding `$fraction`, `$count`, `$resolved` and the
inner product it was taken in, one entry per group node; `$coverage`;
`$reference`, the fraction a pure-noise cell reports; `$layer` and
`$reading`, which mark the record as descriptive; `$queries` and
`$query_labels`, the selected rows of the bank it counted over; and the
`$index`, `$ledger`, `$semantics`, `$normalization` and `$receipt` of
the result it read. `as.data.frame(x, measure = )` returns one measure
in long form.

## Details

Both are descriptive, both live on the latent layer of
`design/conservative-geometry-contract.md` section 6 as
`design/population-form-contract.md` section 6.5 applies it, and neither
carries or admits an error bar. The record says so in `$layer`, in
`$reading`, and on every printed line.

## Why this is a latent-layer record

The layer projects nothing — there is no eigenvalue truncation here and
no moved mass to account for — but `value > threshold` is the same kind
of operation and fails in the same direction. It is a per-participant
sign clamp: it discards the magnitude of a crossvalidated estimate and
keeps the sign, which is the part that carries the noise. **A group node
and query at which nothing reproduces reports a fraction near `0.5`, not
near `0`**, because every participant contributes an independent coin
flip. `$reference` carries that number so the comparison is not made
against zero by habit.

The threshold is applied strictly (`>`) and absolutely, in the ledger's
own units, with no relative tolerance. That is the guard every per-node
fraction in the package uses, and whether these guards should take a
relative tolerance is an open contract decision
(`conservative-geometry-contract.md` section 11.4, gap G3) that this
layer inherits rather than settling privately. A cell that is zero only
to within rounding is counted by its rounding.

## There is no inference here, and none is derivable

A count of participants is not an estimate of a population quantity on
this layer. The participants are not a sample from a declared
superpopulation of transports, and the binomial variance of a
thresholded crossvalidated estimate is not the sampling variance of
anything the plan sealed. This function computes no standard error, no
interval and no p-value, the record refuses to carry a field named like
one, and constructing a prevalence leaves the result's own
`$uncertainty` untouched.
[`population_uncertainty()`](https://bbuchsbaum.github.io/crossform/reference/population_uncertainty.md)
is the separate inferential layer, and it carries its own measured
calibration record.

## What alignment is an inner product in

Participant `i`'s profile at a group node is `v_i = H y_i`, with `H` the
packed queries the driver lowers the bank to and `y_i` the participant's
packed transported form. The scored quantity is \\\langle v_i, \bar
v\_{-i}\rangle\\ against the **leave-one-out** mean of the other
participants. The reference excludes the participant being scored on
purpose: the plain group mean puts the positive self term \\\lVert
v_i\rVert^2/n\\ inside every product, so the fraction would be inflated
by construction and would read `1` at `n = 1`. Excluding it makes every
product a cross-participant one, which is the same device
`population-form-v1` sections 6.2 and 7.1 use.

That inner product is the Frobenius inner product of the transported
forms *restricted to the bank's span* exactly when the packed queries
are orthonormal, and a Euclidean inner product in coordinates the bank
chose otherwise. `$alignment$readout_gram_deviation` measures the
departure and `$alignment$frobenius_equivalent` states the verdict; a
two-condition contrast bank is not orthonormal.

## The sink, and coverage

The sink is excluded, on the same grounds
[`heterogeneity()`](https://bbuchsbaum.github.io/crossform/reference/heterogeneity.md)
excludes it: it is unmapped territory reported in budget units at no
location, so "the fraction of participants whose unmapped mass is
positive" is a statement about transport failure and not about a place.
The exclusion is recorded at `$receipt$prevalence$sink_excluded` and the
per-participant sink budget stays readable at `$receipt$sink_budget`.

A prevalence is uninterpretable without the number of participants
behind it. `$coverage$contributing` counts, per node and query, the
participants whose transported value is finite and not exactly zero — a
proxy, read off the shipped values, for `population-form-v1` section
7.5's `group_node_subject_coverage`, which is a property of the
transport operators a result does not carry. `coverage_floor` marks the
nodes below a declared floor; it has no default because section 14.3
records the threshold itself as an open maintainer decision.

## Refusals

Each is an `effect_capability_refusal` in namespace
`"population_prevalence"` (see
[`catch_refusal()`](https://bbuchsbaum.github.io/crossform/reference/catch_refusal.md)).

- `population_participant_values` — a
  [`materialize_population()`](https://bbuchsbaum.github.io/crossform/reference/materialize_population.md)
  result. The streamed complete-form route retains no participant axis,
  by design.

- `participant_term_decomposition` — a non-`NULL` `term`. Participants
  carry one transported ledger, not one per group-model column.

- `comparable_ledger_orientation` — a `"unit_budget"` population in
  which some participant's signed native total is negative, so its
  values are sign-flipped relative to the others and a sign count would
  mix two questions.

One refusal is carried rather than raised. The form-space reading — the
fraction of participants whose transported *form* has positive Frobenius
inner product with the group mean form — needs per-participant
transported forms, and no population result ships them: the query-bank
route retains each participant's value only in the coordinates the
declared bank names, and the streamed route retains no participant axis.
Capability `participant_form_prevalence` travels on
`$receipt$form_prevalence_refusal` with its reasons and remedies.

## References

`design/population-form-contract.md` (`population-form-v1`), sections
6.5, 7.5 and 8.1; `design/conservative-geometry-contract.md`
(`conservative-geometry-v1`), section 6.

## See also

[`latent_geometry()`](https://bbuchsbaum.github.io/crossform/reference/latent_geometry.md)
for the projection layer this one shares its discipline with,
[`population_uncertainty()`](https://bbuchsbaum.github.io/crossform/reference/population_uncertainty.md)
for the separate inferential layer, and
[`heterogeneity()`](https://bbuchsbaum.github.io/crossform/reference/heterogeneity.md)
for the form-space reading of how participants differ.

Other population transports:
[`anatomical_transport()`](https://bbuchsbaum.github.io/crossform/reference/anatomical_transport.md),
[`estimate_population()`](https://bbuchsbaum.github.io/crossform/reference/estimate_population.md),
[`external_transport()`](https://bbuchsbaum.github.io/crossform/reference/external_transport.md),
[`heterogeneity()`](https://bbuchsbaum.github.io/crossform/reference/heterogeneity.md),
[`location_transport()`](https://bbuchsbaum.github.io/crossform/reference/location_transport.md),
[`materialize_population()`](https://bbuchsbaum.github.io/crossform/reference/materialize_population.md),
[`plan_population()`](https://bbuchsbaum.github.io/crossform/reference/plan_population.md),
[`population_uncertainty()`](https://bbuchsbaum.github.io/crossform/reference/population_uncertainty.md),
[`population_views`](https://bbuchsbaum.github.io/crossform/reference/population_views.md),
[`transport_values()`](https://bbuchsbaum.github.io/crossform/reference/transport_values.md)

## Examples

``` r
# Six participants on different native frames, carrying one planted effect
# every run of every participant shares (face over house) and a phase-shifted
# pattern that reproduces in none of them. Two group nodes: one every
# participant reaches and one only the larger native frames do.
effects <- effect_space(c("face", "house", "tool"), basis_id = "popprev:v1")
subject <- function(id, n, phase) {
  domain <- abstract_domain(n, coordinates = cbind(x = seq_len(n) - 1),
    feature_ids = paste0("f", seq_len(n)), id = id)
  planted <- outer(c(face = 0.8, house = -0.8, tool = 0), rep(1, n))
  values <- function(shift) matrix(cos(seq_len(3 * n) + shift), 3, n,
    dimnames = list(c("face", "house", "tool"), NULL)) + planted
  rel <- relation(list(run1 = values(0), run2 = values(phase)),
    effects = effects, domain = domain)
  plan_geometry(rel, compile_frame(voxelwise(), domain),
    cross_partitions(rel))
}
carrier <- function(n) anatomical_transport(
  native_coords = cbind(seq_len(n) - 1), group_coords = cbind(c(0, 6)),
  semantics = "budget", radius = 2
)
sizes <- c(s01 = 4L, s02 = 5L, s03 = 6L, s04 = 7L, s05 = 8L, s06 = 9L)
phases <- c(s01 = 0.2, s02 = 1.1, s03 = 2.4, s04 = 0.6, s05 = 3, s06 = 1.7)
subjects <- stats::setNames(lapply(names(sizes), function(id)
  subject(id, sizes[[id]], phases[[id]])), names(sizes))
plan <- plan_population(subjects, lapply(sizes, carrier))
fit <- estimate_population(plan,
  rbind(`face-house` = c(1, -1, 0), `face-tool` = c(1, 0, -1)))

shared <- population_prevalence(fit)
shared
#> <effect_population_prevalence>
#>   layer:        latent descriptive (fractions over participants)
#>   ledger:       transported_total (component "total")
#>   participants: 6 (s01, s02, s03, s04 (+2 more))
#>   group nodes:  2 (sink excluded)
#>   queries:      face-house, face-tool
#>   threshold:    value > 0 (ledger units, exact)
#>   sign:         median 0.667 (range 0.5 to 1)
#>   alignment:    median 0.917 (range 0.833 to 1), leave-one-out reference
#>   readout:      2 queries; bank Gram off identity by 3, not Frobenius
#>   coverage:     minimum 5 of 6 participants contributing; no floor declared
#>   frame:        undeclared, conservative
#>   transport:    budget, anatomical, cross-fit not declared
#>   estimand:     population-sha256:7b32ff9246f3...
#>   reading:      latent descriptive layer; not for inference
#>   a cell at which nothing reproduces reports a fraction near 0.5, not near
#>     0: thresholding a signed crossvalidated estimate keeps the sign and
#>     discards the magnitude, and the sign is the noisy part.
#>   descriptive only. No standard error, interval or p-value is attached to a
#>     count of participants, and none follows from one;
#>     population_uncertainty() is the separate inferential layer.
#>   next:         as.data.frame(x, measure), x$coverage, x$sign$resolved

# The planted contrast is carried by every participant; the one that
# reproduces in none of them sits at the 0.5 the print warns about.
shared$sign$fraction
#>         query
#> node     face-house face-tool
#>   group1  1.0000000       0.5
#>   group2  0.8333333       0.5

# And how many participants were there to be counted at each node.
shared$coverage$contributing
#>         query
#> node     face-house face-tool
#>   group1          6         6
#>   group2          5         5

# A term axis does not exist on participant values, and asking says why.
catch_refusal(population_prevalence(fit, term = "(Intercept)"))$capability
#> [1] "participant_term_decomposition"
```
