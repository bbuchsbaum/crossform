# Population form: one participant's ledger, carried to a group

**This layer is experimental.** The API and the returned kinds may
change, and nothing in it is calibrated group inference. What it *is* is
a typed way to carry each participant’s conservative attribution ledger
onto a shared set of group nodes, fit one group model there, and read
the result without any of the three steps quietly changing what was
estimated.

Read
[`vignette("conservative-frames")`](https://bbuchsbaum.github.io/crossform/articles/conservative-frames.md)
first. Everything here operates on the attribution instrument from that
article: a map whose values are shares of one fixed budget, and which
therefore *adds*. Adding is what makes a group-level ledger definable at
all, and it is why this layer refuses a detection map.

Every identity below is asserted in code you can see. If one stopped
holding, this article would stop knitting. The governing document is
[`design/population-form-contract.md`](https://github.com/bbuchsbaum/crossform/blob/main/design/population-form-contract.md)
(`population-form-v1`), and section numbers refer to it.

## 1. The three objects

### 1.1 Six participants, six conservative geometry plans

The fixture is generated, tiny, and deterministic. Six participants on
**different native frame sizes** — that is the case §3 says forces
transport to precede the fit, and a fixture where everyone shared a
frame would never exercise it. Four runs each, because the cross-fitted
subject Gram of section 4 needs two disjoint halves and each half needs
a cross-generalized estimate of its own.

Every participant carries the same planted **consensus** direction in
effect space (face above house, present in every run, so it survives the
cross-partition product). One participant, `s06`, additionally carries a
direction nobody else has. Section 4 is about finding it without being
told.

``` r

pop_effects <- effect_space(c("face", "house", "tool"),
  basis_id = "population-vignette:v1")

pop_subject <- function(id, features, gain = 1, tilt = 0, runs = 4L,
                        normalization = "conservative") {
  domain <- abstract_domain(features,
    coordinates = cbind(x = seq_len(features) - 1),
    feature_ids = paste0("f", seq_len(features)), id = id)
  consensus <- outer(c(0.9, -0.9, 0), rep(1, features))
  odd <- outer(c(0, 0.8, -0.8), rep(1, features))
  block <- function(k) {
    set.seed(1000L * k + sum(as.integer(charToRaw(id))) + features)
    values <- matrix(gain * stats::rnorm(3L * features), 3L, features,
      dimnames = list(c("face", "house", "tool"), NULL))
    values + consensus + tilt * odd
  }
  relation <- relation(
    stats::setNames(lapply(seq_len(runs), block), paste0("run", seq_len(runs))),
    effects = pop_effects, domain = domain
  )
  plan_geometry(relation,
    compile_frame(voxelwise(normalization = normalization), domain),
    cross_partitions(relation))
}

sizes <- c(s01 = 9L, s02 = 11L, s03 = 13L, s04 = 10L, s05 = 12L, s06 = 14L)
gains <- c(s01 = 1, s02 = 1.3, s03 = 0.8, s04 = 1.1, s05 = 0.9, s06 = 1.2)
tilts <- c(s01 = 0, s02 = 0, s03 = 0, s04 = 0, s05 = 0, s06 = 1.6)

subjects <- stats::setNames(lapply(names(sizes), function(id)
  pop_subject(id, sizes[[id]], gains[[id]], tilts[[id]])), names(sizes))
sizes
#> s01 s02 s03 s04 s05 s06 
#>   9  11  13  10  12  14
```

These are ordinary single-participant objects:
[`plan_geometry()`](https://bbuchsbaum.github.io/crossform/reference/plan_geometry.md)
sealed six estimands and read no data.
[`voxelwise()`](https://bbuchsbaum.github.io/crossform/reference/voxelwise.md)
is conservative by default, which is the gate section 1.3 enforces.

### 1.2 A transport is an input, and it is typed

A **transport** is the operator that carries one participant’s native
nodes onto the shared group nodes. `crossform` accepts one and refuses
to learn one: image registration and functional-transport fitting stay
outside the package (§9.2). What it does is check the operator’s shape,
record its provenance, and put both into the estimand’s identity — two
transports are two estimands (§1.5).

The cheapest typed transport is built from coordinates. Group centres at
`x = 0, 5, 11`, a radius of 2, and each native node assigned to its
nearest centre within that radius:

``` r

pop_carrier <- function(features, semantics = "budget", radius = 2) {
  anatomical_transport(
    native_coords = cbind(seq_len(features) - 1),
    group_coords = cbind(c(0, 5, 11)),
    semantics = semantics, radius = radius,
    native_index = paste0("f", seq_len(features))
  )
}
transports <- stats::setNames(lapply(names(sizes), function(id)
  pop_carrier(sizes[[id]])), names(sizes))
transports$s01
#> <effect_location_transport>
#>   nodes:      9 native -> 3 group + sink
#>   semantics:  budget
#>   sink:       mass 1 of 9 rows, 11.1% of territory
#>   provenance: anatomical (cross-fit: none)
#>   built:      nearest group centre within radius 2, ties to the lowest gr...
#>   signature:  sha256:c8881705fe27...
```

Three things in that print are load-bearing.

**The sink is a column, always.** The operator is `n × (m + 1)`, not
`n × m`: group nodes plus one accounting column for native mass that
reached no group node. It is materialized even when it is empty, so
partial coverage is a number you read rather than budget that quietly
went missing.

**`semantics` has no default.** `"budget"` and `"density"` are two
different estimands, and the constructor makes you name which.

**Provenance is required.** `method` is a closed set, and `details` must
say how the operator was built.

The same laws hold for an operator you computed elsewhere and are
declaring. Here is a four-node example small enough to check by hand:
`f2` splits 0.6/0.4, `f3` leaves 30 % of its territory unmapped, `f4` is
unmapped entirely.

``` r

P <- rbind(
  f1 = c(anterior = 1,   posterior = 0),
  f2 = c(anterior = 0.6, posterior = 0.4),
  f3 = c(anterior = 0,   posterior = 0.7),
  f4 = c(anterior = 0,   posterior = 0)
)
carrier <- external_transport(P, semantics = "budget",
  provenance = list(details = "partial-volume warp, atlas-tool 2.1"))
c(native = nrow(carrier$matrix), columns = ncol(carrier$matrix))
#>  native columns 
#>       4       3
stopifnot(
  identical(dim(carrier$matrix), c(4L, 3L)),
  identical(nrow(carrier$group_index), 2L),
  max(abs(Matrix::rowSums(carrier$matrix) - 1)) < 1e-12
)
```

Now carry a **signed** ledger through it — signed because a
crossvalidated conservative ledger is signed by construction, and a
fixture of nonnegative numbers would let a mass-law bug pass.

``` r

ledger <- c(f1 = 1.5, f2 = -0.5, f3 = 2, f4 = 0.25)
carried <- transport_values(carrier, ledger)
round(carried, 4)
#>  anterior posterior    <sink> 
#>      1.20      1.20      0.85

c(
  closes = sum(carried) - sum(ledger),
  without_sink = sum(carried[c("anterior", "posterior")]) - sum(ledger)
)
#>       closes without_sink 
#>         0.00        -0.85
stopifnot(
  abs(sum(carried) - sum(ledger)) <= 1e-12 * sum(abs(ledger)),
  abs(sum(carried[c("anterior", "posterior")]) - sum(ledger)) > 0.8
)
```

Budget semantics preserve the total **exactly, including the sink** —
and lose 0.85 units of signed mass the moment the sink is deleted. That
is why it is required rather than optional.

The same operator read as a **density** is a different estimand: each
group column is divided by the row mass that reached it, so the answer
is a value per unit of declared territory and satisfies no conservation
law at all.

``` r

dense <- external_transport(P, semantics = "density",
  provenance = list(details = "the same operator, read as a density"))
carried_density <- transport_values(dense, ledger)
round(carried_density, 4)
#>  anterior posterior    <sink> 
#>    0.7500    1.0909    0.8500
c(density_gap = sum(carried_density) - sum(ledger))
#> density_gap 
#>  -0.5590909
stopifnot(abs(sum(carried_density) - sum(ledger)) > 1e-6)
```

That gap is arithmetic, not a discovery. Density trades conservation
away; what it keeps is section 2’s commutation, because a declared
row-mass ratio is still a fixed linear map.

**A transport fitted to the response data is refused unless it says
which partitions built it.** An operator that saw the same runs the
geometry is estimated from is circular, and the contract measures what
that circularity buys (§7): the refusal is a gate, not a formality.

``` r

refusal <- catch_refusal(external_transport(P, semantics = "budget",
  provenance = list(method = "functional",
    details = "hyperalignment on the analysis runs")))
c(capability = refusal$capability, reason = refusal$reasons)
#>                          capability                              reason 
#>              "cross_fit_provenance" "cross_fit_partitions_not_declared"
stopifnot(
  inherits(refusal, "effect_capability_refusal"),
  identical(refusal$capability, "cross_fit_provenance"),
  identical(refusal$namespace, "location_transport"),
  identical(refusal$reasons, "cross_fit_partitions_not_declared")
)

honest <- external_transport(P, semantics = "budget",
  provenance = list(method = "functional", details = "hyperalignment",
    cross_fit = c("run1", "run2")))
honest$provenance$cross_fit
#> [1] "run1" "run2"
stopifnot(identical(honest$provenance$cross_fit, c("run1", "run2")))
```

Naming the partitions is what the refusal asks for, and the names travel
with the operator into the estimand’s identity. An anatomical or
external transport never saw the responses, so it owes no such record —
which is why the six carriers above passed without one.

### 1.3 `plan_population()`: the group estimand

``` r

plan <- plan_population(subjects, transports)
plan
#> <effect_population_plan>
#>   subjects:      6 (s01, s02, s03, s04 (+2 more))
#>   group nodes:   3 + sink
#>   sink:          present in 6 of 6 subjects, worst 11.1% of territory
#>   transport:     budget, anatomical
#>   model:         ~1 -> 1 column, rank 1
#>   normalization: none
#>   fit:           OLS (subject-constant weights), transport then fit
#>   estimand:      population-sha256:085f2f7724df...
#>   signature:     sha256:28ffd27995a6...
```

``` r

c(
  semantics = plan$semantics,
  normalization = plan$normalization,
  order = plan$fit$evaluation_order
)
#>            semantics        normalization                order 
#>             "budget"               "none" "transport_then_fit"
plan$subject_index[, c("subject", "measurements", "declared_normalization",
  "conserved", "sink_territory")]
#>   subject measurements declared_normalization conserved sink_territory
#> 1     s01            9           conservative      TRUE     0.11111111
#> 2     s02           11           conservative      TRUE     0.09090909
#> 3     s03           13           conservative      TRUE     0.07692308
#> 4     s04           10           conservative      TRUE     0.10000000
#> 5     s05           12           conservative      TRUE     0.08333333
#> 6     s06           14           conservative      TRUE     0.07142857
stopifnot(
  all(plan$subject_index$conserved),
  identical(plan$semantics, "budget"),
  identical(plan$normalization, "none"),
  identical(plan$fit$evaluation_order, "transport_then_fit"),
  identical(plan$fit$commuting, TRUE),
  identical(nrow(plan$group_index), 3L)
)
```

The evaluation order is **recorded, not chosen**. Participants have
different native frames, so there is no common node axis to fit at
before transporting: fit-then-transport is not definable here at all.
Section 2 is about the orders that *are* definable.

The conservative gate is the plan’s headline refusal. A participant
whose frame reports a density has no budget to partition, and the plan
says so by name:

``` r

loose <- subjects
loose$s01 <- pop_subject("s01", sizes[["s01"]], gains[["s01"]],
  normalization = "local")
gate <- catch_refusal(plan_population(loose, transports))
c(capability = gate$capability, reason = gate$reasons[[1L]])
#>                                 capability 
#>            "conservative_subject_geometry" 
#>                                     reason 
#> "normalization_not_conservative:s01:local"
stopifnot(
  identical(gate$capability, "conservative_subject_geometry"),
  "normalization_not_conservative:s01:local" %in% gate$reasons,
  any(grepl("allow_nonconservative", gate$remedies))
)
```

The escape hatch exists, is named in the remedy, and is *recorded* — it
enters the plan’s scientific identity, because a non-conservative
population estimand is a different estimand and not a relaxed setting on
the same one.

The group model is the last piece. It is one-sided, its rows bind to
participants **by name** rather than by position, and it is factorized
once on the plan:

``` r

covariates <- data.frame(
  age = c(24, 31, 27, 44, 38, 22),
  row.names = c("s02", "s01", "s04", "s03", "s06", "s05")
)
aged <- plan_population(subjects, transports, model = ~ age, data = covariates)
aged$model$matrix
#>     (Intercept) age
#> s01           1  31
#> s02           1  24
#> s03           1  44
#> s04           1  27
#> s05           1  22
#> s06           1  38
#> attr(,"assign")
#> [1] 0 1
stopifnot(
  identical(aged$model$columns, c("(Intercept)", "age")),
  identical(aged$model$rank, 2L),
  identical(rownames(aged$model$matrix), names(sizes)),
  # The shuffled `data` rows followed their participant, not their position.
  identical(unname(aged$model$matrix[, "age"]), c(31, 24, 44, 27, 22, 38)),
  # A different group model is a different estimand.
  !identical(plan$scientific_plan_id, aged$scientific_plan_id)
)
```

The rest of this article uses the intercept-only `plan`, whose single
term is the group mean.

## 2. Model then query is query then model

This is the claim the layer rests on (§3). The query mixes
**experimental** coordinates, the transport mixes **spatial** nodes, and
the group fit mixes **participants**. Three different axes of the same
array, so under the OLS default with subject-constant weights they
commute — and the answer does not depend on the order you ran them in.

[`estimate_population()`](https://bbuchsbaum.github.io/crossform/reference/estimate_population.md)
runs `query → transport → fit`: it contracts each participant’s geometry
against the bank first and carries **one number per node per query**.

``` r

bank <- rbind(`face-house` = c(1, -1, 0), `house-tool` = c(0, 1, -1))
fit <- estimate_population(plan, bank)
dim(fit$coefficients)
#> [1] 4 2 1
dimnames(fit$coefficients)[c("query", "term")]
#> $query
#> [1] "face-house" "house-tool"
#> 
#> $term
#> [1] "(Intercept)"
fit$index
#>     node coord1  sink  units
#> 1 group1      0 FALSE budget
#> 2 group2      5 FALSE budget
#> 3 group3     11 FALSE budget
#> 4 <sink>     NA  TRUE budget
```

### The check that is not a check

Before the real one, the trap. Comparing
[`transport_values()`](https://bbuchsbaum.github.io/crossform/reference/transport_values.md)
applied to a participant’s own
[`contrast_energy()`](https://bbuchsbaum.github.io/crossform/reference/contrast_energy.md)
against `fit$values` looks like a commutation test and is not:
[`contrast_energy()`](https://bbuchsbaum.github.io/crossform/reference/contrast_energy.md)
lowers to exactly the packed operator the bank carries, so both sides
are the *same call on the same numbers*. It is worth recording as what
it is — evidence that the executor’s plumbing does not rescale or
reorder rows — and it is exactly zero, which is the tell.

``` r

one <- estimate_population(plan, rbind(`face-house` = c(1, -1, 0)))
plumbing <- max(vapply(names(plan$subjects), function(id) {
  native <- contrast_energy(plan$subjects[[id]], c(1, -1, 0))$total
  max(abs(transport_values(plan$transport[[id]], native) -
    one$values[, "face-house", id]))
}, numeric(1)))
c(executor_plumbing = plumbing)
#> executor_plumbing 
#>                 0
stopifnot(plumbing == 0)
```

### The real order reversal

[`materialize_population()`](https://bbuchsbaum.github.io/crossform/reference/materialize_population.md)
runs the *other* order. It carries all six packed coordinates of the
complete `3 × 3` form to the group nodes, fits the group model at every
one of them, and hands back a coefficient **form** per node. The
[`rdm()`](https://bbuchsbaum.github.io/crossform/reference/rdm.md) view
then contracts the (face, house) edge out of those coefficients —
**afterwards**. That is `transport → fit → query`.

The face − house contrast *is* that RDM edge. Two genuinely different
sequences of arithmetic, one number.

``` r

form <- materialize_population(plan)
edges <- as.data.frame(rdm(form))
edge <- edges[edges$left == "face" & edges$right == "house", ]
node_ids <- as.character(fit$index$node)

# One row per node: with a richer group model `rdm()` emits one row per node
# per term, and `match()` would silently take the first.
stopifnot(nrow(edge) == length(node_ids), !anyDuplicated(edge$node))

commutation <- max(abs(
  edge$estimate[match(node_ids, edge$node)] -
    fit$coefficients[node_ids, "face-house", "(Intercept)"]
))
c(query_first = fit$basis, transport_first = form$basis)
#>     query_first transport_first 
#>    "query_bank" "complete_form"
c(commutation = commutation)
#>  commutation 
#> 5.329071e-15
stopifnot(
  identical(fit$basis, "query_bank"),
  identical(form$basis, "complete_form"),
  commutation < 1e-12
)
```

`1e-12` is the tolerance §11 states for the commutation of query,
transport and OLS under subject-constant weights. Unlike the plumbing
check this one lands at a rounding scale rather than at exactly zero,
because the two routes really are different arithmetic: one contraction
before transport against six after it.

This is what licenses the reader verbs. A population view is a fixed
linear combination of estimated query columns, not a second execution —
so
[`contrast_energy()`](https://bbuchsbaum.github.io/crossform/reference/contrast_energy.md),
[`rdm()`](https://bbuchsbaum.github.io/crossform/reference/rdm.md),
[`rsa()`](https://bbuchsbaum.github.io/crossform/reference/rsa.md) and
[`contribution()`](https://bbuchsbaum.github.io/crossform/reference/contribution.md)
on a result reopen no participant’s geometry, and the combination is
exact rather than approximate.

Reading the query first is also what makes the route cheap: complete
geometry at a native node is `q(q+1)/2` packed coordinates, a bank of
`K` queries is `K` numbers, and complete geometry is never allocated.

## 3. Conservation survives the transport, sink included

Section 2 of the contract is the certificate: under budget semantics
each participant’s transported total, sink included, equals its native
total exactly.
[`estimate_population()`](https://bbuchsbaum.github.io/crossform/reference/estimate_population.md)
asserts it per participant and per query at fit time, with a tolerance
of `1e-12` times the ledger’s **L1 norm** — not its total, which a
signed ledger can drive near zero.

``` r

native_total <- fit$receipt$native_total
carried_total <- apply(fit$values, c("subject", "query"), sum)
l1 <- apply(abs(fit$values), c("subject", "query"), sum)
budget_deviation <- max(abs(carried_total - native_total) / l1)

fit$receipt$budget$asserted
#> [1] TRUE
c(
  recorded = fit$receipt$budget$max_relative_deviation,
  recomputed = budget_deviation
)
#>     recorded   recomputed 
#> 1.984263e-16 2.059721e-16
stopifnot(
  isTRUE(fit$receipt$budget$asserted),
  identical(fit$receipt$budget$scale, "relative_to_ledger_l1_norm"),
  fit$receipt$budget$max_relative_deviation < 1e-12,
  budget_deviation < 1e-12
)
```

The sink is not a discarded remainder. Drop it and the identity breaks,
by an amount that is a real number about a real failure of coverage:

``` r

group_only <- apply(fit$values[!fit$index$sink, , , drop = FALSE],
  c("subject", "query"), sum)
without_sink <- max(abs(group_only - native_total) / l1)
c(
  with_sink = budget_deviation,
  without_sink = without_sink,
  sink_territory = max(plan$subject_index$sink_territory)
)
#>      with_sink   without_sink sink_territory 
#>   2.059721e-16   3.271382e-01   1.111111e-01
stopifnot(without_sink > 1e-3, all(plan$subject_index$sink_territory > 0))
```

Between 7 % and 11 % of each participant’s native territory here falls
outside every group node’s radius. The sink is also **fitted like any
other row** — its coefficients are estimated rather than dropped —
because covariate-linked sink mass is how a differentially failing
transport becomes visible (§3.3). A transport that works worse in older
participants shows up as an age effect in the sink, and only if you left
the sink in the model.

Note the units discipline in `fit$index`: the sink is reported in
**budget** units under either semantics, so a density population still
has one honest column of unmapped mass.

``` r

fit$index$units
#> [1] "budget" "budget" "budget" "budget"
stopifnot(all(fit$index$units == "budget"))
```

## 4. Heterogeneity: what the participants disagree about

The group fit gives the consensus. The other half of the section 6
decomposition is the scatter: the `N × N` Gram of participants’
deviations from the group fit, in the packed geometry coordinates.
Geometry-space covariance across `N` participants has rank at most
`N − 1`, so the whole spectrum lives in that small matrix and the
`D × D` covariance never has to be formed.

``` r

het <- heterogeneity(plan, estimator = "cross_fit")
c(estimator = het$estimator, space = het$space)
#>     estimator         space 
#>   "cross_fit" "packed_form"
c(residual_df = het$residual_df)
#> residual_df 
#>           5
round(het$spectrum, 3)
#>   mode1   mode2   mode3   mode4   mode5   mode6 
#> 133.865   9.184   0.926   0.000  -6.809 -14.851
stopifnot(
  identical(het$estimator, "cross_fit"),
  identical(het$space, "packed_form"),
  identical(het$residual_df, 5L)
)
```

**The Gram is indefinite, and that is reported rather than repaired.**
The cross-fitted estimator subtracts a cross-participant term to remove
the consensus-squared part, and nothing constrains the difference of two
estimates to be positive semidefinite. The contract measures it in 100 %
of 2 000 replications (§6.2); here it shows up as two genuinely negative
modes.

``` r

c(
  negative_modes = het$latent$negative_modes,
  moved_share = round(het$latent$moved_share, 4),
  n_eff = round(het$latent$n_eff, 3)
)
#> negative_modes    moved_share          n_eff 
#>         2.0000         0.1308         1.1510
stopifnot(
  sum(het$spectrum < -1e-8) == 2L,
  identical(het$latent$negative_modes, 2L),
  het$latent$moved_share > 0
)
```

The effective mode count and the cumulative curve are legal only on the
declared PSD projection, which records the mass it moved — the same
latent-layer discipline
[`latent_geometry()`](https://bbuchsbaum.github.io/crossform/reference/latent_geometry.md)
applies to a single participant’s map.

### The planted participant, found

The loadings are the eigenvectors of that Gram: one number per
participant per mode. Nothing above told the estimator that `s06` is
different.

``` r

round(het$loadings[, 1:2], 3)
#>      mode1  mode2
#> s01 -0.217  0.382
#> s02 -0.153  0.641
#> s03 -0.161 -0.189
#> s04 -0.178 -0.592
#> s05 -0.203 -0.238
#> s06  0.911 -0.004
leading <- het$loadings[, 1L]
others <- leading[names(leading) != "s06"]
names(which.max(abs(leading)))
#> [1] "s06"
c(
  separation = round(max(abs(leading)) / max(abs(others)), 2),
  mode1_share = round(het$latent$cumulative[[1L]], 3)
)
#>  separation mode1_share 
#>        4.20        0.93
stopifnot(
  identical(names(which.max(abs(leading))), "s06"),
  max(abs(leading)) > 4 * max(abs(others)),
  het$latent$cumulative[[1L]] > 0.85
)
```

Mode 1 is one participant against the rest — `s06` at 0.91 with
everybody else between −0.15 and −0.22, the signature of a single odd
participant rather than a gradient — and it holds over 90 % of the
projected heterogeneity mass. Pass `nodes =` to reconstruct the
*geometry* of a mode at named group nodes and see which effect pair it
lives in.

### Why the estimator is not a free choice

The plug-in Gram books every participant’s within-subject sampling noise
as between-participant heterogeneity. The contract measures the
consequence: over 2 000 Monte Carlo replications at `N = 12`, the
plug-in inflates the heterogeneity trace by **+62.7 %** (mean `tr Q^H`
of `+17.61` against a true `+10.82`, predicted bias `+6.72` versus
measured `+6.79`), while the cross-fitted estimator’s bias is within one
Monte Carlo standard error of zero (§6.2). That number is a property of
the estimator, not of this fixture; one draw of a trace estimates
nothing. What *this* fixture can show is the direction, and that the two
Grams disagree about size while agreeing about direction.

``` r

plug <- heterogeneity(plan, estimator = "plug_in")
c(
  cross_fit_trace = round(sum(diag(het$gram)), 3),
  plug_in_trace = round(sum(diag(plug$gram)), 3),
  mode1_agreement = round(abs(sum(het$loadings[, 1L] * plug$loadings[, 1L])), 4)
)
#> cross_fit_trace   plug_in_trace mode1_agreement 
#>        122.3140        157.3450          0.9987
stopifnot(
  sum(diag(plug$gram)) > sum(diag(het$gram)),
  abs(sum(het$loadings[, 1L] * plug$loadings[, 1L])) > 0.99,
  identical(names(which.max(abs(plug$loadings[, 1L]))), "s06")
)
```

Section 6.3 makes the split normative: **loading directions may be read
from either Gram provided the source is named; any eigenvalue, spectrum,
variance-explained figure or effective mode count must come from the
cross-fitted Gram.** The package enforces it rather than advising it —
the plug-in record carries no latent layer, and says why.

``` r

plug_refusal <- plug$receipt$latent_refusal
c(capability = plug_refusal$capability, reasons = plug_refusal$reasons)
#>                                                    capability 
#>                                "plug_in_spectrum_functionals" 
#>                                                      reasons1 
#>                "within_subject_noise_booked_as_heterogeneity" 
#>                                                      reasons2 
#> "plug_in_trace_inflated_62_7_percent_on_the_contract_fixture"
stopifnot(
  is.null(plug$latent),
  identical(plug_refusal$capability, "plug_in_spectrum_functionals"),
  "plug_in_trace_inflated_62_7_percent_on_the_contract_fixture" %in%
    plug_refusal$reasons
)
```

## 5. Two error bars, kept apart — and a count that is neither

### 5.1 The between-subject layer

The scatter of the participants about the group fit. For the
intercept-only model this has a closed form nobody needs a package for,
which is exactly why it is worth checking against one:

``` r

uncertainty <- population_uncertainty(fit)
between <- uncertainty$between

hand_estimate <- apply(fit$values, c("node", "query"), mean)
hand_se <- apply(fit$values, c("node", "query"),
  function(y) stats::sd(y) / sqrt(length(y)))
c(
  residual_df = between$residual_df,
  estimate_gap = max(abs(between$estimate[, , "(Intercept)"] - hand_estimate)),
  se_gap = max(abs(between$se[, , "(Intercept)"] - hand_se))
)
#>  residual_df estimate_gap       se_gap 
#> 5.000000e+00 7.105427e-15 4.440892e-16
stopifnot(
  identical(between$residual_df, 5L),
  max(abs(between$estimate[, , "(Intercept)"] - hand_estimate)) < 1e-12,
  max(abs(between$se[, , "(Intercept)"] - hand_se)) < 1e-12
)
```

**The `t` is labelled uncalibrated, and the label does not move.** The
arithmetic has been checked against a 2 000-replication null simulation
and it is right: under a *correctly specified* group model the nominal
95 % interval covered the null term in 0.9485 of replications at `N = 6`
and 0.9520 at `N = 12`. The arithmetic was never the part in doubt. The
second arm of the same simulation is why the label stays: when each
participant’s transported value carries a variance that depends on the
group covariates — which is what a transport whose quality varies with
age or motion produces — coverage falls to 0.9230 at `N = 6` and
**0.8850 at `N = 24`**. It gets *worse* with more participants, because
the bias is in the standard error and not in the sample size.

``` r

between$calibration
#> [1] "uncalibrated"
c(level = between$level, t_max = round(max(abs(between$t), na.rm = TRUE), 3))
#> level t_max 
#> 0.950 5.321
stopifnot(identical(between$calibration, "uncalibrated"))
```

Report the statistic; do not report a p-value derived from it without an
argument that the section 7.5 transport diagnostics are benign.

### 5.2 The within-subject layer, and where it is refused

The second bar is the sampling variance of one participant’s transported
value. It exists **only where it is exact**. A transported value is a
fixed linear functional `w'z` of that participant’s native query values,
so its variance needs the covariance *between* native nodes — an object
the package refuses (capability `cross_node_sampling_covariance`). The
layer is therefore admitted exactly where `w` has one nonzero entry:
there the cross terms carry weight zero and `Var = w²Var(z)` is exact,
with **no independence assumption anywhere**.

That needs participants with an error channel, so this is a separate
small fixture arriving through the ingestion route. Group nodes sit on
the native grid itself, so a four-node participant maps one-to-one; the
five-node participant sends two native rows into the last group node.

``` r

wu_subject <- function(id, features, seed) {
  set.seed(seed)
  domain <- abstract_domain(features,
    coordinates = cbind(x = seq_len(features) - 1),
    feature_ids = paste0("f", seq_len(features)), id = id)
  runs <- c("run-1", "run-2", "run-3")
  scans <- function(run) paste0(run, "-scan-", 1:6)
  design <- cbind(face = c(1, 0, 0, 1, 0, 0), house = c(0, 1, 0, 0, 1, 0),
    tool = c(0, 0, 1, 0, 0, 1))
  target <- diag(3)
  dimnames(target) <- list(c("face", "house", "tool"), colnames(design))
  relation_fit <- estimate_relation(plan_relation(
    study(observations(
      stats::setNames(lapply(runs, function(run)
        matrix(stats::rnorm(6 * features), 6L, features)), runs),
      stats::setNames(lapply(runs, function(run)
        observation_index(scans(run), run)), runs),
      domain)),
    raw_design_model(stats::setNames(lapply(runs, function(run) {
      value <- design; rownames(value) <- scans(run); value
    }), runs)),
    raw_effect_map(target),
    observation_model("ols", sampling_unit = "scan")
  ))
  list(fit = relation_fit, features = features,
    plan = plan_geometry(relation_fit$relation,
      compile_frame(voxelwise(), domain),
      cross_partitions(relation_fit$relation, independence = "independent")))
}

built <- list(u01 = wu_subject("u01", 4L, 41L),
  u02 = wu_subject("u02", 4L, 52L), u03 = wu_subject("u03", 5L, 63L))
grid_carrier <- function(features) anatomical_transport(
  native_coords = cbind(seq_len(features) - 1), group_coords = cbind(0:3),
  semantics = "budget", native_index = paste0("f", seq_len(features)))

wplan <- plan_population(lapply(built, `[[`, "plan"),
  lapply(built, function(value) grid_carrier(value$features)))
wbank <- rbind(`face-house` = c(1, -1, 0), `face-tool` = c(1, 0, -1))
wfit <- estimate_population(wplan, wbank,
  uncertainty = lapply(built, function(value)
    rdm_sampling_covariance(value$plan, value$fit, target = "null",
      at = seq_len(value$features))))
within <- population_uncertainty(wfit)$within
```

``` r

within$admitted
#>          u01   u02   u03
#> group1  TRUE  TRUE  TRUE
#> group2  TRUE  TRUE  TRUE
#> group3  TRUE  TRUE  TRUE
#> group4  TRUE  TRUE FALSE
#> <sink> FALSE FALSE FALSE
c(scope = within$scope, assumption = within$assumption)
#>                                          scope 
#>             "transported_single_source_column" 
#>                                     assumption 
#> "none: the cross-node terms carry weight zero"
c(admitted_columns = within$admitted_columns)
#> admitted_columns 
#>               11
stopifnot(
  identical(within$scope, "transported_single_source_column"),
  identical(within$assumption, "none: the cross-node terms carry weight zero"),
  # One-to-one for the four-node participants at every group node; `u03`'s
  # last group node collects two native rows and is therefore absent.
  identical(unname(within$admitted[, "u01"]), c(rep(TRUE, 4L), FALSE)),
  identical(unname(within$admitted[, "u03"]), c(rep(TRUE, 3L), FALSE, FALSE)),
  identical(within$admitted_columns, 11L),
  all(is.na(within$variance["<sink>", , ])),
  all(is.na(within$variance["group4", , "u03"])),
  "transported_value_mixes_native_nodes" %in% within$refusal$reasons
)
```

Where the column mixes native rows there is an **absence**, never a
diagonal approximation standing in for the missing cross terms.
Overlapping supports under spatially correlated noise are positively
correlated, so a diagonal sum would under-estimate the variance in a
known direction — worse than nothing.

The six-participant population of sections 1 to 4 admits nothing at all:
every group node there collects several native nodes. **That is the
refusal being visible, not the layer being broken**, and it is the
ordinary situation for a hard anatomical parcellation.

### 5.3 The two are never pooled

``` r

uncertainty$layers
#> [1] "between_subject" "within_subject"
uncertainty$separation
#> [1] "between_subject and within_subject are reported separately and are never pooled"
stopifnot(
  identical(uncertainty$layers, c("between_subject", "within_subject")),
  # There is no field holding their sum, and `as.data.frame()` emits one
  # layer at a time.
  is.null(uncertainty$total),
  !any(c("total", "pooled", "combined") %in% names(uncertainty)),
  identical(
    sort(unique(as.data.frame(uncertainty, layer = "between")$layer)),
    "between_subject")
)
```

They answer different questions and are not summands. The
between-subject residual already contains whatever measurement error
survived into each participant’s transported value, so adding the
within-subject variance would double-count the shared part *and* still
miss the covariance a variance-components model would need.

### 5.4 Prevalence is descriptive, and says so

[`population_prevalence()`](https://bbuchsbaum.github.io/crossform/reference/population_prevalence.md)
counts how many participants stand behind a group value. It is not the
third error bar; it is on the **latent layer**, the same place
[`latent_geometry()`](https://bbuchsbaum.github.io/crossform/reference/latent_geometry.md)
confines fractions to, and it carries no interval.

``` r

before <- fit$uncertainty
prevalence <- population_prevalence(fit, coverage_floor = length(subjects))
round(prevalence$sign$fraction, 3)
#>         query
#> node     face-house house-tool
#>   group1      0.833      0.833
#>   group2      1.000      1.000
#>   group3      0.833      0.667
c(reference = prevalence$reference)
#> reference 
#>       0.5
c(layer = prevalence$layer, reading = prevalence$reading)
#>                                         layer 
#>                          "latent_descriptive" 
#>                                       reading 
#> "latent descriptive layer; not for inference"
# Nothing on the record may look like an error bar, at any depth.
named <- local({
  flatten <- function(x) {
    if (!is.list(x)) return(character())
    c(names(x), unlist(lapply(x, flatten), use.names = FALSE))
  }
  flatten(prevalence[c("sign", "alignment", "coverage")])
})
stopifnot(
  identical(prevalence$layer, "latent_descriptive"),
  identical(prevalence$reading, "latent descriptive layer; not for inference"),
  identical(prevalence$reference, 0.5),
  identical(intersect(named, c("se", "t", "lower", "upper", "p_value", "p",
    "statistic", "level", "confidence")), character(0)),
  identical(prevalence$receipt$prevalence$inference, "none_derivable"),
  # And constructing one does not reach into the inferential layer.
  identical(before, fit$uncertainty)
)
```

**A group node and query at which nothing reproduces reports a fraction
near `0.5`, not near `0`.** Thresholding a signed crossvalidated
estimate discards its magnitude and keeps its sign, and the sign is the
noisy part; every participant contributes an independent coin flip.
`$reference` carries that number so the comparison is not made against
zero out of habit.

Coverage is reported beside it, because a high prevalence at a node only
a few participants reached is a different object from a high prevalence
at a node they all reached:

``` r

prevalence$coverage$minimum
#> group1 group2 group3 
#>      6      6      5
prevalence$coverage$below_floor
#> [1] "group3"
stopifnot(
  identical(prevalence$coverage$floor, 6L),
  identical(prevalence$coverage$below_floor, "group3"),
  prevalence$coverage$minimum[["group3"]] < length(subjects)
)
```

`group3` sits at `x = 11`, out of radius of every native node the
smallest participant has. Five participants stand behind it, not six — a
fact about the transport, available before any of its numbers are
interpreted.

## 6. What this article does not give you

**No calibrated group inference.** Section 5.1 is the whole story: the
standard error is the group OLS’s own and the `t` is arithmetically
correct, but the recorded null simulation shows the interval losing
coverage in exactly the regime a real population fit occupies —
transport quality varying with the group covariates. Nothing in this
package converts that `t` into a defensible p-value, and
[`population_prevalence()`](https://bbuchsbaum.github.io/crossform/reference/population_prevalence.md)
is descriptive rather than an alternative route to one.

**No transport learning, and no registration.** §9.2 is a list of four
things the package refuses and requires as typed input: image
registration of any kind; functional-transport learning (a `P^F` bearing
cross-fit provenance is *accepted and evaluated*, never fitted);
resampling or interpolation of subject images — transport acts on nodes,
not voxels; and any group frame over group features, and therefore any
group-node coherent component. That last one is why the transported
components take their own names — `native_coherent_ledger` is
native-node coherent evidence carried to a group node, not a group-node
common mode, and reading it as one is the error the name exists to
prevent.

``` r

c(total = fit$ledger,
  coherent = estimate_population(plan, bank, component = "coherent")$ledger)
#>                    total                 coherent 
#>      "transported_total" "native_coherent_ledger"
stopifnot(identical(fit$ledger, "transported_total"))
```

**No cross-node sampling covariance**, and therefore no error bar on a
conserved budget. Node estimates of an overlapping frame are strongly
positively correlated; variances do not add. Section 5.2’s carve-out is
the whole of what exists.

**No precision weighting.** `normalization = "precision_weighted"` is in
the closed set and refused at plan construction, because the per-subject
budget variance it would need is exactly the missing object above.

``` r

weighted <- catch_refusal(plan_population(subjects, transports,
  normalization = "precision_weighted"))
weighted$capability
#> [1] "precision_weighted_normalization"
stopifnot(
  identical(weighted$capability, "precision_weighted_normalization"),
  "per_subject_budget_variance_unavailable" %in% weighted$reasons
)
```

## See also

- [`design/population-form-contract.md`](https://github.com/bbuchsbaum/crossform/blob/main/design/population-form-contract.md)
  — the normative document. §1 the transport object, §2 the budget
  certificate, §3 the commutation claim and the pinned evaluation order,
  §4 the normalizations, §5–6 the subject Gram and the measured plug-in
  bias, §6.5 the latent layer, §7 transport diagnostics and
  `η_transport`, §9.2 the four refusals of section 6 above, §11 every
  tolerance asserted here.
- [Haxby 2001
  exemplar](https://github.com/bbuchsbaum/crossform/tree/main/exemplars/haxby2001),
  script `10-population-slice1.R` and **§10 of its README** — this
  article’s identities on six real participants at region-level nodes,
  with committed receipts. Read its note on what an exact zero would
  mean: the executor-plumbing check of section 2 is copied from there.
- [`exemplars/population-slice2/DECISION.md`](https://github.com/bbuchsbaum/crossform/tree/main/exemplars/population-slice2)
  — the hard half, and the record of why OpenNeuro `ds003745` was chosen
  for it: searchlight-level transport in a common normalized space, with
  real sink mass from partial coverage and a functionally-informed
  transport cross-fitted from independent data. Slice 1 is the case
  where transport is trivial; slice 2 is the case this contract was
  written for.
- [`vignette("conservative-frames")`](https://bbuchsbaum.github.io/crossform/articles/conservative-frames.md)
  — the attribution instrument every participant’s ledger comes from,
  and why summing a detection map is the error this layer is built to
  avoid.
- [`vignette("interpreting-results")`](https://bbuchsbaum.github.io/crossform/articles/interpreting-results.md)
  — the interpretive traps that apply to a single participant’s map and
  apply again, unchanged, to a group one.
