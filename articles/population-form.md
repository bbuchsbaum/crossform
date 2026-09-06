# Population form: one participant's ledger, carried to a group

How can participants with different native measurement grids contribute
to one group contrast estimate? Supply a conservative geometry plan for
each participant, a transport to shared group nodes, and a group model.
The result is a group estimate in budget units, with unmapped evidence
retained in a separate sink node.

Read [Conservative
frames](https://bbuchsbaum.github.io/crossform/articles/conservative-frames.md)
first: this guide builds on maps that add to a fixed total. The first
workflow fits two contrasts for six generated participants. Later
sections check conservation, identify a planted participant difference,
and separate between-participant uncertainty from within-participant
measurement error.

**The population API is experimental.** Its uncertainty procedures
condition on the supplied transport; they do not account for learning
that transport or provide simultaneous inference across nodes. Section 5
explains how to read the returned uncertainty and calibration labels.
Derivations and simulation regimes are in the [population
contract](https://github.com/bbuchsbaum/crossform/blob/main/design/population-form-contract.md)
(`population-form-v1`); § references below point there.

## 1. The three objects

### 1.1 Six participants, six conservative geometry plans

The generated participants have different native frame sizes, so their
maps must reach shared nodes before a group model can be fitted. Each
has four runs. This also lets section 4 form two disjoint pairs of runs
for estimating participant heterogeneity.

Every participant carries the same planted **consensus** direction in
effect space (face above house, present in every run, so it survives the
cross-partition product). One participant, `s06`, additionally carries a
direction nobody else has. Section 4 is about finding it without being
told.

Generated participant helper (expand to copy the complete setup)

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
```

Every participant reaches each of the three ordinary group nodes. Some
native territory remains unmapped and will enter the sink. The `tilts`
vector gives `s06` its additional house-versus-tool direction.

``` r

sizes <- c(s01 = 12L, s02 = 14L, s03 = 16L, s04 = 13L, s05 = 15L, s06 = 17L)
gains <- c(s01 = 1, s02 = 1.3, s03 = 0.8, s04 = 1.1, s05 = 0.9, s06 = 1.2)
tilts <- c(s01 = 0, s02 = 0, s03 = 0, s04 = 0, s05 = 0, s06 = 4)

subjects <- stats::setNames(lapply(names(sizes), function(id)
  pop_subject(id, sizes[[id]], gains[[id]], tilts[[id]])), names(sizes))
knitr::kable(data.frame(participant = names(sizes), native_nodes = sizes,
  noise_gain = gains, additional_direction = tilts), row.names = FALSE)
```

| participant | native_nodes | noise_gain | additional_direction |
|:------------|-------------:|-----------:|---------------------:|
| s01         |           12 |        1.0 |                    0 |
| s02         |           14 |        1.3 |                    0 |
| s03         |           16 |        0.8 |                    0 |
| s04         |           13 |        1.1 |                    0 |
| s05         |           15 |        0.9 |                    0 |
| s06         |           17 |        1.2 |                    4 |

These are ordinary single-participant objects:
[`plan_geometry()`](https://bbuchsbaum.github.io/crossform/reference/plan_geometry.md)
sealed six estimands and read no data.
[`voxelwise()`](https://bbuchsbaum.github.io/crossform/reference/voxelwise.md)
is conservative by default, which is the requirement checked by
[`plan_population()`](https://bbuchsbaum.github.io/crossform/reference/plan_population.md)
below.

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
#>   nodes:      12 native -> 3 group + sink
#>   semantics:  budget
#>   sink:       mass 1 of 12 rows, 8.3% of territory
#>   provenance: anatomical, fixed (cross-fit: none)
#>   built:      nearest group centre within radius 2, ties to the lowest gr...
#>   inference:  conditional_on_realized_transport; uncertainty not propagated
#>   signature:  sha256:0264631bb78a...
```

### 1.3 Fit on shared group nodes

``` r

plan <- plan_population(subjects, transports)
plan
#> <effect_population_plan>
#>   subjects:      6 (s01, s02, s03, s04 (+2 more))
#>   group nodes:   3 + sink
#>   sink:          present in 6 of 6 subjects, worst 23.5% of territory
#>   transport:     budget, anatomical
#>   model:         ~1 -> 1 column, rank 1
#>   normalization: none
#>   coverage:      all_planned, operator mass > 0 relative tolerance; node ...
#>   inference:     fixed; conditional on realized transport; uncertainty no...
#>   fit:           OLS (subject-constant weights), transport then fit
#>   estimand:      population-sha256:618aafd36c83...
#>   signature:     sha256:c5362ebf396c...
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
#> 1     s01           12           conservative      TRUE     0.08333333
#> 2     s02           14           conservative      TRUE     0.07142857
#> 3     s03           16           conservative      TRUE     0.18750000
#> 4     s04           13           conservative      TRUE     0.07692308
#> 5     s05           15           conservative      TRUE     0.13333333
#> 6     s06           17           conservative      TRUE     0.23529412
```

Participants have different native frames, so transport first creates
the shared node axis needed by the group model. The plan records this
order. Section 2 checks that a fixed contrast query can be evaluated
before or after the transport and group fit.

### Estimate two group contrasts

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

The coefficient array has axes **group node × query × model term**. Here
the only term is the intercept, so each coefficient is the mean
transported contrast across the six participants. The sink is retained
as a fourth row.

``` r

knitr::kable(fit$coefficients[, , "(Intercept)"], digits = 3,
  caption = "Group mean contrast evidence in budget units, including the sink.")
```

|        | face-house | house-tool |
|:-------|-----------:|-----------:|
| group1 |      7.381 |     14.141 |
| group2 |     13.934 |     29.884 |
| group3 |     11.157 |     28.307 |
|        |      8.628 |     21.983 |

Group mean contrast evidence in budget units, including the sink.
{.table}

### Add a group covariate or check an unsupported input

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
```

The subsequent analyses use the intercept-only `plan`, whose term is the
group mean.

### Transport options: unmapped territory, density, and learned alignment

Three fields explain how a transport will treat the native map.

**The sink is an accounting column.** The operator is `n × (m + 1)`, not
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
```

``` r


honest <- external_transport(P, semantics = "budget",
  provenance = list(method = "functional", details = "hyperalignment",
    fitting_sample = c("session-A", "session-B"),
    cross_fit = c("run1", "run2"),
    cross_fit_folds = c("fold-A", "fold-B")))
honest$provenance$conditioning
#> $source
#> [1] "hyperalignment"
#> 
#> $operator_status
#> [1] "estimated"
#> 
#> $fitting_sample
#> [1] "session-A" "session-B"
#> 
#> $cross_fit_folds
#> [1] "fold-A" "fold-B"
#> 
#> $circularity_control
#> [1] "cross_fit_partitions_declared"
#> 
#> $inference_scope
#> [1] "conditional_on_realized_transport"
#> 
#> $uncertainty_propagated
#> [1] FALSE
#> 
#> $marginal_over_transport
#> [1] FALSE
#> 
#> $excluded_uncertainty
#> [1] "transport_operator_estimation" "cross_fit_fold_assignment"    
#> 
#> $future
#> $future$capability
#> [1] "transport_uncertainty_propagation"
#> 
#> $future$status
#> [1] "not_implemented"
#> 
#> $future$requires
#> [1] "transport_sampling_law"              "joint_transport_response_resampling"
#> [3] "validated_propagation_operator"
```

Naming the partitions is what the refusal asks for, and the fitting
sample, folds, and fixed-versus-estimated status travel with the
operator into the estimand’s identity. This limits circularity; it does
**not** propagate uncertainty from learning the alignment. The sealed
record therefore says `conditional_on_realized_transport`,
`uncertainty_propagated = FALSE`, and `marginal_over_transport = FALSE`.
An anatomical or external transport never saw the responses, so it owes
no cross-fit record — which is why the six carriers above passed without
one — but it still declares the same conditional boundary for its fixed
operator.

## 2. Why querying before or after the group fit agrees

The query combines experimental coordinates; transport combines spatial
nodes; the group model combines participants. Under OLS with
subject-constant weights, these linear operations commute wherever the
relevant axes are shared. Native frames differ here, so transport must
precede the group fit. The query can still be evaluated on either side
of those operations.

### Check the transported values

First compare the bank executor’s participant values with individually
queried and transported contrasts. This checks row alignment and
scaling. Both routes use the same query operation, so agreement alone
does not test a different order of evaluation.

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
```

### Reverse the query order

[`materialize_population()`](https://bbuchsbaum.github.io/crossform/reference/materialize_population.md)
runs the *other* order. It carries all six packed coordinates of the
complete `3 × 3` form to the group nodes, fits the group model at every
one of them, and hands back a coefficient **form** per node. The
[`rdm()`](https://bbuchsbaum.github.io/crossform/reference/rdm.md) view
then contracts the (face, house) edge out of those coefficients —
**afterwards**. That is `transport → fit → query`.

The face − house contrast is that RDM edge, so the two routes should
agree at each group node.

``` r

form <- materialize_population(plan)
edges <- as.data.frame(rdm(form))
edge <- edges[edges$left == "face" & edges$right == "house", ]
node_ids <- as.character(fit$index$node)

# One row per node: with a richer group model `rdm()` emits one row per node
# per term, and `match()` would silently take the first.
```

``` r


commutation <- max(abs(
  edge$estimate[match(node_ids, edge$node)] -
    fit$coefficients[node_ids, "face-house", "(Intercept)"]
))
c(query_first = fit$basis, transport_first = form$basis)
#>     query_first transport_first 
#>    "query_bank" "complete_form"
c(commutation = commutation)
#>  commutation 
#> 3.552714e-15
```

The routes agree to the `1e-12` tolerance in §11. Their order of
arithmetic differs: the first transports a queried value, while the
second transports the complete form and queries the fitted result.

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
#> 1.517733e-16 1.529792e-16
```

Removing the sink loses the evidence assigned to unmapped native
territory. Compare the conservation gap with and without that row:

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
#>   1.529792e-16   4.348050e-01   2.352941e-01
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
```

The diagnostic view keeps those warnings beside their subject
provenance. It summarizes node-wise coverage and sink exposure, reports
descriptive associations with the declared design and outcomes, and
shows what a declared coverage or retained-territory threshold would
remove. The thresholded summary is never substituted for the fitted
result:

``` r

diagnostics <- population_diagnostics(fit,
  minimum_coverage = 0.8,
  minimum_transport_quality = 0.7,
  material_change = 0.2)
c(cells = nrow(diagnostics$cells), warnings = nrow(diagnostics$warnings))
#>    cells warnings 
#>        8        0
unique(diagnostics$sensitivity$target_status)
#> [1] "sensitivity_descriptive_not_primary"
```

The component identity also survives the group model, but only when all
three fits share that exact plan and its cellwise subject sets. The
decomposition object checks comparability first and carries the
cross-component covariance needed for derived coefficient contrasts:

``` r

coherent_fit <- estimate_population(plan, bank, component = "coherent")
configuration_fit <- estimate_population(plan, bank,
  component = "configuration")
decomposition <- population_decomposition(
  fit, coherent_fit, configuration_fit, estimator = "HC3")
c(coefficient_gap = decomposition$max_coefficient_gap,
  covariance_gap = decomposition$direct_derived_total_covariance_gap)
#> coefficient_gap  covariance_gap 
#>    0.000000e+00    9.094947e-13
```

This is an additive estimand law. It does not make coherent and
configurational contributions separate neural mechanisms.

For presentation,
[`population_component_view()`](https://bbuchsbaum.github.io/crossform/reference/population_component_view.md)
maps all three coefficients and their intervals onto one signed
symmetric axis. Its `$data` and `$coverage` tables are the plotting
inputs, so the visual mapping can be audited without comparing raster
files:

``` r

component_view <- population_component_view(decomposition, "(Intercept)",
  query = "face-house")
component_view$axis
#> $type
#> [1] "shared_symmetric"
#> 
#> $limits
#> [1] -21.44838  21.44838
#> 
#> $zero
#> [1] 0
#> 
#> $sign
#> [1] "positive_up_negative_down"
#> 
#> $units
#> [1] "signed transported evidence coefficient"
head(component_view$data[, c("component", "node", "estimate", "lower",
  "upper", "visual_magnitude")])
#>   component   node  estimate    lower     upper visual_magnitude
#> 1     total group1  7.380736 5.676156  9.085316        0.3441162
#> 2     total group2 13.933647 6.418910 21.448384        0.6496362
#> 3     total group3 11.156569 8.160234 14.152904        0.5201590
#> 4     total <sink>  8.627982 2.760327 14.495636        0.4022672
#> 5  coherent group1  7.380736 5.676156  9.085316        0.3441162
#> 6  coherent group2 13.933647 6.418910 21.448384        0.6496362
```

Scale-profile bands select their uncertainty method explicitly. They
retain the same coefficient and `subject_set_id` as the point, annotate
`n`, coverage fraction, and effective sample size, and leave sparse
cells as gaps. These are pointwise bands: simultaneous coverage and maxT
are not implemented or calibrated.

The matched hierarchical certification behind the supported-regime
statement uses 200 paired 24-subject replications. It recovers the
planted component ordering and scale profile while keeping an
informative-coverage arm as a failed marginal target; see
`inst/extdata/certification/population-interpretability-verdicts.csv`.

``` r

scale_profile <- population_scale_profile(decomposition, "(Intercept)",
  query = "face-house", interval = "HC3")
interval_record <- scale_profile$interval
knitr::kable(data.frame(
  field = c("Method", "Coverage", "Simultaneous coverage", "Calibration scope"),
  value = c(interval_record$method, interval_record$semantics,
    interval_record$simultaneous_coverage, interval_record$calibration_scope)
))
```

| field | value |
|:---|:---|
| Method | HC3 |
| Coverage | pointwise |
| Simultaneous coverage | not_available_unimplemented_uncalibrated |
| Calibration scope | Matched simulation regimes in inst/extdata/certification/population-calibration-results.csv; no marginal claim under informative coverage or transport estimation. |

``` r

knitr::kable(head(scale_profile$data[, c("component", "node", "estimate",
  "lower", "upper", "n", "fraction", "subject_set_id", "gap")]), digits = 3)
```

| component | node   | estimate | lower |  upper |   n | fraction | subject_set_id | gap   |
|:----------|:-------|---------:|------:|-------:|----:|---------:|:---------------|:------|
| total     | group1 |    7.381 | 5.676 |  9.085 |   6 |        1 | set1           | FALSE |
| total     | group2 |   13.934 | 6.419 | 21.448 |   6 |        1 | set1           | FALSE |
| total     | group3 |   11.157 | 8.160 | 14.153 |   6 |        1 | set1           | FALSE |
| total     |        |    8.628 | 2.760 | 14.496 |   6 |        1 | set1           | FALSE |
| coherent  | group1 |    7.381 | 5.676 |  9.085 |   6 |        1 | set1           | FALSE |
| coherent  | group2 |   13.934 | 6.419 | 21.448 |   6 |        1 | set1           | FALSE |

An interpretive display can bind the effect rows to coverage, effective
N, sink territory, and retained-territory support by exact node/query
identity. Selections are synchronized and recorded, while each
diagnostic retains its own unit-labelled axis:

``` r

diagnostic_view <- population_diagnostic_view(scale_profile, diagnostics,
  query = "face-house")
names(diagnostic_view$panels)
#> [1] "coverage"          "effective_n"       "sink"             
#> [4] "transport_quality"
diagnostic_view$filters
#> $node
#> [1] "group1" "group2" "group3" "<sink>"
#> 
#> $query
#> [1] "face-house"
#> 
#> $effect_rows_before
#> [1] 12
#> 
#> $effect_rows_after
#> [1] 12
#> 
#> $support_rows_before
#> [1] 8
#> 
#> $support_rows_after
#> [1] 4
#> 
#> $operation
#> [1] "synchronized_exact_key_selection"
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
#>    mode1    mode2    mode3    mode4    mode5    mode6 
#> 2884.697    4.658    0.473    0.000   -8.933  -27.363
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
#>         2.0000         0.0124         1.0040
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
#> s01 -0.173  0.753
#> s02 -0.188 -0.106
#> s03 -0.182  0.091
#> s04 -0.191 -0.096
#> s05 -0.178 -0.636
#> s06  0.913 -0.005
leading <- het$loadings[, 1L]
others <- leading[names(leading) != "s06"]
names(which.max(abs(leading)))
#> [1] "s06"
c(
  separation = round(max(abs(leading)) / max(abs(others)), 2),
  mode1_share = round(het$latent$cumulative[[1L]], 3)
)
#>  separation mode1_share 
#>       4.770       0.998
```

Mode 1 is one participant against the rest — `s06` at 0.91 with
everybody else between −0.15 and −0.22, the signature of a single odd
participant rather than a gradient — and it holds over 90 % of the
projected heterogeneity mass. Pass `nodes =` to reconstruct the
*geometry* of a mode at named group nodes and see which effect pair it
lives in.

[`population_influence()`](https://bbuchsbaum.github.io/crossform/reference/population_influence.md)
links that cross-fitted Gram to leave-one-subject coefficient and
component changes, exact cellwise coverage, and each subject’s transport
provenance. Its default court is bounded; a larger run requires
`mode = "deep"`. The output is descriptive and never an automatic
exclusion rule:

``` r

influence <- population_influence(decomposition, heterogeneity = het)
head(influence$influence[, c("subject", "node", "query", "term",
  "component", "abs_delta", "primary_subject_set_id", "sink_territory")])
#>   subject   node      query        term component  abs_delta
#> 1     s01 group1 face-house (Intercept)     total 0.05827883
#> 2     s02 group1 face-house (Intercept)     total 0.11219769
#> 3     s03 group1 face-house (Intercept)     total 0.23224745
#> 4     s04 group1 face-house (Intercept)     total 0.08525023
#> 5     s05 group1 face-house (Intercept)     total 0.34584411
#> 6     s06 group1 face-house (Intercept)     total 0.49286527
#>   primary_subject_set_id sink_territory
#> 1                   set1     0.08333333
#> 2                   set1     0.07142857
#> 3                   set1     0.18750000
#> 4                   set1     0.07692308
#> 5                   set1     0.13333333
#> 6                   set1     0.23529412
```

### Why the estimator is not a free choice

The plug-in Gram books every participant’s within-subject sampling noise
as between-participant heterogeneity. The contract measures the
consequence: over 2 000 Monte Carlo replications at `N = 12`, the
plug-in inflates the heterogeneity trace by **+62.7 %** (mean `tr Q^H`
of `+17.61` against a true `+10.82`, predicted bias `+6.72` versus
measured `+6.79`), while the cross-fitted estimator’s bias is within one
Monte Carlo standard error of zero (§6.2). The bias mechanism is
general; the 62.7 % magnitude belongs to that simulation’s signal and
noise regime. In this guide’s single fixture, compare the trace
estimates and the leading loading directions:

``` r

plug <- heterogeneity(plan, estimator = "plug_in")
c(
  cross_fit_trace = round(sum(diag(het$gram)), 3),
  plug_in_trace = round(sum(diag(plug$gram)), 3),
  mode1_agreement = round(abs(sum(het$loadings[, 1L] * plug$loadings[, 1L])), 4)
)
#> cross_fit_trace   plug_in_trace mode1_agreement 
#>       2853.5320       3062.0320          0.9999
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
```

## 5. Uncertainty and descriptive prevalence

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
  residual_df = unique(as.numeric(between$residual_df)),
  estimate_gap = max(abs(between$estimate[, , "(Intercept)"] - hand_estimate)),
  se_gap = max(abs(between$se[, , "(Intercept)"] - hand_se))
)
#>  residual_df estimate_gap       se_gap 
#> 5.000000e+00 1.065814e-14 3.552714e-15
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

The expanded certification court now compares classical, HC3, and a
null-imposed wild bootstrap on the same 500 datasets in each of eight
regimes. Its versioned results are in
`inst/extdata/certification/population-calibration-results.csv`. HC3
improves coverage in the declared heteroskedastic arm, but no method is
licensed for a marginal claim under informative coverage, and every
result remains conditional on the realized transport. That bounded
simulation evidence is why the API still reports
`calibration = "uncalibrated"` rather than turning a selected synthetic
success into a universal guarantee.

``` r

between$calibration
#> [1] "uncalibrated"
c(level = between$level, t_max = round(max(abs(between$t), na.rm = TRUE), 3))
#>  level  t_max 
#>  0.950 12.193
```

Report the statistic; do not report a p-value derived from it without an
argument that the section 7.5 transport diagnostics are benign.

The classical interval assumes equal subject-level variance. HC3 is the
leverage-adjusted sandwich sensitivity analysis when that assumption is
not credible. It uses the same exact cellwise subject set and returns
the full coefficient covariance, along with the leverage and assumptions
that produced it:

``` r

hc3 <- population_uncertainty(fit, estimator = "HC3")$between
ratio <- hc3$se / between$se
c(
  estimator = hc3$estimator,
  max_leverage = round(max(hc3$max_leverage, na.rm = TRUE), 3),
  min_se_ratio = round(min(ratio, na.rm = TRUE), 3),
  max_se_ratio = round(max(ratio, na.rm = TRUE), 3)
)
#>    estimator max_leverage min_se_ratio max_se_ratio 
#>        "HC3"      "0.167"      "1.095"      "1.095"
```

Neither estimator repairs an unidentified cell. Rank deficiency,
saturation, or leverage with `1 - h` at the declared numerical tolerance
produces an explicit per-cell refusal in `$status` and `$reason` instead
of an infinite or silently unstable standard error.

The included-versus-excluded uncertainty boundary is concrete here. HC3
includes heteroskedasticity in the observed between-subject residuals,
conditional on the transported values. If those values came from
`honest`, it would exclude variation from estimating the hyperalignment
and from assigning its cross-fitting folds. Cross-fitting protects the
evaluation split; it does not turn the HC3 interval into an interval
marginal over transport learning.

For a resampling diagnostic, the null-imposed wild bootstrap keeps the
participant—not the node—as the randomization unit. The single stored
subject-by-replicate weight matrix is reused across every node and
query, so the joint spatial/readout dependence within a participant is
never broken:

``` r

wild <- population_wild_bootstrap(
  fit, "(Intercept)", null = 0, replicates = 199, seed = 20260821
)
wild_table <- as.data.frame(wild)
wild_table[, c("node", "query", "observed_t", "p_value",
  "monte_carlo_se", "successful_replicates", "status")]
#>     node      query observed_t p_value monte_carlo_se successful_replicates
#> 1 group1 face-house 11.1304765   0.060    0.016792856                   199
#> 2 group2 face-house  4.7663121   0.060    0.016792856                   199
#> 3 group3 face-house  9.5713181   0.060    0.016792856                   199
#> 4 <sink> face-house  3.7798635   0.060    0.016792856                   199
#> 5 group1 house-tool  0.9456873   0.150    0.025248762                   199
#> 6 group2 house-tool  1.0205499   0.085    0.019719914                   199
#> 7 group3 house-tool  0.9748672   0.005    0.004987484                   199
#> 8 <sink> house-tool  1.0285528   0.060    0.016792856                   199
#>      status
#> 1 estimated
#> 2 estimated
#> 3 estimated
#> 4 estimated
#> 5 estimated
#> 6 estimated
#> 7 estimated
#> 8 estimated
```

The tested null and coefficient contrast are fields, not prose. The
plus-one p-value is accompanied by its Monte Carlo standard error, and
the sink remains in the same table—with its own estimate or
refusal—rather than disappearing. `weights = "mammen"` selects the
mean-zero, variance-one two-point Mammen distribution for an
asymmetric-error sensitivity analysis; Rademacher is the default.

The bootstrap inherits exactly the same boundary. Its participant
weights resample group residuals while holding the realized operator
fixed; they do not refit the transport. The future extension point and
its required ingredients are recorded at `wild$conditioning$future`.

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
```

A column that mixes native rows needs their cross-covariances. Omitting
those terms can underestimate or overestimate variance, depending on
their signs and the transport weights. The result therefore reports
unavailable variance for those cells. In the six-participant example,
every group node collects several native nodes, so none has an admitted
within-participant variance.

### 5.3 The two are never pooled

``` r

uncertainty$layers
#> [1] "between_subject" "within_subject"
uncertainty$separation
#> [1] "between_subject and within_subject are reported separately and are never pooled"
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

# Make partial coverage deliberate for this diagnostic: group3 moves to x=16,
# and the explicit policy says that its coefficient targets the participants
# available there. The introductory OLS fit above remains all-planned.
coverage_carrier <- function(features) anatomical_transport(
  native_coords = cbind(seq_len(features) - 1),
  group_coords = cbind(c(0, 5, 16)), semantics = "budget", radius = 2,
  native_index = paste0("f", seq_len(features))
)
coverage_transports <- stats::setNames(lapply(names(sizes), function(id)
  coverage_carrier(sizes[[id]])), names(sizes))
coverage_plan <- plan_population(subjects, coverage_transports,
  coverage_policy = "available_at_node")
coverage_fit <- estimate_population(coverage_plan, bank)
before <- coverage_fit$uncertainty
prevalence <- population_prevalence(coverage_fit,
  coverage_floor = length(subjects))
round(prevalence$sign$fraction, 3)
#>         query
#> node     face-house house-tool
#>   group1          1      0.667
#>   group2          1      0.833
#>   group3          1      0.667
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
```

The returned `$reference` uses `0.5` as a sign-balance comparison. This
is a useful reference under a continuous, symmetric null distribution,
not a universal null expectation: zero mean alone does not imply equal
probabilities of positive and negative estimates. Shared partition
products can produce asymmetric estimate distributions. Neither the
observed positive fraction nor its distance from `0.5` is a prevalence
test.

Coverage is reported beside it, because a high prevalence at a node only
a few participants reached is a different object from a high prevalence
at a node they all reached:

``` r

prevalence$coverage$minimum
#> group1 group2 group3 
#>      6      6      3
prevalence$coverage$below_floor
#> [1] "group3"
```

For this diagnostic, `group3` sits at `x = 16`, out of radius of every
native node in the smaller frames. Its exact available-at-node subject
set is recorded beside the fraction — a fact about the realized
transport, available before any of its numbers are interpreted, and
explicitly a different target from the all-planned fit used earlier.

## 6. Boundaries of the group result

**Calibration depends on the regime.** Classical, HC3, and
wild-bootstrap outputs have the conditional interpretation described in
section 5. They do not establish a marginal population claim under
informative coverage.
[`population_prevalence()`](https://bbuchsbaum.github.io/crossform/reference/population_prevalence.md)
remains descriptive.

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

**No marginal inference over transport estimation.** Functional
cross-fitting limits circularity, but all intervals and resampling
results condition on the realized operator. The transport-conditioning
record names the fitting sample, folds, included and excluded
uncertainty, and the currently unimplemented
`transport_uncertainty_propagation` extension point.

``` r

c(total = fit$ledger,
  coherent = estimate_population(plan, bank, component = "coherent")$ledger)
#>                    total                 coherent 
#>      "transported_total" "native_coherent_ledger"
```

**No cross-node sampling covariance**, and therefore no error bar on a
conserved budget. Overlapping node estimates can have nonzero
covariance; summing their variances omits these terms. Section 5.2
describes the supported single-node transport case.

This absence does **not** block the ordinary unweighted population
estimate. That OLS is performed across participants separately at each
node and query; it does not combine sampling errors across nodes. The
gate begins only when an operation asks for transported within-subject
precision, a conserved-budget variance, or joint spatial inference.

**No precision weighting.** `normalization = "precision_weighted"` is in
the closed set and refused at plan construction, because the per-subject
budget variance it would need is exactly the missing object above.

``` r

weighted <- catch_refusal(plan_population(subjects, transports,
  normalization = "precision_weighted"))
weighted$capability
#> [1] "precision_weighted_normalization"
```

The future admission law is explicit but not implemented: a covariance
must bind the subject, native-node and query indices; declare its error
model; pass finite, symmetry and positive-semidefinite checks; and fit
either a dense-byte or sparse-nonzero compute budget without implicit
densification. The tiny fixture in
`design/oracles/cross-node-covariance.R` proves
`Cov(P' z) = P' Cov(z) P` and shows why `sum(diag(Sigma))` is not a
budget variance. It exports no package function. Optional maxT,
simultaneous-band and multiple-comparison procedures are later consumers
requiring their own calibration; they are not part of the core
estimation gate.

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
