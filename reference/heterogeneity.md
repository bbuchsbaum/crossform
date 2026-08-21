# Decompose a population form into consensus and heterogeneity

`heterogeneity()` builds the subject Gram of
`design/population-form-contract.md` (`population-form-v1`) sections 5
and 6: the \\N \times N\\ matrix of inner products between participants'
deviations from the group fit, its spectrum, the subject loadings on its
modes, and — at group nodes the caller names — the geometry of those
modes.

## Usage

``` r
heterogeneity(x, ...)

# S3 method for class 'effect_population_plan'
heterogeneity(
  x,
  estimator = c("cross_fit", "plug_in"),
  component = c("total", "coherent", "configuration"),
  nodes = NULL,
  modes = NULL,
  partitions = NULL,
  projection = "psd_projection",
  coordinate_tile = NULL,
  ...
)

# S3 method for class 'effect_population_result'
heterogeneity(
  x,
  estimator = c("plug_in", "cross_fit"),
  nodes = NULL,
  modes = NULL,
  projection = "psd_projection",
  ...
)
```

## Arguments

- x:

  An `effect_population_plan` from
  [`plan_population()`](https://bbuchsbaum.github.io/crossform/reference/plan_population.md),
  read through the streamed complete-form route, or an
  `effect_population_result` from
  [`estimate_population()`](https://bbuchsbaum.github.io/crossform/reference/estimate_population.md),
  read in the space of its own query bank.

- ...:

  Passed to the method; every method refuses an unmatched argument
  rather than absorbing it.

- estimator:

  `"cross_fit"` (the default on a plan) or `"plug_in"`. See the section
  above; the choice is part of the record's identity.

- component:

  Which component of each participant's conservative geometry to carry:
  `"total"` (the default), `"coherent"` or `"configuration"`, as in
  [`materialize_population()`](https://bbuchsbaum.github.io/crossform/reference/materialize_population.md).

- nodes:

  Group nodes whose mode geometry is reconstructed, by identifier or
  position. `NULL`, the default, reconstructs none: the Gram, the
  spectrum and the loadings come out of the first pass, and the mode
  forms cost a second one.

- modes:

  Number of leading modes to reconstruct at those nodes. `NULL` takes
  `min(3, residual df)`.

- partitions:

  Optional list of two character vectors naming the partition halves the
  cross-fit reads across. `NULL` interleaves each participant's own
  partitions, which balances the halves against anything drifting across
  a session.

- projection:

  The named PSD projection the latent layer is built with, from
  [`latent_geometry()`](https://bbuchsbaum.github.io/crossform/reference/latent_geometry.md)'s
  closed set. It is validated on either method and enters the record's
  identity, but it does no work on a result: the query-space route is
  plug-in, so that record has no latent layer to project.

- coordinate_tile:

  Optional positive number of packed coordinates held in flight at once,
  as in
  [`materialize_population()`](https://bbuchsbaum.github.io/crossform/reference/materialize_population.md).
  The tile bounds memory and is not part of the estimand.

## Value

An `effect_population_heterogeneity`: a sealed record carrying the
`$gram`, its signed `$spectrum` and subject `$loadings`, the
`$mode_forms` reconstructed at the selected `$nodes`, the `$latent`
projection, and a `$receipt` recording the estimator, the partition
split, the streaming bound and every participant's read.
[`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html) returns
the subject loadings in long form.

`$mode_forms` is `NA` at a `(node, coordinate)` cell the Gram could not
see — a group node reached by no native mass under density semantics —
because zero there is a form and would read as agreement rather than as
absence. `$receipt$unresolved_columns` counts them.

## Details

Geometry-space covariance across `N` participants has rank at most
\\N-1\\, so the whole spectrum lives in an \\N \times N\\ matrix and the
\\D \times D\\ covariance never has to be formed. The Gram is additive
over group nodes, so it accumulates while streaming and the \\N \times m
\times p\\ stack is never built either.

## What the deviations are

At group node `u` the record holds `K_u = D_u D_u^T`, where `D_u` is the
\\N \times p\\ matrix of participant deviations *from the fitted value
under the group model*, in the packed `svec` coordinates of
`effect-form-v1` section 3. Under the default `~ 1` design that is each
participant minus the group mean; under any other design it is what is
left after the covariates. The global Gram is `sum_u K_u`, divided by
the residual degrees of freedom.

The `sqrt(2)` on the packed off-diagonals is load-bearing: only with it
is the Euclidean inner product of two packed rows the Frobenius inner
product of the two forms. A Gram built on a naive `upper.tri` packing is
off by an \\O(1)\\ amount and is not a Gram of geometries at all.

The sink is excluded. It is unmapped territory reported in budget units
under either semantics, so an inner product including it would add a
mass to a density; `$receipt$gram$sink_excluded` records the exclusion
and the sink budget stays readable on `$receipt$sink_budget`.

## The two estimators

The plug-in Gram books every participant's own within-subject sampling
noise as between-participant heterogeneity. Measured over 2 000 Monte
Carlo replications of the contract's fixture, it inflates the
heterogeneity trace by **+62.7 %** (section 6.2). The cross-fitted Gram
takes every inner product across two independent halves of each
participant's partitions, \\\widehat\Gamma\_{ii'} = \tfrac12(\langle
Z_i^{A},Z\_{i'}^{B}\rangle + \langle Z_i^{B},Z\_{i'}^{A}\rangle)\\,
which is unbiased for the Gram of the noiseless participant vectors — on
the diagonal too, because the two halves are independent.

The cross-fitted Gram is **indefinite by construction** (measured in 100
% of the contract's replications), and a single draw of its trace is not
an estimate of the between-participant trace: measured standard
deviation `3.58` against a mean of `10.81`. Both facts are reported and
never repaired.

`estimator = "cross_fit"` re-plans each participant over two disjoint
halves of their partitions, so it needs at least four partitions per
participant and costs a second streaming read. `estimator = "plug_in"`
reads one pass and is admissible for loading *directions*, which section
6.3 measures surviving the inflation; its `$latent` layer is `NULL`, and
the refusal that says why travels on `$receipt$latent_refusal` as a
record rather than an error.

**One departure from section 6.4, stated rather than buried.** Its
literal text refuses "any spectrum, variance-explained figure or
`n_eff`" from the plug-in Gram; `crossform` reports the plug-in
**spectrum** and refuses the other two. The Gram itself is a field of
the record and the validator re-derives its eigenvalues from it, so
withholding them would be theatre rather than discipline — while a
fraction of an inflated whole, or an effective count over noise modes,
is a number that misleads whatever it is labelled. The spectrum is
reported beside the estimator's name, the measured inflation and the
refusal, never alone.

Under `"cross_fit"` the mode *geometry* is reconstructed from the
average of the two halves' deviations — a plug-in reconstruction of a
direction, which section 6.3 admits provided the source is named. It is
named on `$receipt$gram$mode_reconstruction` and on the printed `modes`
line. The eigenvalues beside it still come from the cross-fitted Gram.

## The latent layer

The heterogeneity spectrum is signed. Every functional that treats it as
a nonnegative partition — a fraction, a cumulative curve, an effective
mode count — is defined only on a named projection of it, exactly as
[`latent_geometry()`](https://bbuchsbaum.github.io/crossform/reference/latent_geometry.md)
requires one level down (section 6.5). `$latent` carries the projected
spectrum, its cumulative curve, its `n_eff`, and the mass the projection
moved; `$spectrum` stays signed and carries none of them.

## Refusals

Each is an `effect_capability_refusal` (see
[`catch_refusal()`](https://bbuchsbaum.github.io/crossform/reference/catch_refusal.md)).

- `cross_fitted_subject_gram` — fewer than four partitions, a declared
  split that shares a partition between halves or names one a
  participant does not carry, or a metric schedule that cannot be
  re-planned onto a half.

- `identified_heterogeneity` — a saturated group design, where no
  residual degree of freedom is left for participants to differ in.

- `query_space_cross_fit` — `estimator = "cross_fit"` on an estimated
  result, which carries one partition split and cannot supply a second.

- `complete_form_normalization` — a plan declaring `"unit_budget"`, as
  in
  [`materialize_population()`](https://bbuchsbaum.github.io/crossform/reference/materialize_population.md).

- `selected_group_nodes` — a `nodes` identifier the population does not
  carry.

One further refusal is *recorded* rather than raised, because the record
it belongs to is still worth having: `plug_in_spectrum_functionals` sits
on `$receipt$latent_refusal` of every plug-in record and carries the
reasons and remedies for the absent `$latent` layer.

## References

`design/population-form-contract.md` (`population-form-v1`), sections 5
and 6.

## See also

[`estimate_population()`](https://bbuchsbaum.github.io/crossform/reference/estimate_population.md)
and
[`materialize_population()`](https://bbuchsbaum.github.io/crossform/reference/materialize_population.md)
for the consensus half of the same decomposition, and
[`latent_geometry()`](https://bbuchsbaum.github.io/crossform/reference/latent_geometry.md)
for the single-participant form of the projection discipline.

Other population transports:
[`anatomical_transport()`](https://bbuchsbaum.github.io/crossform/reference/anatomical_transport.md),
[`estimate_population()`](https://bbuchsbaum.github.io/crossform/reference/estimate_population.md),
[`external_transport()`](https://bbuchsbaum.github.io/crossform/reference/external_transport.md),
[`location_transport()`](https://bbuchsbaum.github.io/crossform/reference/location_transport.md),
[`materialize_population()`](https://bbuchsbaum.github.io/crossform/reference/materialize_population.md),
[`plan_population()`](https://bbuchsbaum.github.io/crossform/reference/plan_population.md),
[`population_prevalence()`](https://bbuchsbaum.github.io/crossform/reference/population_prevalence.md),
[`population_uncertainty()`](https://bbuchsbaum.github.io/crossform/reference/population_uncertainty.md),
[`population_views`](https://bbuchsbaum.github.io/crossform/reference/population_views.md),
[`transport_values()`](https://bbuchsbaum.github.io/crossform/reference/transport_values.md)

## Examples

``` r
# Four participants over four runs, so each cross-fit half holds two.
effects <- effect_space(c("face", "house"), basis_id = "heterogeneity:v1")
subject <- function(id, n, gain) {
  domain <- abstract_domain(n, coordinates = cbind(x = seq_len(n) - 1),
    feature_ids = paste0("f", seq_len(n)), id = id)
  run <- function(seed) {
    set.seed(seed)
    matrix(gain * stats::rnorm(2 * n), 2, n,
      dimnames = list(c("face", "house"), NULL))
  }
  rel <- relation(
    list(run1 = run(1), run2 = run(2), run3 = run(3), run4 = run(4)),
    effects = effects, domain = domain
  )
  plan_geometry(rel, compile_frame(voxelwise(), domain),
    cross_partitions(rel))
}
carrier <- function(n) anatomical_transport(
  native_coords = cbind(seq_len(n) - 1), group_coords = cbind(c(0, 3)),
  semantics = "budget", radius = 2
)
sizes <- c(s01 = 5L, s02 = 6L, s03 = 7L, s04 = 8L)
gains <- c(s01 = 1, s02 = 1.4, s03 = 0.7, s04 = 1.1)
subjects <- stats::setNames(lapply(names(sizes), function(id)
  subject(id, sizes[[id]], gains[[id]])), names(sizes))
plan <- plan_population(subjects, lapply(sizes, carrier))

# The cross-fitted subject Gram: one N-by-N matrix, streamed. Naming a
# group node reconstructs that node's mode geometry in a second pass.
split <- heterogeneity(plan, nodes = "group1", modes = 1L)
split
#> <effect_population_heterogeneity>
#>   estimator: cross_fit -- cross-partition inner products (section 6.4)
#>   gram:      4 x 4 over 2 group nodes, 3 packed coordinates; residual df 3
#>   ledger:    transported_total
#>   subjects:  s01, s02, s03, s04
#>   spectrum:  4.585e-01, 4.996e-16, -5.551e-17 (+1 more), 1 negative
#>   n_eff:     1 modes after psd_projection, moved mass 2.504 (84.52%)
#>   cross-fit: interleaved: [run1, run3] x [run2, run4]
#>   modes:     1 at group1 (direction: half average)
#>   estimand:  population-sha256:78c219ca45db...
#>   cross-fitted: indefinite by construction, and one draw of its trace is
#>     not an estimate of the between-subject trace (population-form-v1
#>     section 6.4). Reported as-is.
#>   next:         as.data.frame(x), x$loadings, x$gram

# Signed and routinely indefinite; the nonnegative functionals live on the
# named projection beside it.
split$spectrum
#>         mode1         mode2         mode3         mode4 
#>  4.585164e-01  4.996004e-16 -5.551115e-17 -2.503957e+00 
split$latent$n_eff
#> [1] 1
split$latent$moved_mass
#> [1] 2.503957

# Subject loadings on the leading mode, and that mode's geometry at the
# named group node, unpacked from its `svec` coordinates.
split$loadings[, "mode1"]
#>        s01        s02        s03        s04 
#>  0.6771134  0.2181094 -0.6635208 -0.2317020 
crossform:::.unsvec_symmetric(split$mode_forms["group1", "mode1", ], 2L)
#>           [,1]       [,2]
#> [1,] 0.2773710  0.1594747
#> [2,] 0.1594747 -0.4601125
split$receipt$gram$mode_reconstruction
#> [1] "plug_in_average_of_partition_halves"

# The plug-in Gram costs one pass instead of two and is admissible for
# directions; it carries no `n_eff`, and says why.
plug <- heterogeneity(plan, estimator = "plug_in")
is.null(plug$latent)
#> [1] TRUE
plug$receipt$latent_refusal$capability
#> [1] "plug_in_spectrum_functionals"
```
