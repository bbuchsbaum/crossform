# Read a geometry in model coordinates: the trace split and a structure fit

A model geometry reads one crossvalidated geometry through a
[`model_basis()`](https://bbuchsbaum.github.io/crossform/reference/model_basis.md).
It lowers the plan's relation into model coordinates (`B~ = Q'B`,
composing the basis with each partition's extractor), materializes the
complete compressed form `S_x = Q' G_x Q` at every measurement, and
reads two fixed trace queries from the original plan. From those it
reports two things that live on different layers and are kept apart.

## Usage

``` r
model_geometry(
  x,
  basis,
  structure = c("isotropic", "diagonal", "block", "shared"),
  rank = NULL,
  component = c("total", "coherent", "configuration"),
  warning_fraction = 0.9,
  tolerance = 1e-10,
  max_sweeps = 1000L,
  training = NULL
)
```

## Arguments

- x:

  An `effect_geometry_plan` from
  [`plan_geometry()`](https://bbuchsbaum.github.io/crossform/reference/plan_geometry.md)
  on the **original** relation (condition effects, not model
  coordinates), a complete self form under the implicit identity metric
  or one fixed
  [`neural_metric()`](https://bbuchsbaum.github.io/crossform/reference/neural_metric.md).
  The view re-plans the lowered relation with the same frame, pairing,
  metric and compute policy and executes both; there is no separate
  compute override, because the plan already carries one.

- basis:

  A
  [`model_basis()`](https://bbuchsbaum.github.io/crossform/reference/model_basis.md)
  declared over the plan's effects in the plan's order (build it with
  `conditions = relation$effect_space`).

- structure:

  Which model-side metric to fit; see the table above.

- rank:

  The rank budget: `NULL` for none (block: each model's own dimension;
  shared: the basis dimension), one positive whole number (shared: the
  budget; block: shared by every model), or one per model for block,
  named by model. Budgets are clamped to the dimension they act on and
  the clamp is recorded. Not accepted for isotropic or diagonal.

- component:

  Geometry component to read and fit.

- warning_fraction:

  Span fraction at or above which a shared fit is admitted but labelled
  `near_saturated`. At a span fraction of one the shared fit is refused
  (see Refusals).

- tolerance:

  Positive relative tolerance for the nonnegative and block solvers'
  convergence.

- max_sweeps:

  Positive cap on block coordinate descent sweeps per measurement.

- training:

  `NULL`, or
  [`metric_training_policy()`](https://bbuchsbaum.github.io/crossform/reference/metric_training_policy.md)
  with kind `"exclude_evaluation"`, which additionally cross-fits the
  learned readout of the shared structure: for each evaluation edge of
  the plan's pairing the rank-`s` eigenbasis of the compressed form is
  learned from the pairing's edges disjoint from it (never from a
  product the pairing did not declare), the energy is read on the edge
  itself, and the edges are reduced with the pairing's weights. The
  result, `$cross_fit`, is signed and on the estimation layer: zero in
  expectation under pure noise when the pairing's partitions are
  independent, where the plug-in fitted energy is not. It needs a
  pairing over at least four partitions in which every edge has a
  disjoint edge, and refuses otherwise. Defined for the shared structure
  only in this version. A tie in the training spectrum straddling the
  rank cut would make the readout platform dependent; for continuous
  data that is a measure-zero event.

## Value

An `effect_model_geometry`.

## Details

On the **signed estimation layer**, the trace decomposition
\\\operatorname{tr}(H G_x H) = \operatorname{tr}(P G_x) +
\operatorname{tr}((H - P) G_x)\\: the centered total energy of the
geometry equals its model-addressable energy, which is `tr(S_x)`, plus
its model-orthogonal energy. Both terms are fixed bilinear queries with
vanishing marginals, so both are unbiased and signed, and the identity
is exact at every measurement. It mirrors
`coherent + configuration = total` along the representational axis
instead of the spatial one.

On the **latent layer**, one structure fit of `S_x` to the model family:
the nonnegative model-side metric `C` in the parameterization
`G_x ~ T C T'`, with `T = Q R` the stacked model factor. The four
structures are readings of one form and answer four different questions:

|  |  |  |
|----|----|----|
| structure | `C` | tests |
| `isotropic` | one nonnegative weight per model | the truncated model geometries with nonnegative weights |
| `diagonal` | one nonnegative weight per model coordinate | the model axes with unequal, possibly sparse weights |
| `block` | one PSD rank-limited block per model | an independent learned metric per model (with one model, the shared fit) |
| `shared` | one PSD rank-limited matrix over the joint span | a learned readout of the joint model span, not the supplied geometries |

Only the isotropic structure tests the model geometries as supplied. The
shared fit is
[`latent_geometry()`](https://bbuchsbaum.github.io/crossform/reference/latent_geometry.md)'s
rank-budgeted truncation of `S_x` (`[S_x]_{+,s}`, the closed-form
Frobenius minimizer over PSD forms of rank at most `s`), not a second
implementation of it. Every fit returns `C` and never a factor `W`,
because `W` is identified only up to rotation. The fitted energy `tr(A)`
and the residual `||S_x - A||_F^2` are latent quantities and are never
reported as fractions of a signed denominator.

## Structure

One row per spatial measurement, with the estimation-layer split and the
latent-layer fit kept in separately named fields.

- `$addressable`, `$orthogonal`, `$centered_total`: the signed trace
  decomposition per measurement; `$identity_error` is the largest
  departure of `addressable + orthogonal` from `centered_total`.

- `$fit`: a data frame per measurement with `residual` (Frobenius, in
  model coordinates), `fitted_energy` (`tr(A)`), `iterations`,
  `converged` (block descent needs `max_sweeps` of at least two to
  observe a quiet sweep unless the basis has one model), and for the
  shared structure the two-part moved mass `clipped_negative_mass` and
  `truncated_positive_mass`.

- `$coefficients`: the fitted `C`, a model-coordinate by
  model-coordinate by measurement array, coordinates named `<model>.<k>`
  as in the basis's `$R`. For `isotropic` and `diagonal`, `$weights`
  additionally holds the nonnegative weights per model or per
  coordinate.

- `$cross_fit` (with `training`): `energy` per measurement, the per-edge
  energies `by_edge`, the `edges` with their weights, the partitions
  each readout was trained on and how many pairing edges trained it, the
  `rank`, the policy, and the pairing's declared `independence`. It is
  listed under `$layer$estimation`; the plug-in `$fit$fitted_energy`
  stays on the latent layer, and the two are never the same number. The
  execution `$receipt` carries none of this; its schema is sealed.

- `$structure`, `$hypothesis`, `$tests_supplied_geometries`: which fit
  was made and what it tests.

- `$rank`: the requested and effective budgets and whether any was
  clamped.

- `$df`: `fit`, the free parameters of the structure, beside `free`,
  those of the unconstrained rank-limited PSD form on `n - 1` centered
  directions; `$model_sizes` are the per-model dimensions they were
  computed from.

- `$span_fraction`, `$near_saturated`: the basis's span fraction and the
  warning label.

- `$basis_signature`, `$source_scientific_plan_id`,
  `$lowered_scientific_plan_id`: the identities the reading was built
  from; `$receipt` is derived from the lowered execution.

- `$reading`: the latent reading line, which applies to `$fit`,
  `$coefficients` and `$weights` and not to the trace split.

Any element not listed here is internal and may change.

## Refusals

Each is an `effect_capability_refusal` (see
[`catch_refusal()`](https://bbuchsbaum.github.io/crossform/reference/catch_refusal.md))
in namespace `"model_coordinate"`. A plan on a whitened or learned
metric schedule refuses capability `"model_coordinate_lowering"` with
reason `metric_schedule_not_lowerable`; a shared fit on a basis whose
span fraction is one refuses capability `"model_geometry_test"` with
reason `model_span_saturated`, because a full-span basis reproduces the
best rank-limited PSD approximation of the centered geometry whatever
models produced it; a block fit over a single model is the shared fit
under another name and is gated the same way; a cross-fit on fewer than
four partitions refuses capability `"cross_fitted_model_energy"` with
reason `insufficient_disjoint_edges`.

## See also

[`model_basis()`](https://bbuchsbaum.github.io/crossform/reference/model_basis.md)
for the basis;
[`latent_geometry()`](https://bbuchsbaum.github.io/crossform/reference/latent_geometry.md)
for the rank-budgeted projection the shared fit delegates to;
[`geometry_spectrum()`](https://bbuchsbaum.github.io/crossform/reference/geometry_spectrum.md)
for the signed spectrum of a lowered geometry;
[`rsa()`](https://bbuchsbaum.github.io/crossform/reference/rsa.md) for
the fixed linear reading of model RDMs on the original effects.

Other geometry plans and views:
[`aggregate_first()`](https://bbuchsbaum.github.io/crossform/reference/aggregate_first.md),
[`bilinear_query()`](https://bbuchsbaum.github.io/crossform/reference/bilinear_query.md),
[`coherence_spectrum()`](https://bbuchsbaum.github.io/crossform/reference/coherence_spectrum.md),
[`compute_policy()`](https://bbuchsbaum.github.io/crossform/reference/compute_policy.md),
[`contrast_energy()`](https://bbuchsbaum.github.io/crossform/reference/contrast_energy.md),
[`contribution()`](https://bbuchsbaum.github.io/crossform/reference/contribution.md),
[`crossnobis()`](https://bbuchsbaum.github.io/crossform/reference/crossnobis.md),
[`evaluate_geometry()`](https://bbuchsbaum.github.io/crossform/reference/evaluate_geometry.md),
[`example_fmri_effects()`](https://bbuchsbaum.github.io/crossform/reference/example_fmri_effects.md),
[`geometry_component()`](https://bbuchsbaum.github.io/crossform/reference/geometry_component.md),
[`geometry_spectrum()`](https://bbuchsbaum.github.io/crossform/reference/geometry_spectrum.md),
[`latent_geometry()`](https://bbuchsbaum.github.io/crossform/reference/latent_geometry.md),
[`materialize_geometry()`](https://bbuchsbaum.github.io/crossform/reference/materialize_geometry.md),
[`plan_crossnobis()`](https://bbuchsbaum.github.io/crossform/reference/plan_crossnobis.md),
[`plan_geometry()`](https://bbuchsbaum.github.io/crossform/reference/plan_geometry.md),
[`plot_views`](https://bbuchsbaum.github.io/crossform/reference/plot_views.md),
[`query_geometry()`](https://bbuchsbaum.github.io/crossform/reference/query_geometry.md),
[`rdm()`](https://bbuchsbaum.github.io/crossform/reference/rdm.md),
[`rsa()`](https://bbuchsbaum.github.io/crossform/reference/rsa.md),
[`variation_query()`](https://bbuchsbaum.github.io/crossform/reference/variation_query.md)

## Examples

``` r
example <- example_fmri_effects()
relation <- example$fit$relation
plan <- plan_geometry(relation, example$frame,
  cross_partitions(relation, independence = "independent",
    generalizes_over = "run"))

# The README category model as a basis over the plan's four conditions:
# its Gram has rank one, so the basis spans one of three centered
# directions.
basis <- model_basis(list(category = example$model_rdm),
  distance = "squared_euclidean", conditions = relation$effect_space)
basis$span_fraction
#> [1] 0.3333333

# The trace split at every searchlight, and the isotropic fit: one
# nonnegative weight on the category geometry.
reading <- model_geometry(plan, basis)
reading
#> <effect_model_geometry>
#>   measurements: 280
#>   component:    total
#>   structure:    isotropic,
#>                 tests the truncated model geometries with nonnegative weights
#>   basis:        1 model (category); span fraction 0.333
#>   addressable:  signed tr(S_x), median 0.008763 (range -0.0109 to 4.103),
#>                 orthogonal median -0.003373 (range -0.03502 to 0.02833),
#>                 split closes to 8.9e-16
#>   fit:          fitted energy median 0.008763 (range 0 to 4.103),
#>                 residual median 1.233e-32 (range 0 to 0.0001188),
#>                 df 1 of 3 free
#>   reading:      latent descriptive layer; not for inference,
#>                 applies to fit, coefficients and weights
#>  measurement addressable orthogonal centered_total fit_residual fitted_energy
#>            1   0.0003580 -0.0002863      7.171e-05    1.175e-38      0.000358
#>            2   0.0129223  0.0085732      2.150e-02    3.009e-36      0.012922
#>            3  -0.0026194 -0.0033962     -6.016e-03    6.861e-06      0.000000
#>            4  -0.0006461  0.0043554      3.709e-03    4.174e-07      0.000000
#>            5   0.0130512  0.0024283      1.548e-02    3.009e-36      0.013051
#>            6   0.0039387 -0.0067667     -2.828e-03    0.000e+00      0.003939
#>  converged weight_category
#>          1        0.000358
#>          1        0.012922
#>          1        0.000000
#>          1        0.000000
#>          1        0.013051
#>          1        0.003939
#>   ... 274 more measurements
#>   next:         as.data.frame(x), x$coefficients, x$fit
head(as.data.frame(reading))
#>   measurement   addressable    orthogonal centered_total fit_residual
#> 1           1  0.0003579634 -0.0002862531   7.171038e-05 1.175494e-38
#> 2           2  0.0129223176  0.0085732138   2.149553e-02 3.009266e-36
#> 3           3 -0.0026193537 -0.0033961656  -6.015519e-03 6.861014e-06
#> 4           4 -0.0006460581  0.0043554350   3.709377e-03 4.173911e-07
#> 5           5  0.0130511532  0.0024283483   1.547950e-02 3.009266e-36
#> 6           6  0.0039387093 -0.0067666814  -2.827972e-03 0.000000e+00
#>   fitted_energy converged weight_category
#> 1  0.0003579634         1    0.0003579634
#> 2  0.0129223176         1    0.0129223176
#> 3  0.0000000000         1    0.0000000000
#> 4  0.0000000000         1    0.0000000000
#> 5  0.0130511532         1    0.0130511532
#> 6  0.0039387093         1    0.0039387093
reading$identity_error
#> [1] 8.881784e-16

# The shared fit learns a readout of the model span; here the span is one
# dimensional so it coincides with the isotropic fit's energy.
shared <- model_geometry(plan, basis, structure = "shared")
all.equal(shared$fit$fitted_energy, reading$fit$fitted_energy)
#> [1] TRUE
```
