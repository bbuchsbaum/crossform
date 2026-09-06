# Declare a model family as a basis of the effect axis

A model basis turns one or more model geometries over the relation's
conditions into a single declared value: an orthonormal, centered basis
`Q` of their joint span, and the factor `R` that recovers each model's
own coordinates (`T = Q R`, with `T = [T_1, ..., T_k]` the stacked model
factors). Lowering a relation through it, `B~_r = Q' B_r`, is an
ordinary relation on the model-coordinate effect space, and the ordinary
geometry plan then yields the complete geometry in model coordinates,
`S_x = Q' G_x Q`.
[`geometry_spectrum()`](https://bbuchsbaum.github.io/crossform/reference/geometry_spectrum.md)
and
[`latent_geometry()`](https://bbuchsbaum.github.io/crossform/reference/latent_geometry.md)
read that form as they read any other; two participants lowered through
one basis share one effect-space signature, which is what the population
layer pools on. The basis is built before any neural data is read and
its signature enters every plan downstream.

## Usage

``` r
model_basis(
  models = NULL,
  distance = NULL,
  features = NULL,
  rank = NULL,
  conditions = NULL,
  tolerance = 1e-10,
  negative_share = 0,
  kernels = NULL,
  normalize = c("none", "trace")
)
```

## Arguments

- models:

  One squared Euclidean model RDM over the conditions, or a named list
  of them. Each is condition-by-condition, symmetric, zero on the
  diagonal, and either labelled on both axes with the condition names
  (any order) or unlabelled with the condition count as its dimension.

- distance:

  The dissimilarity `models` declare. Required whenever `models` is
  supplied, and the only admitted value is `"squared_euclidean"`: only a
  squared Euclidean RDM has a centered Gram whose factor reproduces it.
  This is a narrower admission than
  [`rsa()`](https://bbuchsbaum.github.io/crossform/reference/rsa.md)'s,
  which reads any symmetric zero-diagonal matrix for a different
  purpose.

- features:

  One condition-by-feature model matrix, or a named list of them; rows
  are conditions, labelled or not under the same rule. A model given as
  features enters through `H F F' H`, which equals the Gram of its
  squared Euclidean RDM exactly.

- rank:

  `NULL` to keep every model's full effective rank, one positive whole
  number shared by all models, or one per model (named by model or in
  model order). Ranks are clamped to each model's effective rank and the
  clamp is recorded.

- conditions:

  The condition axis the basis is declared over: condition names, an
  [`effect_space()`](https://bbuchsbaum.github.io/crossform/reference/effect_space.md)
  (typically the relation's), or a
  [`condition_space()`](https://bbuchsbaum.github.io/crossform/reference/condition_space.md).
  `NULL` reads the order from the first labelled model. A space lends
  its shared unit and scale to the model coordinates and must carry one
  of each. Supplying a condition space additionally yields
  `$effect_map`, the ingestion-path product.

- tolerance:

  Positive relative tolerance under which an eigenvalue or singular
  value counts as zero, against the largest absolute one.

- negative_share:

  The largest admitted ratio of dropped negative Gram mass to positive
  Gram mass. Above it the model is refused as non-Euclidean. Roots
  within `tolerance` of zero count on neither side, so a Euclidean model
  records a dropped mass of exactly zero and `negative_share = 0` (the
  default) admits it. A larger value explicitly admits replacement by
  the positive part; input-kernel repair mass is recorded separately
  from centered-Gram mass. Numerical negative roots are also recorded
  separately and do not consume this allowance.

- kernels:

  One positive semidefinite condition-by-condition kernel, or a named
  list of them. Both axes follow the same alignment rules as RDMs. The
  input kernel is checked for positive semidefiniteness before
  centering, so centering cannot hide an indefinite constant direction.

- normalize:

  `"none"` preserves each retained model's scale. `"trace"` divides each
  retained kernel by its trace, after rank truncation. Use this to
  declare comparable scales before pooling model kernels. The original
  trace, retained trace and normalization scale are recorded.

## Value

An `effect_model_basis`: a sealed list with the condition-by-basis
matrix `$Q` (orthonormal columns, each summing to zero), the
basis-by-model-coordinate factor `$R` such that `Q %*% R` is the stacked
model factor, the ordered `$conditions`, `$n_conditions`, the basis
`$dimension` `q`, a per-model data frame `$models` (kind, requested and
effective rank, whether the rank was clamped, positive and dropped
negative Gram mass and their ratio), `$columns` naming each model's
columns of `R`, the full Gram `$spectra`, the declared `$distance` and
`$tolerance` and `$negative_share`, the `$span_fraction` `q / (n - 1)`
with `$saturated` true at one, `$span_overlap_rank` (the
model-coordinate count minus union dimension) and `$basis_condition`
(the condition number on the retained union span), `$baseline_invariant`
with the measured `$centering_error`, the model-coordinate
`$effect_space`, the `$extractor` that lowers a beta block (map `t(Q)`),
the `$effect_map` over the supplied condition space or `NULL`, and the
`$signature`. `$normalize` and `$centering` record the scale convention
and the condition-mean centering origin, units and scale. A supplied
effect space retains that space's signature in the centering
declaration, so identical labels cannot silently change their meaning. A
rank cut through a positive tied eigenspace refuses with reason
`ambiguous_rank_cut`; retaining the whole group is allowed.

## Construction

For each model the centered Gram is formed, `K_i = -H D_i H / 2` for a
squared Euclidean RDM, `K_i = H F_i F_i' H` for a feature matrix, or
`K_i = H K_input H` for an explicitly supplied kernel, and truncated at
a relative eigenvalue tolerance. The requested rank is clamped to the
effective rank and the clamp is recorded, never silent: the literal
factor with a rank above the model's true rank has a zero column, and a
thin QR of it returns a direction that is not centered and leaks a
condition-independent baseline into the lowered form. Every
decomposition runs in an explicit orthonormal basis of the centered
subspace, so a poorly conditioned model (a feature scaled far below the
others) cannot leak an uncentered direction either. The stacked factor
is then factored by a thin SVD at the same tolerance. Negative Gram
eigenvalues are dropped and their mass recorded; when it exceeds
`negative_share` of the positive mass the RDM is not a Euclidean
geometry and the constructor refuses. Overlapping model spans are
admitted for invariant form prediction; `$R` can be rectangular. Readers
that require unique redundant coefficients refuse there.
`$basis_condition` reports the condition number on the retained union
span. Centering, `Q'1 = 0`, is asserted and recorded, because it is
exactly the condition under which every reading of the lowered form
cancels additive baselines.

## Alignment

Models are aligned to the declared conditions by name, both axes of an
kernel or RDM and the rows of a feature matrix; an unlabelled model is
admitted only when its dimension equals the condition count. The
`$extractor` product lowers a beta block by **position**:
[`relation()`](https://bbuchsbaum.github.io/crossform/reference/relation.md)
aligns a labelled beta matrix to the basis's `$conditions` by row name
and refuses a mismatch, but an unlabelled block is taken in the basis's
condition order, so build the basis with
`conditions = relation$effect_space` when in doubt.

## What it is not

The basis is not a geometry, and building it subtracts nothing from the
neural side. The rank-limited fit of a geometry to the model span is a
latent projection and lives on the latent layer; the signed compressed
form is the estimation layer.

## Refusals

Each is an `effect_capability_refusal` (see
[`catch_refusal()`](https://bbuchsbaum.github.io/crossform/reference/catch_refusal.md))
in namespace `"model_coordinate"`. A model whose centered Gram carries
more negative mass than `negative_share` admits refuses capability
`"euclidean_model_geometry"` with reason `model_gram_indefinite`. A
positive rank-cut tie refuses `"identified_rank_projection"` with reason
`ambiguous_rank_cut`. Overlap itself is admitted: the fitted form does
not require a unique coefficient for every redundant model coordinate.

## Identity

The signature covers the ordered conditions, the model names, kinds and
declared distance, the tolerance, the requested and effective ranks, and
the unrounded `Q` and `R`. Signs of the computed factors are fixed
deterministically, but a model with tied Gram eigenvalues has a factor
determined only up to a rotation inside the tie, so two platforms whose
LAPACK resolves the tie differently can produce equal spans under
different signatures. Compare declarations by signature on one platform.
Across platforms compare each reconstructed model kernel using its
columns of `Q %*% R`; comparing only `Q %*% t(Q)` establishes equal
spans but discards the model's directional strengths.

## See also

[`effect_extractor()`](https://bbuchsbaum.github.io/crossform/reference/effect_extractor.md)
and
[`effect_map()`](https://bbuchsbaum.github.io/crossform/reference/effect_map.md),
the two forms the basis is lowered through;
[`relation()`](https://bbuchsbaum.github.io/crossform/reference/relation.md)
on the betas path and
[`plan_relation()`](https://bbuchsbaum.github.io/crossform/reference/plan_relation.md)
on the ingestion path;
[`latent_geometry()`](https://bbuchsbaum.github.io/crossform/reference/latent_geometry.md)
for the rank-limited reading of the lowered form;
[`rsa()`](https://bbuchsbaum.github.io/crossform/reference/rsa.md) for
the fixed linear reading of model RDMs that does not lower anything. Use
[`fit_geometry()`](https://bbuchsbaum.github.io/crossform/reference/fit_geometry.md)
and
[`score_geometry()`](https://bbuchsbaum.github.io/crossform/reference/score_geometry.md)
for regularized prediction and independent signed gain;
[`vignette("predictive-geometry")`](https://bbuchsbaum.github.io/crossform/articles/predictive-geometry.md)
gives the complete workflow.

## Examples

``` r
conditions <- c("face", "body", "house", "tool")

# A category model: animate against inanimate. Its Gram has rank one, so a
# requested rank of two is clamped to one and the clamp is on the record.
animate <- c(1, 1, 0, 0)
category <- outer(animate, animate, function(a, b) (a - b)^2)
dimnames(category) <- list(conditions, conditions)
basis <- model_basis(list(category = category),
  distance = "squared_euclidean", rank = 2)
basis
#> model_basis<q = 1 of 3 centered; 1 model; span 0.33>
#>     model kind rank_requested rank_effective rank_clamped dropped_negative_mass
#>  category  rdm              2              1         TRUE                     0
#> baseline_invariant: TRUE (max |Q'1| = 2.2e-16); basis_condition = 1
#> signature: model-basis-sha256:06f8e7882add...
basis$models[, c("rank_requested", "rank_effective", "rank_clamped")]
#>   rank_requested rank_effective rank_clamped
#> 1              2              1         TRUE
max(abs(colSums(basis$Q)))
#> [1] 2.220446e-16

# A second model given as features spans two more centered directions;
# with four conditions the family then saturates the effect axis.
shape <- cbind(elongation = c(0.2, 0.9, 0.1, 0.8), size = c(1, 2, 3, 1))
rownames(shape) <- conditions
both <- model_basis(list(category = category), "squared_euclidean",
  features = list(shape = shape))
both$dimension
#> [1] 3
both$span_fraction
#> [1] 1

# The lowering products: the extractor for betas, and the effect map when
# a condition space is supplied.
both$extractor$map
#>           face       body     house       tool
#> mc1 -0.4606768  0.0265771 0.8056283 -0.3715286
#> mc2 -0.3189599 -0.5816980 0.1725705  0.7280875
#> mc3  0.6603343 -0.6410312 0.2668003 -0.2861034
with_space <- model_basis(list(category = category), "squared_euclidean",
  conditions = condition_space(conditions))
with_space$effect_map
#> effect_map<1 effects x 4 conditions; canonical-amplitude>

# Duplicate models retain separate factors in one common span.
duplicate <- model_basis(
  list(a = category, b = category), distance = "squared_euclidean"
)
duplicate$span_overlap_rank
#> [1] 1
dim(duplicate$R)
#> [1] 1 2
```
