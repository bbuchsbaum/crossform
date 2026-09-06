# Conservative frames: attribution, not detection

How much of a contrast’s total evidence is assigned to each territory? A
**conservative frame** answers this question by partitioning one global
total among its spatial measurements. Its map can be added across
territories.

This guide assumes the `relation → frame → plan → contrast_energy()`
workflow from [the
introduction](https://bbuchsbaum.github.io/crossform/articles/introduction.md).
It uses the same generated fMRI example to compare two normalizations:

| Question | Normalization | Meaning of one node |
|----|----|----|
| Where is evidence locally strong? | `"local"`: normalize rows | Mean evidence density in the node’s support |
| How is the global total distributed? | `"conservative"`: normalize columns | Evidence assigned to the node from a fixed global budget |

A frame has one row per measurement and one column per neural feature.
Column normalization allocates each feature across the nodes that
contain it. Start with the conservation check and territory ledger; the
later sections explain scale, signed estimates, and metric choices. The
outputs show the numerical identities; quiet assertions check their
tolerances during rendering. Their derivations are in the [geometry
contract](https://github.com/bbuchsbaum/crossform/blob/main/design/conservative-geometry-contract.md).
References marked § below point to that contract.

``` r

example <- example_fmri_effects()
domain <- example$domain
pairing <- cross_partitions(
  example$fit$relation,
  independence = "independent",
  generalizes_over = "run"
)
energy <- function(frame) {
  contrast_energy(
    plan_geometry(example$fit$relation, frame, pairing), example$contrast
  )
}
c(features = domain$n_features, conditions = length(example$contrast))
#>   features conditions 
#>        280          4
```

The generated fixture plants two blocks of 19 voxels at opposite ends of
an 8 × 7 × 5 volume, at the same per-voxel amplitude, never overlapping.
They differ only in sign structure: in the **pattern** block the sign
alternates between neighbouring voxels, so a neighbourhood average
nearly cancels; in the **mean** block every voxel moves the same way.
Same effect size, two different effect *geometries*. Section 3 compares
their coherent shares across neighbourhood sizes. These are generated
data, not evidence about a brain; the [Haxby 2001
exemplar](https://github.com/bbuchsbaum/crossform/tree/main/exemplars/haxby2001)
(script `07`) is the real-data companion to this article.

## 1. Two instruments from one fit

The same relation, the same 8 mm neighbourhoods, the same six cross-run
pairs. Only the normalization changes.

``` r

detection <- compile_frame(searchlights(4, "local"), domain)
attribution <- compile_frame(searchlights(4, "conservative"), domain)
global <- compile_frame(whole_brain("none"), domain)
```

[`whole_brain()`](https://bbuchsbaum.github.io/crossform/reference/whole_brain.md)
defaults to `"local"`, which divides by the feature count. `"none"` is
the unnormalized whole-domain operator, and it is the *only* correct
comparator for a conservation claim (§2, precondition 1). Using the
default here is the usual source of a spurious conservation failure.

``` r

detection_mass <- frame_conservation(detection)
attribution_mass <- frame_conservation(attribution)
c(
  detection_conserves = detection_mass$conserved,
  attribution_conserves = attribution_mass$conserved
)
#>   detection_conserves attribution_conserves 
#>                 FALSE                  TRUE
c(
  detection_deviation = detection_mass$max_deviation,
  attribution_deviation = attribution_mass$max_deviation
)
#>   detection_deviation attribution_deviation 
#>          1.500000e-01          2.220446e-16
```

[`frame_conservation()`](https://bbuchsbaum.github.io/crossform/reference/frame_conservation.md)
reports on the **total** component only, and it says so in `$component`.
Section 4 below is why that scope is load-bearing rather than an
implementation detail.

Now the two maps, and the identity that separates them.

``` r

detection_map <- energy(detection)
attribution_map <- energy(attribution)
whole_domain <- energy(global)

identity_gap <- sum(attribution_map$total) - whole_domain$total
c(
  ledger_sum = sum(attribution_map$total),
  whole_domain_total = whole_domain$total,
  identity = identity_gap
)
#>         ledger_sum whole_domain_total           identity 
#>       1.552081e+02       1.552081e+02       5.684342e-14
```

That is claim 2 of the contract: for a column-normalized frame,
`sum_x <H, G_x> = <H, G_Omega>` for every fixed query `H`, exactly. The
280 searchlight values are not 280 measurements of a density; they are
280 pieces of one number.

The detection map has no such identity, and the way it fails is
instructive.

``` r

c(
  detection_sum = sum(detection_map$total),
  whole_domain_total = whole_domain$total,
  ratio = sum(detection_map$total) / whole_domain$total
)
#>      detection_sum whole_domain_total              ratio 
#>         157.987255         155.208072           1.017906
```

Here the detection sum happens to lie close to the global total. Row
normalization does not guarantee this: a feature’s aggregate weight
depends on all the overlapping neighbourhoods that include it. The sum
therefore has no global-budget interpretation.
[`contribution()`](https://bbuchsbaum.github.io/crossform/reference/contribution.md)
checks this requirement before aggregating a map:

``` r

territory <- ifelse(
  domain$feature_ids %in% example$truth$planted_feature_ids, "pattern block",
  ifelse(domain$feature_ids %in% example$truth$mean_feature_ids, "mean block",
    "elsewhere")
)
refused <- catch_refusal(contribution(detection_map, by = territory))
refused$capability
#> [1] "conservative_frame"
```

Read the two maps for what each is. The detection map’s peak (4.103) is
a density, comparable across nodes of different size and correctly
reported as a peak. The attribution map’s peak is a *share of a budget*,
and a large node holds more of it than a small one for the trivial
reason that it covers more features. Never compare them node-for-node.

## 2. The ledger the frame reallocates

Before using the attribution map, understand what it is made of —
because the answer bounds what a multiscale panel is allowed to claim.

Under an identity or diagonal metric the node total is
`G_x = sum_v w_xv b_v^L b_v^R'`. The per-voxel outer products form a
**fixed ledger**, and the frame appears only as a nonnegative linear map
applied after the multivariate content has been fixed per voxel. The
consequence is strong and easy to check: the attribution map at any
radius is a spatially smoothed version of one univariate per-voxel map,
and carries no cross-voxel information whatsoever (§3).

[`voxelwise()`](https://bbuchsbaum.github.io/crossform/reference/voxelwise.md)
is that per-voxel ledger — a point frame, trivially conservative.
Contracting it with the searchlight weights should reproduce the
searchlight map exactly.

``` r

ledger <- energy(compile_frame(voxelwise(), domain))
smoothed <- as.numeric(attribution$weights %*% ledger$total)
c(
  ledger_sums_to_global = sum(ledger$total) - whole_domain$total,
  smoothed_equals_map = max(abs(smoothed - attribution_map$total))
)
#> ledger_sums_to_global   smoothed_equals_map 
#>          8.526513e-14          0.000000e+00
```

The searchlight attribution map *is* the voxel ledger, reallocated.
Nothing between the two steps mixes features. This is scoped to metrics
for which `metric_capabilities()$feature_additive` is `TRUE`; section 6
is the case where the premise fails.

### Per-scale totals are your own weights

Now stack several radii.
[`searchlights()`](https://bbuchsbaum.github.io/crossform/reference/searchlights.md)
takes a vector, which builds an α-weighted **family**: each member is
column-normalized on its own, the weights sum to one, and the stack is
therefore conservative. Every row carries its own provenance.

``` r

radii <- c(3.01, 4.3, 6.1)
family <- compile_frame(searchlights(radii, "conservative"), domain)
head(family$index, 3)
#>      measurement      family node scale center     alpha
#> 1 radius-3.01::1 radius-3.01    1  3.01      1 0.3333333
#> 2 radius-3.01::2 radius-3.01    2  3.01      2 0.3333333
#> 3 radius-3.01::3 radius-3.01    3  3.01      3 0.3333333
```

[`coherence_spectrum()`](https://bbuchsbaum.github.io/crossform/reference/coherence_spectrum.md)
aggregates that family by scale. Look at the `total` column first.

``` r

spectrum <- coherence_spectrum(
  plan_geometry(example$fit$relation, family, pairing), example$contrast
)
as.data.frame(spectrum)[c("scale", "alpha", "total", "coherence_fraction")]
#>   scale     alpha    total coherence_fraction
#> 1  3.01 0.3333333 51.73602          0.2981374
#> 2  4.30 0.3333333 51.73602          0.2112469
#> 3  6.10 0.3333333 51.73602          0.1572821

per_scale_gap <- max(abs(
  spectrum$total - as.data.frame(spectrum)$alpha * whole_domain$total
))
per_scale_gap
#> [1] 2.842171e-14
```

The `total` column is the same number three times, and it is exactly
`alpha_s` times the whole-domain total — **whatever the data say**
(§3.1). It is a picture of the `weights` argument you passed, not of
spatial scale. Change the weights and the column changes with them,
because you changed it.

The plotting API makes this distinction explicit:

``` r

panel <- catch_refusal(plot(spectrum, which = "profile"))
panel$capability
#> [1] "scale_energy_panel"
```

Use the decomposition view to see how the data divide each scale’s
budget. The next section explains why its coherent fraction can vary
with scale even when total energy is fixed.

## 3. The coherence spectrum is the informative object

What the weights do *not* fix is the split of each scale’s fixed budget
into a coherent part — the rank-one component along the node’s own
weighted common mode — and a configuration part. That share depends on
the data through cross-voxel products, and it is the only scale-resolved
quantity in the table above that the weights do not determine.

It is stronger than “not proportional to α”: the share is **exactly
invariant** to α, because total and coherent are both homogeneous of
degree one under a rescaling of a row, so the ratio cancels α (§3.2).
Measure it rather than believe it: rerun the same three radii under a
lopsided weighting.

``` r

skewed <- coherence_spectrum(
  plan_geometry(
    example$fit$relation,
    compile_frame(searchlights(radii, "conservative", weights = c(0.6, 0.3, 0.1)), domain),
    pairing
  ),
  example$contrast
)
c(
  energy_moved = max(abs(skewed$total - spectrum$total)),
  share_moved = max(abs(skewed$coherence_fraction - spectrum$coherence_fraction))
)
#> energy_moved  share_moved 
#> 4.138882e+01 2.775558e-17
```

Changing α moves total energy while leaving the coherent fraction
unchanged within numerical precision. The fraction still depends on the
chosen radii, frame, and metric, which should accompany its
interpretation. The aggregation record names the weight-determined
column as `alpha_fixed` and the invariant column as `alpha_invariant`.

``` r

spectrum_record <- spectrum$metadata$aggregation
c(
  fixed = spectrum_record$alpha_fixed,
  invariant = spectrum_record$alpha_invariant
)
#>                fixed            invariant 
#>              "total" "coherence_fraction"
```

### Two effect geometries, one budget

The share is a function of `(location, scale)`, not a number per scale,
so a location-wise reading is a table rather than a collapse:
`by_location = TRUE` returns one row per centre per radius, and no
reduction over scales is offered. The fixture’s two planted blocks are
the demonstration — same amplitude, same voxel count, disjoint,
differing only in whether the effect keeps one sign across neighbouring
voxels.

``` r

by_location <- as.data.frame(
  coherence_spectrum(
    plan_geometry(example$fit$relation, family, pairing), example$contrast,
    by_location = TRUE
  )
)
by_location$block <- ifelse(
  by_location$center %in% as.character(example$truth$planted_feature_ids),
  "pattern block",
  ifelse(by_location$center %in% as.character(example$truth$mean_feature_ids),
    "mean block", "elsewhere")
)
knitr::kable(tapply(
  by_location$coherence_fraction,
  list(by_location$block, by_location$scale),
  median, na.rm = TRUE
), digits = 3, caption = "Median valid coherent fraction by block and radius (mm).")
```

|               |  3.01 |   4.3 |   6.1 |
|:--------------|------:|------:|------:|
| elsewhere     | 0.170 | 0.100 | 0.070 |
| mean block    | 0.466 | 0.478 | 0.407 |
| pattern block | 0.058 | 0.088 | 0.044 |

Median valid coherent fraction by block and radius (mm). {.table}

Among nodes with valid fractions, the mean block’s median coherent share
is roughly two fifths at each radius. The pattern block’s median is much
smaller: alternating signs largely cancel under neighbourhood averaging.
Its evidence is primarily configuration at these radii. The elsewhere
median decreases as neighbourhoods grow in this fixture; it is a
descriptive result for these data, not a general scale law for noise.

Two disciplines apply to that table. Shares are **masked, never
clamped**, wherever the components are not a nonnegative partition, so
`NA` in this table means “this node’s decomposition does not admit a
fraction”, not zero. A node’s fraction divides its coherent component by
its own total. A territory’s share of summed coherent mass, introduced
next, has a different denominator. That coherent budget is
**frame-relative** (§4): changing the frame changes the common modes, so
their sum is not a conserved global quantity.

``` r

c(
  rows = nrow(by_location),
  valid_shares = sum(!is.na(by_location$coherence_fraction))
)
#>         rows valid_shares 
#>          840          550
```

## 4. Territory ledgers

Because the map is a budget, it adds.
[`contribution()`](https://bbuchsbaum.github.io/crossform/reference/contribution.md)
is that addition and nothing else: group the rows, sum them, and report
what conserves and what does not.

Grouping is **by row**. A searchlight belongs to the territory its
centre falls in, because a node is not divisible; splitting an
overlapping node’s mass between two territories needs a second partition
that can double-count, and it is deliberately not offered. The
measurement labels of a centre-assigned searchlight frame are the centre
feature identifiers, which is what makes the alignment checkable rather
than assumed.

``` r

ledger_by_territory <- contribution(attribution_map, by = territory)
as.data.frame(ledger_by_territory)
#>     measurement n_rows signed  coherent configuration    total
#> 1     elsewhere    242     NA 11.957596      48.20155 60.15914
#> 2    mean block     19     NA 32.261188      14.95224 47.21342
#> 3 pattern block     19     NA  2.054545      45.78096 47.83551
#>   coherence_fraction
#> 1         0.19876606
#> 2         0.68330542
#> 3         0.04295022
```

``` r

territory_gap <- sum(ledger_by_territory$total) - whole_domain$total
territory_gap
#> [1] 5.684342e-14
```

The territories re-add to the whole-domain total exactly. Two readings a
detection map cannot give. First, the 19 nodes centred in the mean block
— 7 % of the 280 — hold about 30 % of the whole contrast budget, and so
do the 19 centred in the pattern block. Second, those two near-equal
budgets are split completely differently: roughly two thirds of the mean
block’s is coherent against about 4 % of the pattern block’s. Equal
amplitude, equal extent, comparable share of the budget, opposite
geometry.

What the aggregate refuses to give is as important.

``` r

ledger_record <- ledger_by_territory$metadata$aggregation
c(budget_exact = ledger_record$budget_exact, masked = ledger_record$masked)
#> budget_exact       masked 
#>      "total"     "signed"
ledger_record$frame_relative
#> [1] TRUE
```

`signed` is masked rather than summed. A contrast view’s signed marginal
is the local weighted *mean* contrast, already divided by the node’s own
frame mass — a density — and adding densities over a territory is
exactly the error a conservative frame exists to avoid. So a territory
ledger carries no direction; read direction at the node level, from the
signed map.

`coherent` and `configuration` add as arithmetic but their sums are
frame-relative, and the object says so rather than leaving it to the
reader.

## 5. Descriptive spectra from a nonnegative projection

Crossvalidated estimates are **signed**. That is the visible cost of the
cross-partition pairing that removes the noise term, and conservation is
a statement about a signed sum: it does not make the summands
nonnegative.

``` r

c(
  negative_nodes = sum(attribution_map$total < 0),
  nodes = length(attribution_map$total)
)
#> negative_nodes          nodes 
#>             77            280
```

So the arithmetic that treats a node’s value as part of a nonnegative
whole — fractions of the total, cumulative-contribution curves,
effective counts, “top-k explains X %” — is invalid on this layer. The
implied shares can leave `[0, 1]`, and clipping to fix that reintroduces
the bias the pairing removes and destroys the conservation identity
along with it.

Those functionals are legal on a **declared** nonnegative projection,
which is a separate, labelled object that records the mass it moved
(§6).

``` r

geometry <- materialize_geometry(
  plan_geometry(example$fit$relation, attribution, pairing)
)
latent <- latent_geometry(geometry)
latent
#> <effect_latent_geometry>
#>   measurements: 280
#>   component:    total
#>   projection:   psd_projection (eigenvalue truncation at zero)
#>   moved_mass:   4.849 moved, 280 of 280 measurements clipped, max share 1
#>   n_eff:        median 1.03 (range 1 to 2.48), 4 measurements masked
#>   reading:      latent descriptive layer; not for inference
#>  measurement n_eff moved_mass moved_share    root1    root2    root3 root4
#>            1 1.565   0.006469      0.2789 0.012775 0.003955 0.000000     0
#>            2 1.679   0.013687      0.2787 0.026583 0.005275 0.003558     0
#>            3 1.460   0.017744      0.3171 0.030720 0.007485 0.000000     0
#>            4 1.577   0.009267      0.3622 0.012383 0.003934 0.000000     0
#>            5 1.996   0.012853      0.3278 0.013790 0.012567 0.000000     0
#>            6 1.000   0.019517      0.7017 0.008298 0.000000 0.000000     0
#>   ... 274 more measurements
#>   next:         x$cumulative, x$projection
```

``` r

peak <- which.max(attribution_map$total)
all(latent$spectrum >= 0)
#> [1] TRUE
c(
  measurements = length(latent$n_eff),
  clipped = sum(latent$moved_mass > 0),
  masked = sum(is.na(latent$n_eff)),
  median_moved_share = round(median(latent$moved_share), 3)
)
#>       measurements            clipped             masked median_moved_share 
#>             280.00             280.00               4.00               0.34
c(
  peak_n_eff = round(latent$n_eff[[peak]], 3),
  peak_moved_share = round(latent$moved_share[[peak]], 4),
  peak_C1 = round(unname(latent$cumulative[peak, 1L]), 4)
)
#>       peak_n_eff peak_moved_share          peak_C1 
#>           1.0060           0.0024           0.9971
```

The projection is not a formality on this fit: every node clipped
something, the median node moved about a third of its absolute mass, and
a handful of nodes had **no** nonnegative mass at all. Those are masked
— `n_eff` and the whole cumulative row are `NA` — rather than reported
as zero, on the same template the coherence fraction uses. At the
strongest node the projected geometry is effectively one-dimensional and
the projection cost almost nothing; that is a property of that node, not
of the map.

None of this is inference, and none of it may be read back onto the
signed estimates of section 1. The object says so itself: its `reading:`
line prints *latent descriptive layer; not for inference*. And the
projection enters identity rather than sitting beside it.

``` r

substring(latent$receipt$task_partition_id,
  nchar(latent$receipt$task_partition_id) - 14L)
#> [1] "+psd_projection"
query_readout <- catch_refusal(latent_geometry(attribution_map))
query_readout$capability
#> [1] "latent_projection_source"
```

A contrast view is refused because its values are already contracted
against fixed weights: there is no spectrum left to truncate, and
clamping such a value at zero would be a *different* projection moving
different mass. The layer takes a named projection from a closed set
rather than substituting one silently.

A rank budget makes a second descriptive choice: retain at most that
many positive modes. The receipt separates negative mass clipped for PSD
admission from positive mass discarded by the rank budget.

``` r

rank_one <- latent_geometry(geometry, rank = 1)
c(
  peak_clipped_negative = rank_one$clipped_negative_mass[[peak]],
  peak_truncated_positive = rank_one$truncated_positive_mass[[peak]],
  peak_moved = rank_one$moved_mass[[peak]]
)
#>   peak_clipped_negative peak_truncated_positive              peak_moved 
#>             0.009944622             0.012430189             0.022374811
```

This does not show that one dimension predicts independent data. For
that question,
[`fit_geometry()`](https://bbuchsbaum.github.io/crossform/reference/fit_geometry.md)
learns a regularized model-supported form and
[`score_geometry()`](https://bbuchsbaum.github.io/crossform/reference/score_geometry.md)
reads it on separate observations. With the same frozen prediction,
coherent and configuration test inner products still add; the
prediction’s squared norm is subtracted once. Separately projected or
fitted latent components need not add to a projection of their total.
See [Predictive
geometry](https://bbuchsbaum.github.io/crossform/articles/predictive-geometry.md)
([`vignette("predictive-geometry")`](https://bbuchsbaum.github.io/crossform/articles/predictive-geometry.md)
offline) for an executable example.

## 6. Which metrics preserve the budget?

Everything above assumed the node metric is the identity. A metric
composes with the frame by symmetric congruence in the square-root
weights, `K_x = D(sqrt(w_x)) Q D(sqrt(w_x))`, and whether conservation
survives depends entirely on whether `Q` is diagonal.

A dense metric is a `p × p` object and the point here is algebraic, so
this section drops to a nine-feature line domain.

``` r

set.seed(20260820L)
p <- 9L
line <- abstract_domain(p, coordinates = cbind(seq_len(p) - 1, 0),
  feature_ids = paste0("v", seq_len(p)), id = "conservative-frames-metric")
block <- function() {
  value <- matrix(stats::rnorm(3 * p), 3, p)
  rownames(value) <- c("face", "house", "tool")
  value
}
line_relation <- relation(list(run1 = block(), run2 = block()), domain = line)
line_pairing <- cross_partitions(line_relation, independence = "independent")
line_contrast <- c(face = 1, house = -1, tool = 0)
line_frame <- compile_frame(searchlights(1.01, "conservative"), line)
line_global <- compile_frame(whole_brain("none"), line)

line_energy <- function(frame, metric = NULL, composition = "native") {
  contrast_energy(
    plan_geometry(line_relation, frame, line_pairing,
      metric = metric, composition = composition),
    line_contrast
  )
}
```

[`metric_capabilities()`](https://bbuchsbaum.github.io/crossform/reference/metric_capabilities.md)
answers, before anything is read, whether the metric keeps the per-voxel
ledger of section 2 intact.

``` r

diagonal <- neural_metric(diag(stats::runif(p, 0.5, 2.5)), line)
raw <- matrix(stats::rnorm(p * p), p, p)
dense <- neural_metric(crossprod(raw) / p + diag(p), line)
c(
  diagonal = metric_capabilities(diagonal)$feature_additive,
  dense = metric_capabilities(dense)$feature_additive
)
#> diagonal    dense 
#>     TRUE    FALSE
```

A diagonal metric folds into the frame weights:
`D(sqrt(w)) D(q) D(sqrt(w))` is `D(w q)`, still diagonal, still
feature-additive, so `sum_x w_xv = 1` still closes the ledger. A dense
metric picks up off-diagonal terms `S_uv = sum_x sqrt(w_xu w_xv)`, which
column normalization forces to one only on the diagonal. Conservation
fails, and it fails as algebra rather than as tolerance.

``` r

gap <- function(metric = NULL, composition = "native") {
  local_sum <- sum(line_energy(line_frame, metric, composition)$total)
  reference <- line_energy(line_global, metric, composition)$total
  (local_sum - reference) / abs(reference)
}
relative_gaps <- c(
  identity = gap(),
  diagonal = gap(diagonal),
  dense_native = gap(dense),
  dense_whitened = gap(dense, "whitened")
)
signif(relative_gaps, 3)
#>       identity       diagonal   dense_native dense_whitened 
#>       0.00e+00      -1.50e-16      -4.32e-02       1.27e-16
```

The dense metric’s measured gap belongs to this fixture. Its size and
sign can vary with the data, while the algebra explains why column
normalization does not guarantee conservation. A small gap on one
dataset would not establish a conservation law.

`composition = "whitened"` places the frame in whitened coordinates,
`Q^(1/2) D(w_x) Q^(1/2)`, whose sum over a conservative frame is `Q`
exactly for any SPD `Q`. It conserves — and it is **a different
estimand, not a bug fix**.

``` r

native_nodes <- line_energy(line_frame, dense)$total
whitened_nodes <- line_energy(line_frame, dense, "whitened")$total
node_gap <- max(abs(native_nodes - whitened_nodes)) / max(abs(native_nodes))
round(node_gap, 3)
#> [1] 0.247
```

Under the native composition a node weights *features* and measures them
in the `Q` geometry. Under the whitened composition it weights *whitened
coordinates*, which are spatially delocalized whenever `Q` is dense —
the node’s support is no longer its support. Here the two disagree by
roughly a quarter of the largest node value. Both are defensible; they
answer different questions, and the package makes you name which.

That naming is enforced through plan identity, so a whitened analysis is
never mistaken for a native one after the fact.

``` r

native_plan <- plan_geometry(line_relation, line_frame, line_pairing,
  metric = dense)
whitened_plan <- plan_geometry(line_relation, line_frame, line_pairing,
  metric = dense, composition = "whitened")
!identical(native_plan$scientific_plan_id, whitened_plan$scientific_plan_id)
#> [1] TRUE
```

The composition and the root convention — `"whitened"` means the
symmetric positive semidefinite root, since conservation holds for *any*
root `R R' = Q` while the node values do not — both enter the plan’s
scientific identity. A conservation certificate is therefore not
sufficient evidence that two whitened analyses computed the same thing.
And the switch must never be applied silently to repair a failed
conservation check: that is a change of estimand, and it belongs in the
record.

## 7. Reporting an attribution map

**No inference.** Conservation is a point-estimate law.
`sum_x theta_x = theta_Omega` says nothing about uncertainty: node
estimates of an overlapping frame can be correlated, variances do not
generally add, and per-node standard errors must never be summed to put
an error bar on a conserved budget (§7.6a). The cross-node sampling
covariance that would license one does not exist yet;
`sampling_covariance(scope = "cross_measurement")` refuses by capability
rather than returning per-node margins that could be mistaken for a
joint block. Until then a conservative attribution map is a point ledger
reported without inference.

**No cross-voxel content in a total.** Section 2 is a caveat as much as
a demonstration: under a feature-additive metric, nothing in a `total`
column at any radius is multivariate across space. The multivariate
content of this package lives in the coherent/configuration split and in
the RDM and RSA views, not in the ledger.

**Advanced tier.** This is an advanced-tier topic and uses advanced-tier
exports:
[`frame_conservation()`](https://bbuchsbaum.github.io/crossform/reference/frame_conservation.md),
[`materialize_geometry()`](https://bbuchsbaum.github.io/crossform/reference/materialize_geometry.md),
[`latent_geometry()`](https://bbuchsbaum.github.io/crossform/reference/latent_geometry.md),
[`neural_metric()`](https://bbuchsbaum.github.io/crossform/reference/neural_metric.md),
and
[`metric_capabilities()`](https://bbuchsbaum.github.io/crossform/reference/metric_capabilities.md).
Everything else here is core:
[`example_fmri_effects()`](https://bbuchsbaum.github.io/crossform/reference/example_fmri_effects.md),
[`compile_frame()`](https://bbuchsbaum.github.io/crossform/reference/compile_frame.md),
[`searchlights()`](https://bbuchsbaum.github.io/crossform/reference/searchlights.md),
[`voxelwise()`](https://bbuchsbaum.github.io/crossform/reference/voxelwise.md),
[`whole_brain()`](https://bbuchsbaum.github.io/crossform/reference/whole_brain.md),
[`abstract_domain()`](https://bbuchsbaum.github.io/crossform/reference/abstract_domain.md),
[`relation()`](https://bbuchsbaum.github.io/crossform/reference/relation.md),
[`cross_partitions()`](https://bbuchsbaum.github.io/crossform/reference/cross_partitions.md),
[`plan_geometry()`](https://bbuchsbaum.github.io/crossform/reference/plan_geometry.md),
[`contrast_energy()`](https://bbuchsbaum.github.io/crossform/reference/contrast_energy.md),
[`contribution()`](https://bbuchsbaum.github.io/crossform/reference/contribution.md),
[`coherence_spectrum()`](https://bbuchsbaum.github.io/crossform/reference/coherence_spectrum.md),
and
[`catch_refusal()`](https://bbuchsbaum.github.io/crossform/reference/catch_refusal.md).

## See also

- [`vignette("interpreting-results")`](https://bbuchsbaum.github.io/crossform/articles/interpreting-results.md)
  — reading a detection map, and the three interpretive traps that make
  a coherent share look like a finding when it is an artefact of the
  frame. Read it before reporting any share from this article.
- [`vignette("introduction")`](https://bbuchsbaum.github.io/crossform/articles/introduction.md)
  — the same fixture as a detection-first workflow, and where the
  planted blocks come from.
- [`design/conservative-geometry-contract.md`](https://github.com/bbuchsbaum/crossform/blob/main/design/conservative-geometry-contract.md)
  — the normative document: the two estimands, the conservation theorem,
  the smoothed-ledger claim and its metric precondition, the
  α-invariance argument, and the numerical tolerances used by this
  guide’s checks.
- [Haxby 2001
  exemplar](https://github.com/bbuchsbaum/crossform/tree/main/exemplars/haxby2001),
  script `07-conservative-geometry.R` and section 07 of its README — the
  same four instruments on one subject’s ventral temporal cortex, with
  committed receipts.
- [Population
  form](https://bbuchsbaum.github.io/crossform/articles/population-form.md)
  — carry each participant’s conservative ledger to shared group nodes
  and estimate a group model. This layer is experimental.
