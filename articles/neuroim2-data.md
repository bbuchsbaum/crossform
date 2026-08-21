# Bring your own neuroim2 data

You already have per-run condition estimates as
[neuroim2](https://github.com/bbuchsbaum/neuroim2) images and a brain
mask. This guide is the shortest correct path from those objects to a
contrast map written back as a `NeuroVol`. It covers four things that
are easy to get wrong: the orientation of the matrices
[`relation()`](https://bbuchsbaum.github.io/crossform/reference/relation.md)
expects, which searchlight provider to use, the difference between a
*measurement* and a *feature*, and how to expand region results back to
voxels.

The `crossform` side is the same as in
[`vignette("introduction", package = "crossform")`](https://bbuchsbaum.github.io/crossform/articles/introduction.md);
only the input and output adapters are new.

``` r

library(neuroim2)
```

## A small example volume

Real data would be read with
[`neuroim2::read_vec()`](https://bbuchsbaum.github.io/neuroim2/reference/read_vec.html)
and
[`neuroim2::read_vol()`](https://bbuchsbaum.github.io/neuroim2/reference/read_vol.html).
So that this vignette is self-contained and fast, the same objects are
generated here: a 10-by-10-by-8 volume at 3 mm spacing with a boxy mask,
and one 4-D image holding 3 runs by 4 conditions with a multivariate
animacy pattern planted in a small central blob.

``` r

set.seed(11)
dims <- c(10L, 10L, 8L)
spacing <- c(3, 3, 3)
conditions <- c("face", "house", "cat", "chair")
n_run <- 3L
n_volume <- n_run * length(conditions)

volume_space <- NeuroSpace(dims, spacing = spacing)
mask_array <- array(0L, dims)
mask_array[2:9, 2:9, 2:7] <- 1L
mask <- NeuroVol(mask_array, volume_space)

# One planted multivariate pattern, present in every run.
blob <- array(0, dims)
blob[5:6, 5:6, 4:5] <- 1
pattern <- array(rnorm(prod(dims)), dims)
animacy <- c(face = 0.5, house = -0.5, cat = 0.5, chair = -0.5)

estimates <- array(rnorm(prod(dims) * n_volume), c(dims, n_volume))
k <- 0L
for (run in seq_len(n_run)) {
  for (condition in conditions) {
    k <- k + 1L
    estimates[, , , k] <- estimates[, , , k] +
      3 * animacy[[condition]] * blob * pattern
  }
}

runs_vec <- DenseNeuroVec(
  estimates,
  NeuroSpace(c(dims, n_volume), spacing = spacing)
)
runs_vec
#> <DenseNeuroVec> [82.2 Kb] 
#> ── Spatial ───────────────────────────────────────────────────────────────────── 
#>   Dimensions    : 10 x 10 x 8 (12 timepoints)
#>   Spacing       : 3 x 3 x 3
#>   Origin        : 0, 0, 0
#>   Orientation   : RAS
#> ── Data ──────────────────────────────────────────────────────────────────────── 
#>   Mean +/- SD   : 0.001 +/- 1.023 (t=1)
#>   Label         : none
```

The volumes are ordered run-major: volumes 1 to 4 are run 1’s four
conditions, volumes 5 to 8 are run 2’s, and so on. Keep that bookkeeping
explicit — `crossform` will check the condition names but cannot check
that you sliced the runs correctly.

## The domain fixes the voxel order

``` r

domain <- neuroim2_volume_domain(mask)
c(
  features = domain$n_features,
  first_voxel_index = domain$feature_ids[1],
  last_voxel_index = domain$feature_ids[domain$n_features]
)
#>          features first_voxel_index  last_voxel_index 
#>               384               112               689
```

[`neuroim2_volume_domain()`](https://bbuchsbaum.github.io/crossform/reference/neuroim2_volume_domain.md)
records the mask’s **stable full-volume indices**
(`domain$feature_ids`), the physical coordinates, the spacing, and a
hash of the full `NeuroSpace`. Those indices are the contract between
every later step: they define which column of a response matrix is which
voxel, and they let
[`as_neurovol()`](https://bbuchsbaum.github.io/crossform/reference/as_neurovol.md)
put results back without interpolation, smoothing, or any
reinterpretation of the coordinate frame.

``` r

identical(domain$feature_ids, which(mask_array != 0))
#> [1] TRUE
```

## Convert images to run-by-condition matrices

[`relation()`](https://bbuchsbaum.github.io/crossform/reference/relation.md)
wants a named list with one matrix per run, each
**condition-by-feature**: rows are the named conditions, columns are the
domain’s features in `domain$feature_ids` order.

`neuroim2::series(x, i)` extracts the time series at a set of voxel
indices and returns a **volume-by-voxel** matrix. That is already the
orientation `crossform` wants, so no transpose is needed — just index
the rows for one run’s volumes and index the columns with the domain’s
feature ids.

``` r

run_matrices <- lapply(seq_len(n_run), function(run) {
  volumes <- ((run - 1L) * length(conditions) + 1L):(run * length(conditions))
  m <- series(runs_vec, domain$feature_ids)[volumes, , drop = FALSE]
  rownames(m) <- conditions
  m
})
names(run_matrices) <- paste0("run", seq_len(n_run))

dim(run_matrices$run1)
#> [1]   4 384
rownames(run_matrices$run1)
#> [1] "face"  "house" "cat"   "chair"
```

If your estimates are a *list of `NeuroVol`s* rather than one
`NeuroVec`, index each volume at the same feature ids and stack the
rows. `vol[ids]` returns the values at those full-volume indices:

``` r

condition_vols <- lapply(seq_len(length(conditions)), function(k) {
  NeuroVol(estimates[, , , k], volume_space)
})

from_vols <- t(vapply(
  condition_vols,
  function(vol) as.numeric(vol[domain$feature_ids]),
  numeric(domain$n_features)
))
rownames(from_vols) <- conditions

isTRUE(all.equal(from_vols, run_matrices$run1))
#> [1] TRUE
```

Both routes go through `domain$feature_ids`, which is what makes them
agree. Do not use `as.array(vol)[mask_array != 0]` unless you are
certain the mask you are subsetting with is the same one that built the
domain.

With the matrices in hand the relation is ordinary `crossform`:

``` r

rel <- relation(
  run_matrices,
  effects = effect_space(conditions, units = "beta"),
  domain = domain
)
```

## Two searchlight providers, one support index

There are two ways to get spheres over this domain, and they are not
redundant.

``` r

frame_neuroim2 <- neuroim2_searchlights(mask, radius = 6, domain = domain)
frame_grid <- compile_frame(searchlights(radius = 6), domain)

c(
  neuroim2_measurements = nrow(frame_neuroim2$index),
  grid_measurements = nrow(frame_grid$index)
)
#> neuroim2_measurements     grid_measurements 
#>                   384                   384
```

- `neuroim2_searchlights(mask, radius)` delegates to
  [`neuroim2::searchlight_indices()`](https://bbuchsbaum.github.io/neuroim2/reference/searchlight_indices.html).
  Use it when the same neighborhoods must be shared with other
  neuroim2-based tools: the sphere membership is literally the one that
  provider returns, and the frame records the upstream provider and
  commit in its specification.
- `compile_frame(searchlights(radius), domain)` uses `crossform`’s own
  grid stencil. It needs no neuroim2 at all, so it also works on a
  [`volume_domain()`](https://bbuchsbaum.github.io/crossform/reference/volume_domain.md)
  built from a plain array — which is what the rest of the documentation
  uses.

Both produce the same kind of object: an additive `effect_frame`
carrying a **support index** (which features each measurement reads) and
a sparse measurement-by-feature weight matrix. On a regular grid with a
Euclidean ball they agree, and you can check it rather than trust it:

``` r

identical(
  as.matrix(frame_neuroim2$weights),
  as.matrix(frame_grid$weights)
)
#> [1] TRUE
```

``` r

frame_neuroim2$specification
#> $kind
#> [1] "neuroim2_searchlights"
#> 
#> $radius
#> [1] 6
#> 
#> $units
#> [1] "mm"
#> 
#> $nonzero
#> [1] TRUE
#> 
#> $upstream_commit
#> [1] "77b1ddb"
range(Matrix::rowSums(frame_neuroim2$weights != 0))
#> [1] 11 33
```

The radius is in millimeters for both, and both clip spheres at the mask
edge, so the sphere sizes vary. Use whichever provider matches the rest
of your pipeline; the plan below does not care which one produced the
frame.

## Plan and read the contrast

``` r

plan <- plan_geometry(
  rel,
  at = frame_neuroim2,
  over = cross_partitions(
    rel,
    independence = "independent",
    generalizes_over = "run"
  )
)

effect <- contrast_energy(plan, animacy)
peak <- which.max(effect$total)
c(
  peak_measurement = peak,
  peak_voxel_index = domain$feature_ids[peak],
  total_energy = effect$total[peak]
)
#> peak_measurement peak_voxel_index     total_energy 
#>       220.000000       445.000000         1.869888
```

The planted blob is recovered. Searchlights whose support touches a
planted voxel carry most of the reproducible energy:

``` r

planted_features <- which(domain$feature_ids %in% which(blob != 0))
membership <- frame_neuroim2$weights != 0
touches_blob <- Matrix::rowSums(membership[, planted_features, drop = FALSE]) > 0

c(
  searchlights_touching_blob = sum(touches_blob),
  mean_energy_touching = mean(effect$total[touches_blob]),
  mean_energy_elsewhere = mean(effect$total[!touches_blob]),
  peak_touches_blob = unname(touches_blob[peak])
)
#> searchlights_touching_blob       mean_energy_touching 
#>                88.00000000                 0.65743727 
#>      mean_energy_elsewhere          peak_touches_blob 
#>                -0.04201012                 1.00000000
```

Away from the blob the estimate is centered near zero and is negative
about as often as not. That is the crossvalidated estimator behaving
correctly, not a bug; see
[`vignette("introduction", package = "crossform")`](https://bbuchsbaum.github.io/crossform/articles/introduction.md)
for why those values are kept rather than truncated.

## Measurements are not features

This is the distinction that decides how you write results back.

- A **feature** is a voxel in the domain, addressed by
  `domain$feature_ids`. There are `domain$n_features` of them.
- A **measurement** is a row of the frame — one searchlight, one region,
  one whole-brain summary. Results (`effect$total`, `rdm()$values`,
  `rsa()$coefficients`) have one row per **measurement**.

[`as_neurovol()`](https://bbuchsbaum.github.io/crossform/reference/as_neurovol.md)
takes one finite value per **feature**. With searchlights the two
coincide: the frame has exactly one measurement per domain feature, in
the same order, because each sphere is labeled by its center voxel.

``` r

c(
  measurements = nrow(frame_neuroim2$index),
  features = domain$n_features,
  index_matches_features =
    identical(frame_neuroim2$index$measurement, domain$feature_ids)
)
#>           measurements               features index_matches_features 
#>                    384                    384                      1
```

So a searchlight map can be written straight back:

``` r

energy_vol <- as_neurovol(
  effect$total,
  mask,
  domain,
  label = "animacy total energy"
)
energy_vol
#> <DenseNeuroVol> [13 Kb] 
#> ── Spatial ───────────────────────────────────────────────────────────────────── 
#>   Dimensions    : 10 x 10 x 8
#>   Spacing       : 3 x 3 x 3 mm
#>   Origin        : 0, 0, 0
#>   Orientation   : RAS
#> ── Data ──────────────────────────────────────────────────────────────────────── 
#>   Range         : [-0.452, 1.870]
#>   NAs           : 416

# Values land at exactly the domain's voxels; everything else is `fill`.
isTRUE(all.equal(as.numeric(energy_vol[domain$feature_ids]), effect$total))
#> [1] TRUE
all(is.na(as.array(energy_vol)[mask_array == 0]))
#> [1] TRUE
```

[`as_neurovol()`](https://bbuchsbaum.github.io/crossform/reference/as_neurovol.md)
re-derives the domain from `mask` and refuses if the geometry does not
match the domain that produced the values, so a mask/result mismatch is
an error rather than a silently misplaced map. It also requires every
value to be finite: it will not guess what a missing measurement means.

## Regions from an atlas, and expanding them back

[`regions()`](https://bbuchsbaum.github.io/crossform/reference/regions.md)
takes one label per **feature** — a vector as long as
`domain$feature_ids`, not an image. Read an atlas `NeuroVol` at the
domain’s voxel indices to get it:

``` r

atlas_array <- array(0L, dims)
atlas_array[2:4, 2:9, 2:7] <- 1L
atlas_array[5:7, 2:9, 2:7] <- 2L
atlas_array[8:9, 2:9, 2:7] <- 3L
atlas_array[mask_array == 0L] <- 0L
atlas <- NeuroVol(atlas_array, volume_space)

region_labels <- as.integer(atlas[domain$feature_ids])
table(region_labels)
#> region_labels
#>   1   2   3 
#> 144 144  96

frame_regions <- compile_frame(regions(region_labels), domain)
frame_regions$index
#>   measurement
#> 1           1
#> 2           2
#> 3           3
```

``` r

plan_regions <- plan_geometry(
  rel,
  at = frame_regions,
  over = cross_partitions(
    rel,
    independence = "independent",
    generalizes_over = "run"
  )
)
region_effect <- contrast_energy(plan_regions, animacy)
round(region_effect$total, 3)
#> [1] -0.065  0.408 -0.051
```

Now there are only 3 measurements but still `domain$n_features`
features, so `as_neurovol(region_effect$total, ...)` would be an error.
Expand the region values to voxels first. The general route uses the
frame’s own membership pattern, which works for any additive frame:

``` r

region_membership <- frame_regions$weights != 0
voxel_values <- as.numeric(
  Matrix::crossprod(region_membership, region_effect$total)
)

length(voxel_values)
#> [1] 384
region_vol <- as_neurovol(voxel_values, mask, domain, label = "region energy")
```

Use the membership pattern (`weights != 0`), not `weights` itself: under
the default `"local"` normalization the weights are `1/region_size`, so
`crossprod(weights, values)` would rescale the numbers you are trying to
display.

For disjoint labels the same thing is a lookup, which is easier to read:

``` r

by_lookup <- region_effect$total[
  match(region_labels, frame_regions$index$measurement)
]
isTRUE(all.equal(by_lookup, voxel_values))
#> [1] TRUE
```

One caveat. If your atlas does not cover the whole mask,
[`regions()`](https://bbuchsbaum.github.io/crossform/reference/regions.md)
drops the unlabeled features, and the lookup above yields `NA` there —
which
[`as_neurovol()`](https://bbuchsbaum.github.io/crossform/reference/as_neurovol.md)
refuses. Build the domain from the atlas coverage
(`neuroim2_volume_domain(NeuroVol((atlas_array != 0) * 1, volume_space))`)
so that domain and results describe the same voxels, rather than filling
missing regions with a number that will later be read as data.

## Tidy output

Every result view has an
[`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html) method.
The measurement index is carried through, so for a searchlight frame the
`measurement` column is the center voxel’s full-volume index:

``` r

head(as.data.frame(effect), 3)
#>   measurement     signed    coherent configuration     total coherence_fraction
#> 1         112 0.21241731  0.04510716    0.08493888 0.1300460          0.3468553
#> 2         113 0.10493450 -0.01416982    0.26524800 0.2510782                 NA
#> 3         114 0.01008603 -0.02316588    0.15649282 0.1333269                 NA
```

``` r

distances <- rdm(plan)
head(as.data.frame(distances)[, 1:4], 3)
#>   measurement face - house  face - cat face - chair
#> 1         112   0.08167025  0.02844973   0.21860346
#> 2         113   0.49200015 -0.08265131   0.29447256
#> 3         114   0.02929930 -0.22957438   0.05166341
```

That frame joins directly to anything else keyed by voxel index, and
[`as_neurovol()`](https://bbuchsbaum.github.io/crossform/reference/as_neurovol.md)
remains available for the map itself.

## Where to go next

- [`vignette("introduction", package = "crossform")`](https://bbuchsbaum.github.io/crossform/articles/introduction.md)
  for what the energies, RDMs, and RSA coefficients mean.
- [`vignette("from-observations", package = "crossform")`](https://bbuchsbaum.github.io/crossform/articles/from-observations.md)
  if you have scan-by-feature responses rather than per-run condition
  estimates, and want the residual channel that analytic standard errors
  require.
- [`?neuroim2_volume_domain`](https://bbuchsbaum.github.io/crossform/reference/neuroim2_volume_domain.md),
  [`?neuroim2_searchlights`](https://bbuchsbaum.github.io/crossform/reference/neuroim2_searchlights.md),
  and
  [`?as_neurovol`](https://bbuchsbaum.github.io/crossform/reference/as_neurovol.md)
  for the exact adapter contracts, including what each one refuses.
