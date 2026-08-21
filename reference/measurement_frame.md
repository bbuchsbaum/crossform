# Declare identified neural measurements

`measurement_frame()` preserves the oriented measurement legs needed for
cross-region forms. An additive frame can be adapted as total, coherent,
or coherent/configuration measurements. A named list of fixed matrices
can be used for prespecified measurements. Learned operators are
intentionally not accepted by this convenience constructor because they
require frozen training provenance.

## Usage

``` r
measurement_frame(
  x,
  domain = NULL,
  mode = c("total", "coherent", "coherent_configuration"),
  id = "fixed-measurements:v1",
  units = "arbitrary"
)
```

## Arguments

- x:

  A compiled additive `effect_frame`, or a uniquely named list of fixed
  measurement matrices.

- domain:

  Required for a list of matrices; the exact source `effect_domain`.

- mode:

  For an additive frame, one of `"total"`, `"coherent"`, or
  `"coherent_configuration"`.

- id:

  Stable identity used for custom measurement axes.

- units:

  One unit or one per row of each custom operator.

## Value

An `effect_measurement_frame` with one oriented `$legs` entry per node,
the `$node_ids` naming them, the stacked `$frame_operator`, and
`$coverage`, `$injectivity`, and `$dual` diagnostics used by
[`reconstruct_evidence()`](https://bbuchsbaum.github.io/crossform/reference/reconstruct_evidence.md).

## Details

Version 0.1 exposes this as a deliberately small-node dense path.
Construction fails before dense conversion when its estimated frame and
leg payload exceeds 256 MiB. Brain-scale work should use the
support-local geometry-plan path until matrix-free measurement frames
are qualified.

## See also

[`edge_frame()`](https://bbuchsbaum.github.io/crossform/reference/edge_frame.md)
to request node pairs,
[`measurement_form()`](https://bbuchsbaum.github.io/crossform/reference/measurement_form.md)
to evaluate them, and
[`additive_frame()`](https://bbuchsbaum.github.io/crossform/reference/additive_frame.md)
or
[`compile_frame()`](https://bbuchsbaum.github.io/crossform/reference/compile_frame.md)
for the additive frames this can adapt.

Other neural domains and frames:
[`abstract_domain()`](https://bbuchsbaum.github.io/crossform/reference/abstract_domain.md),
[`additive_frame()`](https://bbuchsbaum.github.io/crossform/reference/additive_frame.md),
[`as_neurovol()`](https://bbuchsbaum.github.io/crossform/reference/as_neurovol.md),
[`compile_frame()`](https://bbuchsbaum.github.io/crossform/reference/compile_frame.md),
[`edge_frame()`](https://bbuchsbaum.github.io/crossform/reference/edge_frame.md),
[`frame_conservation()`](https://bbuchsbaum.github.io/crossform/reference/frame_conservation.md),
[`frame_family()`](https://bbuchsbaum.github.io/crossform/reference/frame_family.md),
[`neuroim2_searchlights()`](https://bbuchsbaum.github.io/crossform/reference/neuroim2_searchlights.md),
[`neuroim2_volume_domain()`](https://bbuchsbaum.github.io/crossform/reference/neuroim2_volume_domain.md),
[`regions()`](https://bbuchsbaum.github.io/crossform/reference/regions.md),
[`searchlights()`](https://bbuchsbaum.github.io/crossform/reference/searchlights.md),
[`volume_domain()`](https://bbuchsbaum.github.io/crossform/reference/volume_domain.md),
[`voxelwise()`](https://bbuchsbaum.github.io/crossform/reference/voxelwise.md),
[`whole_brain()`](https://bbuchsbaum.github.io/crossform/reference/whole_brain.md)

## Examples

``` r
# Two scalar regional measurements, each a fixed oriented row of weights.
native <- abstract_domain(4, id = "demo:native:v1")
nodes <- measurement_frame(
  list(
    anterior = matrix(c(1, 0, 0, 0), 1),
    posterior = matrix(c(0, 1, 0, 0), 1)
  ),
  domain = native, id = "demo:regional-means:v1"
)
nodes$node_ids
#> [1] "anterior"  "posterior"
nodes$legs$anterior$operator
#>      [,1] [,2] [,3] [,4]
#> [1,]    1    0    0    0

# An additive frame can be adapted instead. In coherent/configuration
# mode each measurement keeps its weighted-mean direction and the
# orthogonal remainder as separate components.
additive <- additive_frame(
  matrix(c(1, 2, 1, 0, 0, 1, 2, 1), 2, 4, byrow = TRUE), domain = native
)
decomposed <- measurement_frame(additive, mode = "coherent_configuration")
decomposed$node_ids
#> [1] "measurement_1" "measurement_2"

# An additive frame already names its domain, so supplying one is refused.
refused <- try(
  measurement_frame(additive, domain = native), silent = TRUE
)
conditionMessage(attr(refused, "condition"))
#> [1] "An additive frame already identifies its neural domain."
```
