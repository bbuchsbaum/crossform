# Declare requested neural measurement edges

Every edge is explicit. The constructor never creates all pairs on the
user's behalf. `from` and `to` name neural measurement nodes; the words
`left` and `right` are reserved for experimental relation sides in
[`measurement_form()`](https://bbuchsbaum.github.io/crossform/reference/measurement_form.md).

## Usage

``` r
edge_frame(from, to, frame, to_frame = frame, weight = NULL)
```

## Arguments

- from, to:

  Equal-length vectors of node identifiers.

- frame:

  Measurement frame containing `from` nodes.

- to_frame:

  Measurement frame containing `to` nodes. Defaults to `frame`.

- weight:

  Optional finite edge weights recorded as part of edge identity.

## Value

An `effect_edge_frame` for the `between` argument, holding the
`$from_frame` and `$to_frame` it draws nodes from, the requested
`$edges` table (`left`, `right`, `weight`), and a `$signature`.

## See also

[`measurement_frame()`](https://bbuchsbaum.github.io/crossform/reference/measurement_frame.md)
for the nodes,
[`measurement_form()`](https://bbuchsbaum.github.io/crossform/reference/measurement_form.md)
which consumes this edge set, and
[`reconstruct_evidence()`](https://bbuchsbaum.github.io/crossform/reference/reconstruct_evidence.md),
which needs every directed pair.

Other neural domains and frames:
[`abstract_domain()`](https://bbuchsbaum.github.io/crossform/reference/abstract_domain.md),
[`additive_frame()`](https://bbuchsbaum.github.io/crossform/reference/additive_frame.md),
[`as_neurovol()`](https://bbuchsbaum.github.io/crossform/reference/as_neurovol.md),
[`compile_frame()`](https://bbuchsbaum.github.io/crossform/reference/compile_frame.md),
[`frame_conservation()`](https://bbuchsbaum.github.io/crossform/reference/frame_conservation.md),
[`frame_family()`](https://bbuchsbaum.github.io/crossform/reference/frame_family.md),
[`measurement_frame()`](https://bbuchsbaum.github.io/crossform/reference/measurement_frame.md),
[`neuroim2_searchlights()`](https://bbuchsbaum.github.io/crossform/reference/neuroim2_searchlights.md),
[`neuroim2_volume_domain()`](https://bbuchsbaum.github.io/crossform/reference/neuroim2_volume_domain.md),
[`regions()`](https://bbuchsbaum.github.io/crossform/reference/regions.md),
[`searchlights()`](https://bbuchsbaum.github.io/crossform/reference/searchlights.md),
[`volume_domain()`](https://bbuchsbaum.github.io/crossform/reference/volume_domain.md),
[`voxelwise()`](https://bbuchsbaum.github.io/crossform/reference/voxelwise.md),
[`whole_brain()`](https://bbuchsbaum.github.io/crossform/reference/whole_brain.md)

## Examples

``` r
native <- abstract_domain(4, id = "demo:native:v1")
nodes <- measurement_frame(
  list(
    anterior = matrix(c(1, 0, 0, 0), 1),
    posterior = matrix(c(0, 1, 0, 0), 1)
  ),
  domain = native, id = "demo:regional-means:v1"
)

# Request exactly the edges the question needs. A correlation view also
# needs the two self-pairs, because they supply its denominator.
pairs <- expand.grid(
  from = nodes$node_ids, to = nodes$node_ids, stringsAsFactors = FALSE
)
between <- edge_frame(pairs$from, pairs$to, nodes)
between$edges$edges
#>        left     right weight
#> 1  anterior  anterior      1
#> 2 posterior  anterior      1
#> 3  anterior posterior      1
#> 4 posterior posterior      1

# A seed-to-target edge set is just a shorter list; nothing is created on
# your behalf.
edge_frame("anterior", "posterior", nodes)$edges$edges
#>       left     right weight
#> 1 anterior posterior      1

# Node names are checked against the frame, so a typo cannot become a
# silently missing edge.
refused <- try(edge_frame("anterior", "postrior", nodes), silent = TRUE)
conditionMessage(attr(refused, "condition"))
#> [1] "Measurement edges must be unique explicit ordered frame-node pairs."
```
