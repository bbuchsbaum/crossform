# Build a hard nearest-centre transport from node coordinates

`anatomical_transport()` is the transport a registration or a
parcellation induces: each native node is assigned whole to one group
node, and nothing is estimated from the response data.

## Usage

``` r
anatomical_transport(
  native_coords,
  group_coords,
  semantics,
  radius = NULL,
  native_index = NULL,
  group_index = NULL,
  provenance = list(),
  row_mass = NULL
)
```

## Arguments

- native_coords:

  Native node centres, one row per native node.

- group_coords:

  Group node centres, one row per group node, with the same number of
  coordinate columns as `native_coords`.

- semantics:

  `"budget"` or `"density"`. Required; see
  [`location_transport()`](https://bbuchsbaum.github.io/crossform/reference/location_transport.md).

- radius:

  Optional positive assignment radius. A native node whose closest group
  centre is further than this goes entirely to the sink.

- native_index, group_index:

  Optional index tables or node-identifier vectors, as in
  [`location_transport()`](https://bbuchsbaum.github.io/crossform/reference/location_transport.md).
  The coordinates are appended to them as `coord1`, `coord2`, ...
  columns, so the displacement diagnostics can be recomputed from the
  transport alone.

- provenance:

  Optional extra provenance entries, such as the atlas or warp identity.
  `method` is fixed at `"anatomical"` and `details` records the
  assignment rule unless the caller supplies its own.

- row_mass:

  Optional declared positive territory measure, one entry per native
  node.

## Value

An `effect_location_transport`.

## The assignment rule

Every native node is assigned, with mass exactly one, to the group node
whose centre is closest in Euclidean distance in the coordinate system
both matrices are already expressed in. No registration is performed and
no image is resampled: the transport maps nodes, not voxels. Ties go to
the lowest group-node position, so the operator is a deterministic
function of its inputs. When `radius` is supplied and the closest group
centre is further away than `radius`, the native node is left unassigned
and its entire mass goes to the sink; without `radius` every native node
is assigned and the sink is identically zero but still materialized.

The result is a hard assignment, so every row has entropy zero. A
partial volume warp does not, which is why the population diagnostics
report row entropy rather than assuming it.

## See also

[`location_transport()`](https://bbuchsbaum.github.io/crossform/reference/location_transport.md)
for the general constructor and the field semantics.

Other population transports:
[`estimate_population()`](https://bbuchsbaum.github.io/crossform/reference/estimate_population.md),
[`external_transport()`](https://bbuchsbaum.github.io/crossform/reference/external_transport.md),
[`heterogeneity()`](https://bbuchsbaum.github.io/crossform/reference/heterogeneity.md),
[`location_transport()`](https://bbuchsbaum.github.io/crossform/reference/location_transport.md),
[`materialize_population()`](https://bbuchsbaum.github.io/crossform/reference/materialize_population.md),
[`plan_population()`](https://bbuchsbaum.github.io/crossform/reference/plan_population.md),
[`population_prevalence()`](https://bbuchsbaum.github.io/crossform/reference/population_prevalence.md),
[`population_uncertainty()`](https://bbuchsbaum.github.io/crossform/reference/population_uncertainty.md),
[`population_views`](https://bbuchsbaum.github.io/crossform/reference/population_views.md),
[`transport_values()`](https://bbuchsbaum.github.io/crossform/reference/transport_values.md)

## Examples

``` r
# Five native nodes on a line, two group centres at 0 and 5, radius 2.
transport <- anatomical_transport(
  native_coords = cbind(c(0, 1, 4, 5, 9)),
  group_coords = cbind(c(0, 5)),
  semantics = "budget", radius = 2,
  provenance = list(atlas = "example:line-atlas:v1")
)
as.matrix(transport$matrix)
#>      [,1] [,2] [,3]
#> [1,]    1    0    0
#> [2,]    1    0    0
#> [3,]    0    1    0
#> [4,]    0    1    0
#> [5,]    0    0    1

# The unassigned node is visible as sink territory, not as absence.
transport
#> <effect_location_transport>
#>   nodes:      5 native -> 2 group + sink
#>   semantics:  budget
#>   sink:       mass 1 of 5 rows, 20.0% of territory
#>   provenance: anatomical (cross-fit: none)
#>   built:      nearest group centre within radius 2, ties to the lowest gr...
#>   signature:  sha256:df387cad865b...
```
