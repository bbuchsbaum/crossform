# Admit a transport built outside crossform

`crossform` does not learn transports. `external_transport()` is the
door a transport built elsewhere — by a registration pipeline, a
parcellation tool, or a functional alignment method — comes through: a
thin wrapper over
[`location_transport()`](https://bbuchsbaum.github.io/crossform/reference/location_transport.md)
that takes the group block as it stands, reads node identifiers from its
dimension names when they are there, and requires the provenance that
makes the operator auditable.

## Usage

``` r
external_transport(
  P,
  semantics,
  provenance,
  native_index = NULL,
  group_index = NULL,
  row_mass = NULL,
  tolerance = 1e-09
)
```

## Arguments

- P:

  The `n_native` by `m` nonnegative group block, as a base matrix or a
  `Matrix`. Dimension names, when present, become the default node
  identifiers.

- semantics:

  `"budget"` or `"density"`. Required; see
  [`location_transport()`](https://bbuchsbaum.github.io/crossform/reference/location_transport.md).

- provenance:

  A list carrying at least `details`; `method` defaults to `"external"`.

- native_index, group_index:

  Optional index tables or node-identifier vectors overriding the
  dimension names.

- row_mass:

  Optional declared positive territory measure, one entry per native
  node.

- tolerance:

  Positive row-sum tolerance, as in
  [`location_transport()`](https://bbuchsbaum.github.io/crossform/reference/location_transport.md).

## Value

An `effect_location_transport`.

## Details

Declaring `method = "external"` says only that the operator was built
outside this package. A transport that was fitted to response data is
`"functional"` whichever program fitted it, and must name its
`cross_fit` partitions; `"external"` is not a way around that.

## See also

[`location_transport()`](https://bbuchsbaum.github.io/crossform/reference/location_transport.md),
[`anatomical_transport()`](https://bbuchsbaum.github.io/crossform/reference/anatomical_transport.md).

Other population transports:
[`anatomical_transport()`](https://bbuchsbaum.github.io/crossform/reference/anatomical_transport.md),
[`estimate_population()`](https://bbuchsbaum.github.io/crossform/reference/estimate_population.md),
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
# A partial-volume operator produced by some other program, with its node
# names already on the dimensions.
P <- matrix(c(0.7, 0.2, 0, 0.3, 0.8, 0.5), nrow = 3,
  dimnames = list(c("n1", "n2", "n3"), c("A", "B")))
transport <- external_transport(
  P, semantics = "density",
  provenance = list(details = "partial-volume warp from atlas-tool 2.1")
)
transport
#> <effect_location_transport>
#>   nodes:      3 native -> 2 group + sink
#>   semantics:  density
#>   row mass:   unit (one per native node)
#>   sink:       mass 0.5 of 3 rows, 16.7% of territory
#>   provenance: external (cross-fit: none)
#>   built:      partial-volume warp from atlas-tool 2.1
#>   signature:  sha256:4286fd2f411f...

# Density divides the transported budget by the transported row mass; with
# the default unit row mass that is the transported node count.
transport_values(transport, c(1, 1, 1))
#>      A      B <sink> 
#>    1.0    1.0    0.5 
```
