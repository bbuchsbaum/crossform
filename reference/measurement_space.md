# Name a common measurement space

`measurement_space()` names the shared coordinate system that both legs
of a
[`measurement_bridge()`](https://bbuchsbaum.github.io/crossform/reference/measurement_bridge.md)
map into, so two different neural domains can be related without ever
forming a dense cross-domain operator. Use it before building a bridge
between, for example, two participants' native feature spaces.

## Usage

``` r
measurement_space(n_measurements, id, provenance = list())
```

## Arguments

- n_measurements:

  Positive common-coordinate count.

- id:

  Stable nonempty identity.

- provenance:

  Portable fixed-space provenance.

## Value

An `effect_measurement_space` carrying `$id`, `$n_measurements`,
`$provenance`, and a content-addressed `$signature` that binds the
identity into every bridge built on it.

## See also

[`measurement_bridge()`](https://bbuchsbaum.github.io/crossform/reference/measurement_bridge.md),
which requires both legs to have this many rows.

## Examples

``` r
# A named 2-coordinate common space shared by two native neural domains.
common <- crossform:::measurement_space(
  2, id = "study:common-modes:v1",
  provenance = list(source = "group template modes")
)
c(id = common$id, n = common$n_measurements)
#>                      id                       n 
#> "study:common-modes:v1"                     "2" 

# The signature is content-addressed: the same declaration reproduces it,
# so a stored bridge can be checked against the space it claims.
identical(
  common$signature,
  crossform:::measurement_space(
    2, id = "study:common-modes:v1",
    provenance = list(source = "group template modes")
  )$signature
)
#> [1] TRUE
```
