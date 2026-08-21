# Construct a fixed factorized bridge between neural spaces

The induced cross-space operator is `t(left_leg) %*% right_leg`. Only
the factorized fixed legs are public; an arbitrary dense cross operator
is not.

## Usage

``` r
measurement_bridge(
  left_leg,
  right_leg,
  left_domain,
  right_domain,
  common_space,
  provenance = list(),
  reverse = FALSE
)

reverse_bridge(bridge)
```

## Arguments

- left_leg, right_leg:

  Finite common-coordinate-by-neural-feature matrices.

- left_domain, right_domain:

  Exact neural domain values or references.

- common_space:

  A named
  [`measurement_space()`](https://bbuchsbaum.github.io/crossform/reference/measurement_space.md)
  shared by both legs.

- provenance:

  Portable provenance for the fixed legs.

- reverse:

  One flag. `TRUE` exchanges the two sides of the declaration before the
  bridge is built.

- bridge:

  For `reverse_bridge()`, a factorized bridge whose two sides are to be
  exchanged.

## Value

An `effect_measurement_bridge` holding `$left_leg`, `$right_leg`, the
two `$left_domain`/`$right_domain` references, the shared
`$common_space`, and a `$signature` binding all of them.

## Details

`reverse = TRUE` declares the same two legs with their sides exchanged,
so a bridge built left-to-right can be reused when the same two domains
appear in the opposite order; the induced operator is the transpose of
the forward one, and the common space and provenance are preserved
exactly. `reverse_bridge()` is the internal shorthand that applies that
exchange to an already-built bridge.

## See also

[`measurement_space()`](https://bbuchsbaum.github.io/crossform/reference/measurement_space.md)
for the shared coordinates.

## Examples

``` r
# Two participants measured in different native feature spaces, related
# through two shared modes rather than a dense 3-by-4 operator.
left_domain <- abstract_domain(3, id = "subject01:native")
right_domain <- abstract_domain(4, id = "subject02:native")
common <- crossform:::measurement_space(2, id = "study:common-modes:v1")
bridge <- crossform:::measurement_bridge(
  left_leg = rbind(c(1, 0, 0), c(0, 1, 0)),
  right_leg = rbind(c(1, 0, 0, 0), c(0, 0, 1, 0)),
  left_domain = left_domain, right_domain = right_domain,
  common_space = common
)
dim(bridge$left_leg)
#> [1] 2 3

# The induced cross-space operator is never stored; it is exactly this.
t(bridge$left_leg) %*% bridge$right_leg
#>      [,1] [,2] [,3] [,4]
#> [1,]    1    0    0    0
#> [2,]    0    0    1    0
#> [3,]    0    0    0    0

# Declaring the same legs reversed transposes that operator exactly, and
# `reverse_bridge()` reaches the same value from the built bridge.
reversed <- crossform:::measurement_bridge(
  bridge$left_leg, bridge$right_leg, left_domain, right_domain, common,
  reverse = TRUE
)
identical(reversed, crossform:::reverse_bridge(bridge))
#> [1] TRUE
all.equal(
  t(bridge$left_leg) %*% bridge$right_leg,
  t(t(reversed$left_leg) %*% reversed$right_leg)
)
#> [1] TRUE

# Leg widths must match their declared domains, so a mismatch is caught
# before any neural value is read.
mismatch <- try(
  crossform:::measurement_bridge(
    rbind(c(1, 0, 0), c(0, 1, 0)), rbind(c(1, 0, 0), c(0, 1, 0)),
    left_domain, right_domain, common
  ),
  silent = TRUE
)
conditionMessage(attr(mismatch, "condition"))
#> [1] "Bridge leg widths must match their exact neural domains."
```
