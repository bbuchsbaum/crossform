# Split a neural metric into coherent and configuration components

For SPD `K` and an oriented covector `a`, the exact metric-aware split
is `K_coh = a a^T / (a^T K^-1 a)` and `K_cfg = K - K_coh`. The coherent
amplitude remains `B a`; the metric controls only its normalization.

## Usage

``` r
metric_components(metric, coherent = NULL)
```

## Arguments

- metric:

  A materialized
  [`neural_metric()`](https://bbuchsbaum.github.io/crossform/reference/neural_metric.md)
  or an internal composed frame-metric value.

- coherent:

  A matching
  [`coherent_functional()`](https://bbuchsbaum.github.io/crossform/reference/coherent_functional.md)
  or numeric vector. It may be omitted for a composed frame metric.

## Value

An `effect_metric_components` value.

## Examples

``` r
domain <- abstract_domain(3, id = "component-example")
metric <- neural_metric(diag(c(1, 2, 3)), domain)
mean_functional <- crossform:::coherent_functional(rep(1 / 3, 3), domain)
components <- crossform:::metric_components(metric, mean_functional)
all.equal(components$coherent + components$configuration, metric$value)
#> [1] TRUE
```
