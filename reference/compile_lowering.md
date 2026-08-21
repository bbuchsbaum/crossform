# Compile a frame-query pair to its algebraic lowering

Only an additive diagonal frame paired with a fixed bilinear query
admits the searchlight-collapse lowering. Factor frames, locally
estimated transforms, adaptive queries, and nonlinear readouts remain
distinct work.

## Usage

``` r
compile_lowering(frame, query)
```

## Arguments

- frame:

  An `effect_frame`.

- query:

  An `effect_query`.

## Value

A small compiler-decision value.
