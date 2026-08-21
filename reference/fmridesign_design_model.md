# Compile a semantic design model from fmridesign

This exact-version adapter binds a compiled
[`fmridesign::event_model()`](https://bbuchsbaum.github.io/fmridesign/reference/event_model.html)
to the event facts and observation rows of a study. Scientific condition
effects remain declared in condition space; concrete design coding is a
compilation receipt.

## Usage

``` r
fmridesign_design_model(
  model,
  study,
  basis_id,
  units,
  scale = 1,
  block_map = NULL,
  semantic_map = NULL,
  specification = NULL,
  solver = "auto"
)
```

## Arguments

- model:

  A `fmridesign` event model from the certified package version.

- study:

  A bound
  [`study()`](https://bbuchsbaum.github.io/crossform/reference/study.md)
  containing the exact timed event facts.

- basis_id, units, scale:

  Scientific coordinate declarations passed to
  [`condition_space()`](https://bbuchsbaum.github.io/crossform/reference/condition_space.md).

- block_map:

  Optional named integer map from study partitions to fmridesign blocks.

- semantic_map:

  Optional condition-by-coefficient matrix, or named list per partition,
  for models that cannot be lowered from column metadata.

- specification:

  Optional portable semantic declaration. When omitted, the adapter
  derives one from the pinned fmridesign model.

- solver:

  Numerical route recorded in the design receipt.

## Value

An
[`design_model()`](https://bbuchsbaum.github.io/crossform/reference/design_model.md)
suitable for
[`plan_relation()`](https://bbuchsbaum.github.io/crossform/reference/plan_relation.md),
whose `$condition_space` holds the fmridesign condition names, whose
`$designs` hold one compiled matrix per study partition with study
observation ids as row names, and whose `$compiler` records the pinned
fmridesign version.

## See also

[`design_model()`](https://bbuchsbaum.github.io/crossform/reference/design_model.md)
for the generic contract,
[`plan_relation()`](https://bbuchsbaum.github.io/crossform/reference/plan_relation.md)
for the next step, and
[`compiler_conformance()`](https://bbuchsbaum.github.io/crossform/reference/compiler_conformance.md)
to check the receipt this adapter produced.

Other relation planning and fitting:
[`adapter_version_certificate()`](https://bbuchsbaum.github.io/crossform/reference/adapter_version_certificate.md),
[`compiler_conformance()`](https://bbuchsbaum.github.io/crossform/reference/compiler_conformance.md),
[`effect_extractor()`](https://bbuchsbaum.github.io/crossform/reference/effect_extractor.md),
[`estimate_relation()`](https://bbuchsbaum.github.io/crossform/reference/estimate_relation.md),
[`file_matrix_source()`](https://bbuchsbaum.github.io/crossform/reference/file_matrix_source.md),
[`fmrireg_relation()`](https://bbuchsbaum.github.io/crossform/reference/fmrireg_relation.md),
[`lm_extractor()`](https://bbuchsbaum.github.io/crossform/reference/lm_extractor.md),
[`lm_relation_fit()`](https://bbuchsbaum.github.io/crossform/reference/lm_relation_fit.md),
[`plan_relation()`](https://bbuchsbaum.github.io/crossform/reference/plan_relation.md),
[`relation()`](https://bbuchsbaum.github.io/crossform/reference/relation.md),
[`relation_block()`](https://bbuchsbaum.github.io/crossform/reference/relation_block.md),
[`relation_fit()`](https://bbuchsbaum.github.io/crossform/reference/relation_fit.md),
[`relation_fit_capabilities()`](https://bbuchsbaum.github.io/crossform/reference/relation_fit_capabilities.md),
[`relation_plan_receipts()`](https://bbuchsbaum.github.io/crossform/reference/relation_plan_receipts.md),
[`residual_block()`](https://bbuchsbaum.github.io/crossform/reference/residual_block.md),
[`residual_df()`](https://bbuchsbaum.github.io/crossform/reference/residual_df.md),
[`source_capabilities()`](https://bbuchsbaum.github.io/crossform/reference/source_capabilities.md)

## Examples

``` r
# This adapter is certified against exactly one fmridesign version and
# refuses any other, so the example runs only under that version.
if (requireNamespace("fmridesign", quietly = TRUE) &&
    identical(as.character(utils::packageVersion("fmridesign")), "0.6.0")) {
  set.seed(1)
  domain <- abstract_domain(3L, id = "fmridesign-example")
  index <- observation_index(
    1:20, "run-1", time = seq(0, by = 2, length.out = 20L), units = "seconds"
  )
  events <- data.frame(
    partition = "run-1", event_id = paste0("e", 1:4),
    onset = c(2, 10, 18, 26), duration = 0,
    condition = c("face", "body", "face", "body")
  )
  facts <- study(
    observations(list(`run-1` = matrix(rnorm(60), 20L, 3L)),
      list(`run-1` = index), domain),
    observation_events(events)
  )

  # The event model must be the one compiled from these exact facts; the
  # adapter re-checks onsets, durations, and condition labels against them.
  external <- fmridesign::event_model(
    onset ~ fmridesign::hrf(condition), data = events,
    block = rep(1L, nrow(events)),
    sampling_frame = fmridesign::sampling_frame(
      blocklens = 20L, TR = 2, start_time = 0
    ),
    durations = events$duration
  )
  model <- fmridesign_design_model(
    external, facts, basis_id = "canonical-hrf-amplitude",
    units = "percent-signal-change"
  )
  print(model$condition_space$coordinates)
  print(model$capabilities$coding_invariant)

  # Effects are now named in condition space, so the contrast below does not
  # depend on which columns the HRF compiler emitted.
  weights <- rbind(`face-body` = c(1, -1))
  colnames(weights) <- model$condition_space$coordinates
  plan <- plan_relation(
    facts, model, effect_map(weights, model$condition_space),
    observation_model("ols", sampling_unit = "scan")
  )
  print(all(as.matrix(compiler_conformance(plan)[-1L])))
}
#> [1] "condition.face" "condition.body"
#> [1] TRUE
#> [1] TRUE
```
