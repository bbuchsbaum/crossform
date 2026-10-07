# Execute a relation plan with fmrireg

This installed-consumer adapter independently executes the OLS point
relation using
[`fmrireg::fmri_ols_fit()`](https://bbuchsbaum.github.io/fmrireg/reference/fmri_ols_fit.html).
It deliberately returns no residual error channel; analytic
second-moment uncertainty therefore refuses rather than borrowing
unsupported standard-error output.

## Usage

``` r
fmrireg_relation(x)
```

## Arguments

- x:

  A
  [`plan_relation()`](https://bbuchsbaum.github.io/crossform/reference/plan_relation.md)
  result using an OLS observation model.

## Value

An `effect_relation_fit` with point-relation capabilities only: its
`$error_models` are all `NULL`, so `residual_blocks`,
`effect_covariance`, and `residual_df` are `FALSE` for every partition.
Its `$provenance` records the `relation_plan_id`, the adapter, and the
pinned `adapter_version`.

## See also

[`estimate_relation()`](https://bbuchsbaum.github.io/crossform/reference/estimate_relation.md),
the in-package execution verb that does return a residual channel, and
[`relation_fit_capabilities()`](https://bbuchsbaum.github.io/crossform/reference/relation_fit_capabilities.md)
to see the difference.

Other relation planning and fitting:
[`adapter_version_certificate()`](https://bbuchsbaum.github.io/crossform/reference/adapter_version_certificate.md),
[`compiler_conformance()`](https://bbuchsbaum.github.io/crossform/reference/compiler_conformance.md),
[`effect_extractor()`](https://bbuchsbaum.github.io/crossform/reference/effect_extractor.md),
[`estimate_relation()`](https://bbuchsbaum.github.io/crossform/reference/estimate_relation.md),
[`file_matrix_source()`](https://bbuchsbaum.github.io/crossform/reference/file_matrix_source.md),
[`fmridesign_design_model()`](https://bbuchsbaum.github.io/crossform/reference/fmridesign_design_model.md),
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
# This adapter is certified against exactly one fmrireg version and refuses
# any other, so the example runs only under that version.
if (requireNamespace("fmrireg", quietly = TRUE) &&
    identical(as.character(utils::packageVersion("fmrireg")), "0.1.2")) {
  set.seed(1)
  domain <- abstract_domain(3L, id = "fmrireg-example")
  index <- observation_index(paste0("scan-", 1:4), "run-1")
  facts <- study(observations(
    list(`run-1` = matrix(rnorm(12), 4L, 3L)), list(`run-1` = index), domain
  ))
  design <- cbind(face = c(1, 0, 1, 0), body = c(0, 1, 0, 1))
  rownames(design) <- paste0("scan-", 1:4)
  target <- rbind(`face-body` = c(1, -1))
  colnames(target) <- colnames(design)
  plan <- plan_relation(
    facts, raw_design_model(list(`run-1` = design)), raw_effect_map(target),
    observation_model("ols", sampling_unit = "scan")
  )

  # An independent execution of the same point relation: use it to check
  # parity against the in-package estimator.
  external <- fmrireg_relation(plan)
  print(all.equal(
    relation_block(external, "run-1", 1:3),
    relation_block(estimate_relation(plan), "run-1", 1:3)
  ))

  # Parity is for point estimates only; no residual channel is claimed.
  print(relation_fit_capabilities(external)$residual_blocks)

  # A plan outside the certified OLS slice is refused, not approximated.
  gls <- plan_relation(
    facts, raw_design_model(list(`run-1` = design)), raw_effect_map(target),
    observation_model(
      "fixed_gls", sampling_unit = "scan", whitener = diag(4)
    )
  )
  print(catch_refusal(fmrireg_relation(gls))$capability)
}
```
