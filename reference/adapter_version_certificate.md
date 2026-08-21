# Certify the installed version of an adapter's upstream package

[`vignette("crossform-extending")`](https://bbuchsbaum.github.io/crossform/articles/crossform-extending.md)
places two obligations on an adapter author, and this is the first of
them: certify against the *installed* version of the package you adapt,
and record what you certified against. The obligation is part of the
protocol, so the refusal that discharges it is too — an adapter outside
this package raises the same two refusals, in the same namespace, as the
two adapters shipped here.

## Usage

``` r
adapter_version_certificate(package, supported)
```

## Arguments

- package:

  Name of the upstream package the adapter is certified against.

- supported:

  The one version string the adapter has been tested with.

## Value

The installed version of `package`, as a character string equal to
`supported`. Anything else refuses: `installed_compiler_adapter` when
`package` is not installed, and `supported_compiler_version` when the
installed version is not the certified one. Both refusals carry the
`relation_compiler` namespace and name the certified version in their
remedies.

## Details

Call it first, before touching the upstream package, and record the
returned string in the provenance of whatever you build. A version other
than the certified one is refused rather than attempted: an untested
upstream release that still runs is the failure mode this exists to
prevent.

## See also

[`compiler_conformance()`](https://bbuchsbaum.github.io/crossform/reference/compiler_conformance.md),
which checks the receipt an adapter produced rather than the version
that produced it, and
[`catch_refusal()`](https://bbuchsbaum.github.io/crossform/reference/catch_refusal.md)
to inspect either refusal.

Other relation planning and fitting:
[`compiler_conformance()`](https://bbuchsbaum.github.io/crossform/reference/compiler_conformance.md),
[`effect_extractor()`](https://bbuchsbaum.github.io/crossform/reference/effect_extractor.md),
[`estimate_relation()`](https://bbuchsbaum.github.io/crossform/reference/estimate_relation.md),
[`file_matrix_source()`](https://bbuchsbaum.github.io/crossform/reference/file_matrix_source.md),
[`fmridesign_design_model()`](https://bbuchsbaum.github.io/crossform/reference/fmridesign_design_model.md),
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
# crossform certifies its own fmridesign adapter against exactly one
# version; this is the call that enforces it.
if (requireNamespace("fmridesign", quietly = TRUE)) {
  installed <- as.character(utils::packageVersion("fmridesign"))
  print(identical(
    adapter_version_certificate("fmridesign", installed), installed
  ))
}
#> Registered S3 method overwritten by 'fmridesign':
#>   method               from   
#>   print.sampling_frame fmrihrf
#> [1] TRUE

# An uncertified version is refused, not attempted. The refusal names the
# missing capability, so a caller can branch on the cause.
refusal <- catch_refusal(
  adapter_version_certificate("stats", "0.0.0-never-released")
)
refusal$capability
#> [1] "supported_compiler_version"
refusal$remedies
#> [1] "Use `stats` 0.0.0-never-released."                                      
#> [2] "Add and certify a version-specific adapter before granting conformance."

# A package that is not installed at all refuses differently.
catch_refusal(
  adapter_version_certificate("crossformNotAPackage", "1.0.0")
)$capability
#> [1] "installed_compiler_adapter"
```
