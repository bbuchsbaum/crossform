# Describe a read-only column-major matrix file

The file must contain only consecutive IEEE-754 doubles in R
column-major order. The constructor records an absolute path and
verifies a strong content revision. It never opens a persistent handle.

## Usage

``` r
file_matrix_source(
  path,
  dim,
  offset_bytes = 0,
  endian = .Platform$endian,
  stable_revision = NULL
)
```

## Arguments

- path:

  Existing binary matrix file.

- dim:

  Two positive integers: observations by neural features.

- offset_bytes:

  Nonnegative byte offset before the matrix payload.

- endian:

  Byte order used by the file.

- stable_revision:

  Optional expected `sha256:` content revision. When omitted it is
  computed from the file.

## Value

An `effect_source_descriptor`: a list with `$kind` (`"file_matrix"`),
the integer `$dim`, `$access` (`"reopenable"`), the `$stable_revision`
content hash, and a `$spec` holding the absolute `path`, `offset_bytes`,
and `endian`. Treat it as immutable.

## See also

[`source_capabilities()`](https://bbuchsbaum.github.io/crossform/reference/source_capabilities.md)
for the capability value it implies, and
[`observations()`](https://bbuchsbaum.github.io/crossform/reference/observations.md)
or
[`relation()`](https://bbuchsbaum.github.io/crossform/reference/relation.md),
which accept descriptors in place of in-memory matrices.

Other relation planning and fitting:
[`adapter_version_certificate()`](https://bbuchsbaum.github.io/crossform/reference/adapter_version_certificate.md),
[`compiler_conformance()`](https://bbuchsbaum.github.io/crossform/reference/compiler_conformance.md),
[`effect_extractor()`](https://bbuchsbaum.github.io/crossform/reference/effect_extractor.md),
[`estimate_relation()`](https://bbuchsbaum.github.io/crossform/reference/estimate_relation.md),
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
# A 4-observation by 3-feature matrix written as raw column-major doubles.
path <- tempfile(fileext = ".bin")
writeBin(as.vector(matrix(as.double(1:12), 4L, 3L)), path, size = 8L)

descriptor <- file_matrix_source(path, dim = c(4L, 3L))
descriptor$dim
#> [1] 4 3
descriptor$access
#> [1] "reopenable"

# The content hash is recorded now and rechecked whenever the file is
# reopened, so a silently edited source is caught rather than used.
substr(descriptor$stable_revision, 1, 24)
#> [1] "sha256:defe2a44f5b1a83ad"

# Declaring a revision that no longer matches the bytes is an error.
try(file_matrix_source(
  path, dim = c(4L, 3L), stable_revision = paste0("sha256:", strrep("0", 64))
))
#> Error : The matrix file does not match `stable_revision`.
unlink(path)
```
