# `effect_coupling_result`: one class, seven kinds

Status: decided and implemented
Date: 2026-08-17
Ticket: A7
Scope: `R/coupling-views.R`, `man/effect_coupling_result.Rd`,
`tests/testthat/test-coupling-views.R`

## Decision

**KEEP one class.** `effect_coupling_result` continues to carry all seven
readings, discriminated by `$kind`. The enumeration is now closed, enforced at
construction and at validation, documented as a class contract, and tested by
iterating over every kind.

This is the opposite call to the one the elite audit made about the
label-subclass zoo elsewhere in the package, and it is the same principle
applied to the opposite evidence. The rule the maintainer stated is that *a
`kind` string is fine only where readers are shared*. Here they are shared:
not one reader of this class branches on `kind`. The table below is the
evidence.

## Evidence

Reader inventory for `effect_coupling_result`, complete as of this date:

| Reader | Location | `switch(x$kind` / `if (kind ==` branches |
|---|---|---|
| `format.effect_coupling_result` | `R/print-methods.R:1377` | **0** — prints `x$kind` as a recorded string, the way it prints `$terminology` |
| `print.effect_coupling_result` | `R/print-methods.R:1383` | **0** |
| `.validate_coupling_result` | `R/coupling-views.R` | **0** |
| `as.data.frame` | — | no method exists; all seven kinds fail identically |
| `plot` | — | no method exists; all seven kinds fail identically |
| accessors | — | none; every field is read directly off the list |

Per-kind detail. All seven carry the identical fourteen sealed fields in the
identical order (`kind`, `values`, `edge_index`, `source_plan`,
`source_receipt`, `normalization_axis`, `summary_axis`, `stage_order`,
`regularization`, `units`, `terminology`, `partition_policy`,
`edge_completeness`, `signature`) under the identical signature scheme, so
`.sealed_fields()` is one exact-names check for the whole class:

| `$kind` | `$values` shape | columns / element | `$regularization` | `$partition_policy` | `$units` |
|---|---|---|---|---|---|
| `effect_coupling` | `edge_blocks` | list, one matrix per edge | `NULL` | `NULL` | `NULL` |
| `covariance_coupling` | `edge_blocks` | list, one matrix per edge | `NULL` | `NULL` | `NULL` |
| `pearson_correlation` | `edge_table` | `edge_id`, `correlation` | `NULL` | `NULL` | `correlation` |
| `partitioned_pearson_coupling` | `edge_table` | `edge_id`, `value`, `transform` | `NULL` | **present** | `correlation` or `fisher_z` |
| `canonical_coupling` | `edge_table` | `edge_id`, `mode`, `canonical_correlation` | **present** | `NULL` | `canonical_correlation` |
| `geometry_alignment` | `edge_table` | `edge_id`, `geometry_alignment` | `NULL` | `NULL` | `linear_cka` |
| `gaussian_mutual_information` | `edge_table` | `edge_id`, `information`, `units` | **present** | `NULL` | `nats` or `bits` |

Note the shape column: the only structural variation across seven kinds is
**two-valued**. Splitting into per-kind subclasses would mint seven types to
express a two-way distinction, and would then need a parent for the shared
sealed-field set, the shared signature, the shared validator, and the two
shared readers — i.e. everything except one boolean.

`$units` varies *within* a kind (`partitioned_pearson_coupling` reports
`correlation` or `fisher_z` depending on the edge transform), so it is
recorded, not constrained.

## Contrast: `R/tomography.R` is already split, correctly

`tomography.R` passes a `kind` to `.tomography_signature()`, but that is a
signature-domain salt across **three separate classes**
(`effect_tomography_resource_plan`, `effect_measured_block_form`,
`effect_tomography_result`), each with its own sealed field set and its own
validator. It is never read back as a discriminator. Those three genuinely
differ — different fields, different validators, and only the last has
readers — so they are three classes. Nothing to change there; the two files
are consistent applications of the same rule, not an inconsistency.

## What was implemented

1. `.coupling_kinds` in `R/coupling-views.R`: the closed registry mapping each
   kind to its value shape, its exact `names(values)` when it is a table, and
   whether `regularization` / `partition_policy` must be present.
2. `.coupling_kind_contract()` raises `effect_input_error` on any unlisted
   kind. Previously the constructor accepted any string, so a typo'd or
   invented kind produced a well-formed, signed, unreadable record.
3. `.check_coupling_kind_contract()` enforces the shape, the exact column
   names, and the presence of the fields that record what changed the numbers.
   It runs in `.new_coupling_result()` and again in
   `.validate_coupling_result()`.
4. `.coupling_value_shape()`: the one branch a reader of `$values` needs.
5. `man/effect_coupling_result.Rd` (roxygen in `R/coupling-views.R`) states
   the closed enumeration and what each kind guarantees, under
   `@section Kinds:`.
6. Three tests in `tests/testthat/test-coupling-views.R`: the enumeration is
   closed *and* exhaustive (registry set == reachable set, so a new view
   without a registry entry fails), every kind satisfies the shared reader
   contract, and mismatched shapes/columns/optional fields are refused.

## Two print defects found while testing the shared-reader claim

Both were in `R/print-methods.R`, which this ticket did not own; both were
handed to that file's owner and are **now fixed**. Neither is a
counter-argument to KEEP — neither was a *per-kind* branch; both were one
shared reader mis-reading kind-dependent content, and both were one-line
fixes in one place, which is exactly the cost profile a single class should
have. Iterating print over all seven kinds is what surfaced them: two of the
seven had never been printed anywhere in the suite, because the only prior
print test used `effect_coupling()`, the one kind that triggers neither.

1. `values = paste0(length(x$values), " blocks")`. For the five `edge_table`
   kinds `length()` is the **column count**, so a `canonical_coupling` over 4
   edges with 8 rows printed `values: 3 blocks`. Correct for the two
   `edge_blocks` kinds only. Fixed by branching on `.coupling_value_shape(x)`
   — a sideways call, layer 5 to layer 5, permitted — to print
   `n rows x m columns` for table kinds.

2. `.pf_num(unlist(x$regularization, use.names = FALSE))`. The
   `effect_measurement_regularization` record holds `kind`, `applied`, and
   `signature` alongside the two lambdas, so `unlist()` coerced the whole
   record to character and `.pf_num()`'s `as.numeric()` yielded
   `regularization: NA, 0.05, 0.02, NA, NA` plus two
   "NAs introduced by coercion" warnings. Affected `canonical_coupling` and
   `gaussian_mutual_information`. Fixed by a `.pf_regularization()` helper
   reading `lambda_left` / `lambda_right` by name and keeping the ridge kind
   visible (`ridge 0.05`, or `ridge (left 0.05, right 0.02)` when the sides
   differ). `.pf_num()` itself was additionally hardened to fall back to a
   text rendering rather than a row of `NA`, since it has ~30 call sites and
   a printer should never warn.

Both `unlist()`-on-a-record and `length()`-on-a-data-frame are shape errors
that a type would not have caught either: a per-kind subclass would have had
the same two bugs, in seven places instead of one.

## Deferred, deliberately

No `as.data.frame` or `plot` method was added. Neither exists today, so
neither is a reader that differs — all seven kinds behave identically under
both. Adding one is new public API surface (NAMESPACE, Rd, and a decision
about what a block kind should flatten to) and belongs to its own ticket, not
to a decision about whether the existing class is honest. Coupling and
tomography remain an experimental, small-node tier.
