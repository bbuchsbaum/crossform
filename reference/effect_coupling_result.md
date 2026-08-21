# Coupling results and the readings they may carry

Every view over a completed `effect_measurement_form` returns one
`effect_coupling_result`. There are no subclasses: the seven readings
below are distinguished by the `$kind` string, which is a **closed**
enumeration. A result of any kind carries the same fourteen fields, is
validated the same way, and prints the same way, so code that reads one
kind reads them all; only the shape of `$values` differs, and it differs
two ways rather than seven.

## Details

Every kind carries `$values`, `$edge_index` naming the edges those
values are keyed to, `$source_plan` and `$source_receipt` identifying
the forms it was read from, `$normalization_axis` and `$summary_axis`
stating what was divided out and what was summarized over,
`$stage_order` listing the pipeline that produced it, `$units` and
`$terminology` stating what the numbers may be called,
`$edge_completeness`, and a `$signature` over all of it.
`$regularization` and `$partition_policy` are present exactly when the
kind applied one; they are `NULL` otherwise, never silently defaulted.

## Kinds

`$kind` is one of the following seven values and no others. Constructing
a result with an unlisted kind, or with a `$values` shape that disagrees
with its kind, is a contract error rather than an accepted record.

Two kinds report **blocks**: `$values` is a named list holding one
numeric matrix per edge, in `$edge_index$edge_id` order.

- `"effect_coupling"` – the raw measurement block, with no covariance
  claim attached. `$normalization_axis` is `"none"` and `$units` is
  `NULL`, because nothing has been divided out and nothing is claimed
  about what the numbers mean.

- `"covariance_coupling"` – the same blocks, now certified as
  repeated-sample covariance: the form established repeated variation of
  effective rank above one and symmetric positive self-blocks.

Five kinds report a **table**: `$values` is a data frame keyed by
`edge_id` whose remaining columns are fixed by the kind.

- `"pearson_correlation"` – columns `edge_id`, `correlation`. One signed
  scalar per edge, requiring rank-one measurement axes.

- `"partitioned_pearson_coupling"` – columns `edge_id`, `value`,
  `transform`. A weighted reduction across several source forms; carries
  a `$partition_policy` recording the weights, the placement, and the
  edge transform, and `$units` is `"correlation"` or `"fisher_z"`
  according to that transform.

- `"canonical_coupling"` – columns `edge_id`, `mode`,
  `canonical_correlation`. Several rows per edge, one per canonical mode
  in descending order; carries the `$regularization` whose ridge changed
  the values.

- `"geometry_alignment"` – columns `edge_id`, `geometry_alignment`.
  Static linear CKA/RV-like alignment, normalized over form entries
  rather than over experimental samples.

- `"gaussian_mutual_information"` – columns `edge_id`, `information`,
  `units`. Carries the `$regularization` used for the underlying
  canonical spectrum, and a `$terminology` naming the signature of the
  [`gaussian_covariance_model()`](https://bbuchsbaum.github.io/crossform/reference/gaussian_covariance_model.md)
  declaration it rests on.

## See also

[`effect_coupling()`](https://bbuchsbaum.github.io/crossform/reference/coupling_views.md),
[`covariance_coupling()`](https://bbuchsbaum.github.io/crossform/reference/coupling_views.md),
[`canonical_coupling()`](https://bbuchsbaum.github.io/crossform/reference/coupling_views.md),
[`geometry_alignment()`](https://bbuchsbaum.github.io/crossform/reference/coupling_views.md),
and
[`connectivity()`](https://bbuchsbaum.github.io/crossform/reference/connectivity.md),
which produce these results, and
[`gaussian_covariance_model()`](https://bbuchsbaum.github.io/crossform/reference/gaussian_covariance_model.md)
for the declaration Gaussian information requires.

Other coupling and connectivity views:
[`connectivity()`](https://bbuchsbaum.github.io/crossform/reference/connectivity.md),
[`control_coupling()`](https://bbuchsbaum.github.io/crossform/reference/control_coupling.md),
[`coupling()`](https://bbuchsbaum.github.io/crossform/reference/coupling.md),
[`coupling_contrast()`](https://bbuchsbaum.github.io/crossform/reference/coupling_contrast.md),
[`coupling_views`](https://bbuchsbaum.github.io/crossform/reference/coupling_views.md),
[`match_control()`](https://bbuchsbaum.github.io/crossform/reference/match_control.md),
[`match_coupling()`](https://bbuchsbaum.github.io/crossform/reference/match_coupling.md),
[`measurement_components()`](https://bbuchsbaum.github.io/crossform/reference/measurement_components.md),
[`measurement_form()`](https://bbuchsbaum.github.io/crossform/reference/measurement_form.md),
[`pair_lm_query()`](https://bbuchsbaum.github.io/crossform/reference/pair_lm_query.md),
[`pair_query()`](https://bbuchsbaum.github.io/crossform/reference/pair_query.md)
