# `crossform` source layering

Status: enforced by `tests/testthat/test-architecture.R`
Date: 2026-08-16
Companion to: [crossform-package-design.md](crossform-package-design.md)

## Why this document exists

An architecture audit of the 2026-08-15 baseline found 24 of 61 source files in
one strongly connected component, and 91% of cross-file calls going to
dot-prefixed internals. The immediate cause was four calls from the compute
core up into the presentation layer: `R/kernel.R` reached into `R/views.R` for
`.unsvec_symmetric()` and `.align_contrast()`, and into `R/result.R` for
`.svec_symmetric()`. Each of those is a pure leaf function with no business
being in a view file; together they made the kernel undiscussable without the
views and vice versa.

This document states the layering the package now follows, so that the next
such call is a deliberate decision rather than an accident of where a helper
happened to be typed.

## The rule

Six layers. **A file may call downward or sideways. It may never call upward.**

| # | Layer | Files |
|---|-------|-------|
| 1 | **primitives** | `primitives.R`, `message-helpers.R`, `conditions.R`, `check.R`, `RcppExports.R` |
| 2 | **values** | `domain.R`, `frame.R`, `pairing.R`, `relation.R`, `effect-space.R`, `effect-map.R`, `metric.R`, `metric-learning.R`, `source.R`, `capabilities.R`, `study.R`, `study-facts.R`, `design-model.R`, `observation-model.R`, `extractor.R`, `scope.R`, `support-index.R`, `numerics.R`, `query-structured.R`, `pair-query.R`, `operations.R`, `receipt.R`, `reliability.R`, `validation-memo.R`, `measurement.R`, `measurement-storage.R`, `relation-fit.R`, `residual-statistics.R`, `bridge.R`, `compute-policy.R`, `memory-plan.R`, `crossform-package.R` |
| 3 | **plans** | `geometry-plan.R`, `relation-plan.R`, `crossnobis.R`, `evidence-task.R`, `evidence-sampling.R`, `evidence-sampling-kernel.R`, `evidence-sampling-product.R`, `compiler-conformance.R` |
| 4 | **compiler / execution, and the records execution produces** | `compiler.R`, `execution-driver.R`, `kernel.R`, `task.R`, `storage.R`, `measurement-kernel.R`, `crossnobis-driver.R`, `result.R` |
| 5 | **results / views** | `views.R`, `geometry-entry.R`, `coupling-views.R`, `tomography.R`, `measurement-result.R`, `measurement-decomposition.R`, `format-results.R`, `print-methods.R`, `plot-methods.R` |
| 6 | **adapters and facade** | `adapter-bids.R`, `adapter-fmridesign.R`, `adapter-fmrireg.R`, `neuroim2-adapter.R`, `bridge.R` consumers, `benchmark.R`, `example-data.R`, `evidence-api.R` |

The package deliberately has **no `Collate:` field**. R sources are loaded in
locale-collation order and every definition is resolved at call time, so file
order carries no meaning and adding a `Collate:` would only create a second,
silently divergent statement of dependency. The layering is a *design*
constraint, checked by a test; it is not a load-order constraint.

### Layer 1 in detail

`R/primitives.R` holds the leaves that more than one layer needs and that
nothing in the package should own:

- **Content addressing.** `.sha256_signature()`, `.sha256_string()`,
  `.sha256_file()`, `.is_sha256_signature()`, `.strong_sha256()`. Every
  identity in the package is a SHA-256 over an R value, and `serialize` /
  `serializeVersion` are pinned here in exactly one place. Before this pass the
  same `paste0("sha256:", digest::digest(x, algo = "sha256", serialize = TRUE,
  serializeVersion = 2L))` expression appeared at 102 sites; a change to either
  argument at any one of them would silently have invalidated recorded
  signatures without a single test noticing the inconsistency.
- **Tiling and store I/O.** `.tile_starts()`, `.validate_tile_size()`,
  `.tile_io()`, `.tile_zero_fill()`, `.effect_form_codec_format()`. The
  package has two block stores with genuinely different on-disk formats — the
  dense column-major geometry store in `R/storage.R` and the variable-block,
  manifest-backed measurement store in `R/measurement-storage.R` — and they
  are not interchangeable. What they share is the transfer primitive
  underneath both: a headerless little-endian float64 payload, opened once,
  seeked to an element offset, transferred as one run, closed once.
  `.tile_io()` is that primitive and the only place in the package that names
  the element size, the endianness, and the `offset * 8` arithmetic.
- **Symmetric packing.** `.svec_symmetric()`, `.unsvec_symmetric()`,
  `.physical_query_operator()`, `.physical_query_operators()`.
- **Contrast alignment.** `.align_contrast()`.

`R/conditions.R` holds the condition taxonomy: `.input_error()`,
`.contract_error()`, `.invariant_error()`, and the older
`.capability_refusal()`. The first three build conditions that inherit
`c("<class>", "effect_error", "error", "condition")` and carry `arg`,
`received`, and `expected`.

`R/check.R` holds the argument guards that raise them —
`.check_string()`, `.check_flag()`, `.check_count()`, `.check_number()`,
`.check_matrix()`, `.check_class()` — plus the two halves of the sealed-record
validator prologue, `.sealed_fields()` and `.check_signature()`. Before this
file the same five shape tests were hand-rolled at hundreds of sites, each
with its own phrasing and its own bare `stop()`.

Within layer 1 the direction is `primitives.R` → `check.R` →
`message-helpers.R` → `conditions.R`. `conditions.R` calls nothing in the
package at all; everything else in the package may call any of the four.

### Layer 4 in detail: the compiler is not the executor

Layer 4 is two jobs, and it used to be one file. `R/compiler.R` held task
construction, task identity, query lowering, memory planning, the exported
entry points, *and* a 324-line `.execute_geometry_plan()` with seven inline
closures. Because plans needed the task constructors and the compiler needed
`plan_geometry()`, the plan and compiler layers were mutually recursive and
sat inside one 26-file strongly connected component.

The direction is now **entry → plan → compiler → kernel**, and the four
responsibilities live apart:

- `R/geometry-entry.R` (layer 5) owns the two exported entry points,
  `materialize_geometry()` and `evaluate_geometry()`, including the relation
  compatibility form `f(relation, at, over)` that builds a plan first. Plan
  construction happens on the caller's side of the boundary.
- `R/geometry-plan.R` (layer 3) is pure data: the metric schedule, the plan's
  estimand identity, `plan_geometry()`, and the plan validator. It calls
  nothing in layer 4.
- `R/evidence-task.R` (layer 3) owns effect-task construction and identity —
  `.compile_effect_evidence_task()`, `.effect_task_semantic()`,
  `.effect_task_id()`, `.effect_task_base_id()`, `.validate_compiled_effect_task()`.
  Naming a task is plan work; it moved out of the compiler.
- `R/compiler.R` (layer 4) lowers a plan it is handed into an
  `effect_geometry_execution_plan`: the physical query, the component
  requirements, the tiles and memory plan, the kernel version, the execution
  identity. It constructs no plan and executes nothing, and it has **zero**
  upward edges.
- `R/execution-driver.R` (layer 4) is the runtime. The executor is a sequence
  of named stages — `.open_execution_sources()`, `.execute_node_block()`,
  `.execute_coherent_block()`, `.execute_marginal_block()`,
  `.execution_metadata()`, `.execution_geometry_result()` — driven by
  `.execute_geometry_compute()`, which is 30 lines at nesting depth 2 with no
  inline closures. `on.exit()` is registered in that one frame, so the source
  session still closes exactly once whether the run completes, errors, or is
  interrupted.

Three declaration files moved down to layer 2 at the same time, because that
is what they always were: `R/compute-policy.R` (was `execution.R`),
`R/memory-plan.R` (was `memory.R`), and `R/capabilities.R`, which holds the
relation-source admission check formerly called `.compiler_capabilities()` and
now named `.relation_source_capabilities()`.

Two files joined layer 4 in the pass that emptied the register:

- `R/crossnobis-driver.R` is the crossnobis runtime, split out of
  `R/crossnobis.R` on exactly the precedent above. `crossnobis.R` was a plan
  file that also ran its own plan, reaching up into `.run_geometry_compiler()`
  and `.support_streamed_scheduled_crossnobis()`. What is left in
  `crossnobis.R` is the plan — `noise_precision()`, `plan_crossnobis()`, the
  validator, the print method — and what moved is everything that executes:
  the planned receipt, the learned-metric runtime, and the exported
  `crossnobis()` entry that dispatches between the learned route and the
  ordinary geometry compiler.
- `R/result.R` holds the sealed result records — `effect_form`,
  `effect_geometry`, `effect_view`, `effect_contrast_view` — and their
  validators. It was classified as a view, which made the executor's
  construction of its own output an upward call. It is not a view: it calls
  nothing above layer 3, it renders nothing, and every reader of a result
  (`views.R`, `format-results.R`, `print-methods.R`, `plot-methods.R`) sits
  above it. It is the executor's output contract, so it sits with the
  executor. `.new_effect_contrast_view()` moved into it from `views.R` for the
  same reason; the derivation that feeds it stays in `contrast_energy()`.

The layering permits `kernel.R` to call `result.R` sideways now, but the test
still forbids it: a kernel that builds a result record has stopped being a
numerical primitive. Only the executor constructs results.

## What the test checks

`tests/testthat/test-architecture.R` parses `R/*.R`, builds the internal call
graph by matching every symbol against the map of top-level function
definitions, and asserts:

1. every source file has a declared layer (so a new file cannot be silently
   unclassified);
2. no upward edge exists that is not in the `allowed_upward` register;
3. no register entry is stale — an edge that has been removed must also be
   removed from the register, so the debt list can only shrink;
4. `kernel.R` and `task.R` never call into any results/views file;
5. `primitives.R` calls nothing above layer 1;
6. no values file reaches into the compiler or executor outside the register.

It runs in about half a second. It reads the package **sources**, not the
installed namespace, so it skips with an explicit reason when sources are not
on disk beside the test.

## The register: empty

`allowed_upward` is `character()`. Every upward edge the 2026-08-15 audit
found has been removed rather than tolerated, so there is no standing debt to
summarise here and any new entry is a new violation that has to be argued for
on its own terms.

The last nine went as follows.

**`crossnobis.R` was a plan and its own executor (2 edges).** Split, as its
follow-up asked: the runtime is now `R/crossnobis-driver.R` in layer 4 and the
plan file executes nothing. See *Layer 4 in detail* above.

**A source session validated the task it was handed (1 edge).** Inverted.
`.open_effect_task_source_session()` no longer calls
`.validate_compiled_effect_task()`; both package callers already passed
`validate = FALSE`, so nothing changed at runtime. The validator moved to the
point where the value is *built*, `.as_compiled_effect_task()` in
`evidence-task.R` — the plan layer validating its own record. The two files
also stopped duplicating the session: `source.R` and `measurement-kernel.R`
had 55 near-identical lines of two-sided open/close/read wiring, now one
`.open_two_sided_source_session()` parameterized by class, side noun, and
whether the relations are revalidated.

**The executor constructing its own results (3 edges).** Reclassified and
moved. `result.R` is layer 4 (it is the executor's output contract, not a
view); `.new_effect_contrast_view()` moved into it from `views.R`; and
`.effect_form_codec_format()`, which `storage.R` needed, was a leaf constant
and moved to `primitives.R` where both stores can reach it.

**`evidence-api.R` was two files wearing one name (3 edges).** Split three
ways. `edge_frame()` and its validator moved down to `measurement.R`, the
values file that owns everything they are built from. `coupling()` moved *in*:
it was `R/coupling-plan.R`, classified as a plan, but it takes a compiled plan
and returns a completed measurement form — a public entry that plans and
executes, exactly like the `measurement_form()` it wraps — and being
misclassified was the sole cause of three of these edges. `coupling-plan.R` is
gone. `measurement_frame()` stayed in the facade, deliberately: it adapts an
additive frame through `.measurement_frame_from_additive_decomposition()` in
the layer-5 `measurement-decomposition.R`, so moving it to a values file would
have opened a worse edge than it closed. The
`measurement-decomposition.R -> evidence-api.R` edge was not a call at all —
a local variable named `geometry_alignment` shadowed the exported view of the
same name — and was fixed by renaming the local.

### What left the register earlier

Twelve of the original twenty-one edges went in the previous pass, and two
changed target because the code they named moved:

| Removed edge | How |
|---|---|
| `metric.R -> geometry-plan.R` | `.validate_geometry_metric_schedule()` moved down to `metric.R` |
| `receipt.R`, `coupling-plan.R`, `crossnobis.R`, `geometry-plan.R` `-> execution.R` | `execution.R` became the layer-2 `compute-policy.R` |
| `receipt.R`, `residual-statistics.R`, `crossnobis.R` `-> memory.R` | `memory.R` became the layer-2 `memory-plan.R` |
| `relation.R`, `relation-fit.R` `-> compiler.R` | `.compiler_capabilities()` became `.relation_source_capabilities()` in the layer-2 `capabilities.R` |
| `evidence-task.R -> compiler.R` | task construction and identity moved into `evidence-task.R` |
| `geometry-plan.R -> compiler.R` | the plan stopped calling the compiler; execution-plan lowering moved into `compiler.R` |

| Retargeted edge | Why |
|---|---|
| `source.R -> compiler.R` became `source.R -> evidence-task.R` | the validator moved beside the task value, as its follow-up asked |
| `crossnobis.R -> compiler.R` became `crossnobis.R -> execution-driver.R` | the runner crossnobis calls moved out of the compiler |
| `compiler.R -> result.R`/`views.R` became `execution-driver.R -> result.R`/`views.R` | result construction moved with the executor |

## What remains

No file has an upward edge. The largest strongly connected component is down
from the 26 files it held before the layer-4 split, to 17 after it, to **15**
now: the compiler, both executors, the kernels, the plans, the results, and
the views are all outside it.

What is left is one genuine tangle, in layer 2 and entirely within it:

```
capabilities.R  effect-map.R   effect-space.R  extractor.R  measurement.R
memory-plan.R   metric.R       operations.R    pairing.R    receipt.R
relation-fit.R  relation.R     scope.R         source.R     study-facts.R
```

These fifteen value files are mutually recursive. That is a different problem
from the one this document was written about — it is not a layering violation,
because a file may call sideways, and there is no direction to restore, only a
knot to untie. Untangling it means deciding which of `relation`, `source`,
`metric`, and `receipt` is the more primitive value, and that is a design
question about the vocabulary rather than a refactor. It is the honest next
piece of architecture work, and no pass so far has attempted it.
