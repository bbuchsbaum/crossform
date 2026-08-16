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
| 1 | **primitives** | `primitives.R`, `message-helpers.R`, `conditions.R` |
| 2 | **values** | `domain.R`, `frame.R`, `pairing.R`, `relation.R`, `effect-space.R`, `effect-map.R`, `metric.R`, `metric-learning.R`, `source.R`, `capabilities.R`, `study.R`, `study-facts.R`, `design-model.R`, `observation-model.R`, `extractor.R`, `scope.R`, `support-index.R`, `numerics.R`, `query-structured.R`, `pair-query.R`, `operations.R`, `receipt.R`, `reliability.R`, `validation-memo.R`, `measurement.R`, `measurement-storage.R`, `relation-fit.R`, `residual-statistics.R`, `bridge.R`, `compute-policy.R`, `memory-plan.R`, `crossform-package.R` |
| 3 | **plans** | `geometry-plan.R`, `relation-plan.R`, `crossnobis.R`, `coupling-plan.R`, `evidence-task.R`, `evidence-sampling.R`, `evidence-sampling-kernel.R`, `evidence-sampling-product.R`, `compiler-conformance.R` |
| 4 | **compiler / execution** | `compiler.R`, `execution-driver.R`, `kernel.R`, `task.R`, `storage.R`, `measurement-kernel.R` |
| 5 | **results / views** | `result.R`, `views.R`, `geometry-entry.R`, `coupling-views.R`, `tomography.R`, `measurement-result.R`, `measurement-decomposition.R`, `format-results.R`, `print-methods.R`, `plot-methods.R` |
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
- **Tiling.** `.tile_starts()`, `.validate_tile_size()`.
- **Symmetric packing.** `.svec_symmetric()`, `.unsvec_symmetric()`,
  `.physical_query_operator()`, `.physical_query_operators()`.
- **Contrast alignment.** `.align_contrast()`.

`primitives.R` calls only `message-helpers.R`. `message-helpers.R` and
`conditions.R` call nothing in the package at all.

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

## The register: 9 upward edges that remain

Each is real debt with a named follow-up. They are listed in `allowed_upward`
in the test with the same comments; the summary is:

**`crossnobis.R` is a plan and its own executor (2 edges).** `crossnobis.R`
builds a crossnobis plan *and* runs it, calling `.run_geometry_compiler()` in
`execution-driver.R` and `.support_streamed_scheduled_crossnobis()` in
`kernel.R` from the plan layer. *Follow-up: split the crossnobis executor out
into layer 4, exactly as the geometry executor was split out of the compiler.*

**A source session validates the task it is handed (1 edge).**
`.open_effect_task_source_session()` in `source.R` calls
`.validate_compiled_effect_task()`, which now sits beside the task value it
validates in `evidence-task.R`. *Follow-up: hand the session an
already-validated task so a layer-2 file stops re-checking a plan-layer value.*

**The executor constructing its own results (3 edges).**
`execution-driver.R` calls `effect_form()`/`effect_geometry()`/`effect_view()`
and `.new_effect_contrast_view()`; `storage.R` calls
`.effect_form_codec_format()`. *Follow-up: a small result-builder in layer 4,
or explicit constructor entry points handed to the executor.*

**`evidence-api.R` is two files wearing one name (3 edges).** It is the public
facade *and* the home of `measurement_frame()`, `measurement_form()`,
`edge_frame()` and `geometry_alignment()`. Anything that needs those
constructors must call the facade. *Follow-up: move the constructors down into
the plan and view layers and leave only the facade.*

### What left the register

Twelve of the original twenty-one edges are gone, and two changed target
because the code they named moved:

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

## What is not in the register

`primitives.R`, `kernel.R`, `task.R`, `compiler.R`, `geometry-plan.R`,
`evidence-task.R`, `geometry-entry.R`, `views.R`, and `result.R` have **zero**
upward edges; `execution-driver.R` has only the two registered result-builder
edges. The largest strongly connected component is down from the 26 files it
held before the layer-4 split to 17,
and the compiler, the executor, the kernel, the geometry plan, the results,
and the views are all outside it. What remains is the layer-2 value tangle
around `relation.R`, `source.R`, and `receipt.R`, which this pass did not
touch.
