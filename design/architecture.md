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
| 2 | **values** | `domain.R`, `frame.R`, `pairing.R`, `relation.R`, `effect-space.R`, `effect-map.R`, `metric.R`, `metric-learning.R`, `source.R`, `study.R`, `study-facts.R`, `design-model.R`, `observation-model.R`, `extractor.R`, `scope.R`, `support-index.R`, `numerics.R`, `query-structured.R`, `pair-query.R`, `operations.R`, `receipt.R`, `reliability.R`, `validation-memo.R`, `measurement.R`, `measurement-storage.R`, `relation-fit.R`, `residual-statistics.R`, `bridge.R`, `crossform-package.R` |
| 3 | **plans** | `geometry-plan.R`, `relation-plan.R`, `crossnobis.R`, `coupling-plan.R`, `evidence-task.R`, `evidence-sampling.R`, `evidence-sampling-kernel.R`, `evidence-sampling-product.R`, `compiler-conformance.R` |
| 4 | **compiler / execution** | `compiler.R`, `kernel.R`, `task.R`, `execution.R`, `memory.R`, `storage.R`, `measurement-kernel.R` |
| 5 | **results / views** | `result.R`, `views.R`, `coupling-views.R`, `tomography.R`, `measurement-result.R`, `measurement-decomposition.R`, `format-results.R`, `print-methods.R`, `plot-methods.R` |
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

## The register: 21 upward edges that remain

Each is real debt with a named follow-up. They are listed in `allowed_upward`
in the test with the same comments; the summary is:

**Declaration constructors filed under "execution" (7 edges).**
`compute_policy()` (`execution.R`) and `memory_plan()` (`memory.R`) are pure
declarations — `compute_policy()` explicitly "never starts a worker itself" and
`execution.R` calls nothing but `conditions.R` and `message-helpers.R`. Values
and plans that declare a budget therefore appear to call upward. *Follow-up:
move both constructors into the values layer. Nothing else moves.* This is the
cheapest item on the list and would delete a third of the register.

**`.compiler_capabilities()` (2 edges).** A version/capability record living in
`compiler.R`. *Follow-up: extract to a leaf `R/capabilities.R`.*

**`.validate_compiled_effect_task()` (1 edge).** A value validator living in the
compiler. *Follow-up: move it beside the task value it validates.*

**Plans reaching into the compiler for task construction and task identity
(4 edges: `crossnobis.R`, `evidence-task.R`, `geometry-plan.R`).** The largest
item. `.compile_effect_evidence_task()`, `.effect_task_id()`,
`.effect_task_semantic()` and friends build and name a task; that is plan work,
not execution work. *Follow-up: extract task construction and identity from
`compiler.R` into a plan-layer file, leaving `compiler.R` to lower an
already-identified task.*

**The compiler constructing its own results (3 edges).** `compiler.R` calls
`effect_form()`/`effect_geometry()`/`effect_view()` and
`.new_effect_contrast_view()`; `storage.R` calls `.effect_form_codec_format()`.
*Follow-up: a small result-builder in layer 4, or explicit constructor entry
points handed to the compiler.*

**`.validate_geometry_metric_schedule()` (1 edge).** *Follow-up: move the
schedule validator down to `metric.R`.*

**`evidence-api.R` is two files wearing one name (3 edges).** It is the public
facade *and* the home of `measurement_frame()`, `measurement_form()`,
`edge_frame()` and `geometry_alignment()`. Anything that needs those
constructors must call the facade. *Follow-up: move the constructors down into
the plan and view layers and leave only the facade.*

## What is not in the register

`kernel.R` and `primitives.R` have **zero** upward edges. That is the property
this pass established and the property the test exists to keep.
