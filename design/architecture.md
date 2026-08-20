# `crossform` source layering

Status: enforced by `tests/testthat/test-architecture.R`
Date: 2026-08-17
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
| 2 | **values** | `domain.R`, `frame.R`, `pairing.R`, `relation.R`, `relation-session.R`, `effect-space.R`, `effect-map.R`, `metric.R`, `metric-learning.R`, `source.R`, `capabilities.R`, `study.R`, `study-facts.R`, `design-model.R`, `observation-model.R`, `extractor.R`, `scope.R`, `support-index.R`, `numerics.R`, `query-structured.R`, `pair-query.R`, `operations.R`, `receipt.R`, `reliability.R`, `validation-memo.R`, `measurement.R`, `measurement-storage.R`, `relation-fit.R`, `residual-statistics.R`, `bridge.R`, `compute-policy.R`, `memory-plan.R`, `crossform-package.R` |
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
  `crossnobis.R` is the plan — `noise_precision()`, `plan_crossnobis()`, and
  the metric- and pairing-role checks — and what moved is everything that
  executes: the exported `crossnobis()` entry that reads a contrast off a
  compiled plan. The second runtime this file also carried at the time of the
  split — a planned receipt, a private learned-metric driver, and a dispatch
  on an `effect_crossnobis_plan` class — was retired (B3, 2026-08-20) once the
  learned local metric became an ordinary compiler lowering (B2); there is now
  one runtime, `.run_geometry_compiler()`, on both routes.
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

## Untangling layer 2

No file has an upward edge, and once that was true the remaining structural
question was a different one: not which calls run the wrong way between layers,
but which of the *values* is the more primitive, so that mutual recursion
inside layer 2 has a direction to be given. Two passes answered it. The first
took the receipt out of the cycle; the second wrote the order down and moved
everything that disobeyed it.

### The receipt left the tangle

The fifteen-file component was, in part, an accident of one call. `receipt.R`
entered it through a single outgoing edge — `receipt.R -> memory-plan.R`, for
`memory_plan()`. `.validate_memory_plan_for_receipt()` rebuilt a memory plan
from the plan's own categories and compared the derived totals, which is a real
integrity check: it is what refuses a receipt whose `modeled_workspace_bytes`
was edited after the fact. But it made the receipt a caller of the vocabulary
it is supposed to be a record *of*.

The record — its scalar validation, its byte arithmetic, and the classed object
itself — now lives in `.workspace_plan_record()` in the primitive layer
(`R/primitives.R`). `memory_plan()` is the named-argument face of that record;
the receipt rebuilds through the same primitive. The check is unchanged, with
identical arithmetic and identical errors, but it is now reachable from both
files without joining them.

One edge was worth six files. `capabilities.R`, `measurement.R`,
`memory-plan.R`, `metric.R`, and `scope.R` were in the cycle only by way of
`receipt.R`, so removing it dropped the component from 15 to 9.

The three in-edges that were left — `relation.R`, `relation-fit.R`, and
`study-facts.R` calling `source_capabilities()` and
`.validate_source_capabilities()` — have since been retargeted rather than
kept, because they were pointing at the wrong file. A receipt is where a
capability is finally *written down*; it is not where a capability is
*defined*. Three value files have to state what their sources can do long
before any receipt exists, and making them ask the receipt for the vocabulary
inverted the record and the thing recorded. `source_capabilities()` and its
validator now live in `R/capabilities.R`, the file already named for them and
already holding the admission check `.relation_source_capabilities()` that
consumes them; `receipt.R` calls *down* to the same validator to canonicalize
the sources it stores. `tests/testthat/test-architecture.R` fails on any new
edge from `receipt.R` into the value files, and — because a legitimate
downward edge such as `receipt.R -> capabilities.R` cannot be told from an
illegitimate one by a name list — also fails if anything the receipt calls can
reach the receipt again by any route.

### The value order

Removing the receipt left one tangle, in layer 2 and entirely within it: nine
mutually recursive value files (`effect-map.R`, `effect-space.R`,
`extractor.R`, `operations.R`, `pairing.R`, `relation-fit.R`, `relation.R`,
`source.R`, `study-facts.R`). That was not a layering violation — a file may
call sideways — so there was no direction to restore, only a knot to untie,
and untying it meant first deciding which of `relation`, `source`, and
`effect-space` is the more primitive value.

That decision is recorded here. **Within layer 2 the value vocabulary has an
order, lowest first:**

| # | Band | Files |
|---|------|-------|
| 1 | **spaces and domains** | `effect-space.R`, `domain.R` |
| 2 | **maps and extractors** | `effect-map.R`, `extractor.R` |
| 3 | **sources and what they can do** | `source.R`, `capabilities.R` |
| 4 | **relations** | `relation.R` |
| 5 | **fits, and sessions over a relation's sources** | `relation-fit.R`, `relation-session.R` |
| 6 | **queries, pairings, frames, metrics** | `operations.R`, `pairing.R`, `frame.R`, `query-structured.R`, `pair-query.R`, `metric.R` |
| 7 | **study facts, and what quotes them** | `study-facts.R`, then `study.R`, `design-model.R`, `reliability.R`, `residual-statistics.R` |
| 8 | **records of an execution** | `receipt.R` |

and above them, in layer 3, the evidence tasks and plans that quote all of it.

The order is not enforced numerically — a value file may still call sideways,
and bands 2 and 3 are genuinely incomparable (nothing connects `extractor.R`
and `source.R` in either direction). What it buys is a criterion: when two
value files call each other, the one lower in this table is the one that keeps
the shared function. Each arrow below is the edge that actually exists in
`Rscript benchmarks/call-graph-scc.R`.

**Spaces are below everything.** `effect-space.R` calls nothing in the package
above layer 1 — not one edge — and six files in five bands above it call in:
`effect-map.R -> effect-space.R` (`effect_space`, `.as_effect_space`,
`.validate_effect_space`, `.validate_effect_provenance`,
`.validate_effect_names`), `extractor.R -> effect-space.R` (`effect_space`,
`.as_effect_space`), `relation.R -> effect-space.R` (the first four of those
plus `.same_effect_space`), and `relation-fit.R`, `metric.R`, and
`study-facts.R -> effect-space.R` (`.validate_effect_provenance`). A space is
what every other value is *expressed in*, so nothing it names may be something
built out of it.

**Maps and extractors are below relations.** `relation.R -> extractor.R`
(`effect_extractor`, `.validate_effect_extractor`) and
`relation-fit.R -> extractor.R` (`.compile_lm_estimator`,
`.resolve_observation_whitener`). A relation carries a per-partition extractor;
an extractor knows nothing about partitions. There is no edge back.

**Sources are below relations.** `relation.R -> source.R`
(`.source_descriptor`, `.validate_source_descriptor`,
`.with_source_descriptor`), `relation-fit.R -> source.R`
(`.validate_source_features`), `study-facts.R -> source.R`
(`.validate_source_descriptor`). A relation *is* a set of partition-keyed
sources plus the extractors that read them, so it names sources; a source is a
descriptor, a handle, and a read, and has no idea what a partition is. Since
the session split below, `source.R` calls nothing in the package outside
layer 1.

**Capabilities sit with sources.** `capabilities.R` calls nothing above layer
1, and `relation.R`, `relation-fit.R`, `study-facts.R`, and `receipt.R` all
call into it. What a source can do is a fact about the source.

**Fits are above relations.** `relation-fit.R -> relation.R` (`relation`,
`.validate_relation`, `.relation_family_identity`). A fit is a relation plus an
error channel; the reverse reading — a relation that knows what a fit is — is
what produced the one back-edge in this band, and it is gone (see below).

**Pairings are above relations and above operations.**
`pairing.R -> relation.R` (`.validate_relation`) and
`pairing.R -> operations.R` (`.new_partition_reducer`). A pairing says which
partitions of which relations are compared and how the results are reduced, so
it names both.

**Study facts are above the vocabulary.** `study-facts.R` calls `relation.R`,
`source.R`, `capabilities.R`, and `effect-space.R`, and nothing in bands 1–6
calls `study-facts.R` — the layer-2 files that do (`study.R`,
`design-model.R`, `reliability.R`, `residual-statistics.R`) sit above it in
turn. Facts about a study are stated in the vocabulary; the vocabulary does
not consult the facts.

**The receipt is above everything it records**, which is the previous section.

#### What untying it took

Nine files, seven changes, no behaviour change. Four were functions defined in
the wrong band, each one holding a cycle open by itself:

| Function | From | To | Why |
|---|---|---|---|
| `.validate_effect_names()` | `extractor.R` | `effect-space.R` | The rule for what counts as a legal set of coordinate names belongs to the space, not to the first consumer that happened to need it. Written in `extractor.R`, it made `effect_space()` — the value an extractor is declared *against* — call up into its own consumer. `kernel.R`, `task.R`, `memory-plan.R`, `design-model.R` and `effect-map.R` need the same rule, and none of them is an extractor. |
| `.validate_nonempty_id()` | `effect-map.R` | `check.R` (layer 1) | A two-line wrapper over `.check_string()` with no effect-map content in it, called from five files. It made `study.R`, `study-facts.R`, `design-model.R` and `observation-model.R` appear to depend on effect maps in order to spell a string. |
| `.validate_partition_reducer()` | `pairing.R` | `operations.R` | The record `.new_partition_reducer()` and the only statement of what makes it canonical sat on opposite sides of a two-file cycle: `operations.R` called up to `pairing.R` to check a value `pairing.R` had called down to `operations.R` to build. |
| `source_capabilities()`, `.validate_source_capabilities()` | `receipt.R` | `capabilities.R` | See the previous section. |

One was a file that was two files:

- **`R/relation-session.R`** (new, layer 2, band 5). Two hundred lines —
  `.open_relation_source_session()`, `.close_source_session()`,
  `.open_two_sided_source_session()`, `.open_effect_task_source_session()` —
  moved out of `source.R`. A source session is not source vocabulary: it is
  what a *relation* does with its sources, opening one handle per distinct
  descriptor, reading feature blocks through it, and closing each exactly once.
  While it sat in `source.R` its single `.validate_relation()` call ran against
  the order — the file defining a source reached up to the file defining a
  relation over sources — and that one edge held `source.R` inside the tangle.
  Now `relation-session.R -> source.R` and `relation-session.R -> relation.R`
  both run downward. Nothing moved but the text: the functions, their
  arguments, and their errors are unchanged, and the tests that reach them
  through `crossform:::` did not move either.

One was a genuine dependency inversion, resolved with a hook:

- `relation_block()` accepts an `effect_relation` *or* an
  `effect_relation_fit`, and unwrapped the second by calling
  `.validate_relation_fit()` — `relation.R -> relation-fit.R`, straight against
  band 4 → 5. Knowing how to check a fit is the fit's business. `relation.R`
  now declares `.as_read_relation()`, a generic whose default is the identity,
  and `relation-fit.R` defines `.as_read_relation.effect_relation_fit()`,
  which runs exactly the validation that used to run in `relation.R` and
  returns `x$relation`. The caller sees no difference — the same message, with
  the same condition class, for a fit whose fields are wrong — and the arrow
  now points down. The generic is internal and is only ever called from inside
  the namespace, so it needs no `S3method()` line in `NAMESPACE`; dispatch
  finds the method in the package environment.

The last was not a call at all. Three of the component's file-to-file edges
came from local names shadowing package functions, which the call-graph rule
cannot tell from a reference — the same defect as the `geometry_alignment`
local that the `evidence-api.R` split turned up:

| Apparent edge | Actually |
|---|---|
| `extractor.R -> relation-fit.R` | `.compile_lm_estimator()` assigned locals named `effect_covariance` and `residual_df`, which are also two exported accessors on a fit |
| `extractor.R -> study-facts.R` | `.resolve_observation_whitener()` took a parameter named `observations`, which is also the exported study fact |
| `relation-fit.R -> operations.R` | `.error_capabilities()` assigned a logical named `covariance`, which is also the edge normalizer `covariance()` |

They are now `coordinate_covariance`, `residual_degrees`, `n_observations`,
and `has_covariance`. A fourth shadow — `.open_relation_source_session()`'s
`relation` parameter — needed no rename, because that function moved to
`relation-session.R`, where a dependency on `relation.R` is real and points
down. Renaming the other four is worth doing on its own terms and not only for
the graph: a reader of `effect_covariance <- tcrossprod(...)` cannot tell the
local from the accessor of the same name, and neither can a reader of the
call graph.

## What remains

**Layer 2 contains no cycle at all.** Not a smaller one — none. The value
vocabulary is a directed acyclic graph, and the order above is the topological
sort of it.

The largest strongly connected component anywhere in the package is down from
the 26 files it held before the layer-4 split, to 17 after it, to 15, to 9
after the receipt left, to **1** now: every file in `R/` is in a component of
its own. The last multi-file component outside layer 2 was the sampling trio
(`evidence-sampling.R`, `evidence-sampling-kernel.R`,
`evidence-sampling-product.R`), and it went the same way as three of the nine
above — it was never a call. `.require_sampling_covariance()` read its
capability with `$sampling_covariance`, which puts the name of an entry point
one layer up into the file's syntax tree; extracting with
`[["sampling_covariance"]]` says the same thing and closes the phantom
plan → product → kernel → plan cycle.

### The evidence-sampling triple

The edges among the three files, counted by symbol
(`SCC_FOCUS=… Rscript benchmarks/call-graph-scc.R`). Before:

```
evidence-sampling.R          -> evidence-sampling-product.R    1   phantom
evidence-sampling-product.R  -> evidence-sampling.R            9
evidence-sampling-product.R  -> evidence-sampling-kernel.R     5
evidence-sampling-kernel.R   -> evidence-sampling.R            2
```

After — the same 16 real edges, minus the phantom, and now a DAG:

```
evidence-sampling-product.R  -> evidence-sampling.R            9
evidence-sampling-product.R  -> evidence-sampling-kernel.R     5
evidence-sampling-kernel.R   -> evidence-sampling.R            2
```

Nothing moved. The one edge that closed the loop was never a call, and the
remaining three run one way: product depends on kernel, kernel depends on
plan. That is the data flow read backwards — a plan is compiled in
`evidence-sampling.R`, the kernel builds the covariance form from it, and the
product specializes both to a relation fit — so each file can now be read
knowing only the files below it.

The order is layer 3's internal order, so the layering test cannot see it: all
three files are on the same layer and are free to call sideways. It is held
instead by `test_that("the evidence-sampling triple is a DAG")` in
`tests/testthat/test-architecture.R`, which takes the transitive closure of
the induced three-file subgraph and fails if any of the three reaches itself.
A mutual pair fails it as surely as the full loop, and the failure prints the
induced edges rather than the file names, so it names the call to move.

Note the measurement's one bias, since two of the phantoms above came from it
and 44 like them remain: the extractor counts every symbol a file mentions
that resolves to a definition elsewhere, including the field name in `x$f`. It
therefore over-reports. That is the safe direction — a graph it calls acyclic
is acyclic — but a component it reports is worth reading before it is
believed.

That is the whole file-level dependency structure: six layers, an order inside
layer 2, and no loop anywhere. What the next pass should defend is not a
number but the shape — see `tests/testthat/test-architecture.R`, which now
fails on an upward edge, on a stale register entry, and on any route from the
receipt back to itself.
