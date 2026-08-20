# `crossform` execution model

Status: architecture proposal; review corrections incorporated, not implementation
Date: 2026-08-12; section 16 (learned local metric as a compiler lowering) added 2026-08-17
Companion to: [crossform-package-design.md](crossform-package-design.md)

## Version 0.1 execution contract freeze

The audited baseline is
`01f3763c6aa00ebe50e69fc966d5201a861e4b09`. Before API freeze, execution is
governed by these non-negotiable rules:

- semantic validation compares exact effect-space and ordered feature-domain
  signatures before opening a source;
- the component requirement graph removes unrequested work before planning;
- the compiler searches legal feature, row, and coordinate tiles within
  `workspace_bytes` rather than testing one fixed default;
- the memory model includes the actual dense or sparse frame, resident sources,
  scoped handles, component-dependent local state, output, contraction and
  replacement copies, and active buffers;
- baseline RSS, incremental RSS, absolute peak RSS, and planned package
  workspace are distinct quantities in both benchmarks and receipts;
- each distinct descriptor is admitted and opened once per execution scope,
  reused across feature tasks, and closed with a checked outcome;
- the coordinator advances receipt task and byte counters only when canonical
  reduction commits a task; failure receipts therefore report genuine partial
  progress without exposing a partial scientific result;
- requested BLAS threads and observed BLAS threads are separate fields. Unknown
  observations stay unknown. Version 0.1 does not mutate ambient global thread
  or executor state to manufacture compliance;
- cleanup is successful only when every owned close, detach, and deletion is
  confirmed; and
- optional `neuroim2` integration supplies only domain and neighborhood index
  geometry. It never supplies the scientific kernel or an executor.

These rules supersede older `memory_budget` sketches below. The public policy
name is `workspace_bytes`; a future `rss_limit_bytes` may be admitted only with
portable observation and failure semantics. Full certification requires
component-elision probes, source lifecycle fault injection, measured memory
evidence, installed-artifact tests, conditional neuroim2 parity, rendered
documentation, and `R CMD check`.

## Outcome

`crossform` should not replace rMVPA’s `future` calls with another parallel map. For its additive fixed-bilinear core, it should eliminate the unit of work that made that machinery necessary.

rMVPA schedules regions or searchlights because each region is treated as a separate fit. `crossform` compiles additive diagonal-frame locations with fixed bilinear queries into matrix algebra:

\[
V=W\,\Phi_\Gamma(EY)\,C.
\]

For this admitted core, the execution engine therefore works over coarse feature blocks, not searchlight centers. It streams each neural feature once, forms its cross-generalized geometry once, and lets sparse multiplication distribute that evidence to every voxel, ROI, searchlight, parcel, or scale that uses it.

The execution principle is:

> Compile globally, stream featurewise, reduce deterministically, and make parallel execution an explicit resource policy rather than inherited session state.

`shard` is valuable here. The correction is not to discard shared memory or
supervised workers; it is to prevent them from becoming another scientific
engine. `crossform` should use `shard` as an optional data-plane and executor
adapter beneath the one relation/geometry kernel.

## What the current rMVPA execution model is compensating for

The live rMVPA audit at commit `3b12a8855b06f549e5772ccb9af070639a34752b` found:

- `future`, `future.apply`, and `furrr` are mandatory imports ([DESCRIPTION](https://github.com/bbuchsbaum/rMVPA/blob/3b12a8855b06f549e5772ccb9af070639a34752b/DESCRIPTION)).
- `mvpa_iterate()` serially builds batches of ROI matrices in the main process, then sends those ROI objects through `furrr::future_pmap()`; searchlight and regional analyses use different chunk heuristics ([mvpa_iterate.R](https://github.com/bbuchsbaum/rMVPA/blob/3b12a8855b06f549e5772ccb9af070639a34752b/R/mvpa_iterate.R#L808-L1054)).
- The worker copy of a model spec is manually stripped of its dataset to reduce serialization, while comments identify file-backed proxies as future work ([mvpa_iterate.R](https://github.com/bbuchsbaum/rMVPA/blob/3b12a8855b06f549e5772ccb9af070639a34752b/R/mvpa_iterate.R#L1090-L1114)).
- `run_future.default()` derives behavior from the globally active future plan via `future::nbrOfWorkers()`, constructs closures around model/process functions, and chooses chunks according to the historical analysis type ([mvpa_iterate.R](https://github.com/bbuchsbaum/rMVPA/blob/3b12a8855b06f549e5772ccb9af070639a34752b/R/mvpa_iterate.R#L1124-L1290)).
- An experimental `shard` backend successfully avoids much of the serial ROI extraction and repeated dataset serialization, but does so by duplicating worker execution and dataset-specific extraction logic ([shard_backend.R](https://github.com/bbuchsbaum/rMVPA/blob/3b12a8855b06f549e5772ccb9af070639a34752b/R/shard_backend.R)). The performance mechanism is worth retaining; the duplicate analysis path is not.
- CLI and custom entry points temporarily mutate and restore the process-global future plan ([utils.R](https://github.com/bbuchsbaum/rMVPA/blob/3b12a8855b06f549e5772ccb9af070639a34752b/R/utils.R#L299-L356)).
- Progress depends on `progressr`; remote shared-memory workers intentionally lose the local progress closure. Runtime failures become per-ROI or per-batch tibble rows.
- The execution-related source surface is large: the iterator, shard backend, engine selection, fast searchlight paths, global path, custom workflows, CLI, and utilities together exceed 8,000 physical R lines.

This architecture has been progressively hardened, but it is solving the consequences of per-ROI fitting: extraction, serialization, task overhead, load balancing, result reconciliation, and method-specific fast paths.

`crossform` removes that causal chain.

The current `shard` 0.2.0 source at commit
`233c71186ebac0e2e98eba70f10671450fc5e1be` already provides the operational
pieces that remain useful after that removal: immutable shared inputs,
memory-mapped/POSIX backing, explicit output buffers, supervised persistent
workers, deterministic chunk-ordered reduction, memory/copy diagnostics, and
worker recycling ([README](https://github.com/bbuchsbaum/shard/blob/233c71186ebac0e2e98eba70f10671450fc5e1be/README.md),
[`shard_map()`](https://github.com/bbuchsbaum/shard/blob/233c71186ebac0e2e98eba70f10671450fc5e1be/R/shard_map.R),
[`shard_reduce()`](https://github.com/bbuchsbaum/shard/blob/233c71186ebac0e2e98eba70f10671450fc5e1be/R/shard_reduce.R)).
Those are execution capabilities, not analysis identities.

## 1. Execution is separate from scientific semantics

The scientific plan consists of relation, frame, pairing, and query. Execution adds only resource policy:

```r
g <- materialize_geometry(
  rel,
  at = searchlights(6),
  over = cross_partitions(),
  compute = compute_policy(
    workers = 1,
    memory_budget = 8 * 1024^3
  )
)
```

`compute_policy()` is a small immutable value. It does not alter the session, create a global plan, or change the estimand. Proposed fields:

```text
workers
memory_budget
io_block_features       # NULL means compiler chooses
reduction_microblock    # fixed numerical grouping, initially same as I/O block
measurement_tile_rows   # NULL means compiler chooses
geometry_tile_width     # NULL means compiler chooses
max_inflight            # bounded task count
max_reorder_bytes       # bounded out-of-order results
process_backend         # sequential in v0.1
source_staging          # none, response, or relation
threads_per_worker
```

The policy is immutable, serializable, and contains no callbacks or destination
paths. Changing it must not change the estimand or output ordering. Version 0.1
guarantees agreement within a documented absolute-plus-relative tolerance, not
bitwise invariance across block sizes, BLAS implementations, or platforms.

Coordinator-side reporting and checkpoint destinations are separate execution
attachments:

```r
materialize_geometry(..., compute = policy, reporter = reporter, checkpoint = checkpoint)
```

They are excluded from scientific-plan and numerical-policy identity. Reporter
failure is recorded and disables further reporting; it does not cancel,
retry, or alter a scientifically valid computation unless the user explicitly
requests cancellation through a separate control channel.

## 2. The compiled stages

### Stage 0: validate and lower

Before reading response data, the compiler:

1. validates source/extractor dimensions and experimental coordinate identity;
2. canonicalizes pairing edges and weights;
3. compiles the spatial frame to sparse \(W\);
4. compiles requested linear views to \(C\);
5. chooses full geometry width \(h=q(q+1)/2\) or direct-query width \(k\);
6. chooses an in-memory or block-backed output representation without changing
   the semantic result type;
7. estimates output, contraction-temporary, serialization, reorder,
   checkpoint, source, and per-task memory;
8. chooses feature, measurement-row, and geometry-coordinate tiles satisfying
   the declared memory budget;
9. determines whether each source is coordinator-readable or safely reopenable by workers.

The compiler fails before execution if the conservative plan is impossible. It
reports the estimate and remedies: direct queries, block-backed output, smaller
tiles, fewer future workers, a file-backed source, or a larger budget. Planning
is validated against measured peak memory; it is never described as exact
portable allocation accounting.

### Stage 1: feature-block relation kernel

For canonical, disjoint feature block \(I_t\), a pure task computes:

\[
B_{r,I_t}=E_rY_r[:,I_t]
\]

for every required partition, followed by either:

\[
Z_{I_t}=\Phi_\Gamma(B_{I_t})
\]

for full packed geometry, or

\[
U_{I_t}=\Phi_\Gamma(B_{I_t})C
\]

for direct linear queries.

A task returns only:

```text
block identity
feature interval/indices
packed atoms or direct feature queries
partition relation blocks needed for coherent/marginal summaries
timing, bytes, and diagnostics
```

It never receives a searchlight, ROI object, model specification, output map, or arbitrary callback.

### Stage 2: tiled deterministic sparse contraction

The coordinator consumes completed feature blocks in canonical feature order.
A full \(m\times h\) multiplication is forbidden when it exceeds the temporary
budget. For measurement-row tile \(J_a\), geometry-coordinate tile \(K_b\), and
feature block \(I_t\), it accumulates:

\[
\mathcal G[J_a,K_b]\mathrel{+}=
W[J_a,I_t]Z_{I_t}[:,K_b]
\]

or the direct-query equivalent. Weighted local relations are likewise tiled by
measurement rows:

\[
U_r[J_a,:]\mathrel{+}=W[J_a,I_t]B_{r,I_t}^\top.
\]

An accumulator tile is updated in place where the storage backend permits, then
written once to its destination. The implementation must not express the hot
path as `total <- total + Wb %*% Z`, which can allocate both a full contraction
temporary and a replacement output. Geometry-coordinate tiling may require the
atom kernel to emit only \(K_b\), rather than materializing all \(h\) columns.

The coordinator then derives coherent geometry and the pairing-appropriate
endpoint or left/right marginals from \(U_r\).

This design has four important properties:

- workers never write the final overlapping spatial result; a later executor
  may permit disjoint intermediate atom/query-buffer writes;
- task outputs occupy disjoint feature intervals;
- the coordinator reduces in the same order regardless of completion order or worker count;
- no worker returns a full `measurement × geometry` partial matrix;
- contraction temporaries are bounded by the chosen row/coordinate tile, not
  by \(m\times h\).

The final sparse contraction stays in one process initially. Compiler-generated
row and coordinate tiles are algebraic contraction tiles, not scientific
searchlight jobs. Parallelizing this stage in R would risk oversubscription and
large partial-result transfers.

### Stage 3: assemble and certify

The coordinator creates a complete `effect_geometry` for full materialization,
or an `effect_view` for direct evaluation, only after every task succeeds and
all invariants pass:

- expected feature coverage is exact;
- no block is duplicated or missing;
- total/coherent dimensions and bases agree;
- point-frame and conservation checks pass when applicable;
- direct-query and materialization metadata are complete.

A partial geometry is never silently returned as a normal result, and a direct
view is never classed as complete geometry.

## 3. Parallelism hierarchy

Use parallelism in this order.

### First: participants

Participants are independent, have separate native spaces, and produce separate geometry objects. This is the cleanest outer parallel axis. The core should document that cohort workflows parallelize participants outside `materialize_geometry()` and avoid nested within-participant process pools.

### Second: feature blocks within one participant

When relation extraction or cross-atom computation is expensive, process disjoint feature blocks concurrently. Dynamic scheduling is safe because each block has a canonical identity and reduction order.

### Third: native numerical libraries

Sparse matrix multiplication and dense small-matrix kernels may use optimized compiled code. The executor must avoid multiplying R processes by multithreaded BLAS workers. Default `threads_per_worker = 1`; any other choice is explicit and recorded.

### Not an axis for the additive core: searchlight centers

Additive-frame searchlights are rows of \(W\), not jobs submitted to workers.
No additive searchlight scheduler exists. Locally trained classifiers, locally
estimated covariance, nonlinear normalization, and generic factor frames fall
outside this collapse theorem; later extensions may retain location-dependent
work, but must use a separately declared lowering rather than pretending to be
the additive contraction.

## 3A. Bounded scheduling and backpressure

> Not implemented. `compute_policy()` has no `max_inflight` or
> `max_reorder_bytes` field, and the package ships no scheduler. An unwired
> in-process simulator (`.simulate_bounded_schedule()`) existed until
> 2026-08-16 and was deleted so that test coverage stops overstating delivered
> capability; this section remains the design target.

Dynamic scheduling is allowed only behind a bounded reorder window. Let the
next reducible feature task be \(t\). The coordinator may hold completed tasks
with ids greater than \(t\), but enforces both:

```text
number of dispatched-but-unreduced tasks <= max_inflight
bytes of completed out-of-order results <= max_reorder_bytes
```

Both limits are derived from the aggregate memory budget. Dispatch pauses when
either limit is reached and resumes only after ordered reduction releases
space. A completed task is released immediately after its contribution has
been reduced. Optional checkpointing may spill a validated completed block, but
spill does not permit unbounded dispatch.

This makes canonical ordered reduction and dynamic load balancing compatible
without turning a slow early block into unbounded memory growth.

## 4. Source capability protocol

Parallelism cannot make data locality disappear. Each response source declares capabilities:

```text
seekable
reopenable
concurrent_read_safe
descriptor_serializable
in_memory
```

The compiler selects one access mode.

### Coordinator-read

For ordinary in-memory R objects or non-reopenable connections, the coordinator reads bounded blocks. Version 0.1 executes these locally rather than duplicating a whole dataset into every process. Vectorized algebra is expected to make this path fast.

### Worker-read

For file-backed or otherwise reopenable sources, workers receive only source descriptors and block indices, open their own read-only handles during initialization, and read blocks directly.

### Explicit shared source

A `shard` adapter may expose an in-memory response or materialized relation as
a read-only, worker-attachable source. It implements the same relation-source
protocol; it does not create a dataset subclass, model specification, geometry
class, or separate scientific engine.

There are two useful staging points:

1. **Share the response \(Y_r\).** Workers attach once, read feature blocks, and
   apply \(E_r\) themselves. This avoids copying large raw/runwise matrices and
   retains fully lazy extraction.
2. **Share the compiled relation \(B_r=E_rY_r\).** Pay one extraction pass, then
   reuse the much smaller \(q\times p\) relation across frames, queries, and
   repeated geometry plans. This is a transient compiler cache, not a required
   beta-file stage.

The compiler should choose neither silently in the first release. An explicit
staging declaration records which object was shared, its dimensions, backing,
owner, and lifetime. Later, an evidence-backed `"auto"` policy may choose from
estimated bytes, \(n_r/q\), extraction cost, and expected reuse.

Phase-E2 API sketch, not version 0.1:

```r
compute_policy(
  workers = 8,
  process_backend = "shard",
  source_staging = "relation"
)
```

The exact public spelling can wait; the semantic choices cannot.

There is no automatic “try shared memory, then silently fall back” behavior. The selected access mode appears in the execution receipt.

All published shared inputs are immutable. Workers may attach and detach, but
only the creating coordinator owns unlinking. The owner remains live until all
tasks have joined; a recycled worker must reopen current descriptors rather
than reuse cached handles from an earlier generation.

## 5. Executor ownership and backend boundary

The core owns no global scheduler state.

### Version 0.1

- `workers = 1`: deterministic in-process streaming; always available.
- any other value fails before source access with `workers > 1 is not implemented
  in crossform 0.1`; it never silently ignores the request or inherits an
  ambient pool.
- Phase E2 may admit `workers > 1` after the benchmark gate. The preferred first
  local executor is an optional `shard` adapter, because it already owns
  shared-input transport, bounded outputs, worker supervision, recycling, and
  memory diagnostics.

The `materialize_geometry()` call creates or explicitly receives an execution scope, initializes workers once, and closes only resources it owns. It does not inspect or modify `future::plan()`. A reusable caller-owned `shard` pool may be supported only through an explicit scope object with unambiguous ownership; the default call never attaches to ambient state.

### Backend-neutral task protocol

Internally, every task is a small serializable manifest and is executed by one package function with no discovered globals:

```text
run_geometry_task(plan_ref, task_id)
```

This permits adapters for `shard`, `future`, `mirai`, batch schedulers, or remote systems without changing kernel semantics. An adapter schedules fixed tasks and returns fixed task results; it may not alter relations, frames, pairings, queries, reduction, or errors.

`shard` is the leading local adapter because its contracts match this workload;
it should remain an optional dependency rather than a mandatory ontology.
`future` is not imported by the core and no adapter is promised. It may be
evaluated later only if it can provide explicit executor ownership rather than
inherit the process-global plan.

## 5A. Shard without a shard backend

The intended integration is deliberately asymmetric:

```text
crossform owns                    shard owns
-------------------------------   -------------------------------
relation/frame/pairing/query       shared segments and descriptors
feature-block task manifest       worker pool and supervision
geometry kernel                   task dispatch and recycling
canonical scientific reduction    memory/copy telemetry
effect_geometry result             optional bounded buffers
```

There is no `run_geometry.shard_model`, no `shard_geometry` result, and no
dataset-specific shared extractor. The same `run_geometry_task()` used by the
sequential executor is dispatched over coarse feature blocks. A worker receives
only a task id, immutable compiled metadata, and shared/reopenable source
descriptors.

The first adapter should use `shard_map()` with a small number of coarse blocks,
not thousands of searchlight jobs. Returned relation/atom blocks are bounded
and reduced by `crossform` in canonical feature order. If gathered blocks
would exceed the memory plan, the adapter may write disjoint atom/query rows to
a `shard::buffer()` and let the coordinator scan that buffer in order. Workers
still never write the final overlapping spatial geometry.

`shard_reduce()` is attractive for scalar or naturally associative diagnostics,
but the geometry accumulator should not be forced into its API unless the
combine law and ordered floating-point semantics are exact. Scientific
reduction remains `crossform`'s responsibility.

Single-worker execution bypasses `shard` completely unless the user explicitly
requests a reusable shared staging object. This avoids paying process and
shared-segment overhead when direct streaming is superior.

## 6. Memory model

The user declares one aggregate execution-memory budget. The compiler accounts for:

\[
\begin{aligned}
M_{\mathrm{coordinator}}={}&
M_W+M_{\mathrm{local\ relations}}+M_{\mathrm{overhead}},\\
M_{\mathrm{task}}={}&
M_{Y\mathrm{block}}+M_{B\mathrm{block}}+M_{\mathrm{atom/query\ block}},\\
M_{\mathrm{peak}}\approx{}&
M_{\mathrm{coordinator}}+
M_{\mathrm{shared}}+
n_{\mathrm{active}}M_{\mathrm{task}}+
M_{\mathrm{contraction}}+
M_{\mathrm{serialization}}+
M_{\mathrm{reorder}}+
M_{\mathrm{checkpoint}}+
M_{\mathrm{result\ buffer}}.
\end{aligned}
\]

Here:

- \(M_{\mathrm{contraction}}\) includes the dense row-by-coordinate tile
  produced by sparse multiplication and any replacement/copy required by R;
- \(M_{\mathrm{serialization}}\) covers task payload overlap during send/receive;
- \(M_{\mathrm{reorder}}\) is capped by `max_reorder_bytes`;
- \(M_{\mathrm{checkpoint}}\) covers compression/write buffers, not total disk
  capacity;
- \(M_{\mathrm{result\ buffer}}\) is either the in-memory result or bounded
  pages/cache for block-backed storage.

The planner distinguishes shared mapped bytes from private RSS and reports both.
Peak-memory evidence must aggregate the coordinator and workers; low coordinator
RSS alone is not proof that execution is memory bounded. Block size is selected
so estimated aggregate private memory stays below a conservative fraction of
the budget, while shared-segment size is reported separately. The receipt
records estimates, actual block and tile dimensions, worker peak RSS, hidden
copy/materialization diagnostics, and end-of-run memory return when the
executor supplies them.

There is no user-facing “ROI batch size” whose meaning changes with analysis
type. There are orthogonal I/O-feature, reduction-microblock,
measurement-row, and geometry-coordinate tiles, normally compiler-selected.

The engine does not call `gc()` as a scheduling mechanism. It controls reachability through short-lived block scopes and bounded buffers.

## 7. Determinism and random numbers

The core kernel is deterministic and uses no random numbers. Therefore it does not need `future.seed`, worker RNG streams, or schedule-dependent seeding.

If a randomized frame, pairing, bootstrap, or permutation is needed, it is constructed upstream as an explicit immutable object. Its realized indices/weights and seed provenance enter the plan hash. Execution merely contracts it.

Determinism laws:

1. output is identical across block sizes within numerical tolerance;
2. output is identical across worker counts within the same tolerance;
3. task completion order cannot change reduction order;
4. repeated plans read the same source revisions or fail provenance validation;
5. thread/process counts do not change row ordering or diagnostics.

These are three distinct guarantees:

- **scheduling invariance:** completion order does not determine semantic order;
- **block-partition agreement:** different legal block sizes agree within the
  declared tolerance, but may group floating-point operations differently;
- **cross-platform agreement:** separate tolerance contract conditioned on R,
  BLAS, sparse backend, and precision.

Version 0.1 guarantees the first and tolerance-qualified versions of the latter
two. A stronger later mode may separate large I/O blocks from fixed reduction
microblocks and use a canonical compensated/pairwise tree independent of worker
count and memory policy. Until then, no bitwise-invariance claim is permitted.

## 8. Failure semantics

### Compile-time failures

Dimension, basis, frame, pairing, estimability, unit, and memory problems fail before workers start.

### Runtime failures

A task failure contains:

```text
task id
stage
feature interval
source identity
condition class and message
worker/process identity
timing and retry count
```

Default behavior is fail-fast at task boundary:

- stop dispatching new work;
- cancel/close only owned workers;
- detach worker source handles, then close/unlink coordinator-owned shared segments;
- preserve valid checkpoint shards if checkpointing was requested;
- return neither a successful `effect_geometry` nor `effect_view`.

The raised `crossform_execution_error` carries the partial execution receipt,
cleanup state, and checkpoint manifest reference. A receipt is therefore
available even though no successful geometry exists to contain it.

No failure is converted into thousands of result rows. No backend silently falls back to a different backend. No incomplete field is printed as successful.

Invalid measurement rows discovered from metadata are compact diagnostics and are excluded before execution. Unexpected numerical/runtime failure is exceptional, not a per-location result type.

Shared-memory lifecycle is part of correctness, not best-effort cleanup:

- the coordinator is the sole segment owner and unlink authority;
- workers receive descriptors and attach read-only;
- every dispatch generation has a distinct identity, so a persistent or
  recycled worker cannot reuse a stale handle;
- owner close occurs only after all worker activity has stopped;
- success, error, interrupt, timeout, and retry exhaustion exercise the same
  idempotent cleanup path;
- an execution receipt distinguishes clean release from leaked/uncertain
  release, and the latter is a failed run.

## 9. Checkpointing and resume

> Not implemented. No public entry point accepts a `checkpoint` argument, and
> `receipt$observed$checkpoint` is permanently `enabled = FALSE`. Unwired
> shard read/write helpers existed until 2026-08-16 and were deleted; the
> `checkpoint_buffer` memory-budget category is retained because it is part of
> the hashed `memory_plan()` schema, not because the package checkpoints.

Checkpointing is optional and should arrive after the base kernel.

Each completed feature task writes one atomic, content-addressed shard:

```text
plan-hash/
  manifest.rds
  tasks/
    000001.rds
    000002.rds
    ...
```

Before dispatch, checkpointing preflights estimated total disk need, free space,
destination writability, restrictive permissions, and same-filesystem temporary
placement. Each task artifact contains a schema version, task/reduction
partition, dimensions, precision, and checksum. Temporary writes occur inside
the destination filesystem and are fsynced/closed before atomic rename.

Resume accepts a shard only when scientific plan hash, numerical execution
policy, package/kernel/schema version, strong source revision, exact feature and
coordinate tiles, precision, query/materialization policy, and reduction plan
match. A path plus modification time is not a strong source revision; accepted
sources need content hashes, immutable object/version ids, or an adapter-specific
integrity receipt.

Task shards contain intermediate feature blocks, not partial scientific maps.
They may contain sensitive derived neural data and therefore use owner-only
permissions by default. Collection reruns the identical canonical reduction
plan. “Resume equality” means tolerance-qualified equality to uninterrupted
execution with the same task partition, coordinate tiles, precision, kernel,
source revisions, and reduction plan.

On success, task shards are removed by default; `keep_checkpoint = TRUE` is explicit.

## 10. Progress and observability

Workers do not own UI and do not serialize progress closures. They return structured receipts. The coordinator emits events:

```text
compiled
worker_started
block_started
block_completed
bytes_read
reduction_advanced
checkpoint_written
completed
failed
cancelled
```

Progress is weighted by planned work or processed feature/frame nonzeros, not
raw searchlight count. A base text reporter can be tiny; notebooks and GUIs can
supply a coordinator-side observer that is not part of `compute_policy()` or
its hash. Reporter errors are caught, recorded, and isolated from scientific
execution.

Every completed geometry or direct view records an execution receipt:

```text
scientific plan hash
numeric execution policy hash
executor and source access modes
workers and threads per worker
memory budget and chosen block size
tasks planned/completed/retried
features and bytes processed
stage timings
checkpoint/resume status
reporter identity and reporter failures (non-semantic)
R, package, Matrix, platform, and BLAS identity
completion status
```

This replaces an elaborate logging framework with one durable, inspectable account of what happened.

## 11. Query-aware execution

Query-aware execution has two explicit public contracts. It is never a hidden
memory-policy branch inside `materialize_geometry()`.

### Complete geometry

```r
g <- materialize_geometry(
  rel,
  at = frame,
  over = pairing,
  materialize = "full",
  storage = "memory"       # or "block"
)
```

This returns `effect_geometry`. `storage` changes physical representation, not
semantic completeness: all packed coordinates required by the declared
components exist and may be queried later.

### Direct query

```r
v <- evaluate_geometry(
  rel,
  at = frame,
  over = pairing,
  query = geom_contrast(...)
)
```

This compiles \(U=\Phi_\Gamma(EY)C\) featurewise and returns `effect_view`, not
`effect_geometry`. It records the exact query, requested components, units,
source/relation identity, and execution receipt. It cannot later answer an
unmaterialized query without rerunning the relation.

### One or a few linear queries

Use `evaluate_geometry()`. Compute \(U=\Phi_\Gamma(EY)C\) featurewise and
contract \(WU\). Do not materialize the packed geometry field.

### Many views or nonlinear views

Use `materialize_geometry(materialize = "full")`. Materialize packed geometry once. RDM,
linear RSA, contrast energy, and later views reuse it without rerunning source
extraction.

### Coherent/configuration requested

Retain the bounded partition relation blocks long enough to accumulate \(WB_r^\top\). Configuration is derived after total and coherent geometry complete.

### Total-only requested

Do not calculate local relation marginals or coherent geometry.

The execution plan is thus driven by algebraic requirements rather than method
names, while the return type remains honest about what was computed.

## 12. Why this is faster even before parallelism

The largest improvement is work elimination:

- each feature relation is extracted once per block, not once per overlapping searchlight;
- each cross-Gram atom is formed once, not once per containing region;
- RSA regressions compile once into a query matrix;
- multiple views reuse one geometry;
- result assembly is dense/sparse matrix placement, not tibble/list-column combination;
- no model, dataset, ROI, prediction table, or closure is serialized per task;
- no engine-eligibility ladder chooses among duplicate scientific implementations.

Parallelism should be introduced only after the vectorized sequential kernel is benchmarked. A simple one-process contraction may beat a large process pool because the old parallel workload largely consisted of redundant computation and serialization.

## 13. Thirty execution decisions, critically evaluated

Scores: impact/effort/risk are 1–5; confidence is architectural confidence.

| # | Decision | I | E | R | Conf. | Verdict |
|---:|---|---:|---:|---:|---:|---|
| 1 | Eliminate additive-frame per-searchlight tasks | 5 | 3 | 2 | 99% | Keep with scope |
| 2 | Deterministic sequential default | 5 | 1 | 1 | 99% | Keep |
| 3 | Explicit immutable compute policy | 5 | 2 | 1 | 97% | Keep |
| 4 | No mandatory `future`/`furrr` | 5 | 2 | 1 | 98% | Keep |
| 5 | Canonical disjoint feature-block tasks plus contraction tiles | 5 | 4 | 2 | 94% | Keep |
| 6 | Source capability negotiation | 5 | 3 | 2 | 95% | Keep |
| 7 | Canonical coordinator reduction order with tolerance contract | 5 | 3 | 2 | 94% | Keep |
| 8 | Query-aware direct/full lowering | 5 | 3 | 2 | 95% | Keep |
| 9 | Participant-level outer parallelism | 4 | 2 | 1 | 96% | Keep/document |
| 10 | Prohibit implicit nested parallelism | 5 | 2 | 2 | 96% | Keep |
| 11 | Conservatively planned feature/row/coordinate tiles | 5 | 4 | 3 | 75% | Keep after measured-memory validation |
| 12 | Cost model using bytes, \(q\), edges, and frame nnz | 4 | 3 | 2 | 90% | Keep |
| 13 | Worker-readable file-backed sources | 4 | 3 | 2 | 92% | Keep protocol |
| 14 | Optional `shard` local executor adapter | 5 | 3 | 3 | 90% | Preferred after benchmark |
| 15 | `future` executor adapter | 3 | 3 | 3 | 65% | Defer |
| 16 | `mirai` executor adapter | 3 | 3 | 3 | 65% | Defer |
| 17 | `shard`-backed shared relation source | 5 | 3 | 3 | 90% | Keep optional adapter |
| 18 | Fork-only multicore default | 2 | 2 | 5 | 25% | Reject |
| 19 | Dynamic scheduling with bounded reorder/backpressure | 4 | 3 | 2 | 80% | Keep after simulator |
| 20 | Coordinator-side structured progress events | 4 | 2 | 1 | 96% | Keep |
| 21 | Durable execution receipt | 5 | 3 | 1 | 97% | Keep |
| 22 | Atomic checkpoint/resume | 4 | 4 | 3 | 85% | Phase 2 |
| 23 | Automatic arbitrary-error retries | 2 | 2 | 5 | 20% | Reject |
| 24 | Silent backend fallback | 1 | 1 | 5 | 5% | Reject |
| 25 | Per-location runtime error rows | 2 | 3 | 5 | 10% | Reject |
| 26 | Fail-fast task-boundary semantics | 5 | 2 | 1 | 98% | Keep |
| 27 | Ordinary views on partial geometry | 1 | 2 | 5 | 5% | Reject |
| 28 | Manual `gc()` as memory control | 2 | 2 | 4 | 20% | Reject |
| 29 | Explicit worker/BLAS thread budget | 4 | 3 | 3 | 75% | Keep as requested/recorded; enforcement backend-dependent |
| 30 | Sequential/parallel/worker-count benchmark gates | 5 | 3 | 1 | 98% | Keep |

## 14. Verification gates

Execution correctness is tested against the dense mathematical oracle.

### Test families

- **Contract tests:** result type/completeness, source capabilities, v0.1
  worker rejection before source access, policy hash purity, failure receipts,
  and reporter isolation.
- **Oracle tests:** dense reference versus streamed/tiled, direct versus
  materialized query, in-memory versus block-backed storage, and ordinary
  versus shard-backed sources when Phase E2 exists.
- **Metamorphic tests:** feature permutation, undirected endpoint swap,
  legal block/tile changes, completion-order permutations, conservative frame
  sums, and coherent/configuration decomposition.
- **Adversarial tests:** a deliberately slow first task with fast successors,
  tiny reorder budgets, reporter exceptions, output larger than RAM budget,
  rank-deficient extractors, sparse/dense frame extremes, interruption during
  write, worker recycle, and stale descriptor attempts.
- **Regression tests:** every observed allocation blow-up, shared-handle reuse
  failure, accidental remote execution for `workers = 1`, or partial-result
  misclassification receives a minimal permanent reproducer.
- **Performance guardrails:** fixed representative \((p,q,m)\) regimes record
  wall time, aggregate private RSS, mapped bytes, serialization, contraction
  temporary size, and end-of-run memory return.

### Semantic parity

- sequential dense = sequential streamed;
- direct query = materialize then query;
- one worker = multiple workers;
- every legal block size agrees;
- every legal measurement-row and geometry-coordinate tiling agrees;
- task completion permutations agree;
- resume = uninterrupted execution;
- file-backed worker-read = coordinator-read;
- shard-backed response/relation sources = ordinary-source reference;
- conservative and decomposition laws remain exact to tolerance.

### Resource behavior

- peak resident memory stays within a declared tolerance of the planner estimate;
- increasing workers cannot multiply full-source memory silently;
- workers are closed on success, error, interrupt, and cancellation;
- worker-owned file handles and temporary files are released;
- shared owners outlive all attachments, stale handles are never reused, and
  repeated interrupted executions leave no live segments;
- adversarial completion order never exceeds `max_inflight` or
  `max_reorder_bytes`, and reduced task payloads become unreachable promptly;
- in-memory and block-backed complete geometries agree;
- reporter exceptions leave scientific outputs unchanged and appear only in
  receipt observation details;
- no global future plan or RNG state changes;
- no nested process/thread oversubscription under the default policy.

### Performance gates

- the vectorized one-worker kernel is the baseline, not brute-force rMVPA;
- parallel execution is shipped only if it wins on realistic file-backed inputs after startup/serialization costs;
- shard staging is benchmarked on aggregate worker/coordinator RSS,
  serialization bytes, hidden materialization, cleanup, and wall time—not just
  coordinator RSS;
- a one-worker request must always use the direct local loop;
- direct-query mode must demonstrate the expected memory reduction;
- full geometry at large \(m,h\) must demonstrate contraction-temporary bounds
  under row/coordinate tiling;
- performance claims require numerical parity and identical estimands.

## 15. Staged implementation

### Phase E0: sequential compiler

Finish the semantic contract, then implement query-aware deterministic feature
streaming, row/coordinate-tiled sparse contraction, in-memory and block-backed
complete geometry, and distinct `effect_geometry`/`effect_view` results. Add
execution receipts and conservative memory planning validated against measured
peak memory. No parallel dependency; version 0.1 rejects `workers > 1`.

This phase may already deliver most of the speedup.

### Phase E1: in-process scheduling simulator and task protocol

Not started. A partial simulator was written and tested in isolation, never
wired into the compiler, and deleted on 2026-08-16.

Factor the block kernel into pure task inputs/outputs. Simulate adversarial
completion order, bounded reordering, backpressure, failure receipts, reporter
failure, and cleanup without creating worker processes. Add source capabilities
and worker-readable file-backed sources.

### Phase E2: optional `shard` executor and shared source

Only if end-to-end benchmarks justify it, add a `shard` adapter for immutable response or
relation staging and coarse feature-block dispatch. Reuse `shard` supervision,
recycling, buffers, and telemetry; keep geometry reduction and result assembly
in the one `crossform` kernel. Do not expose a global plan or a second
`shard`-specific analysis path.

#### Phase E2 admission result, 2026-08-12

The first admission benchmark did **not** justify a public or default `shard`
executor. It compared the installed `crossform` 0.0.0.9000 sequential
compiler with installed `shard` 0.2.0 using exactly the same pure feature task,
canonical ordered reducer, simulated relation, frame, and pairing. It measured
the complete cold path, including shared staging and pool creation, in isolated
processes and polled aggregate coordinator-plus-descendant RSS.

All five modes agreed numerically for total and coherent geometry. Normal and
injected-failure cleanup both stopped owned pools and removed owned mmap
backings. Runtime failed the predeclared 1.1x cold-speedup gate:

| Mode | Seconds | Speed relative to sequential |
|---|---:|---:|
| Sequential response | 0.419 | 1.00x |
| Shard response, cold | 1.627 | 0.26x |
| Shard response, reused | 0.780 | 0.54x |
| Shard relation, cold | 1.529 | 0.27x |
| Shard relation, reused | 0.528 | 0.79x |

This result preserves the architectural value of `shard` without admitting an
unearned integration. Shared relation staging reduced mapped source size from
113.2 MB to 12.6 MB, and aggregate sampled process-tree RSS was lower than the
sequential run in this regime, but neither cold nor reused execution was
faster. Version 0.1 therefore remains sequential-only and has no `shard`
dependency or executor code path. A future reconsideration requires a new
benchmark on materially heavier per-feature work, more blocks, or a genuinely
reused caller-owned execution scope; it may not weaken the cold-path gate or
duplicate the scientific kernel.

The reproducible harness is `benchmarks/run-shard-admission.R`; the evidence
record is `benchmark-results/shard-admission.rds` with its human-readable CSV
summary.

### Phase E3: checkpointing

Add disk preflight, restrictive permissions, schema/checksums, strong source
receipts, same-filesystem atomic writes, and exact resume-plan validation.

### Phase E4: external adapters, if contract-compatible

Only then consider `future`, `mirai`, or cluster schedulers through the fixed
task protocol. An adapter that inherits ambient process-global state or cannot
prove ownership does not satisfy the executor contract and is not promised.

#### Phase E4 evaluation result, 2026-08-12

No external executor is admitted for version 0.1.

- **`future` is rejected as a core-owned executor.** Its documented package
  contract asks package authors not to choose or mutate the active plan; a
  backend is selected through process-global plan state. Inheriting that state
  violates crossform's explicit ownership and resource-budget contract, while
  locally replacing it would make crossform the kind of plan-mutating package
  this architecture is intended to avoid. A caller may orchestrate independent
  crossform jobs with `future`, but the core will not treat the ambient plan
  as its executor.
- **`mirai` remains technically compatible but deferred.** Named compute
  profiles and scoped daemon selection can provide explicit isolation, and
  `daemons(0, .compute = ...)` provides a clear owned-cleanup operation.
  However, the immediately preceding `shard` benchmark showed that two-process
  dispatch could not beat the sequential feature compiler on the admitted
  workload. Building a second process adapter before a heavier irreducible
  workload exists would add dependencies and lifecycle code without an earned
  product benefit.
- **Batch schedulers are outside the local 0.1 compiler.** They are useful for
  participant- or plan-level orchestration, where jobs already have disjoint
  outputs and coarse runtimes. They are not an appropriate replacement for
  feature-block execution inside one geometry, and no configured batch backend
  is part of the current certification environment.

Reconsider `mirai` or a batch adapter only when a real workload cannot meet its
resource target sequentially and the adapter can be benchmarked against the
same feature task, ordered reducer, source descriptors, receipt identity, and
cleanup laws. The public core remains executor-neutral rather than advertising
unverified adapters.

## 16. Learned local metric as a compiler lowering

Added 2026-08-17; implemented (B2) 2026-08-20 and completed (B3) 2026-08-20 —
the lowering, both identity hashes, the capability refusals, and the sugar
route below are in the tree, exercised by
`tests/testthat/test-geometry-learned-metric.R`. `effect_crossnobis_plan`,
`.plan_learned_crossnobis()`, `.execute_learned_crossnobis()` and their
identity, memory and receipt helpers are **retired (B3, 2026-08-20)**; the
numerical equality with that driver is pinned as recorded golden values in
`tests/testthat/fixtures/learned-crossnobis-golden.rds`. The paragraph below
describes the fork as it stood before B2/B3 and is kept as the record of what
was collapsed. The IR is not the fork: `effect_evidence_task` is already
sole, built by `plan_geometry()` (`R/geometry-plan.R:324`) and by
`plan_crossnobis()` (`R/crossnobis.R:344`) alike. The fork is at the
**executor**. Fixed-metric crossnobis routes `R/crossnobis-driver.R:248` into
`.run_geometry_compiler()`; learned crossnobis bypasses `R/compiler.R` through a
second plan class (`effect_crossnobis_plan`), a second identity pair
(`.crossnobis_scientific_plan_id`, `.crossnobis_plan_signature`), a second
memory planner, a second planned receipt, and a second driver
(`.execute_learned_crossnobis`, `R/crossnobis-driver.R:40`) — ~370 lines of plan
and driver in front of a 118-line kernel. This section specifies how that path
becomes an ordinary lowering so both can be retired.

### 16.1 Metric schedule kind `learned_local_before_frame`

`.geometry_metric_schedule()` (`R/geometry-plan.R:3`) emits two kinds today; it
gains a third when `plan_geometry(metric = )` receives an `effect_metric_recipe`
rather than an `effect_neural_metric`:

```
kind = "learned_local_before_frame"; frame_composition = "sqrt_weight_congruence"
feature_additive = FALSE; support_dense = TRUE   # a per-support solve is not additive
materialization = "on_demand_local"; scope = "support_local"
lowering = "derive_then_support_streamed_pair_contraction"
metric_signature = NULL; schedule = <effect_frozen_metric_schedule>
```

`$schedule` is exactly today's object from `compile_metric_schedule()`
(`R/metric-learning.R:503`): recipe, canonical `residual_pair_statistics()`,
support index, evaluation pairing, `metric_training_policy()`, one training
record per edge, capability block. Nothing in `R/metric-learning.R` changes; it
stops being reachable only through `plan_crossnobis()`. The `derive_then_*`
lowering string already exists (`R/metric.R:522`) and is unreachable today
because no plan admits a recipe.

`.validate_geometry_metric_schedule()` (`R/metric.R:540`) gains a third branch
calling `.validate_frozen_metric_schedule()` plus the cross-object checks
`.validate_crossnobis_plan()` makes today (`R/crossnobis.R:387-396`):
`schedule$pairing` identical to the plan pairing,
`schedule$support_index$signature` identical to
`frame$support_index$signature`, recipe domain identical to relation domain.
The first is load-bearing at runtime — `R/kernel.R:683-694` re-checks that the
schedule's records line up with `task$ordered_edges`, and the kernel loops over
`schedule$pairing` rows, not the ordered edges, to contract each declared edge
once (self-adjoint shortcut, `R/kernel.R:735-737`). The new lowering string
joins `valid_lowering` in `.validate_geometry_plan()` (`R/geometry-plan.R:382`);
the `lowering` and `materialization` enumerations widen.

Three entry-level consequences, stated rather than discovered:

- **`plan_geometry()` must accept an `effect_relation_fit`.** A recipe needs a
  residual channel, so `.validate_geometry_plan_inputs()`
  (`R/geometry-plan.R:147`) gains a fit-unwrapping preamble calling
  `.require_relation_fit_capability(x, "learned_metric_input")` *before* shape
  validation — the diagnosis order `R/crossnobis.R:309-314` documents. A bare
  relation plus a recipe gets that same capability refusal.
- **A learned schedule requires a support-indexed frame.**
  `.residual_statistics_support_index()` (`R/residual-statistics.R:3-13`)
  refuses a frame without one, so `plan_geometry(voxelwise(), metric = recipe)`
  refuses. Inherited behaviour, new entry point.
- **The plan reads.** `plan_geometry()` promises validation "without reading
  relation blocks" (`R/geometry-plan.R:199-203`). A learned schedule accumulates
  residual sufficient statistics at plan time in one streamed pass, as
  `R/crossnobis.R:340` does. Amend the promise for this kind rather than
  deferring accumulation: deferring re-accumulates per contrast and loses the
  plan reuse the frozen schedule exists to provide. The two refusals that fire
  before any residual read today — `.preflight_metric_training()`'s partition
  shortage and the workspace-budget overflow — must keep firing before it.

### 16.2 The lowering `support_streamed_scheduled_metric_query_contraction`

`.compile_geometry_execution_plan()` (`R/compiler.R:128`) gains one branch and
emits an ordinary `effect_geometry_execution_plan`. The kernel already exists:
`.support_streamed_scheduled_crossnobis()` (`R/kernel.R:666`), today called from
one site, `R/crossnobis-driver.R:85`. Its call site becomes
`.execute_node_block()` (`R/execution-driver.R:311`), beside the fixed-metric
`.support_streamed_metric_contraction()` (`R/kernel.R:422`, called at
`R/execution-driver.R:316`). Branch points:

- `.execution_support_streamed()` (`R/execution-driver.R:307`) matches
  `^support_streamed_metric_` today; it becomes a three-way selector over
  `plan$lowering` — fixed kernel, scheduled kernel, or
  `.streamed_effect_form_contraction()`.
- `.compile_geometry_execution_plan()` computes `learned <-
  identical(plan$metric_schedule$kind, "learned_local_before_frame")` beside its
  existing `explicit_metric`/`support_streamed` flags, emitting the lowering
  string above when `learned && query_fused`.
- `.geometry_kernel_version()` (`R/compiler.R:113`) gains a leading
  `if (isTRUE(learned)) "support-streamed-scheduled-metric-v1"`, preserving the
  kernel id asserted at `test-crossnobis-learned.R:58`.
- `.validate_geometry_execution_plan()` (`R/compiler.R:299`) re-derives both
  strings through the same widened helpers, so validator and compiler cannot
  drift.

**Stage sequence** is the geometry executor's, unchanged — which is the point.
`.execute_geometry_compute()` (`R/execution-driver.R:535`) opens the source
session, records `source_admission`, runs `.execute_node_block()` and records
`feature_tasks`, then skips coherent and marginal blocks because the scheduled
requirements demand neither. The private `support_tasks` stage label disappears;
the fixed support-streamed route already reports its node loop as
`feature_tasks`.

**Tiling.** The kernel streams one support at a time and holds at most two local
covariance matrices live: `feature_block = max(support_size)`, `row_tile = 1L`,
`coordinate_tile = 1L` (`.support_metric_memory_plan()`'s `min(output_width,
64L)` already yields this at `output_width == 1`), `task_count =
plan$measurements` — the branch `R/compiler.R:284` already takes.

**Memory.** `.support_metric_memory_plan()` (`R/memory-plan.R:495`) refuses a
non-fixed schedule at line 498. Its learned branch must charge what
`.learned_crossnobis_memory_plan()` charges and the fixed branch does not:
retained pair coordinates plus one cross-product vector per partition, the
support-index CSR membership and canonical symmetric pair-pattern CSC slots, and
`2*k*k` contraction bytes for the live covariance pair. Copying the fixed
branch's `object.size(schedule$metric$value)` would silently drop the whole
residual-statistics payload, the dominant resident term at 52k features. The
`fits_budget` refusal keeps its current message shape.

**Admissible reads.** The kernel evaluates a *rank-one signed* contrast only: it
forms `c %*% B` per endpoint and contracts one scalar per edge. No coherent
component, no endpoint marginals, no packed form. The lowering admits
`component = "total"` with `storage = "memory"` and refuses `"coherent"`,
`"configuration"`, `"contrast"`, full materialization, block storage and
rectangular plans — one capability `"scheduled_metric_component"` in namespace
`"geometry_views"`, remedy naming `crossnobis()`. New refusal surfaces, not
regressions: no route to them exists today.

**The contrast hint.** The kernel needs the signed vector `c`, not packed
`svec(cc^T)`. The execution plan carries `signed_query` for `component =
"total"` when the schedule is learned, guarded by the outer-product equality
check `R/compiler.R:165-178` already applies to contrast plans, so the hint
cannot alter the estimand. It enters `.geometry_execution_signature()` (already
listed there) but is excluded from `.geometry_view_scientific_id()` for this
component: `total` of `cc^T` is one estimand whether or not the executor was
handed `c`, and folding the sign in would make `c` and `-c` name different
estimands.

**Return shape.** The kernel returns `values` as a bare numeric vector plus
`metric_receipts` and `endpoints_read`; `.execute_node_block()` must wrap it as
`list(value = matrix(values, ncol = 1L), diagnostics = , metric_receipts = )`
so `.execution_component_values()` and `.execution_geometry_result()` consume it
unchanged.

### 16.3 Identity

**`scientific_plan_id`** — `.geometry_plan_scientific_id()`
(`R/geometry-plan.R:65`), prefix `geometry-sha256:`, digests exactly:
(1) `schema_version = 1L`; (2) `evidence_task` = `task$task_id`; (3) `frame` =
`.additive_frame_signature(frame)`; (4) `metric_schedule` =
`metric_schedule$signature`; (5) `component` (`"full"` for a plan);
(6) `signed_query` (`NULL` for a plan).

The metric is part of the estimand through (4). The learned schedule's semantic
digest is its own fields minus `$schedule` and `$signature`, plus the frozen
schedule signature; `.metric_schedule_signature()` (`R/metric-learning.R:427`)
in turn digests `role`, `recipe_specification`, `recipe$signature` (estimator
kind, shrinkage, both variance floors, spectral floor, randomness, seed, bound
domain), `statistics$signature`, `support_index$signature`,
`.metric_pairing_identity(pairing)`, `training_policy$signature` (kind,
`includes_evaluation_residuals`, assumption, justification text), the per-edge
`training_signature` vector (training-partition assignment, atomic statistic
signatures, source revisions, residual revisions), and the capability block.

Nothing in `.crossnobis_scientific_plan_id()` is lost. Its `relation_fit` field
re-enters through `statistics$signature`, which carries the fit signature; its
`pairing = .metric_pairing_identity(over)` re-enters through (4) — which matters,
because the evidence-task semantic (`.evidence_task_general_semantic()`,
`R/evidence-task.R`) digests the ordered partition products but **not** the
pairing's declared `independence`, `estimate`,
`self_pairs`, `directed` or `generalizes_over`. Those declarations stay
estimand-bearing for a learned plan; the pre-existing asymmetry on fixed-metric
geometry plans is out of scope here.

For a view, `.geometry_view_scientific_id()` (`R/geometry-plan.R:102`) digests
`schema_version`, `role = "geometry_view"`, `parent` (the plan estimand id),
`component`, `.query_identity_semantic(query)`, `signed_query` — `NULL` for a
learned total, per §16.2.

**Execution `signature`** — `.geometry_execution_signature()` (`R/compiler.R:92`)
digests exactly: `schema_version`, `parent` (plan `$signature`),
`scientific_plan_id`, `task` (`task$task_id`), `storage`, `storage_path`,
`component`, `signed_query`, `requirements`, `output_width`, `feature_block`,
`row_tile`, `coordinate_tile`, `memory` (unclassed `memory_plan`), `lowering`,
`kernel_version`. The schedule reaches it transitively through
`scientific_plan_id` and directly through `lowering`; tiles and kernel version
are explicit.

`residual_workspace_bytes` — a cache-capacity knob that never changes the
canonical numerical tile shape, defaulting to `compute$workspace_bytes` or
512 MiB — becomes a `plan_geometry()` argument admitted only with a recipe,
recorded in an `$execution_hints` slot digested by `.geometry_plan_signature()`
(`R/geometry-plan.R:133`) as a fourth component beside `scientific_plan_id`,
`compute` and `dense_payload_bytes`. It must not reach
`.geometry_plan_scientific_id()`. Folding it into `compute_policy()` would be
tidier but would move every existing plan signature and invalidate every
recorded certification digest.

### 16.4 Capabilities, not class checks

`.metric_schedule_capabilities()` (`R/metric-learning.R:414`) already sets
`calibration_requires_metric_uncertainty = TRUE` for every non-identity recipe.
That flag is the whole refusal and must reach the sampling layer without a class
check.

Do **not** add a `$capabilities` field to `effect_metric_schedule`: its semantic
digest is `unclass(x[!names(x) %in% c("metric", "signature")])`
(`R/metric.R:579`), so a new stored field moves every existing geometry
`scientific_plan_id`. Derive it instead. A predicate
`.metric_schedule_requires_metric_uncertainty(schedule)` returns
`isTRUE(schedule$schedule$capabilities$calibration_requires_metric_uncertainty)`
for the learned kind, `isTRUE(schedule$metric$capabilities$learned_frozen)` for
a fixed metric materialized from a learned handle, `FALSE` otherwise.

`.sampling_evidence_descriptor()` (`R/evidence-sampling.R:100`) then drops its
`effect_crossnobis_plan` branch (`:116-131`) and sets `metric_status <-
if (.metric_schedule_requires_metric_uncertainty(x$metric_schedule)) "learned"
else "fixed"` in the one remaining geometry branch. The refusal at
`R/evidence-sampling-product.R:569-584` (capability `fixed_metric_sampling_law`,
reason `learned_metric_law_not_admitted`) then fires for learned geometry plans
unchanged, and the second refusal at `:585-599` (`crossnobis_plan_not_routed`)
is deleted with the class it tests for. `sampling_covariance()` and
`plan_sampling()` need no other change: both reach the descriptor before
compiling anything.

### 16.5 `crossnobis()` as a view, `plan_crossnobis()` as sugar

`crossnobis()` (`R/crossnobis-driver.R:224`) keeps one argument type,
`effect_geometry_plan`, and validates rather than dispatches. It requires either
a `fixed_metric_before_frame` schedule whose metric carries
`$provenance$metric_role == "noise_precision"`, or a `learned_local_before_frame`
schedule whose recipe is not `identity` — the two claims
`.crossnobis_plan_metric()` (`R/crossnobis.R:68`) and `R/crossnobis.R:319-324`
make today, under the one capability `"declared_noise_metric"`. It refuses the
implicit identity metric, a fixed metric not built by `noise_precision()`, an
`identity_metric()` recipe, and a pairing that is not independent,
cross-partition and self-product-free. `missing(weights)` is now checked on both
routes; the learned route skips it today (`R/crossnobis-driver.R:225`).

One behaviour moves: `plan_crossnobis()` enforces the pairing contract at plan
time (`R/crossnobis.R:316`), `plan_geometry()` does not, so the pairing refusal
fires at view time — where the fixed-metric route already puts it
(`R/crossnobis-driver.R:241`). `.preflight_metric_training()` still refuses a
training-partition shortage at plan time, so a learned plan still cannot be
built on a pairing unusable under `exclude_evaluation`.

`plan_crossnobis(x, at, over, metric = shrinkage_precision(), training, compute,
residual_workspace_bytes)` becomes a thin wrapper over `plan_geometry(...)` with
the same defaults. Keep it exported: it names the intent,
`exemplars/haxby2001/05-crossnobis-uncertainty.R:178` and
`vignettes/from-rmvpa.Rmd:76,377` use it, and retiring the verb buys nothing.
What is retired is the *class*.

The view's public fields are unchanged — `$values`, `$contrast`, `$estimand`,
`$metric`, `$pairing`, `$index`, `$receipt` — with `$metric` still the frozen
schedule signature on the learned route and the neural metric signature on the
fixed route. Signed negative finite estimates are still retained.

`.execution_metadata()` (`R/execution-driver.R:426`) gains a learned branch
carrying what `R/crossnobis-driver.R:106-141` carries today: recipe signature
and kind, training policy, per-edge records (evaluation endpoints, training
partitions, training signature), `local_metric_storage =
"none_derived_on_demand"`, `retained_factor_table = FALSE`,
`calibration_requires_metric_uncertainty`, and the kernel's `metric_receipts`.
Kernel diagnostics (`pair_atoms_materialized`, `pair_frame_materialized`,
`metric_factor_table_retained`, `metric_handles_derived`,
`max_local_covariance_bytes`) move from a flat `$metadata$diagnostics` to
`$metadata$diagnostics$total`, where every other geometry route puts them.
`$metadata$source_session` keeps its named per-partition `read_count`, which is
what proves the executor reads evaluation endpoints only.

### 16.6 Migration

| Symbol / class / field | Disposition |
|---|---|
| `plan_crossnobis()`, `crossnobis()` (exported) | **Kept**: sugar over `plan_geometry()`; validating view over `effect_geometry_plan` only |
| `noise_precision()`, `shrinkage_precision()`, `diagonal_precision()`, `identity_metric()`, `metric_training_policy()`, `metric_capabilities()`, `residual_pair_statistics()` | **Unchanged** |
| `effect_crossnobis_plan` (class) | **Retired** — done (B3, 2026-08-20); `plan_crossnobis()` returns `effect_geometry_plan` |
| `print.`/`format.effect_crossnobis_plan` (`R/format-results.R:126`) | **Retired** — done (B3, 2026-08-20), both `S3method()` lines dropped from `NAMESPACE`; `print.effect_geometry_plan` carries the recipe-kind and training-policy lines |
| `.validate_crossnobis_plan`, `.crossnobis_scientific_plan_id`, `.crossnobis_plan_signature`, `.learned_crossnobis_memory_plan` | **Retired** — done (B3, 2026-08-20); geometry equivalents plus the learned branch in `.support_metric_memory_plan()` |
| `.planned_crossnobis_receipt`, `.execute_learned_crossnobis` | **Retired** — done (B3, 2026-08-20); `.planned_execution_receipt()` and `.execute_geometry_plan()` cover both |
| `R/crossnobis.R`, `R/crossnobis-driver.R` | Done (B3, 2026-08-20): `crossnobis.R` (251 lines) keeps `noise_precision()`, the metric- and pairing-role checks, and the sugar; `crossnobis-driver.R` (136 lines, 71 of them roxygen) keeps `crossnobis()` and may still fold into `R/views.R`. `design/architecture.md` updated; both files keep their layer-3 / layer-4 entries in `test-architecture.R`, so the layer map is unchanged |
| `.support_streamed_scheduled_crossnobis()` (`R/kernel.R:666`) | **Kept unchanged**; call site moves to `.execute_node_block()` |
| `compile_metric_schedule()`, `.metric_schedule_provider()`, `materialize_metric()`, `effect_frozen_metric_schedule` | **Unchanged**; reached through the schedule kind |
| receipt `$kernel_version = "support-streamed-scheduled-metric-v1"` | **Kept** |
| receipt `$scientific_plan_id` prefix | **Renamed** `crossnobis-sha256:` → `geometry-sha256:`; every learned plan and view id changes once |
| receipt `$task_partition_id = "ascending-supports-one-live-node"` | **Renamed** to the geometry form (`R/execution-driver.R:42`) |
| execution stage `"support_tasks"` | **Renamed** `"feature_tasks"` |
| `metadata$execution_plan$materialization = "direct_crossnobis_contrast"` | **Renamed** `"direct_total"` |
| `.sampling_evidence_descriptor()` crossnobis branch; `crossnobis_plan_not_routed` refusal | **Retired** — done (B3, 2026-08-20); `.metric_schedule_requires_metric_uncertainty()` replaces both, and the surviving `learned_metric_law_not_admitted` refusal is covered by `test-geometry-learned-metric.R` |

Must pass **unchanged**: all of `test-crossnobis-known.R` (the fixed route is
untouched); all of `test-metric-learning.R` (recipes, policies, schedule
compilation, providers, oracle agreement, leakage refusals — none of it names
the plan class); and in `test-crossnobis-learned.R`, oracle agreement (`:49-56`,
`:122-135`), receipt status and kernel version (`:57-59`), lowering string
(`:60`), zero-residual-read metric receipts (`:69-73`), per-partition source
read counts (`:76-80`), identity separation (`:107-121`), plan-time refusal
ordering (`:142-191`), recipe hyperparameters (`:193-200`), the
scientific-vs-execution identity split (`:202-226`), and the
no-revalidation-per-node bound (`:228-257`).

Needs **rewriting** (all done, B2 and B3, 2026-08-20):
`test-crossnobis-learned.R:54` (class assertion →
`"effect_geometry_plan"`), `:62-68` (diagnostics path moves under
`$diagnostics$total`), `:107-109` (same assertions, read off the geometry plan);
`test-print-methods.R:677` (printed-class inventory loses
`effect_crossnobis_plan`); `test-architecture.R:37,50` (layer map — no change
was needed in the end, both files survive in their existing layers). The B2
old-versus-new equality test is now a golden-value pin against
`fixtures/learned-crossnobis-golden.rds`, recorded from the retired driver's
last run before deletion. New tests
are owed for the component refusals of §16.2, the sampling refusal reached by
capability on a learned geometry plan, and
`identical(plan_crossnobis(...)$scientific_plan_id, plan_geometry(..., metric =
recipe)$scientific_plan_id)`.

### 16.7 Open risks

1. **Per-tile training leakage.** Training partitions are assigned once per
   evaluation edge, before tiling; `feature_block = max(support_size)` makes a
   tile a support, so a support-dependent training set would be a different
   estimator per node. *Resolution:* the record is frozen in
   `compile_metric_schedule()` and read-only in the kernel; assert in the
   learned validator that `record$training_partitions` is independent of
   `support_index`, and keep `.metric_schedule_provider()` the only path to a
   handle.
2. **Evaluation-edge cross-fitting.** `exclude_evaluation` trains on partitions
   disjoint from both endpoints; `all_partitions_residual_orthogonality` does
   not, and buys admission with a recorded justification, not proof. Moving to a
   geometry plan widens the surface, since `plan_geometry()` accepts pairings
   `plan_crossnobis()` refused. *Resolution:* keep
   `.preflight_metric_training()` at plan time and the independence refusal at
   view time, so no crossnobis reading can come from a dependent pairing.
3. **Memory of per-feature precision blocks.** The 52k-feature gate (59.98 s /
   903 MB) is dominated by retained residual pair statistics, not local
   precision — handles are derived and discarded per node, no factor table is
   kept. The risk is the learned memory branch inheriting the fixed branch's
   metric-object accounting. *Resolution:* port the payload formulas from
   `.learned_crossnobis_memory_plan()` verbatim and assert a learned plan's
   planned workspace exceeds a fixed-metric plan's on the same frame.
4. **Identity churn.** Every learned scientific plan id and view id changes
   prefix and content. *Resolution:* take the one-time break, bump the schedule
   schema version, re-record certification artifacts in the same commit; no
   compatibility shim.
5. **Conditional plan-time reads.** One entry point now sometimes reads and
   sometimes does not — the kind of conditional promise this design elsewhere
   refuses. *Resolution:* state it in `plan_geometry()`'s `@section Structure`,
   print it on the plan, and record the accumulation in `$execution_hints` so it
   is visible in the plan signature, not only in the receipt.

### 16.8 Metric schedule kind `whitened_metric_before_frame` (D6)

`.geometry_metric_schedule()` emits a fourth kind when
`plan_geometry(metric = <fixed metric>, composition = "whitened")` is compiled.
It is the mirror image of §16.1: where the learned kind is the *least* additive
schedule, this one is the identity lowering wearing a metric's name.

```
kind = "whitened_metric_before_frame"; composition = "whitened"
root = "symmetric_psd_root"; frame_composition = "sqrt_weight_congruence"
feature_additive = TRUE; support_dense = FALSE   # additive in whitened coordinates
materialization = "whitened_effect_coordinates"; scope = "domain_operator"
lowering = "additive_contraction"; metric_signature = <metric>
```

The estimand is `K_x = Q^(1/2) D(w_x) Q^(1/2)`, which is dense on the whole
domain even though `D(w_x)` is supported on one node. It therefore cannot go
through `.compose_frame_metric()`, whose contract requires the metric's support
to equal the node's — the reason
`design/conservative-geometry-contract.md` §5.2.2 budgeted D6 for a
relation-level transform rather than a metric-schedule variant. The seam is
`plan_geometry()`: it whitens the effect coordinates once, `B~ = B Q^(1/2)`,
compiles the task on `B~`, and lets the ordinary implicit-identity pipeline run,
because `B~ D(w_x) B~ᵀ = B Q^(1/2) D(w_x) Q^(1/2) Bᵀ` exactly. Nothing in
`R/compiler.R`, `R/execution-driver.R`, or `R/kernel.R` learns a new branch;
`.metric_additive_frame()` refuses the kind outright, since folding `diag(Q)`
into the frame would compose the metric a second time.

Three consequences, stated rather than discovered:

- **`composition` is estimand-bearing and the root is part of it.**
  `Σ_x R D(w_x) Rᵀ = Q` holds for any `R` with `RRᵀ = Q`, so a conservation
  certificate cannot tell two roots apart while the node values differ by
  double-digit percentages (contract §5.2.1). Both `composition` and `root`
  enter the schedule's semantic digest and so the `$scientific_plan_id`. They
  are **absent** from the other three kinds' field lists rather than present
  with a `"native"` value, because the digest is the whole field list and a new
  field would have moved every geometry plan identity in existence.
- **It is the second conditional plan-time read** (§16.7 risk 5, now two
  instances). The source is read one feature block at a time and accumulated,
  so the input stays bounded, but the output cannot: a global congruence has no
  blockwise output. The resident cost — one whitened effect-by-feature matrix
  per partition, plus the symmetric root — is recorded in `$execution_hints`,
  enforced against a declared `compute_policy(workspace_bytes = )` before the
  first read, and reaches `$signature` without touching
  `$scientific_plan_id`, because a cost is not an estimand.
- **Admitted for fixed positive-definite domain-wide metrics only.** A learned
  recipe gets a capability refusal (`whitened_metric_composition`): a
  per-support operator has no single global root, so there is no one set of
  whitened coordinates every node could share. A support-local fixed metric is
  refused for the same reason at smaller scale, and a positive-semidefinite but
  singular metric is refused by naming the offending eigenvalue rather than
  truncated.

## 17. Identity schema consolidation (B6)

Two identity schemas used to coexist in `R/evidence-task.R`. Every evidence
task carried an `$identity_schema` field naming which of them had produced its
`$task_id`:

- `evidence-pairing-v1` — the boundary-typed semantic built by
  `.evidence_task_general_semantic()`: four identified spaces, the ordered
  partition products and their expansion, **both** boundary signatures, the
  stage plan signature, the materialization signature.
- `effect-form-v1` — a flatter legacy semantic (`.effect_task_semantic()`) used
  only by the bridged effect-form route: relation ids, the two effect spaces,
  the ordered edges, the bridge signature, the three operations, the query
  identity. It could not mention either boundary signature, because it predated
  the boundary-typed IR. It existed so that ids recorded before that IR stayed
  byte-stable, and it borrowed its name from the scientific contract
  `design/effect-form-contract.md`, which is a different object and is
  untouched by this change.

The compatibility window is now closed. `evidence-pairing-v1` is the only
identity schema; `.effect_task_semantic()` is deleted; `.new_evidence_task()`
no longer takes a schema or a legacy compatibility semantic and stamps the one
schema unconditionally, so a forged `$identity_schema` fails the rebuild in
`.validate_evidence_task()`. Two gates that were phrased as schema checks are
now phrased structurally: `.as_compiled_effect_task()` admits an effect-form
task with an open experimental boundary and a **bridge**-closed neural boundary
(the closure kind matters — a query-closed effect form has no `$bridge` to
project), and `.validate_compiled_effect_task()` recomputes the recorded id
from the same boundary-typed parts its evidence task was named from, via the
shared `.effect_form_evidence_parts()`.

### The id migration

Consolidation renames every task on the bridged effect-form route, and
everything derived from those names. It renames nothing else.
`tests/testthat/test-identity-schema.R` and
`tests/testthat/fixtures/identity-schema.rds` record both sides: eleven
representative identities, taken under the two-schema tree (`$before`) and the
consolidated tree (`$after`).

| identity | moved? |
| --- | --- |
| `effect_form_complete`, `effect_form_complete_adapter` | yes |
| `effect_form_pair_query`, `effect_form_pair_query_base` | yes |
| `effect_form_physical_query`, `effect_form_reversed` | yes |
| `effect_task_plan_id` | yes |
| `geometry_plan_fixed_metric`, `geometry_plan_learned_metric` | yes |
| `measurement_form` | **no** |
| `effect_form_neural_query` | **no** |

The two that do not move are the routes that were already native to
`evidence-pairing-v1`: a measurement form (closed experimental boundary, open
neural boundary) and an effect form whose neural boundary closes with a fixed
query rather than a bridge. Retiring the duplicate naming rule did not rename
the tasks that never used it.

Because `scientific_plan_id` digests `task$task_id`, every `geometry-sha256:`
and `crossnobis-sha256:` id recorded in `inst/extdata/certification/` is stale
after this change, as it already was after the rest of this program's `R/`
edits. Re-recording those artifacts is ticket B8. No test compares a
certification artifact's recorded id to a freshly computed one — the artifact
tests match the id *format* and re-derive verdicts from recorded measurements —
so the staleness is a provenance debt, not a failing gate.

### Materialization kinds

The materialization enum lost its third arm in the same pass. `scalar_field`
was legislated for the `closed/closed` boundary pair — both boundaries closed
by a fixed query, so the task materializes one scalar per frame node — and was
never constructed anywhere: the query-fused geometry route reaches those
numbers by keeping the experimental boundary open and carrying the query in the
materialization projection, which is why every such task is spelled as a
`query_only` `effect_form`. Two kinds remain, `effect_form` (`open/closed`) and
`measurement_form` (`closed/open`), and `.validate_evidence_boundary_combination()`
now refuses `closed/closed` for both. The comment at the removal site in
`R/evidence-task.R` lists what a future scalar-field materialization must
re-introduce: the enum value, the necessarily-`query_only` invariant, the
`closed/closed` arm, a reversal rule for a task whose experimental *and* neural
queries both transpose, and an executor that admits the kind.

## Final recommendation

The top three execution improvements are:

1. **Scope the collapse exactly.** Additive diagonal-frame searchlights with
   fixed bilinear queries are sparse matrix rows, not scheduled fits. Learned,
   locally estimated, and nonlinear methods are not covered by that theorem.
2. **Close the contraction-memory hole.** Plan and measure dense contraction
   temporaries, R copies, serialization, reordering, checkpoints, and output;
   tile rows and geometry coordinates or use block-backed complete geometry.
3. **Make scheduling bounded and numerically honest.** Enforce inflight/reorder
   backpressure, reduce canonically, promise tolerance-qualified agreement, and
   separate full geometry from direct-query results.

The package should ship its sequential compiler before any parallel backend.
After those contracts are executable, `shard` should be the preferred first parallel experiment, not because the
package needs a new engine, but because `shard` already solves the remaining R
systems problem: several supervised processes reading one large immutable
object without several full copies. The admission gate is end-to-end numerical,
memory, cleanup, and runtime evidence against the sequential compiler.
