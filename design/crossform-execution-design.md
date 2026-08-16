# `crossform` execution model

Status: architecture proposal; review corrections incorporated, not implementation
Date: 2026-08-12
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
