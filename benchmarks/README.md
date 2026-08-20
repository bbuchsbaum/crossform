# Benchmark and certification harness

**Re-certifying after a source edit:** `benchmarks/RECERTIFY.md` is the exact
command sequence — temp library install, runners in order, promotion, and the
verifying test run — with the traps that cost time to rediscover. This file
names the map-scale verbs those gates certify and explains what each gate
asserts; that one explains how to make the recorded evidence bind again.

## Provenance binding

A recorded benchmark artifact is evidence for exactly one source tree. Every
runner that persists an artifact records a `provenance` list into it
(`benchmarks/provenance.R`):

| Field | Meaning |
|---|---|
| `source_digest` | SHA-256 over the sorted per-file SHA-256 of `R/*.R` |
| `git_commit`, `git_dirty` | HEAD at record time, and whether `R/` was dirty |
| `package_version` | `packageVersion("crossform")`, or DESCRIPTION when the package is not yet installed |
| `r_version`, `platform`, `blas`, `lapack`, `lapack_version` | numerical environment |
| `runner`, `recorded_at` | which script wrote it, when |

`tests/testthat/helper-certification.R` refuses to read a recorded gate as a
boolean unless that provenance still binds it, in one of two tiers:

1. **`source_tree`** - a source checkout is reachable from the test directory,
   so the recorded `source_digest` must equal the current digest of `R/*.R`.
   Any change to package source invalidates recorded certification until the
   runner is re-run. This is the tier that applies to `testthat::test_local()`.
2. **`package_version`** - under `R CMD check` the package sources are not
   shipped (`R/` holds a lazy-load database, not `.R` files) and `.git` is
   absent, so only the recorded package version can be compared. Binding is
   version-level, and the tier is reported so it is never mistaken for the
   stronger one.

A non-binding artifact makes its test **skip with a loud message** naming the
runner that re-certifies it, rather than pass silently. Hard failure was
rejected deliberately: it would turn every ordinary source edit into a red
suite, which trains people to disable the check. The skip reasons are
greppable (`CERTIFICATION UNBOUND`, `CERTIFICATION STALE`,
`CERTIFICATION ABSENT`), and a pass at the weaker tier announces itself with a
`CERTIFICATION TIER package_version` message so that a version-level binding
is never read as a source-level one.

`benchmarks/run-measurement-benchmarks.R` and
`benchmarks/run-support-index-benchmark.R` persist nothing - they print to the
console - so they carry no provenance and certify nothing.

## Where artifacts live

Runners write to `benchmark-results/`, which is git-tracked for the two large
validation records but never shipped (`.Rbuildignore`). The small gate
artifacts the test suite reads live in `inst/extdata/certification/`, so that
certification **runs** under `R CMD check` instead of skipping. Promotion is
one script:

```sh
Rscript benchmarks/promote-artifacts.R benchmark-results .
```

It refuses to promote an artifact that carries no provenance, and refuses
anything above a 64 KiB shipped cap. The two large validation records
(`sampling-covariance-validation.rds`, 1.0 MB;
`learned-metric-policy-validation.rds`, 239 KB) are never promoted; their
tests read them from `benchmark-results/` when a local run has produced them
and skip with a message when it has not.

Tests resolve an artifact in this order: a fresh unpromoted run in
`benchmark-results/`, then `system.file("extdata", "certification", ...)`,
then the checkout's `inst/extdata/certification/`.

## Map-scale admission coverage

A recorded gate is evidence for a **map-scale compute verb**: an export that,
given a compiled plan and neural values, produces a scientific result at many
spatial measurements. Constructors, printers, inspectors, and adapters do not
count. A verb is certified only when its row's artifact is in the promote
table, is under the 64 KiB shipped cap, and still binds.

The machine-readable coverage list is `benchmarks/admission-coverage.R`.
`promote-artifacts.R` derives the shipped set from it.
`tests/testthat/test-admission-coverage.R` fails if a certified row is missing
from the promote table, or if a shipped `.rds` has no listed role.

| Verb | Artifact | Role |
|---|---|---|
| `rdm()`, `contrast_energy()` | `public-map-scale-gate.rds` | certified |
| query-first `rdm()`, `rsa()`, `contrast_energy()` | `query-first-scale-gate.rds` | certified |
| `evaluate_geometry()` | query-first gate | covered |
| `materialize_geometry()` | public-map and query-first late/comparator paths | covered |
| `crossnobis()` | `crossnobis-scale-gate.rds` | certified |
| `rdm_sampling_covariance()` | `sampling-covariance-scale.rds` | certified |
| `plan_relation()`, `estimate_relation()`, `fmrireg_relation()` | `first-moment-vertical-slice.rds` | certified |
| sequential memory contraction | `small-dense-memory-cold.rds`, `medium-sparse-memory-cold.rds`, `medium-sparse-block-cold.rds`, `medium-sparse-memory-warm.rds` | certified |
| `measurement_form()` | none yet (`run-measurement-profile.R` is print-only) | **gap** |
| shard executor | `shard-admission.rds` | refused (unbound) |

Local-only records stay in `benchmark-results/` and are never promoted:

| Record | Why it stays local |
|---|---|
| `sampling-covariance-validation.rds` (1.0 MB) | 10,000-rep Monte Carlo; exceeds the shipped cap |
| `learned-metric-policy-validation.rds` (239 KB) | 500-rep statistical recovery; exceeds the shipped cap |

Out of scope: `compile_frame()`, `plan_geometry()`, `relation()`, `pairing()`,
print and plot methods, `sampling_capabilities()`, neuroim2 adapters, and
every other constructor or inspector. They do not need a 576-node gate.

The sections below say what each retained gate asserts.

## Public geometry map gate

Run the matched implicit-identity and explicit fixed-diagonal map paths with:

```sh
Rscript benchmarks/run-public-map-scale-gate.R . benchmark-results
```

The isolated 576-feature, 576-center fixture uses 12 conditions, eight
partitions, and 66 RDM coordinates. Both public plans are warmed once and then
timed three times in alternating order. Sixteen spatial nodes are also checked
against a direct loop over condition contrasts, partition pairs, and
searchlight weights. The gate requires `1e-12` oracle parity, `1e-12` matched-
path parity, at most 60 seconds per complete map, an explicit-to-implicit
median runtime ratio no greater than five, and no more than 1 GiB incremental
peak RSS in the isolated measured phase.

The same isolated worker also fits an error-bearing 12-condition relation and
evaluates the factorized RDM-variance diagonal over a 576-center frame with a
mean support of about 42 features. Twenty-five evenly spaced nodes are warmed
and repeated three times before one complete sweep. The sampling arm must be
bit-identical across those reads, finite and nonnegative, no slower than one
second per probe node, and complete within ten minutes. The broad budget is a
cliff detector; the recorded runtime is the evidence.

The intentionally generous absolute thresholds are regression shields, not
marketing claims. The relative threshold specifically prevents a validated
domain-wide diagonal metric from falling back to the historical minutes-long
path while the matched implicit metric remains subsecond. Metric construction
is recorded separately from repeated map execution.

`tests/testthat/test-public-map-scale.R` invokes the same runner when
`CROSSFORM_RUN_SCALE_TESTS=true`. The scheduled and manually dispatchable
scale workflow sets that flag, which also activates the two 50k topology tests.

## Memory benchmark harness

Run from the repository root:

```sh
Rscript benchmarks/run-memory-benchmarks.R . /tmp/crossform-memory-results
```

Each fixed-seed scenario runs the actual sparse additive-frame path in a fresh
child R process: relation blocks, featurewise cross-Gram atoms, sparse
contraction, retained local relations, and coherent geometry. A ready-file
handshake starts resident-memory polling only after fixture construction (and,
for warm scenarios, warm-up) has completed. This prevents setup peaks from
being mislabeled as contraction memory while preserving cold and warm kernel
regimes explicitly.

The resulting RDS records measured child-process peak RSS, `Rprofmem`
allocation counts and sizes, every named live kernel temporary, total and
coherent output-storage bytes, and the conservative memory plan.

The plan predicts **crossform-owned workspace**, not process RSS: its
`prediction_kind` is `crossform_owned_workspace_upper_bound`, and its named
categories (frame, resident source, source handles and blocks, relation and
atom blocks, local state, output, contraction, replacement copy, serialization
overlap, reorder and checkpoint buffers) are summed and multiplied by a 1.25
safety factor. It deliberately makes no claim about interpreter, BLAS,
allocator-arena, or native pages, so it must not be compared against process
RSS. `tests/testthat/test-benchmark.R` therefore checks the plan against what
it actually bounds - the measured live temporaries and durable local relations
- and checks the recorded RSS separately as an OS measurement of an isolated
child polled only after the fixture signalled ready.

Version 0.1 has no worker processes, so worker peak RSS is explicitly zero;
future executor benchmarks must aggregate coordinator and worker peaks rather
than reporting coordinator memory alone.

Representative dimensions and frame density are versioned by
`.memory_benchmark_scenarios()`, and every scenario it declares must have a
recorded artifact.

`memory-benchmark-summary.csv` is the compact table: OS peak RSS, incremental
peak above the signalled pre-contraction baseline, the planned workspace,
allocation totals, and the largest measured live temporary.

## Shard admission benchmark

The optional executor gate is separate from the sequential memory harness. It
requires installed development versions of `crossform` and `shard` in an
isolated library, then runs:

```sh
Rscript benchmarks/run-shard-admission.R \
  . benchmark-results /path/to/isolated/library
```

The harness compares sequential response extraction, shared-response staging,
and shared-relation staging. Cold modes include staging and pool startup; warm
modes measure execution with those resources already available. Every mode
uses `crossform`'s canonical feature task and ordered reducer. The parent
process samples the full child process tree, so the memory record includes the
coordinator and all workers rather than coordinator RSS alone.

Admission requires total/coherent numerical parity, normal and injected-error
cleanup, and at least 1.1x speedup on a cold end-to-end mode. The 2026-08-12 run
passed parity and cleanup but failed runtime, so the adapter was not admitted.
`shard-admission.rds` is the complete evidence record and
`shard-admission-summary.csv` is its compact table.

## Learned-metric policy validation

Run the paired statistical validation in an installed package environment:

```sh
Rscript benchmarks/run-learned-metric-policy-validation.R 500 benchmark-results
```

The harness executes the public `lm_relation_fit() -> plan_crossnobis() ->
crossnobis()` path under two residual-training policies. It reports evaluator
error relative to each realized learned metric separately from error relative
to the population Mahalanobis target. A Gaussian iid arm tests the declared
GLM orthogonality condition; a stationary AR(1) arm deliberately fits the
wrong observation metric. Policy comparison is an equivalence analysis with a
predeclared 0.005 margin, not a failed difference test. The estimated
shrinkage target is recorded, and no confidence interval for scientific data
is implied because uncertainty in the learned metric is not calibrated.

## Learned-crossnobis brain-scale gate

Run the named 52,416-feature volume fixture in an isolated child process:

```sh
Rscript benchmarks/run-crossnobis-scale-gate.R . benchmark-results
```

The parent begins OS RSS polling only after the child has built its raw-data
fixture, fitted relation envelope, and sparse searchlight frame. The gate then
includes residual pair accumulation, an on-demand dense local shrinkage metric,
one independent crossnobis contrast, and exact NeuroVol back-mapping. It records
support and union-pair topology, local factorization work, residual and effect
source reads, planned workspace, incremental peak RSS, stage timings, and
whether any pair frame, pair atoms, or node-by-edge factor table was retained.
The initial admission limits are 4 GiB incremental RSS and 30 minutes for the
complete analysis.

Since the union pair graph became lazy, an ordinary searchlight compile leaves
`pair_pattern` NULL and the crossnobis path never needs it, so the recorded
`structural_bytes` is the support index the analysis actually used and no
longer includes the union graph. The `union_pair_stored_nnz` and `union_degree`
topology figures come from a receipt-only materialization taken after the
measured analysis, and `pair_pattern_route` records which route the analysis
ran.

This fixture qualifies execution and storage only. Its training-only arm has
30 residual degrees of freedom and support sizes up to 33, so the declared
shrinkage estimator is load-bearing. Statistical recovery and residual-reuse
policy claims belong to the separate 500-replication validation above.

## Query-first scale gate (Gate 5)

Certify that selected RDM edges, the full RDM, and fixed linear RSA
coefficients execute at q = 100 conditions over 1,080 searchlights without
materializing the packed geometry field or any dense pair-query matrix:

```sh
Rscript benchmarks/run-query-first-scale.R . benchmark-results
```

The worker validates three things before timing begins: an independent
raw-beta cross-partition oracle over probe nodes and selected pairs, exact
agreement between `rdm(plan, pairs = )` columns and the corresponding full-RDM
columns, and route-stable scientific identity between the fused and
materialize-then-project executions. The gate then requires the fused full
RDM to be at least as fast as materialize-then-project and selected pairs to be
no slower than the full fused sweep. Three additional fresh workers each run
exactly one of selected-RDM, full-fused-RDM, and RSA execution. Each worker
resets R's high-water heap counter after constructing the common fixture; the
largest one-call heap increment must remain below 512 MiB. This is explicitly
an R-heap claim, not an OS RSS claim. The materialized comparator remains in
the timing court and is excluded from the query-first memory claim. The
receipt separately records the allocations avoided by structured
execution: a ~200 MB (191 MiB) dense packed query matrix and an ~87 MB
(83 MiB) two-component geometry field. The recorded artifact reports 100
selected pairs in 0.13 s, the full fused RDM in 1.59 s, and the materialized
route in 2.08 s (fused/materialized ratio 0.76 against a 1.2 ceiling — the
native packed-form kernel made the materialized comparator fast, so this
margin is now thin and worth watching), a 277 MB (264 MiB) maximum
fresh-worker incremental R heap, and oracle error 4.4e-16.

## Call-graph components

`benchmarks/call-graph-scc.R` is analysis rather than a gate: it records no
artifact, binds no provenance, and nothing skips on it. It rebuilds the
file-level call graph of `R/` using the same extraction rule as
`tests/testthat/test-architecture.R` — top-level `name <- function(...)`
definitions, then every symbol a file references that resolves to a definition
in another file — and reports the strongly connected components with Tarjan.

```sh
Rscript benchmarks/call-graph-scc.R
SCC_FOCUS=receipt.R Rscript benchmarks/call-graph-scc.R
```

It prints the largest component's size and membership, every component of more
than one file, and the edges into and out of `SCC_FOCUS` (default `receipt.R`).
This is what produces the numbers `design/architecture.md` quotes under "What
remains", so re-run it before changing them.
