# Memory benchmark harness

## Public geometry map gate

Run the matched implicit-identity and explicit fixed-diagonal map paths with:

```sh
Rscript benchmarks/run-public-map-scale-gate.R . benchmark-results
```

The isolated 576-feature, 576-centre fixture uses 12 conditions, eight
partitions, and 66 RDM coordinates. Both public plans are warmed once and then
timed three times in alternating order. Sixteen spatial nodes are also checked
against a direct loop over condition contrasts, partition pairs, and
searchlight weights. The gate requires `1e-12` oracle parity, `1e-12` matched-
path parity, at most 60 seconds per complete map, an explicit-to-implicit
median runtime ratio no greater than five, and no more than 1 GiB incremental
peak RSS in the isolated measured phase.

The same isolated worker also fits an error-bearing 12-condition relation and
evaluates the factorized RDM-variance diagonal over a 576-centre frame with a
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
coherent output-storage bytes, and the conservative memory plan. The plan
includes a named 64 MiB process-level runtime reserve for measured R,
Matrix/BLAS initialization, allocator arenas, and native pages that cannot be
assigned to one scientific array. This is an explicit benchmark-calibrated
category, not an unnamed multiplier. Version 0.1 has no worker processes, so
worker peak RSS is explicitly zero; future executor benchmarks must aggregate
coordinator and worker peaks rather than reporting coordinator memory alone.

The harness is evidence, not a timing assertion in the fast test suite.
Representative dimensions and frame density are versioned by
`.memory_benchmark_scenarios()`.

`summary.csv` compares the conservative plan with measured incremental peak
RSS above the signaled pre-contraction baseline. Every scenario must report
`plan_covers_incremental_peak = TRUE`; a false value blocks both closure of the
sequential compiler epic and every later parallel-executor gate.

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
selected pairs in 0.27 s, the full fused RDM in 5.16 s, and the materialized
route in 13.40 s (fused/materialized ratio 0.39), a 291 MB (277 MiB) maximum
fresh-worker incremental R heap, and oracle error 4.4e-16.
