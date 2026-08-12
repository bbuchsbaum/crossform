# Memory benchmark harness

Run from the repository root:

```sh
Rscript benchmarks/run-memory-benchmarks.R . /tmp/effectagram-memory-results
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
requires installed development versions of `effectagram` and `shard` in an
isolated library, then runs:

```sh
Rscript benchmarks/run-shard-admission.R \
  . benchmark-results /path/to/isolated/library
```

The harness compares sequential response extraction, shared-response staging,
and shared-relation staging. Cold modes include staging and pool startup; warm
modes measure execution with those resources already available. Every mode
uses `effectagram`'s canonical feature task and ordered reducer. The parent
process samples the full child process tree, so the memory record includes the
coordinator and all workers rather than coordinator RSS alone.

Admission requires total/coherent numerical parity, normal and injected-error
cleanup, and at least 1.1x speedup on a cold end-to-end mode. The 2026-08-12 run
passed parity and cleanup but failed runtime, so the adapter was not admitted.
`shard-admission.rds` is the complete evidence record and
`shard-admission-summary.csv` is its compact table.
