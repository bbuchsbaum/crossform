# Predictive geometry performance evidence

The benchmark compares declared numerical work, with no universal speedup
threshold. Run it from the repository root:

```sh
Rscript benchmarks/predictive-geometry/run-performance.R . benchmark-results/predictive-geometry-performance all
```

Each primary case runs in a fresh child process. Workers execute serially.
They load the package, construct metadata and admit their workspace before
the parent starts monitoring execution. A readiness handshake ensures RSS
polling is active before neural reads. The measured region includes fitting,
scoring, requested persistence, and bounded extraction of the comparison
outputs. Reference checks run afterwards, outside that region. Hardware,
available memory and host load are recorded; an idle machine is not assumed.

## Matched work

| Case | Conditions | Model rank | Partitions | Neural features | Measurements | Local support | Response rank |
|---|---:|---:|---:|---:|---:|---:|---:|
| Tiny | 12 | 4 | 4 | 64 | 32 | 8 | 2 |
| Medium | 64 | 12 | 8 | 2,048 | 2,000 | 32 | 6 |
| Large | 200 | 20 | 8 | 20,000 | 10,000 | 32 | 8 |

All cases use deterministic bounded matrix readers, a trace-normalized
anisotropic model, penalty 0.005, independent training/evaluation partition
sets, nonuniform within-set edge weights, the identity neural metric, total
geometry and local frame normalization. Feature blocks are at most 128;
prediction/readout blocks are at most 64 measurements. The configuration and
generator are retained with the artifact.

Tiny and medium compare public memory and block execution with an original-
effect-space complete-form baseline. The latter uses the same executor,
partitions, frame, metric and feature blocks, then compresses the full forms
and applies an independent base-R fit and score. A pilot caught unused
parent partitions being read by this baseline; it now subsets the original
relation to the same actual endpoints before execution. Read counts must
agree across matched routes.

The baseline returns the same numerical prediction, signed/shifted spectra,
diagnostics, mode amplitudes/evidence/gain and total gain. It omits the
predictive record's validation and durable storage machinery. Timings are
therefore labelled by route; they are not a claim that all surrounding
software responsibilities cost the same. Reported retained sizes also name
their units: R object bytes, durable file bytes, or plain numerical-reference
bytes.

For tiny and medium, every measurement is additionally checked against raw
weighted partition products constructed independently of the executor and
packed codec. The large block case checks 64 prespecified locations against
that dense oracle. It does not allocate a full 200-condition geometry at
10,000 measurements merely to make a timing comparison.

The existing descriptive reader is compared only on overlapping outputs:
signed model-space trace and the unpenalized shared PSD form. Its projector
energy is not equated with predictive gain, and no runtime comparison across
those different statistical tasks is presented.

## Memory and reuse

Every route is admitted against a 512 MiB workspace budget before neural
reads. The dense baseline additionally reserves overlap with its retained
first full form. Public-route receipts include predictive buffers and the
existing compiler's workspace model. The large blocked request initially
admits 145,792,120 planned bytes.

The parent polls OS RSS; workers record baseline and final RSS. `Rprofmem`
records cumulative R allocations and the largest single allocation without
loading its complete log into memory. Native work is additionally observed
through compiler receipts, kernel counters and RSS. Planned buffer bytes,
cumulative allocation and incremental RSS are different quantities: R heap
retention, JIT/runtime allocation and allocator pages contribute to RSS. The
512 MiB admission is not a retroactively chosen RSS threshold.

Instrumentation records actual native-kernel calls, eigen calls, fitted
spectral paths, materializations, source calls/columns, maximum read width
and requested packed width. Predictive materialization is guarded against
the original effect dimension. Source geometry uses `r(r+1)/2` coordinates,
not `n(n+1)/2`.

A separate warm comparison on tiny and medium runs six fixed recipes:
three response ranks at each of two penalties. Cached and uncached runs
must agree on every candidate gain and sampled prediction form. Counters
require two cached materializations versus twelve uncached ones, and two
spectral paths per measurement versus six. These timings run in the order
uncached then cached and are labelled accordingly; they do not establish a
general performance ratio.

## Artifacts and gates

`benchmark-results/predictive-geometry-performance/` contains individual
source-bound case records, worker logs, a combined receipt and a summary CSV.
Every record binds the R source, all four harness files, configuration,
platform and execution observations. The runner refuses source/harness drift
and fails on numerical or matched-read disagreement. The relative-error gate
is fixed at `1e-8`; no wall-time threshold is inferred from this first run.

`test-predictive-performance.R` checks metadata-only admission for all three
sizes, the tiny public/full/raw comparison, existing-reader overlap and
fixed-recipe reuse. Later edits to production R invalidate the measurement's
source binding and require a rerun for final certification.

## Initial measured baseline

The 2026-09-05 run binds source
`sha256:10bfc4506bcf9c110f8ac10a0bfba425fb0f3eb4c5975e608bd1cd095690bf12`.
All seven routes passed. Maximum independent dense relative error was
`3.24e-14`. Tiny and medium routes also agreed on all requested numerical
outputs and exact source-read counts.

| Case / route | Seconds | Planned MiB | Incremental RSS MiB | Largest R allocation MiB | Packed coordinates |
|---|---:|---:|---:|---:|---:|
| Tiny / memory | 0.856 | 0.54 | 20.94 | 0.06 | 10 |
| Tiny / block | 0.910 | 0.52 | 15.83 | 0.06 | 10 |
| Tiny / dense numerical reference | 0.589 | 0.37 | 9.89 | 0.06 | 78 |
| Medium / memory | 5.503 | 22.67 | 222.81 | 2.20 | 78 |
| Medium / block | 5.984 | 16.83 | 231.98 | 2.20 | 78 |
| Medium / dense numerical reference | 4.716 | 328.66 | 1191.86 | 31.74 | 2080 |
| Large / block | 65.949 | 139.04 | 1273.69 | 16.02 | 210 |

The dense numerical baseline is faster on these tiny/medium tasks; the
predictive route also constructs and validates its richer records. Its
advantage demonstrated here is reduced geometry width and substantially
smaller medium-case peak allocation/RSS, not a universal latency improvement.
The large run retained approximately 22.24 MB of durable prediction/score
files and read 160,000 feature columns across the eight sources, with no
read wider than 128 features.

The observed large-process RSS increase is about 1.24 GiB despite a 139 MiB
buffer plan. These measurements do **not** support treating `workspace_bytes`
as a hard process-memory cap or assuming that 512 MiB of free RAM suffices.
They expose R/runtime overhead alongside the declared numerical workspace;
both quantities must remain in future comparisons.

The six-candidate medium reuse comparison reduced materializations from
12 to 2 and fitted spectral paths from 12,000 to 4,000. Its instrumented warm
time changed from 10.718 to 9.479 seconds, with identical candidate gains and
prediction witnesses. The cached run reserved 149,139,320 bytes. The tiny
comparison likewise changed 12 to 2 materializations and 192 to 64 paths.
