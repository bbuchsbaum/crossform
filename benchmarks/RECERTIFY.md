# How to re-certify after editing `R/`

Every recorded benchmark artifact is evidence for exactly one source tree. Each
runner stamps a `provenance` block into the artifact it writes, including a
`source_digest`: the SHA-256 over the sorted per-file SHA-256 of `R/*.R`
(`benchmarks/provenance.R`). `tests/testthat/helper-certification.R` refuses to
read a recorded gate as a boolean unless that digest still equals the current
one, so **any** edit to any file under `R/` — even a comment — turns all eight
certification tests into loud `CERTIFICATION STALE` skips until the runners are
re-run and the artifacts re-promoted.

This file is the exact sequence that clears those skips. See
`benchmarks/README.md` § "Map-scale admission coverage" for which public
verbs the gates certify, and the later sections there for what each gate
asserts. The machine-readable coverage list is
`benchmarks/admission-coverage.R`. Recertify is re-running that table on a
frozen `R/` digest, not inventing new gates.

## Before you start: freeze `R/`

The whole sequence takes about 12 minutes and the digest must not move for any
of it. If `R/` changes partway through, the artifacts written before the change
bind to a digest that no longer exists, promotion ships evidence for a source
state that is gone, and the tests skip `STALE` again. Nothing warns you except
the skip.

Check that the tree is quiet and remember the digest:

```sh
cd /path/to/crossform
git status --porcelain -- R          # expect no output
Rscript -e 'source("benchmarks/provenance.R"); cat(.crossform_source_tree_digest("."), "\n")'
```

Run the same command again at the end. If the two digests differ, the run is
void — do it again on a quiet tree.

Do not run anything else heavy on the machine while the gates execute. The
public map gate and the crossnobis gate measure wall-clock time and resident
memory; under heavy CPU contention the public map sampling sweep has been
observed at 889 s against its 600 s budget, and passed at 54 s on the same
source once the machine was idle. A timing gate failure on a loaded machine is
not evidence of a regression — re-run it alone before believing it.

## 1. Install the package into a temp library

Only `run-learned-metric-policy-validation.R` needs an installed copy; it calls
`library(crossform)`. Every other runner loads the source tree directly with
`pkgload::load_all()` or `devtools::load_all()`.

Install from an out-of-tree copy so that `R CMD INSTALL` does not disturb the
`src/*.o` and `src/*.so` objects that `load_all()` is using in your checkout:

```sh
LIB=$(mktemp -d)/lib && mkdir -p "$LIB"
SRC=$(mktemp -d)
rsync -a --exclude '.git' --exclude 'exemplars' --exclude 'benchmark-results' \
  --exclude 'dist' --exclude 'doc' --exclude '*.tar.gz' \
  --exclude 'src/*.o' --exclude 'src/*.so' ./ "$SRC/crossform/"
(cd "$SRC" && R CMD INSTALL --library="$LIB" --no-docs crossform)
```

**Exclude `exemplars/`.** It holds the downloaded Haxby and population
datasets and is now about 11 GB; without the exclusion this copy takes many
minutes and fills the temp volume before `R CMD INSTALL` ever starts. The
excluded directories are build products and data, none of which the install
needs — with them out, the copy is roughly 9 MB and finishes instantly.

## 2. Run every persisting runner, in this order

Run them serially, from the package root. Approximate runtimes are from an
idle Apple silicon laptop; treat them as order-of-magnitude.

```sh
Rscript benchmarks/run-memory-benchmarks.R . benchmark-results              #  ~15 s
Rscript benchmarks/run-sampling-covariance-scale.R                          #   ~6 s
Rscript benchmarks/run-sampling-covariance-validation.R                     #  ~10 s
Rscript benchmarks/run-population-null-coverage.R 2000 benchmark-results     #  ~32 s
Rscript benchmarks/run-first-moment-vertical-slice.R . benchmark-results    #  ~18 s
Rscript benchmarks/run-public-map-scale-gate.R . benchmark-results          #  ~68 s
Rscript benchmarks/run-query-first-scale.R . benchmark-results              #  ~40 s
Rscript benchmarks/run-crossnobis-scale-gate.R . benchmark-results          #  ~64 s
R_LIBS="$LIB" Rscript benchmarks/run-learned-metric-policy-validation.R 500 benchmark-results   # ~260 s
```

Notes that cost real time to rediscover:

- `run-sampling-covariance-scale.R` and `run-sampling-covariance-validation.R`
  take no repository argument. They use `getwd()`, so run them from the root.
- **Give `run-sampling-covariance-validation.R` no repetition argument.** Its
  default is 10,000 Monte Carlo repetitions and the test bands are calibrated
  to that. At 1,000 repetitions the run is honest but too noisy, and it fails
  its own recorded bands (a null variance ratio of 1.09 against a 0.95-1.05
  band, and a null covariance relative error of 0.066 against a 0.05 limit).
  The argument exists for debugging, not for certification.
- The `500` for the learned-metric validation is pinned:
  `test-certification-artifacts.R` asserts
  `contract$replications == 500L`.
- **Give `run-population-null-coverage.R` the `2000`.** The runner refuses
  fewer than 500 replications, and the coverage bands
  `test-population-uncertainty.R` ratchets were measured at 2,000 (Monte Carlo
  standard error 0.005). At 500 the misspecification arm's floor is inside the
  noise. Its second argument is the results directory; without one it writes to
  `benchmark-results` anyway.
- Each runner exits nonzero when its gate fails. **Do not lower a threshold to
  make a gate pass** — that is the one move that makes the whole apparatus
  worthless. Record the failure and fix the code.
- **A runner that aborts is usually harness drift, not a gate failure.** The
  runners reach into result and plan internals, and a refactor under `R/` can
  move a field without any test noticing, because the only reader is the
  runner. The 2026-08-20 re-certification hit five such reads at once (ticket
  B3 turned `plan_crossnobis()` into an `effect_geometry_plan`). The
  signatures to recognise:
  - `Error in !result$work$… : invalid argument type` — the field read as
    `NULL`; `!NULL` is an error in R 4.5, not `logical(0)`.
  - `arguments imply differing number of rows: 1, 0` from the summary
    `data.frame()` — some recorded field is `NULL`.
  - `values must be length 1, but FUN(X[[1]]) result is length 0` from a
    `vapply` over plans — a plan-level field is gone.
  - A package-side refusal such as "Frozen metric-schedule fields are missing
    or noncanonical" — the runner passed a wrapper where a sealed object is
    wanted (`metric_schedule` vs its inner `metric_schedule$schedule`).

  Repair the read, never the assertion, and prove the repair by checking that
  a structural number the refactor should not have touched still reproduces
  the previous record exactly — residual reads and planned workspace bytes are
  the two that pinned the 2026-08-20 run.
- **Silent `NULL`s are the quieter version of the same failure.** A field that
  moves but is never dereferenced records `NULL` into the shipped artifact and
  nothing fails. After promoting, walk each artifact for zero-length leaves.
  Known-legitimate ones as of 2026-08-20: the four memory scenarios'
  `plan$budget_bytes` and `plan$measured_*`/`*_rss_bytes` (no budget is set for
  those scenarios), the query-first `identity$fused_lowering` and
  `fused_materialization` (the additive-query-fused view does not populate
  `execution_plan`), and everything under the unbound `shard-admission.rds`.
  Anything else is drift.
- `run-shard-admission.R` is deliberately excluded. It needs an installed
  `shard` at a pinned version in an isolated library, and its recorded
  artifact carries no provenance at all, so its test always skips
  `CERTIFICATION UNBOUND`. That is the designed state, not a defect.
- `run-measurement-benchmarks.R`, `run-measurement-profile.R`, and
  `run-support-index-benchmark.R` persist nothing and certify nothing; they
  only print. `measurement_form()` is a documented coverage gap until it has
  a compact promoted gate.

## 3. Promote the small gate artifacts

Runners write to `benchmark-results/`, which is never shipped. The artifacts
the test suite reads under `R CMD check` live in
`inst/extdata/certification/`:

```sh
Rscript benchmarks/promote-artifacts.R benchmark-results .
```

Promotion refuses any artifact with no source digest and anything over the
64 KiB shipped cap, and exits nonzero if it refused something. The two large
validation records are never promoted; their tests read them from
`benchmark-results/`.

One trap: promotion copies a `*-summary.csv` whenever the file is present,
without checking provenance, because the summaries are human-readable receipts
rather than test evidence. If a gate failed and left a stale summary from an
earlier run in `benchmark-results/`, that summary is promoted next to an
unrelated `.rds` and the shipped receipt then disagrees with the artifact
beside it. **Delete the summary of any gate that did not produce an `.rds` in
this run** before promoting.

## 4. Verify

```sh
Rscript -e 'testthat::test_local(".", filter = "certification|scale|vertical|benchmark", reporter = "summary")'
Rscript -e 'source("benchmarks/provenance.R"); cat(.crossform_source_tree_digest("."), "\n")'
```

Expect zero failures. Exactly these skips are correct and expected:

| Skip | Why |
|---|---|
| `the recorded executor admission record refuses the adapter` | `CERTIFICATION UNBOUND` — shard admission is excluded by design |
| `the public map-scale benchmark passes its declared gate` | opt-in, needs `CROSSFORM_RUN_SCALE_TESTS=true` |
| `the query-first scale benchmark passes its declared gate` | opt-in, needs `CROSSFORM_RUN_SCALE_TESTS=true` |

Any remaining `CERTIFICATION STALE` skip names the runner that re-certifies it
in the skip message. If one appears, that runner either failed or ran against a
different source state than the one you are testing.

To also execute the two scale gates live rather than reading their receipts:

```sh
CROSSFORM_RUN_SCALE_TESTS=true Rscript -e 'testthat::test_local(".", filter = "scale", reporter = "summary")'
```

That flag additionally activates the two 50k topology tests in
`test-support-index.R`. The scheduled scale workflow sets it.
