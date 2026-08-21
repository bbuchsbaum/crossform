# Statistical and brain-scale validation receipt

Date: 2026-08-15

Umbrella: `bd-01KZXA4C1RJDV567TS9X175ND0`  
Final evidence epic: `bd-01KZXA62QZQMH79RS7KW3PVCGT`

## Errata / reconciliation (2026-08-15)

This receipt is preserved as issued. Two of its statements do not survive a
reading against the repository, and are corrected here rather than silently
edited above.

**Export count.** "Public surface and scale regimes" states that the artifact
exports 78 functions. `NAMESPACE` exports 105 functions and registers 207 S3
methods — 101 `print`, 95 `format`, 6 `as.data.frame`, 5 `plot` (source:
`NAMESPACE`). The semantic grouping in that section is still accurate; the
count is not.

*Erratum to the erratum (2026-08-17).* This paragraph as first issued said "43
S3 methods". That figure was itself stale; the registered total is 207. Both
counts are now recorded, with their derivation, in "Recorded certification
metrics" below.

**Public-data analysis.** "Explicitly unrun or deferred evidence" states that
no public-data matched-estimand analysis was run. The Haxby 2001 exemplar
(`exemplars/haxby2001`, smoke report dated 2026-08-13) had already been run.
What it does and does not show:

- *Shows implementation parity on a matched estimand.* On condition-level
  crossvalidated squared Euclidean distance under the identity metric, averaged
  over cross-run partition pairs, `crossform`'s `rdm()` agrees with rMVPA's own
  crossnobis estimator to a maximum absolute difference of `8.88e-16` at all
  577 VT searchlight centers, and with an independent reference loop written
  from the definition to `1.33e-15`. Refitting the raw responses to recover the
  error channel reproduces the point RDM to `4.44e-16`.
- *Does not show matched-estimand scientific validation.* This is numerical
  agreement among three implementations of one estimand, on one subject, at the
  smoke tier. There is no ground truth, no external criterion, no group or
  population inference, and no calibration of the analytic standard errors the
  exemplar's uncertainty arm reports. The estimand the exemplar originally
  specified — a `1 - Pearson` correlation-distance RDM scored by a Spearman
  correlation — is expressible by neither package, so the compared estimand is
  the nearest one both systems express exactly, not the one first intended. No
  matched-estimator speed advantage is claimed.

The remainder of that section stands: no hosted CI, cross-platform, CRAN,
R-universe, or downstream-consumer court was run, and no interval, LD-t,
bootstrap, permutation, or group-level calibration is implemented.

## Verdict

The evidence-pairing architecture now has one complete local statistical
vertical slice:

```text
raw responses
  -> relation_fit with residual error channel
  -> residual pair sufficient statistics
  -> provenance-frozen, on-demand local precision recipe
  -> compiled evidence plan
  -> signed query-first crossnobis
  -> exact NeuroVol index reconstruction
```

This is locally algebraically established, numerically verified,
statistically validated in the stated simulation regimes, and scale-qualified
on the named 52,416-feature fixture. Those four statuses are deliberately
separate. This report does not claim hosted CI, cross-platform behavior, CRAN
acceptance, calibrated intervals, LD-t, or arbitrary real-data validity.

## Exact artifact and reproducibility boundary

| Property | Value |
|---|---|
| Artifact | `crossform_0.0.0.9000.tar.gz` |
| SHA-256 | `b18e47ca1933caf671ed65a79f64d56f53b26edac5f1535eff4af6f9f236dbfa` |
| Payload size | 433,291 bytes |
| R | 4.5.1 (2025-06-13) |
| Platform | `aarch64-apple-darwin20` |
| Host | macOS Sonoma 14.3 |
| Locale | `en_US.UTF-8` |
| Exact local check | `R CMD check --no-manual`, CRAN incoming disabled |
| Result | `Status: OK` |

The artifact was built from the staged workspace immediately before
publication. Its hash identifies the exact validated payload, and the
publication commit records the same package-payload paths. This establishes a
commit-reproducible local receipt; it is not a release, hosted-CI, or
cross-platform certification.

The same artifact also passed `R CMD check --as-cran --no-manual` with no
errors or warnings and one expected incoming note: the development version
`0.0.0.9000` contains large components. Remote incoming checks were disabled;
the deterministic local package court above is the final artifact's
authoritative result, and incoming/hosted evidence remains qualified.

Workspace-only `.claude`, `.tmp`, `.mote`, benchmark artifacts, and this report
are excluded from the package payload. The exact tarball contains none of
those paths.

## Package rename integrity

The 2026-08-15 rename to `crossform` is present in package metadata, namespace
loading, generated help, vignettes, examples, benchmark drivers, scale-test
controls, and the pkgdown site. All three vignettes render against an installed
`crossform` package. The pkgdown site completes its reference, article, search,
redirect, and integrity stages when the articles are rendered in-process; the
installed pkgdown version's isolated article subprocess fails in its own error
wrapper before reporting the underlying condition.

Previously hashed scientific tokens are not rewritten. Both the legacy and
public design-protocol spellings lower to the same canonical semantic token,
while compiler-route and execution receipts carry `crossform`. The regression
suite verifies that this preserves `design_model_id` and `relation_plan_id`
while changing route and design-receipt identities. Remaining `effectagram`
strings in active source are limited to those frozen tokens, their compatibility
test, and an explicit note about historical exemplar receipts.

## Four evidence levels

### Algebraically established

The governing observable is

\[
\mathscr E_{LR}(H,K)
=
\operatorname{tr}(H^\top B_L K B_R^\top)
=
\langle H, B_L K B_R^\top\rangle_F
=
\langle B_L^\top H B_R, K\rangle_F.
\]

Direct matrix-law tests cover the forward/adjoint identity, rectangular and
self forms, second-order vectorization, reversal, explicit bridges,
decomposition lifting, normalization order, and dual-frame reconstruction.
The independent helpers in `helper-effect-form-laws.R`,
`helper-evidence-pairing-laws.R`, `helper-measurement-form-oracle.R`, and
`helper-residual-statistics.R` implement their reference calculations in base
R and call no package kernels.

The neural noise model enters the existing neural operator slot:

\[
K=I \quad\text{for Euclidean geometry},
\qquad
K=\widehat\Sigma^{-1} \quad\text{for noise-normalized geometry}.
\]

The relation remains the pure experimental-neural map. `relation_fit()` adds a
separate capability-bearing error channel with residual blocks, residual
degrees of freedom, effect-coordinate covariance, estimator identity, and
source provenance.

### Numerically verified

The final fast source court ran 49 test files and 362 tests:

```text
PASS 1822
FAIL 0
WARN 0
SKIP 2
```

> Superseded by the 2026-08-16 court recorded under "Certification provenance
> binding" below: `FAIL 0 | WARN 0 | SKIP 10 | PASS 4112`. The block above is
> the 2026-08-15 source-checkout run and is preserved as issued; 1,822 is not
> the current expectation count.

The two skips are the declared 50k topology tier. Running the same support-index
court with `CROSSFORM_RUN_SCALE_TESTS=true` produced:

```text
PASS 47
FAIL 0
WARN 0
SKIP 0
```

Product-to-oracle tests establish:

- known-metric crossnobis equals a direct evidence-pairing implementation;
- learned-metric crossnobis equals independently accumulated residual
  covariance, shrinkage, local solve, and edge reduction;
- raw-response and precomputed-effect paths agree;
- query-first and complete materialization agree on small fixtures;
- additive, support-streamed, and measurement-form lowerings preserve their
  declared values within the numerical contract;
- signed and indefinite crossvalidated estimates are retained;
- source revisions, workspace refusal, plan identity, and training provenance
  fail closed when forged or incompatible.

The 64-replication known-metric recovery guard pins per-replication reference
standard deviations and derives Monte Carlo standard-error thresholds from
the active replication count. Its bias and variance gates therefore do not
silently change when the test duration changes.

### Statistically validated

`benchmark-results/learned-metric-policy-validation.rds` records the exported
`lm_relation_fit() -> plan_crossnobis() -> crossnobis()` path over 500 paired
replications in each of two regimes. The policy contrast is all-run residual
reuse minus training-only estimation, with a predeclared equivalence margin of
plus or minus 0.005.

| Regime | Quantity | Mean difference | 95% Monte Carlo interval | Equivalent? |
|---|---|---:|---:|---|
| Correct Gaussian iid GLM, null | estimate | 0.0002205 | [-0.0000014, 0.0004424] | yes |
| Correct Gaussian iid GLM, signal | estimate | -0.0042779 | [-0.0049124, -0.0036434] | yes |
| Unmodelled AR(1), null | estimate | 0.0000164 | [-0.0000710, 0.0001038] | yes |
| Unmodelled AR(1), signal | estimate | -0.0093603 | [-0.0114013, -0.0073193] | no |

For the misspecified AR(1) signal arm, the conditional learned target changed
by -0.0089104, while evaluator error changed by -0.0004499 and remained within
the equivalence margin. The material difference is therefore attributable to
the metric learned under the two provenance policies, not to a second
evaluation engine. Training-only remains the conservative default; all-run
reuse remains explicit and justification-bearing.

This simulation validates those named regimes and estimators only. The
separable GLM error model is a declared capability, not a universal noise law.
No interval or LD-t is exported because calibration has not yet propagated
uncertainty in the estimated metric.

### Scale-qualified

`inst/extdata/certification/crossnobis-scale-gate.rds` (and its shipped
`crossnobis-scale-gate-summary.csv`) records one exact brain-scale execution.
The values below are the 2026-08-20 re-record described under "Re-certification
after WS-B" below, `provenance$recorded_at = 2026-08-20 19:51:26 UTC`, fixture
`crossnobis-scale-52k:v1`:

| Quantity | Measured value | Recorded field |
|---|---:|---|
| Active volume features | 52,416 | `fixture$features` |
| Searchlight radius / spacing | 6.1 mm / 3 mm | `fixture$radius_mm`, `fixture$spacing_mm` |
| Support size | min 11, median 33, mean 31.16, max 33 | `support$min/median/mean/max` |
| Union-pair stored nonzeros | 4,365,964 | `support$union_pair_stored_nnz` |
| Local metric derivations | 52,416 | `work$local_metric_derivations` |
| Factorization work units | 1,653,408,544 | `work$factorization_units` |
| Residual source reads | 11,163 | `reads$residual_total` |
| Evaluation source reads | 104,832 | `reads$evaluation_total` |
| Planned owned workspace | 281,589,160 bytes | `memory$planned_workspace_bytes` |
| Incremental peak RSS | 570,572,800 bytes | `memory$incremental_peak_rss_bytes` |
| OS peak RSS | 1,199,013,888 bytes | `memory$os_peak_rss_bytes` |
| Plan plus residual statistics | 41.365 s | `timing$plan_and_residual_statistics_seconds` |
| Crossnobis evaluation | 57.859 s | `timing$crossnobis_seconds` |
| Exact NeuroVol mapping | 0.035 s | `timing$output_mapping_seconds` |
| Complete measured analysis | 101.461 s | `timing$analysis_seconds` |

The structural rows are unchanged, which is the load-bearing invariance: the
WS-B refactors did not change how much work the brain-scale execution does or
how much workspace it plans. All fourteen structural quantities — features,
the four support statistics, union-pair nonzeros, structural bytes, local
metric derivations, dense metric entries, factorization units, both read
counts, planned workspace, and output values — compare identical between the
2026-08-17 and 2026-08-20 records; the 2026-08-15 correspondence is as
recorded by the previous erratum.

The timing and RSS rows supersede two earlier records, preserved here rather
than silently replaced:

| Row | 2026-08-20 (current) | 2026-08-17 | 2026-08-15 as issued |
|---|---:|---:|---:|
| Incremental peak RSS | 570,572,800 | 903,479,296 | 948,568,064 |
| Plan plus residual statistics | 41.365 s | 18.021 s | 64.348 s |
| Crossnobis evaluation | 57.859 s | 41.069 s | 57.790 s |
| Exact NeuroVol mapping | 0.035 s | 0.049 s | 0.084 s |
| Complete measured analysis | 101.461 s | 59.982 s | 123.498 s |

Same machine (Apple silicon, aarch64-apple-darwin20, R 4.5.1, Accelerate BLAS)
across all three. **The 2026-08-20 timings were measured under heavy load** —
the host carried a 15-minute load average near 87 from unrelated concurrent R
jobs, which `benchmarks/RECERTIFY.md` warns invalidates timing comparisons.
They are recorded because they are what this run measured, and the gate passed
its 30-minute limit with three orders of magnitude of headroom, but they are
not a clean-machine regression signal and the 101.461 s should not be read as a
69% slowdown against 2026-08-17. Peak RSS fell by 37% on the same comparison,
which contention does not explain and which is consistent with the WS-B
executor work. A clean-machine re-record is the only way to settle the timing
question.

The gate passed its 4 GiB incremental-RSS and 30-minute limits. It retained no
pair-atom field, `p^2` pair frame, or node-by-edge factor table. Compact output
was mapped back to the exact full-volume indices with the original NeuroSpace,
and negative finite estimates were preserved.

This is execution evidence, not statistical recovery evidence. The fixture has
30 training residual degrees of freedom and supports up to 33 features, so the
declared shrinkage estimator is load-bearing for invertibility. The separate
500-replication study above carries the estimator-policy evidence.

The exact same 52,416-feature workload initially took 3,476.43 s because every
node revalidated the full frame and support graph. Canonical boundary
validation plus trusted CSR node access removed that work, without changing the
estimand, covariance recipe, solve, or reduction order. Against the 2026-08-20
re-record above, complete analysis is 101.461 s (34.3 times faster) and the
evaluation stage 57.859 s (59.1 times faster). The 2026-08-17 re-record read
59.982 s (58.0 times faster) and 41.069 s (83.2 times faster); as issued on
2026-08-15 the same comparison read 123.498 s (28.1 times faster) and 57.790 s
(59.1 times faster). The spread across these three records is dominated by host
load, not by the source tree, so the order-of-magnitude claim is the durable
one and the specific multiplier is not.
The 3,476.43 s and 3,417.17 s pre-optimization figures are the only
numbers in this section not held by a shipped artifact: they were measured
once, before the gate runner existed, and no `.rds` records them. Oracle suites
remained green. This is the package's performance constitution in executable
form: remove mathematical and administrative work; never buy speed by changing
the scientific answer.

## Statistical estimands and physical lowerings

| Public result | Estimand | Executed lowering(s) |
|---|---|---|
| `materialize_geometry()` / `evaluate_geometry()` | symmetric cross-generalized effect form or fixed query under explicit neural metric | additive form/query-fused contraction for diagonal metrics; support-streamed metric form/query contraction for non-diagonal metrics |
| `crossnobis()` with fixed `noise_precision()` | signed crossvalidated squared Mahalanobis contrast | geometry query through the same fixed-metric compiler |
| `plan_crossnobis()` + `crossnobis()` | signed crossvalidated squared Mahalanobis contrast with provenance-frozen learned local precision | `support_streamed_scheduled_metric_query_contraction` |
| `measurement_form()` | requested neural self/cross blocks under a fixed experimental query | `pull_h`, `forward_k`, low-rank `factorized_h`, rank-one `scalar_stack`, or requested multivariate blocks |
| coupling/connectivity views | raw effect coupling, or capability-gated covariance/correlation/CCA/alignment/information summaries | views over the materialized requested blocks; no separate connectivity engine |
| `reconstruct_evidence()` | exact or explicitly projected neural evidence operator | Parseval or canonical-dual reconstruction within dense resource gates |

The compiler chooses only mathematically equivalent physical routes. Scientific
plan identity records the estimand; execution receipts record the actual route,
tiles, source revisions, workspace, and materialization.

## Public surface and scale regimes

The artifact exports 78 functions. *(Superseded: `NAMESPACE` exports 105
functions and registers 207 S3 methods — see the erratum above and "Recorded
certification metrics" below. The semantic grouping that follows is unchanged.)*
They fall into these semantic groups:

- relation and error sources: `relation()`, `relation_fit()`,
  `lm_relation_fit()`, block/capability readers;
- experimental queries and partition pairings: effect spaces, pair queries,
  contrast/control builders, variation queries, and reducers;
- neural domains, frames, metrics, metric recipes, and bridges;
- query-first geometry plans and derived contrast/RDM/RSA/spectrum views;
- learned and fixed crossnobis;
- explicit-edge measurement forms, coupling views, and tomography;
- neuroim2 domain, searchlight, and exact `as_neurovol()` adapters;
- compute, source, and numerical-contract controls.

The supported regimes are intentionally unequal:

- brain-scale identity/diagonal geometry uses the shared sparse additive path;
- brain-scale dense local precision uses support-streamed on-demand solves and
  is qualified only for the named support/fold regime above;
- complete geometry materialization is preflighted, while query-first output is
  the primary large-`q` path;
- measurement-form connectivity is a small-node dense path gated at 256 MiB;
- tomography defaults to a 512 MiB reconstruction ceiling;
- execution is sequential and deterministic in version 0.1.

Print and `format()` methods summarize plans, fits, geometry, and primary views
without hidden source reads. `as.data.frame()` is lossless for already
materialized query, contrast, RDM, RSA, spectrum, and crossnobis views. Complete
geometry still requires an explicit component read, so printing never triggers
large materialization. Runnable Rd examples cover query-first geometry, fixed
crossnobis, and the raw-response learned-metric journey.

## Subtraction receipt

The architectural reset removed the unreachable public edge-operation surface:

- exports removed: `correlation()`, `cosine()`, `covariance()`, `fisher_z()`,
  and `rank_edges()`;
- evaluator removed: `.evaluate_edge_operation()`;
- tests removed or shortened where they exercised only that unreachable
  evaluator;
- retained because they are executable: `inner_product()`,
  `aggregate_first()`, and `reduce_partitions()`.

The deletion was scientifically correct, but the pre-deletion files were never
committed and are not recoverable from Git. This receipt reconstructs the
deleted surface and rationale from the tracker review; it cannot reproduce the
exact removed text. That provenance loss remains a process defect and is one
reason this report refuses to call the current workspace a release commit.

## Explicitly unrun or deferred evidence

- No hosted GitHub Actions workflow or hosted scale tier was run.
- No Windows, Linux, Intel macOS, CRAN, R-universe, or downstream-consumer
  court was run.
- The final CRAN incoming metadata court was network-qualified; the complete
  local package court is clean.
- No public-data matched-estimand analysis was run.
- No LD-t, confidence interval, bootstrap, permutation, or group-level
  calibration is implemented.
- No classifier-accuracy, nonlinear learner, dynamic informational-
  connectivity, distributed executor, or arbitrary plugin registry is claimed.
- Full brain-scale connectivity/tomography is not claimed; those exports are
  explicitly resource-gated small-node facilities.
- Surface output follows the exact index contract conceptually, but this slice
  added and scale-qualified only the NeuroVol output adapter.

## Closure decision

The local implementation objective is achieved: the universal evaluator now
supports a statistically meaningful, query-first, noise-normalized searchlight
journey at realistic volume scale, with exact results, bounded storage,
independent oracles, generative evidence, and honest capability boundaries.

The implementation epic can close. Release certification cannot: it remains
blocked on a user-authorized commit and, later, hosted/cross-platform and
external scientific evidence.

## Certification provenance binding (2026-08-16)

This receipt's evidence sections cite recorded `.rds` artifacts. Until today
the test suite read those artifacts as bare booleans
(`expect_true(artifact$gate$passed)`) with nothing tying an artifact to the
source that produced it, and `benchmark-results/` was `.Rbuildignore`d, so
under `R CMD check` every certification block skipped silently. A stale
artifact passed forever, and the tier the project calls certification was the
tier that did not run where it mattered. Eight of the twelve committed
artifacts, including both large validation records, were read by no test at
all.

**Binding.** Every runner that persists an artifact now records a
`provenance` list into it (`benchmarks/provenance.R`): an aggregate SHA-256
over the sorted per-file digests of `R/*.R`, the git commit and whether `R/`
was dirty, the package version, and the R/platform/BLAS/LAPACK environment.
`tests/testthat/helper-certification.R` refuses to read a recorded gate unless
that provenance still binds it, in two tiers. From a source checkout the
recorded source digest must equal the current one. Under `R CMD check` the
package sources are not shipped and `.git` is absent, so only the package
version can be compared; that weaker tier is named rather than presented as
the stronger one. A non-binding artifact **skips with a loud, greppable
message** naming the runner that re-certifies it. Hard failure was rejected:
it would redden the suite on every ordinary source edit, which teaches people
to switch the check off.

**Shipping.** The small gate artifacts moved from `benchmark-results/` to
`inst/extdata/certification/`, so certification now runs under check instead
of skipping. `benchmarks/promote-artifacts.R` is the only supported path from
a fresh run into the package; it refuses an artifact without provenance and
anything above a 64 KiB cap. The two large validation records stay out of the
tarball, and their tests skip with a message when they are absent.

**Re-recorded today**, each passing its own gate as recorded, on R 4.5.1,
aarch64-apple-darwin20, Accelerate BLAS: the four memory scenarios, the
first-moment vertical slice, the sampling-covariance scale record, the
brain-scale crossnobis gate, the public map gate, the query-first gate, and
the sampling-covariance validation. Two artifacts remain unbound and their
tests skip with instructions: `shard-admission.rds` (needs an isolated
library with `shard` and `crossform` installed; the recorded verdict remains
"not admitted") and `learned-metric-policy-validation.rds`.

The re-recorded artifacts stamp `provenance$recorded_at` between
`2026-08-17 04:39:40 UTC` and `2026-08-17 05:02:07 UTC`; this section is dated
by the local working day on which they were run. The stamps, not the heading,
are what the tests bind against.

**Two recorded claims did not survive re-execution and were corrected rather
than re-asserted.** The memory harness documentation described a conservative
plan carrying a named 64 MiB process-level runtime reserve, and required every
scenario to report `plan_covers_incremental_peak = TRUE`. The current
`memory_plan()` emits no such reserve and no such column: its `prediction_kind`
is `crossform_owned_workspace_upper_bound`, a claim about crossform-owned
workspace and explicitly not about process RSS. The committed artifacts were
schema 2 against a schema 3 emitter. `benchmarks/README.md` now describes what
the code does, and the tests check the plan against what it actually bounds
while checking measured RSS separately.

**Court.** `testthat::test_local()`: FAIL 0. `R CMD check --as-cran
--no-manual` on the built tarball: tests report `FAIL 0 | WARN 0 | SKIP 10 |
PASS 4112`, with the certification blocks running rather than skipping. The
remaining skips are the two unbound or absent validation records, the
`CROSSFORM_RUN_SCALE_TESTS` tier, and three `On CRAN` blocks. This 4,112
supersedes the `PASS 1822` block under "Numerically verified" above, which is
the 2026-08-15 source-checkout run.

## Re-certification after WS-B (2026-08-20)

Every artifact above was re-recorded against the post-WS-B tree. The frozen
source is `git rev-parse HEAD = 01e78b565b46cac0d724e1541aae7ed159cd5686`
(branch `elite-pass`), source digest
`sha256:346aba47b81e768edc36e3fbbc51fddfff8220e3bc8522c778e31cf02a917ede`,
superseding `sha256:67f21604d0f0…` (2026-08-17) and `sha256:2ab5640f8e25…`
(the 2026-08-20 population-null-coverage record made mid-WS-B). `R/` and
`src/` were `git`-clean for the whole sequence and the digest was re-read
afterwards unchanged, so all twelve artifacts bind to one source state.

All nine persisting runners were re-run in the `benchmarks/RECERTIFY.md` order
and each passed its own gate. `benchmarks/promote-artifacts.R` promoted nine
`.rds` artifacts and ten summary CSVs with no refusals; the two
`shard-admission` files were skipped as not present, which is the designed
state for a deliberately excluded runner. The E8 population-null-coverage
record ships as its summary CSV only, per its `local` role in
`benchmarks/admission-coverage.R`.

**Three runners had to be repaired before they could gate at all.** Ticket B3
(`09e18c4`, "Retire effect_crossnobis_plan and the private learned driver")
made `plan_crossnobis()` return an `effect_geometry_plan`, which moved six
fields the benchmark harness read. Nothing in `R/` was touched to fix this;
the staleness was entirely in `benchmarks/`:

| Stale read | Current location | Effect if unrepaired |
|---|---|---|
| `diagnostics$pair_atoms_materialized` (and `pair_frame_`, `metric_factor_table_retained`) | `diagnostics$total$…` | `!NULL` aborts the crossnobis gate with "invalid argument type" |
| `plan$memory$planned_workspace_bytes`, `…$budget_bytes` | `value$receipt$memory$…` | summary `data.frame()` aborts on a zero-length column |
| `plan$kernel_version` | `value$receipt$kernel_version` | same abort; also aborts the learned-metric runner |
| `plan$metric_schedule$statistics` | `plan$metric_schedule$schedule$statistics` | residual reads silently recorded as `0` |
| `materialize_metric(plan$metric_schedule, …)` | `…$metric_schedule$schedule` | package refuses the unsealed object |

That the repairs are faithful rather than convenient is checkable: the rebuilt
crossnobis artifact reproduces the 2026-08-17 record's residual reads (11,163)
and planned workspace (281,589,160 bytes) exactly, and the query-first record
reproduces its planned workspace (149,400,790 bytes) exactly. No threshold was
lowered and no gate assertion was removed; the one plan-level check that could
not be restored — equality of `kernel_version` across policy plans, which the
user-facing plan no longer carries — is subsumed by the per-result receipt
check against a literal kernel name that the same runner already performs.

**What did change scientifically.** The 52k crossnobis fixture now records
`execution$lowering = derive_then_support_streamed_pair_contraction` where the
2026-08-17 record read `support_streamed_scheduled_metric_query_contraction`.
The executed kernel is unchanged (`execution$kernel_version =
support-streamed-scheduled-metric-v1` in both), so this is a plan-level
relabelling that follows from the class change, not a different execution
path. The plan identity prefix moved with it: `crossnobis-sha256:` no longer
exists anywhere in `R/`, and crossnobis plans now carry `geometry-sha256:`
ids like every other geometry plan.

## Recorded certification metrics (2026-08-20)

Every current number this receipt asserts, with the artifact or file that holds
it. Superseded values are listed newest-first so each correction stays visible
rather than being silently edited away.

| Metric | Recorded value | Superseded values | Source |
|---|---:|---:|---|
| Exported functions | 112 | 105, 78 | `NAMESPACE` (`grep -c '^export(' NAMESPACE`) |
| Registered S3 methods | 80 | 207, 43 | `NAMESPACE` (`grep -c '^S3method(' NAMESPACE`) |
| — `print` methods | 28 | 101 | `NAMESPACE` (`S3method(print,`) |
| — `format` methods | 22 | 95 | `NAMESPACE` (`S3method(format,`) |
| — `as.data.frame` methods | 12 | 6 | `NAMESPACE` (`S3method(as.data.frame,`) |
| — `plot` methods | 5 | 5 | `NAMESPACE` (`S3method(plot,`) |
| Distinct classes with any S3 method | 31 | 103 | `NAMESPACE` (distinct second argument of `S3method(`) |
| Brain-scale complete analysis | 101.461 s | 59.982 s, 123.498 s | `inst/extdata/certification/crossnobis-scale-gate-summary.csv`, `analysis_seconds` |
| Brain-scale incremental peak RSS | 570,572,800 bytes | 903,479,296, 948,568,064 | same CSV, `incremental_peak_rss_bytes` |
| Query-first selected 100 pairs | 0.142 s | 0.130 s | `query-first-scale-gate-summary.csv`, `selected_median_seconds` |
| Query-first full fused RDM | 1.936 s | 1.586 s | same CSV, `full_fused_median_seconds` |
| Query-first materialize-then-project | 2.734 s | 2.079 s | same CSV, `materialized_median_seconds` |
| Query-first fused/materialized ratio | 0.708 | 0.763 | same CSV, `fused_to_materialized_ratio` |
| Query-first incremental R heap | 257,949,696 bytes | 276,719,200 | same CSV, `incremental_peak_r_heap_bytes` |
| Public map sampling full sweep | 101.864 s | — | `public-map-scale-gate-summary.csv`, `sampling_full_sweep_seconds` |
| Test expectations, source checkout | 8,353 over 125 files | 1,822 (2026-08-15) | `devtools::test()` on 2026-08-20; 8,346 pass, 6 skip, 1 fail |
| Test expectations under `R CMD check` | 4,112 | 1,822 | the 2026-08-16 court above (check log; not a shipped artifact) |

The 80 registrations cover 31 distinct classes. `print` and `format` are
registered for the same class in all 22 `format` cases; 6 classes have `print`
alone and none has `format` alone, so 28 classes carry their own printer. The
drop from 207 registrations over 103 classes is the WS-A print/format
consolidation, not a loss of coverage: the printer surface was merged onto
shared classes rather than deleted. Exports rose from 105 to 112 over the same
period. These are surface-size facts about the package, not certification
passes or failures.

The five query-first rows and the brain-scale timing row were all re-measured
on 2026-08-20 under heavy host load; see the contention caveat under
"Scale-qualified". The vignette quotes in `vignettes/novelty.Rmd` were updated
to match these shipped CSV values exactly.

The 8,353 figure is a source-checkout `devtools::test()` run and is not
comparable to the 4,112 `R CMD check` figure above, which was not re-run on
2026-08-20; the check court counts fewer expectations because the
`CROSSFORM_RUN_SCALE_TESTS` tier and several source-only blocks do not execute
there. Both are recorded so neither is mistaken for the other.

**Court after re-certification (2026-08-20).** `devtools::test()`: 8,346 pass,
6 skip, 1 fail. Every `CERTIFICATION STALE` skip is gone — the eight
certification blocks now execute against artifacts bound to the current
digest, which was the point of the exercise. The six remaining skips are all
accounted for: `shard-admission` (`CERTIFICATION UNBOUND`, by design), the
four `CROSSFORM_RUN_SCALE_TESTS` opt-ins (public map gate, query-first gate,
and the two 50k topology tests), and the population slice 2 cross-fit axis
whose exemplar CSV had not yet been produced when this run started.

Two plan-identity assertions were carried through to B3's new scheme.
`test-certification-artifacts.R` matched `^crossnobis-sha256:[[:xdigit:]]{64}$`
at two sites; B3 retired that prefix, so both artifacts correctly record
`geometry-sha256:` ids, and both sites now match that instead. This is the
format match `crossform-execution-design.md` § "Retiring the duplicate naming
rule" anticipates — the artifact tests match the id *format* and re-derive
verdicts from recorded measurements, never comparing a recorded id to a
freshly computed one — so a strong prefixed SHA-256 identity is still pinned
and nothing was loosened. With that carried through, every certification block
passes.

One unrelated failure remains: `test-population-slice2.R:95` is an exemplar
receipt assertion from concurrent slice 2 work. That file reads no
certification artifact.

Three caveats on traceability. The S3 and export counts are counted from
`NAMESPACE` in the working tree, not from a recorded artifact, and will drift
with any export change; the commands above reproduce them. The 4,112
expectation count comes from a check run, and no shipped artifact carries it —
it is the one recorded number here that cannot be re-derived without running
`R CMD check`. The 3,476.43 s / 3,417.17 s pre-optimization timings quoted
under "Scale-qualified" predate the gate runner and are likewise unheld by any
artifact.

## Hot-path epic closure (2026-08-21)

The legacy native-kernel epic was reconstructed in the current Mote store as
`bd-01M0J3FYZ6Q20K1SR1R9DMD8SQ` because its original records were not present
after the tracker-store replacement. The reconstruction audited the landed
commits before adding new evidence: `d1504bd` introduced the fused pair-query
kernel, `2bd7ba3` added the sampling, packed/coherent, and topology work, and
`f6baa14`, `03e1d33`, and `8fdbb17` supplied the prior re-certification and
admission ledger. No legacy child was credited from prose alone.

Two remaining evidence gaps are now closed:

- `measurement-profile.rds` promotes the no-Rcpp decision for
  `measurement_form()`. Across the scalar and requested-multivariate routes,
  the maximum independent-oracle error is `7.10543e-15`, plan identity is
  stable, and the classified R-loop share is zero. The recorded decision is
  `no_rcpp_keep_blas`: another native kernel is not admitted under the declared
  15 percent R-loop / projected 1.25x speedup rule.
- `native-pair-allocation.rds` records a five-repetition cumulative-allocation
  court for the fused pair-query kernel. The native route allocated a median
  4,638,576 bytes against 286,853,280 bytes for the retained two-pass R oracle,
  a ratio of 0.0161706 against the maximum 0.70. Its maximum numerical error is
  `9.43690e-16`; the median runtimes were 0.030 and 0.094 seconds,
  respectively. This receipt concerns cumulative R allocation, not peak heap
  or process RSS.

The source tree remained bound to the certification digest
`sha256:346aba47b81e768edc36e3fbbc51fddfff8220e3bc8522c778e31cf02a917ede`.
The complete source-checkout test suite passed with the six declared skips
(one unbound executor, one absent local-only population artifact, and four
opt-in scale/topology courts). A compiler-neutral `R CMD check` using Apple
clang passed with zero errors and zero warnings. Its single remaining note
combines new-submission metadata with URLs for the not-yet-published site and
files not yet on the hosted `main` branch; it is a release-state note, not a
hot-path or package-correctness failure.
