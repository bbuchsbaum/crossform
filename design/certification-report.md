# Statistical and brain-scale validation receipt

Date: 2026-08-15

Umbrella: `bd-01KZXA4C1RJDV567TS9X175ND0`  
Final evidence epic: `bd-01KZXA62QZQMH79RS7KW3PVCGT`

## Errata / reconciliation (2026-08-15)

This receipt is preserved as issued. Two of its statements do not survive a
reading against the repository, and are corrected here rather than silently
edited above.

**Export count.** "Public surface and scale regimes" states that the artifact
exports 78 functions. `NAMESPACE` exports 105 functions and registers 43 S3
methods. The semantic grouping in that section is still accurate; the count is
not.

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

`benchmark-results/crossnobis-scale-gate.rds` records one exact brain-scale
execution:

| Quantity | Measured value |
|---|---:|
| Active volume features | 52,416 |
| Searchlight radius / spacing | 6.1 mm / 3 mm |
| Support size | min 11, median 33, mean 31.16, max 33 |
| Union-pair stored nonzeros | 4,365,964 |
| Local metric derivations | 52,416 |
| Factorization work units | 1,653,408,544 |
| Residual source reads | 11,163 |
| Evaluation source reads | 104,832 |
| Planned owned workspace | 281,589,160 bytes |
| Incremental peak RSS | 948,568,064 bytes |
| Plan plus residual statistics | 64.348 s |
| Crossnobis evaluation | 57.790 s |
| Exact NeuroVol mapping | 0.084 s |
| Complete measured analysis | 123.498 s |

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
validation plus trusted CSR node access reduced complete analysis to 123.498 s
(28.1 times faster) and the evaluation stage from 3,417.17 s to 57.790 s
(59.1 times faster), without changing the estimand, covariance recipe, solve,
or reduction order. Oracle suites remained green. This is the package's
performance constitution in executable form: remove mathematical and
administrative work; never buy speed by changing the scientific answer.

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

The artifact exports 78 functions. They fall into these semantic groups:

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
`CROSSFORM_RUN_SCALE_TESTS` tier, and three `On CRAN` blocks.
