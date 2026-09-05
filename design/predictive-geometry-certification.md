# Predictive geometry implementation and certification

2026-09-05. G01–G22 implementation and associated T01–T62 verification are
complete. This is local source and package evidence, not hosted CI or
publication evidence.

The public workflow is `model_basis()` → `fit_geometry()` → `score_geometry()`.
Models retain their eigenvalue preferences through the inverse-kernel penalty.
Training uses signed reduced geometry; scoring reads a frozen fitted form on
independent declared observation origins. The existing `rsa()` and descriptive
`model_geometry()` estimands remain separate.

**Source identity and scope**

- Branch: `elite-pass`; baseline HEAD `1d93c857ec240154541e4c6283fa5542912e0f8a`.
- Final R source digest:
  `sha256:10bfc4506bcf9c110f8ac10a0bfba425fb0f3eb4c5975e608bd1cd095690bf12`.
  The digest remained unchanged through performance, mutation, integration and
  final statistical recertification. The working tree is intentionally uncommitted.
- Baseline archive:
  `/private/tmp/crossform-predictive-baseline-dz2i2luh/baseline.tar.gz`.
  Pre-existing unrelated edits were preserved. No commit, push or publication
  was performed.
- All 22 work packages and all 62 test obligations are covered by the
  [implementation state](../.planning/2026-09-04-predictive-geometry-implementation.json)
  and [requirement evidence ledger](predictive-geometry-test-evidence.md).
  The [contract](predictive-geometry-contract.md) defines the admitted estimand;
  the [migration map](predictive-geometry-migration.md) records old-contract decisions.

**Delivered behavior**

The typed model basis admits explicit kernels, features or declared squared
Euclidean distances, aligns condition identities, normalizes retained kernels
after truncation, admits overlapping spans and restricts zero-weight mixtures
to their positive pooled support. Fitting solves the exact regularized PSD
rank problem in reduced coordinates. Zero response rank, ambiguous positive
tie cuts, saturation and numerical negative mass have explicit contracts.

Fit and score records retain invariant forms through bounded factors and sealed
provenance. Memory/block execution, actual source support, independent origin
admission, frozen mode/group evidence, weighted edge cross-fitting, nested
selection and isolated rank-path caches are tested. The public surface remains
small; cross-fit/selection orchestration is private in this version.

**Verification**

| Evidence | Observed result |
|---|---|
| Full source suite, `NOT_CRAN=true` | 11,851 passing assertions; 0 failures/errors/warnings; 7 documented skips |
| Added final infrastructure regressions | 3 oracle-binding assertions and 15 frozen-count recovery assertions pass |
| Predictive tests included in the full suite | 1,411 passing assertions |
| Base-R mathematical oracle | 14 fixtures, 16 wrong alternatives, dense partition/minimum-norm/completed-square checks pass |
| Fixed-parameter calibration | 87,000 independent datasets, 12 bias/precision gates and 4 generator-moment gates pass |
| Complete selection calibration | 2,000 null + 5,000 signal datasets; both risk-identity bias/precision gates pass |
| Production mutations | All 20 killed by named primary tests; baseline/no-op controls pass and every changed function executes |
| Performance | Seven serial routes pass dense form/spectrum/mode/gain parity; large blocked form error below 3.3e-14 |
| Final compact/binding/coverage tests | Pass, including rejection of an edited oracle at unchanged R source and refusal of stale promotion |
| Existing certification | Nine serial runners pass; summaries agree with artifact fields |
| R/Python external parity | Pinned environment passes; regenerated after an LF-only run-receipt writer fix |
| Documentation | Executable example, rendered predictive guide and reference pages pass; package-check vignette/example checks pass |
| R CMD check, final Apple compiler | 0 errors, 0 warnings, 1 incoming-review note; installed tests 9,538 pass, 0 fail/warn, 119 designed skips |
| Diff hygiene | `git diff --check` passes |

The final package check uses Apple clang 15 with R's configured C++17 standard
through `/private/tmp/crossform-g22-Apple-Makevars`; no user/global compiler
configuration was changed. The initial completed check had 0 errors, 1 warning
and 1 note: the warning came from an R-header pragma under the user's Homebrew
clang 20. An intermediate isolated build omitted R's default C++ standard and
was corrected before the final completed check. No production-source change
was required.

The remaining CRAN incoming note covers new-submission/development-version
metadata, optional packages available from the declared R-universe rather than
mainstream repositories, and web URL checks reporting unavailable/redirected
hosted pages. The package is not claimed to have a completely clean CRAN
incoming review or a published documentation site. Local documentation,
examples, installation, code checks and vignette rebuilds pass. The final
installed run deliberately skips 119 source-only/optional checks; these are
not counted as source-bound evidence.

The seven full-suite skips are the designed unbound shard refusal, the opt-in
48-seed matched-simulation rebuild, public-map and query-first live reruns, the
review-bundle rebuild and two 50k topology opt-ins. The public-map/query-first
runners were executed separately in this certification batch. None of these
skips is counted as a passing experiment, and no `CERTIFICATION STALE` skip
remains.

The fixed-split calibration uses pilot-derived counts 20,000/27,000/20,000/20,000
for null/aligned/outside-model/heterogeneous-metric arms and the same fixed
rank/penalty budget for supplied, isotropic and mismatched model kernels.
Production seeds are 2026090403 (fixed parameters) and 2026090404 (selection),
with the declared per-arm offsets. Pilot seeds/counts are separate and frozen.
All risk checks satisfy both the five-MCSE criterion and the prespecified
2% reference-scale equivalence bound. Deliberately leaky controls detect
optimism; their scores are not reported as valid model evidence.

`benchmarks/predictive-geometry/restore-counts.R` restores only the frozen count
contract from the shipped receipt when local pilot records are absent. It
verifies model-independent harness/configuration hashes and never re-pilots
from the observed production results. Existing conflicting counts refuse.

**Compact and complete records**

The shipped receipt is
[`predictive-geometry-validation.rds`](../inst/extdata/certification/predictive-geometry-validation.rds).
It retains statistical summaries/configurations/seeds/counts, paired model
comparisons, measured performance and reuse counters, mutation positive controls
and hashes, the oracle output, and all relevant source/harness/test bindings.
Its final serialized size is 27,476 bytes, below the 64 KiB shipment cap. Full Monte Carlo arrays and mutation
bodies/logs remain under `benchmark-results/` and can be regenerated with
[RECERTIFY.md](../benchmarks/RECERTIFY.md).

The bare-R binding gate checks 11 bound shipped receipts (including the external
CSV) against the current tree. Only `shard-admission.rds` is explicitly permitted
to remain unbound. Predictive oracle, test, calibration and promotion-tool edits
also invalidate the compact receipt. Installed tests explicitly report the
weaker package-version binding when a source checkout is unavailable; local
source binding is separate evidence.

The serial batch below records elapsed wall time, not universal timing claims.
The external parity pipeline was rerun after correcting its receipt line
endings; its final log is listed below.

| Job | Exit | Seconds |
|---|---:|---:|
| memory | 0 | 10.48 |
| sampling-scale | 0 | 6.98 |
| sampling-validation | 0 | 8.17 |
| population | 0 | 27.35 |
| first-moment | 0 | 15.16 |
| public-map | 0 | 55.56 |
| query-first | 0 | 25.69 |
| crossnobis | 0 | 50.10 |
| learned-metric | 0 | 222.02 |
| predictive-sampling | 0 | 62.43 |
| predictive-selection | 0 | 34.87 |
| external-parity | 0 | 43.95 |

**Performance and scientific limits**

The large case has 200 conditions, model rank 20, 8 partitions, 20,000 features
and 10,000 measurements. Block execution took 65.949 seconds in its recorded
worker, used packed width 210, planned 145,792,120 numerical-buffer bytes and
observed 1,335,558,144 incremental process-RSS bytes. The 512 MiB workspace
setting is a logical numerical-buffer admission limit; it is not a hard
process-memory cap. R heap/runtime overhead matters.

The medium six-candidate reuse comparison reduces source materializations from
12 to 2 and fitted spectral paths from 12,000 to 4,000 with unchanged outputs.
The dense numerical baseline is faster on these tiny/medium fixtures; the
measured benefits are reduced geometry width, memory and repeated-fit work,
not a universal speedup. Full details and matched-baseline limits are in
[predictive-geometry-performance.md](predictive-geometry-performance.md).

This version evaluates the same conditions/effects and spatial target on
independent declared origins, with an undirected self form and implicit identity
or fixed SPD neural metric. Fixed-metric complete execution is memory-only;
block execution uses the admitted identity route. Learned metrics, generic GLS,
condition transfer and fitted-eigenvalue inference are not admitted. Origin
metadata records an independence assumption; it cannot discover hidden physical
dependence from numerical values. Compression does not identify a whole-geometry
residual or explained fraction.

**Local logs and follow-ups**

- Full suite: `/private/tmp/crossform-g22-full-suite.log` and
  `/private/tmp/crossform-g22-full-suite-results.rds`.
- Final focused courts: `/private/tmp/crossform-g22-certification-court-final.log`,
  `/private/tmp/crossform-g22-integration-laws.log` and
  `/private/tmp/crossform-g22-external-parity-court.log`.
- Binding/promotion: `/private/tmp/crossform-g22-binding-final.log`,
  `/private/tmp/crossform-g22-promotion.log`,
  `/private/tmp/crossform-g22-summary-consistency.log`.
- Serial runner receipts: `/private/tmp/crossform-g22-recertify-results.json`;
  final parity log `/private/tmp/crossform-g22-external-parity-final.log`.
- Package checks: `/private/tmp/crossform-g22-check.log` and
  `/private/tmp/crossform-g22-check-final.log`.

The existing nonidentity-extractor measurement/coupling defect
`bd-01M1Q62WMQVY1S210CH8ZH9J2H` is outside the admitted predictive execution path.
The old Haxby conservative-geometry artifact still needs its absent data cache
for regeneration. Neither is counted as completed by this feature. Hosted
publication and a Git commit are subsequent actions.
