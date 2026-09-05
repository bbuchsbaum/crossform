# Predictive geometry: test specification

Date: 2026-09-04. Status: proposed acceptance tests; not implemented production
coverage. Companion [work packages](2026-09-04-predictive-geometry-work-packages.md).
The previous coder's passing suite is a baseline for its own estimands. It
does not establish correctness of this regularized predictor or its score.

The test bar is: **a plausible scientific error must fail a specific test**.
Counts, line coverage and golden output alone do not meet that bar.

**Objects and exact reference cases**

Let `S` be the signed compressed geometry; `J = U D U'` is the positive
support of the pooled model kernel; `A = U [U' S U - lambda D^-1]_(+,s) U'`.
The score is `gain(A, S_test) = 2 * sum(A * S_test) - sum(A * A)`.
All tests distinguish this from projector energy and from trace energy.

The executable [fixture file](2026-09-04-predictive-geometry-test-fixtures.R)
contains 14 base-R reference cases and deliberately wrong alternatives.
It does not call crossform. Its [output](2026-09-04-predictive-geometry-test-fixtures.txt)
validates the arithmetic and sensitivity of selected assertions, not a future
implementation. Embed its small matrices into an explicit orthonormal basis
of the centered condition space when exercising public package entry points.

| Fixture | Exact content and use |
|---|---|
| F01 | `S=diag(6,3,-1)`, `J=diag(.5,.25,.25)`, `lambda=.5`; shifted roots `(5,1,-3)`, rank-2 `A=diag(5,1,0)`, objective `10`. Rank path: zero, `diag(5,0,0)`, `A`, `A`. |
| F02 | `S=diag(2,7)`, `J=diag(1,0)`, `lambda=1`; fit `diag(1,0)`. A pseudoinverse-only shortcut incorrectly fits `diag(1,7)`. |
| F03 | `S=[[1,2],[2,1]]`, `J=diag(.8,.2)`, `lambda=.2`; shifted matrix `[[.75,2],[2,0]]`. Rank-one PSD solution calculated from the quadratic formula and spectral projector, without `eigen()`. |
| F04 | Same full `G=[[1,2],[2,1]]`, compress to first coordinate: signed `S=1`, while premature full PSD clipping gives `1.5`. |
| F05 | `A=diag(2,1)`, `S_test=diag(1,-1)`; inner product `1`, squared norm `5`, gain `-3`. |
| F06 | `A=[[1,.5],[.5,1]]`, `S_test=[[0,1],[1,0]]`; inner product `1`, squared norm `2.5`, gain `-.5`. Catches packed off-diagonal errors. |
| F07 | Amplitudes `(2,1)`, mode evidence `(2,0)`; gains `(4,-1)`. Rank-1 gain `4` exceeds rank-2 gain `3`. Also amplitude `3`, evidence `1` gives gain `-3`. |
| F08 | F05 prediction; coherent test `diag(1,0)`, configuration test `diag(0,-1)`. Evidence adds, total gain `-3`; summing two component gains incorrectly gives `-8`. |
| F09 | Fixed prediction and target with test noise drawn equiprobably from `{N,-N}`. Enumerate the expectation exactly. Under a zero target, F05 prediction gives expected gain `-5`. |
| F10 | Model eigenvalues `(8,2)`, keep rank 1. Unit retained trace gives `diag(1,0)`; normalization by the original trace incorrectly gives `diag(.8,0)`. |
| F11 | Full forms `diag(1,0)` and `diag(1,4)` have identical first-coordinate compression. Prediction `diag(1,0)` has full residual squared `0` versus `16`. |
| F12 | `S1=diag(1,-2)`, `S2=diag(-2,1)`: sum of separate PSD fits is identity, but PSD fit of the sum is zero. |
| F13 | Full-span `S=3I`, kernels `diag(.8,.2)` and `diag(.2,.8)`, `lambda=.4`: fits `diag(2.5,1)` and `diag(1,2.5)`. Both give `3I` when unregularized at rank 2. |
| F14 | Shifted form `2I`, rank 1: two equal optima `diag(2,0)` and `diag(0,2)` have test gains `0` and `-4` against `diag(1,0)`. A rank-cut tie can change the prediction. |

**Reference independence and tolerance policy**

- The independent court must not import `.latent_rank_psd_form()`, package
  packing helpers, kernel pooling, or package fixture constructors. Derive
  small answers analytically and dense cross-partition products directly.
- A second spectral implementation is useful differential evidence, but
  shares possible conceptual errors. Require F01–F14, metamorphic laws and
  dense feature-space checks as independent anchors. Numerical optimization
  of the factor objective is supplementary; nonconvergence is not evidence
  that the spectral result is wrong or right.
- For ordinary well-conditioned fixtures use `abs_error <= 1e-12 +
  1e-10 * reference_scale`; use matrix Frobenius norms with documented scale.
  Exact rank, dimensions, IDs, selected support and refusals use exact checks.
- For stress cases report residual, condition estimate and spectral gaps.
  Set the admissible error before running the test, using backward residuals
  and gap/conditioning where appropriate. Do not enlarge tolerances until
  a result passes. Near-rank thresholds test the declared classification;
  do not demand stable eigenvectors where the mathematical problem is unstable.
- All outputs admitted as numeric must be finite. Rank and PSD tests use
  numerical tolerances separate from scientific rank truncation thresholds.
  Boundary ties use the explicit contract policy, not arbitrary LAPACK signs.
- Seeds are fixed before running: reserve `2026090401` for generated kernels,
  `2026090402` for pullback fixtures, `2026090403` for sampling calibration,
  `2026090404` for selection calibration and `2026090405` for benchmarks.
  Use independent per-replicate streams; worker count cannot change the data.
  Failure output includes the seed, minimal matrices, policies and source SHA.

Tiers below: **F** = fast mandatory PR; **D** = dense differential PR;
**S** = scheduled stochastic/stress; **B** = isolated benchmark/release.
Times are budgets to measure, not claims about tests that do not yet exist.

**Kernel admission and coordinate meaning**

Suggested file: `tests/testthat/test-predictive-kernel.R`.

| ID | Tier | Required assertion | Error it must expose |
|---|---|---|---|
| T01 | F | A small named feature matrix, its squared Euclidean RDM, and its centered Gram give the same retained **kernel**, fitted form and score after identical normalization. Compare more than spans. | Silently discarding model eigenvalues; confusing distances with squared distances. |
| T02 | F | Independently permute labelled rows and columns, then align. Duplicate, missing or one-sided labels refuse before neural reads. Positional admission requires an explicit same-size coordinate declaration. | Correct shapes with incorrect condition identity. |
| T03 | F | Nonfinite, nonsquare, materially asymmetric or materially indefinite kernels refuse with structured reasons. Numerical-scale negative roots can be repaired only under the recorded policy. All-zero centered models refuse. | Silent model replacement, NaN propagation, asymmetric eigen input. |
| T04 | F | Effect coordinates, shared units/scales, centering origin and distance convention are checked. A nonorthogonal coordinate scaling changes the declared Frobenius problem; do not demand false invariance. | Treating every same-size effect array as the same estimand. |
| T05 | F | F10: truncation followed by retained-trace normalization gives unit retained trace, with raw and retained scales recorded. | Wrong normalization order. |
| T06 | F | Independently rescale each model by positive factors from `1e-8` to `1e8`: with retained-trace normalization and unchanged numerical rank, pooled kernels and predictions are unchanged. | Model units covertly changing regularization. |
| T07 | F | Above-rank requests clamp and record requested/effective ranks; invalid model ranks refuse. Zero response rank is valid. A zero model is not silently replaced by a zero predictor. | Confusing model-rank admission with response-rank zero. |
| T08 | F | Duplicate one normalized model and split its weight: pooled kernel, fitted form and score remain equal. Rectangular `R` is admitted. Include overlapping, independent and near-overlapping spans. | Requiring identified redundant coefficients or counting shared directions twice. |
| T09 | F | F02 and weights `(1,0)`/`(0,1)`: positive pooled support controls the fit. All-zero, negative, nonfinite or malformed named weights refuse; valid simplex vertices work. | Fitting null-space directions; silently repairing invalid weights. |
| T10 | D | Add a condition-independent neural baseline and a constant offset to feature rows. Centered model predictions and compressed forms remain unchanged. Retain the existing `1e-4` feature-scale centering regression. | Baseline leakage through poorly conditioned eigenvectors or excess QR columns. |
| T11 | F | A model truncation splitting a positive tied eigenspace follows the declared refusal policy; retaining the full tied group is allowed. Rotate that retained group and recover the same kernel. | Sign fixing mistaken for full eigenbasis identification. |

**Estimator and rank path**

Suggested files: `test-predictive-fit-numerics.R` and a base-R court under
`design/oracles/predictive-geometry.R`.

| ID | Tier | Required assertion | Error it must expose |
|---|---|---|---|
| T12 | F | F01 gives exact shifted roots, fitted form, effective rank and objective; distinguish original negative mass `1` from shifted negative mass `3`. | Penalty sign/factor mistakes; misleading neural moved-mass labels. |
| T13 | F | F03 agrees with the analytic 2x2 formula; rotate the coordinate system and obtain the same lifted form. | Diagonal penalization in an arbitrary basis; assuming source and kernel commute. |
| T14 | F | F04 fits from signed compression, yielding `1` before any model penalty, not `1.5`. | Premature full-neural PSD projection. |
| T15 | F | F03 differs from the result of clipping `S` before subtracting its anisotropic penalty; require the correct analytic answer. | Premature compressed-neural PSD projection. |
| T16 | F | F01's complete rank path, zero and all-negative sources, `s=0`, `s=r`, and oversized budgets follow declared clamping/zero rules. Zero predictor has exactly zero score. | Empty-index bugs, negative-root retention, fabricated modes. |
| T17 | F | At fixed kernel/penalty, optimal objective is nonincreasing as rank budget expands. At fixed rank, optimal penalized value is nondecreasing in penalty and `tr(J^+ A)` is nonincreasing. Explicitly do not assert monotone test gain. | Wrong spectral ordering; testing a false property of held-out performance. |
| T18 | D | Orthogonal coordinate changes and effect permutations preserve the lifted prediction. Scaling `(S,lambda)` by `c>0` scales `A` by `c` and jointly scaled test gain by `c^2`; scaling unnormalized `(J,lambda)` by `c` preserves `A`. | Wrong units, orientation-dependent results, penalty-scale confounding. |
| T19 | D | For arbitrary admissible candidate `A`, verify the completed-square equality including its constant and compare full versus compressed objective differences. | Claiming an identified absolute full-space residual from compression. |
| T20 | D | For redundant `T`, independently compute a minimum-norm factor and verify `T W W' T'=F` and `||W||^2=tr(K^+ F)`. Tiny multi-start factor optimization cannot beat the analytic optimum and should reach it on a pinned nondegenerate case. | Wrong regularizer equivalence; treating redundant coordinates as unique. |
| T21 | F | F13: zero-penalty full-span fits coincide; positive-penalty fits differ. Generic baseline and model-specific interpretation are labelled separately. | Blanket full-span refusal or geometry-blind fitting. |
| T22 | F | F14 refuses an ambiguous positive rank-cut tie. Ties wholly retained yield invariant fitted forms and grouped score; individual mode orientation is not asserted unique. Include near-tie cases on both sides of the tolerance. | Platform-dependent prediction silently presented as unique. |
| T23 | S | Generate centered kernels with controlled condition numbers and gaps, varying `n`, admitted rank, scale and overlap. Check finiteness, symmetry, PSD, rank/support residuals and dense-reference objective. Log any deliberately refused unsafe cases. | Numerical overflow, unstable inverse use, unjustified universal tolerances. |

**Lowering, training support and compatibility**

Suggested files: `test-predictive-lowering.R`,
`test-predictive-training-support.R` and `test-predictive-score-contract.R`.

| ID | Tier | Required assertion | Error it must expose |
|---|---|---|---|
| T24 | D | Explicit dense `sym(B_a M_x B_b')` over a small nonuniform pairing matches the package's lowered forms for total/coherent/configuration. Include nonidentity/non-square extractors, fixed admitted metrics and regional/searchlight supports. | Wrong extractor composition, pairing normalization or component lowering. |
| T25 | F | Instrument read/materialization boundaries: fitting requests only reduced form width, and no original-space trace bank. A deliberately forbidden full-form path raises if entered. | Numerically correct results obtained by an unintended quadratic allocation. |
| T26 | F | A parent relation carries unused partitions whose value readers raise. Training and its actual dependency record use only declared training endpoints; metadata inspection is allowed. | Fitting all parent partitions, or mistaking stored-but-unused data for training use. |
| T27 | F | Declare raw observation origins, transformations and derived partitions. Overlapping actual origins refuse even if endpoint names differ; distinct disjoint origins admit under the required assumption. | Treating disjoint names as disjoint observations. |
| T28 | F | Copied/renamed/repacked sources retain known origin ancestry and refuse as independent evaluation. Conversely, two explicitly independent origins with identical numeric values remain distinguishable. | Using content hashes as either a complete independence proof or an automatic dependence proof. |
| T29 | F | Changed training origins, source revisions, transforms, selected hyperparameters or frozen prediction fields invalidate the relevant fit record. Execution knobs remain separate from scientific semantics. | Stale or tampered model admitted by a cached identity. |
| T30 | F | Both within-test cross-product independence and train/test independence must be declared/admitted. Three partitions and a four-partition star pairing expose the insufficiency of count-only checks. | One correct independence axis masking failure of the other. |
| T31 | F | Different training/test pairings are allowed for a matching target. Changed units, effect meaning, frame/measurement identities, metric or component refuse as appropriate. Same dimensions and different plan IDs are not the criterion. | Either rejecting valid held-out plans or accepting scientifically different targets. |
| T32 | F | Change only evaluation outcomes while holding origins and declarations fixed: the learned form and all selections are unchanged, while scores can change. Compare prediction semantics, not necessarily the parent family's hash. | Evaluation data entering fitting through preprocessing, caching or parameter selection. |

**Frozen scoring and result records**

Suggested files: `test-predictive-score-algebra.R` and
`test-predictive-result.R`.

| ID | Tier | Required assertion | Error it must expose |
|---|---|---|---|
| T33 | F | F05 returns inner product `1`, norm squared `5`, gain `-3`, with the factor of two and scalar offset explicit. | Returning alignment/projector energy under a predictive-gain name. |
| T34 | F | F06 matches dense, packed and single-operator query routes. | Missing sqrt(2) off-diagonal scaling or double-counting it. |
| T35 | F | F07's overestimated mode gives negative gain despite positive evidence. Negative signed test evidence remains negative. | Clipping evidence or assuming every replicated mode is beneficial. |
| T36 | F | F07 preserves both mode gains `(4,-1)` and total `3`; no score-time pruning or rank reselection occurs. | Selecting modes on the final test set. |
| T37 | F | F08 evidence terms add and the prediction cost appears once; F12 shows why independently fitted components cannot be claimed to conserve. | Double-subtracted score offsets or nonlinear conservation claims. |
| T38 | F | Enumerate F09's two-point noise distribution and match the conditional risk identity exactly, for zero and nonzero targets. | Mistaking a Monte Carlo coincidence for the score's algebraic law. |
| T39 | F | F11 demonstrates that equal compressed data can have different full residuals. Output uses `prediction_norm_sq`, subspace residual and declared units; no unidentifiable explained fraction is emitted. | Trace/Frobenius confusion or fabricated denominators. |
| T40 | F | Permute evaluation rows and realign by stable measurement IDs; scores match the original order after joining. Duplicate or missing IDs refuse instead of positional pairing. | Applying a searchlight's learned operator to a different searchlight. |
| T41 | F | Serialization/reopening retains prediction, model declaration, actual support and score. Mutate rank, amplitudes, dimensions, offsets and index mapping; validators reject inconsistencies. | Plausible but internally contradictory fitted objects. |
| T42 | F | Byte/content checks show signed inputs are unchanged; scoring does not call the fitting or selection routines. Fitted forms have no raw-estimation inheritance that enables unsupported inference methods. | Silent mutation, refitting, or false capability inheritance. |

**Storage, statistical calibration and validation orchestration**

Suggested files: `test-predictive-storage.R`, `test-predictive-folds.R`,
`test-predictive-selection.R`; scheduled runner
`benchmarks/run-predictive-geometry-validation.R`.

| ID | Tier | Required assertion | Error it must expose |
|---|---|---|---|
| T43 | F | A guarded block store accepts only bounded row reads. Fit/score complete without reading all rows; retained predictions use reduced factors, not a dense redundant-coordinate array. | An ostensibly blocked API that eagerly materializes its store. |
| T44 | D | Memory/block and every admitted worker route agree; vary row/feature block sizes. Inject read errors and interrupted writes: incomplete fits/scores cannot carry complete receipts or be reopened as complete. | Storage-dependent answers and false successful completion. |
| T45 | S | Independent Gaussian partition cross-products under a zero target: across independent datasets, `mean(gain + ||A||^2)` matches zero with a prespecified MCSE criterion. Compare projector energy's different null expectation. | Null gain incorrectly centered at zero; shared training/test noise. |
| T46 | S | Known planted effect matrices yield dense `G*`; test `gain - (||G*||^2 - ||G*-F||^2)` across signal, noise, anisotropic features, heterogeneous partition variance and admitted nonuniform pairings. | A score correct only for homoskedastic identity-metric toy data. |
| T47 | S | Negative-control simulations reuse the fitting geometry or introduce known train/test dependence and show the anticipated optimistic bias. Known lineage violations refuse; hidden dependence is labelled a violated assumption, not claimed detectable from values. | A validation story that cannot distinguish leakage from honest prediction. |
| T48 | F | The edge schedule uses only the declared pairing, with exact orientation/weight handling. Include complete, subset, star and disconnected pairings; unsupported folds refuse before execution. | Inventing undeclared training edges or accepting partitions merely by count. |
| T49 | D | Manual per-fold fits and scores match cross-fit results, including unequal weights and the weighted sum of per-fit squared norms. It is generally not the squared norm of the average fit. | Training once on all data, or applying the wrong aggregated offset. |
| T50 | D | Eight partition fixture: outer train six/test two, inner train four/validation two. Trace every fit, candidate selection and final refit; none uses outer-test observations. Invalid nested schedules fail before data reads. | Reusing the final test as inner validation; assuming four partitions always suffice. |
| T51 | F | Candidate tuples and fold weights are explicit; aggregate validation gain selects a candidate. Exact score ties use a pinned rule: lower rank, then larger penalty, then stable canonical candidate ID. Changing presentation order does not select a different model. | Hidden tuning defaults or iteration-order-dependent selection. |
| T52 | S | Repeat the full inner-selection/outer-score pipeline under null and planted targets. Calibrate against each outer fit's conditional target risk. A deliberately test-selected best candidate/mode is a negative control, not the reported estimator. | Correct fixed-parameter tests masking biased automatic selection. |
| T53 | F | Read/eigensolver counters show reuse for compatible penalty/rank paths. Changes to source revisions, folds, effects, component, normalization or kernel invalidate the correct cache entries. Same-parent held-out value changes cannot alter training products. | Stale reuse across folds or recomputing neural geometry for every candidate. |

**Scientific usefulness, product contract and release evidence**

Suggested files: `test-predictive-example.R`, existing architecture/API tests,
`test-predictive-certification.R`; isolated benchmark runner
`benchmarks/run-predictive-geometry-scale.R`.

| ID | Tier | Required assertion | Error it must expose |
|---|---|---|---|
| T54 | D/S | An analytic same-span example isolates kernel geometry at a common positive penalty. A prespecified simulation compares zero, same-span isotropic shrinkage, supplied kernels and a mismatched kernel using identical selection budgets. Report paired gain differences and MCSE, without claiming universal model superiority. | Advertising generic rank/shrinkage benefit as model-specific explanation. |
| T55 | F | Unsupported condition transfer, covariance-weighted loss, unadmitted learned metrics and generic eigenvalue inference refuse with a precise capability reason. Supported fixed-split scores still work. | Quietly broadening the claimed estimand or solving GLS with the spectral shortcut. |
| T56 | B | Benchmark tiny, medium and large declared workloads with read counts, actual packed width, incremental RSS, wall time, retained bytes and matched numerical parity. Run heavy runners serially. | An unmeasured performance claim or hidden full-form allocation. |
| T57 | F/B | Under a known workspace limit, select admitted blocks or refuse before expensive reads. Requested unsupported worker/storage paths refuse rather than silently switching statistical semantics. | Budget bypass, allocation failure masquerading as numerical failure. |
| T58 | F | Register new files; the layer graph remains acyclic with no upward calls. Fit records and their validators can be consumed by scoring without a fit/score mutual dependency. | Architectural complexity concealed by a short public API. |
| T59 | F | A migration map connects each changed old contract assertion to its replacement. Preserve source-lowering and signed/latent laws. Revise overlap/saturation tests where scientifically justified; do not freeze the old export count or four-structure taxonomy as goals. | Letting partial implementation dictate the new design, or deleting useful guarantees wholesale. |
| T60 | F | The minimal example executes and prints the gain's meaning, units, rank and training/evaluation scope. Fitted amplitudes and replicated evidence are separate. Public exports, generated docs and API tiers agree. | An elegant internal design with a misleading or unusable user workflow. |
| T61 | B | Fresh validation/performance receipts bind to the final source, oracle/test harness, seeds, configurations and runtime. Existing certification gates are re-run on the frozen source; no unexplained STALE/UNBOUND result is accepted as passing evidence. | Shipping old certification as evidence for a new estimator. |
| T62 | F/S | Deliberately wrong variants each fail their named primary test; record the kill matrix. Apply mutations in disposable test-local variants, not shared source files. Every survived critical mutant triggers a test improvement. | Many green tests that would also accept a scientifically wrong algorithm. |

**Mutation kill matrix**

The exact fixture script already distinguishes 16 wrong alternatives. The
following production-path mutations must additionally be exercised after the
corresponding code exists; do not call the fixture check a production mutation run.

| Mutation | Required detecting test |
|---|---|
| Replace `D^-1` by `D`, or add rather than subtract the penalty | T12 |
| Fit the null space after subtracting `J^+` | T09 |
| Clip the full or compressed neural form before fitting | T14, T15 |
| Normalize by pre-truncation trace | T05 |
| Drop sqrt(2) packing weights | T34 |
| Omit the factor two or use trace instead of norm squared | T33 |
| Omit or subtract the prediction norm twice | T33, T37 |
| Pick the best modes using test evidence | T36 |
| Treat endpoint labels/content hashes as complete origin identity | T27, T28 |
| Train on all parent partitions or all outer folds | T26, T32, T50 |
| Use positional measurement matching | T40 |
| Reuse cached training data after a fold/source change | T53 |
| Fit full geometry while returning only a reduced result | T25, T43, T56 |

**Simulation precision and CI placement**

Do not make PR success depend on a noisy expectation estimated from 30
replicates. F09 is an exact deterministic expectation test. Keep small
simulation smoke runs for execution, not calibration claims.

For scheduled T45/T46, initially prespecify 20,000 independent datasets per
small-matrix arm, processed serially in batches; for full nested-selection
T52, prespecify 2,000 independent datasets with the candidate grid and
splits recorded. The unit for MCSE is a generated dataset after any dependent
within-dataset folds are aggregated. These are proposed budgets, not existing
certification requirements. Pin the final arm list and MCSE decision rule
before producing release receipts; never select a favorable seed or keep
adding replicates only until significance disappears.

For identity calibration, define the squared-geometry reference scale
`b = ||G*||_F^2 + E[||G_test - G*||_F^2]` from the known generator, not from
the fitted result. Validate that generator moment independently. If `b=0`,
use the exact zero fixture. Otherwise use the proposed practical margin
`delta = 0.02 * b` and record the normalized errors as well as raw units.
Require both `abs(mean(error)) <= 5 * MCSE + numerical_tolerance` and
`abs(mean(error)) + 5 * MCSE <= delta`, with a familywise rule frozen for the
arm list. The first criterion checks systematic bias; the second requires
enough precision for the declared equivalence margin.

Use an independent, separately seeded pilot only to select the final fixed
replication count before the certification run, targeting `5 * MCSE <=
delta / 2`; retain the pilot and planning calculation. The counts above are
initial budgets. If an arm needs more, increase its declared count before
freezing the production run rather than silently weakening the margin.
Insufficient precision at the frozen count is **not a pass**. Negative
controls have prespecified detectable effects and do not share the
calibrated-estimator pass criterion. Never use per-fold standard errors as
the Monte Carlo uncertainty of independent datasets.

PR: F cases and bounded D cases, targeting under two minutes incremental
runtime after package load. Scheduled: generated stress and stochastic
calibration, with seeds and minimal failures retained. Release: all admitted
routes, full R tests/R CMD check, and isolated benchmarks/recertification.
Do not require a worker backend that the feature does not support; require
a refusal test for it. Do not mistake absence of an optional data cache for
a failed numerical test.

Benchmark proposals: tiny `(n=12,r=4,partitions=4,features=64,nodes=32)`;
medium `(64,12,8,2048,2000)`; large `(200,20,8,20000,10000)` with bounded
deterministic sources and local supports. Full-form comparison is required
where it fits the declared budget; the large arm compares sampled dense
references plus allocation/read invariants instead of forcing a multi-GB
reference allocation. Record hardware and competing load. Set an initial
512 MiB incremental workspace budget for the large blocked arm, then prove
admission from a byte model before running; revise the workload/design if
its necessary buffers exceed that budget. It is not a retroactive RSS ceiling.

**Done means**

Every supported behavior has an oracle/contract assertion and a failure mode;
every admitted execution route has parity and budget evidence; every
prediction claim has an explicit training/evaluation boundary; and every
reported result is tied to the code and test configuration that produced it.
The old 10,433-pass count remains historical evidence for the old tree.
