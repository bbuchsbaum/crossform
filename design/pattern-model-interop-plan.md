# Pattern-model interoperability and agent-state implementation plan

Status: proposed implementation, documentation-only revision.  
Contract: `scientific-state-v1`; date: 2026-09-12.  
Bases inspected: crossform `b7cffa9` / rMVPA `7c93f8d`.

## 1. Preserve the package boundary

rMVPA already has pattern fitting, weighted training, regional prediction,
confirmation and loading-level group code at the inspected commit. Crossform's
`elite-pass` has `model_basis`, `fit_geometry`, `score_geometry`, partitioned
relations, error capabilities and population forms. Do not schedule their
reinvention. The implementation gap is a supported interchange and inspectable
state, not another universal estimator.

Crossform remains independent of rMVPA. The optional rMVPA adapter calls public
crossform constructors and queries. No private-class construction, monkey
patching, hidden package installation or new mandatory numerical dependency.

## 2. Mathematical interchange

For the frozen forward relation B_hat = T A', a neural measurement K_j gives
F_j = T M_j T', M_j = A' K_j A. For an experimental query H,
<H,F_j> = <T' H T,M_j>. Read feature blocks as T A[J,]' or contract in rank-sized
coordinates; neither p-by-p covariance nor full condition geometry is required.

With K=Psi^-1, M is the decoder's G. This shared numerical core does not make
quadratic geometry queries equivalent to nonlinear posterior prediction.

Predictions must bind the full raw-target transformation (training means,
scales, whitening, null-space handling and C), identified effect rows, original
neural units, measurement IDs/normalization, metric composition and all training
origins. A learned model evaluated on a new reference feature table is not a
data-independent `model_basis` merely because that table is fixed.

For fixed F and independent unbiased signed G_test, score gain as
2<F,G_test>-||F||^2. Retain negative gain and prediction cost. Under no signal,
nonzero F has negative expected gain. Do not turn this into a p-value or an
explained fraction. New-condition generalization is not licensed by a
same-condition score; cross-kernel/new-condition support needs its own contract.

## 3. Three adapters, implemented in order

### I0. Metadata bridge and status inspection

Add package-local read-only state adapters over existing objects and validators.
Use `agent-protocol.schema.json`; the schema is not a required public R class.
Return source/runtime version, quantities, dimensions/units, training and exposure
ancestry, usable operations, cost/retention limits and explicit refusals.

Owning areas: rMVPA `R/pattern_result.R`, `R/validate_analysis.R`; crossform
relation/plan/fit/score validators and receipt code. Resolve actual symbols before
editing. Do not duplicate the compiler or CV runner. Existing console and
machine-readable views must derive from the same state.

Exit: inspection performs no model fit or source-image read; unknown capabilities
stay unknown; representative existing objects and proposed/refused actions validate
against the schema. Legacy objects get an explicit limited-metadata view rather
than invented independent origins.

### I1. Public frozen geometry-prediction contract

Proposed name: `geometry_prediction()`. Stabilize it only after checking naming
and migration costs. It describes an externally frozen PSD prediction, not a
`fit_geometry` optimization result. Required payload: target identity, a
bounded factor/provider read interface, effect/measurement axes, numerical
rank/tolerance, training/selection dependencies, implementation identity,
read/write policy and integrity checks.

PSD must be guaranteed by a validated factor representation or checked within
budget. Symmetry, finiteness and dimensions are verified without unbounded dense
materialization. Provider code is registered trusted package code, not executable
text read from a manifest. Preserve the existing estimator diagnostics only on
`effect_geometry_fit`; external predictors do not acquire those diagnostics.

Refactor scoring through this small protocol while keeping existing public
`score_geometry(fit_geometry(...), plan)` behavior. In rMVPA, proposed
`as_crossform_prediction()` supplies a factorized predicted relation/form.
Crossform does not call rMVPA to refit anything during scoring.

Exit: existing scores unchanged within declared tolerance; external and native
frozen forms match independent dense oracles; all used inputs/metrics/effect
names match; overlapping origins, changed scales, incompatible normalization and
outcome-selected predictors refuse. No undocumented private-field dependence.

### I2. Confirmation relation export

Proposed rMVPA adapter: `as_crossform_confirmation()`. Prefer exporting retained,
partition-resolved confirmation inputs or explicitly supported effect estimates.
A pooled confirmation result alone cannot reconstruct independent partitions or
a residual channel that was not retained. Missing data yields a capability
refusal or point-only relation, not fabricated residual degrees of freedom.

Freeze task coordinates on discovery rows. In independent confirmation partition
b, estimate X_b=T_b A_b'+E_b with appropriate nuisance/error treatment. Export
B_b=A_b' on the common frozen task axis. Queries of A_a' K_j A_b then measure
reproduced task-linked effects. For a fixed task covariance Phi, a voxel query
can read A_a[v,] Phi A_b[v,]'. Its signed evidence is not the fitted nonnegative
signal-SD map.

Use `lm_relation_fit` only when its separable observation/error assumptions are
actually met. rMVPA's CR1 and restricted wild-bootstrap covariance objects do not
automatically qualify as crossform's separable-GLM residual model. Declare the
capability translation explicitly and retain narrower rMVPA inference otherwise.

Exit: correct feature/target/partition mappings; independent error oracle where
admitted; refusal for pooled-only input, incompatible covariance, missing
preprocessing ancestry and reused confirmation origins. No change to the existing
loading-level group estimand.

## 4. Metric and locality extension: do not force it into I1

Initially use identity or an admitted fixed, bounded metric for geometry export.
Crossform's matrix-based `neural_metric()` is not a license to densify a
whole-brain structured precision. A later operator extension must distinguish
`restrict_metric` from `marginal_precision`: they compute different operators.

An ROI-only predictor uses (Psi_RR)^-1. Re-estimating a local covariance creates
an adapted predictor. Frame restriction of a fixed Q uses its declared native or
whitened composition and may not mean either predictor. Each operation records
support, units, training data and applicability to additive conservation.

Gate any extension on dense small-matrix parity, singular/rank-zero behavior,
outside-ROI poisoning tests, non-diagonal counterexamples and measured memory.
Current small-node `measurement_frame` and population metric gates remain in
force; a numerical wrapper cannot promote an unsupported operator.

## 5. Common state and receipts

The adapter negotiates `scientific-state-v1` and the geometry-prediction protocol
independently. Data meaning, selection state and implementation compatibility are
separate checks. Unsupported versions return a structured refusal before reads.

Publish one manifest of adapter conformance cases, not a second scientific claim
registry. Map evidence to the existing crossform registry only when a named claim
is being advanced. rMVPA benchmark receipts remain owned by rMVPA. This design
revision does not change their evidence class or recertify their source.

Adapters preserve observation weighting versus scoring weighting, fold-resolved
versus pooled prediction ledgers, selected rank versus rank support, and filtered
versus estimated-zero features. Missing uncertainty stays unavailable.

## 6. Acceptance matrix

| Gate | Required independent check |
|---|---|
| Factorization | Direct B K B' equals T(A'KA)T' and fixed-query contraction |
| Rotation | Invertible coordinate transport preserves whole-form predictions; singular transforms refuse |
| Cost | Same scientific plan under two tile sizes; bounded reads; no p-by-p allocation |
| Locality | (Psi_RR)^-1 differs from cropped precision on a planted correlated case; local predictions match the former |
| Evidence | Negative null geometry kept; nonzero null-prediction gain penalized; PSD descriptions never relabeled unbiased |
| Components | Coherent/configuration evidence adds before one cost; arbitrary oblique gains retain interactions |
| Identity | Permuted named axes align; changed units/basis or unseen target names refuse |
| Independence | Aliases and upstream normalization dependencies expose overlaps; CV train means are not independent partitions |
| Group | Common experimental space is necessary; loading pooling and geometry pooling remain distinct |
| Recovery | Restart/interruption preserves finished artifacts, RNG and outcome exposure; plan changes fork |

Every gate has a positive fixture, a negative fixture and a mutation-sensitive
oracle where numerical equivalence is claimed. Wire only the relevant subset
into each PR. A missing R dependency is a recorded unrun gate, not a pass.

## 7. Release and scope

Order: metadata-only I0; frozen prediction I1; confirmation I2; then measured
operator/retention extensions. Each PR preserves the existing API and has a
source-pinned receipt. Crossform code changes follow its existing certification
binding procedure; documentation-only changes do not earn new numerical status.

TV envelopes, graph-local precision, sequential supported-rank tests, arbitrary
oblique inferential attribution and joint hierarchical fitting remain separate
research/implementation decisions. This integration is not a reason to reopen
the accepted `pattern_model` name or move its numerical core.
