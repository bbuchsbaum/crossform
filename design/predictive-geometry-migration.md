# Predictive geometry migration

Decision: G01, 2026-09-04. Contract: `predictive-geometry-v1`.

The supported predictive workflow is `model_basis()` -> `fit_geometry()` ->
`score_geometry()`. The model declaration is independent of neural outcomes;
the fit is a biased latent prediction; the score is independent signed
evidence. The public surface follows these roles, not an export-count target.

| Existing behavior | Decision | Destination and test obligations |
|---|---|---|
| `model_basis()` centered factor construction | Generalize | Retain the explicit centered frame, rank clamp, axes and lowering products. Add kernel inputs and declared normalization; T01–T11. |
| Constructor overlap refusal | Replace | Admit redundant factors and report union/overlap rank. A fitted form does not require identified redundant coefficients; T08/T09. |
| Permissive default model negative mass (1%) | Replace | Default `negative_share=0` tolerates numerical roots only. A nonzero declared limit explicitly admits positive-part model replacement and records its cost; T03. |
| `model_geometry()` trace decomposition | Retain, advanced diagnostic | Its exact signed addressable/orthogonal trace decomposition is useful and distinct from Frobenius predictive gain; T37/T39/T59. |
| `model_geometry()` isotropic/diagonal/block/shared coefficients | Retain, advanced hypothesis-specific diagnostic | These represent different constraints. Require identified model coordinates in this reader until it has a separately specified redundant-coefficient convention. It is absent from the primary predictive workflow. |
| `model_geometry(training=...)` | Retain its projector-energy meaning | Do not rename it predictive gain. New cross-fitting calls the fitted-form scorer and retains amplitudes and norm offsets; T33/T45/T49. |
| Shared full-span saturation refusal | Keep for the existing descriptive model-test claim | New prediction permits full-span fits, labels zero-penalty results model-agnostic, and compares regularized kernels with same-span isotropic baselines; T21/T54. |
| `relation(extract=basis)`, `plan_relation(effects=basis)` | Retain | Typed lowering and linear-stage restrictions remain; add actual training ancestry where prediction needs it; T24–T32. |
| `latent_geometry(rank=...)` | Retain descriptive semantics | Reuse its numerical PSD/rank primitive after independent tests. Predictive fit/score records never inherit raw geometry or descriptive projection receipts; T12–T23/T42. |
| `rsa()` | Retain unchanged | Fixed unrestricted linear regression in RDM observation space remains a different scientific objective from nonnegative additive fitting and from predictive Frobenius loss. |
| Old model basis serialized records | Version explicitly | New basis signatures bind normalization, PSD policy, centering and complete model records. Old-schema values must be rebuilt from model inputs, not silently reinterpreted. |
| Old recertification artifacts | Preserve as baseline until recertification | They certify the original source digest only. G22 replaces them with evidence bound to the implemented source, and separately adds predictive evidence. |

Callers of the retained advanced reader continue to receive the same estimands
on admitted nonoverlapping input. Its overlap guard moves to the reader rather
than preventing compilation of useful overlapping models. No existing
statistical result is silently renamed, and no claim is retained solely
because an earlier test asserted it.

The initial tree is recoverable from the archive recorded in
`.planning/2026-09-04-predictive-geometry-implementation.json`. No commit or
publication is part of this implementation request.

G19 keeps 123 governed exports: the prior 121 plus `fit_geometry()` and
`score_geometry()`. This count follows the role decisions above. The primary
reference section groups the basis compiler with those two verbs; the
descriptive reader has a separate section. Cross-fit scheduling, nested
selection and training-local caches are currently private orchestration
contracts, not additional public fitting frameworks. The public fixed-split
workflow is executable in `vignette("predictive-geometry")` and the standalone
exemplar. Empty and repeated score-row selections preserve their documented
summary/mode/group schemas and ordering.
