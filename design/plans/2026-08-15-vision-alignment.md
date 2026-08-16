# Vision-alignment plan — 2026-08-15

Goal: make the codebase correctly, coherently, and beautifully implement the
closure-calculus vision recorded in `novelty.md` (2×2 boundary closures ×
spatial frames, coherent/configuration decomposition, query-first execution,
estimand-bound generalization, executable capability contract, evidence
ledger with gates).

Basis: four independent full audits (public API surface, execution layer,
capability/refusal contract, documentation/gate artifacts), each verified by
reading source and executing the package. File:line references below come
from those audits.

## Verdict

The self-form spine is real and well-built: `plan_geometry()` →
`contrast_energy()`/`rdm()`/`rsa()` are provably views of one object; `contrast_energy()`
is genuinely query-first (one fused pass, 27× cheaper than materialization
at one edge); plan/receipt separation exists with the right shape of test;
the precomputed-beta uncertainty refusal is exemplary; `novelty.md` and
`README.md` already carry the converged framing with an honest ledger.

But four of the vision's five headline claims are currently contradicted —
not merely unfinished — in specific, measured ways:

1. **Query-first**: `rdm()`/`rsa()` eagerly build an O(q⁴) dense query
   matrix (191 MB at q=100 before any data is read) and the query path is
   3.5–3.9× *slower* than full materialization for the flagship RDM query.
   There is no public selected-pairs API; `pair_query()` is refused by the
   compiler. Memory admission omits the query matrix (unsound by up to 20×).
2. **Generalization in estimand identity**: cross-run vs cross-session plans
   get identical `scientific_plan_id`s when partition labels are generic
   (`p1..pK`) — Γ identity currently rides on free-text strings. Meanwhile
   the *same* estimand gets *different* ids depending on route
   (`rdm(plan)` vs `rdm(materialize_geometry(plan))`). Identity is both too coarse and
   too fine.
3. **Executable refusals**: all capability refusals are unclassed
   `stop(<string>)` (1059 stops; zero classed capability conditions).
   Failure-gallery cases 1 and 3 (correlation-distance, univariate removal)
   do not exist as refusals — users get `unused argument (...)`. The
   learned-metric admission is a class check; `metric_status` is hard-coded
   `"fixed"` for every geometry plan (`evidence-sampling.R:113`); the actual
   block fires by accident with a *wrong* message.
4. **Closure calculus**: rectangular forms are internal-only (kernel fully
   supports them; no exported constructor), coupling is a disconnected
   pipeline with its own vocabulary and a 256 MiB ceiling (so coupling ×
   searchlight is impossible), and conservation Σ_x G_x = G_Ω holds only for
   `total` (coherent residual 2.49 on a 3×3 grid) while `searchlights()`
   defaults to the non-conservative `"local"` normalization.

Docs are the *most* aligned layer; code must catch up to its own claims.

---

## Phase 1 — Integrity (nothing may contradict the contract)

Cheap relative to stakes; do first. Every item is a place where current
behavior falsifies a stated guarantee.

1. **Fix learned-metric admission (scientific-integrity bug).**
   - Derive `metric_status` from `metric$estimation`/`capabilities$learned_frozen`
     instead of hard-coding `"fixed"` (`R/evidence-sampling.R:113`).
   - Replace the `inherits(x, "effect_crossnobis_plan")` type test with a
     capability check (`R/evidence-sampling-product.R:371-378`).
   - Preserve `frozen`/`training_signature` through support restriction
     (`R/evidence-sampling-product.R:86-96`) so the misleading
     "requires frozen training provenance" message can't fire on a metric
     that has exactly that provenance.
   - Read `calibration_requires_metric_uncertainty` somewhere, or delete it.
   - Bind `neural_metric(estimation="learned_frozen")` provenance to a real
     training object, not a 64-hex regex (`R/metric.R:155-162`, `R/source.R:3-6`).
   - Tests for each refusal.

2. **Typed generalization axis.** `pairing(..., generalizes_over =
   "run"|"session"|"task"|...)` (and via `cross_partitions()`), bound into
   `.effect_task_semantic` and `.metric_pairing_identity` so cross-run vs
   cross-session yield distinct ids regardless of label strings. Test with
   identical generic labels. Also reconsider `independence = "independent"`
   as a silent default (`R/pairing.R:17`) — the declaration that licenses
   crossnobis unbiasedness should be explicit on the paths that consume it.

3. **Route-stable `scientific_plan_id`.** Same estimand + query must hash
   identically whether executed fused (`rdm(plan)`) or projected
   (`rdm(materialize_geometry(plan))`). Derive the id from the scientific request
   (task + frame + metric + pairing + query), never from the execution
   route's recompiled component. Add the id-equality assertions the tests
   currently omit (`test-geometry-plan.R:150-152`).

4. **Classed refusals.** One `effect_capability_refusal` condition class
   (pattern already exists at `R/tomography.R:236-243`) carrying
   `capability`, `namespace`, `reasons` (all of them, not
   `unavailable_reasons[[1L]]`), and `remedies`. Route every `.require_*`
   through it. Give the eight cryptic sampling reasons real sentences +
   remedies and tests (currently zero test hits for six of eleven reasons).

5. **Make failure-gallery cases 1 and 3 exist.** `rdm(normalize=)` and
   `contrast_energy(remove_univariate=)` accept-and-refuse stubs that name the
   capability (`guaranteed_psd` — currently computed but gating nothing)
   and point to `correlation-distance-policy.md` and
   `component = "coherent"|"configuration"|"total"` respectively.

6. **Honest invariance + identity tests.**
   - Replace the exact-integer tiling-invariance fixture with a
     floating-point one and a declared tolerance; cover the query-fused path
     (`test-effect-form-certification.R:109-131`; measured drift ~5e-15
     under `block_features` changes).
   - Bind the q=120 sampling-covariance scale artifact to a q=120 plan:
     `.sampling_covariance_from_components` must check contrast width
     against the plan's effect space (`benchmarks/run-sampling-covariance-scale.R:49-68`,
     `R/evidence-sampling-kernel.R:79-98`).
   - Test exact recomposition (total = coherent + configuration) at the
     public level for `rdm()` and `rsa()` (currently only `contrast_energy()`;
     `component="coherent"` never passed to rdm/rsa in any test).

7. **Conservation honesty.**
   - Decide the `searchlights()` default: either `"conservative"` by
     default, or document loudly that the default breaks Σ_x G_x = G_Ω.
   - Document + test that conservation holds for `total` only (coherent is
     locally defined and does not sum; measured residual 2.49 on a toy grid).
   - Export a `frame_conservation()` diagnostic.
   - Warn in `geometry_spectrum()` docs that `component=` is accepted but
     eigenvalues are not additive.

8. **Kill remaining silent fallbacks.** `pmin(1, pmax(0, svd(...)$d))`
   without a bound check and the silent pseudo-inverse in
   `R/measurement-decomposition.R:263-308` (its sibling in
   `coupling-views.R:459-462` errors — make them consistent);
   `reconstruct_evidence(allow_projection = TRUE)` default; message on
   invalid `coherence_fraction`. Also fix `plan_crossnobis` check order so
   a missing residual channel is diagnosed as such, not as a
   metric-training-partition problem (`R/crossnobis.R:227-253`), and give
   unequal-length-runs a legible capability error in `lm_relation_fit`.

## Phase 2 — Query-first for real

The claim must become true before Gate 5 can exist.

1. **Factored edge queries.** RDM edge queries are rank-1
   (`tcrossprod(e_i - e_j)`); represent queries as sparse/factored operators
   and push the factorization into `.effect_form_feature_task`
   (`R/task.R:150-160`, `R/kernel.R:976-996`) instead of the dense q×q
   matmul per view. This removes both the O(q⁴) eager matrix
   (`R/views.R:137`) and the measured 3.5–3.9× regression.
2. **Public selected-pairs API.** `rdm(plan, pairs = ...)` and a compiler
   path that accepts `pair_query()` (`.compiler_query`,
   `R/compiler.R:229-246`). "The RDM is optional" must be callable.
3. **Sound memory admission.** Add the live query matrix (and the
   validation-memo retention, `R/validation-memo.R:27-41`) to
   `.compiler_memory_plan`'s categories; add one gate comparing
   `planned_workspace_bytes` to a measured peak.
4. **Cheap wins.** Single-pass `contrast_energy()` on materialized geometry
   (currently two full store reads, `R/views.R:51-52`); guard the double
   `source_session$close()`; decide the fate of tested-but-unwired
   checkpoint/scheduler/tiled-contraction code (wire or move to attic —
   test coverage currently overstates delivered capability).
   *Resolved 2026-08-16: `R/scheduler.R` and `R/checkpoint.R` deleted with
   their tests; `.tiled_contraction()` retained (it has compiler callers).*
5. **Gate 5 benchmark.** `benchmarks/run-query-first-scale.R`: large-q
   trial-level, selected edges + RSA coefficients, materialization receipts
   (`pair_atoms_materialized: FALSE`), recorded result artifact + opt-in test;
   commit the artifact together with the implementation before claim promotion.

## Phase 3 — One calculus, one vocabulary (beauty; breaking changes OK pre-release)

1. **Public rectangular plans.** `plan_geometry(x, right = ..., at, over)`
   (kernel already supports `codec="rectangular"` with the decomposition,
   tested only via `:::`). Ordered/pair queries on the plan path; pair-set
   `rdm()`/cross-axis `rsa()`; covariates via `pair_lm_query` on a spatial
   rectangular plan. This is closure (c)'s public realization and the
   prerequisite for Gate 3.
2. **Coupling reachable from the plan.** A `coupling(plan, between, by)`
   verb (or shared plan object) so the adjoint closure reuses the same
   frame/pairing/query vocabulary instead of a parallel pipeline; a staged
   plan for lifting the 256 MiB small-node ceiling (searchlight-resolved
   coupling is currently impossible).
3. **Naming unification.** One public noun for the evidence pairing
   (`evidence_*`); resolve "frame" (3 meanings), `geometry_alignment` vs
   `geometry_*` collision, "form" vs "geometry" for the same object.
   Rename or absorb: `connectivity(view=)` duplicating three coupling
   exports; `crossnobis()` gaining the decomposition it currently drops
   (`R/crossnobis.R:524`) or becoming an explicit thin alias of
   `contrast_energy()$total` with a metric gate.
4. **Prune dead exports.** `inner_product`, `measurement_bridge`,
   `reverse_bridge`, `measurement_space`, `pairing_marginals`,
   `additive_frame`, metric-plumbing quintet — unexport, wire up, or
   document as advanced with a real consumer. Split the seven
   one-bit relation-fit capability flags into independently earned bits or
   collapse them honestly.
5. **Capability introspection.** Export `sampling_capabilities(plan, fit)`
   returning the full reason table; extend `print.effect_geometry_plan` to
   show metric status, pairing/generalization axis, and covariance
   availability. A scientist should be able to ask, not provoke.

## Phase 4 — Evidence and documentation (describe the fixed reality)

1. **Failure gallery** (`failure-gallery.md` or vignette): the five agreed
   cases, executable, with verbatim classed refusals — three already have
   verbatim evidence in `exemplars/haxby2001/results/smoke-report.md`,
   two become possible after Phase 1.4–1.5. Link from novelty.md + README.
2. **Evidence ledger upgrade** in `novelty.md`: the four labels
   (established algebraically / implemented / demonstrated / prospective),
   each row linked to its machine-checkable artifact (test file, law
   helper, exemplar script, benchmark RDS); labels inside the 2×2 table
   cells (mirror in README); restore the Diedrichsen, Provost &
   Zareamoghaddam (2016) citation.
3. **Gate 2 completion**: show `effect$signed` beside the three energies in
   `introduction.Rmd` + README; add the energy-vs-signed caution and the
   w-weighted orthogonality qualifier; optionally a Haxby `06` decomposition
   script as the real-data flagship.
4. **Gate 1**: rsatoolbox parity via pinned-version fixture export
   (no Python in CI), `exemplars/haxby2001/06-rsatoolbox-parity`; also
   resolve the recorded rMVPA shadow-install hazard.
5. **Gate 4**: overlapping-searchlight double-counting vs conservative
   frame demo with an interpretive punchline (uses Phase 1.7).
6. **Gate 3**: one real ER/cross-task rectangular exemplar (needs Phase 3.1
   and a suitable dataset — longest pole; can be scouted in parallel).
7. **Small fixes**: retitle the `evidence-pairing` cross-references in
   README:171-173 and introduction.Rmd:270-272 (the vignette is a
   coupling/measurement vignette, not the promised rectangular one) or
   split it; verify `novelty.html`/`correlation-distance-policy.html`
   render in the pkgdown build; align `mission.md`/`vision.md` with the
   calculus framing.

## Sequencing rationale

Integrity first because those items are cheap and every one is a live
contradiction of the package's own contract philosophy (the learned-metric
admission hole is the kind of silent substitution the architecture exists
to prevent). Performance second because "query-first" is currently a
measured regression and Gate 5 cannot be honest until it isn't. API shape
third because renames/unexports are free now and expensive after adoption.
Docs last because they must describe the repaired reality — they are
already the most aligned layer and would otherwise drift twice.

Verification bar per phase: full test suite + the relevant law helpers;
R CMD check clean; for Phase 2, recorded benchmark artifacts with receipts,
committed before claim promotion; for Phase 4, pkgdown build + exemplar reruns.

## Execution record (2026-08-15)

This record closes the four-phase Vision-alignment implementation epic. It
does not close the broader novelty epic: Gates 1, 3, and 4 remain open as
listed below.

- **Phase 1 — complete.** All eight items; R CMD check Status: OK.
- **Phase 2 — complete.** Structured rank-1 queries end to end; Gate 5
  landed (q = 100 / 1,080 searchlights: selected-100 0.23 s, fused 4.6 s vs
  materialized 18.3 s, oracle 4.4e-16, < 512 MiB); `rdm(pairs = )` public.
- **Phase 3 — complete except deferred renames.** Public rectangular plans
  (`plan_geometry(right = )` + `pair_query()` execution + rectangular
  materialization), `coupling()` adjoint closure from the plan vocabulary,
  five dead exports internalized (bridge quartet retained per the frozen
  vocabulary certification), `sampling_capabilities()` + enriched plan
  print, crossnobis documented and tested as the named total of the
  metric-carrying contrast. Pure renames recorded as a 0.2 surface pass on
  the Phase 3 bead.
- **Phase 4 — documentation current.** Four-label artifact-linked ledger,
  closure-table cell labels, failure-gallery.md (Gate 6 landed), Gate 2
  landed via the one-plan vignette family, Diedrichsen et al. (2016)
  citation restored, README/vignette cross-references corrected, vision.md
  aligned. Open gates tracked on the novelty epic: rsatoolbox parity
  (Gate 1), real rectangular exemplar (Gate 3), operational conservation
  demo (Gate 4).
