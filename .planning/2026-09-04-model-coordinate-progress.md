# Model-coordinate epic — implementation progress (2026-09-04)

Epic `bd-01M1P2X7BB4TP1D2JXA2HAK3EA`; index `.planning/2026-09-04-mote-index-model-coordinate.md`;
contract `design/model-coordinate-geometry-contract.md`; oracle `design/oracles/model-coordinate-geometry.R`
(now O1–O8); package witness `design/oracles/model-coordinate-package-route.R`.
Branch `elite-pass`. **Nothing is committed**; everything below is in the working tree.
Certification artifacts are STALE by design until M9 (test-certification-artifacts.R skips, does not fail).

## Review status (2026-09-04, later)

- M2: REJECT → all fixes applied (deflated construction in a basis V of 1-perp; validator re-derives fields;
  name alignment under `extract`; etc.); contract §1.1(d)/§1.3/§1.4/§15.1a amended; tests green (246).
- M4: ACCEPT WITH AMENDMENTS → all applied (rank/symmetry validation on `.latent_rank_psd_form`, binding-row
  check, budgeted formats, docs); tests green. Nit for M9: regenerate `exemplars/haxby2001/results/conservative-geometry.rds`.
- M3: reviewer failed (session limit) before reporting; needs a new review together with M5.

## Later the same day: M6, M7, M8 landed in the working tree

- M6: `model_geometry(training = metric_training_policy("exclude_evaluation"))` → `$cross_fit` (estimation layer),
  `.model_geometry_cross_fit()`; refusal `cross_fitted_model_energy`/`insufficient_disjoint_edges`; tests incl.
  exactness vs the original plan's executor and a 30-rep Monte Carlo null. mote M6 begun (`rv-01M1Q59X9KD30MZS5AT1CVTY5F`).
- M7: `tests/testthat/test-model-coordinate-contract.R` (README-fixture anchors: 1.15 trap, R1 1e-12, R2 bank,
  peak 144 split 4.10545 = 4.10297 + 0.00248, latent clipped 79/280; vocabulary and refusal names in the
  contract and terminology; layers apart).
- M8: `design/terminology.md` (seven rows + geometry/basis sentence), `design/relation-to-prior-work.md`
  (new section + Cook 2018, Jozwik 2016, Kaniuth & Hebart 2022 refs), `vignettes/novelty.Rmd` (ledger row),
  `README.md` (advanced-tier paragraph), `_pkgdown.yml` (section with both exports).
- Full suite (NOT_CRAN=true, `testthat::test_local`) after M2–M8: **10187 passed, 0 failed, 0 errors, 16 skipped**
  (certification STALE/UNBOUND by design until M9, opt-in scale gates, CRAN-only snapshots).
- M3+M5 review: both ACCEPT WITH AMENDMENTS; all applied and re-tested (lowering 82, model-geometry 209,
  contract 53 expectations pass). M3 and M5 closed in mote. Pre-existing measurement-path defect filed as
  `bd-01M1Q62WMQVY1S210CH8ZH9J2H`.
- M6+M7 review: M6 ACCEPT WITH AMENDMENTS (all applied: pairing-driven edges and floor, refusals for star/subset/
  self-pair pairings, invariant checks, bookkeeping validation, wording, contract §7.2/§14), M7 ACCEPT. M6, M7, M8
  closed in mote. **R/ frozen for M9 from here.** (agent `m3-m5-reviewer`), an M6 review, `mote done` for M2–M8, M9
  (re-run `benchmarks/run-crossnobis-scale-gate.R`, `run-shard-admission.R`, `run-sampling-covariance-validation.R`,
  `run-learned-metric-policy-validation.R`, then `benchmarks/promote-artifacts.R`), full test suite, commit. The Haxby exemplar rds
  (`exemplars/haxby2001/results/conservative-geometry.rds`, old-schema latent receipt) needs the data cache,
  which is absent here (`bash exemplars/data.sh fetch` first, then `Rscript exemplars/haxby2001/07-conservative-geometry.R`);
  recorded on M9 as a follow-up, not blocking (tests read the CSV receipts only).

## M9 done (2026-09-04, evening)

Re-certified on frozen digest `sha256:24afa8ce392d…` (before == after): nine runners passed, promotion clean,
certification court 540/0/4 designed skips, empty-leaf walk clean; recorded in `design/certification-report.md`.
M2–M9 all closed in mote; M10–M12 remain deferred under the open epic. Final full suite on the certified tree:
**10433 passed, 0 failed, 0 errors, 7 skipped** (all opt-in or designed; no STALE). **Nothing is committed** — the tree (R/, tests, man, design, README,
vignette, inst/extdata/certification, benchmark-results) is ready for the maintainer to commit.

## Done in the working tree

| ticket | state | files | tests |
|---|---|---|---|
| M2 `model_basis()` | code + tests green; fresh-context review **running** (agent `m2-reviewer`) | `R/model-basis.R` (new), `man/model_basis.Rd`, `NAMESPACE`, `tests/testthat/test-architecture.R` (layer_of), `tests/testthat/test-api-surface.R` (list + count 120), `design/api-tiers.md` (row + totals), `_pkgdown.yml` (new section) | `tests/testthat/test-model-basis.R` (new, 134 expectations pass) |
| M4 `latent_geometry(rank=)` | code + tests green (snapshots unchanged); review **running** (`m4-reviewer`) | `R/latent.R` (rank budget, two-part mass, `.latent_rank_psd_form()`, schema_version 2 identities), `man/latent_geometry.Rd` | appended to `tests/testthat/test-latent-geometry.R` (4 new blocks; one pre-existing forgery fixture extended) |
| M3 lowering entries + gate | code + tests green; review **running** (`m3-reviewer`) | `R/relation.R` (`extract = basis`), `R/relation-plan.R` (`effects = basis`), `R/evidence-task.R` (`.model_coordinate_lowering_gate()` in `.new_evidence_task()`), `R/model-basis.R` (`.is_model_coordinate_space()`) | `tests/testthat/test-model-coordinate-lowering.R` (new, 72 expectations pass) |
| M5 `model_geometry()` | code + tests green; registered (layer_of 5L, api-surface 121, api-tiers row, pkgdown); Rd generated; contract §4.2/§13/G3 amended with O8; review **pending** (combined M3+M5 reviewer to be spawned; the first M3 reviewer died at its session limit) | `R/model-geometry.R` (new), `R/relation.R` (`.relation_lowered()`), `man/model_geometry.Rd` | `tests/testthat/test-model-geometry.R` (new) |

Oracle O8 (block coordinate descent: monotone to 4.5e-13, fixed point 3.6e-10, nesting shared ≤ block ≤ diagonal)
was added and runs clean; contract §4.2 / §13 still say "unmeasured" and must be amended when M5 lands
(also update the oracle header comment to list O8).

## Decisions taken (record in mote if not already)

- M2: centering failure is `effect_invariant_error` (package-bug class), tolerance scales with `basis_condition`
  (backward-stability bound); products `$extractor` always, `$effect_map` only with a `condition_space()`;
  Q coordinates named `mc1..mcq`, R columns `<model>.<k>`; no `provenance` arg; `negative_share` default 0.01;
  `distance` required with RDMs (closed set `"squared_euclidean"`), features via `features =`; saturation
  refusal helper `.model_basis_refuse_saturated()` for M5; near-saturation warning fraction lives on the fit (M5).
- M4: `rank = NA_integer_` when unbudgeted; ALL latent identities changed (schema 1→2); `moved_mass` stays the
  sum; as.data.frame/print add the two parts only under a budget (snapshots byte-identical).
- M3: gate in `.new_evidence_task()` (only place both relations meet the stage plan); `covariance()` admitted
  per contract wording; conveniences `relation(extract = basis)` and `plan_relation(effects = basis)`.
- M5 design: `model_geometry(x = ORIGINAL plan, basis, structure, rank, component, warning_fraction, ...)`
  lowers via composed extractors (`.relation_lowered()`), re-plans with the same frame/pairing/metric/compute,
  materializes the lowered form, reads `tr(HGH)` and `tr((H-P)G)` as a two-column packed bank on the original
  plan, fits per measurement (isotropic/diagonal NNLS in-house Lawson–Hanson; block = exact block descent using
  `.latent_rank_psd_form`; shared = `.latent_rank_psd_form(S, rank)`), returns `C` (array), `fit` df, `weights`,
  `df$fit/free`, latent reading line on fit fields only.

## Next steps, in order

1. Collect the three review reports (agents `m2-reviewer`, `m4-reviewer`, `m3-reviewer`; their final reports
   arrive as task notifications). Apply amendments, rerun the affected test files, then `mote done` M2, M4, M3
   with a note each (M2 `rv-01M1P4EXWNH22FZXGBJ1SH1DC3`, M4 `rv-01M1P56JXP141EBR6ZJJA010BH`,
   M3 `rv-01M1P5Y2HFB18S7T5ZPKANAMHP` are the reservation ids from `mote begin`).
2. Finish M5: add `.relation_lowered(x, basis)` to `R/relation.R` (copy relation, replace `$extractors` with
   `effect_extractor(t(Q) %*% e$map, effects = basis$effect_space, estimator = "model_basis", diagnostics = ...)`,
   `$effect_space`, `$effects`; `.validate_relation()`); write `tests/testthat/test-model-geometry.R` (trace split
   identity on README fixture vs oracle-style algebra; NNLS KKT; diagonal g from S equals g from G; block descent
   monotone + nesting; shared == `latent_geometry(lowered, rank)` energy; saturation refusal; learned-metric refusal;
   df values; print/as.data.frame; forgery); register file (layer 5) and export; api-tiers row; pkgdown; document();
   amend contract §4.2/§13 with O8 numbers; fresh-context review; `mote done` M5.
3. M6 cross-fit (`model_geometry(..., training = )`, four-partition floor, refusal `insufficient_disjoint_edges`),
   M7 contract test file, M8 vocabulary/novelty/prior-work/README, M9 re-certify (`benchmarks/run-*.R` then
   `benchmarks/promote-artifacts.R`), then commit.
4. Housekeeping each `devtools::document()`: `git checkout -- man/population_wild_bootstrap.Rd` (roxygen
   trailing-space churn from a newer roxygen2 than RoxygenNote 7.3.3); `man/fmrireg_relation.Rd` change is a
   legitimate sync with the last commit and should be kept.

## How to run things

```sh
Rscript -e 'suppressPackageStartupMessages(pkgload::load_all(".", quiet=TRUE)); testthat::test_file("tests/testthat/test-model-basis.R", reporter="summary")'
NOT_CRAN=true Rscript -e '...test_file("tests/testthat/test-latent-geometry.R"...)'   # snapshots need NOT_CRAN
Rscript design/oracles/model-coordinate-geometry.R
Rscript design/oracles/model-coordinate-package-route.R
Rscript -e 'suppressMessages(devtools::document())'
```
