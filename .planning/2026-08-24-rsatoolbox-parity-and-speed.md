# Reproducing rsatoolbox, and the scaling comparison

**Date:** 2026-08-24 · **Branch:** `elite-pass` (clean at `0217add`)

## The question

A researcher with a working rsatoolbox analysis asks: *can I run this in
crossform, get the same numbers, and get them faster?*

## The finding

**Crossform already reproduces the standard RSA workflow end to end. Nobody has
ever shown it.** Every step below was verified by execution on 2026-08-24, not
inferred from the source:

| # | Step | rsatoolbox | crossform route | Verified |
|---|---|---|---|---|
| 1 | Patterns from data | `Dataset` | `relation()`, `lm_relation_fit()` | betas = per-run OLS at **6.7e-16** (existing) |
| 2 | Noise precision | `prec_from_residuals(method="full")` | `noise_precision()` | **2.3e-14** (existing) |
| 3a | `calc_rdm("crossnobis")` | | `rdm()` + `cross_partitions()` | **3.8e-15** (existing) |
| 3b | `calc_rdm("euclidean")` | | `rdm()` + biased self-pairing | **2.2e-16** (new) |
| 3c | `calc_rdm("correlation")` | | `G` recovered downstream from a PSD self form | **4.9e-15** (new) |
| 4 | Model RDMs | `ModelFixed` | `models=` argument | equivalent |
| 5 | `compare(method=)` | | base R on `rdm()$values` | **7.2e-16** across cosine, corr, spearman, kendall, tau-a (new) |
| 6 | Inference | `eval_bootstrap`, noise ceilings | — | out of scope by standing declaration |

**So the work is evidence, not features.** This plan adds **zero exports** and
amends **zero contracts**. That is not a constraint I worked around; it is what
the founding documents require, and the verification shows it was never
actually in conflict with the goal.

Two points of care, both load-bearing:

- **Step 5 is downstream by design, not by evasion.** `rdm()$values` is a
  measurement × pair matrix already in `np.triu_indices` order. The similarity
  statistic is then `cor(d, model, method = ...)` — a function of two numeric
  vectors that knows nothing about brains. `design/api-tiers.md:305` records
  that **no exported function accepts a normalizer** (`cosine()`/`correlation()`
  are private "for exactly this reason", and `inner_product()` was demoted to
  enforce it), and `design/common-geometry-equivalence.md:156` excludes
  "Spearman or Kendall RSA, rank transforms, correlations or cosine
  normalization of the observed RDM" from the equivalence theorem. Exporting a
  similarity function would contradict both. Demonstrating the statistic
  downstream contradicts neither, and satisfies *"scalarization occurs as late
  as practical"*.
- **Step 3c does not reopen the correlation-distance refusal.**
  `rdm(normalize=)` refuses because *cross-generalized* diagonals can be zero or
  negative. The downstream route uses a **within-sample** self form, which is
  PSD with strictly positive diagonals — the case
  `vignettes/correlation-distance-policy.Rmd` already classifies as a legitimate
  "disciplined nonlinear view". `G_ij = (G_ii + G_jj − d_ij)/2` from `rdm()` and
  `contrast_energy()`, then `1 − G_ij/√(G_ii G_jj)`. **Reproduction is not
  endorsement**: the package's recommendation to prefer crossnobis for
  cross-generalized inference is unchanged.

---

## Track A — The demonstration  ✅ COMPLETE 2026-08-24

**Claim:** every rsatoolbox quantity in the standard workflow is shown to agree,
under a version-pinned environment, with the conventions recorded rather than
conceded.

The exemplar already has the right architecture — R → pinned Python → R, CSV as
the only channel, no reticulate. This is new rows in an existing table.

1. **Distances.** Euclidean and Mahalanobis via
   `pairing(r, r, self_pairs = "allow_biased", independence = "not_independent")`;
   crossnobis already landed. Correlation distance via the PSD self-form route,
   which must **check** the positivity precondition and refuse otherwise, exactly
   as the policy's required contract specifies.
2. **Noise estimators.** `prec_from_residuals(method="diag")` against
   `diagonal_precision()`. Shrinkage as **conditional parity**: rsatoolbox
   derives its coefficient analytically, crossform takes λ as a declaration and
   refuses to tune it on evaluation data — feed crossform rsatoolbox's λ and show
   the precisions coincide. That is a policy difference, not a missing feature,
   and should be labelled as one.
3. **Comparison statistics.** `06-similarity.R` computing cosine, Pearson,
   Spearman, Kendall tau-b, tau-a, and rho-a from `rdm()$values`, against
   `compare()` on the same vectors. **The fixture must contain deliberate ties** —
   verified: without ties, `tau-a` equals `kendall` and `rho-a` equals
   `spearman`, so four of the six rows would silently test nothing.
   Whitened `corr_cov` / `cosine_cov` are the same pattern with Σ from
   `sampling_covariance()` — and are where crossform is *ahead*, since
   rsatoolbox takes `sigma_k` as an input while crossform derives it.
4. **The coverage ledger.** `coverage.csv`: one disposition per rsatoolbox
   capability — `parity`, `parity_conditional`, `reproducible_downstream`,
   `refused_by_design`, `out_of_scope` — each with a pointer to its evidence row
   or its contract. Drift-tested so it cannot rot away from reality.

**And make the table worth trusting.** Four holes currently undercut anything
built on it: the ratchet never fires under `R CMD check` (`^exemplars$` is
Rbuildignored, so `skip_if` swallows it) and no CI workflow runs Python;
`check-certification-binding.R:32` globs `\.rds$`, leaving the parity receipt as
the one certification artifact with no source-digest binding; the manifest
hashes no crossform source and no git SHA; and `requirements.txt` leaves `scipy`
unpinned (it drifted 1.18.0 → 1.18.1 during verification). Fix all four, add a
CI job running the pipeline from a source checkout, and correct the README
sentence claiming a CI enforcement that does not exist.

**Gate A — met.** Agreement table 9 → **24 rows**, worst **5.11e-15**; 23 rows
at `1e-10` and one at `1e-8` (the documented `fit_regress` route). Three noise
rows were *tightened* from `1e-8` to `1e-10`; none loosened. **Zero new
exports** — `git diff` over `NAMESPACE`, `R/`, `src/` is empty. Full suite
**9683 pass / 0 fail / 0 warn**. `coverage.csv` gives 28 capabilities one
disposition each, drift-guarded. Evidence chain closed: environment fully
pinned (25 packages), manifest schema 2 binds the crossform source digest and
git commit through `benchmarks/provenance.R` across 31 artifacts,
`check-certification-binding.R` now covers the CSV receipt and was
negative-tested in both failure modes, `run-receipt.csv` carries execution
provenance, and the README's false CI claim is corrected.

*One item unverified:* `.github/workflows/rsatoolbox-parity.yaml` is written
and YAML-valid but has never executed — it needs a push to confirm.

Findings worth carrying forward:

- **Correlation distance is not frame-composable.** Its centering *and* its
  normalizer both depend on the support, so each support needs its own plan,
  unlike Euclidean and Mahalanobis which are fixed bilinear queries and come
  out of one multi-support frame in a single pass. A measured consequence of
  being a nonlinear view.
- **`cosine_cov` is linear CKA** on the double-centred second-moment matrices
  `G = -0.5 H D H`; `corr_cov` is the same on mean-centred vectors.
- **The tied fixture is load-bearing.** `tau-a` separates from Kendall tau-b
  by 0.243 and `rho-a` from Spearman by 0.212 *only* because the model RDMs
  are binary. On untied data four of the eight comparison rows would test the
  same number twice; a test now asserts the separation is nonzero.

---

## Track B — The scaling comparison

**Claim:** on map-scale workloads, crossform produces the same numbers at lower
cost — and where it does not, that is published too.

`vignettes/novelty.Rmd` carries the row *"Cross-package speed advantage — Not
demonstrated"* together with its bar: *matched estimands, independent parity,
warm-up, repeated timings, hardware, and map-scale budgets.*

**Frame it as scaling, not as a race.** `design/vision.md`'s development track 3
says comparative exemplars must "compare matched estimands, **scaling**, and
sampling rules — **not merely similar-looking outputs or runtime**." So the
deliverable is how cost grows with searchlights, conditions, and queried pairs.
A speed ratio is one readout of that, not the headline.

**Be honest about where an advantage can come from.** crossform's kernels are
single-threaded by design (`src/Makevars` sets no OpenMP) against multithreaded
Accelerate numpy, and the package's own profiling recorded `no_rcpp_keep_blas`
with an R-loop share of 0.000 — it is BLAS-bound, not loop-bound. There is no
micro-optimization headroom. Any advantage is architectural:

- **Map scale.** rsatoolbox has no searchlight engine; the idiom is a Python
  loop calling `calc_rdm` per searchlight, re-restricting noise each time.
- **Query-first.** Computing k of q(q−1)/2 pairs when only k are wanted, where
  rsatoolbox must materialize the full RDM.

Both sides already link Apple Accelerate, so the BLAS objection is answered
before it is raised.

**The measurement is the hard part.** The repo's own gates show a **1.3–2×
run-to-run envelope** on identical fixtures (n=3–5, medians only, no thread or
load control; `RECERTIFY.md` documents 16× degradation under contention). So the
harness adds three things the current one lacks:

- **Thread control**, pinned on both sides, with 1-thread *and* all-core
  configurations both reported. Hiding crossform's single-threadedness would
  make the comparison dishonest.
- **Dispersion.** ≥15 repetitions; median, IQR, MAD; and a standing rule that
  **no ratio is reported whose interval crosses 1**.
- **A matched-estimand interlock.** Every timed pair must also agree numerically
  to tolerance *in the same run*, or the timing row is void — this is
  `design/vision.md`'s *"speed never licenses a change in the estimand"* made
  executable.

Provenance gains CPU model, cores, RAM, thread env vars, and load average, none
of which `provenance.R` records today.

**Gate B:** a provenance-bound artifact ≤64 KiB, promoted and ratcheted; the
ledger row moves to what the evidence supports; the axes where crossform loses
are published beside the ones where it wins.

---

## Track C — Publish

`vignettes/reproducing-rsatoolbox.Rmd` carrying the workflow table, the coverage
ledger, the conventions resolved, and the boundaries with their reasons. Update
`novelty.Rmd`'s ledger rows and gate 1 text. README section. One
end-of-program re-certification per `benchmarks/RECERTIFY.md`, after all `R/`
edits, never per ticket.

---

## What stays out, and why

- **Correlation distance as a crossform view.** `rdm(normalize=)` stays refused.
  Reproducing rsatoolbox's within-sample estimator downstream does not license a
  cross-generalized normalizer, which is what the refusal is about.
- **Inference** — noise ceilings, condition bootstrap, the two-factor
  subject × condition bootstrap — stays absent, recorded as `absent_no_claim`.
  `design/relation-to-prior-work.md:61`: "crossform supplies no inference layer,
  and says so."
- **`poisson`, `poisson_cv`, unbalanced designs, RDM movies:** out of scope.
- **Any new export.** The mission measures progress "not by the number of
  exported functions"; nothing here needs one.

## Sequencing

Track A and Track B's harness start immediately and in parallel — A touches only
the exemplar, B needs only existing verbs. Neither blocks the other, because
Track A adds no `R/` code. Track C follows both.

## Risks

| Risk | Handling |
|---|---|
| Measurement envelope swallows the speedup | Dispersion + the interval-clear-of-1 rule; claim only what clears it |
| Speed claim fails on single small ROI | Expected, and published as a loss |
| Untied fixture makes 4 similarity rows vacuous | Verified failure mode; fixture carries deliberate ties |
| Correlation route hits a non-positive diagonal | Precondition checked and refused, per the policy's own contract |

## Board

Umbrella `bd-01M0V36K7PWFN95R5JQW9V82N9`; tracks A/B/C under it. `mote ready`.
