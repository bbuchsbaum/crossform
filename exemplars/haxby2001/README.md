# Exemplar: Haxby 2001, matched-estimand searchlight RSA

Status: scaffold under construction
Date: 2026-08-14
Vision track: 3 (comparative exemplars) — compare matched estimands, scaling,
and sampling rules on public data; not merely similar-looking outputs.

## Purpose

This is the first public-data, matched-estimand comparison between a
`crossform` workflow and its `rMVPA` counterpart. It exists to earn — or
honestly fail to earn — one row of the replacement map: *searchlight
representational similarity analysis with cross-run generalization*.

The exemplar succeeds when either

1. both systems, given identical inputs and a matched estimand, agree
   numerically within a stated tolerance at every searchlight center, or
2. they disagree, and the disagreement is traced to a named semantic
   difference (pairing rule, distance convention, normalization, tie handling,
   center inclusion) that is documented here rather than left as an anecdote.

A silent numerical discrepancy is a failure of the exemplar, not a footnote.

## Design principle: isolate the analysis semantics

Both systems consume **identical per-run condition-mean matrices** computed
once by a shared preparation script. No GLM machinery, HRF model, or
preprocessing choice is allowed to differ between arms; the comparison
isolates the pairing, distance, and RSA semantics that the two packages
actually own. Richer beta estimation (fmrireg) is a planned v2 arm, applied
identically to both systems when added.

## Dataset

Haxby et al. (2001), subject 1, from the PyMVPA distribution
(`http://data.pymvpa.org/datasets/haxby2001/`): 12 runs, block design, 8
stimulus categories (face, house, cat, bottle, scissors, shoe, chair,
scrambled), with the ventral-temporal (VT) mask shipped alongside the BOLD
series. This is the canonical public dataset for category-structure RSA and
is small enough (~300 MB) to download and analyze in one sitting.

## Estimand (stated once, used twice)

For each searchlight sphere (radius 3 voxels, centers restricted to the VT
mask for the smoke tier; whole brain for the full tier):

- **Inputs:** per-run condition-mean patterns, one 8 × P matrix per run,
  where P is the number of in-sphere voxels; each run's patterns are computed
  from within-run z-scored time series averaged over HRF-shifted block
  volumes by `01-prepare-data.R`.
- **Effect estimate:** the cross-run second-moment structure over the 8
  conditions, restricted to **off-diagonal run pairs only** (an effect must
  reproduce across runs; within-run structure never contributes).
- **RDM:** correlation distance (1 − Pearson r) between condition patterns,
  computed from the cross-run structure.
- **Report:** Spearman correlation between the lower triangle of the
  empirical RDM and the binary animate/inanimate + category model RDM
  defined in `models.R`.

In `rMVPA` this is an `rsa_design`/`rsa_model` searchlight with run-wise
cross-validation structure. In `crossform` it is a relation over the eight
condition effects with `cross_partitions()` over runs, a searchlight frame,
and `rdm()`/`rsa()` views. Where a system cannot express exactly this
estimand, the nearest expressible estimand is documented and the residual
difference becomes part of the comparison, not noise.

## Beyond parity: the error-channel arm

Parity is the floor, not the point. A third arm runs `crossform`'s
crossnobis RDM with the identified error channel from `lm_relation_fit()`,
and transports sampling covariance to the RSA coefficient — something
`rMVPA` does not expose. This arm demonstrates the scientific surplus of the
new system on the same data, clearly labeled as a *different, better-defined
estimand* rather than a comparison row.

## Scripts

| script | role |
|---|---|
| `00-download.R` | fetch and checksum subject 1 tarball into `data/` (gitignored) |
| `01-prepare-data.R` | one shared preparation: per-run condition means, VT + whole-brain masks, saved as `.rds` |
| `02-crossform-searchlight.R` | the crossform arm |
| `03-rmvpa-searchlight.R` | the rMVPA arm |
| `04-compare.R` | center-by-center agreement report, tolerance gates, divergence taxonomy |
| `05-crossnobis-uncertainty.R` | the error-channel arm (crossform only) |
| `06-coherent-configuration.R` | the coherent/configuration arm (crossform only) |
| `models.R` | model RDMs and shared constants (radius, tolerance, seeds) |

Run order is numeric. Every script is idempotent and writes only under
`data/` (raw + prepared inputs, gitignored) and `results/` (small summary
`.rds` and the comparison report, committed as evidence).

## Tiers

- **Smoke:** VT-mask centers only (~500–600 voxels), minutes not hours.
  This tier gates commits to this directory.
- **Full:** whole-brain centers, both arms, wall-clock and peak memory
  recorded per the benchmark conventions in `benchmarks/README.md`.

## Current status

Updated 2026-08-13. **Smoke tier implemented and run end to end.** All seven
scripts exist and are idempotent; `results/smoke-comparison.rds` and
`results/smoke-report.md` are the evidence. Read the report for the full
argument — this section records the outcome and the one place where reality
contradicted the estimand written above.

The committed smoke report predates the package rename and therefore retains
`effectagram` in its historical receipt text. Rerunning the current scripts
writes `crossform`-labeled artifacts; the earlier evidence is not rewritten.

### Smoke tier result: matched estimand agrees, README estimand is not expressible

Both arms consume the identical prepared object, and their searchlight
geometry was verified rather than assumed: `neuroim2::searchlight_indices`
(crossform) and `neuroim2::searchlight` (rMVPA) return **bit-identical
sphere membership at all 577 VT centers** at radius 11.25 mm.

| comparison | result |
|---|---|
| matched estimand, **crossform vs rMVPA's own crossnobis estimator** | max abs difference **8.88e-16**, 0/577 centers above the 1e-8 gate — PASS |
| matched estimand, crossform vs an independent reference loop | max abs difference 1.33e-15, 0/577 — PASS |
| shared Spearman score | max abs difference **0** (exact), map r = 1 — PASS |
| native estimands (crossform vs rMVPA `rsa_model`) | r ≈ 0.54, fully accounted for by five named semantic differences |

### The estimand section above is not achievable as written

The "Estimand" section specifies a **correlation-distance (1 − Pearson r) RDM
scored by a Spearman correlation**. Implementation established that neither
system can express that combination:

- **crossform** cannot produce a correlation-distance RDM at all. Correlation
  distance rescales each pattern to unit norm inside the sphere, which is
  nonlinear in the patterns and therefore outside the bilinear core.
  (Per-sphere mean-centering *would* be expressible as a non-diagonal metric;
  the unit-norm denominator is not.) crossform's `rsa()` is also OLS, not a
  rank correlation.
- **rMVPA's `rsa_model`** couples the two choices: `distmethod` sets both the
  neural RDM metric and the second-order correlation method, and `regtype` is
  inert for `"pearson"`/`"spearman"` — verified on this data, max abs
  difference exactly 0 across 577 centers between the two `regtype` values.
  (`vector_rsa_model` does decouple them, but it averages 96 per-trial row
  correlations rather than correlating over a condition RDM's lower triangle,
  so it is a different estimand again.)

The nearest estimand **both** systems express exactly, and what the smoke tier
therefore compares, is a condition-level cross-run RDM using the
**cross-validated squared Euclidean distance** over the 66 off-diagonal run
pairs, plus a Spearman statistic supplied by the exemplar's own
`spearman_rsa()` applied identically to both arms. The residual difference from
the text above is the distance convention, and it is a hard expressibility
boundary rather than a configuration choice.

Per the design principle in this README, the Purpose/Design/Estimand sections
are left unedited; this note records where they and the implementation diverge.

### Replacement-map claim

**Earned for the matched estimand, with two conditions attached.** Both sides
of the headline comparison are package estimators — crossform's `rdm()` over
a cross-run geometry against rMVPA's
`compute_crossvalidated_means_sl(estimation_method = "crossnobis")` +
`compute_crossnobis_distances_sl()` — and they agree to 8.88e-16 at every
center.

The conditions:

1. The rMVPA side is **not** `rsa_model`, the function this README originally
   named for the arm. `rsa_model` cannot express the estimand (it never
   aggregates to condition level), and neither can `vector_rsa_model` (trial
   level, returns one scalar) or `contrast_rsa_model` (computes the right
   `G_cv` internally but exposes only contrast betas). The crossnobis helper is
   the only route, and `compute_crossnobis_distances_sl` is **not exported**,
   so `03` calls it with `:::`.
2. It requires a current rMVPA. The 2026-04 build refuses crossnobis without a
   whitening matrix and, given one, derives folds from leave-one-run-out
   *training* means; those overlap across folds, so noise does not cancel and
   the distances come out positively biased. `03` therefore checks the
   **semantics, not the version**: it asserts that the returned fold estimates
   equal this run's per-run condition means and stops with an explanatory
   error otherwise. On the build used here that assertion passes at max abs
   difference 0, and crossnobis on white noise averages −0.011 (unbiased),
   confirming the folds are independent.

The row this earns is the crossnobis/cross-validated-distance estimand, not
correlation-distance searchlight RSA — no system in this comparison computes
the latter the way the Estimand section describes.

### Cost

The parity tier (scripts 02–04) is genuinely cheap: about 20 s total, of which
crossform's `rdm()` over 577 centers is 0.48 s and rMVPA's crossnobis RDM is
0.40 s. The error-channel arm (05) is heavier, but the old estimate in this
section was stale. A 2026-08-14 rerun of 25 evenly spaced nodes on the current
tree took 16.2 s total (0.65 s/node). A linear projection is therefore about
78 s for the 120-node default and 6.3 min for `UNCERTAINTY_NODES=all`; those two
numbers are projections, not measured full-sweep gates. The separately
versioned public map gate measures a complete 576-node, 66-distance covariance
sweep at 183.3 s on its synthetic scale fixture. Scripts 02–04 remain the
quick parity tier; script 05 is an explicitly timed extended gate.

The committed `smoke-uncertainty.rds` predates the validation-memo optimization
and still records 9.03 s/node. Keep it as historical result evidence, not a
current runtime claim; rerunning script 05 writes the current per-node timings.

### Error-channel arm: ran end to end

`05` refits the same effects from raw volumes with `lm_relation_fit()` —
reproducing `02`'s RDM to 4.44e-16, so the point estimand is unchanged — and
then does what rMVPA cannot. Over 120 evenly spaced VT centers: RDM standard
errors (median 0.0207), and an **exact** transport of the RDM sampling
covariance to the RSA coefficient (median estimate 0.0372, median SE 0.0121,
median |z| 3.47, 72.5 % of centers above |z| = 1.96). The transport is exact
rather than approximate because `rsa()` is a linear functional of the
distances, verified to 1.1e-16. Learned-metric crossnobis (face − house) over
all 577 centers: median 0.300, positive at 99.8 %.

All three capability **refusals** fired and are recorded verbatim in the
report: no analytic covariance on a learned-metric plan, none from precomputed
effects with no error channel, and no inferred `target`. The second one is why
this arm had to refit from raw volumes rather than reuse `01`'s means.

### Environment hazard, needs a maintainer decision

A subagent installed rMVPA from source into the user library during this
session (`~/Library/R/arm64/4.5/library/rMVPA`, built 2026-08-14 00:59 UTC).
It now shadows the pre-existing framework build (2026-04-25). This was not
requested; it was left in place rather than deleted, and needs an explicit
accept-or-undo. It is load-bearing for the results: the crossnobis arm only
works on the newer build, and `rsa_model` with `distmethod="spearman"` differs
between builds by 2.97e-4 (`"pearson"` agrees to 1.11e-16).

### Coherent/configuration arm: the central claim on real data

Added 2026-08-16. `06-coherent-configuration.R` runs the decomposition that
had previously only been shown on a generated fixture. Same prepared object,
same VT searchlight frame, same cross-run pairing as `02`; the only new thing
is what is read off the plan. Two contrasts:

- **face − house**, weights `+1 / −1`;
- **animate − inanimate**, `(cat, face)` at `+1/2` each against
  `(bottle, chair, house, scissors, shoe)` at `−1/5` each. `scrambledpix` is a
  visual control rather than an object category and carries weight 0. The
  animacy coding is `models.R`'s `ANIMATE`, not a new one.

The partition is exact where it must be: over 577 searchlights,
`max |total − (coherent + configuration)|` is **5.55e-17** for face − house and
**1.39e-17** for animate − inanimate.

| face − house | value |
|---|---|
| peak searchlight (position 573, centre voxel 77745) | signed **−1.0694**, coherent **1.1390**, configuration **0.0902**, total **1.2292**, coherent share **0.927** |
| coherent share over VT (536 of 577 searchlights with a valid fraction) | median **0.529**, IQR **[0.278, 0.773]**; above 0.5 at 281 of 536 |
| total > 0 | **576 / 577** (99.8 %) |
| configuration > 0 | 553 / 577 (95.8 %) |
| coherent > 0 | 559 / 577 (96.9 %) |
| signed < 0 | 568 / 577 |

| animate − inanimate | value |
|---|---|
| peak searchlight (position 278, centre voxel 67417) | signed **−0.3786**, coherent **0.1420**, configuration **0.0556**, total **0.1977**, coherent share **0.719** |
| coherent share over VT (503 of 577 valid) | median **0.442**, IQR **[0.200, 0.659]** |
| total > 0 | 575 / 577 (99.7 %) |
| configuration > 0 | 554 / 577 (96.0 %) |
| coherent > 0 | 526 / 577 (91.2 %) |

Read at a single whole-VT region (`regions()` over all 577 VT voxels), from
the same relation and pairing:

| contrast | signed | coherent | configuration | total | coherent share |
|---|---|---|---|---|---|
| face − house | −0.46292 | 0.21220 | 0.22371 | 0.43592 | 0.4868 |
| animate − inanimate | −0.16541 | 0.02662 | 0.06287 | 0.08949 | 0.2975 |

Three things are worth reading out of that.

1. **A substantial share of the reproducible energy is coherent.** In this
   subject's VT, the median searchlight puts about half of the crossvalidated
   face/house energy in its own weighted common spatial mode. An analysis that
   demeans each searchlight and reports the remainder would discard it; one
   that reports only the total would not know it was there.
2. **The signed marginal carries the direction, and both directions are
   negative.** For face − house, `signed` is negative at 568 of 577
   searchlights and at the whole-VT reading, so the common mode runs *house
   above face*; animate − inanimate is negative at 512 of 577 and at the
   whole-VT reading, *inanimate above animate*. Energies are squared
   quantities and cannot say this; keeping the first-moment marginal in the
   same plan is what makes the coherent part interpretable rather than merely
   large.
3. **The coherent part is the crossvalidated version of that mean
   difference.** Across searchlights, `coherent` correlates with `signed²` at
   r = 0.9999 (face − house) and 0.9994 (animate − inanimate), with a median
   ratio of 0.964 and 0.942 — just below one, which is what a bias-free
   cross-partition product should be relative to a squared average that still
   contains its own noise.

Cost: 0.41 s for face − house and 0.49 s for animate − inanimate over 577
searchlights, 0.16 s for both whole-VT region reads, 1.8 s of analysis, and
3.4 s of wall clock including R startup and reading the mask. Evidence:
`results/coherent-configuration.rds` and `results/coherent-configuration.png`.

**Caveats, which are not small.** One subject. A block design. Per-run
condition means, not GLM betas, computed from within-run z-scored time series
with a 2 TR haemodynamic shift. An identity metric, not a whitened or learned
one. And **no inference of any kind**: no permutation test, no group model, no
p-values, no correction. Every number above is a description of one subject's
ventral temporal cortex, and the claim it supports is that the decomposition
runs exactly and informatively on real data — not anything about faces,
houses, or animacy.

### Not done

- Full (whole-brain) tier: out of scope here, and not run.
- Wall-clock and peak-memory instrumentation per `benchmarks/README.md`: the
  timings recorded are plain elapsed times from single runs, not benchmarks.
- A v2 arm with richer beta estimation (fmrireg) applied identically to both
  systems.
