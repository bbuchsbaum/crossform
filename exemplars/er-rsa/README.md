# Rectangular encoding-retrieval RSA

This exemplar closes claim-promotion gate 3 of
[`vignettes/novelty.Rmd`](../../vignettes/novelty.Rmd): *distinct unequal
axes, match/control coding, pair-space covariates, missing items, and
explicit operation order*.

It runs one rectangular encoding-by-retrieval analysis end to end. The left
axis is 36 encoding effects, the right axis is 30 retrieval probes, and the
two are not the same set: 12 studied items are never retrieved and 6
retrieval probes were never studied. Encoding-retrieval couplings are read
with axis-bound pair queries, the item-specific reinstatement contrast is
compiled two ways (a normalized coupling difference and a nuisance-adjusted
weighted least-squares coefficient), and an encoding-trial covariate enters
as a pair-space regression term. The estimand pairs runs across study-test
cycles only.

Three planted regional structures are recovered, qualitatively and
quantitatively, against a closed-form ground truth.

## Dataset decision

**This is a designed simulation with recovered ground truth, not an
empirical dataset.**

A public encoding-retrieval fMRI dataset with per-trial betas — the input
this analysis actually needs — could not be fetched or verified in this
environment. Raw BIDS sets of encoding-retrieval memory experiments exist,
but turning one into per-trial encoding and retrieval betas requires running
first-level GLMs on data this repository cannot download, and shipping an
unverified derivative would put an empirical claim behind a pipeline nobody
here has checked. Ledger discipline says an unverifiable dataset is worse
than an honest simulation, so the demonstration is a simulation whose truth
is known exactly.

What that buys, and what it costs:

- **Buys.** Every readout has an exact planted value, computed in closed
  form from the noiseless per-run patterns (`er_planted_geometry()` in
  `00-common.R`), with no Monte Carlo. Recovery is therefore a numerical
  claim, not an impression: 30 of 30 cross-run readouts land inside their
  95% across-subject interval, and the largest bias is 0.0093 (1.2 standard
  errors).
- **Costs.** Nothing here is evidence about human memory, about any brain
  region, or about how the estimator behaves under real fMRI noise
  (autocorrelation, motion, physiological structure, subject
  heterogeneity). The simulation's noise is spatially correlated but
  temporally independent, so no whitener is supplied and none is claimed.

The upgrade that would retire this boundary is named in the gate: a public
encoding-retrieval dataset with per-trial betas, analysed with the same
scripts.

## How to run

```sh
cd exemplars/er-rsa
Rscript 01-simulate.R    # design, truth, data          (~4 s)
Rscript 02-analyze.R     # the crossform analysis       (~15 s)
Rscript 03-recover.R     # planted vs estimated, tests  (~1 s)
```

Each script is re-runnable and deterministic; running them in order
reproduces every number below. `01-simulate.R` and `02-analyze.R` write
`data/` (git-ignored, ~8 MB); everything committed is in `results/`.
`ER_SUBJECTS=3 Rscript 01-simulate.R` shrinks the run for a smoke test — it
changes the numbers, so the committed results are always produced with the
default 12.

`02-analyze.R` loads crossform from the source tree with
`pkgload::load_all()` when run inside the repository and from the installed
package otherwise. `03-recover.R` calls no crossform function at all: it
reads the committed CSV and does its own statistics, so the recovery check
is external to the thing being checked.

Recorded run: R 4.5.1, aarch64-apple-darwin20, 2026-08-20.

## The design

Three study-test cycles over one 36-item list. Each cycle is one scanning
run holding that cycle's encoding block and then its retrieval block, so a
cycle's encoding and retrieval coefficients come out of the *same* GLM.

| | |
| --- | --- |
| encoding effects (left axis) | 36 (`enc_item01` … `enc_item36`), 12 per category |
| retrieval effects (right axis) | 30 (24 old probes + 6 lures), 10 per category |
| encoded but never retrieved | 12 items |
| probes never encoded | 6 lures |
| study-test cycles (runs) | 3, each 300 TRs at TR = 2 s |
| trials per run | 66 (36 encoding + 30 retrieval), jittered ISI |
| design columns per run | 70 = 66 trial regressors + baseline + 2 drift + motion |
| residual df per run | 230 |
| voxels | 90 (regionA 40, regionB 30, regionC 20) |
| simulated subjects | 12, differing only in noise realization |
| trial covariate | encoding study duration, 1.5–4.0 s, centered at 2.75 s |
| cross-run edges | 6 directed (all ordered cycle pairs `k != l`) |

Repeated study-test cycles are not decoration. A crossform relation
estimates its whole shared effect space in every partition, so each run must
contain every encoding item and every retrieval probe. A single-cycle design
with items split across runs is a different partitioning problem and is out
of scope here.

**Nuisance terms.** Run baseline, two Legendre drift polynomials, and a
motion trace are columns of the same design in both fits and belong to
neither effect space. The simulation puts real drift, baseline, and motion
signal in the data (drift amplitude sd 4, motion sd 2, against effect
amplitudes near 1), so those columns are load-bearing rather than
ornamental.

**Planted structure.** Item patterns, category patterns, and a same-run item
state are drawn once and held fixed across subjects:

| region | what is planted |
| --- | --- |
| regionA | item-specific reinstatement, scaled by the encoding study duration (gain = 0.55 + 0.18 × centered duration), *plus* a category component |
| regionB | category structure only, in both phases |
| regionC | no task structure at all |

regionA deliberately carries category structure too. That is what makes the
naive control set wrong rather than merely imprecise.

**Same-run item state.** Each cycle adds an item-specific pattern (sd 0.5,
all voxels) to that item's encoding trial *and* to its retrieval trial in
the same run. It stands in for anything that couples an item's two trials
inside one run. Cross-cycle pairing removes it; same-run pairing does not.

## Planted pair values, before any estimation

Mean planted coupling per cell class, computed in closed form
(`results/planted-cell-summary.csv`):

| region | edges | match (24) | same-category control (336) | other-category control (720) |
| --- | --- | ---: | ---: | ---: |
| regionA | cross-run | 0.8808 | 0.3531 | −0.0362 |
| regionB | cross-run | 1.4370 | 1.4394 | −0.0642 |
| regionC | cross-run | 0.0062 | −0.0024 | 0.0000 |
| regionA | same-run | 1.1249 | 0.3547 | −0.0369 |
| regionB | same-run | 1.6932 | 1.4377 | −0.0640 |
| regionC | same-run | 0.2537 | 0.0026 | 0.0006 |

Two facts to hold on to. In regionB the match cells and the same-category
controls are *identical* (1.4370 vs 1.4394): there is no item-specific
reinstatement there, only category structure. And in regionC the same-run
match cells sit at 0.2537 while the cross-run ones sit at 0.0062: the
same-run comparator invents an effect out of nothing.

## What crossform estimates

Two relations over one neural domain, then one rectangular plan:

```r
encoding  <- lm_relation_fit(sources, designs, encoding_target, domain = domain)
retrieval <- lm_relation_fit(sources, designs, retrieval_target, domain = domain)

over <- pairing(                       # 6 directed cross-cycle edges
  cross_edges$left, cross_edges$right,
  directed = TRUE, independence = "independent", generalizes_over = "run"
)
plan <- plan_geometry(encoding$relation, frame, over, right = retrieval$relation)
view <- evaluate_geometry(plan, query = contrast_category)
```

Both relations use the same partition names, so `pairing()`'s default
`self_pairs = "forbid"` is what excludes same-run encoding-retrieval
products, and `generalizes_over = "run"` is bound into plan identity. The
same-run comparator is a *different* plan with a different declared estimate
(`self_product_biased` rather than `cross_generalized`); it exists in this
exemplar only to show what the estimand refuses to include.

## Results

All numbers are means over 12 simulated subjects, cross-run plan, with the
planted value beside them (`results/planted-vs-estimated.csv`).

### Match versus control coupling

| region | match coupling | control, all eligible | control, category-matched | item-specific contrast | planted | t |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| regionA | 0.8828 | 0.0890 | 0.3556 | **0.5272** | 0.5277 | 62.1 |
| regionB | 1.4415 | 0.4134 | 1.4487 | **−0.0072** | −0.0024 | −1.1 |
| regionC | 0.0065 | 0.0010 | 0.0003 | **0.0062** | 0.0086 | 1.2 |

The item-specific contrast is `coupling_contrast()` over the
category-matched eligible set. It separates the three regions exactly as
planted: a large effect in regionA, nothing in regionB, nothing in regionC.

### The control set is a scientific choice, not a formality

| region | contrast vs all eligible controls | planted | contrast vs category-matched controls | planted |
| --- | ---: | ---: | ---: | ---: |
| regionA | 0.7937 | 0.7931 | 0.5272 | 0.5277 |
| regionB | **1.0281** | 1.0228 | −0.0072 | −0.0024 |
| regionC | 0.0054 | 0.0070 | 0.0062 | 0.0086 |

Matched pairs are always same-category, so an unrestricted control set
compares within-category matches against mostly across-category controls.
Under that control set regionB — which has no item-specific structure
whatsoever — reports a reinstatement effect of 1.03 (t = 104). Restricting
the eligible set to same-category pairs, through `match_coupling(…,
eligible = )`, removes it. Both numbers are recovered to their planted
values, so this is a property of the estimand, not of the estimator.

### The trial covariate

`pair_lm_query()` over all 1080 eligible pairs, with `match`, `match_depth`
(the centered encoding study duration on matched cells), `same_category`,
and both item nuisance families:

| region | `match` | planted | `match_depth` | planted |
| --- | ---: | ---: | ---: | ---: |
| regionA | 0.5329 | 0.5341 | **0.2206** | 0.2154 |
| regionB | −0.0080 | −0.0042 | 0.0544 | 0.0602 |
| regionC | 0.0083 | 0.0085 | 0.0001 | −0.0044 |

The regression recovers both the item-specific reinstatement (0.5329, next
to `coupling_contrast()`'s 0.5272 and `match_control()`'s 0.5166 on the
restricted set) and the study-duration slope, while adjusting for the
category confound over the *full* eligible set rather than discarding cells.

regionB's `match_depth` of 0.054 is not a failure: its planted value is
0.060. The realized same-run item state leaves a small item-level residue in
the planted truth, and the estimator reports it. That is the point of
comparing against a closed-form planted value rather than against zero.

### Cross-run generalization has a measurable consequence

Same estimand, same queries, same data — only the pairing changes:

| region | cross-cycle edges | planted | same-run edges | planted |
| --- | ---: | ---: | ---: | ---: |
| regionA | 0.5272 | 0.5277 | 0.7740 | 0.7702 |
| regionB | −0.0072 | −0.0024 | 0.2710 | 0.2555 |
| regionC | 0.0062 | 0.0086 | **0.2560** | 0.2511 |

Under same-run pairing, regionC — which contains no task structure at all —
reports item-specific reinstatement of 0.256 (t = 67), and regionA is
inflated by 47%. Both plans recover their own planted values; they are
simply different quantities, and crossform makes them different objects.

### Recovery, execution, and boundaries

- **Recovery.** 30 of 30 cross-run readouts have their planted value inside
  the 95% across-subject interval. Largest bias 0.0093, largest |bias|/se
  1.2 (`results/planted-vs-estimated.csv`).
- **Route.** Query-first `evaluate_geometry()` and materialize-then-project
  `query_geometry()` agree to `6.66e-16` over all ten queries, and the
  materialized rectangular form satisfies `total = coherent + configuration`
  to `2.22e-16` (`results/route-check.csv`).
- **Verdicts.** All eight qualitative claims pass their stated test
  (`results/recovery-verdicts.csv`), and
  [`tests/testthat/test-er-rsa-exemplar.R`](../../tests/testthat/test-er-rsa-exemplar.R)
  ratchets the record.

Three boundaries were hit and recorded in `results/refusals.csv`:

| situation | capability | outcome |
| --- | --- | --- |
| a fixed noise metric on the rectangular plan | `rectangular_fixed_metric` | refused; this analysis runs on the implicit identity metric |
| `rdm()` on a rectangular cross-axis plan | `symmetric_self_form` | refused; 36 × 30 is not a square RDM |
| `match_control()` with both nuisance families on category-blocked cells | rank deficiency | classed input error; the exemplar drops the retrieval family explicitly |

The first is the live limitation. Rectangular plans refuse a fixed noise
metric, so there is no crossnobis-style noise normalization here and no
analytic sampling covariance for a rectangular readout; the uncertainty
reported above is the across-subject spread over 12 independent noise
realizations. Admitting a fixed metric on rectangular plans is ticket B7.

The third is a real identifiability fact, not a crossform quirk: on a
category-blocked eligible set the encoding and retrieval level shifts trade
off exactly within each block, so the two nuisance families are not jointly
estimable. crossform refuses rather than silently dropping a column.

## What each of the five pair-query exports contributed

The novelty ledger keeps these exports only if a real analysis needs them
(`design/api-tiers.md`, maintainer decision 3). All five are load-bearing
here.

| export | what it did that nothing else did |
| --- | --- |
| `match_coupling()` | declared the encoding-retrieval correspondence over unequal axes with 12 unmatched left items and 6 unmatched right items, and carried the `eligible` restriction that turns the naive 1.03 in regionB into −0.01 |
| `control_coupling()` | named the comparison cells, in both variants: the 1056/336 non-matched controls, and `include_matches = TRUE` for the eligible-pair baseline level (regionA 0.107, regionB 0.436) |
| `coupling_contrast()` | compiled the headline item-specific readout as one fixed operator; it is also the constructor whose balance diagnostics reveal the problem below |
| `match_control()` | the same comparison as a nuisance-adjusted coefficient (0.5166 against `coupling_contrast()`'s 0.5272), additive-baseline invariant where the coupling difference is not, and the constructor whose nuisance flags had to be set explicitly when both families turned out not to be jointly identified |
| `pair_lm_query()` | the covariate model: `match_depth` (0.2206 vs planted 0.2154) alongside `match`, `same_category`, and 64 item nuisance columns, none of which the coupling constructors can express |

The clearest evidence that `match_control()` and `pair_lm_query()` are not
redundant with `coupling_contrast()` is in `results/query-diagnostics.csv`.
On this rectangular design — unequal axes, missing matches — the normalized
coupling difference is **not** additive-baseline invariant: its row and
column marginals are nonzero, because a retrieved item contributes `1/24` of
the matched mass and `29/1056` of the control mass, which do not cancel. The
regression forms have exactly zero marginals (`rank == columns`, 68 of 68),
so a per-item or per-probe level shift cannot masquerade as reinstatement.
On the square balanced example in `?coupling_contrast` the two agree; on a
real rectangular design they do not, and the difference is visible in the
compiled operator before any data is read.

## Files

| file | what it is |
| --- | --- |
| `00-common.R` | design constants, the simulation, and the closed-form planted geometry |
| `01-simulate.R` | builds design, truth, and data; writes the ground-truth tables |
| `02-analyze.R` | the crossform analysis: relations, plan, ten queries, two pairings, 12 subjects |
| `03-recover.R` | planted vs estimated, the match/control tests, the verdicts |
| `results/design-items.csv`, `design-probes.csv`, `design-summary.csv` | the design as data |
| `results/planted-cell-summary.csv` | planted pair values by region and cell class |
| `results/query-diagnostics.csv` | what each compiled operator is: cells, rank, balance, claim |
| `results/estimates-by-subject.csv` | every subject × plan × query × region estimate, with its planted value |
| `results/planted-vs-estimated.csv` | the recovery table |
| `results/coupling-levels.csv` | the headline match/control levels |
| `results/route-check.csv` | query-first vs materialized, and rectangular recomposition |
| `results/recovery-verdicts.csv` | the eight qualitative claims and their tests |

## Boundary statement

This is a designed simulation with recovered ground truth. It demonstrates
that crossform's rectangular machinery — unequal directed axes, missing
matches, eligibility restriction, match/control coupling, pair-space
regression with a trial covariate and item nuisance families, and an
estimand that generalizes across runs — estimates what it claims to
estimate, on a design where the answer is known. It is not evidence about
memory, about any brain region, or about behaviour under real fMRI noise.
No empirical encoding-retrieval result is claimed here, and none should be
read out of it.
