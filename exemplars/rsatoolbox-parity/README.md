# `rsatoolbox` parity, and what crossform reads from the same fit

This exemplar closes claim-promotion gate 1 of
[`vignettes/novelty.Rmd`](../../vignettes/novelty.Rmd): *reproduce a fixed
crossnobis plus linear-RSA result with a version-pinned environment and an
independent oracle*.

It has two arms.

1. **Parity.** On one shared synthetic fixture, crossform's fixed-metric
   crossnobis RDM and its linear RSA coefficients agree with Python
   `rsatoolbox` 0.3.2 to **3.8e-15** and **8.6e-16** in maximum absolute
   difference. Those are floating-point reordering differences, not
   estimator differences.
2. **Strict extension.** From the *same* relation fit, the *same* fixed noise
   precision, and the *same* cross-run pairing, crossform additionally
   returns the signed contrast energy, the exact coherent/configuration/total
   partition of it, the same partition of every RDM entry, and analytic
   sampling covariance for the RDM and for a linear RSA coefficient, under a
   declared and refusable model.

Nothing here claims that `rsatoolbox` is wrong or slow. The parity arm exists
so that the extension arm is anchored: the extra quantities come out of an
object that reproduces the standard answer exactly.

The scope is deliberately narrow: the external result supports the fixed
crossnobis and fixed linear-RSA rows described below, not every statistic
called RSA. The exact mapping into the common-geometry theorem is recorded in
[`design/common-geometry-equivalence.md`](../../design/common-geometry-equivalence.md#8-external-parity-binding).

## How to run

R first, then Python, then R.

```sh
cd exemplars/rsatoolbox-parity

# 1. Environment (once). Any Python 3.12 with the pinned wheels will do.
uv venv -p 3.12 rsaenv
uv pip install --python rsaenv/bin/python -r requirements.txt

# 2. The four steps.
Rscript 01-fixture.R          # build the fixture, fit it, export CSVs
rsaenv/bin/python 02-rsatoolbox.py
Rscript 03-compare.R          # agreement table -> results/agreement.csv
Rscript 04-extension.R        # the strict extension
Rscript 05-manifest.R         # bind sources, environment, and outputs

# or, equivalently
RSA_PYTHON=rsaenv/bin/python ./run-all.sh
```

`01-fixture.R` loads crossform from the source tree with `pkgload::load_all()`
when run inside the repository, and from the installed package otherwise.
There is no `reticulate` anywhere: the only channel between the two languages
is the CSV directory `results/`.

Recorded run: R 4.5, CPython 3.12.11, `rsatoolbox` 0.3.2, `numpy` 2.5.2,
`scipy` 1.18.0, darwin/arm64, 2026-08-17.

To revalidate the recorded external output against a changed crossform source
without rebuilding Python, run `Rscript 01-fixture.R`, `Rscript 03-compare.R`,
and `Rscript 05-manifest.R`. This regenerates the deterministic crossform arm,
rechecks every matched value against the pinned external CSVs, and refreshes
the source binding. A full external regeneration uses `run-all.sh` and the
pinned environment above.

`results/parity-manifest.csv` is the drift diagnostic. It records byte sizes
and MD5 digests for the fixture and comparison sources, Python implementation,
environment lock, algebraic claim, machine-readable fixture contract, and
recorded outputs. `tests/testthat/test-rsatoolbox-parity.R` recomputes every
entry, so editing a producer or an output without regenerating and reviewing
the parity evidence fails CI/certification. A digest mismatch identifies the
exact stale path; a semantic assertion then diagnoses version, tolerance,
metric, centering, partition, ordering, or objective drift.

## The fixture

Deterministic (`set.seed(20260817)`), and small enough to inspect by hand:

| | |
| --- | --- |
| conditions | 6 (`cond1_faceA` … `cond6_toolB`) |
| runs (cross-validation folds) | 4 |
| voxels | 40 |
| observations per run | 48 (8 trials per condition, indicator design) |
| residual df | 168 = 4 × (48 − 6) |
| residual covariance | non-spherical: AR-like correlation × heterogeneous voxel scale + a low-rank bump; condition number **182** |
| measurements | one whole-brain (40 voxels) + three regions (16 / 14 / 10 voxels) |

The condition names are chosen so that alphabetical order — what
`rsatoolbox`'s `sort_by()` and `np.unique()` impose — is the same as
crossform's declared effect order. The residual covariance is deliberately
far from spherical: an identity metric would give visibly different numbers,
so the noise metric is load-bearing rather than decorative.

The betas exported to Python are crossform's relation blocks, checked against
plain per-run OLS on the same design (max abs diff **6.7e-16**) before they
leave R. Both packages therefore contract the same estimates.

## Conventions, and how each was resolved

A parity claim is mostly a claim about conventions. Four had to be matched,
and all four were resolved exactly — none needed a loosened tolerance.

**1. The noise precision.** `crossform::noise_precision()` is a *constructor
for a fixed metric*, not an estimator: it takes the precision the analyst
supplies and records that it was fixed before effect evaluation. So the
estimator is the exemplar's, applied once in R — pooled within-run residual
covariance, `sum_r E_r' E_r / nu` with `nu = 168` — and the identical 40 × 40
matrix is handed to `rsatoolbox` as `noise=`. To show that this is not a
private convention, `02-rsatoolbox.py` also recomputes it with
`rsatoolbox.data.noise.prec_from_residuals(residuals, dof=168, method="full")`
and gets the same matrix to **2.3e-14** (relative 3.6e-15).

**2. The cross-validation folding.** `calc_rdm_crossnobis` crossvalidates
*leave-one-fold-out*: for each fold it contracts that fold's condition means
against the mean of the remaining folds. `crossform::cross_partitions()`
instead declares the C(P, 2) = 6 unordered run pairs with uniform weight.
These are the same estimator on a balanced design: LOO is
`(1/P) Σ_f (1/(P−1)) Σ_{g≠f}`, which is the uniform mean over ordered pairs,
and the fixed precision is symmetric, so ordered and unordered means coincide.
`02-rsatoolbox.py` computes the explicit all-pairs mean in plain numpy as a
third oracle rather than asserting the identity: LOO vs all-pairs agree to
**3.3e-16**.

**3. The division by channel count.** `_calc_rdm_crossnobis_single` divides by
`n_channels`. crossform gets the same division from the *frame*, not from the
distance: `whole_brain()` and `regions()` with `normalization = "local"` have
row sums of one, so each measurement's weights are `1/|support|` and the
support-streamed kernel contracts `sqrt(w_i) K_ij sqrt(w_j)`. Same number,
arrived at from opposite directions. For a region, both sides use the
precision *submatrix* `K[support, support]` — not the inverse of the
covariance submatrix — which is what crossform's kernel contracts.

**4. The RSA objective.** `crossform::rsa()` is ordinary least squares in RDM
space: the models are vectorised over the row-major upper triangle
(`utils::combn(seq_len(q), 2)`, which is `np.triu_indices(q, 1)`), an
intercept column is added by default, and the coefficient map is compiled
into the geometry query. The like-for-like comparator is therefore
`numpy.linalg.lstsq` on the vectorised RDMs. `rsatoolbox`'s own
`ModelWeighted` + `fit_regress` is a *different objective*: with
`method="cosine"` it first divides the data vector by its RMS
(`pool_rdm(..., "cosine")`), so its `theta` is the no-intercept OLS
coefficient divided by `sqrt(mean(d²))`, and by default it renormalises
`theta` to unit norm as well. That is a rescaling, not a different fit, and
`02-rsatoolbox.py` records the scale factor and undoes it. Both routes appear
in the agreement table.

## Agreement table

`results/agreement.csv`, verbatim (`03-compare.R` fails the run if any row
exceeds its tolerance):

| quantity | comparator | n | max abs diff | max rel diff | tolerance |
| --- | --- | ---: | ---: | ---: | ---: |
| `crossnobis_rdm` | `rsatoolbox::calc_rdm_crossnobis` | 60 | 3.77e-15 | 6.70e-14 | 1e-10 |
| `crossnobis_rdm` | explicit all-pairs numpy oracle | 60 | 3.77e-15 | 6.69e-14 | 1e-10 |
| `crossnobis_rdm[roiA]` | `rsatoolbox::calc_rdm_crossnobis` | 15 | 3.77e-15 | 9.19e-15 | 1e-10 |
| `crossnobis_rdm[roiB]` | `rsatoolbox::calc_rdm_crossnobis` | 15 | 1.03e-15 | 7.23e-15 | 1e-10 |
| `crossnobis_rdm[roiC]` | `rsatoolbox::calc_rdm_crossnobis` | 15 | 9.44e-16 | 6.70e-14 | 1e-10 |
| `crossnobis_rdm[whole_brain]` | `rsatoolbox::calc_rdm_crossnobis` | 15 | 7.77e-16 | 4.92e-14 | 1e-10 |
| `linear_rsa_coefficients` | numpy least squares on vectorised RDMs | 20 | 8.60e-16 | 1.11e-12 | 1e-10 |
| `linear_rsa_coefficients[intercept]` | numpy least squares on vectorised RDMs | 12 | 4.94e-16 | 1.11e-12 | 1e-10 |
| `linear_rsa_coefficients[no intercept]` | `ModelWeighted` + `fit_regress(cosine)`, rescaled | 8 | 8.88e-16 | 1.82e-14 | 1e-8 |

Plus, from `results/rsatoolbox-meta.csv`:

| check | value |
| --- | ---: |
| `prec_from_residuals(method="full", dof=168)` vs the R pooled precision | 2.31e-14 |
| LOO folding vs uniform C(4,2) pairing | 3.33e-16 |

**No convention had to be conceded.** The declared tolerance is `1e-10` and
the worst observed difference across all 60 RDM entries and 20 RSA
coefficients is **3.8e-15**, which is the accumulation order of the two
implementations' inner products, not a difference in estimand. The one
looser row (`1e-8`) is the `fit_regress` route, and it is looser only because
its coefficient passes through an extra division and multiplication by the
data RMS; it in fact lands at 8.9e-16 too.

The largest relative difference, `1.1e-12`, belongs to a near-zero quantity
(the roiC RSA intercept, `−9.3e-05`); its absolute difference is 4.9e-16.

## The strict extension

Every number below comes from the same `plan_geometry(metric =
noise_precision(...))` object that produced the parity RDM. Nothing is
refitted. Contrast: `mean(face) − mean(house)`.

`results/extension-table.csv`, whole-brain measurement:

| quantity | `rsatoolbox` | crossform, from the same fit |
| --- | --- | --- |
| crossnobis RDM (15 pairs) | yes (`calc_rdm_crossnobis`) | yes, agrees to 3.77e-15 |
| linear RSA coefficient (`category`) | yes (`fit_regress` / external OLS) | yes, one compiled query, agrees to 8.60e-16 |
| signed contrast energy, face − house | no — RDM entries are squared and unsigned | **+0.230130** |
| coherent energy (common spatial mode) | no | 0.008804 |
| configuration energy (pattern remainder) | no | 0.198452 |
| total energy = coherent + configuration | no | 0.207256 |
| coherence fraction | no | 0.0425 |
| coherent/configuration split of *every* RDM entry | no | yes, e.g. `cond1_faceA` vs `cond2_faceB`: 0.006279 + (−0.002832) |
| analytic SE of one RDM entry | bootstrap / subject-level inference | 0.029394 (plugin) |
| analytic SE of the RSA coefficient | bootstrap / subject-level inference | 0.026697 (plugin) / 0.008601 (null) |
| refusal when the estimand is not admitted | no typed refusal channel | 3 of 3 provoked refusals are classed conditions |

Three internal identities are asserted, not assumed, in `04-extension.R`:

* `crossnobis(plan, w)` equals `contrast_energy(plan, w)$total` **exactly**
  (difference `0`) — one estimand, two named views;
* `coherent + configuration − total` = **0** for the contrast and ≤ 2.8e-17
  for every RDM entry — the partition is exact, not approximate;
* every RDM entry re-derived as `crossnobis(plan, e_i − e_j)` reproduces
  `rdm(plan)` to ≤ 2.2e-16 — the RDM really is one linear functional of the
  same evidence, not a separately computed product.

Across the four measurements (`results/extension.csv`):

| measurement | signed | coherent | configuration | total | coherence fraction |
| --- | ---: | ---: | ---: | ---: | ---: |
| whole_brain | +0.2301 | 0.008804 | 0.19845 | 0.20726 | 0.0425 |
| roiA | +0.6546 | 0.219423 | 0.30787 | 0.52730 | 0.4161 |
| roiB | +0.0270 | −0.002776 | −0.02074 | −0.02352 | *not reported* |
| roiC | −0.1646 | 0.003122 | 0.02025 | 0.02337 | 0.1336 |

Two things in that table are the point of the arm. `roiC` has a **negative**
signed contrast alongside a positive energy: the squared distance cannot say
which way the effect goes, and the RDM entry alone would report the region as
weakly "different" with no direction. And `roiB`'s coherence fraction is
withheld rather than printed, because its raw cross-generalized components do
not form a nonnegative partition — the ratio would be an artefact, so
`$coherence_fraction_valid` is `FALSE` there.

Analytic uncertainty (`results/extension-uncertainty.csv`), from
`rdm_sampling_covariance()` and its exact transport through the fixed linear
RSA readout:

| measurement | RSA `category` | SE (null) | SE (plugin) | ν | P_eff |
| --- | ---: | ---: | ---: | ---: | ---: |
| whole_brain | 0.1173 | 0.00860 | 0.02670 | 168 | 52.81 |
| roiA | 0.4536 | 0.02201 | 0.10345 | 168 | 14.64 |
| roiB | 0.0782 | 0.02559 | 0.07319 | 168 | 12.20 |
| roiC | 0.0431 | 0.02698 | 0.04936 | 168 | 9.76 |

The RSA coefficient is a *fixed linear functional* of the 15 distances, so
its variance is an exact transport of the RDM sampling covariance rather than
a resampling estimate; `sampling_covariance(cv, "transport", query = ...)`
performs it without materialising the 15 × 15 matrix. `plugin` and `null` are
different declared calibration policies, never chosen implicitly, and the
residual channel behind the quadratic noise term is reported (ν degrees of
freedom against P_eff effective directions) rather than assumed.

### Boundaries

* One simulated subject, no group-level inference. `rsatoolbox`'s inferential
  machinery (subject and condition bootstrap, noise ceilings, model
  comparison, as in Schütt et al. 2023) is out of scope here and is not
  claimed to be contained in the bilinear core.
* The analytic law is admitted only under the declared model: equal
  partitions, one common **fixed** metric, and a separable plug-in error
  channel. `04-extension.R` provokes and records the refusals that mark that
  edge (`results/extension-refusals.csv`): a learned metric or a missing
  `target` gets a classed `effect_capability_refusal`, not a plausible
  number.
* Correlation distance is outside the fixed-linear subset and is not compared;
  see `vignette("correlation-distance-policy")`.
* This is a numerical-agreement exemplar. It records no timing and makes no
  cross-package speed claim.

## Files

| file | role |
| --- | --- |
| `00-common.R` | fixture constants, model RDMs, the pooled precision estimator |
| `01-fixture.R` | build + fit + export (R → CSV) |
| `02-rsatoolbox.py` | the external arm (CSV → rsatoolbox → CSV) |
| `03-compare.R` | the agreement table; non-zero exit if any row fails |
| `04-extension.R` | the strict extension, identities, and refusals |
| `requirements.txt` | `rsatoolbox==0.3.2`, `numpy==2.5.2` |
| `run-all.sh` | the four steps in order |
| `results/agreement.csv` | the ratchet read by `tests/testthat/test-rsatoolbox-parity.R` |
| `results/extension-table.csv` | the two-column comparison above |
| `results/extension-uncertainty.csv` | RDM SEs and the transported RSA SE |
| `results/extension-refusals.csv` | the three refusals, verbatim |

The cross-language exchange payload (`betas.csv`, `residuals.csv`,
`precision.csv`, `covariance.csv`, `fixture.rds`) is regenerated
deterministically by `01-fixture.R` and is not versioned; see
`results/.gitignore`.
