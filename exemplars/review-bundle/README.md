# crossform external review bundle

This directory builds a self-contained evidence bundle for a reviewer outside
the project. The tarball it produces is also shipped with this file as its
top-level `README.md`, so what you are reading is either the build recipe (in
the source tree) or the front matter of the bundle itself (in the tarball).

## Who this is for

A second lab, a statistical reviewer, or a methods referee who wants to check
crossform's numerical claims without cloning the repository, installing the
package, fetching 300 MB of fMRI data, or trusting a narrative.

It contains the contract documents that say what the package promises, the
oracles that check those promises against independent implementations, the
recorded outputs of every exemplar analysis, the shipped certification
artifacts, and a manifest that hashes all of it and names the source commit it
was cut from. It contains no fetched data and no compiled code.

Three of the six oracles run in **about a second each on a bare R installation
with no packages at all** (measured from an extracted bundle: 1.25 s, 0.46 s,
0.56 s including R startup). That is the cheapest possible entry point: it
exercises the population-geometry algebra with the package deliberately not
loaded.

```sh
tar -xzf crossform-review-bundle-<date>.tar.gz
cd crossform-review-bundle-<date>
shasum -a 256 -c SHA256SUMS          # or: Rscript build-bundle.R --verify=...
Rscript design/oracles/population-geometry-split.R
```

## Building it (from the source tree)

```sh
Rscript exemplars/review-bundle/build-bundle.R
Rscript exemplars/review-bundle/build-bundle.R --verify=dist/crossform-review-bundle-<date>.tar.gz
```

The bundle lands in `dist/` (git-ignored; the script writes the `.gitignore`).
It reads the working tree and writes nothing else. `--help` lists the options;
the ones that matter are `--force` to overwrite and `--max-file-mb=N` to change
the per-file inclusion cap.

## What is in it

| path | what it is |
| --- | --- |
| `MANIFEST.txt` | every payload file with its SHA-256 and byte count, plus the source commit, branch, working-tree state, and `R/` digest the bundle was cut from |
| `SHA256SUMS` | the same hashes as plain `shasum -a 256 -c` input |
| `build-bundle.R` | the script that cut this bundle, shipped so the bundle describes its own assembly; `--verify` re-checks the manifest and needs no checkout |
| `RUNNING-TESTS.md` | one page: re-running the fast test suite from a source checkout, and how to read the certification skips |
| `design/*.md` | the six contract documents |
| `design/oracles/*.R` | six standalone oracle scripts |
| `exemplars/rsatoolbox-parity/` | scripts, pinned `requirements.txt`, and the recorded agreement and extension CSVs |
| `exemplars/er-rsa/` | scripts and the recorded recovery CSVs |
| `exemplars/haxby2001/` | scripts, README, and the committed smoke receipts (`.rds`, one `.png`) |
| `exemplars/population-slice2/` | the dataset decision record, the pinned manifest, and the fetch script |
| `certification/artifacts/` | the shipped gate records and per-gate summary CSVs |
| `certification/BINDING.txt` | per artifact: which `R/` digest it binds to, and whether that was current when the bundle was cut |
| `certification/` | the binding helper, the digest function, the promotion script, and the re-certification instructions |

Total is a few megabytes. Nothing under `exemplars/*/data/` is included at any
size: the Haxby subject-1 tarball (~300 MB, SHA-256 recorded in
`exemplars/haxby2001/00-download.R`) and the OpenNeuro ds003745 subset (7.21 GB
across 392 files, every one md5-pinned in
`exemplars/population-slice2/manifest.csv`) are fetched, never versioned. Both
have checksummed, idempotent fetch scripts in the bundle.

---

# The claims this bundle evidences

Seven claims, each with the file that carries it. Numbers are quoted from the
recorded artifacts, not recomputed here.

### 1. Fixed-metric crossnobis and linear RSA agree with an independent Python implementation to floating-point reordering

`exemplars/rsatoolbox-parity/`, `results/agreement.csv`.

On one shared synthetic fixture (6 conditions, 4 runs, 40 voxels, non-spherical
residual covariance with condition number 182), crossform's fixed-metric
crossnobis RDM agrees with Python `rsatoolbox` 0.3.2 to **3.77e-15** maximum
absolute difference across all 60 RDM entries, and its linear RSA coefficients
agree with `numpy.linalg.lstsq` on the vectorised RDMs to **8.60e-16** across 20
coefficients. The declared tolerance is 1e-10. `03-compare.R` exits non-zero if
any row of the agreement table exceeds its tolerance.

Four conventions had to be matched — the noise precision, the cross-validation
folding, the division by channel count, and the RSA objective — and all four
were resolved exactly rather than by loosening a tolerance. Two of them are
cross-checked by third oracles inside `02-rsatoolbox.py`:
`prec_from_residuals(method="full", dof=168)` reproduces the R pooled precision
to **2.31e-14**, and leave-one-fold-out folding matches the explicit uniform
C(4,2) pairing to **3.33e-16**.

`exemplars/rsatoolbox-parity/README.md` documents each convention and how it
was resolved.

### 2. From the same fit, crossform returns quantities the RDM cannot express, under a declared and refusable model

`exemplars/rsatoolbox-parity/results/extension-table.csv`,
`extension.csv`, `extension-uncertainty.csv`, `extension-refusals.csv`.

Nothing is refitted between claim 1 and claim 2: the same
`plan_geometry(metric = noise_precision(...))` object yields the signed contrast
energy (whole brain, face − house: **+0.230130**), an exact coherent /
configuration / total partition of it (0.008804 + 0.198452 = 0.207256), the same
partition of every RDM entry, and analytic sampling covariance for both the RDM
and a linear RSA coefficient.

Three identities are asserted in `04-extension.R` rather than assumed:
`crossnobis(plan, w)` equals `contrast_energy(plan, w)$total` **exactly**
(difference 0); `coherent + configuration − total` is **0** for the contrast and
≤ 2.8e-17 for every RDM entry; and each RDM entry re-derived as
`crossnobis(plan, e_i − e_j)` reproduces `rdm(plan)` to ≤ 2.2e-16.

Two entries in that table are boundaries, not results. Region `roiC` has a
negative signed contrast (−0.1646) alongside a positive energy (0.02337) — the
squared distance cannot say which way the effect goes. Region `roiB`'s coherence
fraction is **withheld rather than printed**, because its cross-generalized
components do not form a nonnegative partition there and the ratio would be an
artefact; `$coherence_fraction_valid` is `FALSE`. Three provoked capability
refusals (a learned metric, a missing `target`) return classed
`effect_capability_refusal` conditions instead of plausible numbers, recorded
verbatim in `extension-refusals.csv`.

### 3. Conservation, composition, and transport hold — or fail — exactly where the contracts say

`design/oracles/`, six scripts. All six exit 0 and print `DONE`. **Read the
numbers, not the exit status**: only `conservative-multiscale-ledger.R` gates
(`stopifnot`, plus a printed `PASS (<= 1e-12)`); the other five print measured
errors without asserting them, so a regression would print a bad number and
still exit 0.

*Conservation and composition* (`conservative-metric-composition.R`,
`conservative-multiscale-ledger.R`). Under a fixed neural metric, `sum_x G_x =
G_Omega` holds to **1.11e-16** for an identity metric and **2.63e-16** for a
diagonal one, and **fails at 2.38e-01 relative — a signed trace error of
+21.2 %** for a dense metric under the native composition
`D(sqrt(w)) Q D(sqrt(w))`. The failure is not noise: the oracle predicts it in
closed form as `sum_x G_x = B_L (S ∘ Q) B_R^T` and the prediction matches to
**8.88e-16**. The whitened composition `Q^(1/2) D(w) Q^(1/2)` restores
conservation to **1.53e-15**. Over a 12-draw sweep the native signed trace error
ranges **[−34.5 %, +188.4 %]**. A separate result with contract consequences:
symmetric-PSD and Cholesky roots both conserve, yet differ per-node by **15.7 %**
of the largest node value — so the choice of root is part of the estimand, not
an implementation detail. The multiscale ledger oracle is the one hard gate:
scale-wise conservation `sum_{x in scale s} G_x = alpha_s * G_Omega` passes at
worst absolute **8.88e-16** against its 1e-12 bound, while the *coherent*
component deliberately does not conserve (deviations 7.410 / 2.703 / 1.271 by
scale), which is why the coherent share is the only scale-resolved quantity.

*Transport* (`conservative-transport-readiness.R`,
`population-transport-contract.R`, `population-transport-diagnostics.R`,
`population-geometry-split.R`). The `sqrt(2)` symmetric packing is load-bearing:
the correct codec is Frobenius-consistent to **3.11e-15** while the naive
packing is wrong at **O(1)** (2.96). The N×N subject Gram reproduces the D×D
covariance spectrum to **1.07e-14** at rank N−1, so the large matrix is never
formed. Row-stochastic transport with a sink preserves the signed budget exactly
(`|group + sink − native|` = **8.88e-16**, worst **7.11e-15** over 500 random
transports); *without* a sink column the same transport silently loses **48.1 %**
of the budget. Query, transport, and fit commute under a subject-constant weight
operator (three evaluation orders agree to **8.88e-16**) and demonstrably do not
under per-node or per-coordinate weights (**25.2 %** and **25.5 %** of the
largest coefficient). The plug-in consensus/heterogeneity split is biased by
**+62.7 %** over 2000 Monte Carlo replications (`tr(Q^H)` bias +6.787 against a
predicted +6.720); the cross-fitted split has bias **+0.061** (MC se 0.082) but
is indefinite in **100 %** of replications, which is why PSD projection is a
recorded decision rather than a cleanup. The transport diagnostic
`eta_transport` is **+0.050286** honestly and **+0.158552** circularly — a
**3.15×** inflation — and under a null transport it is negative in **89.5 %** of
200 draws, so it needs its coverage diagnostics to be interpretable at all.

The three `population-*.R` oracles are **pure base R** and state so in their
headers; they load neither crossform nor any other package. The three
`conservative-*.R` oracles call `pkgload::load_all()` and exercise the real
package internals.

`design/conservative-geometry-contract.md` and
`design/population-form-contract.md` are the contracts these oracles check.

### 4. The rectangular encoding-retrieval machinery recovers a known ground truth on a designed simulation

`exemplars/er-rsa/`, `results/planted-vs-estimated.csv`,
`results/recovery-verdicts.csv`, `results/route-check.csv`.

36 encoding effects × 30 retrieval probes over 3 study-test cycles, 12 simulated
subjects, with unequal axes (12 items encoded but never retrieved, 6 lures never
encoded). Every readout has an exact planted value computed in closed form from
the noiseless per-run patterns — no Monte Carlo truth. **30 of 30 cross-run
readouts land inside their 95 % across-subject interval; the largest bias is
0.0093, or 1.2 standard errors.** Query-first `evaluate_geometry()` and
materialize-then-project `query_geometry()` agree to **6.66e-16** over all ten
queries, and the materialized rectangular form satisfies `total = coherent +
configuration` to **2.22e-16**.

Two results in that exemplar are about estimand choice rather than estimator
quality, and both are recovered to their planted values, which is what makes
them properties of the estimand:

- With an unrestricted control set, regionB — which was planted with **no**
  item-specific structure at all — reports a reinstatement effect of **1.03**
  (t = 104). Restricting the eligible set to same-category pairs drops it to
  **−0.01**. Both match their planted values (1.0228 and −0.0024).
- Under same-run rather than cross-cycle pairing, regionC — planted with **no
  task structure whatsoever** — reports item-specific reinstatement of **0.256**
  (t = 67), and regionA inflates by **47 %**. crossform makes these two
  different plan objects with different declared estimates
  (`cross_generalized` vs `self_product_biased`) rather than one number with a
  footnote.

Three boundaries were hit and recorded in `results/refusals.csv`, the live one
being that rectangular plans refuse a fixed noise metric — so there is no
crossnobis-style noise normalization and no analytic sampling covariance for a
rectangular readout, and the uncertainty above is across-subject spread over 12
noise realizations.

**This is a simulation.** See "What is not claimed" below.

### 5. Certification gates are recorded results bound to a source digest — and the shipped records are currently STALE

`certification/artifacts/`, `certification/BINDING.txt`,
`certification/helper-certification.R`, `design/certification-report.md`.

A recorded benchmark artifact is evidence for exactly one source tree. Every
runner stamps a `provenance` block carrying `source_digest`: the SHA-256 over
the sorted per-file SHA-256 of `R/*.R` (`certification/provenance.R`).
`certification/helper-certification.R` refuses to read a recorded gate as a
boolean unless that digest still equals the current one. Any edit to any file
under `R/` — including a comment — turns the certification tests into loud,
greppable `CERTIFICATION STALE` skips rather than silent passes. Under
`R CMD check` the sources are not shipped, so only the weaker package-version
tier is available, and the helper prints `CERTIFICATION TIER package_version`
rather than letting that pass be read as the stronger one.

**Honest status.** The shipped artifacts were recorded on 2026-08-17 against
`R/` digest `67f21604d0f0`. The package has been under active revision since,
so as of this bundle they are **stale against the current `R/` digest**, and the
end-of-program re-certification has not been run. `certification/BINDING.txt`,
generated at build time, records the exact comparison for your copy — which
digest each artifact binds to, which digest the tree was at, and the resulting
`BOUND` / `STALE` / `UNBOUND` state per artifact. Treat the recorded gate
numbers as evidence about the source state named in `BINDING.txt`, not about
whatever `R/` you are holding. `certification/RECERTIFY.md` is the ~12-minute
sequence that clears it.

Two artifacts are unbound by construction and skip with instructions even on a
matching tree: `shard-admission.rds` (whose recorded verdict is itself **"not
admitted"**) and `learned-metric-policy-validation.rds`.

The recorded numbers, from `design/certification-report.md`: the brain-scale
gate on the 52,416-feature fixture completed in **59.982 s** at an incremental
peak RSS of **903,479,296 bytes**, inside its 30-minute and 4 GiB limits. The
statistical validation — 500 paired replications × 2 regimes, predeclared
equivalence margin ±0.005 — found **3 of 4 arms equivalent and one not**: the
unmodelled-AR(1)-with-signal arm has a mean difference of **−0.0094**, 95 % MC
interval [−0.0114, −0.0073], outside the margin. That failure is reported, not
buried; training-only residual estimation remains the conservative default
because of it.

### 6. On real public data, three implementations of one estimand agree to machine precision

`exemplars/haxby2001/`, `results/smoke-comparison.rds`,
`results/smoke-report.md`, `results/coherent-configuration.rds`.

Haxby et al. (2001) subject 1 from the PyMVPA distribution — 12 runs, 8
categories, VT mask, 577 searchlight centers at radius 11.25 mm. Both arms
consume **identical per-run condition-mean matrices** from one shared prep
script, so no preprocessing difference separates them. Searchlight geometry is
verified rather than assumed: `neuroim2::searchlight_indices` and
`neuroim2::searchlight` give **bit-identical sphere membership at all 577
centers**.

Matched estimand, crossform `rdm()` against rMVPA's own crossnobis estimator:
maximum absolute difference **8.88e-16**, 0 of 577 centers above the 1e-8 gate.
Against an independent reference loop: **1.33e-15**. A refit through
`lm_relation_fit()` reproduces the same RDM to **4.44e-16**. The shared Spearman
score agrees **exactly** (max abs diff 0).

The *native* estimands of the two packages correlate at only **r ≈ 0.54**, and
that gap is fully accounted for by five **named** semantic differences rather
than absorbed into a tolerance. The exemplar's own success condition is that any
disagreement be traced to a named difference: "a silent numerical discrepancy is
a failure of the exemplar, not a footnote."

The coherent/configuration decomposition runs on this data with an exact
partition (`max|total − (coherent + configuration)|` = **5.55e-17** for face −
house across 577 searchlights) and is informative rather than degenerate: the
coherent share has median **0.529**, IQR [0.278, 0.773].

### 7. Population transport has a contract, a diagnostics slice, and a pinned dataset — but no result yet

`exemplars/population-slice2/DECISION.md`, `manifest.csv`, `fetch.sh`,
`design/population-form-contract.md`.

The dataset decision is included because it is itself a piece of evidence about
process. The inherited manifest was **rejected and regenerated**: it pointed at
ds000233, and **109 of its 235 fMRIPrep derivative rows were fabricated** — that
dataset's `derivatives/` contains only `mriqc`, the recorded URLs return HTTP
500/404, and the manifest claimed 26.49 GB across 235 files where OpenNeuro
reports 4.34 GB across 455. The replacement is OpenNeuro **ds003745** snapshot
2.1.1 (CC0, fMRIPrep 21.0.2), a deliberate 6 + 6 younger/older split, **392
files / 7.21 GB**, every file md5-pinned from S3 ETags (**392/392 coverage**).

The property that justified it: all 84 run masks sit on **one** voxel grid, so
the anatomical transport operator P^A is cheap to build and its row-stochasticity
trivially certifiable — while **31 %** of the union mask is covered for some
subjects and not others, so the sink still carries real subject-varying mass.

`fetch.sh` is checksum-gated and idempotent, promotes only after both size and
md5 match, and handles interruption without leaving `.part` files. The
recommended starter subset is **118.53 MiB** — enough to build and certify P^A
with no BOLD data at all.

**No population-level result is claimed.** This is a pinned dataset, a fetch
tool, and a contract. See below.

---

# How to run each layer

| layer | command | needs | time |
| --- | --- | --- | --- |
| verify the bundle | `shasum -a 256 -c SHA256SUMS` | nothing | < 1 s |
| population oracles (3) | `Rscript design/oracles/population-geometry-split.R` etc. | **bare R, no packages** | 1.25 s, 0.46 s, 0.56 s |
| conservative oracles (3) | `Rscript design/oracles/conservative-metric-composition.R` etc. | source checkout + `pkgload` | ~3.2 s each (mostly `load_all`) |
| fast test suite | `Rscript -e 'devtools::test()'` | source checkout, compiler | minutes — see `RUNNING-TESTS.md` |
| scale gates | `CROSSFORM_RUN_SCALE_TESTS=true Rscript -e 'devtools::test()'` | idle machine, ~4 GiB | opt-in; timing-sensitive |
| ER-RSA | `Rscript 01-simulate.R; 02-analyze.R; 03-recover.R` | package; writes ~8 MB to `data/` | ~4 s + ~15 s + ~1 s |
| rsatoolbox parity | `RSA_PYTHON=rsaenv/bin/python ./run-all.sh` | package **+ pinned Python** | ~1 min after setup |
| Haxby | `00-download.R` … `06-…R` | package, rMVPA, **~300 MB fetch** | ~2 min prep, ~20 s parity |
| re-certification | `certification/RECERTIFY.md` | idle machine | ~12 min |
| population slice 2 | `./fetch.sh --tier core` | network, up to 7.21 GB | dataset-dependent |

**The oracles are the cheapest real check** and the only layer that needs
nothing installed. The three `population-*.R` scripts run on a bare R with no
packages; the three `conservative-*.R` scripts need a crossform source checkout
because they exercise the package's own internals, which is the point.

## Reproducing the pinned Python environment (parity arm only)

Any Python 3.12 with the pinned wheels will do:

```sh
cd exemplars/rsatoolbox-parity
uv venv -p 3.12 rsaenv
uv pip install --python rsaenv/bin/python -r requirements.txt

Rscript 01-fixture.R            # build + fit + export CSVs
rsaenv/bin/python 02-rsatoolbox.py
Rscript 03-compare.R            # non-zero exit if any agreement row fails
Rscript 04-extension.R
```

There is no `reticulate` anywhere. The only channel between the two languages
is the CSV directory `results/`. The recorded run is CPython 3.12.11,
`rsatoolbox` 0.3.2, `numpy` 2.5.2, `scipy` 1.18.0, R 4.5, darwin/arm64,
2026-08-17.

Note that `results/` in this bundle already holds the recorded outputs. Re-running
overwrites them; the cross-language exchange payload (`betas.csv`,
`residuals.csv`, `precision.csv`, `covariance.csv`) is regenerated
deterministically and was never versioned, so it is absent here by design.

---

# What is NOT claimed

Stated plainly, because a bundle of agreeing numbers invites over-reading.

1. **No group-level or population inference is calibrated.** There is no LD-t,
   no bootstrap, no permutation test, no confidence-interval coverage study, no
   multiple-comparison correction, and no calibration of the analytic standard
   errors against a resampling reference. The analytic sampling covariance in
   claim 2 is an exact transport of a declared model, not a validated coverage
   claim. Claim 7's transport machinery has a contract and diagnostics; it has
   no result.

2. **The ER-RSA evidence is simulation-only.** A public encoding-retrieval fMRI
   dataset with per-trial betas could not be fetched and verified, so the
   demonstration is a designed simulation whose truth is known exactly. It is
   evidence that the rectangular machinery estimates what it claims to estimate
   on a design where the answer is known. It is **not** evidence about human
   memory, about any brain region, or about behaviour under real fMRI noise —
   the simulated noise is spatially correlated but temporally independent, so no
   whitener is supplied and none is claimed.

3. **The Haxby evidence is implementation parity, not scientific validation.**
   One subject, one hemisphere-scale mask, a block design, per-run condition
   means rather than GLM betas, an identity metric, at the smoke tier (577 VT
   centers; the whole-brain tier was not run). There is no ground truth, no
   external criterion, and **no inference of any kind** — no permutation test,
   no group model, no p-values, no correction. It shows that three
   implementations of one estimand agree; it says nothing about faces, houses,
   or animacy. The replacement-map row it earns is the cross-validated-distance
   estimand, not correlation-distance searchlight RSA, which **neither** package
   computes the way the exemplar's own estimand section describes.

4. **Every recorded run is macOS / Apple silicon.** R 4.5.x on
   `aarch64-apple-darwin20` with Accelerate BLAS. There is no Windows, Linux,
   Intel macOS, CRAN, or R-universe court, and no hosted CI. Numbers at the
   1e-15 level are BLAS-dependent; a different LAPACK will move the last digits.
   Each artifact records its own BLAS and LAPACK in its provenance block.

5. **The certification artifacts are stale.** See claim 5 and
   `certification/BINDING.txt`. They are recorded results for a named earlier
   source state. Re-certification against the final tree is pending.

6. **One certified arm did not pass.** The unmodelled-AR(1)-with-signal
   equivalence arm is *not* equivalent at the predeclared ±0.005 margin. It is
   in the report and in this list rather than only in the report.

7. **`shard-admission` records "not admitted".** That gate's recorded verdict is
   a refusal, and it is shipped as such.

8. **No speed claim is made anywhere.** The timings in these exemplars are plain
   elapsed times from single runs, without memory instrumentation, and in
   several cases the two arms compute different estimands. The only
   instrumented, budgeted measurement is the brain-scale certification gate, and
   it is a claim about crossform against its own declared budget — not against
   any other package.

9. **This bundle is not a release.** The package version is `0.0.0.9000`. The
   certification report declines to call the workspace a release commit.

## Pointers, if you want to go further

- `design/architecture.md` — how the pieces fit.
- `design/api-tiers.md` — which exports are load-bearing and why; the ER-RSA
  exemplar exists partly to justify five pair-query exports against this
  document.
- `design/crossform-execution-design.md` — the execution model behind the
  query-first / materialize-then-project equivalence checked in claim 4.
- `design/certification-report.md` — the full receipt, including its own errata
  and its recorded process defects.
- `certification/benchmarks-README.md` — what each gate asserts.
