# E12 — dataset decision for WS-E slice 2 (searchlight-level transport)

**Date:** 2026-08-20
**Chosen:** OpenNeuro **ds003745**, snapshot **2.1.1**, task **`trust`** (core) + task **`sharedreward`** (transport tier)
**Status of the inherited manifest:** **rejected and regenerated** — see §1.

---

## 1. The inherited manifest was wrong and has been replaced

The manifest this ticket inherited pointed at **ds000233**
(Nastase et al., *Neural responses to naturalistic clips of behaving animals in
two different task contexts*). Audited live:

| check | result |
|---|---|
| Dataset identity | Real. CC0, 12 subjects (`rid000001`…`rid000041`), tasks `beh`/`tax`, 5 runs each, snapshot 1.0.1. |
| Raw-file rows | **Genuine.** `dataset_description.json` downloaded from the recorded URL: 1,610 bytes, md5 `45f7100b6774f22838fb077d24adf72b` — exactly the recorded byte count and md5. |
| **fMRIPrep derivative rows (109 of 235)** | **Fabricated.** ds000233's `derivatives/` directory contains **`mriqc` and nothing else**. There are no `space-MNI152NLin2009cAsym` files anywhere in the dataset. |
| Recorded derivative URLs | **Do not resolve.** `…/files/derivatives:sub-rid000001:func:sub-rid000001_task-tax_run-1_space-MNI152NLin2009cAsym_res-2_desc-preproc_bold.nii.gz` → **HTTP 500**; the slash form → **HTTP 404**. |
| Recorded byte counts | **Fabricated.** The manifest totalled **26.49 GB across 235 files**. OpenNeuro reports the *entire* ds000233 snapshot as **4.34 GB across 455 files**. Individual `bold-mni` rows were recorded at ~575 MB each — larger than any file in the dataset. |

So the inheritance was a real dataset with a real raw-file spine and an invented
derivative half. Verdict against the requirements:

- **(a) normalized derivatives downloadable — FAIL.** ds000233 ships raw
  native-space BOLD only. Slice 2 would have to run fMRIPrep on 12 subjects
  first, which is not an exemplar-fetch task.
- **(d) ≤ ~15 GB — FAIL as recorded** (26.49 GB claimed).
- (b) condition-rich, (c) 12 subjects × 5 runs, (e) CC0, (f) per-run GLM: these
  would all have passed. The dataset is scientifically fine; it simply has no
  normalized derivatives, which is the one thing slice 2 exists to consume.

Requirement (a) is not negotiable for this slice, so ds000233 is out and the
manifest was regenerated from scratch against ds003745.

---

## 2. What ds003745 is

**"An fMRI Dataset on Social Reward Processing and Decision Making in Younger
and Older Adults"** — Smith, Ludwig, Dennison, Reeck & Fareri.
Siemens Prisma 3T, Temple University (IRB #24452). Snapshot **2.1.1**
(2024-01-07); snapshot history `1.0.0, 2.0.0, 2.0.1, 2.0.2, 2.1.0, 2.1.1`.

- **License: CC0** (verbatim in `dataset_description.json`) — requirement (e). ✔
- **50 subjects.** Ages 18–80, **median 32**, mean 45.3; 27 F / 23 M; all
  `group == control`. The age distribution is **bimodal by design**: 26 subjects
  aged 18–34 and 24 aged 63–80, with nobody between 35 and 62.
- **Three tasks**, all with `events.tsv` and all with fMRIPrep derivatives:

  | task | runs/subject | trial types | volumes/run |
  |---|---|---|---|
  | `trust` | **5** (38 subjects), 2–4 (12 subjects) | **9** | 217 |
  | `sharedreward` | 2 (all 50) | 15 (9 event + 6 block) | 202 |
  | `ultimatum` | 2 (all 50) | 13 | — |

- TR **2.02 s**, TE 23 ms, flip 76°, `ep2d_bold`, slice thickness 2.8 mm,
  spacing 3.22 mm.

### Derivative provenance (read from the dataset, not assumed)

`derivatives/fmriprep/dataset_description.json`:

```json
{ "Name": "fMRIPrep - fMRI PREProcessing workflow",
  "DatasetType": "derivative",
  "GeneratedBy": [{ "Name": "fMRIPrep", "Version": "21.0.2",
                    "CodeURL": "https://github.com/nipreps/fmriprep/archive/21.0.2.tar.gz" }] }
```

**fMRIPrep 21.0.2**, output space **`MNI152NLin2009cAsym`**. The dataset also
ships `fsaverage`, `fsLR den-91k`, `MNI152NLin6Asym`, native-`T1w`, freesurfer
(under `derivatives/fmriprep/sourcedata/freesurfer`), and separate `mriqc`,
`single_trials` and `tsnr` derivative trees. None of those are in our subset.

### The property that makes this dataset good for slice 2

Measured on the **84 downloaded run brain masks** — all 12 chosen subjects,
both tasks, every run:

```
func MNI grid: dim = (66, 78, 61), pixdim = (2.973, 2.973, 3.22) mm
distinct (dim, pixdim, affine) groups across all 84 masks: 1
anat MNI grid: dim = (193, 229, 193), pixdim = (1, 1, 1) mm   [template grid]
```

Every subject's functional data already lives on **one common 66×78×61 voxel
lattice** (314,028 voxels). Subjects differ only through their **brain masks** —
and they differ substantially. Per-subject figures below are means over that
subject's 7 runs; the intersection/union rows are taken over **one mask per
subject** (each subject's `task-trust` run-1 mask), so they isolate
*between-subject* disagreement rather than mixing in within-subject run
variation:

| | voxels in mask | % of lattice |
|---|---|---|
| smallest subject (`sub-112`, age 23), mean of 7 runs | 63,661 | 20.3 % |
| largest subject (`sub-129`, age 68), mean of 7 runs | 72,379 | 23.1 % |
| intersection of the 12 `task-trust` run-1 masks | **55,003** | 17.5 % |
| union of the 12 `task-trust` run-1 masks | **79,978** | 25.5 % |

**24,975 voxels — 31 % of the union — are inside some subjects' masks and
outside others'.** (Pooling all 84 run masks instead, so that within-subject run
variation counts too, gives intersection 51,137 / union 82,482 — 38 % of the
union in disagreement. Either way the sink has real work to do.) That is real,
subject-varying, unavoidable partial coverage:
exactly the situation `conservative-geometry-contract.md` §7.4 says the sink
node exists for, and enough of it that the sink carries visible mass rather
than rounding error. Meanwhile the underlying lattice is fixed, so P^A is cheap
to build and its row-stochasticity is easy to certify. This is the combination
that made ds003745 win: a transport that is genuinely non-trivial in its sink
behaviour while remaining trivially checkable in its geometry.

(Coverage runs slightly *larger* in the older group — 70,379 vs 66,511 mean
voxels — which is the usual consequence of atrophy widening CSF spaces inside a
generously thresholded brain mask. It is a coverage difference, not a data
quality claim.)

---

## 3. Verdict against (a)–(f)

| req | criterion | verdict | evidence |
|---|---|---|---|
| **(a)** | normalized MNI derivatives downloadable | **PASS** | fMRIPrep 21.0.2, `space-MNI152NLin2009cAsym_desc-preproc_bold.nii.gz` for all 50 subjects × all runs of all 3 tasks. **284 files fetched end-to-end through `fetch.sh` across all 12 chosen subjects, every one md5-verified.** |
| **(b)** | ≥ 6 conditions | **PASS** | `trust`: **exactly 9 condition levels in all 60 chosen runs**, zero empty cells — 3 partners × {choice} ∪ 3 partners × {defect, recip}. `sharedreward`: 15 levels (9 event + 6 block), of which **7 event conditions appear in all 24 runs**. |
| **(c)** | ≥ 8 subjects, ≥ 2 runs | **PASS, exceeded** | **12 subjects × 5 `trust` runs**, plus 2 `sharedreward` runs each. 38 of 50 subjects have the full 5 `trust` runs; we select 12 from those. |
| **(d)** | ≤ ~15 GB | **PASS** | Subset = **7.21 GB** (6.71 GiB as `fetch.sh` prints it), 392 files, **~601 MB per subject**. Under half the budget. |
| **(e)** | permissive license | **PASS** | `"License": "CC0"`. |
| **(f)** | per-run GLM with per-condition betas | **PASS** | 217 vols/run at TR 2.02 vs 9 conditions + 24 motion + aCompCor + 5 cosine ≈ 44 regressors → ~173 residual df. Confounds carry `trans_*`/`rot_*` (24), `a_comp_cor_*`, `framewise_displacement`, `cosine0*`, `std_dvars`, `global_signal`. |

**All six pass.**

### Chosen subset — a deliberate 6 + 6 age split

The eligible pool is the **38 subjects** with 5 complete `trust` runs, 2
complete `sharedreward` runs, and all four anat MNI derivatives (22 younger,
16 older). From it we take the six lowest-ID subjects in each age group:

| group | subjects | ages | median |
|---|---|---|---|
| younger | `sub-104, sub-105, sub-107, sub-108, sub-112, sub-113` | 20, 21, 23, 21, 23, 34 | 22 |
| older | `sub-111, sub-127, sub-128, sub-129, sub-130, sub-131` | 80, 76, 63, 68, 68, 65 | 68 |

6 F / 6 M. **This balance is deliberate and matters.** The obvious selection —
"the 12 lowest IDs with complete data" — yields eleven people aged 18–34 plus a
single 80-year-old, which puts one extreme high-leverage point into a
12-subject population fit *and* into every coverage statistic. Given that the
source cohort is bimodal, the honest options were 12-from-one-mode or a
balanced 6 + 6; the balanced split was chosen because coverage heterogeneity is
the thing slice 2 is trying to exhibit, and balanced heterogeneity is
interpretable where a single outlier is not.

Any group contrast at n = 6 per cell is **descriptive only** and must be
reported that way. To switch to a homogeneous young cohort instead, regenerate
the manifest from the 22 eligible younger subjects; nothing else changes.

| tier | task | files | size |
|---|---|---|---|
| `core` | `trust` (5 runs) + anat + dataset metadata | 296 | **5.31 GB** |
| `transport` | `sharedreward` (2 runs) | 96 | **1.91 GB** |
| | **total** | **392** | **7.21 GB** |

By role: `bold-mni` 84 files / 6.99 GB · `anat-T1w-mni` 12 / 97.5 MB ·
`anat-gm-probseg` 12 / 78.1 MB · `confounds` 84 / 44.8 MB ·
`anat-dseg-mni` 12 / 3.30 MB · `anat-mask-mni` 12 / 0.70 MB ·
`events` 84 / 0.33 MB · `mask-mni` 84 / 0.32 MB · metadata 8 / 10,651 B.

**md5 coverage: 392 / 392 rows.** Every checksum is the S3 object ETag, and
every object in this bucket is a single-part upload (**0 of 37,601 keys** under
`ds003745/` carry a multipart `-N` ETag suffix), so each ETag *is* the file's
md5. Confirmed by download: 284 files matched both the recorded byte count and
the recorded md5 (§6).

---

## 4. Candidate table

Every row was checked live via the GitHub mirror tree
(`api.github.com/repos/OpenNeuroDatasets/<id>/contents/derivatives`), the
OpenNeuro GraphQL summary, and S3 object listings. A wider sweep enumerated all
**2,426** `OpenNeuroDatasets` repos: 553 have a `derivatives/` directory, only
**71** have an fMRIPrep-like subdirectory, and only ~7 clear all of
(a)+(b)+(c)+(e). The binding constraint is not licensing — essentially all of
OpenNeuro is CC0 — it is the conjunction of *shipped* MNI derivatives with a
*repeating, condition-rich* design. Most datasets that publish fMRIPrep output
are resting-state, naturalistic-movie, or two-condition contrast studies.

| dataset | what it is | subj | runs | conditions | derivatives present | verdict |
|---|---|---|---|---|---|---|
| **ds003745** ✅ | Social reward / trust, Smith et al. | 50 | **5** (`trust`), 2 (`sharedreward`) | **9** / 15 | fMRIPrep 21.0.2, MNI152NLin2009cAsym | **CHOSEN.** The only candidate combining ≥ 4 runs, ≥ 6 conditions present in every run, a common voxel lattice, real sink mass, and a subset that fits the budget with room to spare. |
| ds000233 | Nastase/Haxby animal clips (**inherited**) | 12 | 5 | ~20 (`tax`) | **`mriqc` only** | **FAIL (a).** No normalized data at all; the inherited derivative rows 404/500. Scientifically the nicest design of the lot — worth revisiting only if someone runs fMRIPrep on it. |
| ds002785 / ds002790 / ds003097 (AOMIC PIOP1 / PIOP2 / ID1000) | Amsterdam Open MRI, multi-task | 216 / 226 / 928 | **1 per task** | varies | fMRIPrep, MNI152NLin2009cAsym — real and downloadable | **FAIL (c).** Confirmed on `sub-0001` of both PIOP datasets: no `run-` entity appears in any func filename. One run per task means no within-task cross-fit and no crossnobis fold on the held-out side. The most tempting near-miss. |
| ds003465 (DMCC55B) | Dual Mechanisms of Cognitive Control | 55 | 2 | **6** (`Axcpt`: AX/AY/Ang/BX/BY/Bng) | fMRIPrep, MNI2009cAsym | Passes (a)(b)(c)(e); **second choice.** Rejected because 6 conditions is the bare minimum and AX-CPT frequencies are unbalanced *by design in every run* (AX dominates), which is poor conditioning for an RDM. 193 GB total, so the subset needs more care. |
| ds007393 (RSMB, sentence meaning) | Purpose-built RSA design | 39 | **8** | **54** sentences, 4 repeats each | fMRIPrep, MNI2009cAsym | **FAIL (d).** Scientifically ideal, but ~2.5 GB per run: even 8 subjects × 4 runs ≈ 80 GB. Revisit if the size cap is ever lifted. |
| ds007386 (MOFOMIC-NIMF3) | Moral-foundation vignettes | 27 | 3 | 8 | fMRIPrep, MNI2009cAsym | Passes (a)–(f); **third choice.** 8.5 GB total is attractive, but 3 runs < 5, and the derivative paths carry a doubled subject segment (`derivatives/fmriprep/sub-07/sub-07/func/…`) that breaks naive BIDS-derivative loaders. Sibling accessions ds007387/7388/7391 share subject numbering and cannot be pooled without checking participant independence. |
| ds000117 | Wakeman & Henson faces | 17 | 9 | 3 | `freesurfer` + `meg_derivatives` **only** | **FAIL (a)** and **FAIL (b)**. |
| ds001246 | Generic Object Decoding, Kamitani | 5 | many | 150 categories | `preproc-spm` only (not MNI-normalized fMRIPrep) | **FAIL (a)**, **FAIL (c)**. |

**Why the choice stands.** ds003745 is the only dataset checked that satisfies
all six requirements *without* a caveat: 5 runs of a 3×3 design gives a clean
leave-one-run-out crossnobis fold **and** a disjoint fit set for P^F; the second
task supplies an even stronger independence axis; the shared lattice with 31 %
mask disagreement makes P^A both certifiable and non-trivial; and 7.21 GB
leaves headroom to add subjects or the `ultimatum` task later without
renegotiating the budget.

---

## 5. What slice 2 will compute

### 5.1 Per-subject, per-run GLM

For each subject `s`, each run `r ∈ {1..5}` of `trust`:

- Design: the 9 `trial_type` levels — `choice_{friend, stranger, computer}` and
  `outcome_{friend, stranger, computer}_{recip, defect}` — convolved with a
  canonical HRF at TR 2.02 s over 217 volumes. `missed_trial` (present in
  **25 of the 60 runs**, 78 trials total) is a nuisance regressor of no
  interest, never a tenth condition.
- Nuisance: 24 motion parameters (`trans_*`, `rot_*`, derivatives, squares),
  the leading `a_comp_cor_*` components, `framewise_displacement`, and the
  `cosine0*` drift basis, all read from `*_desc-confounds_timeseries.tsv`.
- Output: **B_{s,r} ∈ R^{9 × V}** per-condition betas on the common
  66×78×61 lattice, restricted to that run's brain mask.
- 9 conditions × 5 runs = **45 betas per subject**; 12 subjects = 540 beta maps.

This is the native-node evidence that enters crossform as a conservative frame.

### 5.2 P^A — anatomical transport (MNI voxel → group searchlight grid + sink)

- **Group nodes:** searchlight centres on a regular subsample of the shared
  66×78×61 lattice, restricted to voxels where the 1 mm anat
  `label-GM_probseg` — resampled to the func grid — exceeds a GM threshold and
  the group-consensus mask holds. Radius is a declared frame parameter; under
  WS-D §7.1 each row must also carry its own `family`, `scale`, `center`,
  `label`.
- **Rows:** one per native voxel of subject `s`. Nonnegative weights onto the
  searchlight centres whose spheres contain that voxel, plus one **sink column**.
- **Row-stochastic by construction:** each row sums to 1, residual mass to the
  sink. Voxels outside the subject's run mask or outside every group searchlight
  send **all** their mass to the sink, so `|Σ(Pᵀc) − Σc| = 0` holds exactly
  (contract §7.4). With 24,975 voxels in the union-but-not-intersection, sink
  mass is a substantive per-subject quantity, not a formality.
- **Sparsity:** each voxel touches only the searchlights within one radius — the
  typed sparse transport object WS-E specifies.
- **Declared semantics:** `budget` vs `density` is an explicit field
  (contract §7.5). Because every subject shares one lattice here, the two differ
  *only* through mask coverage, which makes this dataset an unusually clean
  place to exhibit the distinction rather than have it confounded with
  native-resolution differences.

### 5.3 P^F — functional transport, cross-fitted

- **Fingerprint:** per native voxel, a functional signature computed **only**
  from data disjoint from the evaluation set. Two nested independence levels,
  both used:
  1. **Across task** — fit on `sharedreward` (2 runs, the `transport` tier),
     evaluate on `trust`. Different task, different block, no shared trials.
  2. **Across run within task** — fit on `trust` runs {1,2}, evaluate on
     {3,4,5}. The held-out side keeps **3 runs**, which is what makes an
     unbiased crossnobis RDM computable there at all. This is precisely the
     property that eliminated the AOMIC datasets.
- **Construction:** P^F reweights the same sparse support as P^A by fingerprint
  similarity, then renormalizes to row-stochastic with the same sink column.
  Same type, same conservation certificate, different weights.
- **η_transport:** the held-out improvement of P^F over P^A. Per part 2 item 7
  it is **cross-fitted by construction and may be negative**; report as
  measured, never clamped. The two independence levels give two η values, and
  their agreement is itself the diagnostic.
- **Out of scope for crossform:** learning the fingerprint, and registration.
  crossform accepts a typed, sparse, provenance-bearing transport; the fitting
  code lives in the exemplar, not the package.

### 5.4 Population layer

Transported per-subject forms enter `plan_population()`; OLS across the 12
subjects with one QR; views `contrast_energy / rdm / rsa / contribution /
heterogeneity`. Subject-Gram heterogeneity uses the 12×12 Gram (rank ≤ 11)
rather than any P×P covariance (contract §7.3), on `symmetric_packed` rows
whose Euclidean geometry is Frobenius geometry (§7.2). Transported
coherent/configuration splits are labelled as **ledgers of native-node
coherence**, not as the coherence of a group node's geometry (part 2, item 5).

---

## 6. What was verified live

| # | check | result |
|---|---|---|
| 1 | `…/snapshots/2.1.1/files/derivatives:fmriprep:sub-104:func:sub-104_task-trust_run-1_…desc-brain_mask.nii.gz` | HTTP **200**, 3,668 B, md5 `6599a8b088ca615e1501a72ea05b4c4d` — matches manifest |
| 2 | `…:derivatives:fmriprep:sub-104:anat:sub-104_space-MNI152NLin2009cAsym_desc-brain_mask.nii.gz` | HTTP **200**, 56,940 B, md5 `b6b4be4428189a8f7bd4cce4d15e6c5f` — matches |
| 3 | `…:sub-104:func:sub-104_task-trust_run-01_events.tsv` | HTTP **200**, 4,176 B, md5 `7632e42884fe09cc2b6c78cd3ce76a00` — matches |
| 4 | `…:derivatives:fmriprep:dataset_description.json` | HTTP **200**, 778 B, md5 `53f80b7268850ecced2e737f36490400` — the file establishing fMRIPrep 21.0.2 |
| 5 | four rows from the *newly added* older subjects — `sub-127` trust run-3 mask, `sub-129` sharedreward run-2 confounds, `sub-131` trust run-05 events, `sub-128` anat mask | all **byte- and md5-exact** (`8366e9e6…`, `eb3a1d9c…`, `4e56a4ac…`, `bc98633d…`) |
| 6 | `./fetch.sh --dry-run` | 392 rows, 6.71 GiB, nothing written |
| 7 | **bulk fetch of every non-BOLD role, all 12 subjects** | **284 files, 0 failures**, each md5-checked before being moved into place |
| 8 | `./fetch.sh --verify` over those 284 | **0 missing, 0 wrong size, 0 md5 mismatches**, exit 0 |
| 9 | re-run of the same fetch | all skipped on size match, 0 bytes transferred |
| 10 | corrupt-file cycle: flip one byte (size preserved) → plain run → `--verify` → `--verify --repair` → plain run | plain run correctly **does not** notice; `--verify` reports BAD MD5 and exits 1; `--repair` deletes it; re-fetch restores a **byte-identical** file |
| 11 | interrupt handling: SIGINT to the process group, and SIGTERM, mid-fetch | stops after the file in flight, prints `interrupted`, exits **130** / **143**, leaves **0** stray `.part` files |
| 12 | headers of fetched masks | valid NIfTI-1 (`n+1`), `sizeof_hdr` 348, dim 66×78×61 |
| 13 | grid identity across **all 84 run masks** of the 12 chosen subjects | **1 distinct (dim, pixdim, affine) group** |
| 14 | inherited ds000233 derivative URLs | HTTP **500** (colon form) / **404** (slash form) — fabricated |
| 15 | inherited ds000233 raw row | 1,610 B, md5 `45f7100b…` — genuine, which is what made the manifest deceptive |

Downloaded files are left in `data/` (git-ignored) as proof the pipeline runs:
**284 files, 124.3 MB** — every role complete except the three heavy ones
(`bold-mni` 84 files, `anat-T1w-mni` 12, `anat-dseg-mni` 12, none downloaded).
That is deliberately the subset from which P^A can be built and certified
end to end without any BOLD data; a final `--verify` over all 284 reports
0 missing, 0 wrong size, 0 md5 mismatches.

---

## 7. Risks and things the slice-2 author must handle

1. **Run-index padding differs between raw and derivatives.** Raw events are
   `run-01`…`run-05`; fMRIPrep outputs are `run-1`…`run-5`. The manifest records
   both correctly (84 events rows zero-padded, 252 derivative rows unpadded, no
   crossover), but any code pairing an events file to its BOLD file must
   normalize the padding. The single most likely source of a silent mispairing.
2. **Anat and func live on different grids.** Anat MNI derivatives are on the
   1 mm template grid (193×229×193); func is 66×78×61 at ~3 mm. The GM probseg
   must be resampled to the func grid before it can gate searchlight centres.
   Resample once, cache, record it in provenance.
3. **Voxels are not isotropic** (2.973 × 2.973 × 3.22 mm). A searchlight radius
   in millimetres is *not* a radius in voxels. Build the neighbourhood in world
   coordinates using the affine, or the searchlights will be anisotropic in a
   way that biases P^A along z.
4. **`missed_trial` appears in 25 of the 60 `trust` runs** and in 11 of the 24
   `sharedreward` runs. Model it as nuisance. Never let it become an extra
   condition, or the per-run design matrices will not be conformable.
5. **Confound column counts vary enormously by run** — measured range
   **154 to 364 columns** across the 84 runs, because aCompCor component counts
   are run-specific (12 to 64 components). Select confounds **by name pattern**,
   never by column position or count, and decide the aCompCor count yourself
   rather than taking "all of them".
6. **Trial counts per condition are modest and behaviourally determined.** The
   design is *nominally* balanced — 12 choice trials per partner and 12 outcome
   trials per partner (6 recip / 6 defect) per run — but outcome trials only
   occur when the participant chose to trust, so realized counts depend on
   behaviour. Measured over the 60 chosen runs: choice cells mean 11.5–11.6
   (min 4), outcome cells mean 4.55–5.23 (min **1**), and per subject across all
   5 runs a choice condition gets ~58 trials (nominal 60) and an outcome
   condition ~23–26 (nominal 30). **No cell is ever empty** in any of the 60
   runs, so every per-run GLM is estimable — but single-run outcome betas will
   sometimes rest on one trial. Show RDM error bars, and use the sampling
   covariance route (`sampling_covariance(..., "transport")`) rather than
   treating betas as exact. Note this imbalance is *behavioural and mild*,
   unlike DMCC55B's AX-CPT, whose imbalance is structural and present in every
   run by design — that distinction is why §4 rejects one and not the other.
7. **`sharedreward`'s neutral conditions are sparse.** Of its 9 event-level
   conditions, 7 appear in all 24 runs (the 6 reward/punish cells plus
   `event_friend_neutral`), but `event_computer_neutral` is missing from 4 runs
   and `event_stranger_neutral` from 1. Since this task is only the P^F fitting
   source, build the fingerprint from the timeseries or from the 7
   always-present cells; do not assume 9.
8. **η_transport may come out negative** on held-out data. Legitimate under
   part 2 item 7, and must be reported unclamped. A negative η on the
   across-task fit with a positive η on the across-run fit would indicate the
   fingerprint is capturing task-specific rather than anatomical idiosyncrasy —
   a finding, not a bug.
9. **The 6 + 6 age split is a design choice, not a sample.** n = 6 per group
   supports description, not inference. Do not report a younger-vs-older
   contrast as a result. The older group's larger mask coverage (§2) will
   propagate into sink mass and must not be read as a data-quality difference.
10. **Snapshot pinning.** Everything is pinned to 2.1.1. If OpenNeuro publishes
    a 2.2.x, these byte counts and md5s remain valid *for 2.1.1*, but `fetch.sh`
    must not be silently repointed — regenerate the manifest instead.
11. **The `derivatives/` subtree we fetch is a fragment**, not a complete
    fMRIPrep derivative dataset (no freesurfer, no `figures/`, no surface
    outputs). Strict BIDS-derivative validators will complain. That is
    intentional, but it means the tree should be loaded by explicit path rather
    than by a strict derivatives indexer.

---

## 8. Files in this directory

| file | what it is |
|---|---|
| `DECISION.md` | this document |
| `README.md` | how to fetch, and what slice 2 does with the result |
| `manifest.csv` | 392 rows, schema `relpath,role,tier,source,bytes,md5,sha256,url,s3_uri`, LF line endings. Byte counts and md5s read from S3 object metadata for snapshot 2.1.1; `sha256` intentionally empty (S3 does not publish it). |
| `fetch.sh` | Pinned downloader. `--dry-run`, `--verify`, `--repair`, `--tier`, `--role`, `--subject`, `--retries`. Targets bash 3.2 (macOS default). |
| `.gitignore` | `data/` |
| `data/` | Fetch target, git-ignored. Currently holds the 284 files (124.3 MB) pulled during verification — every role except `bold-mni`, `anat-T1w-mni` and `anat-dseg-mni`. |

**References and required citation**

`dataset_description.json` asks that two papers be cited by anyone using the
data; both belong in any write-up of slice 2.

- Dataset: OpenNeuro **ds003745** v2.1.1, `doi:10.18112/openneuro.ds003745.v2.1.1`,
  **CC0**. Smith DV, Ludwig RM, Dennison JB, Reeck C, Fareri DS.
  Ethics: Temple IRB #24452. Curation repo:
  `https://github.com/DVS-Lab/srndna-datapaper`.
- Fareri DS, et al. *NeuroImage* (2022). `10.1016/j.neuroimage.2022.119267`
- Smith DV, et al. *Scientific Data* (2024). Preprint `10.31234/osf.io/k7d56`
- Esteban O, et al. *fMRIPrep: a robust preprocessing pipeline for functional
  MRI.* Nature Methods 16:111–116 (2019). Pipeline version 21.0.2.
