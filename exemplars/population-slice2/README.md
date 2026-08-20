# population-slice2 — searchlight-level transport exemplar

This directory pins the public dataset that **WS-E slice 2** runs on, and gives
you a reproducible way to fetch it. It contains no data and no analysis code —
only a manifest, a downloader, and the record of why this dataset.

Slice 1 does population form at *region* level on the six Haxby subjects, where
transport is trivial. Slice 2 is the hard half: **searchlight-level transport**
between subjects in a common normalized space, with sink mass for partial
coverage, and a functionally-informed transport cross-fitted from independent
data.

- **Dataset:** OpenNeuro [`ds003745`](https://openneuro.org/datasets/ds003745),
  snapshot **2.1.1** — *An fMRI Dataset on Social Reward Processing and Decision
  Making in Younger and Older Adults* (Smith, Ludwig, Dennison, Reeck & Fareri).
  **CC0.**
- **Derivatives:** fMRIPrep **21.0.2**, space `MNI152NLin2009cAsym`.
- **Subset:** 12 subjects (6 younger, 6 older) × 5 runs of task `trust`
  (9 conditions) plus 2 runs of task `sharedreward` — **392 files, 7.21 GB**.

Read [`DECISION.md`](DECISION.md) for how this dataset was chosen, the
alternatives that were rejected, everything verified live, and the risk list.
**Skim §7 before writing analysis code** — the run-index padding mismatch
between raw and derivative filenames (item 1), the anisotropic voxels (item 3),
and the 154–364 column swing in the confounds files (item 5) will each bite
silently.

## Fetching

```sh
./fetch.sh --dry-run      # list what would be fetched, print the total
./fetch.sh                # fetch all 392 files (7.21 GB) into ./data
./fetch.sh --verify       # re-check sizes and md5s of what is on disk
./fetch.sh --verify --repair   # ...and delete whatever fails, so a plain
                               #    run afterwards re-fetches it
```

Everything lands under `data/`, which is git-ignored, with BIDS relpaths
preserved — so `data/` is a partial BIDS tree with a
`data/derivatives/fmriprep/` subtree beside it.

Re-running is cheap: a file already on disk at the manifest's byte count is
skipped without a network request. Downloads go to a `.part` file and are
promoted only after both the size **and** the md5 match, so an interrupted run
never leaves a truncated NIfTI in place. Note the corollary: a plain run trusts
the byte count, so a file that is corrupt *at the right size* is only caught by
`--verify`.

Useful subsets while developing:

```sh
./fetch.sh --tier core                      # task-trust + anat + metadata (4.94 GiB)
./fetch.sh --subject sub-104,sub-127        # one younger, one older       (1.11 GiB)
./fetch.sh --role events,confounds,mask-mni # design + nuisance, no BOLD   (43.37 MiB)
```

`--role events,confounds,mask-mni,anat-mask-mni,anat-gm-probseg` is the one to
start with: 118.53 MiB, and it is everything needed to **build and certify
P^A** — the searchlight support, the GM gating, and the per-subject mask
coverage that determines sink mass — before a single BOLD volume is downloaded.
Fitting P^F and computing η_transport does need the timeseries, so that step
waits for the full fetch.

(`fetch.sh` reports binary units; `DECISION.md` quotes decimal GB. The 392-file
total is 7.21 GB there and 6.71 GiB here.)

## The manifest

`manifest.csv`, one row per file, LF line endings:

| column | meaning |
|---|---|
| `relpath` | path under `data/`, identical to the path inside the snapshot |
| `role` | `bold-mni`, `mask-mni`, `confounds`, `events`, `anat-T1w-mni`, `anat-mask-mni`, `anat-gm-probseg`, `anat-dseg-mni`, `dataset-metadata`, `derivative-metadata` |
| `tier` | `core` = task `trust` (the RSA/crossnobis target) plus anat and all metadata. `transport` = task `sharedreward`, the independent source for the functional fingerprint. |
| `source` | `raw` (BIDS raw) or `fmriprep` (derivative) |
| `bytes`, `md5` | read from S3 object metadata for snapshot 2.1.1; all 392 rows carry an md5 |
| `sha256` | intentionally empty — S3 does not publish it |
| `url`, `s3_uri` | pinned snapshot URL, and the S3 equivalent |

To change the subset — different subjects, or the `ultimatum` task — regenerate
the manifest against the snapshot rather than hand-editing byte counts.
`DECISION.md` §6 records what a correct manifest was checked against; the
manifest previously in this directory had invented derivative rows and invented
byte counts, which is what §1 is about.

## What slice 2 computes

Full sketch in [`DECISION.md`](DECISION.md) §5. In short:

1. **Per-run GLM.** For each subject and each of the 5 `trust` runs, fit the 9
   conditions (3 partners × {choice} ∪ 3 partners × {defect, recip}) against
   217 volumes at TR 2.02 s, with fMRIPrep confounds as nuisance and
   `missed_trial` as a regressor of no interest. Gives `B ∈ R^{9 × V}` per run
   on the shared 66×78×61 MNI lattice.
2. **P^A — anatomical transport.** Sparse, row-stochastic, MNI voxel → group
   searchlight centre, **with a sink column**. Voxels outside a subject's run
   mask or outside every searchlight put all their mass in the sink, so each
   subject's budget closes exactly and lost coverage stays visible in the
   ledger instead of vanishing from it.
3. **P^F — functional transport.** The same sparse support reweighted by a
   functional fingerprint, renormalized, same sink. Cross-fitted two ways:
   across task (fit on `sharedreward`, evaluate on `trust`) and across run
   (fit on `trust` runs 1–2, evaluate on 3–5 — the held-out side keeps 3 runs,
   which is what makes an unbiased crossnobis RDM computable there).
4. **η_transport**, the held-out improvement of P^F over P^A. Cross-fitted by
   construction, **may be negative, reported unclamped**.
5. **Population layer.** `plan_population()`, OLS over the 12 subjects with one
   QR, and the subject-Gram heterogeneity trick (12×12, rank ≤ 11) instead of
   any P×P covariance.

### Why this dataset makes the transport honest

Every subject's functional data is already on **one common 66×78×61 lattice**
(2.973 × 2.973 × 3.22 mm) — verified across all 84 run masks of the 12 chosen
subjects, a single distinct grid. Subjects differ only through their brain
masks, and they differ a lot: taking one mask per subject (each subject's
`task-trust` run-1 mask), the intersection is 55,003 voxels and the union is
79,978 — so **24,975 voxels, 31 % of the union, are covered for some subjects
and not others**. Pooling all 84 run masks pushes that to 38 %.

That is the point. The geometry is fixed, so P^A is cheap to build and its
row-stochasticity is trivially certifiable — while the sink still carries real,
substantial, subject-varying mass. It is the partial-coverage case
[`design/conservative-geometry-contract.md`](../../design/conservative-geometry-contract.md)
§7.4 exists for, rather than a permutation dressed up as a transport.

## Contract references

- [`design/conservative-geometry-contract.md`](../../design/conservative-geometry-contract.md)
  **§7 "Transport readiness"** — the seven things a population layer needs.
  §7.1 per-row frame metadata (missing, a WS-D deliverable), §7.2 the
  Frobenius-consistent packed codec, §7.3 the subject-Gram eigen trick,
  §7.4 row-stochastic transport with a **required** sink, §7.5 budget-vs-density
  as a declared field.
- [`.planning/2026-08-17-feedback-assessment.md`](../../.planning/2026-08-17-feedback-assessment.md)
  **Part 2** — corrections this slice must encode. Especially item (4)
  transport semantics, item (5) transported coherent/configuration is a ledger
  of *native-node* coherence and must be labelled so, and item (7)
  η_transport must be cross-fitted and may be negative. Item (10) is the ticket
  this directory answers.

## Citation

CC0, so no permission is needed, but `dataset_description.json` asks that
Fareri et al. (2022, *NeuroImage*, `10.1016/j.neuroimage.2022.119267`) and
Smith et al. (2024, *Scientific Data*) be cited. Cite fMRIPrep 21.0.2
(Esteban et al., 2019) for the derivatives and the dataset DOI
`10.18112/openneuro.ds003745.v2.1.1` for the data.
