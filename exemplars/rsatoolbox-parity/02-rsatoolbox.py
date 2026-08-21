#!/usr/bin/env python
"""02-rsatoolbox.py -- the external arm.

Reads the CSVs written by 01-fixture.R and recomputes, with the pinned
``rsatoolbox`` 0.3.2, the two quantities that live in the common fixed-linear
subset:

  1. the crossvalidated Mahalanobis (crossnobis) RDM, via
     ``rsatoolbox.rdm.calc_rdm_crossnobis``;
  2. a linear RSA readout of that RDM.

No reticulate, no shared process: the only channel between the two languages
is the CSV directory.

Convention notes, recorded here because they are the whole substance of a
parity claim:

* ``calc_rdm_crossnobis`` crossvalidates **leave-one-fold-out**: for each fold
  it contracts that fold's condition means against the mean of the remaining
  folds.  ``crossform``'s ``cross_partitions()`` instead declares the
  C(P, 2) unordered run pairs with uniform weight.  With a balanced design
  (one estimate per condition per fold) the two are algebraically identical --
  LOO averages (1/P) sum_f (1/(P-1)) sum_{g != f}, which is the uniform mean
  over ordered pairs, and the fixed precision is symmetric, so ordered and
  unordered means coincide.  The script computes the explicit all-pairs mean
  as a third, independent oracle so the identity is checked rather than
  asserted.
* ``_calc_rdm_crossnobis_single`` divides by ``n_channels``.  ``crossform``
  gets the same division from the *frame*: a ``whole_brain``/``regions`` frame
  with ``normalization = "local"`` has row sums of one, so each measurement's
  weights are 1/|support|.  The two normalisations are the same number
  arrived at from different directions.
* ``noise=`` takes the precision restricted to the measurement's support,
  which is what crossform's support-streamed kernel contracts
  (``K[support, support]``, not the inverse of the covariance submatrix).
* ``rsatoolbox``'s ``fit_regress`` optimises a *normalised* cosine objective:
  ``pool_rdm(..., 'cosine')`` divides the data vector by its RMS, so its
  theta equals the no-intercept OLS coefficients divided by
  ``sqrt(mean(d**2))``.  That scale factor is recorded and undone here; the
  unrescaled comparator is plain least squares on the vectorised RDMs, which
  is exactly ``crossform::rsa()``'s definition.
"""

from __future__ import annotations

import csv
import os
import sys
from importlib import metadata

import numpy as np
from rsatoolbox.data import Dataset
from rsatoolbox.data.noise import prec_from_residuals
from rsatoolbox.model import ModelWeighted
from rsatoolbox.model.fitter import fit_regress
from rsatoolbox.rdm import RDMs
from rsatoolbox.rdm.calc import calc_rdm_crossnobis

HERE = os.path.dirname(os.path.abspath(__file__))
RESULTS = os.path.join(HERE, "results")


def read_csv(name):
    with open(os.path.join(RESULTS, name), newline="") as handle:
        return list(csv.DictReader(handle))


def numeric_block(rows, prefix="V"):
    columns = [k for k in rows[0] if k.startswith(prefix)]
    columns.sort(key=lambda k: int(k[len(prefix) :]))
    return np.array([[float(r[c]) for c in columns] for r in rows])


def main() -> int:
    meta = {r["key"]: float(r["value"]) for r in read_csv("fixture-meta.csv")}
    n_voxels = int(meta["n_voxels"])
    n_runs = int(meta["n_runs"])
    residual_df = int(meta["residual_df"])

    beta_rows = read_csv("betas.csv")
    betas = numeric_block(beta_rows)
    conds = np.array([r["condition"] for r in beta_rows])
    runs = np.array([r["run"] for r in beta_rows])
    assert betas.shape[1] == n_voxels

    precision = np.array(
        [[float(v) for v in r.values()] for r in read_csv("precision.csv")]
    )
    assert precision.shape == (n_voxels, n_voxels)

    residual_rows = read_csv("residuals.csv")
    residuals = numeric_block(residual_rows)

    region_rows = read_csv("regions.csv")
    regions = np.array([r["region"] for r in region_rows])

    model_rows = read_csv("model-rdms.csv")
    pair_left = [r["left"] for r in model_rows]
    pair_right = [r["right"] for r in model_rows]
    model_vectors = {
        "category": np.array([float(r["category"]) for r in model_rows]),
        "animacy": np.array([float(r["animacy"]) for r in model_rows]),
    }

    conditions = sorted(set(conds))
    n_cond = len(conditions)
    # The pair order crossform reports must be the row-major upper triangle,
    # or every element-wise comparison downstream is meaningless.
    triu = np.triu_indices(n_cond, 1)
    expected_left = [conditions[i] for i in triu[0]]
    expected_right = [conditions[j] for j in triu[1]]
    assert expected_left == pair_left, (expected_left, pair_left)
    assert expected_right == pair_right, (expected_right, pair_right)

    # ---- The precision estimator, checked against rsatoolbox's own --------
    rsatoolbox_precision = prec_from_residuals(
        residuals, dof=residual_df, method="full"
    )
    precision_gap = float(np.max(np.abs(rsatoolbox_precision - precision)))
    precision_rel = precision_gap / float(np.max(np.abs(precision)))

    # ---- Measurement supports --------------------------------------------
    supports = {"whole_brain": np.arange(n_voxels)}
    for label in sorted(set(regions)):
        supports[label] = np.flatnonzero(regions == label)

    rdm_records = []
    rsa_records = []
    oracle_gap = 0.0
    for name, support in supports.items():
        noise = precision[np.ix_(support, support)]
        dataset = Dataset(
            measurements=betas[:, support],
            obs_descriptors={"conds": conds, "run": runs},
            channel_descriptors={"voxel": support.astype(str)},
        )
        rdms = calc_rdm_crossnobis(
            dataset, descriptor="conds", cv_descriptor="run", noise=noise
        )
        values = np.asarray(rdms.get_vectors()).ravel()
        assert list(rdms.pattern_descriptors["conds"]) == conditions

        # Independent oracle: uniform mean over the C(P, 2) unordered run
        # pairs of (delta_i - delta_j)' K (delta_i - delta_j) / |support|.
        fold_names = sorted(set(runs))
        patterns = {
            f: np.array(
                [betas[(conds == c) & (runs == f)][0][support] for c in conditions]
            )
            for f in fold_names
        }
        accumulator = np.zeros(len(triu[0]))
        n_edges = 0
        for a in range(n_runs):
            for b in range(a + 1, n_runs):
                first = patterns[fold_names[a]]
                second = patterns[fold_names[b]]
                kernel = first @ noise @ second.T
                square = (
                    np.expand_dims(np.diag(kernel), 0)
                    + np.expand_dims(np.diag(kernel), 1)
                    - kernel
                    - kernel.T
                )
                accumulator += square[triu] / len(support)
                n_edges += 1
        oracle = accumulator / n_edges
        oracle_gap = max(oracle_gap, float(np.max(np.abs(oracle - values))))

        for k, value in enumerate(values):
            rdm_records.append(
                {
                    "measurement": name,
                    "pair": k + 1,
                    "left": pair_left[k],
                    "right": pair_right[k],
                    "rsatoolbox": repr(float(value)),
                    "allpairs_oracle": repr(float(oracle[k])),
                }
            )

        # ---- Linear RSA ---------------------------------------------------
        design = np.column_stack(
            [np.ones_like(values), model_vectors["category"], model_vectors["animacy"]]
        )
        beta, *_ = np.linalg.lstsq(design, values, rcond=None)
        for term, value in zip(("(Intercept)", "category", "animacy"), beta):
            rsa_records.append(
                {
                    "measurement": name,
                    "term": term,
                    "rsatoolbox": repr(float(value)),
                    "route": "numpy_lstsq",
                }
            )

        bare_design = np.column_stack(
            [model_vectors["category"], model_vectors["animacy"]]
        )
        bare_beta, *_ = np.linalg.lstsq(bare_design, values, rcond=None)
        for term, value in zip(("category", "animacy"), bare_beta):
            rsa_records.append(
                {
                    "measurement": name,
                    "term": "nointercept:" + term,
                    "rsatoolbox": repr(float(value)),
                    "route": "numpy_lstsq",
                }
            )

        # rsatoolbox's own weighted-model fitter, with its cosine
        # normalisation undone. `pool_rdm(..., 'cosine')` scales the data
        # vector by 1 / sqrt(mean(d ** 2)); `normalize=False` keeps theta on
        # that scale, so multiplying back recovers the OLS coefficients.
        model = ModelWeighted(
            "weighted",
            RDMs(
                np.stack([model_vectors["category"], model_vectors["animacy"]]),
                dissimilarity_measure="crossnobis",
                pattern_descriptors={"conds": np.array(conditions)},
            ),
        )
        theta = fit_regress(model, rdms, method="cosine", normalize=False)
        cosine_scale = float(np.sqrt(np.mean(values**2)))
        for term, value in zip(("category", "animacy"), theta):
            rsa_records.append(
                {
                    "measurement": name,
                    "term": "nointercept:" + term,
                    "rsatoolbox": repr(float(value) * cosine_scale),
                    "route": "fit_regress_cosine_rescaled",
                }
            )
        rsa_records.append(
            {
                "measurement": name,
                "term": "cosine_scale",
                "rsatoolbox": repr(cosine_scale),
                "route": "diagnostic",
            }
        )

    with open(os.path.join(RESULTS, "rsatoolbox-rdm.csv"), "w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "measurement",
                "pair",
                "left",
                "right",
                "rsatoolbox",
                "allpairs_oracle",
            ],
        )
        writer.writeheader()
        writer.writerows(rdm_records)

    with open(os.path.join(RESULTS, "rsatoolbox-rsa.csv"), "w", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=["measurement", "term", "rsatoolbox", "route"]
        )
        writer.writeheader()
        writer.writerows(rsa_records)

    with open(os.path.join(RESULTS, "rsatoolbox-meta.csv"), "w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["key", "value"])
        writer.writerow(["python", sys.version.split()[0]])
        writer.writerow(["rsatoolbox", metadata.version("rsatoolbox")])
        writer.writerow(["numpy", metadata.version("numpy")])
        writer.writerow(["scipy", metadata.version("scipy")])
        writer.writerow(["prec_from_residuals_max_abs_diff", repr(precision_gap)])
        writer.writerow(["prec_from_residuals_rel_diff", repr(precision_rel)])
        writer.writerow(["loo_vs_allpairs_max_abs_diff", repr(oracle_gap)])
        writer.writerow(["cv_scheme", "leave_one_fold_out"])
        writer.writerow(["distance_normalization", "divide_by_n_channels"])

    print(
        f"rsatoolbox {metadata.version('rsatoolbox')} / "
        f"numpy {metadata.version('numpy')} / python "
        f"{sys.version.split()[0]}"
    )
    print(
        f"prec_from_residuals(method='full', dof={residual_df}) vs the "
        f"R pooled precision: max abs diff = {precision_gap:.3e} "
        f"(relative {precision_rel:.3e})"
    )
    print(
        "calc_rdm_crossnobis (leave-one-run-out) vs explicit all-pairs "
        f"mean: max abs diff = {oracle_gap:.3e}"
    )
    print(f"wrote {len(rdm_records)} RDM rows and {len(rsa_records)} RSA rows")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
