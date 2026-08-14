# Smoke-tier comparison report

Date: 2026-08-13
Tier: smoke (VT mask, 577 centres)
Data: Haxby 2001 subject 1, `data.pymvpa.org`, `subj1-2010.01.14.tar.gz`,
SHA-256 `3c14bd7fad6c869e5d9b81739e24e21bb9feb6a7abef27682bbf131b6b4bec5c`
Versions: effectagram 0.0.0.9000 (source tree), rMVPA 0.1.2, neuroim2 0.19.0, R 4.5.1

## Verdict

**effectagram and rMVPA compute the same matched estimand to 8.88e-16 at every
one of 577 searchlight centres, estimator against estimator.** The shared
second-order score is bit-identical (max absolute difference exactly 0). This
is the exemplar's condition 1, met.

Getting there required rejecting the obvious pairing. The estimand the README
specifies — a 1 − Pearson correlation-distance RDM scored by a Spearman
correlation — is expressible by **neither** package. And rMVPA's `rsa_model`,
the function the README names for this arm, cannot express the matched
estimand either: it never aggregates to the condition level. The match is
between effectagram's `rdm()` over a cross-run geometry and rMVPA's
*crossnobis* machinery, which is a different rMVPA entry point than the one
originally planned.

Five named semantic differences, each established by reading source and
confirmed by running code, account for why the two packages' **`rsa_model`-vs-
`rdm()` native outputs** correlate only r ≈ 0.54. That is the exemplar's
condition 2, also met. No number in either arm was adjusted to produce any
agreement reported here.

## What was compared

Both arms consume the identical object written by `01-prepare-data.R`: twelve
8 × 577 per-run condition-mean matrices, built once from within-run z-scored
time series averaged over block volumes shifted by 2 TR (5 s at TR = 2.5 s),
9 volumes per condition per run. No preprocessing, HRF, or averaging choice
differs between arms.

Searchlight geometry is identical, and this was verified rather than assumed:
at radius 11.25 mm (the README's "3 voxels" along the 3.75 mm in-plane axes),
`neuroim2::searchlight_indices` (effectagram's route) and
`neuroim2::searchlight` (rMVPA's route) return **bit-identical sphere
membership at all 577 centres**, in the same centre order. Sphere sizes: min 6,
median 45, max 88 voxels.

Model RDM: binary animate (`cat`, `face`) vs inanimate over the eight stimulus
categories.

## The matched estimand

Cross-validated squared Euclidean distance under the identity metric, averaged
over the 66 off-diagonal run pairs, with local (1/P) normalisation:

```
d2(i,j) = mean over run pairs r<s of (1/P) * sum_v (b_i^r[v] - b_j^r[v]) * (b_i^s[v] - b_j^s[v])
```

Three independent implementations of it were run, and all three agree.

| # | comparison | max abs difference | centres above 1e-8 gate |
|---|---|---|---|
| 1a | effectagram `rdm()` vs an independent reference loop | 1.33e-15 | 0 / 577 |
| **1b** | **effectagram `rdm()` vs rMVPA's own crossnobis estimator** | **8.88e-16** | **0 / 577** |
| 2 | shared `spearman_rsa()` applied to each arm's own RDM | **0** (exact) | 0 / 577 |

Scale for context: max |d2| = 1.229, so 1b is a relative difference of about
7e-16 — floating-point noise, not agreement within a tolerance.

**1a** is effectagram's compiled sparse frame contraction against a plain
nested loop written from the definition above and fed by rMVPA's
`get_searchlight` spheres. The loop depends on neither package's estimator, so
it arbitrates between them.

**1b** is the comparison that earns a replacement-map row, because both sides
are package estimators:

- effectagram: `relation()` → `cross_partitions()` → `plan_geometry()` →
  `rdm()`, 0.48 s for all 577 centres.
- rMVPA: `compute_crossvalidated_means_sl(estimation_method = "crossnobis")`
  → `rMVPA:::compute_crossnobis_distances_sl(U, P_voxels = ncol(X))`,
  0.4 s for all 577 centres.

The two arrive at the same numbers by visibly different routes: effectagram
contracts a sparse additive frame against cross-run feature atoms, while rMVPA
forms `G_cv = (tcrossprod(U_sum) − tcrossprod(U_stacked)) / (M(M−1))` and reads
`(G_ii + G_jj − 2G_ij)/P`. rMVPA averages over *ordered* partition pairs where
effectagram averages over *unordered* ones; these coincide because the summand
is symmetric in (r, s). Pair naming and ordering also coincide — rMVPA's
`"<right>_vs_<left>"` lower-triangle order is `combn()` order — and `03`
asserts that rather than trusting it.

**Neither `rsa_model` nor `vector_rsa_model` can express this estimand.**
`rsa_model` builds one pooled 96 × 96 RDM and keeps its between-run cell pairs;
`vector_rsa_model` works at the trial level and returns a single scalar.
`contrast_rsa_model` computes the correct `G_cv` internally but only exposes
contrast betas. The crossnobis helper is the only route, and
`compute_crossnobis_distances_sl` is **not exported**, so `03` calls it with
`:::`.

## Native estimands — reported, deliberately not gated

rMVPA's `rsa_model` and effectagram's cross-run geometry compute **different
things**. Both were run natively and compared as maps:

| rMVPA `distmethod` | r vs effectagram Spearman | ρ vs effectagram Spearman | r vs effectagram OLS `rsa()` | max abs diff | rMVPA median |
|---|---|---|---|---|---|
| pearson | 0.540 | 0.512 | 0.510 | 0.665 | 0.093 |
| spearman | 0.562 | 0.521 | 0.518 | 0.656 | 0.091 |

A map correlation near 0.54 is the expected consequence of the five
differences below, not an unexplained discrepancy. Both arms recover the same
qualitative VT animacy structure; they do not compute the same statistic.
effectagram's Spearman RSA over the matched estimand: median 0.152, range
[−0.206, 0.786].

## Divergence taxonomy

Each entry was established by reading the packages' source, not inferred from
the size of a discrepancy. All five affect the *native* comparison only; none
affects the matched estimand, which agrees exactly.

**1. `aggregation-level`.** `rsa_model` builds one pooled 96 × 96 RDM per sphere
(12 runs × 8 conditions as 96 independent rows) and keeps its 4224 between-run
cell pairs. effectagram aggregates over run pairs *first*, yielding 28
condition-level values. This is not a reparameterisation: a correlation over
4224 cell pairs is not a function of the 28 cell means, because its denominator
uses the full within-cell spread. The two cannot coincide even on identical
patterns.

**2. `distance-convention`.** `rsa_model` uses 1 − correlation between
individual single-run patterns. effectagram uses the cross-validated squared
Euclidean distance under the identity metric. Correlation distance rescales
each pattern to unit norm *inside the sphere*, which is nonlinear in the
patterns and therefore **not expressible in effectagram's bilinear core** at
all. Per-sphere mean-centering alone would be expressible (as a non-diagonal
metric `K = I − 11'/P`); the unit-norm denominator is not.

**3. `second-order-coupling`.** In rMVPA, `distmethod` sets **both** the neural
RDM metric **and** the second-order correlation method; `regtype` is inert for
`"pearson"`/`"spearman"`. Verified on this exemplar's own data by running all
four combinations over the 577 VT centres:

| comparison | max abs difference over 577 centres |
|---|---|
| `distmethod="pearson"`, `regtype` pearson vs spearman | **0** (exact) |
| `distmethod="spearman"`, `regtype` pearson vs spearman | **0** (exact) |
| `regtype` fixed, `distmethod` pearson vs spearman | 0.0282 |

So a Pearson (1 − r) RDM scored by a Spearman second-order correlation — the
README's literal estimand — **cannot be expressed in a single `rsa_model`
call**. The cause is `R/rsa_model.R:270-275`, where the second-order
correlation is taken with `method = obj$distmethod`. The RSA vignette's comment
describing `regtype` as the model-vs-neural comparison method is therefore
incorrect.

The coupling is a property of `rsa_model`, not of rMVPA. `vector_rsa_model`
*does* decouple them —
`vector_rsa_model(distfun = cordist(method = "pearson"), rsa_simfun = "spearman")`
sets the neural distance and the second-order statistic independently. It still
does not produce the README's estimand, because it correlates *row i* of the
96 × 96 distance matrix against row i of the expanded model, restricted to
other-run columns, and averages those 96 per-trial scores
(`R/vector_rsa_model.R:238-264`). That is a mean of per-trial row correlations,
not a correlation over the lower triangle of a condition RDM. Its output is a
single `rsa_score` scalar per centre.

**4. `run-pair-treatment`.** Both arms exclude within-run information, but at
different levels. `rsa_model` drops within-run *cell* pairs from a pooled
vector (4224 of 4560 kept, dropping 12 × 28 = 336). effectagram averages over
the 66 unordered *run* pairs. These agree on what is excluded, not on how what
remains is weighted.

**5. `effectagram-second-order`.** effectagram's native `rsa()` is an ordinary
least-squares regression in RDM space, not a rank correlation, so effectagram
cannot natively produce a Spearman RSA either. The Spearman statistic used in
comparison [2] is supplied by the exemplar's own shared `spearman_rsa()`,
applied identically to both arms. That `rsa()` is exactly the OLS functional
`(X'X)^{-1}X'd` was verified to 1.1e-16, which is what licenses the covariance
transport in `05`.

### Consequence for the README's estimand

The README's stated estimand is not the nearest expressible one in either
system. The nearest estimand **both** systems express exactly is the matched
one above: a condition-level cross-run RDM under cross-validated squared
Euclidean distance, plus a rank statistic applied identically to both arms.
The residual difference from the README's text is the distance convention, and
it is a hard expressibility boundary of effectagram's bilinear core, not a
configuration choice.

## Error-channel arm (05) — the surplus, and its declared limits

Not a comparison row: a different, better-defined estimand that rMVPA does not
expose. The precomputed condition means cannot support it — averaging discards
the residuals — so `05` refits the *same* effects from the raw labelled volumes
with `lm_relation_fit()`. Because the design is a disjoint indicator matrix,
the refitted coefficients are the condition means, and the refit reproduces
`02`'s RDM to **4.44e-16**. The point estimand is unchanged; what is added is
the residual channel that averaging had thrown away.

Run over 120 evenly spaced VT centres (of 577).

**Distances with analytic standard errors.** Fixed identity noise precision,
`rdm_sampling_covariance(target = "plugin")`:

| quantity | value |
|---|---|
| RDM standard error | median 0.0207, range [0.0074, 0.0788] |
| distance point estimates | range [−0.0758, 1.1686] |
| distances with \|d\|/SE > 2 | 76.4 % |

**Transport to the RSA coefficient.** effectagram's `rsa()` is exactly the OLS
functional `(X'X)^{-1}X'd` (verified to 1.1e-16), so the RSA coefficient is a
*linear* functional of the 28 distances and its variance is an exact transport
of the RDM covariance — `sampling_covariance(cv, operation = "transport")`,
no delta-method approximation:

| quantity | value |
|---|---|
| animacy RSA coefficient | median 0.0372 |
| analytic standard error | median 0.0121 |
| \|z\| | median 3.47, max 11.61 |
| centres with \|z\| > 1.96 | 72.5 % |

Note the boundary this exposes: the Spearman score used for the parity
comparison is **not** a linear functional, so it gets no transport. The error
bar is available for the OLS coefficient and not for the rank statistic — a
real constraint, not an oversight.

**Learned-metric crossnobis.** `plan_crossnobis(metric = shrinkage_precision(0.2))`
with a face − house contrast, all 577 centres in 152 s: median 0.300, range
[−0.135, 0.739], positive at 99.8 % of centres — the expected VT face/house
separation.

**Three capability refusals, captured verbatim.** These are as much the point
of the arm as the numbers: the package declines to calibrate what it cannot
calibrate rather than returning a plausible number.

1. Analytic covariance on the learned-metric crossnobis plan:
   > Analytic RDM sampling covariance is unavailable because this plan's neural
   > metric was learned. Version 0.1 does not propagate metric-estimation
   > uncertainty; use a geometry plan with one common fixed metric or retain
   > this result as a signed point estimate.

2. Analytic covariance from precomputed effects (the `02`/`03` relation):
   > Sampling covariance is unavailable because this evidence plan has only a
   > precomputed relation and no error channel. Refit raw observations with
   > `lm_relation_fit()` or supply a validated, identity-bound external error
   > channel; beta matrices alone cannot recover residual uncertainty.

3. `target` omitted:
   > Choose `target = "plugin"` for the partition-mean signal policy or
   > `target = "null"` for the fixed-zero null; the signal-dependent covariance
   > target is not inferred.

Refusal 2 is the one that justifies the whole arm's structure: it is exactly
why `05` had to refit from raw volumes instead of reusing `01`'s means.

**Cost, measured rather than guessed.** 9.03 s/node mean over 120 nodes
(1084 s total). An earlier run looked like it was degrading badly (marginal
cost rising 7.0 → 10.8 s/node), which prompted instrumenting the loop. With
the machine quiet the trend is mild — first 10 nodes 7.73 s/node, last 10
8.93 s/node, regression slope 0.0048 s per node — and RSS oscillates between
819 MB and 1466 MB with a slope of 0.6 MB/node and no net growth. **There is no
memory leak**; the earlier apparent degradation was contention from unrelated
long-running R processes on the machine, and the initial reading was wrong.
The genuine cost point stands: one node per call at ~9 s makes a full 577-node
sweep about 90 minutes, so `05` defaults to 120 nodes.

## Environment hazards found

**A subagent installed rMVPA from source during this session.** There are now
two rMVPA 0.1.2 builds, and the newer one shadows the older:

```
/Users/bbuchsbaum/Library/R/arm64/4.5/library/rMVPA          built 2026-08-14 00:59 UTC  <- loads
/Library/Frameworks/R.framework/Versions/4.5-arm64/.../rMVPA built 2026-04-25 17:21 UTC
```

This was not requested and is flagged for the maintainer to accept or undo; it
was left in place rather than deleted. It matters for the results in two ways:

1. **The crossnobis arm requires the newer build.** The April build refuses
   crossnobis without a whitening matrix, and even when given `W = diag(P)` it
   derives folds from leave-one-run-out *training* means rather than per-run
   means. Those overlap across folds, so cross-fold inner products do not
   cancel noise and the distances come out positively biased — the kind of
   silent wrongness this exemplar exists to catch.

   `03` therefore checks the **semantics rather than the version**: it asserts
   that `compute_crossvalidated_means_sl(estimation_method = "crossnobis")`
   returns exactly this run's per-run condition means, and stops with an
   explanatory error if not. On the build used here that assertion passes at
   max abs difference **0** (against leave-one-out means the gap is 2.14), and
   crossnobis over pure white noise averages **−0.011**, i.e. unbiased around
   zero as an independent-fold estimator must be. That is the direct evidence
   that comparison [1b] is measuring what it claims.
2. **`rsa_model` results are nearly but not exactly build-stable.** Across the
   two builds: `distmethod="pearson"` agrees to 1.11e-16, but
   `distmethod="spearman"` differs by **2.97e-4**. Small, but a real
   behavioural change; the committed numbers come from the newer build.

**A latent geometry bug in this exemplar, caught by the newer rMVPA.** `03`
originally built the 96-volume dataset with `NeuroSpace(c(dim(mask), 96))`,
which silently drops spacing, origin, and orientation. The April build accepted
it; the current build correctly rejects it with `Spatial geometry mismatch
between `train_data` and `mask``. Fixed to
`neuroim2::add_dim(neuroim2::space(vt_mask), n_obs)`. The earlier results were
not affected — spheres come from the mask's own geometry and extraction is
index-based — but the construction was wrong and is now correct.

**rMVPA shard backend, intermittent.** One run of three successive
`run_searchlight` calls in a single session aborted with

```
ERROR [2026-08-13 20:49:12] No successful results to combine (schema combiner). 577 ROIs failed.
Error in combine_schema_standard(model_spec, good_results, bad_results) :
  No valid results for standard searchlight: all ROIs failed to process
```

It did not reproduce: a later session ran four `run_searchlight` calls back to
back cleanly, and driving the same model through `mvpa_iterate` directly
reported `error = FALSE` for every ROI. Each call logs `shard backend:
preparing shared memory for dataset`, so a shared-memory resource not released
between calls is the natural suspect, but this was not pinned down. Recorded
because `03` runs two searchlights in one session.

**effectagram: no defects found.** No incorrect result and no crash was
observed anywhere in this exemplar. Three API frictions worth noting, none of
them wrong behaviour:

- `rdm_sampling_covariance()` evaluates one frame node per call (`at = 1L`), so
  covering a searchlight map is an explicit loop at ~9 s/node.
- `rdm(plan)` with a domain-wide fixed metric costs about 3 minutes over 577
  centres, against 0.48 s for the same view without a metric. The metric here
  is `diag(577)`, which `R/metric.R:171-173` does detect as native-diagonal, so
  the gap is larger than the diagonal fast path would suggest.
- `residual_df(fit)` requires a `partition` argument that its usage was not
  obvious about.

## Timings (smoke tier, sequential, single core)

| step | wall clock |
|---|---|
| effectagram `rdm()`, 577 centres | 0.48 s |
| effectagram `rsa()`, 577 centres | 0.84 s |
| rMVPA crossnobis RDM, 577 centres | 0.40 s |
| rMVPA `rsa_model` searchlight, `distmethod="pearson"` | 3.9 s |
| rMVPA `rsa_model` searchlight, `distmethod="spearman"` | 7.1 s |
| independent reference loop over 577 spheres | 3.5 s |

Not a benchmark claim — no memory instrumentation, single runs, and the
`rsa_model` rows compute a different estimand. Recorded only to show the parity
tier is cheap (about 20 s end to end).

## Reproducing

```sh
cd exemplars/haxby2001
Rscript 00-download.R            # ~300 MB, checksum-gated, idempotent
Rscript 01-prepare-data.R        # ~2 min
Rscript 02-effectagram-searchlight.R
Rscript 03-rmvpa-searchlight.R
Rscript 04-compare.R             # exits non-zero if a gate fails
Rscript 05-crossnobis-uncertainty.R   # ~22 min at the default 120 nodes
```

`05` defaults to 120 evenly spaced VT centres.
`UNCERTAINTY_NODES=all Rscript 05-crossnobis-uncertainty.R` runs all 577
(about 90 minutes). Scripts 02–04 are the commit gate and take about 20 s
together once `01` has run.
