# Null coverage of the between-subject uncertainty layer

Recorded evidence for `population_uncertainty()`'s between-subject standard
error, its residual degrees of freedom, and the `t` statistic it reports.

```sh
Rscript benchmarks/run-population-null-coverage.R 2000 benchmark-results
Rscript benchmarks/promote-artifacts.R benchmark-results .
```

Artifacts: `benchmark-results/population-null-coverage.rds` (local, carries
provenance) and `inst/extdata/certification/population-null-coverage-summary.csv`
(committed receipt). `tests/testthat/test-population-uncertainty.R` reads both
and skips loudly when either is absent.

## What the simulation exercises, and what it does not

It exercises the **group layer** and nothing else: `.population_ols()`, the
unscaled covariance `.population_between_uncertainty()` records at fit time,
and the standard error, `t` and interval `.population_between_statistics()`
computes. Those are the shipped functions, called directly rather than
re-implemented.

It does **not** run the subject geometry kernels. The per-subject transported
query values `z_{i,u,k}` are simulated directly. That is deliberate and
sufficient: the object under test is the sampling behaviour of a standard error
computed from `N` numbers per group node and query, and query, transport and
fit commute (`population-form-v1` §3), so the group layer sees exactly the
array the script hands it. Running a geometry plan 2,000 times would measure
the compiler.

The shortcut is bound to the public verb by a **wiring check** the runner
refuses to proceed without: a real six-participant population is fitted through
`plan_population()` and `estimate_population()`, the fast path is run on the
transported values that fit produced, and the two standard errors must agree.
Measured `max |SE(fast) − SE(population_uncertainty)| = 0` (bit-identical), and
against base R's own `lm()` on the same response, `3.5e-18`.

## The generative model

Four group nodes (three plus the sink), two queries, group model
`~ null_x + signal_x` with both covariates redrawn each replication.

| quantity | value |
|---|---|
| null term | `null_x`, true coefficient **0** |
| signal term | `signal_x`, true coefficient **0.6** |
| between-subject sd | 0.8 |
| within-subject sd | 0.6 |
| across-node correlation | AR(1), ρ = 0.5 over group nodes |
| replications | 2,000 per cell |
| nominal level | 95%, two-sided |

Every distributional summary is computed from **one designated cell**
(`group1`, first query) per replication, so the 2,000 draws are independent,
the Monte Carlo standard error is exact, and the KS statistic is valid.
Coverage pooled over all node-query cells is reported beside it and is labelled
dependent, because cells within one replication share participants.

The `heteroskedastic` arm keeps everything above and makes each participant's
noise scale a function of the covariates the group model regresses on,
`sd_i ∝ exp(0.7 (null_x_i + signal_x_i)/√2)`, rescaled so the mean variance
matches the `gaussian` arm. The two arms differ in shape, not in magnitude.
This is the transport-heterogeneity analogue: a participant whose warp is poor
carries a noisier transported value, and warp quality is not independent of
age, motion or head size.

## Measured (2026-08-20, 2,000 replications per cell)

### `gaussian` — the group model is correct

| N | residual df | null coverage | MCSE | null rejection | signal coverage | power | KS vs t_df | KS p |
|---|---|---|---|---|---|---|---|---|
| 6 | 3 | **0.9485** | 0.0049 | 0.0515 | 0.9565 | 0.142 | 0.0209 | 0.344 |
| 8 | 5 | **0.9500** | 0.0049 | 0.0500 | 0.9555 | 0.213 | 0.0145 | 0.796 |
| 12 | 9 | **0.9520** | 0.0048 | 0.0480 | 0.9480 | 0.374 | 0.0151 | 0.753 |
| 24 | 21 | **0.9480** | 0.0050 | 0.0520 | 0.9510 | 0.747 | 0.0201 | 0.393 |

Every coverage figure is within one Monte Carlo standard error of 0.95, and no
KS test comes close to rejecting `t_df`. The null `t`'s sample standard
deviation tracks the reference `sqrt(df/(df−2))` (1.689 against 1.732 at
`df = 3`; 1.063 against 1.051 at `df = 21`).

**This is the arm where the answer was never in doubt.** It says the
arithmetic is right — the pivoted-QR unscaled covariance, the residual df, the
`qt` quantile — and it says nothing about a real study.

### `heteroskedastic` — the participants' noise is linked to the covariates

| N | residual df | null coverage | MCSE | null rejection | signal coverage | power | KS vs t_df | KS p |
|---|---|---|---|---|---|---|---|---|
| 6 | 3 | **0.9230** | 0.0060 | 0.0770 | 0.9255 | 0.206 | 0.0407 | 2.7e-03 |
| 8 | 5 | **0.9180** | 0.0061 | 0.0820 | 0.9280 | 0.303 | 0.0371 | 8.2e-03 |
| 12 | 9 | **0.9095** | 0.0064 | 0.0905 | 0.9160 | 0.470 | 0.0503 | 8.1e-05 |
| 24 | 21 | **0.8850** | 0.0071 | 0.1150 | 0.8760 | 0.746 | 0.0599 | 1.2e-06 |

Coverage degrades, and **it degrades further as the study grows**: the bias
lives in the standard error, not in the sample size, so more participants buy
a tighter interval around the wrong width. At `N = 24` a nominal 5% test
rejects a true null **11.5%** of the time — a 2.3× inflation of the false
positive rate — and the null `t` is distinguishable from `t_df` at
`p = 1.2e-6`. The null `t`'s sample standard deviation is 1.308 against a
reference 1.051 at `df = 21`.

The misspecification here is mild by the standards of a real population study:
one smooth link between noise scale and two covariates, no outlying
participant, no node-dependent heterogeneity, no violated normality. A real
transport varies in displacement, entropy and subject coverage across
participants (`population-form-v1` §7.5), and a group node's residual scatter
is then a mixture over participants who were not measuring the same thing.

## The conclusion the shipped label rests on

`population_uncertainty()` reports `calibration = "uncalibrated"` in every row
of `$between`, on every printed line, and in the Rd. The `gaussian` arm is the
reason the estimator is shipped at all; the `heteroskedastic` arm is the reason
the label does not move. Report the statistic. A p-value derived from it needs
an argument that §7.5's transport diagnostics are benign for the study at hand,
and this package cannot supply that argument.

## Ratchet

`tests/testthat/test-population-uncertainty.R` holds the recorded numbers to
declared bounds:

* `gaussian` null coverage in `[0.93, 0.97]` at every `N`, and every KS
  statistic below `0.03`.
* `heteroskedastic` null coverage below `0.94` at every `N` and below `0.90`
  at `N = 24`, with the `N = 24` rejection rate above `0.09`.

The second set is a floor on a *failure*: it fails if the misspecification arm
stops demonstrating the anticonservatism the documentation cites. Deleting the
demonstration is the failure mode it exists to catch.
