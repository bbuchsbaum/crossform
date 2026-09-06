# Population calibration certification

Status: normative matched-simulation contract

Contract version: `population-calibration-v1`

Date: 2026-08-21

This court measures the conditional population estimators currently shipped by
Crossform. It does not convert conditional-on-realized-transport inference into
marginal inference over transport learning, coverage selection, or a source
population.

The target is the zero slope in a two-column population model with 24 planned
subjects. Each of 500 replications is evaluated by classical OLS, HC3, and a
399-draw null-imposed Rademacher wild bootstrap on the identical simulated
dataset. The recorded outputs are bias, empirical and reported standard error,
95% interval coverage with binomial Monte Carlo standard error, rejection rate
with Monte Carlo standard error, and failure rate.

Eight declared regimes cover Gaussian homoskedastic errors, covariate-linked
heteroskedasticity, heavy tails plus an influential observation, unequal
transport quality under fixed and cross-fitted labels, independent coverage,
and informative coverage under fixed and cross-fitted labels. Transport quality
is observed conditioning metadata; the simulation does not propagate the
uncertainty of estimating a transport.

The prespecified gates require all methods on every scenario, Gaussian coverage
within 0.04 of 0.95, an HC3 two-MCSE lower coverage bound of at least 0.88 in
regimes whose marginal target is identified, at least 0.02 HC3 improvement over
classical OLS under covariate-linked heteroskedasticity, both transport labels,
and at most 0.02 numerical failure. Informative coverage is deliberately a
refusal gate: none of the implemented intervals is licensed for an unconditional
marginal population claim there, even if its conditional interval happens to
cover in a finite run.

Changing the generator, seeds, replication counts, bootstrap distribution,
target, scenario inventory, or thresholds requires regenerating the results and
checksum manifest and updating the evidence ledger. The replicate artifact
retains `dataset_id`; exactly three method rows per identifier prove paired
comparison rather than independent method-specific simulations.
