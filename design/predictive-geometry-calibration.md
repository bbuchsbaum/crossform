# Predictive geometry calibration protocol

The fixed-fit court uses seed 2026090403 and a separate precision pilot with
seed 2026090413. Four prespecified generators have six conditions, five neural
features, six independent partition estimates, and a two-dimensional model
span: pure noise, aligned signal, signal outside the model span, and
heterogeneous partition noise with an anisotropic fixed SPD metric. Neural
pairings have three nonuniformly weighted edges on each disjoint three-run
subset. The primary unit is a generated dataset, including all its edges.

Every dataset evaluates the same fixed rank-two/penalty-0.15 procedure under
the supplied kernel, a same-span isotropic kernel and a same-span mismatched
kernel. The zero predictor is exact. This comparison separates kernel
preference from generic rank restriction and shrinkage; no comparator is
selected on test outcomes and no universal superiority claim is required.

For row-independent Gaussian partition errors with feature covariance `V`,
the exact dense target is `G* = B M B'`. Let `d_a` be the total incident test
edge weight of endpoint `a`. The estimator noise second moment is

```
(q+1)/2 * tr(B M V M B') * sum_a d_a^2 sigma_a^2
+ q*(q+1)/2 * tr(M V M V) * sum_edges w_ab^2 sigma_a^2 sigma_b^2.
```

The moment formula is checked independently by finite noise enumeration and
by a separately reported empirical generator-moment check. Define reference
scale `b = ||G*||^2 + noise_second_moment` and practical margin `delta=.02*b`.
The primary error is `gain - (||G*||^2 - ||G*-F||^2)`. Each of 12 generator /
predictor checks must satisfy both `abs(mean(error)) <= 5*MCSE + 1e-10` and
`abs(mean(error)) + 5*MCSE <= delta`. The four generator-moment checks use five
MCSE. The corresponding CLT Bonferroni family bound is below `1e-5`.

The 2,500-dataset pilot selects only sample size. It takes the largest error
SD among the three predictors, targets `5*MCSE <= delta/2`, adds 50% variance
headroom, rounds upward to 1,000 datasets, and uses at least 20,000. Final
counts and all pilot output are frozen before production. Each replicate
uses its own L'Ecuyer stream; extra computations cannot alter later datasets.
An imprecise final run fails; counts and tolerances are not adjusted afterward.

A negative control evaluates each prediction on its own training geometry.
In the null arm, its lower five-MCSE bias bound must exceed `.005*b`. It is
labelled leakage, not an admitted score. Known overlapping origins are also
refused by public-path tests; hidden dependence is a violated assumption that
cannot be diagnosed reliably from numeric values alone.

The runner records production and harness digests, generator details, streams,
all per-dataset values, paired gain comparisons and runtime. Final source
certification reruns this court after source changes. The small testthat
generator/public-path checks establish execution parity, not Monte Carlo
calibration by themselves.

The independent pilot selected fixed production counts before the run: null 20,000; aligned 27,000; outside-model 20,000; heterogeneous-metric 20,000. These counts, the config hash and the complete pilot are recorded in `benchmark-results/predictive-geometry-sampling-pilot.rds`. They are not adjusted using production results.
