# Evidence-sampling and calibration contract

Status: normative statistical architecture contract

Contract version: `evidence-sampling-v1`

Date: 2026-08-13

Last corrected: 2026-08-15 (Section 4, signal term; Section 9, plug-in bias)

This contract governs sampling covariance and within-participant calibration
for estimated evidence pairings. It extends `evidence-pairing-v1`; it does not
alter the evidence estimand defined there. The first admitted analytic
specialization is the fixed-metric, equal-partition covariance of
crossvalidated Mahalanobis distances derived by Diedrichsen, Provost, and
Zareamoghaddam (2016).

This contract does not by itself introduce a public standard-error, interval,
or test API. Production support requires product--oracle equivalence,
generative calibration, and scale evidence in addition to the laws below.

## 1. Estimated evidence and its sampling law

For a relation estimate and an evidence query \(a=(H_a,K_a)\), write

\[
\widehat e_a
=
\widehat{\mathscr E}(H_a,K_a).
\]

The evidence pairing defines what was estimated. Its sampling covariance is a
separate object:

\[
\boxed{
\mathcal V_{ab}
=
\operatorname{Cov}
\left[
\widehat{\mathscr E}(H_a,K_a),
\widehat{\mathscr E}(H_b,K_b)
\right].
}
\]

For a finite family of queries, \(V=[\mathcal V_{ab}]\) is a covariance form
over estimated evidence. It is the sampling law of the evidence vector, not a
replacement for that vector and not another scientific estimand.

This yields the architectural sequence

```text
relation fit
  -> evidence plan
  -> estimated evidence
  -> sampling-covariance plan
  -> calibrated view
```

The evidence and calibration plans may share sources, frames, queries, and
compiled sufficient statistics. Their semantic identities remain distinct.

## 2. A second moment of the same observable

The evidence pairing is the irreducible completed observable. Sampling
covariance arises only after the relation becomes random through estimation.
Conceptually, it is a second moment on the tensor square of the evidence
transport:

\[
h_a=\operatorname{vec}(H_a),
\qquad
p_a=\operatorname{vec}\{\widehat{\mathcal P}(K_a)\},
\]

and therefore

\[
\boxed{
\operatorname{Cov}(\widehat e_a,\widehat e_b)
=
h_a^\top\operatorname{Cov}(p_a,p_b)h_b
=
\left\langle
h_ah_b^\top,
\operatorname{Cov}(p_a,p_b)
\right\rangle_F.
}
\]

This expression is explanatory. The compiler must not materialize a dense
fourth-order tensor or Kronecker lift by default. Admitted statistical models
must derive factored, structured, streamed, or direct-query evaluations.

The statistical completion is therefore:

> The evidence pairing defines the effect. The sampling-covariance transport
> defines uncertainty among estimated effects.

## 3. Signal, error information, and capabilities

A pure `relation` carries experimental--neural values and supports point
evidence. It does not acquire an error model from its matrix shape.

A `relation_fit` may additionally provide:

- an estimated relation;
- effect-coordinate covariance information;
- residual blocks or residual sufficient statistics;
- residual degrees of freedom;
- estimator and observation-model identities; and
- sampling, independence, and provenance declarations.

These fields form an error channel around the relation. They do not become the
mathematical identity of the relation itself. A precomputed relation may become
calibratable only when an independently valid, identity-bound external error
channel accompanies it. Betas alone cannot recreate residual information that
was discarded upstream.

Every calibration task records construction capabilities separately from
finite-sample diagnostics. Initial capabilities include:

```text
sampling_covariance: available | unavailable
error_channel: relation_fit | external | absent
error_model: separable_glm | declared_alternative
metric_role: metric | bridge
metric_status: fixed | learned
metric_uncertainty: not_applicable | ignored | propagated
partition_model: equal | heterogeneous
sampling_axis: time | trial | run | subject | ...
calibration_target: point_estimate | null | linear_hypothesis
spatial_covariance: unavailable | local_marginals | modeled
```

An unavailable capability must produce an explanatory refusal naming the
missing information and a viable entry path. It must not be reported as a
malformed-object error.

## 4. First admitted analytic specialization

Let \(K\) experimental conditions generate

\[
D=K(K-1)/2
\]

crossvalidated distance estimates, averaged over \(M\) independent
partitions. Under the equal-partition matrix-normal construction and a common
fixed neural metric, Diedrichsen et al. derive the joint sampling covariance

\[
\boxed{
V_{rs}
=
\frac{4}{M}
\Xi_{rs}
\frac{\mu_r\Sigma_R\mu_s^\top}{P^2}
+
\frac{2}{M(M-1)}
\Xi_{rs}^2
\frac{\operatorname{tr}(\Sigma_R^2)}{P^2}.
}
\]

Here:

- \(\mu_r\) is the whitened contrast pattern of distance \(r\), so that
  \(\Delta_{rs}=\mu_r\mu_s^\top/P\) contains the true distance cross-products;
- \(\Xi\) contains design/effect-coordinate covariance cross-products;
- \(\Sigma_R\) is residual covariance after the declared neural metric has
  been applied; and
- \(P\) is the distance normalization dimension under that convention.

The residual covariance enters the two terms differently. The
signal-independent term sees it only through the scalar
\(\operatorname{tr}(\Sigma_R^2)\). The signal-dependent term sees it as a
*metric on the whitened contrast patterns*: estimation error along a
high-variance residual axis is amplified in proportion to the signal lying on
that same axis. Replacing \(\mu_r\Sigma_R\mu_s^\top\) by a scalar multiple of
\(\mu_r\mu_s^\top\) is a different law. The two agree only in the isotropic
case \(\Sigma_R=\{\operatorname{tr}(\Sigma_R^2)/P\}I\), which includes the
common presentation in which the data are pre-whitened so that
\(\Sigma_R=I\) and distances are normalized by \(P\) features.

This is one \(D\times D\) law. Its diagonal contains variances of individual
distance estimates; its off-diagonal contains covariance between different
distances. The two terms are respectively signal-dependent and
signal-independent.

The estimator in their Eq. 10,

\[
\frac1M\sum_m
\widehat\delta_m^\top
\left(\frac1{M-1}\sum_{n\ne m}\widehat\delta_n\right),
\]

is algebraically the same equal-weight estimator as averaging every unordered
partition pair once:

\[
\boxed{
\frac{1}{\binom M2}
\sum_{m<n}
\widehat\delta_m^\top\widehat\delta_n.
}
\]

Consequently the fixed-metric, equal-partition `cross_partitions()` estimator
is within this specialization without an estimator conversion.

### Correction, 2026-08-15: the signal term is not isotropic

Between the first implementation and 2026-08-15 the compiler evaluated the
signal-dependent term with an isotropic surrogate for the whitened residual
covariance,

\[
\mu_r\Sigma_R\mu_s^\top
\;\longrightarrow\;
\frac{\operatorname{tr}(\Sigma_R^2)}{P}\,\mu_r\mu_s^\top,
\]

that is, it substituted \(\{\operatorname{tr}(\Sigma_R^2)/P\}I\) for
\(\Sigma_R\) inside the first term only. The signal-independent term was
always correct, and so was `target = "null"`, whose signal term vanishes.

The surrogate is exact when \(\Sigma_R\) is proportional to the identity
*with that same scale*, and wrong otherwise. Monte Carlo against the exact
sampling distribution of the estimator, with an AR(1) residual covariance
(\(\rho=0.75\)) over six features and five partitions, gives
surrogate-to-empirical variance ratios of

```text
whole_brain("none")  + identity metric    18.2  8.1  9.6
whole_brain("none")  + exact precision      5.9  5.8  5.7   (= P)
whole_brain("local") + identity metric      3.3  1.4  1.6
whole_brain("local") + exact precision      1.0  1.0  1.0   (coincidence)
```

The last row is the configuration in which frame normalization makes
\(\Sigma_R=I/P\) and the surrogate scale matches, and it was also the only
configuration the product-oracle test exercised. The defect was therefore
invisible to a test suite in which both the implementation and its oracle
carried the same substitution. The corrected law, and the surrogate's failure
against Monte Carlo in the three anisotropic configurations, are both pinned
in `tests/testthat/test-evidence-sampling-nonspherical.R`.

Two consequences are normative:

1. An oracle that shares an algebraic shortcut with the implementation is not
   an independent court. Every law in the sampling court must be checked
   against either the sampling distribution of the estimator itself or a
   scalar-loop transcription of the equation, not only against a second
   expression of the same convenience.
2. A statistical law must be validated in a configuration where its terms are
   distinguishable. A single spherical fixture cannot separate
   \(\mu_r\Sigma_R\mu_s^\top\) from \(\mu_r\mu_s^\top\).

Reference: Diedrichsen, J., Provost, S., and Zareamoghaddam, H. (2016),
[*On the distribution of cross-validated Mahalanobis
distances*](https://arxiv.org/abs/1607.01371), especially Eqs. 10, 13, and 35
and Sections 3.4 and 5.1.

### Correction, 2026-08-16: the noise term is a quadratic functional of an estimate

The law above is written in terms of the true \(\Sigma_R\). The compiler does
not have it. It substitutes the partition-pooled sample residual covariance
\(S_R\) with

\[
\nu=\sum_m\operatorname{df}_m
\]

degrees of freedom. That substitution is *not* neutral, because the
signal-independent term depends on \(\Sigma_R\) through the quadratic
functional \(\operatorname{tr}(\Sigma_R^2)\). For
\(S_R\sim W_P(\nu,\Sigma_R)/\nu\),

\[
E\operatorname{tr}(S_R^2)
=
\frac{\nu+1}{\nu}\operatorname{tr}(\Sigma_R^2)
+
\frac{\operatorname{tr}(\Sigma_R)^2}{\nu},
\]

so a plug-in \(\operatorname{tr}(S_R^2)\) overstates the functional, and the
reported *standard error* by

\[
\sqrt{1+\frac{1+P_{\mathrm{eff}}}{\nu}},
\qquad
P_{\mathrm{eff}}
=
\frac{\operatorname{tr}(\Sigma_R)^2}{\operatorname{tr}(\Sigma_R^2)} .
\]

\(P_{\mathrm{eff}}\) is the participation ratio of the whitened residual
covariance: the number of residual directions the support actually spends
variance on, which equals the support size only for spherical \(\Sigma_R\).
The consequence is not a rounding artifact. At \(\nu=168\) (six runs of a
32-trial, four-condition design), the reported standard error was too large by
a factor of

```text
P_eff =   1  (voxelwise)     1.006
P_eff =  50  (radius ~3 SL)  1.142
P_eff = 120                  1.311
P_eff = 800  (large ROI)     2.40
```

Voxelwise frames were therefore never materially wrong; searchlight and
region frames were, in the conservative direction.

The corrected estimator is the Wishart-unbiased form of the same functional,

\[
\boxed{
\widehat{\operatorname{tr}}(\Sigma_R^2)
=
\frac{\nu^2}{(\nu-1)(\nu+2)}
\left(
\operatorname{tr}(S_R^2)
-
\frac{\operatorname{tr}(S_R)^2}{\nu}
\right),
}
\]

clamped at zero, and it is applied **only** where the residual covariance is a
plug-in. Three points are normative:

1. The signal term \(\mu_r\Sigma_R\mu_s^\top\) is *linear* in \(\Sigma_R\) and
   receives no correction. Plugging an unbiased \(S_R\) into a linear
   functional is unbiased; plugging it into a quadratic one is not. Any future
   term of the law must be classified this way before an estimate is
   substituted into it.
2. A caller that supplies the true \(\Sigma_R\) — an oracle, a known-covariance
   test, a simulation with a declared truth — must *not* receive the
   correction. The distinction is carried explicitly by a `residual_df`
   argument whose `NULL` default means "this is the covariance itself".
3. \(\nu\) buys information about at most \(\nu\) residual directions. When
   \(\nu<P_{\mathrm{eff}}\) the corrected estimator's own sampling error is of
   the order of the quantity and its clamp at zero converts an unusable
   estimate into a confidently small standard error. The compiler refuses with
   capability `sufficient_residual_df` rather than reporting one, and both
   \(\nu\) and \(P_{\mathrm{eff}}\) are reported on every result.

The defect survived the court in §11 because the court's own oracles
recomputed \(\operatorname{tr}(S_R^2)\) from the same plug-in, and because the
Monte Carlo tests supplied the *true* \(\Sigma_R\) to both sides. That is the
same failure mode as the 2026-08-15 correction in a different disguise: an
oracle sharing an assumption with the implementation. §11 now requires an
end-to-end Monte Carlo test in which \(\Sigma_R\) is re-estimated from
residuals on every replication.

Reference for the estimator: Srivastava, M. S. (2005), "Some tests concerning
the covariance matrix in high dimensional data", *Journal of the Japan
Statistical Society* 35(2), 251–272.

### Independent observations within a partition

Both the effect covariance \(\Xi\) and the residual degrees of freedom \(\nu\)
are computed as if the observations within a partition were independent given
the design. fMRI residuals are not. Without an observation whitener the
reported standard error can err in either direction, and by a large factor.
Monte Carlo under AR(1) errors (\(\rho=0.75\), 32 trials, four conditions, six
runs, 50 features, corrected noise term), reporting the ratio of the true
spread of the estimator to the reported standard error:

```text
randomly interleaved order              0.50
blocked order                           5.10
blocked order, correct whitener         1.03
randomly interleaved, correct whitener  1.00
```

Supplying \(L\) with \(L^\top L=\Sigma_t^{-1}\) to `lm_relation_fit()`
restores calibration, because the whitened problem satisfies the assumption.
crossform applies and records the \(L\) it is given; it cannot verify that
\(L\) matches the data's autocorrelation. This is a documented assumption of
the specialization, not a defect, but it must be stated wherever the law is
offered.

## 5. Partition endpoints are independent; pair products are not

`cross_partitions()` declares that distinct partition estimates may be used in
unbiased cross-products. It does not declare that the resulting edge products
are independent observations.

Each partition participates in multiple unordered pairs. Edge products that
share an endpoint therefore share estimation noise. A reducer weight or edge
count has no authority to manufacture sampling independence.

The governing rule is:

> Generalization edges state where evidence must reproduce. They are
> computational contributions to an estimator, not sampling replicates for an
> ordinary spread-across-edges standard error.

In particular, `sd(edge_values) / sqrt(number_of_edges)` is not an admitted
standard-error estimator for the all-pairs mean. Its error is signal-dependent
and can worsen as the number of partitions increases.

For the actual sample-SD construction, the exact signal-dominated
underestimation factor is \(\sqrt{M+1}\), not \(\sqrt{M-1}\). The latter arises
under a different comparison that treats marginal edge variance as known and
divides it by the nominal pair count. Documentation and tests must name which
naive estimator is being analyzed.

## 6. Metric scope and uncertainty

The first analytic specialization requires a fixed common metric. This
includes Euclidean geometry \(K=I\) and a declared fixed precision metric.

A metric learned from residual data is random. Conditioning on its realized
value can validate a point-evidence evaluator, but an interval or test that
ignores its estimation uncertainty generally does not have the same coverage
claim. The task must record

```text
metric_status: learned
metric_uncertainty: propagated | ignored
```

and no calibrated-interval capability may be granted when propagation is
required but absent. Metric construction and calibration remain distinct
stages while sharing the same uncertainty source. LD-t and related procedures
are possible later specializations; they are not implied by this contract.

An edge- or location-indexed schedule \(K_{x,e}\) additionally changes the
sampling law because different evidence contributions can use different
estimated operators. Diedrichsen Eq. 13 must not be applied unchanged to that
case without a derivation establishing equivalence or accounting for operator
randomness and dependence.

## 7. Unequal partitions and bridges

The separable equal-partition GLM model is a declared first capability, not the
foundation of the architecture. When effect-coordinate or residual covariance
differs by partition, the compiler requires a heterogeneous-partition sampling
law, such as the appropriate general expressions in Diedrichsen et al., rather
than silently substituting the equal-partition formula.

For distinct neural spaces, \(K\) is a bridge rather than a self-space metric.
Its sampling law must state whether the bridge is fixed or estimated and how
uncertainty is transported from both sides. A PSD metric assumption cannot be
inferred for a rectangular bridge.

## 8. Local calibration is not spatial calibration

For a spatial frame, the first implementation may produce one experimental
sampling-covariance form \(V_x\) per measurement \(x\). The local support,
normalization, metric, residual covariance, and effective dimension are part of
that task identity.

This does not establish

\[
\operatorname{Cov}(\widehat e_x,\widehat e_y)
\]

between locations. Overlapping searchlights share features and residual noise,
so spatial multiplicity or field inference requires an additional spatial
sampling model or valid resampling construction. Local Eq. 13 variances must
not be relabelled as a calibrated whole-map random field.

Voxel, ROI, searchlight, and whole-brain uncertainty still use one frame-aware
calibration architecture; they differ in the measurement and residual
sufficient statistics supplied to it, not through method-specific inference
engines.

## 9. Calibration targets

Let \(\widehat d\) be an evidence vector with sampling covariance \(V\).
Admitted exact linear operations include:

\[
\operatorname{Var}(a^\top\widehat d)=a^\top Va
\]

and covariance transport \(AVA^\top\) for a declared linear view \(A\).

A variance estimate must state the parameter value or null under which
signal-dependent terms are evaluated. Replacing an unknown nonnegative true
distance by \(\max(\widehat d,0)\) is a named plug-in convention. It is not the
same as evaluating the covariance under the null \(d=0\), and it must not be
presented as the unique interpretation of Diedrichsen et al.'s Section 5.1.

Two distinct plug-ins occur in this law and they must not be conflated:

- the **signal** plug-in, which is a *calibration target* the caller chooses
  (`target = "plugin"` or `"null"`); and
- the **residual-covariance** plug-in, which is unavoidable — \(\Sigma_R\) is
  never known — and is therefore not a choice but a correction obligation. Its
  quadratic contribution carries the finite-sample correction of §4, so
  `target = "null"` is exact on the variance scale rather than merely
  centered: the signal term vanishes and the surviving term is an unbiased
  estimate of the true variance, not an upper bound on it.

Every result must report \(\nu\) and \(P_{\mathrm{eff}}\) so a reader can see
how much residual information stands behind the second term, and must refuse
when \(\nu<P_{\mathrm{eff}}\).

The plug-in policy implemented here, `partition_mean_plugin`, substitutes the
partition mean of the *estimates* \(\bar B\) for the unknown signal. Since
\(\operatorname{Cov}(\widehat\mu_r,\widehat\mu_s)=\Xi_{rs}\Sigma_R/M\),

\[
E\!\left[\widehat\mu_r\Sigma_R\widehat\mu_s^\top\right]
=
\mu_r\Sigma_R\mu_s^\top
+
\frac{\Xi_{rs}}{M}\operatorname{tr}(\Sigma_R^2),
\]

so the plug-in covariance is biased upward by
\(4\Xi_{rs}^2\operatorname{tr}(\Sigma_R^2)/(M^2P^2)\). This is a policy, not a
defect: the inflation is \(O(M^{-2})\), is largest when noise dominates the
true distances, and makes plug-in intervals mildly conservative rather than
anticonservative. It must be disclosed wherever the policy is offered, and
`target = "null"` remains the exact choice for calibrating a test of no
effect. Since 2026-08-16 this is the *only* remaining bias of the plug-in
policy; the separate noise-trace inflation that used to accompany it is
corrected in §4.

An analytic z-test, confidence interval, contrast test, LD-t, bootstrap,
permutation distribution, and population resampling are different calibration
procedures. They may reuse \(V\) but require their own assumptions, targets,
and validation evidence.

Subject-level resampling estimates population uncertainty. It is not a way to
recover discarded within-participant residual information for a single
precomputed-beta relation.

## 10. Queryable and exact execution

For \(D\) distances, dense \(V\) contains \(D^2\) entries. Retaining the
sampling law therefore means preserving exact queryability, not eagerly
allocating every covariance.

An `evidence_sampling_plan` should admit operations such as

```text
diagonal()
selected_entries(index_pairs)
apply(a)
quadratic_form(a)
transport(A)
materialize()
```

The compiler may exploit symmetry, Hadamard structure, low-rank design terms,
sparsity, direct linear transport, and matrix-free action. `materialize()` is
explicit and size-preflighted.

The performance law is the same as for evidence geometry:

> Remove algebraic and administrative work while preserving the exact declared
> estimator. Approximation, stochastic trace evaluation, or reduced precision
> creates a separately named estimator with an explicit error contract.

## 11. Required executable law court

Before an analytic sampling specialization becomes public, independent tests
must cover:

1. equivalence of the Eq. 10 and all-unordered-pairs estimators;
2. Eq. 13 against an independent direct or closed-form oracle, including a
   scalar-loop transcription that shares no algebraic shortcut with the
   implementation;
3. Eq. 13 against the Monte Carlo sampling distribution of the estimator
   itself, under a non-spherical whitened residual covariance and at more
   than one frame normalization and metric, so that the signal term's
   dependence on \(\Sigma_R\) is identifiable;
4. the complete \(D\times D\) covariance, not only its diagonal;
5. exported product-path agreement with an oracle importing no package
   internals;
6. null centering and signal-dependent variance in generative simulations;
7. coverage under the declared plug-in or null policy, and the stated
   direction and size of any plug-in bias;
8. refusal when the error channel is absent;
9. refusal or qualification for learned metric uncertainty;
10. unequal-partition negative fixtures;
11. direct-query versus full-materialization equivalence;
12. runtime and peak-memory gates at realistic \(q\), frame, and feature
    sizes; and
13. **an end-to-end Monte Carlo test of the exported path in which the
    residual covariance is re-estimated from residuals on every replication**
    (added 2026-08-16). The fit must be recomputed per replication, not
    reused; the configuration must reach \(P_{\mathrm{eff}}/\nu\ge 0.3\), the
    regime in which a plug-in quadratic functional is visibly biased; the
    acceptance band must be pre-registered in the test file in units of the
    Monte Carlo standard error of a ratio of standard deviations,
    \(1/\sqrt{2(N-1)}\); and the test must assert both that the corrected
    estimator falls inside that band **and** that the uncorrected one falls
    outside it, so the correction is pinned rather than merely present.

    Requirements 2 and 3 are not substitutes for 13. An oracle that is handed
    the true \(\Sigma_R\), or that recomputes \(\operatorname{tr}(S_R^2)\)
    from the same plug-in the implementation uses, cannot see a bias that
    lives in the substitution itself. Where a test oracle is handed a plug-in
    residual covariance it must apply the §4 correction independently, and
    where it is handed the truth it must not.

The court is implemented in `tests/testthat/test-evidence-sampling-laws.R`
(1, 2, 4), `-nonspherical.R` (2, 3, 7), `-generative.R` (6, 7),
`-product.R` (5, 8, 9, 11), `-df-correction.R` (13), `-scale.R` (12), and
`test-capability-refusals.R` (10).

Verification claims use four levels:

```text
algebraically established
numerically verified
statistically validated
scale-qualified
```

Reaching one level does not imply the next.

## 12. Initial scope boundary

`evidence-sampling-v1` authorizes implementation planning for fixed-metric,
equal-partition analytic covariance and exact linear transports. It does not
authorize:

- treating partition-pair rows as independent samples;
- calibrating a pure precomputed relation without an error channel;
- silently applying the equal-partition formula to heterogeneous partitions;
- ignoring learned-metric randomness while claiming calibrated coverage;
- inferring cross-location covariance from local variances;
- silently choosing a plug-in distance or null;
- presenting subject resampling as recovered first-level uncertainty; or
- claiming a calibrated public API before product, generative, and scale gates
  pass.

Those boundaries are part of the statistical semantics, not temporary
documentation caveats.
