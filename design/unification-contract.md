# Cross-generalized geometry unification contract

Status: normative scientific contract

Contract version: `unification-v1`

Date: 2026-08-20

This document fixes the claim and vocabulary used by the population-estimation
and interpretive-validation program. It is normative for PE-A and for every
later program artifact that describes what `crossform` unifies. The executable
court is `design/oracles/unification-vocabulary.R`; it uses only base matrix
algebra. `tests/testthat/test-unification-contract.R` compares that independent
court with the public package route.

## 1. The claim

`crossform` provides one typed, cross-generalized experimental-neural geometry
from which **fixed linear queries of first and second moments** can be read.
The geometry keeps the experimental axes, spatial frame and scale, neural
metric, partition pairing and reduction, and—at population level—the realized
transport in the estimand rather than treating them as interchangeable
execution details.

This is a claim of unified estimands and interpretation. It is not a claim that
all representational statistics are linear queries, and it is not a claim to
unify preprocessing, effect estimation, transport learning, or every
downstream inferential procedure.

For partition-specific representations

\[
B_r:\mathcal N\longrightarrow\mathcal E,
\qquad B_r\in\mathbb R^{q\times p},
\]

a declared ordered pairing/reducer \(\Gamma\), and a fixed measurement metric
\(K_x\in\mathbb R^{p\times p}\), the common second-moment geometry at node
\(x\) is

\[
G_x^{\mathrm{total}}
=\sum_{r,s}\Gamma_{rs}B_rK_xB_s^\top.
\]

When every nonzero edge joins relation estimates with uncorrelated estimation
errors, those estimates are mean-zero around one shared target, the weights
sum to one, and \(K_x\) is fixed relative to the paired data, this is a
noise-unbiased estimate of the corresponding cross-generalized geometry.
`independence = "independent"` records that assumption; it does not establish
it. Negative crossvalidated estimates are valid and are not clipped.

The complete estimand is therefore not merely the numeric matrix \(G_x\). It
is the tuple

\[
\mathcal T=(\mathcal E_L,\mathcal E_R,\mathcal N,
            \mathcal F,\mathcal S,K,\Gamma,R,P),
\]

where \(\mathcal F\) is the frame, \(\mathcal S\) its scale/family identity,
\(K\) the metric schedule, \(\Gamma\) the partition pairing and reducer,
\(R\) the declared result component/query, and \(P\) the realized transport
when the result crosses subjects. A change to any member changes the scientific
question and must change the plan identity.

## 2. Normative vocabulary

Each term below has one definition, one operational reading, declared units or
normalization, and an executable witness. The HTML comments are machine-read
by the contract test; removing a term breaks the test rather than silently
changing the public vocabulary.

<!-- contract-term: representation -->
### 2.1 Representation

**Mathematical definition.** A representation is a typed map \(B_r\) from a
named neural feature space to an ordered experimental effect space for one
partition \(r\). In the implemented orientation it is a \(q\times p\) matrix:
rows are effect coordinates and columns are neural features. A study supplies
a partition-indexed family \(\{B_r\}\), not one unqualified pattern matrix.

**Operational interpretation.** It is the effect-by-feature content carried by
an `effect_relation`, usually constructed by `relation()` from fitted effects
or by `plan_relation()` followed by `estimate_relation()` from observations.
`effect_space()` and the neural `effect_domain` make equal-looking but
differently identified axes non-interchangeable.

**Units and normalization.** Entries retain the fitted-effect units declared
by the relation. Crossform does not silently standardize them. Partition labels,
effect order, feature identity, units, and provenance belong to the object.

**Executable witness.** `unification_oracle$representation` contains two
named `2 x 4` maps; the package comparison constructs the same maps with
`relation()`.

<!-- contract-term: total_geometry -->
### 2.2 Total geometry

**Mathematical definition.** The total geometry at node \(x\) is the full
cross-partition second moment

\[
G_x^{\mathrm{total}}=\sum_{r,s}\Gamma_{rs}B_rK_xB_s^\top.
\]

“Common geometry” names the typed family of these forms. “Total” names its
undecomposed component at a particular frame node; it does not mean a sum over
nodes. For a symmetric self form it is stored in Frobenius-isometric `svec`
coordinates, but packing is a physical codec, not another estimand.

**Operational interpretation.** `plan_geometry()` binds relation, frame,
metric, and pairing before values are read. `materialize_geometry()` stores the
complete form; `geometry_component(x, "total")` reads it; named or custom
queries may instead use `evaluate_geometry()` without materializing it.

**Units and normalization.** For self forms, units are squared effect units
times the declared neural-metric units. Built-in frame weights and family
weights are dimensionless. `normalization = "local"`, `"none"`, and
`"conservative"` define different node quantities and remain in plan identity.

**Executable witness.** `unification_oracle$geometry[[i]]$total` constructs the
form directly for every family row; `unification_oracle$whole$total` is the
whole-domain comparator.

<!-- contract-term: coherent_component -->
### 2.3 Coherent component

**Mathematical definition.** The coherent component is the part of total
geometry carried by the measurement's weighted common neural mode. For the
canonical additive diagonal node with nonnegative weight vector \(w_x\) and
mass \(a_x=\mathbf1^\top w_x\),

\[
G_x^{\mathrm{coherent}}
=\sum_{r,s}\Gamma_{rs}
  \frac{(B_rw_x)(B_sw_x)^\top}{a_x}.
\]

The metric-general implementation obtains the corresponding fixed coherent
operator from the declared measurement/metric decomposition; it does not
replace a dense metric with the diagonal formula above.

**Operational interpretation.** It is reproducible energy carried by a node's
own weighted mean pattern—the univariate/common-mode part of the same
geometry, not a second fit. Read it with
`geometry_component(x, "coherent")`, `query_geometry(..., "coherent")`, or
the `$coherent` column of `contrast_energy()`.

**Units and normalization.** It has exactly the units and node normalization
of total geometry. It is frame-relative: coherent components do not generally
sum to the whole-domain coherent component.

**Executable witness.** `unification_oracle$geometry[[i]]$coherent` evaluates
the formula above, including an alternating node whose weighted mean cancels.

<!-- contract-term: configurational_component -->
### 2.4 Configurational component

**Mathematical definition.** The canonical noun used by the API is
**configuration**; “configurational” is its adjective. The component is the
exact complement

\[
G_x^{\mathrm{configuration}}
=G_x^{\mathrm{total}}-G_x^{\mathrm{coherent}}.
\]

For the additive diagonal case its neural operator is
\(D(w_x)-w_xw_x^\top/a_x\), the weighted spatial variation orthogonal to the
common mode.

**Operational interpretation.** It is reproducible spatial pattern beyond the
node's weighted mean. It is not “noise,” a demeaned reanalysis, or the residual
of a fitted population model. Read it with
`geometry_component(x, "configuration")`,
`query_geometry(..., "configuration")`, or `$configuration` from
`contrast_energy()`.

**Units and normalization.** It has the same units and normalization as total
and coherent geometry. The identity `coherent + configuration = total` holds
coefficient by coefficient and survives every fixed linear query. Signed
crossvalidated components need not be nonnegative.

**Executable witness.** Every entry of
`unification_oracle$geometry[[i]]$configuration` is constructed as the exact
complement, and the oracle fails if recomposition exceeds `1e-12`.

<!-- contract-term: frame_family -->
### 2.5 Frame family

**Mathematical definition.** A frame family is an ordered set of node-weight
operators \(W_s\) over one neural domain, with each member conservative on its
own, combined as

\[
W_{\mathrm{family}}=
\begin{bmatrix}\alpha_1W_1\\ \cdots\\ \alpha_SW_S\end{bmatrix},
\qquad \alpha_s>0,\quad\sum_s\alpha_s=1.
\]

Each \(W_s\) has unit column mass. Under the feature-additive identity or
diagonal metric used by the executable witness, the family conserves overall
and each scale carries exactly \(\alpha_s\) of the whole-domain **total**
geometry. A dense spatial metric mixes features, so that scale-budget law does
not apply; the metric schedule remains part of the estimand. The law never
extends to coherent components.

**Operational interpretation.** `frame_family()` combines already compiled
conservative frames; multiscale `searchlights()` is shorthand for the same
construction. Each row retains `measurement`, `family`, `node`, `scale`,
`center`, and `alpha` metadata. `frame_conservation()` checks the total-budget
law.

**Units and normalization.** `alpha` is dimensionless and is never silently
renormalized. Every member and the stack use `normalization = "conservative"`.
Changing a scale, member, α, domain, or normalization changes the estimand.

**Executable witness.** `unification_oracle$frame_family` stacks two four-node
searchlight frames at radii 0.5 and 1.01, with α = 0.5 each, and fails unless
every feature's stacked mass is exactly one.

<!-- contract-term: scale_profile -->
### 2.6 Scale profile

**Mathematical definition.** For fixed query \(H\), a scale profile is the
ordered mapping

\[
s\longmapsto(C_s,Q_s,T_s,\phi_s),\qquad
C_s=\sum_{x\in s}\langle H,G_x^{\mathrm{coherent}}\rangle_F,
\]

\[
Q_s=\sum_{x\in s}\langle H,G_x^{\mathrm{configuration}}\rangle_F,
\quad T_s=C_s+Q_s,
\quad \phi_s=C_s/T_s
\]

where \(\phi_s\) is reported only when the aggregated components form a finite,
positive, nonnegative partition. Under a conservative family and a
feature-additive identity or diagonal metric,
\(T_s=\alpha_s\langle H,G_\Omega^{\mathrm{total}}\rangle_F\); total energy by
scale therefore pictures α, while the coherent share is α-invariant and is the
interpretable scale-resolved quantity.

**Operational interpretation.** `coherence_spectrum()` returns this profile
from a geometry plan or contrast view. With `by_location = TRUE`, the profile
at one location is a filter over the returned `(scale, center)` table; no
undeclared collapse over location or scale is performed.

**Units and normalization.** \(C_s,Q_s,T_s\) have query-output units;
\(\phi_s\) and α are dimensionless. The result is explicitly tied to its frame
family and fixed query.

**Executable witness.** `unification_oracle$scale_profile` aggregates both
scales under the identity metric, verifies `T_s = alpha_s * whole_query`, and
recovers different coherent shares for the fine and coarse frames.

<!-- contract-term: query -->
### 2.7 Query

**Mathematical definition.** A second-moment query is a finite matrix \(H\)
fixed independently of the local geometry, evaluated as the Frobenius
functional

\[
q_H(G)=\langle H,G\rangle_F=\operatorname{tr}(H^\top G).
\]

A first-moment query is a fixed vector \(c\) applied to the
pairing-appropriate local marginal. For an undirected pairing,

\[
m_x(c)=c^\top\left\{\frac12\sum_e\omega_e
\frac{(B_{l(e)}+B_{r(e)})w_x}{a_x}\right\}.
\]

`H = cc^T` gives contrast energy; `H = (e_i-e_j)(e_i-e_j)^T` gives one RDM
entry; a fixed linear map of RDM entries gives fixed linear RSA. The detailed
equivalence assumptions and proofs belong to PE-B, not to this vocabulary
contract.

**Operational interpretation.** `bilinear_query()` describes symmetric
self-form operators, `pair_query()` binds rectangular operators to two effect
spaces, and `evaluate_geometry()`/`query_geometry()` execute them. Named
constructors `contrast_energy()`, `rdm()`, `rsa()`, and `crossnobis()` compile
common queries. The query, component, effect axes, and parent plan identity are
retained in the view identity.

**Units and normalization.** A dimensionless \(H\) preserves geometry units;
a query with scientific coefficients carries their units into the result.
RDM and RSA normalizations are part of their declared operators, not generic
properties of geometry.

**Executable witness.** `unification_oracle$query` fixes \(c=(1,-1)\) and
`H = tcrossprod(c)` before evaluation. `first_query` and `second_query` verify
the first- and second-moment formulas independently.

## 3. Public object and function map

| Contract concept | Public objects and constructors | Public readout or audit |
|---|---|---|
| representation | `effect_space()`, `relation()`, `plan_relation()`, `estimate_relation()` | relation print/identity and `source_capabilities()` |
| total geometry | `plan_geometry()`, `materialize_geometry()` | `geometry_component(..., "total")`, `evaluate_geometry()` |
| coherent component | the fixed measurement/metric decomposition inside a geometry plan | `geometry_component(..., "coherent")`, `contrast_energy()` |
| configuration component | exact complement of the same plan | `geometry_component(..., "configuration")`, `contrast_energy()` |
| frame family | `frame_family()`, multiscale `searchlights()`, `compile_frame()` | `frame_conservation()` and the frame `$index` |
| scale profile | a conservative frame-family plan plus a fixed contrast | `coherence_spectrum()` |
| query | `bilinear_query()`, `pair_query()` and named query constructors | `evaluate_geometry()`, `query_geometry()`, `contrast_energy()`, `rdm()`, `rsa()`, `crossnobis()` |
| population transport | `location_transport()`, `anatomical_transport()`, `external_transport()` | `transport_values()`, population receipt diagnostics |
| population estimand | `plan_population()` | `estimate_population()`, `materialize_population()` and population views |

The transport row is part of the unification claim even though transport is
not one of the seven new vocabulary terms. At population level, each subject's
native geometry ledger is carried by a declared `location_transport`. Its
operator content, semantics, row mass, group index, provenance, and signature
enter `plan_population()` identity. Crossform consumes that realized
transport; it does not learn registration or alignment.

## 4. One executable route

Run the independent court from the repository root:

```sh
Rscript design/oracles/unification-vocabulary.R
```

The package court builds the same route through public APIs:

```r
domain <- abstract_domain(
  4L, coordinates = cbind(0:3, 0), feature_ids = paste0("v", 1:4),
  id = "unification-v1"
)
B <- list(
  run1 = rbind(a = c(2, -2, 1, 1), b = rep(0, 4)),
  run2 = rbind(a = c(1, -3, 1, 1), b = rep(0, 4))
)
relation <- relation(B, domain = domain)
family <- compile_frame(
  searchlights(
    c(0.5, 1.01), "conservative", weights = c(0.5, 0.5)
  ),
  domain
)
plan <- plan_geometry(
  relation, family,
  cross_partitions(
    relation, independence = "independent", generalizes_over = "run"
  )
)
contrast <- c(a = 1, b = -1)
effect <- contrast_energy(plan, contrast)
profile <- coherence_spectrum(plan, contrast)
geometry <- materialize_geometry(plan)
stopifnot(
  isTRUE(all.equal(
    geometry_component(geometry, "total"),
    geometry_component(geometry, "coherent") +
      geometry_component(geometry, "configuration")
  )),
  isTRUE(all.equal(effect$total, effect$coherent + effect$configuration))
)
```

This route witnesses all seven terms in one small example. It does not provide
empirical evidence: the matrices are a hand-computable algebra fixture.

## 5. Boundaries and evidence labels

The unification claim covers only fixed linear first- and second-moment
readouts of the declared geometry. In particular:

- nonlinear rank correlations, locally data-adaptive queries, classification,
  feature selection, and learned prediction are not special-case fixed linear
  geometry queries;
- preprocessing, registration, masking, universal HRF modeling, effect
  extraction, and transport learning remain upstream;
- a point-estimate theorem is not an uncertainty or coverage theorem;
- coherent/configuration arithmetic is exact even when either crossvalidated
  estimate is negative, but a coherence fraction is defined only for a
  nonnegative partition; and
- a prospective protocol is not a completed empirical result or an
  independent replication.

Claims governed by this contract must use one or more of the eight evidence
classes defined in `design/evidence-status-ledger.md`: `algebraic_theorem`,
`internal_oracle`, `external_parity`, `matched_simulation`,
`existing_illustration`, `prospective_protocol`,
`completed_real_data_result`, and `independent_replication`. The versioned
machine-readable rows live in
`inst/extdata/certification/evidence-status-ledger.csv`. A prospective protocol
cannot be called demonstrated or replicated, and neither the Haxby nor
ds003745 regression artifact qualifies as a completed prospective result.

## 6. Implementation audit and change control

The contract was checked against the current public implementation on
2026-08-20:

- `R/pairing.R` and `R/execution-driver.R` implement pairing-appropriate first
  moments; undirected pairings retain the orientation-invariant endpoint
  marginal rather than an arbitrary left/right orientation;
- `R/measurement-decomposition.R`, `R/result.R`, and `R/views.R` expose the
  exact total/coherent/configuration split and mask, rather than clamp,
  invalid coherence fractions;
- `R/scope.R`, `R/query-structured.R`, and `R/geometry-entry.R` admit fixed
  pair/bilinear queries and keep nonlinear or adaptive readouts outside this
  lowering;
- `R/frame.R` and `R/views.R` enforce conservative member frames, α weights,
  per-row family metadata, and scale-profile semantics; and
- `R/population-plan.R` includes subject plan identities and transport
  signatures in the population scientific identity.

`tests/testthat/test-unification-contract.R` is the executable review gate. It
fails if a vocabulary marker or required operational mapping disappears, if
the independent formulas cease to match public output, if component
recomposition or scale conservation fails, or if the boundary excluding all
representational statistics is removed. Changes to these definitions require
a new contract version and coordinated updates to public prose; changing an R
class label or physical codec alone does not change the mathematical term.
