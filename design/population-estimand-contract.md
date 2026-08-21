# Population estimand, coverage, and realized-transport contract

Status: normative scientific contract

Contract version: `population-estimand-v1`

Date: 2026-08-20

This contract fixes what a population coefficient means when subjects have
different native frames, realized transports, or group-node coverage. It
refines `design/population-form-contract.md` without changing that contract's
transport algebra. PE-D1 implements the node-specific subject-set record
specified here; later uncertainty work must consume that record rather than
reconstructing availability from fitted values.

The independent finite-sample court is
`design/oracles/population-estimand-targets.R`. It uses only base R.
`tests/testthat/test-population-estimand-contract.R` pins this prose and checks
the current package's plan- and receipt-level provenance.

## 1. Conditioning statement

All population point estimates and uncertainty procedures are conditional on
the following realized objects:

\[
\mathcal C=(I_0,X,\{\mathcal G_i\},\{P_i\},\{\mu_i\},J,
            \mathcal A,\mathcal N,\mathcal Q),
\]

where \(I_0\) is the planned subject set, \(X\) its fixed design matrix,
\(\mathcal G_i\) the subject geometry plans, \(P_i\) the realized transport
operators including their sink columns, \(\mu_i\) their declared row masses,
\(J\) the common group-node index, \(\mathcal A\) the realized availability
sets, \(\mathcal N\) the population normalization, and \(\mathcal Q\) the
fixed query or complete-form coordinate.

The procedures do **not** propagate uncertainty from estimating registration,
alignment, transport, frame membership, coverage, the query, or a learned
subject metric. A functional transport must carry cross-fit provenance, but
cross-fitting prevents circular evaluation; it does not turn the realized
transport into a fixed population constant. Any future procedure that
propagates transport or selection uncertainty is a different inferential
contract and must receive a different scientific identity.

The subject population remains the one named by the sampling protocol. A
finite-sample coefficient is an estimate over the planned or available
participants below; interpreting it as a superpopulation coefficient further
requires the declared participant-sampling assumptions. Transport does not
convert a convenience sample into a probability sample.

## 2. Response and availability

For subject \(i\), ordinary group node \(j\), and query/coordinate \(h\), let

\[
y_{ijh}=\mathcal N_{ih}\left[
  \sum_x P_{i,xj}\langle H_h,G_{i,x}\rangle_F
\right]
\]

for budget semantics, with the declared transported-mass divisor added under
density semantics. Let

\[
m_{ij}=\sum_xP_{i,xj}\mu_{i,x}
\]

be the subject's transported native territory at node \(j\). The ordinary-node
availability indicator is

\[
A_{ijh}=1\{m_{ij}>0\}\,1\{y_{ijh}\text{ is finite}\}\,
          1\{\mathcal N_{ih}\text{ is admitted}\}.
\]

The implementation applies its declared numerical tolerance to `m > 0`; it
must not infer availability from whether `y` happens to equal zero.

This distinction is load-bearing. Under budget semantics, a node with no
transported mass numerically receives zero; under density semantics it receives
`NA`. Both encode the same absence and therefore produce the same \(A_{ijh}\).
A true observed zero with positive mass remains available.

The sink is different. Its value is lost budget, and zero means the observed
fact that no budget was lost. Every planned subject is therefore available for
the sink fit when its native query and normalization are admitted. Positive
sink mass is reported as a diagnostic; it is not the sink's availability rule.

## 3. The two targets

<!-- population-target: all_planned -->
### 3.1 All-planned-subject target

The default policy is `coverage_policy = "all_planned"`. For cell \((j,h)\),
the finite-sample target is the OLS projection over the exact planned set
\(I_0\):

\[
\widehat\beta^{\mathrm{planned}}_{jh}
=\arg\min_b\sum_{i\in I_0}(y_{ijh}-X_i b)^2.
\]

It is returned only when every planned subject is available and the requested
coefficient or contrast is estimable. If any \(A_{ijh}=0\), the coefficient is
unresolved with reason `planned_subject_unavailable`; zero filling, implicit
case deletion, and mean imputation are forbidden. This policy keeps one subject
population across nodes at the cost of honest missing results where the
realized transport does not cover everyone.

The corresponding superpopulation parameter solves
\(E[X(Y_{jh}-X^\top\beta)\mid\mathcal C]=0\) for the protocol's subject
population. The conditioning includes the realized transports; it does not
license a claim about transport-averaged anatomy.

<!-- population-target: available_at_node -->
### 3.2 Available-at-node target

The explicit alternative is `coverage_policy = "available_at_node"`. Define
\(I_{jh}=\{i\in I_0:A_{ijh}=1\}\). Its finite-sample target is

\[
\widehat\beta^{\mathrm{available}}_{jh}
=\arg\min_b\sum_{i\in I_{jh}}(y_{ijh}-X_i b)^2.
\]

This target may use a different subject set at every node and query. It is not
an approximation to the all-planned target and must never inherit its label or
scientific identity. Its superpopulation interpretation is the coefficient
among subjects available at that cell,
\(E[X(Y_{jh}-X^\top\beta)\mid A_{jh}=1,\mathcal C]=0\). If availability is
associated with the design or outcome, this selected-population target may
differ materially from the planned-population target.

The policy enters `plan_population()` scientific identity. A result records
the policy again; a view cannot switch it.

### 3.3 Rank and coefficient-specific estimability

For either policy, \(X_{jh}\) means the rows of the design admitted to the
cell. Record \(r_{jh}=\operatorname{rank}(X_{jh})\) and residual degrees of
freedom \(n_{jh}-r_{jh}\). A saturated, empty, or rank-deficient cell is not
collapsed into generic missingness:

- `empty_subject_set`: \(n_{jh}=0\);
- `rank_deficient`: \(r_{jh}<\operatorname{ncol}(X)\);
- `coefficient_not_estimable`: the requested coefficient/contrast is outside
  the row space of the admitted design; and
- `planned_subject_unavailable`: the all-planned policy was not observed.

A rank-deficient design may still support an estimable contrast. The result
therefore records estimability per reported coefficient, not just one Boolean
for the cell.

## 4. Coverage quantities

<!-- population-diagnostic: node_sample_size -->
### 4.1 Node-specific sample size

\[
n_{jh}=\sum_{i\in I_0}A_{ijh},\qquad
f_{jh}=n_{jh}/|I_0|.
\]

Both count and fraction are reported. Operator coverage
\(C^{P}_{ij}=1\{m_{ij}>0\}\) is also retained separately from analytic
availability \(A_{ijh}\), because query normalization or non-finite inputs can
remove a covered subject.

<!-- population-diagnostic: effective_sample_size -->
### 4.2 Effective sample size

The inferential effective sample size is defined by the actual regression
weights \(w_{ijh}\):

\[
n^{\mathrm{eff}}_{jh}
=\frac{(\sum_{i\in I_{jh}}w_{ijh})^2}
       {\sum_{i\in I_{jh}}w_{ijh}^2}.
\]

Current OLS uses subject-constant unit weights, so
\(n^{\mathrm{eff}}_{jh}=n_{jh}\). Transport mass is not silently used as a
precision weight. The separate descriptive mass-concentration diagnostic is

\[
n^{\mathrm{mass}}_j=(\sum_{i:m_{ij}>0}m_{ij})^2/
                     \sum_{i:m_{ij}>0}m_{ij}^2.
\]

It describes whether a node's territory is dominated by a few subjects; it is
not the degrees of freedom of the OLS fit.

<!-- population-diagnostic: sink_mass -->
### 4.3 Sink mass

Two quantities must not be conflated:

\[
s_i^{\mathrm{territory}}=
\frac{\sum_x\mu_{i,x}P_{i,x\perp}}{\sum_x\mu_{i,x}},
\qquad
s_{ih}^{\mathrm{budget}}=
\sum_xP_{i,x\perp}\langle H_h,G_{i,x}\rangle_F.
\]

Sink territory is data-free and nonnegative. Sink budget is data-dependent,
signed, and query-specific. Report both with the transport signature.

<!-- population-diagnostic: transport_quality -->
### 4.4 Transport quality

“Transport quality” is a diagnostic vector, not one favorable score:

\[
Q_i=(\text{provenance},\text{cross-fit status},\text{displacement},
     \text{entropy},\text{perplexity},s_i^{\mathrm{territory}},
     \text{all-sink rows},\{m_{ij}\}_j).
\]

Displacement, entropy, and perplexity use the definitions in
`population-form-v1` section 7.5. Report their named components and units; do
not average them into a scalar without a separately frozen rule. At node level,
report coverage count/fraction, `mass_n_eff`, and summaries of the contributing
subjects' transport diagnostics. A high alignment statistic cannot override
poor sink or coverage diagnostics.

The default interpretive-view warning thresholds are a planned-subject coverage
fraction below `0.8`, mean retained territory below `0.7` (equivalently, high
sink territory), or removal of at least `0.2` of a cell's primary contributors
under the declared transport-quality sensitivity threshold. A protocol may
declare different thresholds before inspection, but the values and every
resulting filter must travel in the diagnostic-view receipt. Warnings do not
exclude subjects or change the primary estimand.

## 5. Covariates associated with coverage or transport

Informative coverage is diagnosed, not assumed away.

1. For every model covariate and prespecified transport-quality component,
   report its association with \(A_{ijh}\) or the relevant continuous mass
   across nodes. Continuous covariates use a declared standardized mean
   difference or availability model coefficient; categorical covariates use
   level-specific coverage. The diagnostic must identify the node/query and
   actual subject set.
2. A covariate associated with availability may be included in `model =` if it
   was measured independently of the evaluated response and the model was
   declared before looking at local results. Its inclusion changes the
   conditional coefficient but does not turn an available-at-node target into
   the all-planned target.
3. A response-derived or non-cross-fitted functional transport diagnostic may
   not be used as an adjustment covariate in the same evaluation. That would
   condition on an adaptively selected descendant of the outcome.
4. Missing-indicator adjustment, zero filling, and nodewise automatic covariate
   selection are forbidden. If a covariate is absent for one planned subject,
   the plan is invalid; it is not another local coverage policy.
5. Inverse-probability or transport-quality weighting requires a frozen
   selection/weight model, positivity diagnostics, weight truncation policy,
   and uncertainty propagation. It is not admitted by `population-estimand-v1`.
6. Report rank and coefficient estimability after subsetting. Coverage tied to
   a factor level can remove that level and make one coefficient nonestimable
   even when the intercept remains estimable.

These diagnostics are warnings about external interpretation, not automatic
permission to choose the target that gives more non-missing maps.

## 6. Required provenance schema

No coefficient may be returned without recoverable knowledge of its inputs.
Every population result must carry:

| field | requirement |
|---|---|
| `planned_subjects` | ordered identifiers for \(I_0\) |
| `subject_plan_id` | one geometry scientific-plan id per planned subject |
| `transport_signature` | one signature per planned subject, bound to that identifier |
| `coverage_policy` | `all_planned` or `available_at_node` |
| `operator_coverage` | subject-by-node logical coverage, including an explicit sink rule |
| `analytic_availability` | subject-by-node-by-query/coordinate admission, or a lossless equivalent |
| `subject_set_id` | one dictionary key per reported cell |
| `subject_set_dictionary` | exact subject identifiers for every key; a hash alone is insufficient |
| `n`, `n_eff`, `mass_n_eff` | cell count, regression effective size, and separate mass concentration |
| `design_rank`, `residual_df` | values for the exact admitted design |
| `coefficient_estimable` | per-cell, per-coefficient/contrast status and reason |
| `transport_quality` | named diagnostics with units and transport identity |
| `conditioning` | explicit exclusion of transport- and coverage-estimation uncertainty |

Dictionary encoding is permitted so repeated subject sets are stored once.
Exact identifiers must remain recoverable after serialization; counts and
digests alone do not satisfy the contract.

The current `plan_population()` already seals planned subject plan ids,
transport signatures, model data, normalization, and fit order into its
scientific identity. Current population receipts preserve the per-subject plan
and transport table. PE-D1 extends that provenance to the cell-specific fields
above and distinguishes operator absence from a numeric zero.

## 7. Views and aggregation

A view may change presentation or apply a fixed linear query; it may not
silently change the population target.

- A cellwise view preserves `coverage_policy`, `subject_set_id`, rank,
  estimability, and all linked dictionary entries exactly.
- A fixed linear query over complete-form coordinates uses the intersection of
  the coordinate availability needed for that output, then applies the plan's
  declared policy. It cannot reuse one input coordinate's count for another.
- A spatial or scale aggregation is valid when it is computed subject-first
  from the underlying transported values and then refit. For ordinary
  unweighted OLS, summing already-fitted coefficients is exactly the same
  operation when every source cell has the identical subject set and design;
  implementations may use that linear identity but must check the sets and
  preserve their common id. Combining coefficients estimated from different
  subject sets does not define a common coefficient and is refused.
- If an aggregation cannot recover subject-level values, it remains a
  descriptive summary of heterogeneous targets and must expose all source set
  ids. It cannot be labelled a population coefficient at the aggregate node.
- The sink remains a separately labelled audit row and is never folded into an
  ordinary spatial aggregate.

Thus no view can turn `available_at_node` into `all_planned`, or vice versa,
merely because the plotted arrays have the same dimensions.

## 8. Population component conservation

When total, coherent, and configuration ledgers are transported and fitted
under the same plan, query bank, model, normalization, and exact cellwise
subject sets, linearity requires

\[
\beta_{\mathrm{total}} = \beta_{\mathrm{coherent}} +
\beta_{\mathrm{configuration}}
\]

for every estimable coefficient, node, and query. `population_decomposition()`
checks the subject-value and coefficient identities and exposes the joint
cross-component coefficient covariance needed for derived contrasts. Separate
componentwise normalization is nonadditive and therefore refuses, as do
different plans, models, queries, or subject-set identities. Coherent and
configuration remain estimand components, not separate biological mechanisms.

## 9. Executable example and implementation gate

From the repository root:

```sh
Rscript design/oracles/population-estimand-targets.R
```

The fixture has five planned subjects and three nodes with coverage `5/4/4`.
The all-planned policy returns a coefficient only at the fully covered node;
the available policy fits exact, recoverable four-subject sets at the other
two. Availability is deliberately associated with the covariate, so the
selected fits differ for a visible reason. The oracle also separates OLS
`n_eff` from transported-mass `mass_n_eff`.

`tests/testthat/test-population-estimand-contract.R` fails if either target,
any required diagnostic, the covariate rules, the conditioning statement, or
the recoverable-set schema disappears. It also checks that current public
plans and result receipts bind each planned subject to the same stable
geometry-plan and transport identities. PE-D1 adds executable variable-
coverage package tests against the same oracle after the output schema lands.
