# Future cross-node sampling covariance contract

Status: `cross-node-covariance-v0`, not implemented.

This is an admission contract, not an inference feature. Current unweighted
node-wise population OLS, classical covariance, HC3, and the subject-level
wild bootstrap remain valid without it because each cell is fitted across
participants. They do not combine sampling errors across group nodes. A
cross-node covariance becomes necessary only when a method combines node
errors—for example transported within-subject precision, the variance of a
conserved budget, joint spatial inference, or a later maxT procedure.

## Required object

A future covariance object must bind all of the following into its identity:

- participant and source result identifiers;
- the ordered native-node index and ordered query index;
- a declared error model and fitting provenance;
- covariance scope `joint_native_node_query`;
- representation `dense_symmetric` or `sparse_symmetric`;
- numerical tolerance and compute budget;
- conditioning and any estimated-covariance uncertainty.

For one query with native vector `z`, covariance `Sigma`, and row-stochastic
transport `P` including the sink, the defining law is

```
Cov(P' z) = P' Sigma P.
```

The conserved-budget variance is `1' Sigma 1`, not `sum(diag(Sigma))`.
Because `P 1 = 1`, the full transported covariance preserves that variance
when the sink is retained. The executable, non-package oracle in
`design/oracles/cross-node-covariance.R` checks this law and demonstrates that
the diagonal shortcut is wrong on a three-node PSD fixture.

## Dense and sparse admission

Dense input must be finite, square, index-aligned, symmetric within its
declared tolerance, positive semidefinite within a scale-aware tolerance, and
small enough that its explicit bytes fit the compute policy. Sparse input must
carry a canonical symmetric sparse representation, finite stored entries,
matching indices, a bounded nonzero count, and a sparse PSD certificate or
solver check. Production validation must not densify a sparse object.

Both routes must reject:

- asymmetric or indefinite matrices;
- missing, duplicated, or reordered node/query identities;
- an undeclared error model;
- a dense allocation or sparse nonzero count outside the compute budget;
- implicit assembly from marginal per-node blocks.

The tiny oracle may densify its three-by-three sparse fixture solely to prove
representation equivalence. That is not a production scaling strategy.

## Precision-weighting gate

`normalization = "precision_weighted"` remains refused. Lifting it requires
both this validated covariance object and a provenance-bearing precision input
whose weights equal the appropriate budget variances under the same query and
normalization. Marginal `K`-by-`K` blocks at individual nodes are insufficient.

## Optional later joint inference

maxT, spatial multiple-comparison correction, simultaneous bands, cluster
procedures, and covariance regularization are later methods decisions. They
may consume an admitted cross-node covariance, but they are not prerequisites
for unweighted point estimation and they are not silently included by this
contract. Their calibration and multiplicity guarantees require separate
oracles and evidence.
