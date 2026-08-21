# Common-geometry equivalence theorem

Status: algebraic theorem plus executable exact oracle. This document proves a
closed equivalence claim; it does not claim that every representational
statistic is a geometry query.

## 1. One statistic

At one declared measurement let

- `B_r` be the `q x d` effect-coefficient matrix from partition `r`;
- `H` be a fixed `q x q` effect-side operator;
- `K` be a fixed symmetric `d x d` neural metric; and
- `Gamma` be the `R x R` partition-pair weight matrix.

Define

```
E(H, K) = sum_{r,s} Gamma[r,s] tr(H' B_r K B_s')
        = <H, G_K>_F,
G_K     = sum_{r,s} Gamma[r,s] B_r K B_s'.
```

The equality is only linearity of trace. For the symmetric self geometry used
by contrast, RDM and RSA views, `K = K'` and `Gamma = Gamma'`; therefore
`G_K = G_K'`, and `H` is fixed and symmetric. Algebraically the formula allows
any finite `Gamma`. For the crossvalidated interpretation, nonzero weights
join declared independent partitions, self-pairs are absent, and their scale
is the declared averaging normalization (normally the weights sum to one over
the ordered edges). Changing that normalization rescales every special case.

No effect centering is needed for the general statistic. Centering enters the
special cases through zero-sum contrasts. `K` must be positive semidefinite for
a nonnegative population quadratic form; a finite crossvalidated estimate can
still be negative. A strict Mahalanobis metric additionally needs `K` positive
definite on the pattern-difference subspace. The equalities themselves need
neither positive definiteness nor full rank.

## 2. Contrast energy

For a fixed contrast `c` in `R^q`, set `H_c = c c'`. Then

```
E(H_c, K)
  = tr(c c' G_K)
  = c' G_K c
  = sum_{r,s} Gamma[r,s] (c' B_r) K (c' B_s)'.
```

This is the familiar cross-partition inner product of the contrast patterns.
With `K = I` it is Euclidean contrast energy. A zero-sum contrast (`1'c = 0`)
annihilates a common shift of all effect rows; an uncentred contrast is still a
valid geometry query, but it is not recoverable from pairwise distances.

The rank-one identity `rank(c c') = 1` holds for nonzero `c`; it is a
computational structure, not an inferential assumption.

## 3. Pairwise RDM and crossnobis

For effects `i != j`, let `u_ij = e_i - e_j`. The unordered pairwise RDM
operator is

```
D(G)[i,j] = u_ij' G u_ij = <u_ij u_ij', G>_F.
```

Thus its entry is exactly `E(u_ij u_ij', K)`. Because `1'u_ij = 0`, the RDM
is invariant to the common-effect origin. With a declared fixed noise
precision `K = Sigma_noise^{-1}`, the same value is the familiar crossnobis

```
sum_{r,s} Gamma[r,s]
  (B_r[i,] - B_r[j,]) K (B_s[i,] - B_s[j,])'.
```

Calling it crossnobis adds assumptions, not new arithmetic: evaluated
partition pairs must be independent; the precision must have the declared
noise role; and an estimated precision must be frozen from data independent
of the evaluated effect products. Singular PSD precision yields a
pseudodistance; strict Mahalanobis separation requires positive rank on the
contrast-pattern subspace. Any factor such as division by feature count must
be declared in `K`, the frame, or `Gamma`; it is not an extra hidden factor in
the equivalence.

## 4. Fixed linear RSA through the RDM adjoint

Enumerate the `p = q(q-1)/2` unordered pairs and regard `D` as a linear map
from symmetric geometry to `R^p`. Its Frobenius adjoint is

```
D*(a) = sum_{i<j} a[i,j] u_ij u_ij',
```

because `a' D(G) = <D*(a), G>_F`.

Let `X` be a fixed `p x h` RSA design whose columns are vectorized symmetric,
zero-diagonal model RDMs, optional nuisance RDMs, and an optional intercept.
When `X` has full column rank, ordinary least squares is the fixed linear map

```
A    = (X'X)^-1 X',
beta = A D(G_K).
```

For term `l`, let `a_l'` be row `l` of `A`. Then

```
beta_l = a_l' D(G_K)
       = <D*(a_l), G_K>_F
       = E(D*(a_l), K).
```

This is why fixed multiple-regression RSA can be compiled as one geometry
query. The RDM need not be materialized. Full column rank of `X` is required
to identify the requested coefficient; `D` itself need not be inverted.
Including an intercept changes `A` and hence the estimand, while model scaling
changes coefficient units. Both must be declared. The models, nuisance terms,
and intercept choice must be fixed before the observed geometry is inspected.

For any centred contrast there is also an exact RDM-coordinate identity:

```
c c' = -sum_{i<j} c_i c_j u_ij u_ij',  when sum_i c_i = 0.
```

It states precisely which contrast energies can be recovered from distances.

## 5. Hand-computable court

The executable oracle `design/oracles/common-geometry-equivalence.R` uses two
partitions, three effects, two features,

```
Gamma[1,2] = Gamma[2,1] = 1/2,  K = diag(2, 1),
B_1 = [1 0; 0 1; 1 1],          B_2 = [2 0; 0 2; 1 3].
```

It gives

```
G_K = [4 0 3; 0 2 2.5; 3 2.5 5]
D(G_K) = (ab=6, ac=3, bc=2).
```

The direct partition formula gives the same three values. Regressing those
distances on an intercept and model `(1,0,1)` gives `(intercept=3, model=1)`;
the adjoint operators give those coefficients directly as
`E(D*(a_intercept),K)` and `E(D*(a_model),K)`.

## 6. Boundary of the claim

The equivalence includes fixed bilinear contrast energy, signed pairwise
crossvalidated squared distances, crossnobis under a declared fixed or honestly
cross-fitted precision, and fixed *linear* RSA with a full-rank design.

It explicitly excludes Spearman or Kendall RSA, rank transforms, correlations
or cosine normalization of the observed RDM, data-adaptive model selection,
models fitted to the evaluated geometry, locally re-estimated metrics without
independent training, classifiers with learned nonlinear decision rules,
kernels or embeddings learned from the evaluation data, and unsupported named
estimators. Those operations are nonlinear or data-adaptive; resemblance to
`E(H,K)` is not an equivalence proof.

## 7. Numerical policy

The executable courts compare floating-point values with a scale- and
conditioning-aware bound

```
atol + rtol * operations * kappa_supported * max(1, abs(reference)),
```

where `kappa_supported` is computed from singular values above the declared
rank tolerance and capped before it can turn every comparison into a pass.
Exact hand examples still use exact equality. No BLAS-derived digest is a
scientific oracle.

Near-singular PSD metrics remain valid for the total bilinear geometry and
crossnobis pseudodistance. The coherent/configuration decomposition has a
stronger SPD requirement and refuses deterministically rather than falling
back to an arbitrary inverse. Nonfinite blocks, incompatible axes, singular
RSA designs, and correlation-normalized signed RDMs likewise fail before an
equivalence comparison. There is no numerical warning that silently changes
the estimand and no fallback normalization.

## 8. External parity binding

The versioned `exemplars/rsatoolbox-parity` fixture exercises one supported
external slice of this theorem. Its `crossform::rdm()` values instantiate
Section 3 with a fixed pooled-residual precision, zero-sum pair-difference
operators, four independent run partitions, uniform weight over the six
unordered cross-run pairs, and local frame normalization by support size.
Its `crossform::rsa()` values instantiate Section 4 with a fixed full-rank OLS
design over the row-major upper triangle of that same RDM.

The matched external estimator is `rsatoolbox` 0.3.2
`calc_rdm_crossnobis`: its leave-one-run-out average equals the uniform
unordered-pair average on this balanced fixture, and its division by channel
count equals the declared local frame normalization. The fixed linear RSA
comparison is ordinary least squares on the matched RDM vector. The additional
`ModelWeighted` plus `fit_regress(cosine)` row is compared only after its
documented RMS scaling is undone.

This external evidence supports those mapped crossnobis and fixed-linear-RSA
cases only. It does not validate rank RSA, arbitrary similarity objectives,
adaptive models, different precision estimators, unbalanced cross-validation,
or every statistic named RSA. `results/parity-manifest.csv` binds the pinned
environment, producing sources, convention metadata, algebraic claim, and
recorded outputs; the test court rejects a missing, edited, or stale binding.
