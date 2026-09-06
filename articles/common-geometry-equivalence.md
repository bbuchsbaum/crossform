# Contrast, crossnobis, and linear RSA as one geometry

Why can one geometry plan return a contrast energy, an RDM and an RSA
coefficient? Each is a fixed weighted reading of the same effect
geometry. This guide derives that relationship, checks a small example,
and compares the package outputs with independent calculations.

Read the
[introduction](https://bbuchsbaum.github.io/crossform/articles/introduction.md)
first for the analysis workflow. Here the goal is to understand the
algebra. The result covers fixed linear queries; nonlinear fitting and
its independent evaluation are treated at the end.

The normative derivation and exact independent implementation are
[`design/common-geometry-equivalence.md`](https://github.com/bbuchsbaum/crossform/blob/main/design/common-geometry-equivalence.md)
and
[`design/oracles/common-geometry-equivalence.R`](https://github.com/bbuchsbaum/crossform/blob/main/design/oracles/common-geometry-equivalence.R).
The calculations below reproduce their central identities without
sourcing package internals or the oracle file.

## The theorem

Let $`B_r`$ be the effects-by-features matrix estimated in partition
$`r`$, $`H`$ a fixed effect-side operator, $`K`$ a fixed neural metric,
and $`\Gamma`$ the declared partition-pair weights. Define

``` math
G_K = \sum_{r,s} \Gamma_{rs} B_r K B_s^\top,
\qquad
\mathcal E(H,K) = \langle H,G_K\rangle_F.
```

The Frobenius inner product, $`\langle H,G_K\rangle_F`$, multiplies
matching matrix entries and sums them. The operator $`H`$ therefore
specifies which aspects of the geometry to read. Linearity of trace
gives an equivalent expression in terms of the partition effects:

``` math
\mathcal E(H,K)
= \sum_{r,s}\Gamma_{rs}
  \operatorname{tr}(H^\top B_r K B_s^\top).
```

Three substitutions give the familiar estimators.

- A contrast $`c`$ uses $`H=cc^\top`$, so the output is $`c^\top G_Kc`$.
- An RDM entry for effects $`i,j`$ uses $`u_{ij}=e_i-e_j`$ and
  $`H=u_{ij}u_{ij}^\top`$. With a declared noise precision for $`K`$,
  this is crossnobis.
- A fixed OLS RSA coefficient uses the row $`a_\ell^\top`$ of
  $`(X^\top X)^{-1}X^\top`$ and the RDM adjoint
  $`H=\mathcal D^*(a_\ell)=\sum_{i<j}a_{\ell,ij}u_{ij}u_{ij}^\top`$.

The equalities concern the estimand. Crossnobis additionally needs
independent evaluated partitions and an honestly fixed or cross-fitted
precision. RSA additionally needs a fixed, full-rank design and a linear
objective.

## A hand-sized exact case

The fixture has two partitions, three effects, two features,
off-diagonal partition weights, and a non-identity metric. Every
displayed value is computed in this article.

``` r

effect_labels <- c("a", "b", "c")
blocks <- list(
  run1 = matrix(c(1, 0, 0, 1, 1, 1), 3L, 2L, byrow = TRUE,
                dimnames = list(effect_labels, NULL)),
  run2 = matrix(c(2, 0, 0, 2, 1, 3), 3L, 2L, byrow = TRUE,
                dimnames = list(effect_labels, NULL))
)
gamma <- matrix(c(0, 0.5, 0.5, 0), 2L, 2L,
                dimnames = list(names(blocks), names(blocks)))
metric <- diag(c(2, 1))
geometry <- oracle_geometry(blocks, gamma, metric)
knitr::kable(geometry, digits = 3, caption = "Generated common geometry")
```

|     |   a |   b |   c |
|:----|----:|----:|----:|
| a   |   4 | 0.0 | 3.0 |
| b   |   0 | 2.0 | 2.5 |
| c   |   3 | 2.5 | 5.0 |

Generated common geometry {.table}

The pairwise squared distances are generated from $`u^\top G_Ku`$. The
direct court below independently expands every partition, effect, and
feature index.

``` r

pairs <- pair_matrix(length(effect_labels))
pair_labels <- apply(pairs, 1L, function(x) {
  paste(effect_labels[x], collapse = " - ")
})
rdm_values <- rdm_from_geometry(geometry)
direct_values <- vapply(seq_len(nrow(pairs)), function(k) {
  u <- numeric(length(effect_labels))
  u[pairs[k, ]] <- c(1, -1)
  oracle_energy(tcrossprod(u), blocks, gamma, metric)
}, numeric(1))
hand_rdm <- data.frame(
  pair = pair_labels,
  from_geometry = rdm_values,
  direct_index_oracle = direct_values,
  absolute_error = abs(rdm_values - direct_values)
)
knitr::kable(hand_rdm, digits = 15)
```

| pair  | from_geometry | direct_index_oracle | absolute_error |
|:------|--------------:|--------------------:|---------------:|
| a - b |             6 |                   6 |              0 |
| a - c |             3 |                   3 |              0 |
| b - c |             2 |                   2 |              0 |

Now fit a fixed linear RSA model to the generated RDM. The same
coefficients are recovered by compiling each coefficient row back into
an $`H`$ operator.

``` r

rsa_design <- cbind(`(Intercept)` = 1, category = c(1, 0, 1))
coefficient_map <- solve(crossprod(rsa_design), t(rsa_design))
rsa_from_rdm <- drop(coefficient_map %*% rdm_values)
rsa_from_adjoint <- vapply(seq_len(nrow(coefficient_map)), function(term) {
  operator <- matrix(0, length(effect_labels), length(effect_labels))
  for (k in seq_len(nrow(pairs))) {
    u <- numeric(length(effect_labels))
    u[pairs[k, ]] <- c(1, -1)
    operator <- operator + coefficient_map[term, k] * tcrossprod(u)
  }
  oracle_energy(operator, blocks, gamma, metric)
}, numeric(1))
hand_rsa <- data.frame(
  term = colnames(rsa_design),
  from_materialized_rdm = rsa_from_rdm,
  from_compiled_adjoint = rsa_from_adjoint,
  absolute_error = abs(rsa_from_rdm - rsa_from_adjoint)
)
knitr::kable(hand_rsa, digits = 15)
```

|  | term | from_materialized_rdm | from_compiled_adjoint | absolute_error |
|:---|:---|---:|---:|---:|
| (Intercept) | (Intercept) | 3 | 3 | 0 |
| category | category | 1 | 1 | 0 |

The public
[`crossnobis()`](https://bbuchsbaum.github.io/crossform/reference/crossnobis.md),
[`rdm()`](https://bbuchsbaum.github.io/crossform/reference/rdm.md), and
[`rsa()`](https://bbuchsbaum.github.io/crossform/reference/rsa.md)
results also agree with these direct calculations to the declared
tolerance of `1e-12`; this is checked when the article runs.

## Randomized production-to-oracle evidence

The same comparison can be repeated with different effect estimates,
positive definite metrics, contrasts and model RDMs. The following
summary compares six generated cases against the explicit matrix
definition above. The full calculation remains in this vignette’s R
source.

| Cases | Largest absolute error | Tolerance |
|------:|-----------------------:|----------:|
|     6 |                2.8e-14 |     5e-10 |

The package test court extends this article with rank-deficient PSD
metrics, unequal partition weights, near-singular metrics, nonfinite
refusals, permutation equivariance, scaling laws, and a 50-seed deep
mode. See the [`test-common-geometry-*`
files](https://github.com/bbuchsbaum/crossform/tree/main/tests/testthat)
and their independent oracle link above.

## External differential evidence

The external comparison uses a deterministic six-condition, four-run
fixture, a fixed non-spherical pooled-residual precision, local support
normalization, and pinned Python packages. It compares the same RDM
vector and fixed OLS RSA coefficients to `rsatoolbox`. The table is read
from the generated package certification artifact; no result is typed
into this article.

| Quantity             | Entries | Largest absolute error | Tolerance |
|:---------------------|--------:|-----------------------:|----------:|
| Crossnobis distances |      60 |                  4e-15 |     1e-10 |
| OLS RSA coefficients |      20 |                  1e-15 |     1e-10 |

Fixed-estimand parity against rsatoolbox 0.3.2 {.table}

The full fixture, reproduction instructions, convention mapping, and
source/artifact manifest live in
[`exemplars/rsatoolbox-parity`](https://github.com/bbuchsbaum/crossform/tree/main/exemplars/rsatoolbox-parity).
The manifest makes an edited producer, claim, environment lock, or
result a test failure until the external evidence is regenerated or
explicitly revalidated.

## Estimator claim table

This generated table separates an algebraic mapping from the extra
assumptions that license a familiar estimator name.

| Estimator | Mapping and requirements | Scope |
|:---|:---|:---|
| Contrast energy | `H = c c^T`: fixed contrast, metric and partition weights; centering is needed for translation invariance. | Equivalent special case |
| Crossvalidated squared-distance RDM | `H = u_ij u_ij^T`: zero-sum pair differences, fixed metric, declared normalization and independent cross-partition estimates. | Equivalent special case |
| Crossnobis RDM | Same pair operator, with fixed or honestly cross-fitted noise precision and independent cross-partition estimates. | Equivalent special case |
| Multiple-regression RSA coefficient | `H = D*(a_l)`: fixed full-rank OLS design, pair order, intercept and model scaling; inherits the source geometry’s metric and partition weights. | Equivalent special case |
| Correlation/cosine RDM or RSA | Observed-pattern normalization is nonlinear; there is no fixed H in general. | Outside this theorem |
| Spearman/Kendall RSA | Ranking is nonlinear and depends on the tie policy. | Outside this theorem |
| Adaptive classifier, kernel, or embedding | The operator is learned from evaluated data; the training rule changes the estimand. | Outside this theorem |

“Equivalent special case” means unification of the named fixed estimand
under the stated assumptions. It does not mean that Crossform implements
every estimator family, preprocessing convention, uncertainty procedure,
null law, or population generalization target associated with that name.

## Learning a form, then freezing its query

[`fit_geometry()`](https://bbuchsbaum.github.io/crossform/reference/fit_geometry.md)
adds a nonlinear training step: it learns a regularized, rank-limited
PSD form $`F`$ from signed training geometry and a declared model
kernel. That optimization is outside the fixed-query theorem above. Once
training and selection are complete, however, independent evaluation
returns to its algebra:

``` math
\Delta(F;G_{\mathrm{test}})
= 2\langle F,G_{\mathrm{test}}\rangle_F-\|F\|_F^2.
```

For frozen $`F`$, this is a fixed linear geometry query plus a known
offset. If
$`\mathbb E[G_{\mathrm{test}}\mid\mathrm{training}]=G^\star`$, then

``` math
\mathbb E[\Delta\mid\mathrm{training}]
=\|G^\star\|_F^2-\|G^\star-F\|_F^2.
```

Conditional unbiasedness requires evaluation observations independent of
all training and selection dependencies. Marking a learned query as
fixed does not establish that fact.
[`score_geometry()`](https://bbuchsbaum.github.io/crossform/reference/score_geometry.md)
checks declared observation origins and a matching evaluation target;
the sampling assumption remains part of the analysis design.

This extension preserves linear evidence decomposition: the inner
products with coherent and configuration forms add under the same frozen
$`F`$, with the prediction cost subtracted once. The nonlinear fits
themselves do not obey that conservation law. The [predictive geometry
guide](https://bbuchsbaum.github.io/crossform/articles/predictive-geometry.md)
([`vignette("predictive-geometry")`](https://bbuchsbaum.github.io/crossform/articles/predictive-geometry.md)
offline) executes the public fit/score workflow and its mode-wise gain
identity. It uses Frobenius loss on geometry, not the fixed RDM
regression loss used in the RSA equivalence above.

## Limits and evidence status

The theorem is algebraic. The hand, randomized, adversarial, and
external courts are computational evidence. None supplies population
coverage, calibrates a null distribution, or proves validity after
data-adaptive model or metric selection. External parity validates only
the mapped fixed-linear case, not all of `rsatoolbox` and not all RSA.

For the broader package evidence boundary, see the
[`certification report`](https://github.com/bbuchsbaum/crossform/blob/main/design/certification-report.md).
For the vocabulary separating theorem, numerical verification,
statistical validation, empirical demonstration, and interpretation, see
the
[`unification contract`](https://github.com/bbuchsbaum/crossform/blob/main/design/unification-contract.md).
