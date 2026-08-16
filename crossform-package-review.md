# `crossform` package design review

Status: design guidance, not implementation
Date: 2026-08-11

Implementation baseline: commit
`01f3763c6aa00ebe50e69fc966d5201a861e4b09` (2026-08-12). The interim
implementation audit used aggregate source SHA-256
`13fe403a2ad6deec5b70c769dd39d226c993254cba4b4a206b5e3036e6c71e4c`.
Its eight remediation decisions are frozen in the opening sections of
`crossform-package-design.md` and `crossform-execution-design.md`; those
sections supersede any incompatible provisional guidance in this review.

## Recommendation

Build the package, but freeze it as a small effect-geometry compiler—not as the full Unified Effect Geometry platform imagined in the earlier ledger.

The durable architecture is sound:

\[
\text{relation} + \text{frame} + \text{pairing}
\longrightarrow
\text{geometry}
\longrightarrow
\text{views}
\]

Part 7 of `searchlight-conversation-ledger.md` and the distilled design should be authoritative. Parts 1–6 are valuable rationale and future research, not the version 0.1 specification. This preserves the strongest idea: voxel, ROI, searchlight, and whole-brain analyses differ by a spatial frame rather than requiring separate engines.

## Corrections needed before implementation

### 1. Fix undirected pairing marginals

The current edge representation stores each unordered pair once, but the proposed left/right marginals depend on the arbitrary edge orientation. With edges `(1, 2)`, `(1, 3)`, and `(2, 3)`, partition 1 appears only on the left and partition 3 only on the right. The marginals therefore do not merely differ through sampling noise.

For an undirected edge of weight \(\omega\), internally expand it into both directions with weight \(\omega/2\), or store the single orientation-invariant endpoint marginal:

\[
\bar B^{U}_j
=
\sum_e\omega_e
\frac{(B_{l(e)}+B_{r(e)})w_j}{2a_j}.
\]

Undirected pairings should expose one `endpoint` marginal; directed pairings should expose `left` and `right`. Amend the marginal law to require exact orientation invariance. This affects the pairing rules and marginal definition in `crossform-package-design.md`, beginning around lines 76 and 189.

### 2. Add an explicit experimental-space contract

Effect names alone cannot establish that two partitions use the same basis, scaling, or units. Add a small `effect_space` value containing:

- coordinate labels;
- a semantic `basis_id`;
- per-coordinate units and scaling;
- optional contrast/design provenance.

RDM and RSA views should require commensurate coordinates. Directed cross-domain pairing in version 0.1 should still require the same experimental space on both sides. Direction describes endpoint roles, not an asymmetric cross-covariance that survives symmetrization.

### 3. Qualify “unbiased cross estimate”

A cross-product is noise-unbiased only under independent, zero-mean estimation errors, fixed extractors and metrics, and a common latent target. The package can validate distinct partition identifiers but cannot prove independence.

Use `estimate = "cross_product"` and print “noise-unbiased under declared independence.” For version 0.1:

- require nonnegative edge weights summing to one;
- reject mixtures of self and off-diagonal edges;
- label self-products explicitly as noise-biased;
- reserve signed pairing weights for a later, separately specified estimand.

### 4. Resolve full-geometry versus direct-query execution

A packed geometry is not always small. At 100,000 measurements:

- with \(q=50\), one packed component is about 1.02 GB;
- with \(q=100\), one packed component is about 4.04 GB.

Storing both total and coherent roughly doubles that before marginals and R copies. This makes direct-query execution foundational, not merely an optimization.

Keep `geometry()` semantically explicit: either materialize full geometry after a memory preflight, or execute a declared query and return a non-durable view. Never silently change the output type. The tension appears between the packed representation and direct-query execution described in `crossform-package-design.md` around lines 135–158 and 471–477.

### 5. Separate localization weights from metric and conservation semantics

“Conservative columns sum to one” conserves the identity feature metric. More generally, if the global diagonal metric is \(D(g)\), conservation requires:

\[
\sum_j w_{jv}=g_v.
\]

Record the target global metric explicitly. Require users to choose `normalization`; do not default it. Also make searchlight distance units explicit—for example, millimetres rather than an ambiguous `radius = 6`.

### 6. Narrow the version 0.1 views

The design alternately includes and defers information and effective rank. Resolve this by omitting both from version 0.1. Metadata can record a claimed noise metric, but it cannot certify the upstream model strongly enough to make generic information values safe.

Include:

- contrast bundles;
- signed squared dissimilarities;
- fixed linear multiple-regression RSA;
- signed eigenvalues, clearly labelled as finite-sample cross-estimate eigenvalues.

Defer information, effective rank, coherence ratios as primary statistics, latent PSD reconstruction, and Bures/population summaries. This matches the final recommendation in `crossform-package-design.md` despite the inconsistency in its current scope list around lines 502–529.

## Recommended version 0.1 boundary

### Include

- validated dense and block-backed relations;
- identity and explicit numeric extractor matrices;
- additive nonnegative sparse frames only;
- abstract and volume domains;
- voxel, searchlight, region, and global frame constructors;
- normalized off-diagonal pairings;
- packed total/coherent geometry and derived configuration;
- dense reference and streamed implementations;
- contrast, signed squared-RDM, and fixed linear RSA views;
- a conditional `neuroim2` adapter under `Suggests`.

### Defer

- factor frames and dense local metrics;
- `lm_extractor()` until the explicit-extractor kernel is stable;
- information and effective rank;
- learned metrics, classifiers, calibration, and permutations;
- surface support if it complicates the domain protocol;
- group models, geometry transport, Bures methods, and workflow infrastructure.

`lm_extractor()` can still enter version 0.1 late if it passes a separate numerical gate: pivoted stable solve, explicit whitening-factor semantics, estimability diagnostics, and failure by default for nonestimable requested effects.

## API guidance

Avoid exporting `frame()` and `spectrum()`: they collide with `graphics::frame` and `stats::spectrum`. A clearer initial surface would be:

```r
rel <- effect_relation(
  sources = list(run1 = B1, run2 = B2, run3 = B3),
  effects = effect_space(
    c("face", "house", "object"),
    basis_id = "condition-means-v1"
  ),
  domain = volume_domain(mask)
)

at <- searchlights(
  domain(rel),
  radius = 6,
  units = "mm",
  normalization = "local"
)

over <- cross_partitions(rel)

g <- geometry(rel, at = at, over = over)
out <- geom_contrast(g, c(face = 1, house = -1))
```

Accept numeric extractor matrices directly for convenience, but immediately canonicalize them into an internal extractor record with dimensions, effect-space identity, and provenance. Use `normalization` consistently; the current design example also uses `normalize`.

Squared RDMs should inherit their scale from the frame metric. A locally normalized frame gives mean weighted feature energy; an unnormalized or conservative frame gives total weighted energy. `rdm()` should not silently rescale it again.

## Recommended build order

### 1. Freeze the mathematical contract

Resolve undirected marginals, edge normalization, experimental-space identity, metric units, and materialization behavior in writing.

### 2. Build an independent dense oracle

Use explicit loops and unpacked matrices so it does not share implementation machinery with the optimized kernel. Add exact tiny examples plus randomized laws, including the corrected marginal law.

### 3. Add sparse streaming and resource preflight

Verify equivalence across block sizes and feature permutations to scale-aware numerical tolerances. Make stochastic null-centering tests seed-fixed and judge them against their simulated standard error rather than an arbitrary zero tolerance.

### 4. Stabilize the public API, then add volume frames and views

The release vignette should calculate voxel, ROI, searchlight, and global results from one relation, demonstrate total equals coherent plus configuration, demonstrate conservative local equals global, and preserve a deliberately negative crossvalidated dissimilarity or eigenvalue.

## Claim discipline

The package should claim verified algebraic and computational equivalence for its fixed bilinear kernel—not equivalence to every complete RSA, cvMANOVA, classifier, or inferential workflow. The verification ledger in `searchlight-conversation-ledger.md`, beginning around line 2459, correctly lists those broader claims as still requiring independent checking.

The three immediate priorities are therefore:

1. correct and freeze the contracts;
2. implement the dense oracle and mathematical laws;
3. build the streamed compiler before adapters and broader scientific views.
