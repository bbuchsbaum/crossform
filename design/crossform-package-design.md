# `crossform`: clean-room package design

Status: living package design
Date: 2026-08-11; package name changed to `crossform` 2026-08-15
Canonical package name: `crossform`. An exact-name CRAN scan on 2026-08-15
found no collision. That check is not GitHub, domain, or trademark clearance.

## Version 0.1 contract freeze

The implementation baseline is commit
`01f3763c6aa00ebe50e69fc966d5201a861e4b09` (2026-08-12). The interim audit
reviewed the same package sources under aggregate SHA-256
`13fe403a2ad6deec5b70c769dd39d226c993254cba4b4a206b5e3036e6c71e4c`.
The following decisions are binding before the public API freezes:

1. Experimental identity is an `effect_space`, not a character vector. It
   contains ordered coordinate identifiers, a semantic `basis_id`,
   per-coordinate units and scaling, optional design/contrast provenance, and a
   stable signature. Named precomputed partitions are aligned to this order;
   partial names, missing or extra coordinates, incompatible bases, and
   incompatible units fail before source access.
2. Neural identity is a compact immutable domain reference containing ordered
   `feature_ids`, feature count, coordinate units/geometry identity, and a stable
   signature. Relations and frames carry the reference and must agree exactly;
   a shared string and feature count are insufficient.
3. Direct execution is component-aware. Total-only work does not allocate local
   relations or calculate coherent geometry or marginals. Coherent-only work
   retains local relations but not total atoms. Configuration calculates total
   and coherent. Complete geometry calculates both and retains
   pairing-appropriate marginals.
4. `workspace_bytes` means the conservative peak of crossform-owned live
   objects and temporaries. Baseline, incremental, and absolute process RSS are
   separate observations. Any hard aggregate RSS limit is an explicit optional
   policy and is enforced only where RSS can be observed reliably.
5. File-backed execution opens each distinct immutable source once per
   execution scope, validates its revision on admission, reuses the owned
   read-only handle across blocks, and closes it exactly once. Content hashing
   is never repeated for every block.
6. Execution receipts distinguish requested policy from observed execution.
   They record actual task progress, blocks, bytes, timings, access modes,
   runtime identities, memory observations, reporter failures, checkpoint
   state, and checked cleanup outcomes. An unobservable BLAS thread count is
   recorded as unknown, never invented as one.
7. Volumetric searchlight geometry is a conditional `neuroim2` adapter under
   `Suggests`. It consumes `neuroim2::searchlight_indices()` and maps stable
   full-volume indices through the domain reference. It does not import ROI
   iteration, data extraction, result types, or parallel state.
8. The supported public constructor is `compile_frame()`, avoiding the
   `graphics::frame` collision. Factor frames, nonlinear queries, and compiler
   decision objects remain internal until an executable public workflow earns
   them. Every exported consumer reconstructs or validates its value objects at
   the boundary.

These contracts supersede looser examples or field sketches later in this
document where they conflict. The API-freeze gate is the installed source
artifact plus semantic, adversarial, resource, adapter, documentation, and
package checks; a green source-tree unit suite alone is not certification.

## Outcome

Build a small R package whose only numerical responsibility is to compile lazy partitioned brain–experiment relations into cross-generalized experimental geometries. Voxelwise, searchlight, ROI, and global analyses are spatial frames. Contrasts, squared-distance RSA, MANOVA-like subspaces, and spectra are views. Cross-validation is a partition-pair relation. Scalarization occurs last.

The irreducible pipeline is:

\[
\boxed{
\texttt{relation}
\xrightarrow[\texttt{pairing}]{\texttt{frame}}
\texttt{geometry}
\xrightarrow{\texttt{query}}
\texttt{view}
}
\]

For the additive fast path, the compiler normal form is

\[
\boxed{V=W\,\Phi_\Gamma(EY)\,C.}
\]

The package does not implement a method catalog. It implements this algebra.

## Design constraints

- R package, clean-room implementation, no rMVPA compatibility layer.
- Fewer than 3,000 executable R lines for version 0.1; tests and documentation should be larger than production code.
- `Matrix` is the only mandatory non-base dependency.
- `shard` is the preferred suggested execution dependency after the sequential
  kernel is proven; it is never required to construct or query geometry.
- No raw-data workflow platform, formula DSL, classifier registry, model registry, CLI, scheduler abstraction, or registration system.
- One complete geometry result type, `effect_geometry`, plus the explicitly
  partial query result `effect_view`; neither varies by legacy method or spatial scope.
- One block-streamed numerical engine, with at most two algebraic lowerings: additive diagonal frames and fixed factor frames.
- Crossvalidated estimates remain signed/indefinite when sampling fluctuation makes them so.
- Learned operations must be trained externally, frozen, and evaluated on independent partitions.

## Evidence from current rMVPA

This proposal is grounded in a read-only audit of rMVPA commit `3b12a8855b06f549e5772ccb9af070639a34752b`.

- The current package collates 82 R files, imports 23 packages, exports about 190 symbols, and contains roughly 49,000 physical lines under `R/` ([DESCRIPTION](https://github.com/bbuchsbaum/rMVPA/blob/3b12a8855b06f549e5772ccb9af070639a34752b/DESCRIPTION)).
- Its documented mental model is dataset → design → model specification → regional/searchlight engine → result, with RSA and cross-decoding as first-class method families ([README](https://github.com/bbuchsbaum/rMVPA/blob/3b12a8855b06f549e5772ccb9af070639a34752b/README.md#the-mental-model)).
- The public workflow selects `searchlight`, `regional`, or `global` mode and dispatches to three execution functions ([workflow_api.R](https://github.com/bbuchsbaum/rMVPA/blob/3b12a8855b06f549e5772ccb9af070639a34752b/R/workflow_api.R)).
- Searchlight execution selects among legacy, SWIFT, dual-LDA, and naive cross-decoding policies through eligibility and fallback rules ([searchlight_engine.R](https://github.com/bbuchsbaum/rMVPA/blob/3b12a8855b06f549e5772ccb9af070639a34752b/R/searchlight_engine.R)).
- A standardized ROI result is converted back to a one-row tibble required by existing method-specific combiners ([data_roi_result.R](https://github.com/bbuchsbaum/rMVPA/blob/3b12a8855b06f549e5772ccb9af070639a34752b/R/data_roi_result.R)).
- The general image-dataset constructor rejects one-feature inputs for a Feature-RSA-specific reason ([dataset.R](https://github.com/bbuchsbaum/rMVPA/blob/3b12a8855b06f549e5772ccb9af070639a34752b/R/dataset.R#L263)).

These facts identify the architectural pressure; they do not prescribe `crossform` internals. No source is to be copied.

## 1. Canonical mathematics

### 1.1 Orientation

Use one internal coordinate convention:

\[
Y_r\in\mathbb R^{n_r\times p},
\qquad
E_r\in\mathbb R^{q\times n_r},
\qquad
B_r=E_rY_r\in\mathbb R^{q\times p}.
\]

- rows of \(Y_r\): observations/time points;
- columns of \(Y_r\): native neural features;
- rows of \(B_r\): common experimental coordinates;
- columns of \(B_r\): native neural features.

The public relation abstraction is orientation-free; adapters enforce this convention once.

### 1.2 Pairing semantics

Do not make a dense \(\Gamma\) matrix the primary API. Use a canonical edge table

```text
left | right | weight
```

with these rules:

1. Weights are finite and normally sum to one.
2. Every edge contributes a symmetrized cross-product.
3. `cross_partitions()` creates each unordered off-diagonal pair once and marks
   the pairing as undirected; stored edge orientation has no scientific meaning.
4. Directed cross-domain pairings retain distinct left and right roles.
5. Self-pairs are permitted only with `bias = "self-product"` recorded.
6. Duplicate edges are combined during compilation.

For edge set \(\mathcal E\), feature atom \(v\) is

\[
H_v
=\sum_{e\in\mathcal E}\omega_e
\operatorname{sym}
\left(b_{l(e),v}b_{r(e),v}^\top\right),
\]

where \(\operatorname{sym}(X)=(X+X^\top)/2\).

This edge definition avoids double-counting ambiguity and makes units predictable.

### 1.3 General spatial measurement

For spatial element \(j\) with fixed neural metric \(M_j\succeq0\),

\[
G_j
=\sum_{e\in\mathcal E}\omega_e
\operatorname{sym}
\left(B_{l(e)}M_jB_{r(e)}^\top\right).
\]

The package initially supports two representations.

**Additive frame**

\[
M_j=D(w_j),\qquad w_j\ge0.
\]

This supports point, searchlight, ROI, parcel, smoothing, and global measurements through one sparse matrix \(W\).

**Factor frame**

\[
M_j=L_j^\top L_j.
\]

This covers a modest number of fixed dense/low-rank measurements. It is an extension path after the additive kernel is proven.

### 1.4 Packed geometry

Let \(h=q(q+1)/2\). Use an isometric symmetric vectorization:

\[
\operatorname{tr}(A^\top G)
=\operatorname{svec}(A)^\top\operatorname{svec}(G)
\]

for symmetric \(A,G\), with \(\sqrt2\)-scaled off-diagonals.

Stack \(\operatorname{svec}(H_v)^\top\) into \(Z\in\mathbb R^{p\times h}\). For additive frame \(W\in\mathbb R^{m\times p}\),

\[
\mathcal G^{\mathrm{total}}=WZ.
\]

Linear query operators form \(C\in\mathbb R^{h\times k}\), giving

\[
V=WZC.
\]

With a declared direct query, compile as \(W(ZC)\) and return `effect_view`;
with complete geometry, materialize \(WZ\) and return `effect_geometry`. The
queried scalar semantics agree, but the durable result capabilities are
deliberately different.

### 1.5 Coherent/configuration decomposition

For frame row \(w_j\) with mass \(a_j=\mathbf1^\top w_j>0\),

\[
D(w_j)
=\frac{w_jw_j^\top}{a_j}
+\left[D(w_j)-\frac{w_jw_j^\top}{a_j}\right].
\]

Both terms are PSD. Let \(u_{jr}=B_rw_j\). Then

\[
G_j^{\mathrm{coherent}}
=\frac1{a_j}\sum_e\omega_e
\operatorname{sym}(u_{j,l(e)}u_{j,r(e)}^\top),
\]

and

\[
G_j^{\mathrm{configuration}}
=G_j^{\mathrm{total}}-G_j^{\mathrm{coherent}}.
\]

Configuration is derived, not stored. A point frame has zero configuration.

The word “coherent” means coherent under the effective nonnegative metric weights. If diagonal noise precision is folded into the weights, this is a precision-weighted common mode, not necessarily an arithmetic regional mean. Metadata and printing must say which.

### 1.6 Pairing marginals and signed effects

A pairing graph does not generally determine one signed effect. For an
undirected pairing stored once per unordered edge, left/right marginals depend
on arbitrary storage orientation and are therefore invalid. Store the
orientation-invariant endpoint marginal:

\[
\bar B^{\mathrm{endpoint}}_j
=\frac12\sum_e\omega_e
\frac{(B_{l(e)}+B_{r(e)})w_j}{a_j}.
\]

For a genuinely directed pairing, such as encoding on the left and retrieval
on the right, store separate role-specific summaries:

\[
\bar B^{L}_j
=\sum_e\omega_e\frac{B_{l(e)}w_j}{a_j},
\qquad
\bar B^{R}_j
=\sum_e\omega_e\frac{B_{r(e)}w_j}{a_j}.
\]

Views report `endpoint` for undirected pairings and `left`/`right` for directed
pairings. Do not expose left/right summaries for an undirected edge table and do
not store an ambiguous single `mean`.

### 1.7 Frame normalization

Support exactly three declared states:

- `normalization = "none"`: raw supplied weights;
- `normalization = "local"`: every nonempty row sums to one;
- `normalization = "conservative"`: every included feature column sums to one.

Only conservative frames guarantee

\[
\sum_jG_j^{\mathrm{total}}=G^{\mathrm{global}}.
\]

This conservation applies to geometry and linear views, not nonlinear information, accuracy, ranks, ratios, or thresholded maps.

## 2. Durable representations

### 2.1 `effect_relation`

```r
structure(list(
  sources      = sources,
  extractors   = extractors,
  effects      = effect_names,
  partitions   = partition_table,
  n_features   = p,
  domain       = domain,
  capabilities = list(
    seekable = TRUE,
    reopenable = FALSE,
    concurrent_read_safe = FALSE,
    descriptor_serializable = FALSE,
    in_memory = TRUE
  ),
  units        = units,
  provenance   = provenance
), class = "effect_relation")
```

Only computational protocol:

```r
relation_block(x, partition, features)
```

Return shape: `q × length(features)`.

A source needs a block reader returning observation × feature data and an
explicit capability declaration used by the execution compiler. An extractor
is normally a numeric matrix, though a later developer interface may permit a
linear operator with `apply()` and dimension metadata.

### 2.2 `effect_geometry`

```r
structure(list(
  marginals   = pairing_appropriate_marginals,
  total       = total_svec,
  coherent    = coherent_svec,
  effects     = effect_names,
  index       = frame_index,
  frame       = compact_frame_metadata,
  pairing     = compact_pairing_metadata,
  metric      = metric_metadata,
  storage     = storage_metadata,
  receipt     = execution_receipt,
  diagnostics = diagnostics,
  provenance  = provenance
), class = "effect_geometry")
```

Shapes:

- undirected endpoint marginal, or directed left/right marginals:
  `measurement × effect`;
- total/coherent: `measurement × packed-symmetric-coordinate`.

`total` and `coherent` may be in-memory matrices or block-backed symmetric
fields behind one read-only matrix-like protocol. Storage changes access cost,
not semantic completeness.

Derived configuration is `total - coherent`.

The geometry type records estimator status:

- `estimate = "cross"`: unbiased/signed symmetric estimate; may be indefinite;
- `estimate = "self"`: self-product estimate with positive noise bias;
- `estimate = "latent_psd"`: modeled PSD object, never silently substituted for a cross estimate.

Frames, pairings, domains, and queries are small validated values, not durable result families.

### 2.3 `effect_view`

A direct compiled query does not masquerade as complete geometry:

```r
structure(list(
  values     = values,
  query      = compiled_query,
  components = component_metadata,
  index      = frame_index,
  units      = units,
  receipt    = execution_receipt,
  provenance = provenance
), class = "effect_view")
```

An `effect_view` can be rendered or combined under its declared query, but
cannot answer a new geometry query without re-executing the relation.

## 3. Small value objects

### `effect_frame`

Principally a `Matrix::dgCMatrix` plus index, row mass, normalization, domain identity, representation, and effective metric description.

### `effect_pairing`

A data frame with `left`, `right`, `weight`, plus endpoint semantics, normalization, self-pair/bias status, and declared independence.

### `geom_query`

For a linear view, principally a packed query matrix plus labels, units, required metric properties, and component requirements. Nonlinear views receive unpacked matrices only when necessary.

### `compute_policy`

An immutable serializable value containing numeric execution choices only:
memory budget, I/O block size, reduction microblock, row/coordinate tile sizes,
inflight/reorder bounds, worker count, source staging, and thread request. It
contains no reporter closure or checkpoint path.

The scientific plan hash, numeric execution-policy hash, and non-semantic
reporter/checkpoint details are recorded separately in the execution receipt.

### `feature_domain`

Metadata and two capabilities:

```text
build_frame(specification)
render(values)
```

The kernel sees feature count, block reader, and compiled frame only.

## 4. Public API

Version 0.1 exports only functions that lead to an implemented workflow.

```r
# Relations
relation()
relation_block()
effect_extractor()
lm_extractor()
effect_space()
file_matrix_source()
source_capabilities()

# Domains and frames
abstract_domain()
volume_domain()
compile_frame()
voxelwise()
searchlights()
regions()
whole_brain()
additive_frame()
neuroim2_volume_domain()
neuroim2_searchlights()

# Pairings
pairing()
cross_partitions()
pairing_marginals()

# Computation
compute_policy()
materialize_geometry()
evaluate_geometry()

# Views
contrast_energy()
geometry_component()
query_geometry()
rdm()
rsa()
geometry_spectrum()

# Numerical contracts
numerical_contract()
numerical_agreement()
```

Factor frames, nonlinear queries, cross-domain pairings, information summaries,
and presentation adapters remain internal or deferred until a supported public
workflow and contract exist. Surface support remains a suggested adapter or
version 0.2 concern unless it can reuse the exact domain protocol without
expanding the core.

### Raw response example

```r
domain <- volume_domain(mask)
rel <- relation(
  bold_runs,
  extract = lm_extractor(
    design  = X_runs,
    effects = C,
    whiten  = L_runs
  ),
  domain = domain
)

at <- compile_frame(searchlights(radius = 6, normalization = "local"), domain)
over <- cross_partitions(rel)

g <- materialize_geometry(
  rel,
  at = at,
  over = over,
  storage = "memory",
  compute = compute_policy(workers = 1)
)

fh <- contrast_energy(g, c(face = 1, house = -1))
```

The contrast bundle contains pairing-appropriate signed effects, coherent
energy, configuration energy, total energy, and an optional guarded coherence
ratio.

For an undirected `cross_partitions()` pairing, the signed member is the
orientation-invariant endpoint effect. For a directed pairing, the bundle
contains distinct left/right signed effects.

### Direct-query example

```r
fh_direct <- evaluate_geometry(
  rel,
  at = at,
  over = over,
  query = bilinear_query(tcrossprod(c(face = 1, house = -1))),
  compute = compute_policy(workers = 1)
)
```

This returns `effect_view`. It is cheaper when only that contrast is needed,
but it is not a partially populated `effect_geometry`.

### Precomputed effects example

```r
rel <- relation(
  list(run1 = B1, run2 = B2, run3 = B3),
  effects = c("face", "house", "object"),
  domain = volume_domain(mask)
)
```

Identity extractors make this equivalent to raw-source execution through the corresponding \(E_rY_r\).

### Deferred cross-domain example

Cross-domain pairing is not part of the version 0.1 public API. A later
contract may admit syntax such as the following only after it can preserve
explicit domain identity and role-specific marginals:

```r
g_x <- materialize_geometry(
  rel,
  at = regions(atlas),
  over = cross_domains(
    domain = "phase",
    from = "encoding",
    to = "retrieval",
    match = "run"
  )
)
```

Signed views report encoding and retrieval marginals separately; geometry reports their reproducible cross-relation.

## 5. Views

### Contrast

For \(c\), compile \(A_c=cc^\top\). Return a compact view with:

- `signed_left`, `signed_right`, and `signed_mean` only when pairing symmetry makes the latter meaningful;
- `coherent`;
- `configuration`;
- `total`;
- `coherence_fraction` only with guarded denominator and explicit descriptive status.

### Squared-distance RDM

Compile the fixed linear map from packed \(G\) to squared distances. Preserve scaling (sum versus mean over features) in metadata.

### Multiple-regression RSA

Compile one explicit RDM vectorization, intercept policy, weights, nuisance RDMs, and regression solve into a query matrix. Test against explicit RDM construction and the same regression. Spearman and normalized-correlation RSA are separate nonlinear views.

### Spectrum

For a cross estimate, return signed eigenvalues and label them sampling estimates. Do not call them variances or information roots. PSD-only summaries require `latent_psd` or a tolerance-qualified population object.

### Information

Refuse unless all are true:

- geometry is PSD or a declared latent PSD estimate;
- neural metric is noise-normalized in documented units;
- design scaling/prior is defined;
- the requested Gaussian information model is declared.

Then compute \(\frac12\log\det(I+G)\) with stable eigen/log1p numerics.

## 6. Compiler and execution

The detailed resource, scheduling, failure, checkpoint, and observability design
is in [crossform-execution-design.md](crossform-execution-design.md). Its
governing decision is that additive-frame searchlights with fixed bilinear
queries are rows of sparse \(W\), not jobs submitted to workers. Locally
trained, locally estimated, nonlinear, or generic factor-frame extensions are
outside this exact lowering.

Execution is controlled by an explicit immutable compute policy; the package
does not inspect or mutate a global `future` plan. The version-0.1 default is
deterministic in-process streaming, and `workers != 1` is rejected before
source access. After that reference is verified, the
preferred first parallel adapter is `shard`: it can publish response or compiled
relation matrices once as immutable shared inputs and supervise coarse
feature-block tasks without per-worker full copies. `shard` remains beneath the
relation protocol and geometry kernel; it does not introduce a second model,
engine, or result path.

The execution pipeline is:

1. compile/validate semantics and memory requirements;
2. compute relation and cross-atom blocks over canonical feature intervals;
3. reduce blocks in canonical order through sparse frame multiplication;
4. construct a result only after completeness and algebraic laws are verified;
5. attach a durable execution receipt.

Undirected pairings retain one orientation-invariant endpoint marginal;
directed pairings retain separate left and right relation marginals. Runtime
failures stop at task boundaries and never become apparently valid per-location
result rows. Parallel backends are replaceable task schedulers, not scientific
engines.

### Compile-time validation

Fail once before streaming when:

- sources and extractors disagree in observation count;
- partitions disagree in experimental coordinates;
- frame/domain feature identities mismatch;
- frame rows have zero mass (unless explicitly dropped);
- pairing nodes do not exist;
- weights are nonfinite or normalization invalid;
- a query has wrong design basis;
- information is requested without valid metric units.

Record compact per-measure diagnostics: `valid`, `support_size`, `mass`, `rank_hint`, `reason`.

### Additive block kernel

Pseudocode:

```r
for (features in feature_blocks(rel)) {
  B <- relation_blocks(rel, features)        # list of q × block
  Z <- cross_atoms_svec(B, over)             # block × h

  for (rows in measurement_tiles(plan)) {
    Wb <- W[rows, features, drop = FALSE]
    for (coords in geometry_tiles(plan)) {
      tmp <- Wb %*% Z[, coords, drop = FALSE]
      add_tile(total_store, rows, coords, tmp)
    }

    for (r in partitions(rel)) {
      tmp_rel <- Wb %*% t(B[[r]])
      add_tile(local_store[[r]], rows, NULL, tmp_rel)
    }
  }
}

coherent <- cross_coherent(local, row_mass(W), over)
marginals <- pairing_marginals(local, row_mass(W), over)
new_geometry(marginals, total_store, coherent, receipt = receipt, ...)
```

`add_tile()` must avoid replacing a full output matrix. If the complete
\(m\times h\) field exceeds the in-memory result budget, `total_store` and the
coherent destination are block-backed. No neighborhood extraction loop,
prediction table, or per-view result combiner is needed.

### Determinism and memory

- Stream feature-major blocks and eliminate \(p\) early.
- Tile measurement rows and packed-geometry coordinates so sparse contraction
  never creates an unbudgeted `m × h` temporary.
- Materialize complete `m × h` geometry in memory or block-backed storage;
  direct queries return a distinct `effect_view` with `m × k` values.
- Use deterministic edge and task order, with tolerance-qualified agreement
  across block/tile choices.
- Avoid implicit parallelism in version 0.1; allow callers to parallelize participants.
- Never construct a \(p\times p\) covariance matrix in the additive path.

## 7. Correctness laws

These are release gates, not optional tests.

1. **Raw–effect equivalence:** `materialize_geometry(Y,E) == materialize_geometry(EY,I)`.
2. **Dense–streaming equivalence:** all block sizes equal dense reference.
3. **Frame equivalence:** sparse contraction equals explicit measurement loops.
4. **Query equivalence:** `effect_view` from direct evaluation equals the same
   query of complete `effect_geometry` within the declared tolerance.
5. **Trace isometry:** `sum(A * G) == crossprod(svec(A), svec(G))`.
6. **RSA equivalence:** compiled RSA equals explicit distances plus identical regression.
7. **Decomposition:** total equals coherent plus configuration.
8. **Point law:** point-frame configuration is zero.
9. **Conservation:** conservative local total geometries sum to global geometry.
10. **Feature permutation invariance:** jointly permuting relation and frame features changes nothing.
11. **Experimental covariance:** reparameterizing design and transforming queries gives the same scientific scalar.
12. **Pair orientation invariance:** swapping an undirected pair’s endpoints
    changes neither symmetrized geometry nor its endpoint marginal.
13. **Cross-null centering:** independent zero-effect partitions center at zero in simulation.
14. **No hidden PSD repair:** cross estimates preserve negative eigenvalues.
15. **No hidden transport:** incompatible domains/bases fail without explicit correspondence.
16. **Marginal semantics:** undirected pairings expose only the
    orientation-invariant endpoint marginal; directed pairings preserve role-
    specific left/right summaries.
17. **Contraction tiling:** every legal feature/row/coordinate tiling agrees
    with the dense oracle within tolerance.
18. **Storage equivalence:** in-memory and block-backed complete geometries
    answer every supported query identically within tolerance.
19. **Result honesty:** direct execution returns `effect_view`; no partial
    geometry is classed as `effect_geometry`.
20. **Version-scope enforcement:** version 0.1 rejects `workers != 1` before
    opening sources or allocating outputs.

Property-based randomized tests should exercise dimensions, sparsity, weight normalizations, block sizes, and rank-deficient extractors.

## 8. Version 0.1 scope

### Include

- dense/list/function-backed sources;
- explicit source capability metadata;
- identity and explicit matrix extractors;
- stable QR-based `lm_extractor()` with estimability diagnostics;
- abstract and volume domains;
- point, searchlight, ROI, and global additive frames;
- local and conservative normalization;
- arbitrary normalized edge pairings plus cross-partition/cross-domain constructors;
- packed total and coherent geometry in memory or block-backed storage, derived
  configuration, and pairing-appropriate marginals;
- `compute_policy(workers = 1)`, complete `materialize_geometry()`, direct
  `evaluate_geometry()`, `effect_view`, and execution receipts;
- contrast, squared-RDM, linear multiple-regression RSA;
- signed spectrum and effective-rank only when valid;
- dense and streamed kernels;
- `as_neuro()` as presentation adapter.

### Defer

- generic factor frames;
- all `workers > 1` execution until the vectorized sequential compiler is benchmarked;
- optional `shard` source/executor adapter;
- checkpoint/resume and external executor adapters;
- learned metrics and local covariance fitting;
- information summaries beyond a guarded experimental function;
- surface adapter if it expands the kernel;
- Bures or hierarchical population modeling;
- native geometry transport;
- nonlinear feature lifts;
- calibration/permutation engines;
- distributed execution and GUI.

### Exclude constitutionally

- classifier zoo and registry;
- model/result class per method;
- searchlight/regional/global modes;
- arbitrary per-location callback;
- fold objects carrying neural data;
- destructive mean removal;
- implicit image registration;
- raw-time-series design DSL;
- internal CLI/config/scheduler/save platform;
- backward compatibility with rMVPA.

## 9. Source layout and budget

```text
R/
  relation.R       # sources, relation validation, block protocol
  extractor.R      # explicit extractors, stable lm_extractor
  domain.R         # abstract/volume domain and rendering protocol
  frame.R          # sparse additive frames and normalization
  pairing.R        # canonical edge semantics
  symmetric.R      # svec/unsvec and packed queries
  kernel.R         # dense reference and tiled streamed contraction
  storage.R        # in-memory/block-backed complete fields
  execution.R      # immutable policy, receipt, reporter boundary
  geometry.R       # complete geometry and direct-view result contracts
  views.R          # contrast, RDM, RSA, guarded spectra/information
  neuroim2.R       # suggested adapter if kept outside mandatory core
```

Target executable R lines:

| Area | Budget |
|---|---:|
| Relation and sources | 280 |
| Extractors | 200 |
| Domains and frames | 450 |
| Pairings | 150 |
| Packed algebra | 180 |
| Tiled kernel | 520 |
| Geometry and storage | 300 |
| Views | 400 |
| Execution policy and receipts | 220 |
| Adapter | 180 |
| Total | 2,880 |

Do not optimize to the line count at the expense of explicit contracts. The budget is a pressure against new ontologies, not a code-golf target.

## 10. Thirty candidate decisions, critically evaluated

Scores: impact/effort/risk are 1–5; confidence is architectural confidence.

| # | Candidate | I | E | R | Conf. | Evidence | Verdict |
|---:|---|---:|---:|---:|---:|---|---|
| 1 | Lazy relation as primitive | 5 | 3 | 2 | 95% | Raw/beta equivalence and upstream estimator diversity | Keep |
| 2 | One additive sparse-frame kernel | 5 | 3 | 2 | 98% | Exact Gram pooling/collapse | Keep |
| 3 | Canonical normalized edge pairing | 5 | 2 | 2 | 95% | Removes fold/data copying and weight ambiguity | Keep |
| 4 | Packed `svec` geometry IR | 5 | 2 | 2 | 95% | Trace-isometric linear queries and compact fields | Keep |
| 5 | Store total and coherent; derive configuration | 5 | 2 | 2 | 98% | Exact frame decomposition | Keep |
| 6 | Endpoint marginal for undirected; left/right for directed | 5 | 2 | 2 | 98% | Prevents arbitrary stored-edge orientation from changing signed summaries | Keep |
| 7 | Compile squared-distance RSA | 5 | 2 | 1 | 98% | Fixed linear transform of geometry | Keep |
| 8 | Direct-query lowering returning `effect_view` | 4 | 3 | 2 | 95% | Same scalar semantics without falsely claiming complete geometry | Keep after reference |
| 9 | Conservative frame normalization | 4 | 2 | 1 | 98% | Exact local/global conservation | Keep |
| 10 | QR-based raw GLM extractor | 4 | 3 | 3 | 90% | Enables raw sources without owning GLM design | Keep narrowly |
| 11 | Function/file-backed source protocol | 4 | 3 | 2 | 90% | Out-of-core requirement | Keep |
| 12 | Dense reference implementation | 5 | 1 | 1 | 99% | Oracle for every optimized path | Keep first |
| 13 | Property/law test suite | 5 | 3 | 1 | 99% | Mathematics defines correctness | Keep first |
| 14 | Metric/unit validity metadata | 5 | 2 | 2 | 98% | Prevents false information interpretation | Keep |
| 15 | Explicit indefinite cross-estimate type | 5 | 2 | 2 | 98% | Crossvalidation legitimately yields negative eigenvalues | Keep |
| 16 | Abstract plus volume domains in v0.1 | 4 | 3 | 2 | 90% | Matrix tests plus practical fMRI entry | Keep |
| 17 | Surface domain in v0.1 | 3 | 3 | 3 | 65% | Useful but risks adapter-driven delay | Maybe |
| 18 | Generic factor frames in v0.1 | 4 | 4 | 4 | 65% | Needed for dense metrics, not for core proof | Defer |
| 19 | Gaussian information view in v0.1 | 3 | 3 | 4 | 60% | Attractive but units/PSD easily misused | Experimental/defer |
| 20 | Effective-rank view on raw cross estimates | 2 | 2 | 4 | 45% | Signed eigenvalues break standard interpretation | Reject by default |
| 21 | Learned shrinkage metric in core | 3 | 4 | 4 | 55% | Reintroduces training/hyperparameters | Defer to extension |
| 22 | Formula/event/HRF DSL | 2 | 5 | 5 | 20% | Recreates a design package | Reject |
| 23 | Classifier registry | 2 | 5 | 5 | 10% | Recreates algorithm ontology | Reject |
| 24 | Arbitrary per-location callback | 2 | 2 | 5 | 15% | Breaks semantics and optimization | Reject |
| 25 | Built-in permutation framework | 3 | 5 | 4 | 45% | Calibration is orthogonal and large | Defer/separate |
| 26 | Geometry-level linear group model in v0.1 | 3 | 4 | 3 | 60% | Geometry output is already group-ready | Defer |
| 27 | Native location transport in v0.1 | 3 | 4 | 4 | 55% | Important UEG goal, not kernel proof | Defer |
| 28 | Bures population model | 4 | 5 | 5 | 45% | Requires latent PSD and substantial inference work | Strategic later bet |
| 29 | Optional `shard` source/executor adapter | 5 | 3 | 3 | 90% | Directly addresses R worker copies while preserving one kernel | Keep after sequential benchmark |
| 30 | CLI/YAML/workbench in core | 2 | 5 | 5 | 15% | Would reproduce workflow platform | Reject |

## 11. Recommended build order

### First: finish the semantic contract and executable algebraic reference

Fix the additive theorem scope; define undirected endpoint versus directed
left/right marginals; specify complete `effect_geometry` versus direct
`effect_view`; and freeze version 0.1 as sequential-only. Then implement dense
relations, canonical pairings, `svec`, additive frames, and the 20 laws.

Success criterion: randomized dense, brute-force region, and streamed variants agree to tolerance; zero-effect cross estimates center at zero; conservative sums close numerically.

### Second: truly memory-bounded lazy streaming compiler

Add source capabilities, identity and matrix extractors, sparse `W`, tiled
feature/measurement/geometry contraction, block-backed complete geometry,
deterministic accumulation, and direct-query lowering.

Success criterion: all tilings and storage backends agree with the dense oracle;
conservative memory plans include contraction/copy/reorder/serialization
temporaries and are validated against measured peak memory.

### Third: in-process execution simulator

Exercise arbitrary completion orders, bounded inflight/reorder windows,
backpressure, failure receipts, reporter isolation, and cleanup without worker
processes.

Success criterion: adversarial schedules remain memory bounded and produce the
same tolerance-qualified result as canonical sequential execution.

### Fourth: scientific views and one neuroimaging adapter

Add contrast bundles, squared RDMs, compiled linear RSA, strict metric gating, volume frame compilation, and `as_neuro()`.

Success criterion: one vignette reproduces brute-force voxel, ROI, searchlight, and global calculations from the same relation and demonstrates total = coherent + configuration plus conservative local = global.

### Fifth: shard-backed execution, if earned

Add an optional `shard` adapter that can stage either \(Y_r\) or
\(B_r=E_rY_r\), dispatch the existing pure feature-block task, and return
bounded blocks for canonical reduction. Do not fork the geometry code.

Success criterion: ordinary, file-backed, and shard-backed sources agree with
the dense oracle; one-worker execution bypasses the pool; multiworker runs
reduce aggregate private memory/serialization or wall time on realistic data;
and success, interrupt, worker recycling, and failure all release shared
resources without stale-handle reuse.

## 12. Unknowns to resolve before scaffolding

1. Should `relation()` accept extractor matrices directly, or require `effect_extractor()` wrappers for identity/provenance?
2. Should cross-domain pairing normalize globally or separately by left/right node degrees?
3. What exact weighted squared-distance convention is canonical: total feature energy or mean feature energy?
4. Should the default frame normalization be `local` or require explicit choice? Requiring explicit choice is safer for the first release.
5. Is `neuroim2` a Suggests adapter or a separate companion package? `Matrix`-only core is cleaner.
6. What tolerance and rank policy should `lm_extractor()` use for nonestimable contrasts?
7. Which metadata are sufficient to certify noise-normalized information units without pretending to validate the upstream noise model?
8. The public package is named `crossform`. The mathematical intermediates
   retain precise names such as `effect_relation`, `effect_geometry`, and
   `effect_view`; the package name does not replace those domain terms.

## Final recommendation

Proceed with the relation-first design from Part 7, with the review corrections
now incorporated:

1. scope sparse-frame collapse to additive fixed-bilinear analyses;
2. use an endpoint marginal for undirected pairings and left/right only for
   directed pairings;
3. distinguish complete geometry from direct-query views;
4. tile feature, measurement, and packed-coordinate axes and permit
   block-backed complete output;
5. make version 0.1 sequential-only with conservative measured-memory planning;
6. defer information, dense metrics, transport, population geometry, and
   parallel adapters until the additive cross-Gram laws are executable.

The first milestone should be a mathematical reference kernel, not an R package skeleton. If that kernel stays small and every intended v0.1 analysis compiles through it, the package thesis has been demonstrated before API and adapter work begin.
