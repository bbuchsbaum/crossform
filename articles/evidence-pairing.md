# Measure effects within and between neural nodes

Use this guide when your result must relate two neural measurements to
each other, rather than report one number at each. Two regions, or two
multivariate populations, are the usual cases.

A **measurement** is the same spatial unit the other guides report at: a
region, a searchlight, a voxel, or the whole brain. This guide also
calls a measurement a **node** whenever it is one endpoint of a pair,
because that is how the functions name it (`node_ids`, `edge_frame`).
The two words denote the same object.

The inputs are repeated experimental observations, a fixed measurement
operator for each node, a list of requested node pairs, and an
experimental query: a fixed matrix that weights time points, conditions,
or other experimental coordinates.
[`measurement_form()`](https://bbuchsbaum.github.io/crossform/reference/measurement_form.md)
returns one matrix block for every requested pair. A diagonal block
describes one node; an off-diagonal block relates two nodes.

The later sections apply several views to those blocks.
[`connectivity()`](https://bbuchsbaum.github.io/crossform/reference/connectivity.md)
can report an ordinary scalar correlation,
[`canonical_coupling()`](https://bbuchsbaum.github.io/crossform/reference/coupling_views.md)
reports a multivariate correlation spectrum, and
[`geometry_alignment()`](https://bbuchsbaum.github.io/crossform/reference/coupling_views.md)
compares two within-node geometries. Each view checks the assumptions it
needs. Use
[`effect_coupling()`](https://bbuchsbaum.github.io/crossform/reference/coupling_views.md)
when you need the matrix block itself without claiming that it is
covariance or connectivity.

## Can two regional signals recover ordinary correlation?

The hidden setup creates two sessions. In each matrix, eight repeated
time points are rows and four native neural features are columns. The
first two features are negatively related, which gives this example a
known correlation sign.

First declare two scalar measurements. Each one selects a single
oriented feature here; a regional mean would use one row of fixed
weights over several features.

``` r

nodes <- measurement_frame(
  list(
    anterior = matrix(c(1, 0, 0, 0), 1),
    posterior = matrix(c(0, 1, 0, 0), 1)
  ),
  domain = native,
  id = "demo:regional-means:v1"
)

node_pairs <- expand.grid(
  from = c("anterior", "posterior"),
  to = c("anterior", "posterior"),
  stringsAsFactors = FALSE
)
between_nodes <- edge_frame(node_pairs$from, node_pairs$to, nodes)
```

[`edge_frame()`](https://bbuchsbaum.github.io/crossform/reference/edge_frame.md)
computes only the pairs listed in `node_pairs`. This example requests
the two self-pairs and the two directed cross-pairs. Correlation needs
the self-pairs to obtain the two variances in its denominator. For a
seed map or selected network, list the scientific edges you need and the
self-pairs required by the chosen normalization.

The matrix `center / 7` centers the eight time points and divides their
sum of products by seven.
[`variation_query()`](https://bbuchsbaum.github.io/crossform/reference/variation_query.md)
records that operation on the time axis. The `joint_covariance`
construction is valid here because both sides of each product use the
same session and the same fitted relation.

``` r

center <- diag(8) - matrix(1 / 8, 8, 8)
sample_covariance <- variation_query(
  center / 7,
  sample_space,
  sampling_axis = "time",
  construction = "joint_covariance",
  provenance = list(estimator = "centered within session")
)
self_products <- pairing(
  signals$partitions,
  signals$partitions,
  directed = TRUE,
  self_pairs = "allow_biased",
  independence = "not_independent"
)

measured <- measurement_form(
  left = signals,
  between = between_nodes,
  by = sample_covariance,
  over = self_products
)
```

`self_pairs = "allow_biased"` permits each session to multiply by
itself; this is ordinary within-session covariance, not a crossvalidated
estimator. `independence = "not_independent"` prevents the result from
claiming that the two sides were independently estimated.
[`measurement_form()`](https://bbuchsbaum.github.io/crossform/reference/measurement_form.md)
averages the two session products and returns the four requested blocks.

Because both node measurements are scalar and oriented, the correlation
view returns ordinary signed Pearson correlations.

``` r

functional <- connectivity(measured, view = "correlation")
data.frame(
  from = measured$block_index$left,
  to = measured$block_index$right,
  correlation = functional$values$correlation
)
#>        from        to correlation
#> 1  anterior  anterior   1.0000000
#> 2 posterior  anterior  -0.9912024
#> 3  anterior posterior  -0.9912024
#> 4 posterior posterior   1.0000000
```

The negative cross-node value matches the generated data. Its sign is
defined because each scalar measurement has a fixed orientation. A
multivariate node can be rotated without changing its represented
subspace, so individual cross-block entries are not invariant there. The
next section uses summaries that do not depend on an arbitrary basis
rotation.

## What changes when each node retains several modes?

Use two fixed feature selectors as multivariate nodes, then request the
same four blocks under the same experimental query.

``` r

populations <- measurement_frame(
  list(
    anterior = diag(4)[1:2, , drop = FALSE],
    posterior = diag(4)[3:4, , drop = FALSE]
  ),
  domain = native,
  id = "demo:multivariate-populations:v1"
)
population_pairs <- edge_frame(
  node_pairs$from, node_pairs$to, populations
)
population_form <- measurement_form(
  signals, population_pairs, sample_covariance, self_products
)
```

Canonical coupling uses the two self-covariance blocks to scale the
cross-node block, then reports its singular values. These canonical
correlations describe the shared modes in descending order.
`ridge = 0.05` stabilizes the two inverse covariance calculations and is
recorded because changing it changes the reported values.

``` r

canonical <- canonical_coupling(population_form, ridge = 0.05)
canonical$values[
  canonical$values$edge_id == cross_population_edge,
  c("mode", "canonical_correlation")
]
#>   mode canonical_correlation
#> 5    1            0.93961123
#> 6    2            0.06163478
```

[`geometry_alignment()`](https://bbuchsbaum.github.io/crossform/reference/coupling_views.md)
answers a different question: do the two populations induce similar
pairwise geometry over the repeated observations? It computes a linear
CKA/RV-like alignment. Under the joint-covariance construction used
here, the result lies between zero and one.

``` r

alignment <- geometry_alignment(population_form)
alignment$values$geometry_alignment[
  alignment$values$edge_id == cross_population_edge
]
#> [1] 0.9413427
```

Gaussian mutual information treats the normalized canonical modes as
arising from a joint Gaussian model. The call therefore requires an
explicit model declaration and a choice of information units.

``` r

gaussian_model <- gaussian_covariance_model(
  list(assumption = "joint Gaussian time observations")
)
information <- connectivity(
  population_form,
  view = "gaussian_information",
  ridge = 0.05,
  model = gaussian_model,
  units = "bits"
)
information$values[
  information$values$edge_id == cross_population_edge,
  c("information", "units")
]
#>   information units
#> 3     1.54965  bits
```

## Why can one contrast not estimate connectivity?

A contrast vector $`c`$ defines the rank-one query $`H=cc^\top`$. This
query can measure whether two nodes express the same contrast, but it
leaves only one experimental direction. Its neural block is an outer
product, so a normalized correlation would equal one whenever both node
effects are nonzero. The repeated variation needed to estimate
connectivity has been removed.

``` r

effect_direction <- c(-1, -1, -1, -1, 1, 1, 1, 1)
rank_one_query <- variation_query(
  tcrossprod(effect_direction),
  sample_space,
  sampling_axis = "time",
  construction = "joint_covariance"
)
rank_one_form <- measurement_form(
  signals, between_nodes, rank_one_query, self_products
)
rank_one_effect <- effect_coupling(rank_one_form)
c(
  effective_rank = rank_one_form$diagnostics$experimental_effective_rank,
  edge_blocks = length(rank_one_effect$values)
)
#> effective_rank    edge_blocks 
#>              1              4
```

[`effect_coupling()`](https://bbuchsbaum.github.io/crossform/reference/coupling_views.md)
returns the block because it makes no covariance claim.
`connectivity(..., view = "correlation")` stops with an error because
the experimental query has effective rank one. A matrix with two neural
endpoints is therefore not sufficient evidence for connectivity.

## How do coherent and configuration modes cross on an edge?

An additive measurement has a positively oriented weighted-mean
direction, called its coherent mode. The orthogonal remainder is its
configuration subspace. When both nodes use this decomposition, their
cross-node block has four components: coherent-to-coherent,
coherent-to-configuration, configuration-to-coherent, and
configuration-to-configuration.

``` r

additive <- additive_frame(
  matrix(c(
    1, 2, 1, 0,
    0, 1, 2, 1
  ), 2, 4, byrow = TRUE),
  domain = native
)
decomposed_nodes <- measurement_frame(
  additive, mode = "coherent_configuration"
)
decomposed_edge <- edge_frame(
  decomposed_nodes$node_ids[1],
  decomposed_nodes$node_ids[2],
  decomposed_nodes
)
decomposed_form <- measurement_form(
  signals, decomposed_edge, sample_covariance, self_products
)
components <- measurement_components(decomposed_form, edge = 1)
components[, c(
  "left_component", "right_component", "raw_entries_meaningful",
  "frobenius_strength"
)]
#>   left_component right_component raw_entries_meaningful frobenius_strength
#> 1       coherent        coherent                   TRUE          0.4762249
#> 2  configuration        coherent                  FALSE          1.7337607
#> 3       coherent   configuration                  FALSE          0.3089121
#> 4  configuration   configuration                  FALSE          1.6400838
```

The coherent direction has a fixed orientation because the measurement
weights define its positive direction. The configuration subspace is
fixed, but any orthonormal basis within it is interchangeable.
[`measurement_components()`](https://bbuchsbaum.github.io/crossform/reference/measurement_components.md)
therefore reports rotation-invariant Frobenius strengths for components
that involve configuration. Frobenius strength is the square root of the
sum of the block’s squared entries. The function does not present
individual configuration-to-configuration entries as uniquely defined
quantities.

## How do encoding and retrieval use different axes?

[`measurement_form()`](https://bbuchsbaum.github.io/crossform/reference/measurement_form.md)
accepts separate left and right relations. In this example, the left
relation has two encoding conditions over three neural features. The
right relation has three retrieval conditions over four neural features.
[`pair_query()`](https://bbuchsbaum.github.io/crossform/reference/pair_query.md)
therefore uses a rectangular 2 by 3 experimental query, while
[`edge_frame()`](https://bbuchsbaum.github.io/crossform/reference/edge_frame.md)
relates a scalar encoding seed to a two-mode retrieval target.

``` r

encoding_nodes <- measurement_frame(
  list(seed = matrix(c(1, 0, 1), 1)),
  encoding_domain,
  id = "demo:encoding-node:v1"
)
retrieval_nodes <- measurement_frame(
  list(target = matrix(c(1, 0, 0, 1, 0, 1, 1, 0), 2, 4)),
  retrieval_domain,
  id = "demo:retrieval-node:v1"
)
cross_node <- edge_frame(
  "seed", "target", encoding_nodes, to_frame = retrieval_nodes
)
cross_query <- pair_query(
  matrix(c(1, -0.5, 0, 0, 0.5, -1), 2, 3, byrow = TRUE),
  encoding_space,
  retrieval_space
)
cross_form <- measurement_form(
  left = encoding,
  right = retrieval,
  between = cross_node,
  by = cross_query,
  over = pairing("encoding_run", "retrieval_run", directed = TRUE)
)
effect_coupling(cross_form)$values[[1]]
#>      [,1] [,2]
#> [1,]    6   -3
```

The returned matrix is 1 by 2. Its rows belong to the encoding seed and
its columns belong to the retrieval target. No square RDM is constructed
because the two experimental axes and the two neural measurements are
different.

## Did the node and edge blocks retain all neural evidence?

[`reconstruct_evidence()`](https://bbuchsbaum.github.io/crossform/reference/reconstruct_evidence.md)
tests whether the measured blocks retain the complete global neural
operator or only a projection of it. Lossless reconstruction requires
every directed node pair and a measurement frame with full column rank.
This example uses a Parseval frame: stacking the two node operators
gives the identity, so the local blocks can be added back without a
correction matrix.

``` r

parseval_nodes <- measurement_frame(
  list(
    first_half = diag(4)[1:2, , drop = FALSE],
    second_half = diag(4)[3:4, , drop = FALSE]
  ),
  native,
  id = "demo:parseval:v1"
)
parseval_pairs <- expand.grid(
  from = parseval_nodes$node_ids,
  to = parseval_nodes$node_ids,
  stringsAsFactors = FALSE
)
parseval_edges <- edge_frame(
  parseval_pairs$from, parseval_pairs$to, parseval_nodes
)
parseval_form <- measurement_form(
  signals, parseval_edges, sample_covariance, self_products
)
reconstructed <- reconstruct_evidence(parseval_form, parseval_edges)
c(
  method = reconstructed$method,
  status = reconstructed$status,
  lossless = reconstructed$lossless
)
#>                           method                           status 
#>                       "parseval" "exact_algebraic_reconstruction" 
#>                         lossless 
#>                           "TRUE"
```

A full-rank non-Parseval frame reconstructs through its canonical dual,
a correction that reverses overlapping or non-orthogonal frame weights.
A rank-deficient frame returns `status = "projected_reconstruction"` and
`lossless = FALSE`, which tells you that the measurements discarded
neural directions. The function refuses a lossless claim when edges are
missing, frame bases disagree, the reconstruction is too
ill-conditioned, or the declared workspace is too small. Diagonal node
blocks alone are insufficient: two global operators can share the same
diagonal blocks and differ between nodes.

## Which interpretation is permitted?

| View | Minimum interpretation contract |
|----|----|
| [`effect_coupling()`](https://bbuchsbaum.github.io/crossform/reference/coupling_views.md) | Any supported fixed experimental query; no covariance claim |
| [`covariance_coupling()`](https://bbuchsbaum.github.io/crossform/reference/coupling_views.md) | Repeated variation, effective rank above one, a joint-covariance construction, and valid self-blocks |
| `connectivity(..., "correlation")` | The covariance requirements plus scalar, oriented measurement axes |
| [`canonical_coupling()`](https://bbuchsbaum.github.io/crossform/reference/coupling_views.md) | The covariance requirements plus explicit ridge regularization |
| [`geometry_alignment()`](https://bbuchsbaum.github.io/crossform/reference/coupling_views.md) | Joint covariance and nonzero self-geometries; static normalization over form entries |
| Gaussian information | Canonical coupling plus an explicit joint Gaussian model and information units |
| [`reconstruct_evidence()`](https://bbuchsbaum.github.io/crossform/reference/reconstruct_evidence.md) | Every directed edge, matching frame bases, sufficient workspace, and acceptable frame rank and conditioning |

Each result records which axis was normalized. Pattern correlation
normalizes neural features before it produces an experimental form.
Functional correlation normalizes repeated experimental samples before
it interprets a neural form. The two operations both use the word
correlation, but they act on different axes and answer different
questions.

The measurement workflow first averages the requested partition products
and then applies a connectivity normalization. Normalizing or
Fisher-transforming each partition pair before averaging would define a
different estimator, so the package does not substitute that order. The
result identity records the source order, view order, regularization,
units, and sampling axis.

## Is geometry alignment informational connectivity?

No. [`geometry_alignment()`](https://bbuchsbaum.github.io/crossform/reference/coupling_views.md)
compares static Gram matrices and is a rotation-invariant multivariate
generalization of squared scalar correlation. Informational-connectivity
analyses often correlate dynamic discriminability series estimated
across time or trials instead. Estimating those series requires training
and cross-fitting stages that this fixed-form package does not
implement.

Choose the function by the question:

- use
  [`effect_coupling()`](https://bbuchsbaum.github.io/crossform/reference/coupling_views.md)
  for an arbitrary fixed query;
- use
  [`connectivity()`](https://bbuchsbaum.github.io/crossform/reference/connectivity.md)
  only when repeated variation and covariance capabilities are present;
- call static CKA/RV-like comparison
  [`geometry_alignment()`](https://bbuchsbaum.github.io/crossform/reference/coupling_views.md);
- keep dynamic informational-connectivity estimators in an adjacent
  system until their training, cross-fitting, and inference contracts
  are explicit.

For contrasts, RDMs, RSA, and geometry spectra within one experimental
space, continue with
[`vignette("introduction", package = "crossform")`](https://bbuchsbaum.github.io/crossform/articles/introduction.md).
Use the measurement-form functions in this guide only when the returned
object must retain explicit neural node and edge blocks.
