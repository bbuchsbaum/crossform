// Grouped-partner Gram kernel for support-streamed pair differences
//
// The R route in `.support_streamed_metric_contraction()` evaluates, at every
// spatial measurement and for every ordered partition edge, one row-difference
// product per requested effect pair. When the query asks for most of the pairs
// that is the wrong association: the same numbers are the entries of one small
// effect-by-effect matrix per measurement.
//
// With support s of size n, local metric K (already composed with the frame
// weights by the caller), partition blocks L_a (q by n), and the half-edge
// weight w_ab of the unordered partition pair {a, b}:
//
//   T   = sum over ordered edges of w_e L_{a_e} K L_{b_e}'
//       = sum_{a<b} w_ab ( L_a K L_b' + L_b K L_a' )
//   value(i, j) = T_ii + T_jj - T_ij - T_ji
//
// The reduction groups the inner sum by its left partition:
//
//   H   = sum_{a<b} w_ab L_a K L_b'  =  sum_{a=1}^{P-1} L_a K C_a',
//   C_a = sum_{b>a} w_ab L_b,        T = H + H'
//
// which is P-1 products of rank n instead of P(P-1)/2, and P-1 metric products
// instead of P(P-1). The grouping is a re-association of the estimand's own
// sum: nothing is subtracted from a larger quantity of the same sign, so it
// introduces no cancellation the estimand does not already carry. The partner
// sums C_a do combine signed partition blocks, but that is the same summation
// the edge sum they stand for performs. It is deliberately not the collapse
// T = w(S K S' - sum_a L_a K L_a'), which is cheaper still but subtracts a
// within-partition term from a larger between-partition one, and would need an
// accuracy argument this reduction does not.
//
// With no declared metric the local operator is not absent, it is diagonal:
// the frame weights ARE the operator, K = diag(w), because composing the
// identity metric with the frame gives K_ij = delta_ij r_i r_j = delta_ij w_i.
// (It is not the rank-one r r'.) `diagonal` selects that reading, and then
// L_a K is a column scaling rather than a product, which removes the whole
// first stage. The weight multiplies once, on the left factor, which is the
// association the additive route it replaces already uses: it accumulates
// sum_v w_v (dl_v)(dr_v), one weight per feature per term.
//
// T = H + H' holds because a crossform metric is stored exactly symmetric:
// `neural_metric()` returns `(value + t(value))/2` after admitting it within
// its declared tolerance (R/metric.R, `.canonical_symmetric_metric()`), and the
// frame weights compose as K_ij = value_ij r_i r_j, so the composed local
// metric is symmetric bit for bit. Nothing here factors K, so the route is
// metric-agnostic: an indefinite or singular local metric is contracted, not
// refused, exactly as the reference loop does.
//
// The C_a are built once for the whole feature axis, so the per-measurement
// work is a gather and four matrix products. Accumulation order is fixed
// (declared edge order, then partitions ascending, then pairs), there is no
// fast-math and no parallel reassociation, so the result is reproducible for a
// given BLAS.
//
// The caller owns every admission decision; this translation unit assumes its
// inputs are already checked.

#ifndef USE_FC_LEN_T
# define USE_FC_LEN_T
#endif

#include <cstring>
#include <vector>
#include <Rcpp.h>
#include <R_ext/BLAS.h>

#ifndef FCONE
# define FCONE
#endif

namespace {

const char kNoTranspose = 'N';
const char kTranspose = 'T';
const char kLower = 'L';
const double kOne = 1.0;
const double kZero = 0.0;

// Two layouts of the same numbers. `sweep` keeps every partition's slab
// contiguous down one feature column, so a support gathers with one memcpy per
// feature and the metric product is a single call over all partitions at once.
// `paired` splits the gather into per-partition strips so the collapse can be
// one symmetric rank-2k update of inner dimension (P-1)n, which runs closer to
// peak once q is large enough to pay for the wider gather.
const int kSweep = 1;
const int kPaired = 2;

}  // namespace

// [[Rcpp::export(name = ".support_metric_gram_pairs_cpp", rng = false)]]
Rcpp::List support_metric_gram_pairs_cpp(
    const Rcpp::IntegerVector& row_ptr,
    const Rcpp::IntegerVector& column_index,
    const Rcpp::NumericVector& frame_weight,
    const Rcpp::NumericMatrix& metric_value,
    const Rcpp::IntegerVector& metric_row,
    const Rcpp::List& blocks,
    const Rcpp::IntegerVector& edge_left,
    const Rcpp::IntegerVector& edge_right,
    const Rcpp::NumericVector& edge_weight,
    const Rcpp::IntegerVector& pair_left,
    const Rcpp::IntegerVector& pair_right,
    int path,
    bool diagonal) {
  const int n_nodes = row_ptr.size() - 1;
  const int n_partitions = blocks.size();
  const int n_pairs = pair_left.size();
  const int n_edges = edge_left.size();
  const int n_metric = metric_value.nrow();
  if (n_nodes < 1 || n_partitions < 2 || n_pairs < 1 || n_edges < 1) {
    Rcpp::stop("The Gram route requires nodes, partitions, edges and pairs.");
  }

  // The coerced views are retained, not just their data pointers: a block that
  // arrives as anything but a real matrix is converted here and must stay
  // alive while the staging matrix is built from it.
  std::vector<Rcpp::NumericMatrix> held;
  held.reserve(static_cast<std::size_t>(n_partitions));
  std::vector<const double*> block(n_partitions);
  int q = 0;
  int n_features = 0;
  for (int a = 0; a < n_partitions; ++a) {
    held.push_back(Rcpp::as<Rcpp::NumericMatrix>(blocks[a]));
    const Rcpp::NumericMatrix& value = held.back();
    if (a == 0) {
      q = value.nrow();
      n_features = value.ncol();
    } else if (value.nrow() != q || value.ncol() != n_features) {
      Rcpp::stop("Partition blocks must share one effect-by-feature shape.");
    }
    block[a] = &value[0];
  }
  if (q < 2 || n_features < 1) {
    Rcpp::stop("Partition blocks must have at least two effects.");
  }

  if (!diagonal) {
    if (metric_row.size() != n_features || n_metric < 1 ||
        metric_value.ncol() != n_metric) {
      Rcpp::stop("The metric row map must cover every declared feature.");
    }
    for (int f = 0; f < n_features; ++f) {
      if (metric_row[f] < 0 || metric_row[f] >= n_metric) {
        Rcpp::stop("A metric row index falls outside the metric operator.");
      }
    }
  }
  if (column_index.size() != frame_weight.size() ||
      row_ptr[0] != 0 || row_ptr[n_nodes] != column_index.size()) {
    Rcpp::stop("Compressed frame weights are not canonically shaped.");
  }
  for (R_xlen_t k = 0; k < column_index.size(); ++k) {
    if (column_index[k] < 0 || column_index[k] >= n_features) {
      Rcpp::stop("A frame support falls outside the relation feature axis.");
    }
  }
  if (edge_right.size() != n_edges || edge_weight.size() != n_edges) {
    Rcpp::stop("Edge endpoints and weights must be aligned.");
  }
  for (int e = 0; e < n_edges; ++e) {
    if (edge_left[e] < 0 || edge_right[e] < 0 ||
        edge_left[e] >= n_partitions || edge_right[e] >= n_partitions ||
        edge_left[e] >= edge_right[e]) {
      Rcpp::stop("Edges must name two ascending partitions of the family.");
    }
  }
  if (pair_right.size() != n_pairs) {
    Rcpp::stop("Pair index vectors must be aligned.");
  }
  for (int p = 0; p < n_pairs; ++p) {
    if (pair_left[p] < 0 || pair_left[p] >= q ||
        pair_right[p] < 0 || pair_right[p] >= q) {
      Rcpp::stop("A pair names an effect outside the declared space.");
    }
  }
  if (path != kSweep && path != kPaired) {
    Rcpp::stop("The requested Gram layout is not one of the two available.");
  }

  int max_support = 0;
  for (int node = 0; node < n_nodes; ++node) {
    const int m = row_ptr[node + 1] - row_ptr[node];
    if (m < 1) {
      Rcpp::stop("Every spatial measurement must have a nonempty support.");
    }
    if (m > max_support) max_support = m;
  }

  // The staging matrix: [L_0 .. L_{A-1} ; C_0 .. C_{A-1}] down each feature
  // column, so one support feature is one contiguous run of 2 A q doubles.
  // Only a partition that is the smaller endpoint of some edge needs an L
  // slab, and the largest such index is P - 2, so A = P - 1 covers both.
  const int left_slots = n_partitions - 1;
  const int stage_rows = 2 * left_slots * q;
  const std::size_t left_offset = static_cast<std::size_t>(left_slots) * q;
  std::vector<double> stage(
    static_cast<std::size_t>(stage_rows) * n_features, 0.0
  );
  for (int f = 0; f < n_features; ++f) {
    double* target = &stage[static_cast<std::size_t>(f) * stage_rows];
    const std::size_t column = static_cast<std::size_t>(f) * q;
    for (int a = 0; a < left_slots; ++a) {
      std::memcpy(target + static_cast<std::size_t>(a) * q, block[a] + column,
        sizeof(double) * static_cast<std::size_t>(q));
    }
    for (int e = 0; e < n_edges; ++e) {
      const double weight = edge_weight[e];
      double* partner_slab = target + left_offset +
        static_cast<std::size_t>(edge_left[e]) * q;
      const double* source = block[edge_right[e]] + column;
      for (int i = 0; i < q; ++i) partner_slab[i] += weight * source[i];
    }
  }

  // Offsets into the staged per-measurement matrix, resolved once. The sweep
  // layout leaves a full matrix and reads both off-diagonal entries; the
  // paired layout leaves the symmetric combination in one triangle.
  const std::size_t gram_size = static_cast<std::size_t>(q) * q;
  std::vector<std::size_t> offset_left(n_pairs);
  std::vector<std::size_t> offset_right(n_pairs);
  std::vector<std::size_t> offset_cross(n_pairs);
  std::vector<std::size_t> offset_mirror(n_pairs);
  for (int p = 0; p < n_pairs; ++p) {
    const int i = pair_left[p];
    const int j = pair_right[p];
    const int lower = (i <= j) ? i : j;
    const int upper = (i <= j) ? j : i;
    offset_left[p] = static_cast<std::size_t>(i) * q + i;
    offset_right[p] = static_cast<std::size_t>(j) * q + j;
    offset_cross[p] = static_cast<std::size_t>(j) * q + i;
    offset_mirror[p] = (path == kPaired)
      ? static_cast<std::size_t>(lower) * q + upper
      : static_cast<std::size_t>(i) * q + j;
  }

  // The output is measurement-major, so one measurement's pairs land one
  // column apart. Filling a block of measurements before writing turns that
  // scatter into contiguous runs; the block is sized to keep the staged
  // matrices inside a mid-level cache. Nothing about the arithmetic changes:
  // each output entry still comes from one measurement's own matrix.
  const std::size_t block_budget = 1u << 20;
  int node_block = static_cast<int>(block_budget / (gram_size * sizeof(double)));
  if (node_block < 1) node_block = 1;
  if (node_block > 64) node_block = 64;
  if (node_block > n_nodes) node_block = n_nodes;

  Rcpp::NumericMatrix out(n_nodes, n_pairs);
  Rcpp::NumericVector mass(n_nodes);
  const int whitened_rows = left_slots * q;
  const std::size_t whitened_size =
    static_cast<std::size_t>(whitened_rows) * max_support;
  std::vector<double> gathered(
    static_cast<std::size_t>(stage_rows) * max_support
  );
  std::vector<double> metric_block(
    diagonal ? 0u : static_cast<std::size_t>(max_support) * max_support
  );
  std::vector<double> whitened(diagonal ? 0u : whitened_size);
  std::vector<double> partner(path == kPaired ? whitened_size : 0u);
  std::vector<double> staged(static_cast<std::size_t>(node_block) * gram_size);
  std::vector<double> root_weight(max_support);
  std::vector<int> support_column(column_index.size());
  for (R_xlen_t k = 0; k < column_index.size(); ++k) {
    support_column[k] = column_index[k];
  }
  std::vector<int> metric_of_feature(diagonal ? 0 : n_features);
  for (int f = 0; f < static_cast<int>(metric_of_feature.size()); ++f) {
    metric_of_feature[f] = metric_row[f];
  }

  const double* metric_data = diagonal ? NULL : &metric_value[0];
  double* out_data = &out[0];

  for (int node = 0; node < n_nodes; ++node) {
    double* gram = &staged[
      static_cast<std::size_t>(node % node_block) * gram_size
    ];
    const int start = row_ptr[node];
    const int m = row_ptr[node + 1] - start;
    const int* support = &support_column[start];
    const double* weight = &frame_weight[start];

    double node_mass = 0.0;
    for (int i = 0; i < m; ++i) node_mass += weight[i];
    mass[node] = node_mass;

    if (!diagonal) {
      for (int i = 0; i < m; ++i) root_weight[i] = std::sqrt(weight[i]);
      // K = metric[support, support] composed with the frame weight roots.
      for (int j = 0; j < m; ++j) {
        const std::size_t column =
          static_cast<std::size_t>(metric_of_feature[support[j]]) * n_metric;
        const double scale_j = root_weight[j];
        double* target = &metric_block[static_cast<std::size_t>(j) * m];
        for (int i = 0; i < m; ++i) {
          target[i] = metric_data[column + metric_of_feature[support[i]]] *
            root_weight[i] * scale_j;
        }
      }
    }

    if (path == kSweep) {
      for (int j = 0; j < m; ++j) {
        std::memcpy(
          &gathered[static_cast<std::size_t>(j) * stage_rows],
          &stage[static_cast<std::size_t>(support[j]) * stage_rows],
          sizeof(double) * static_cast<std::size_t>(stage_rows)
        );
      }
      const double* left_factor;
      int left_stride;
      if (diagonal) {
        // L_a K is a column scaling. It is applied in place to the gathered
        // left slabs, which the collapse then reads at the staging stride.
        for (int j = 0; j < m; ++j) {
          double* column = &gathered[static_cast<std::size_t>(j) * stage_rows];
          const double scale = weight[j];
          for (int i = 0; i < whitened_rows; ++i) column[i] *= scale;
        }
        left_factor = gathered.data();
        left_stride = stage_rows;
      } else {
        // Y = [L_0; ..; L_{A-1}] K, one product over every partition at once.
        F77_CALL(dgemm)(&kNoTranspose, &kNoTranspose, &whitened_rows, &m, &m,
          &kOne, gathered.data(), &stage_rows, metric_block.data(), &m,
          &kZero, whitened.data(), &whitened_rows FCONE FCONE);
        left_factor = whitened.data();
        left_stride = whitened_rows;
      }
      // H = sum_a Y_a C_a'.
      for (int a = 0; a < left_slots; ++a) {
        const double beta = (a == 0) ? kZero : kOne;
        F77_CALL(dgemm)(&kNoTranspose, &kTranspose, &q, &q, &m, &kOne,
          left_factor + static_cast<std::size_t>(a) * q, &left_stride,
          &gathered[left_offset + static_cast<std::size_t>(a) * q],
          &stage_rows, &beta, gram, &q FCONE FCONE);
      }
    } else {
      // One strip per partition, so the collapse is a single symmetric
      // rank-2k update of inner dimension A n.
      for (int u = 0; u < m; ++u) {
        const double* source =
          &stage[static_cast<std::size_t>(support[u]) * stage_rows];
        for (int a = 0; a < left_slots; ++a) {
          const std::size_t strip = (static_cast<std::size_t>(a) * m + u) * q;
          std::memcpy(&gathered[strip],
            source + static_cast<std::size_t>(a) * q,
            sizeof(double) * static_cast<std::size_t>(q));
          std::memcpy(&partner[strip],
            source + left_offset + static_cast<std::size_t>(a) * q,
            sizeof(double) * static_cast<std::size_t>(q));
        }
      }
      const double* left_factor;
      if (diagonal) {
        for (int a = 0; a < left_slots; ++a) {
          for (int u = 0; u < m; ++u) {
            double* column = &gathered[(static_cast<std::size_t>(a) * m + u) * q];
            const double scale = weight[u];
            for (int i = 0; i < q; ++i) column[i] *= scale;
          }
        }
        left_factor = gathered.data();
      } else {
        for (int a = 0; a < left_slots; ++a) {
          const std::size_t strip = static_cast<std::size_t>(a) * m * q;
          F77_CALL(dgemm)(&kNoTranspose, &kNoTranspose, &q, &m, &m, &kOne,
            &gathered[strip], &q, metric_block.data(), &m,
            &kZero, &whitened[strip], &q FCONE FCONE);
        }
        left_factor = whitened.data();
      }
      const int inner = left_slots * m;
      F77_CALL(dsyr2k)(&kLower, &kNoTranspose, &q, &inner, &kOne,
        left_factor, &q, partner.data(), &q, &kZero, gram, &q
        FCONE FCONE);
    }

    // Flush a full block of staged matrices, or the tail at the end.
    const bool flush = (node % node_block == node_block - 1) ||
      (node == n_nodes - 1);
    if (!flush) continue;
    const int first = node - (node % node_block);
    const int count = node - first + 1;
    for (int p = 0; p < n_pairs; ++p) {
      const std::size_t left = offset_left[p];
      const std::size_t right = offset_right[p];
      const std::size_t cross = offset_cross[p];
      const std::size_t mirror = offset_mirror[p];
      double* column = out_data + static_cast<std::size_t>(p) * n_nodes + first;
      if (path == kPaired) {
        // The staged matrix is already T = H + H'.
        for (int b = 0; b < count; ++b) {
          const double* value =
            &staged[static_cast<std::size_t>(b) * gram_size];
          const double cross_value = value[mirror];
          column[b] = value[left] + value[right] - cross_value - cross_value;
        }
      } else {
        // The staged matrix is H alone, and T = H + H' contributes each of
        // the four entries exactly twice.
        for (int b = 0; b < count; ++b) {
          const double* value =
            &staged[static_cast<std::size_t>(b) * gram_size];
          column[b] = 2 * (value[left] + value[right] -
            value[cross] - value[mirror]);
        }
      }
    }
  }

  return Rcpp::List::create(
    Rcpp::Named("value") = out,
    Rcpp::Named("mass") = mass,
    Rcpp::Named("max_support") = max_support,
    Rcpp::Named("node_block") = node_block,
    Rcpp::Named("stage_rows") = stage_rows
  );
}
