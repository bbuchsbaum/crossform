// Fused structured pair-difference kernel
//
// Evaluates sum_e w_e * (u^T L_e) * (u^T R_e) for u = e_i - e_j without
// allocating the temporary difference matrices dl, dr, or dl*dr. Pair
// products accumulate in declared edge order, then pair order; optional RSA
// coefficients contract afterwards with one BLAS product. There is no
// fast-math and no parallel reassociation.

#include <algorithm>
#include <vector>
#include <Rcpp.h>

namespace {

struct BlockView {
  const double* data;
  int n_effects;
  int n_features;
};

BlockView as_block(SEXP value, const char* what, int expected_effects,
                   int expected_features) {
  if (TYPEOF(value) != REALSXP || !Rf_isMatrix(value)) {
    Rcpp::stop("%s must be a real matrix.", what);
  }
  const int n_effects = Rf_nrows(value);
  const int n_features = Rf_ncols(value);
  if (expected_effects > 0 &&
      (n_effects != expected_effects || n_features != expected_features)) {
    Rcpp::stop("%s dimensions are inconsistent across relation blocks.", what);
  }
  BlockView block;
  block.data = REAL(value);
  block.n_effects = n_effects;
  block.n_features = n_features;
  return block;
}

std::vector<BlockView> as_blocks(const Rcpp::List& blocks, const char* what) {
  const int n_blocks = blocks.size();
  if (n_blocks < 1) {
    Rcpp::stop("%s must contain at least one relation block.", what);
  }
  std::vector<BlockView> views;
  views.reserve(static_cast<size_t>(n_blocks));
  int n_effects = 0;
  int n_features = 0;
  for (int index = 0; index < n_blocks; ++index) {
    const BlockView view = as_block(
      blocks[index], what, n_effects, n_features
    );
    if (index == 0) {
      n_effects = view.n_effects;
      n_features = view.n_features;
    }
    views.push_back(view);
  }
  return views;
}

void check_index_range(const Rcpp::IntegerVector& values, int lower, int upper,
                       const char* what) {
  for (int index = 0; index < values.size(); ++index) {
    const int value = values[index];
    if (value < lower || value > upper) {
      Rcpp::stop("%s contains an index outside the declared range.", what);
    }
  }
}

}  // namespace

// [[Rcpp::export(name = ".fused_pair_difference_atoms_cpp", rng = false)]]
Rcpp::NumericMatrix fused_pair_difference_atoms_cpp(
    const Rcpp::List& left_blocks,
    const Rcpp::List& right_blocks,
    const Rcpp::IntegerVector& left_index,
    const Rcpp::IntegerVector& right_index,
    const Rcpp::NumericVector& edge_weight,
    const Rcpp::IntegerVector& pair_left,
    const Rcpp::IntegerVector& pair_right) {
  const std::vector<BlockView> left = as_blocks(left_blocks, "left_blocks");
  const std::vector<BlockView> right = as_blocks(right_blocks, "right_blocks");
  const int n_effects = left[0].n_effects;
  const int n_features = left[0].n_features;
  if (right[0].n_effects != n_effects || right[0].n_features != n_features) {
    Rcpp::stop("Left and right relation blocks must share one shape.");
  }

  const int n_edges = left_index.size();
  if (right_index.size() != n_edges || edge_weight.size() != n_edges ||
      n_edges < 1) {
    Rcpp::stop("Edge indices and weights must be nonempty and aligned.");
  }
  check_index_range(left_index, 1, static_cast<int>(left.size()),
    "left_index");
  check_index_range(right_index, 1, static_cast<int>(right.size()),
    "right_index");

  const int n_pairs = pair_left.size();
  if (pair_right.size() != n_pairs || n_pairs < 1) {
    Rcpp::stop("Pair index vectors must be nonempty and aligned.");
  }
  check_index_range(pair_left, 1, n_effects, "pair_left");
  check_index_range(pair_right, 1, n_effects, "pair_right");

  Rcpp::NumericMatrix pair_atoms(n_features, n_pairs);
  std::fill(pair_atoms.begin(), pair_atoms.end(), 0.0);

  for (int edge = 0; edge < n_edges; ++edge) {
    const BlockView& L = left[left_index[edge] - 1];
    const BlockView& R = right[right_index[edge] - 1];
    const double weight = edge_weight[edge];
    for (int pair = 0; pair < n_pairs; ++pair) {
      const int i = pair_left[pair] - 1;
      const int j = pair_right[pair] - 1;
      const double* left_i = L.data + i;
      const double* left_j = L.data + j;
      const double* right_i = R.data + i;
      const double* right_j = R.data + j;
      double* atom_column = &pair_atoms(0, pair);
      const int stride = L.n_effects;
      for (int feature = 0; feature < n_features; ++feature) {
        const R_xlen_t offset = static_cast<R_xlen_t>(feature) * stride;
        atom_column[feature] += weight *
          (left_i[offset] - left_j[offset]) *
          (right_i[offset] - right_j[offset]);
      }
    }
  }
  return pair_atoms;
}
