// Packed and coherent full-form materialization
//
// Replaces the R-level coordinate-by-partition-edge accumulation used for
// unqueried packed atoms and coherent first-moment forms. Products accumulate
// in declared edge order, then coordinate order. There is no fast-math and
// no parallel reassociation.

#include <cmath>
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

int packed_width(int q_left, int q_right, bool packed) {
  if (packed) {
    if (q_left != q_right) {
      Rcpp::stop("Symmetric-packed atoms require a square effect form.");
    }
    return q_left * (q_left + 1) / 2;
  }
  return q_left * q_right;
}

}  // namespace

// [[Rcpp::export(name = ".packed_effect_form_atoms_cpp", rng = false)]]
Rcpp::NumericMatrix packed_effect_form_atoms_cpp(
    const Rcpp::List& left_blocks,
    const Rcpp::List& right_blocks,
    const Rcpp::IntegerVector& left_index,
    const Rcpp::IntegerVector& right_index,
    const Rcpp::NumericVector& edge_weight,
    bool packed) {
  const std::vector<BlockView> left = as_blocks(left_blocks, "left_blocks");
  const std::vector<BlockView> right = as_blocks(right_blocks, "right_blocks");
  const int q_left = left[0].n_effects;
  const int q_right = right[0].n_effects;
  const int n_features = left[0].n_features;
  if (right[0].n_features != n_features) {
    Rcpp::stop("Left and right relation blocks must share one feature axis.");
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

  const int n_coords = packed_width(q_left, q_right, packed);
  Rcpp::NumericMatrix atoms(n_features, n_coords);
  std::fill(atoms.begin(), atoms.end(), 0.0);
  const double off_diag = std::sqrt(2.0);

  int coordinate = 0;
  for (int column = 0; column < q_right; ++column) {
    const int row_start = packed ? column : 0;
    for (int row = row_start; row < q_left; ++row) {
      double* dest = &atoms(0, coordinate);
      for (int edge = 0; edge < n_edges; ++edge) {
        const BlockView& L = left[left_index[edge] - 1];
        const BlockView& R = right[right_index[edge] - 1];
        const double weight = edge_weight[edge];
        const double* left_row = L.data + row;
        const double* right_col = R.data + column;
        const int left_stride = L.n_effects;
        const int right_stride = R.n_effects;
        for (int feature = 0; feature < n_features; ++feature) {
          dest[feature] += weight *
            left_row[static_cast<R_xlen_t>(feature) * left_stride] *
            right_col[static_cast<R_xlen_t>(feature) * right_stride];
        }
      }
      if (packed && row != column) {
        for (int feature = 0; feature < n_features; ++feature) {
          dest[feature] *= off_diag;
        }
      }
      ++coordinate;
    }
  }
  return atoms;
}

// [[Rcpp::export(name = ".coherent_effect_form_atoms_cpp", rng = false)]]
Rcpp::NumericMatrix coherent_effect_form_atoms_cpp(
    const Rcpp::NumericVector& left_first,
    const Rcpp::NumericVector& right_first,
    const Rcpp::IntegerVector& left_index,
    const Rcpp::IntegerVector& right_index,
    const Rcpp::NumericVector& edge_weight,
    const Rcpp::NumericVector& mass,
    int row_start,
    int n_rows,
    bool packed) {
  Rcpp::IntegerVector left_dims = left_first.attr("dim");
  Rcpp::IntegerVector right_dims = right_first.attr("dim");
  if (left_dims.size() != 3 || right_dims.size() != 3) {
    Rcpp::stop("First-moment arrays must be measurement-by-effect-by-partition.");
  }
  const int n_meas = left_dims[0];
  const int q_left = left_dims[1];
  const int n_left_parts = left_dims[2];
  const int q_right = right_dims[1];
  const int n_right_parts = right_dims[2];
  if (right_dims[0] != n_meas) {
    Rcpp::stop("Left and right first moments must share a measurement axis.");
  }
  if (row_start < 1 || n_rows < 1 || row_start + n_rows - 1 > n_meas) {
    Rcpp::stop("Coherent tile bounds fall outside the measurement axis.");
  }
  if (mass.size() != n_meas) {
    Rcpp::stop("Coherent mass must have one value per measurement.");
  }

  const int n_edges = left_index.size();
  if (right_index.size() != n_edges || edge_weight.size() != n_edges ||
      n_edges < 1) {
    Rcpp::stop("Edge indices and weights must be nonempty and aligned.");
  }
  check_index_range(left_index, 1, n_left_parts, "left_index");
  check_index_range(right_index, 1, n_right_parts, "right_index");

  const int n_coords = packed_width(q_left, q_right, packed);
  Rcpp::NumericMatrix tile(n_rows, n_coords);
  std::fill(tile.begin(), tile.end(), 0.0);
  const double* left_ptr = REAL(left_first);
  const double* right_ptr = REAL(right_first);
  const int left_plane = n_meas * q_left;
  const int right_plane = n_meas * q_right;
  const int row0 = row_start - 1;
  const double off_diag = std::sqrt(2.0);

  int coordinate = 0;
  for (int column = 0; column < q_right; ++column) {
    const int form_row0 = packed ? column : 0;
    for (int row = form_row0; row < q_left; ++row) {
      double* dest = &tile(0, coordinate);
      for (int edge = 0; edge < n_edges; ++edge) {
        const int left_part = left_index[edge] - 1;
        const int right_part = right_index[edge] - 1;
        const double weight = edge_weight[edge];
        const double* L = left_ptr + row0 +
          static_cast<R_xlen_t>(row) * n_meas +
          static_cast<R_xlen_t>(left_part) * left_plane;
        const double* R = right_ptr + row0 +
          static_cast<R_xlen_t>(column) * n_meas +
          static_cast<R_xlen_t>(right_part) * right_plane;
        for (int local = 0; local < n_rows; ++local) {
          dest[local] += weight * L[local] * R[local];
        }
      }
      for (int local = 0; local < n_rows; ++local) {
        dest[local] /= mass[row0 + local];
        if (packed && row != column) {
          dest[local] *= off_diag;
        }
      }
      ++coordinate;
    }
  }
  return tile;
}
