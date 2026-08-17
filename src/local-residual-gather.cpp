// Gather a dense local residual covariance from the canonical upper CSC
// pair graph. Values accumulate in support-column order, then the upper
// triangle of that column. There is no fast-math and no parallel
// reassociation.

#include <Rcpp.h>

// [[Rcpp::export(name = ".local_residual_gather_cpp", rng = false)]]
Rcpp::NumericMatrix local_residual_gather_cpp(
    const Rcpp::IntegerVector& column_ptr,
    const Rcpp::IntegerVector& row_index,
    const Rcpp::NumericVector& values,
    const Rcpp::IntegerVector& support) {
  const int n_features = column_ptr.size() - 1;
  const int n_pairs = row_index.size();
  if (n_features < 1 || values.size() != n_pairs ||
      column_ptr[0] != 0 || column_ptr[n_features] != n_pairs) {
    Rcpp::stop("Local covariance extraction requires the canonical upper pair graph.");
  }
  const int dimension = support.size();
  if (dimension < 1) {
    Rcpp::stop("A local residual support must contain at least one feature.");
  }
  for (int index = 0; index < dimension; ++index) {
    if (support[index] < 1 || support[index] > n_features) {
      Rcpp::stop("A support index falls outside the pair graph.");
    }
  }

  Rcpp::NumericMatrix value(dimension, dimension);
  std::fill(value.begin(), value.end(), 0.0);
  for (int column = 1; column < dimension; ++column) {
    if (support[column] <= support[column - 1]) {
      Rcpp::stop("A local residual support must be strictly increasing.");
    }
  }
  for (int column = 0; column < dimension; ++column) {
    const int global_column = support[column];
    const int first = column_ptr[global_column - 1];
    const int last = column_ptr[global_column];
    if (first < 0 || last < first || last > n_pairs) {
      Rcpp::stop("The support pair graph has invalid canonical coordinates.");
    }
    int slot = first;
    for (int local_row = 0; local_row <= column; ++local_row) {
      const int desired = support[local_row] - 1;
      while (slot < last && row_index[slot] < desired) {
        ++slot;
      }
      if (slot >= last || row_index[slot] != desired) {
        Rcpp::stop("A support pair is absent from its declared union pair graph.");
      }
      const double pair_value = values[slot];
      value(local_row, column) = pair_value;
      value(column, local_row) = pair_value;
    }
  }
  return value;
}
