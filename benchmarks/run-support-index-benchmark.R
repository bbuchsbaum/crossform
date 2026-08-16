# Realistic support-topology benchmark. Run from the package root with:
# Rscript benchmarks/run-support-index-benchmark.R

pkgload::load_all(quiet = TRUE)

dims <- c(48L, 48L, 24L)
spacing <- c(2, 2, 2)
radius <- 6.01
domain <- volume_domain(
  array(TRUE, dims), spacing = spacing, id = "support-benchmark-55k"
)

elapsed <- system.time({
  index <- crossform:::.euclidean_support_index(domain, radius)
})
evaluation_edges <- 8L
preflight <- crossform:::.support_index_preflight(
  index, 4 * 1024^3, evaluation_edges = evaluation_edges
)

record <- data.frame(
  features = domain$n_features,
  radius_mm = radius,
  support_min = index$cost$support_size[["min"]],
  support_median = index$cost$support_size[["median"]],
  support_mean = index$cost$support_size[["mean"]],
  support_max = index$cost$support_size[["max"]],
  union_degree_mean = index$cost$union_degree[["mean"]],
  union_degree_max = index$cost$union_degree[["max"]],
  memberships = index$cost$support_memberships,
  pair_pattern_nnz = index$cost$pair_pattern_nnz,
  structural_mib = index$cost$estimated_structural_bytes / 1024^2,
  evaluation_edges = preflight$evaluation_edges,
  dense_schedule_per_edge_gib =
    preflight$materialized_dense_metric_bytes_per_edge / 1024^3,
  dense_schedule_total_gib =
    preflight$materialized_dense_metric_bytes_total / 1024^3,
  one_factorization_pass_per_edge =
    preflight$one_dense_factorization_pass_units_per_edge,
  one_factorization_pass_total =
    preflight$one_dense_factorization_pass_units_total,
  factorization_passes = preflight$factorization_passes,
  metric_schedule_storage = preflight$metric_schedule_storage,
  dense_lowering = preflight$dense_metric_lowering,
  elapsed_seconds = unname(elapsed[["elapsed"]]),
  stringsAsFactors = FALSE
)

coordinates <- as.matrix(expand.grid(
  x = seq_len(37L), y = seq_len(37L), z = seq_len(37L)
))
coordinate_domain <- abstract_domain(
  nrow(coordinates), coordinates = coordinates,
  id = "support-benchmark-50k-coordinates"
)
coordinate_elapsed <- system.time({
  coordinate_index <- crossform:::.euclidean_support_index(
    coordinate_domain, 1.1
  )
})

print(record, row.names = FALSE)
cat("\nIrregular-coordinate provider\n")
print(data.frame(
  features = coordinate_domain$n_features,
  lookup = coordinate_index$construction$lookup,
  support_mean = coordinate_index$cost$support_size[["mean"]],
  pair_pattern_nnz = coordinate_index$cost$pair_pattern_nnz,
  structural_mib =
    coordinate_index$cost$estimated_structural_bytes / 1024^2,
  elapsed_seconds = unname(coordinate_elapsed[["elapsed"]]),
  stringsAsFactors = FALSE
), row.names = FALSE)
