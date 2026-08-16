#!/usr/bin/env Rscript

# Exact, query-first scale qualification for evidence-sampling-v1.
# The benchmark uses a realistic 60k-feature volume support index and the
# local covariance shape induced by 120 experimental effects (7,140 RDM
# coordinates) on one 100-feature-scale searchlight. No dense D-by-D
# covariance is formed.

suppressPackageStartupMessages(devtools::load_all(quiet = TRUE))
repo <- normalizePath(getwd(), mustWork = TRUE)
source(file.path(repo, "benchmarks", "provenance.R"), local = TRUE)
provenance <- crossform_benchmark_provenance(
  repo, "run-sampling-covariance-scale.R"
)

measure <- function(expression) {
  garbage <- gc()
  before <- sum(garbage[, 2L]) * 1024^2
  elapsed <- system.time(value <- force(expression))[["elapsed"]]
  after <- sum(gc()[, 2L]) * 1024^2
  list(
    value = value,
    elapsed_seconds = unname(elapsed),
    resident_delta_bytes = max(0, after - before)
  )
}

set.seed(81360)
volume <- measure({
  domain <- volume_domain(
    array(TRUE, c(50L, 40L, 30L)), spacing = c(2, 2, 2),
    id = "sampling-scale-60k"
  )
  frame <- compile_frame(searchlights(6), domain)
  list(domain = domain, frame = frame)
})
domain <- volume$value$domain
frame <- volume$value$frame

conditions <- 120L
features <- as.integer(round(frame$support_index$cost$support_size[["mean"]]))
partitions <- 8L
contrasts <- crossform:::.sampling_distance_contrasts(
  paste0("condition", seq_len(conditions))
)$matrix
dimension <- nrow(contrasts)
signal_patterns <- matrix(rnorm(conditions * features), conditions, features)
effect_covariance <- diag(0.05, conditions)
residual_covariance <- diag(features)

# The scale fixture binds the covariance artifact to a real q = 120 plan:
# `.sampling_covariance_from_components()` checks the contrasts' declared
# experimental dimension and effect names against the bound plan's effect
# space, so the identity this artifact records is the identity it measured.
effect_labels <- paste0("condition", seq_len(conditions))
observations <- 2L * conditions
design <- kronecker(rep(1, 2L), diag(conditions))
colnames(design) <- effect_labels
effect_map <- diag(conditions)
dimnames(effect_map) <- list(effect_labels, effect_labels)
fixture_features <- 4L
fixture_domain <- abstract_domain(
  fixture_features, id = "sampling-scale-plan"
)
sources <- stats::setNames(lapply(seq_len(partitions), function(index) {
  matrix(
    rnorm(observations * fixture_features), observations, fixture_features
  )
}), paste0("run", seq_len(partitions)))
fit <- lm_relation_fit(
  sources, design, effect_map, sampling_unit = "trial",
  domain = fixture_domain
)
evidence <- plan_geometry(
  fit$relation, compile_frame(whole_brain(), fixture_domain),
  cross_partitions(fit$relation, independence = "independent")
)
plan <- crossform:::.compile_evidence_sampling_plan(
  evidence, fit,
  target = crossform:::.sampling_target(
    "point_estimate", policy = "partition_mean_plugin"
  )
)

construction <- measure(crossform:::.sampling_covariance_from_components(
  plan, contrasts, signal_patterns, effect_covariance,
  residual_covariance, labels = rownames(contrasts)
))
covariance <- construction$value

selected_index <- unique(as.integer(round(
  seq.int(1L, dimension, length.out = 100L)
)))
selected_pairs <- cbind(
  selected_index,
  rev(selected_index)
)
vector <- rnorm(dimension)
transport_map <- matrix(rnorm(8L * dimension), 8L, dimension)

diagonal <- measure(crossform:::.sampling_covariance_diagonal(covariance))
selected <- measure(crossform:::.sampling_covariance_entries(
  covariance, selected_pairs[, 1L], selected_pairs[, 2L]
))
quadratic <- measure(crossform:::.sampling_covariance_quadratic(
  covariance, vector
))
transport <- measure(crossform:::.sampling_covariance_transport(
  covariance, transport_map
))

dense_payload_bytes <- 8 * as.double(dimension)^2
dense_refusal <- tryCatch({
  crossform:::.sampling_covariance_materialize(
    covariance, max_bytes = 64 * 1024^2
  )
  FALSE
}, error = function(error) {
  grepl("materialization budget", conditionMessage(error), fixed = TRUE)
})

operations <- data.frame(
  operation = c(
    "volume_frame_compile", "factor_construction", "diagonal",
    "selected_100", "quadratic_form", "transport_8"
  ),
  elapsed_seconds = c(
    volume$elapsed_seconds, construction$elapsed_seconds,
    diagonal$elapsed_seconds, selected$elapsed_seconds,
    quadratic$elapsed_seconds, transport$elapsed_seconds
  ),
  resident_delta_bytes = c(
    volume$resident_delta_bytes, construction$resident_delta_bytes,
    diagonal$resident_delta_bytes, selected$resident_delta_bytes,
    quadratic$resident_delta_bytes, transport$resident_delta_bytes
  ),
  result_bytes = c(
    as.numeric(object.size(frame)), as.numeric(object.size(covariance)),
    as.numeric(object.size(diagonal$value)),
    as.numeric(object.size(selected$value)),
    as.numeric(object.size(quadratic$value)),
    as.numeric(object.size(transport$value))
  ),
  stringsAsFactors = FALSE
)

artifact <- list(
  schema_version = 1L,
  provenance = provenance,
  numerical_contract = "exact_factorized_hadamard_gram",
  generated_at = format(Sys.time(), tz = "UTC", usetz = TRUE),
  dimensions = list(
    domain_features = domain$n_features,
    frame_nodes = nrow(frame$weights),
    membership_nnz = as.double(Matrix::nnzero(frame$weights)),
    pair_pattern_nnz = frame$support_index$cost$pair_pattern_nnz,
    mean_support = frame$support_index$cost$support_size[["mean"]],
    max_support = frame$support_index$cost$support_size[["max"]],
    experimental_conditions = conditions,
    evidence_dimension = dimension,
    local_features = features,
    partitions = partitions,
    dense_covariance_payload_bytes = dense_payload_bytes
  ),
  factor_shapes = list(
    signal = dim(covariance$signal_factor),
    xi = dim(covariance$xi_factor)
  ),
  dense_materialization_refused_at_64_mib = dense_refusal,
  operations = operations,
  checks = list(
    diagonal_finite_nonnegative = all(is.finite(diagonal$value)) &&
      all(diagonal$value >= 0),
    selected_finite = all(is.finite(selected$value)),
    quadratic_finite_nonnegative = is.finite(quadratic$value) &&
      quadratic$value >= 0,
    transport_symmetric = max(abs(transport$value - t(transport$value))) <
      1e-10,
    no_dense_covariance_field = !any(c("delta", "xi", "covariance") %in%
      names(covariance))
  )
)

dir.create("benchmark-results", showWarnings = FALSE, recursive = TRUE)
saveRDS(artifact,
  file.path("benchmark-results", "sampling-covariance-scale.rds")
)
utils::write.csv(operations,
  file.path("benchmark-results", "sampling-covariance-scale-summary.csv"),
  row.names = FALSE
)

print(operations, row.names = FALSE)
print(artifact$dimensions)
print(artifact$checks)
cat("dense materialization refused at 64 MiB:", dense_refusal, "\n")
