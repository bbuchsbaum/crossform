# Independent numerical oracles for evidence-pairing-v1. These helpers do not
# call production crossform functions.

pair_oracle_matrix <- function(value, name) {
  if (!is.matrix(value) || !is.numeric(value) || any(!is.finite(value))) {
    stop(name, " must be one finite numeric matrix")
  }
  value
}

pair_oracle_space <- function(id, coordinates,
                              kind = c("experimental", "neural", "measurement")) {
  kind <- match.arg(kind)
  if (!is.character(id) || length(id) != 1L || is.na(id) || !nzchar(id) ||
      !is.character(coordinates) || anyNA(coordinates) ||
      anyDuplicated(coordinates)) {
    stop("space requires one id and unique ordered coordinates")
  }
  list(id = id, coordinates = coordinates, kind = kind)
}

pair_oracle_relation <- function(value, experimental_space, neural_space) {
  value <- pair_oracle_matrix(value, "relation")
  if (!identical(experimental_space$kind, "experimental") ||
      !identical(neural_space$kind, "neural") ||
      nrow(value) != length(experimental_space$coordinates) ||
      ncol(value) != length(neural_space$coordinates)) {
    stop("relation dimensions do not match its identified boundaries")
  }
  list(
    value = value,
    experimental_space = experimental_space,
    neural_space = neural_space
  )
}

pair_oracle_bound_query <- function(value, left_space, right_space,
                                    role = c("effect", "variation")) {
  value <- pair_oracle_matrix(value, "pair query")
  role <- match.arg(role)
  if (nrow(value) != length(left_space$coordinates) ||
      ncol(value) != length(right_space$coordinates)) {
    stop("query dimensions do not match its identified boundaries")
  }
  list(
    value = value,
    left_space = left_space,
    right_space = right_space,
    role = role
  )
}

pair_oracle_assert_boundaries <- function(left_relation, right_relation,
                                          experimental_query, neural_query) {
  if (!identical(experimental_query$left_space,
        left_relation$experimental_space) ||
      !identical(experimental_query$right_space,
        right_relation$experimental_space)) {
    stop("experimental query axis identities differ from relation boundaries")
  }
  if (!identical(neural_query$left_space, left_relation$neural_space) ||
      !identical(neural_query$right_space, right_relation$neural_space)) {
    stop("neural query axis identities differ from relation boundaries")
  }
  invisible(TRUE)
}

pair_oracle_forward <- function(left, right, neural_query) {
  left <- pair_oracle_matrix(left, "left relation")
  right <- pair_oracle_matrix(right, "right relation")
  neural_query <- pair_oracle_matrix(neural_query, "neural query")
  if (ncol(left) != nrow(neural_query) ||
      ncol(right) != ncol(neural_query)) {
    stop("forward transport dimensions do not conform")
  }
  left %*% neural_query %*% t(right)
}

pair_oracle_adjoint <- function(left, right, experimental_query) {
  left <- pair_oracle_matrix(left, "left relation")
  right <- pair_oracle_matrix(right, "right relation")
  experimental_query <- pair_oracle_matrix(
    experimental_query, "experimental query"
  )
  if (nrow(left) != nrow(experimental_query) ||
      nrow(right) != ncol(experimental_query)) {
    stop("adjoint transport dimensions do not conform")
  }
  t(left) %*% experimental_query %*% right
}

pair_oracle_evidence <- function(left, right, experimental_query,
                                 neural_query) {
  sum(experimental_query * pair_oracle_forward(left, right, neural_query))
}

pair_oracle_evidence_adjoint <- function(left, right, experimental_query,
                                         neural_query) {
  sum(pair_oracle_adjoint(left, right, experimental_query) * neural_query)
}

pair_oracle_typed_evidence <- function(left_relation, right_relation,
                                       experimental_query, neural_query) {
  pair_oracle_assert_boundaries(
    left_relation, right_relation, experimental_query, neural_query
  )
  pair_oracle_evidence(
    left_relation$value,
    right_relation$value,
    experimental_query$value,
    neural_query$value
  )
}

pair_oracle_second_order_lift <- function(left, right, materialize = FALSE,
                                          max_entries = 10000L) {
  left <- pair_oracle_matrix(left, "left relation")
  right <- pair_oracle_matrix(right, "right relation")
  if (!isTRUE(materialize)) {
    stop("the second-order Kronecker lift is conceptual by default")
  }
  entries <- as.double(nrow(left)) * ncol(left) * nrow(right) * ncol(right)
  if (!is.numeric(max_entries) || length(max_entries) != 1L ||
      !is.finite(max_entries) || max_entries < 1 || entries > max_entries) {
    stop("explicit diagnostic lift exceeds its bounded size guard")
  }
  kronecker(right, left)
}

pair_oracle_measurement_form <- function(left, right, experimental_query,
                                         left_leg, right_leg) {
  left_leg <- pair_oracle_matrix(left_leg, "left measurement leg")
  right_leg <- pair_oracle_matrix(right_leg, "right measurement leg")
  neural_form <- pair_oracle_adjoint(left, right, experimental_query)
  if (ncol(left_leg) != nrow(neural_form) ||
      ncol(right_leg) != ncol(neural_form)) {
    stop("measurement legs do not match native neural boundaries")
  }
  left_leg %*% neural_form %*% t(right_leg)
}

pair_oracle_measured_effect_form <- function(left, right, left_leg, right_leg,
                                             measurement_query) {
  left_leg <- pair_oracle_matrix(left_leg, "left measurement leg")
  right_leg <- pair_oracle_matrix(right_leg, "right measurement leg")
  measurement_query <- pair_oracle_matrix(
    measurement_query, "measurement query"
  )
  if (nrow(measurement_query) != nrow(left_leg) ||
      ncol(measurement_query) != nrow(right_leg)) {
    stop("measurement query axes do not match measurement legs")
  }
  neural_query <- t(left_leg) %*% measurement_query %*% right_leg
  pair_oracle_forward(left, right, neural_query)
}

pair_oracle_factorized_evidence <- function(left, right, u, v, r, s) {
  left <- pair_oracle_matrix(left, "left relation")
  right <- pair_oracle_matrix(right, "right relation")
  u <- pair_oracle_matrix(u, "left experimental factor")
  v <- pair_oracle_matrix(v, "right experimental factor")
  r <- pair_oracle_matrix(r, "left neural factor")
  s <- pair_oracle_matrix(s, "right neural factor")
  if (nrow(u) != nrow(left) || nrow(v) != nrow(right) ||
      nrow(r) != ncol(left) || nrow(s) != ncol(right) ||
      ncol(u) != ncol(v) || ncol(r) != ncol(s)) {
    stop("factor dimensions do not match relation boundaries")
  }
  sum((t(u) %*% left %*% r) * (t(v) %*% right %*% s))
}

pair_oracle_query_capabilities <- function(
    role = c("effect", "variation"),
    construction = c("arbitrary", "outer_product", "centering", "covariance"),
    sampling_axis = NULL,
    effective_rank = NA_integer_,
    observed_min_eigenvalue = NA_real_) {
  role <- match.arg(role)
  construction <- match.arg(construction)
  if (!is.numeric(effective_rank) || length(effective_rank) != 1L ||
      (!is.na(effective_rank) && effective_rank < 0)) {
    stop("effective_rank must be one nonnegative diagnostic")
  }
  if (!is.numeric(observed_min_eigenvalue) ||
      length(observed_min_eigenvalue) != 1L) {
    stop("observed eigenvalue must be one diagnostic")
  }
  if (role == "variation" &&
      (!is.character(sampling_axis) || length(sampling_axis) != 1L ||
        is.na(sampling_axis) || !nzchar(sampling_axis))) {
    stop("variation role requires declared sampling provenance")
  }
  if (construction %in% c("centering", "covariance") && role != "variation") {
    stop("sampling constructions require the variation role")
  }
  if (construction == "outer_product" && role != "effect") {
    stop("an outer-product contrast is an effect construction")
  }
  list(
    guarantees = list(
      role = role,
      positive_query = construction %in%
        c("outer_product", "centering", "covariance"),
      joint_covariance = identical(construction, "covariance"),
      sampling_axis = sampling_axis
    ),
    diagnostics = list(
      effective_rank = as.integer(effective_rank),
      observed_min_eigenvalue = observed_min_eigenvalue
    ),
    provenance = list(construction = construction)
  )
}

pair_oracle_require_view <- function(capabilities,
                                     view = c("effect_coupling",
                                              "covariance_coupling",
                                              "canonical_coupling"),
                                     positive_self_blocks = FALSE,
                                     regularization_policy = NULL) {
  view <- match.arg(view)
  if (view == "effect_coupling") return(invisible(TRUE))
  guarantees <- capabilities$guarantees
  diagnostics <- capabilities$diagnostics
  if (!identical(guarantees$role, "variation") ||
      !isTRUE(guarantees$positive_query) ||
      is.null(guarantees$sampling_axis) ||
      is.na(diagnostics$effective_rank) || diagnostics$effective_rank <= 1L) {
    stop("normalized coupling requires valid repeated variation of rank above one")
  }
  if (view == "canonical_coupling" &&
      (!isTRUE(guarantees$joint_covariance) ||
        !isTRUE(positive_self_blocks) ||
        !is.character(regularization_policy) ||
        length(regularization_policy) != 1L ||
        is.na(regularization_policy) || !nzchar(regularization_policy))) {
    stop("canonical coupling requires joint covariance and regularized self-blocks")
  }
  invisible(TRUE)
}

pair_oracle_inverse_sqrt <- function(value, tolerance = 1e-12) {
  value <- pair_oracle_matrix(value, "self block")
  if (nrow(value) != ncol(value) ||
      max(abs(value - t(value))) > tolerance) {
    stop("self block must be symmetric")
  }
  decomposition <- eigen((value + t(value)) / 2, symmetric = TRUE)
  if (min(decomposition$values) < -tolerance) {
    stop("self block is not positive semidefinite")
  }
  scale <- numeric(length(decomposition$values))
  retained <- decomposition$values > tolerance
  scale[retained] <- 1 / sqrt(decomposition$values[retained])
  decomposition$vectors %*% diag(scale, nrow = length(scale)) %*%
    t(decomposition$vectors)
}

pair_oracle_canonical_spectrum <- function(cross, left_self, right_self,
                                           tolerance = 1e-12) {
  normalized <- pair_oracle_inverse_sqrt(left_self, tolerance) %*% cross %*%
    pair_oracle_inverse_sqrt(right_self, tolerance)
  svd(normalized, nu = 0L, nv = 0L)$d
}

pair_oracle_geometry_alignment <- function(cross, left_self, right_self) {
  denominator <- sqrt(sum(left_self^2) * sum(right_self^2))
  if (!is.finite(denominator) || denominator == 0) {
    stop("geometry alignment requires nonzero finite self-blocks")
  }
  sum(cross^2) / denominator
}

pair_oracle_stage <- function(operation,
                              axis = c("neural_features",
                                       "experimental_samples",
                                       "partition_pairs", "form_entries"),
                              placement = c("within_pair", "after_aggregation")) {
  axis <- match.arg(axis)
  placement <- match.arg(placement)
  if (!is.character(operation) || length(operation) != 1L ||
      is.na(operation) || !nzchar(operation)) {
    stop("operation requires one explicit name")
  }
  list(operation = operation, axis = axis, placement = placement)
}

pair_oracle_stage_identity <- function(stages) {
  if (!is.list(stages) || !length(stages) ||
      any(!vapply(stages, is.list, logical(1))) ||
      any(!vapply(stages, function(stage) {
        identical(sort(names(stage)), sort(c("operation", "axis", "placement")))
      }, logical(1)))) {
    stop("every ordered stage must identify operation, axis, and placement")
  }
  vapply(stages, function(stage) {
    paste(stage$operation, stage$axis, stage$placement, sep = "@")
  }, character(1))
}

pair_oracle_scalar_normalization <- function(cross, left_self, right_self,
                                             weight,
                                             policy = c("within_then_reduce",
                                                        "reduce_then_normalize")) {
  policy <- match.arg(policy)
  if (length(cross) != length(left_self) ||
      length(cross) != length(right_self) ||
      length(cross) != length(weight) ||
      any(left_self <= 0) || any(right_self <= 0)) {
    stop("scalar normalization requires matching positive self moments")
  }
  if (policy == "within_then_reduce") {
    return(sum(weight * cross / sqrt(left_self * right_self)))
  }
  sum(weight * cross) /
    sqrt(sum(weight * left_self) * sum(weight * right_self))
}
