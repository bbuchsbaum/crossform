# Capability-gated views over completed measurement forms ------------------

.coupling_result_signature <- function(fields) {
  paste0("sha256:", digest::digest(c(
    list(schema_version = 1L), fields
  ), algo = "sha256", serialize = TRUE))
}

.new_coupling_result <- function(kind, values, x, normalization_axis,
                                 summary_axis, stage_order,
                                 regularization = NULL, units = NULL,
                                 terminology = NULL,
                                 partition_policy = NULL,
                                 source_forms = list(x)) {
  .validate_measurement_form(x, probe = FALSE)
  if (!is.list(source_forms) || length(source_forms) < 1L) {
    stop("Coupling results require at least one completed source form.",
      call. = FALSE)
  }
  lapply(source_forms, .validate_measurement_form, probe = FALSE)
  if (!is.null(partition_policy)) {
    .validate_coupling_partition_policy(partition_policy)
  }
  if (!is.character(kind) || length(kind) != 1L || is.na(kind) ||
      !nzchar(kind) || !is.character(normalization_axis) ||
      length(normalization_axis) != 1L || is.na(normalization_axis) ||
      !normalization_axis %in% c("none", "experimental_samples",
        "neural_features", "partition_pairs", "form_entries") ||
      !is.character(summary_axis) || length(summary_axis) != 1L ||
      is.na(summary_axis) ||
      !summary_axis %in% c("measurement_coordinates", "form_entries",
        "canonical_modes") || !is.character(stage_order) ||
      length(stage_order) < 1L || anyNA(stage_order)) {
    stop("Coupling-view identity or axes are invalid.", call. = FALSE)
  }
  valid_values <- if (is.data.frame(values)) {
    nrow(values) >= 1L && "edge_id" %in% names(values) &&
      all(values$edge_id %in% x$block_index$edge_id) &&
      all(vapply(values, function(value) {
        if (is.numeric(value)) all(is.finite(value)) else
          all(!is.na(value))
      }, logical(1)))
  } else if (is.list(values)) {
    length(values) == nrow(x$block_index) &&
      identical(names(values), x$block_index$edge_id) &&
      all(vapply(values, function(value) {
        is.matrix(value) && is.numeric(value) && all(is.finite(value))
      }, logical(1)))
  } else {
    FALSE
  }
  if (!valid_values) {
    stop("Coupling-view values do not match their completed edge set.",
      call. = FALSE)
  }
  fields <- list(
    kind = kind,
    values = values,
    edge_index = x$block_index,
    source_plan = vapply(source_forms, function(form) {
      form$plan$signature
    }, character(1)),
    source_receipt = vapply(source_forms, function(form) {
      form$receipt$signature
    }, character(1)),
    normalization_axis = normalization_axis,
    summary_axis = summary_axis,
    stage_order = stage_order,
    regularization = regularization,
    units = units,
    terminology = terminology,
    partition_policy = partition_policy,
    edge_completeness = x$edge_completeness
  )
  structure(c(fields, list(
    signature = .coupling_result_signature(fields)
  )), class = "effect_coupling_result")
}

.validate_coupling_result <- function(x) {
  expected <- c("kind", "values", "edge_index", "source_plan",
    "source_receipt", "normalization_axis", "summary_axis", "stage_order",
    "regularization", "units", "terminology", "partition_policy",
    "edge_completeness", "signature")
  if (!inherits(x, "effect_coupling_result") || !is.list(x) ||
      !identical(names(x), expected)) {
    stop("Coupling result is missing or noncanonical.", call. = FALSE)
  }
  fields <- x[setdiff(names(x), "signature")]
  if (!identical(x$signature, .coupling_result_signature(fields))) {
    stop("Coupling-result identity is inconsistent.", call. = FALSE)
  }
  if (!is.character(x$source_plan) || length(x$source_plan) < 1L ||
      anyNA(x$source_plan) || !is.character(x$source_receipt) ||
      length(x$source_receipt) != length(x$source_plan) ||
      anyNA(x$source_receipt)) {
    stop("Coupling-result source identities are inconsistent.", call. = FALSE)
  }
  if (!is.null(x$partition_policy)) {
    .validate_coupling_partition_policy(x$partition_policy)
  }
  invisible(x)
}

.coupling_partition_policy <- function(
    weights, placement = c("within_partition_pair",
                           "after_partition_aggregation"),
    transform = NULL) {
  placement <- match.arg(placement)
  if (!is.numeric(weights) || length(weights) < 1L || anyNA(weights) ||
      any(!is.finite(weights)) || any(weights < 0) || sum(weights) <= 0) {
    stop("Partitioned coupling requires finite nonnegative positive-mass weights.",
      call. = FALSE)
  }
  weights <- as.numeric(weights / sum(weights))
  if (is.null(transform)) transform <- .identity_edge_transform()
  transform <- .validate_edge_transform(transform)
  if (!transform$kind %in% c("identity", "fisher_z")) {
    stop("Partitioned Pearson coupling supports identity or Fisher transforms.",
      call. = FALSE)
  }
  semantic <- list(
    schema_version = 1L,
    placement = placement,
    weights = weights,
    transform = unclass(transform)
  )
  structure(c(semantic[-1L], list(
    signature = paste0("sha256:", digest::digest(
      semantic, algo = "sha256", serialize = TRUE
    ))
  )), class = "effect_coupling_partition_policy")
}

.validate_coupling_partition_policy <- function(x) {
  expected <- c("placement", "weights", "transform", "signature")
  if (!inherits(x, "effect_coupling_partition_policy") || !is.list(x) ||
      !identical(names(x), expected)) {
    stop("Coupling partition policy is missing or noncanonical.",
      call. = FALSE)
  }
  transform <- structure(x$transform, class = "effect_edge_transform")
  rebuilt <- .coupling_partition_policy(x$weights, x$placement, transform)
  if (!identical(x, rebuilt)) {
    stop("Coupling partition-policy identity is inconsistent.",
      call. = FALSE)
  }
  rebuilt
}

.coupling_blocks <- function(x) {
  stats::setNames(lapply(seq_len(nrow(x$block_index)), function(edge) {
    .measurement_block(x, edge)
  }), x$block_index$edge_id)
}

.coupling_self_blocks <- function(x, edge, tolerance = 1e-10,
                                  ridge_left = 0, ridge_right = 0,
                                  require_positive = TRUE) {
  .validate_measurement_form(x, probe = FALSE)
  position <- .measurement_edge_position(x$block_index, edge)
  row <- x$block_index[position, , drop = FALSE]
  if (!identical(x$left_frame$signature, x$right_frame$signature)) {
    stop("Normalized coupling requires one compatible self-measurement frame.",
      call. = FALSE)
  }
  locate <- function(left, right) {
    which(x$block_index$left == left & x$block_index$right == right)
  }
  left_position <- locate(row$left[[1L]], row$left[[1L]])
  right_position <- locate(row$right[[1L]], row$right[[1L]])
  if (length(left_position) != 1L || length(right_position) != 1L) {
    stop("Normalized coupling requires explicit left and right self-blocks.",
      call. = FALSE)
  }
  cross <- .measurement_block(x, position)
  left_self <- .measurement_block(x, left_position)
  right_self <- .measurement_block(x, right_position)
  if (max(abs(left_self - t(left_self))) > tolerance ||
      max(abs(right_self - t(right_self))) > tolerance) {
    stop("Coupling self-blocks must be symmetric.", call. = FALSE)
  }
  left_eigen <- eigen((left_self + t(left_self)) / 2,
    symmetric = TRUE, only.values = TRUE)$values
  right_eigen <- eigen((right_self + t(right_self)) / 2,
    symmetric = TRUE, only.values = TRUE)$values
  if (require_positive &&
      (min(left_eigen) < -tolerance || min(right_eigen) < -tolerance)) {
    stop("Indefinite self-blocks cannot normalize coupling.", call. = FALSE)
  }
  list(
    edge_id = row$edge_id[[1L]],
    cross = cross,
    left_self = left_self,
    right_self = right_self,
    left_regularized = left_self + diag(ridge_left, nrow(left_self)),
    right_regularized = right_self + diag(ridge_right, nrow(right_self)),
    left_eigenvalues = left_eigen,
    right_eigenvalues = right_eigen
  )
}

.require_nondegenerate_variation <- function(x) {
  if (!isTRUE(x$capabilities$repeated_variation) ||
      is.null(x$capabilities$sampling_axis)) {
    stop("Normalized connectivity requires certified repeated variation.",
      call. = FALSE)
  }
  if (x$diagnostics$experimental_effective_rank <= 1L) {
    stop(paste0(
      "Normalized connectivity requires effective sampling rank above one; ",
      "a rank-one effect direction is degenerate."
    ), call. = FALSE)
  }
  invisible(x)
}

.effect_coupling <- function(x) {
  .require_measurement_view_capabilities(x, "effect_coupling")
  .new_coupling_result(
    "effect_coupling",
    .coupling_blocks(x),
    x,
    normalization_axis = "none",
    summary_axis = "measurement_coordinates",
    stage_order = c(x$plan$stages$order,
      "effect_coupling_materialization"),
    terminology = "algebraic_effect_coupling_no_covariance_claim"
  )
}

.covariance_coupling <- function(x, tolerance = 1e-10) {
  .require_nondegenerate_variation(x)
  .require_measurement_view_capabilities(
    x, "canonical_coupling", positive_self_blocks = TRUE
  )
  for (edge in seq_len(nrow(x$block_index))) {
    .coupling_self_blocks(x, edge, tolerance = tolerance)
  }
  .new_coupling_result(
    "covariance_coupling",
    .coupling_blocks(x),
    x,
    normalization_axis = "none",
    summary_axis = "measurement_coordinates",
    stage_order = c(x$plan$stages$order, "joint_covariance_validation",
      "covariance_materialization"),
    terminology = "repeated_sample_covariance_coupling"
  )
}

.pearson_coupling <- function(x, tolerance = 1e-10) {
  .require_nondegenerate_variation(x)
  .require_measurement_view_capabilities(
    x, "canonical_coupling", positive_self_blocks = TRUE
  )
  rows <- lapply(seq_len(nrow(x$block_index)), function(edge) {
    blocks <- .coupling_self_blocks(x, edge, tolerance = tolerance)
    if (!identical(dim(blocks$cross), c(1L, 1L))) {
      stop("Pearson coupling requires rank-one measurement axes.",
        call. = FALSE)
    }
    denominator <- sqrt(blocks$left_self[[1L]] * blocks$right_self[[1L]])
    if (!is.finite(denominator) || denominator <= tolerance) {
      stop("Pearson coupling requires strictly positive scalar self-variance.",
        call. = FALSE)
    }
    value <- blocks$cross[[1L]] / denominator
    if (abs(value) > 1 + 10 * tolerance) {
      stop("Scalar covariance blocks do not define a valid correlation.",
        call. = FALSE)
    }
    data.frame(edge_id = blocks$edge_id,
      correlation = max(-1, min(1, value)), stringsAsFactors = FALSE)
  })
  .new_coupling_result(
    "pearson_correlation",
    do.call(rbind, rows),
    x,
    normalization_axis = "experimental_samples",
    summary_axis = "measurement_coordinates",
    stage_order = c(x$plan$stages$order, "joint_covariance_validation",
      "experimental_sample_normalization", "signed_scalar_summary"),
    units = "correlation",
    terminology = "signed_functional_connectivity"
  )
}

.validate_partitioned_coupling_forms <- function(forms) {
  if (!is.list(forms) || length(forms) < 1L) {
    stop("Partitioned coupling requires a nonempty list of measurement forms.",
      call. = FALSE)
  }
  lapply(forms, .validate_measurement_form, probe = FALSE)
  reference <- forms[[1L]]
  compatible <- vapply(forms, function(form) {
    identical(form$block_index, reference$block_index) &&
      identical(form$left_frame$signature, reference$left_frame$signature) &&
      identical(form$right_frame$signature, reference$right_frame$signature) &&
      identical(form$plan$experimental_query,
        reference$plan$experimental_query) &&
      identical(form$capabilities$sampling_axis,
        reference$capabilities$sampling_axis)
  }, logical(1))
  if (!all(compatible)) {
    stop(paste0(
      "Partitioned coupling forms must share exact measurement edges, frames, ",
      "experimental query, and sampling axis."
    ), call. = FALSE)
  }
  lapply(forms, function(form) {
    .require_nondegenerate_variation(form)
    .require_measurement_view_capabilities(
      form, "canonical_coupling", positive_self_blocks = TRUE
    )
  })
  forms
}

.scalar_correlations_from_blocks <- function(form, blocks,
                                              tolerance = 1e-10) {
  if (!is.list(blocks) || length(blocks) != nrow(form$block_index)) {
    stop("Scalar coupling blocks are absent or incomplete.", call. = FALSE)
  }
  vapply(seq_len(nrow(form$block_index)), function(edge) {
    row <- form$block_index[edge, , drop = FALSE]
    left_self <- which(form$block_index$left == row$left[[1L]] &
      form$block_index$right == row$left[[1L]])
    right_self <- which(form$block_index$left == row$right[[1L]] &
      form$block_index$right == row$right[[1L]])
    if (length(left_self) != 1L || length(right_self) != 1L ||
        !identical(dim(blocks[[edge]]), c(1L, 1L)) ||
        !identical(dim(blocks[[left_self]]), c(1L, 1L)) ||
        !identical(dim(blocks[[right_self]]), c(1L, 1L))) {
      stop(paste0(
        "Partitioned Pearson coupling requires rank-one axes and explicit ",
        "self-blocks."
      ), call. = FALSE)
    }
    denominator <- sqrt(
      blocks[[left_self]][[1L]] * blocks[[right_self]][[1L]]
    )
    if (!is.finite(denominator) || denominator <= tolerance) {
      stop("Partitioned Pearson coupling requires positive self-variance.",
        call. = FALSE)
    }
    value <- blocks[[edge]][[1L]] / denominator
    if (!is.finite(value) || abs(value) > 1 + 10 * tolerance) {
      stop("Partitioned covariance blocks do not define valid correlations.",
        call. = FALSE)
    }
    max(-1, min(1, value))
  }, numeric(1))
}

.partitioned_pearson_coupling <- function(
    forms, weights = rep(1, length(forms)),
    placement = c("within_partition_pair", "after_partition_aggregation"),
    transform = NULL, tolerance = 1e-10) {
  forms <- .validate_partitioned_coupling_forms(forms)
  policy <- .coupling_partition_policy(
    weights, match.arg(placement), transform
  )
  if (length(policy$weights) != length(forms)) {
    stop("Partition coupling weights must match the source forms.",
      call. = FALSE)
  }
  transform <- structure(policy$transform, class = "effect_edge_transform")
  reference <- forms[[1L]]
  if (policy$placement == "within_partition_pair") {
    correlations <- vapply(forms, function(form) {
      result <- .pearson_coupling(form, tolerance)
      result$values$correlation[
        match(reference$block_index$edge_id, result$values$edge_id)
      ]
    }, numeric(nrow(reference$block_index)))
    transformed <- apply(correlations, 2L, function(value) {
      drop(.apply_edge_transform(matrix(value, ncol = 1L), transform))
    })
    if (is.null(dim(transformed))) {
      transformed <- matrix(transformed,
        nrow = nrow(reference$block_index), ncol = length(forms))
    }
    values <- drop(transformed %*% policy$weights)
    stage_order <- c(
      "partition_pair_covariance",
      "joint_covariance_validation",
      "experimental_sample_normalization",
      "edge_transform",
      "partition_reduction"
    )
  } else {
    blocks_by_form <- lapply(forms, .coupling_blocks)
    aggregate <- lapply(seq_len(nrow(reference$block_index)), function(edge) {
      Reduce(`+`, Map(function(blocks, weight) {
        blocks[[edge]] * weight
      }, blocks_by_form, policy$weights))
    })
    correlations <- .scalar_correlations_from_blocks(
      reference, aggregate, tolerance
    )
    values <- drop(.apply_edge_transform(
      matrix(correlations, ncol = 1L), transform
    ))
    stage_order <- c(
      "partition_pair_covariance",
      "partition_reduction",
      "joint_covariance_validation",
      "experimental_sample_normalization",
      "edge_transform"
    )
  }
  value_table <- data.frame(
    edge_id = reference$block_index$edge_id,
    value = values,
    transform = transform$kind,
    stringsAsFactors = FALSE
  )
  .new_coupling_result(
    "partitioned_pearson_coupling",
    value_table,
    reference,
    normalization_axis = "experimental_samples",
    summary_axis = "measurement_coordinates",
    stage_order = stage_order,
    units = if (transform$kind == "identity") "correlation" else "fisher_z",
    terminology = "ordered_signed_functional_connectivity",
    partition_policy = policy,
    source_forms = forms
  )
}

.canonical_coupling <- function(x, regularization, tolerance = 1e-10) {
  .require_nondegenerate_variation(x)
  .require_measurement_view_capabilities(
    x, "canonical_coupling", positive_self_blocks = TRUE
  )
  regularization <- .validate_measurement_regularization(regularization)
  if (regularization$kind != "ridge" || regularization$applied) {
    stop(paste0(
      "Canonical coupling requires an explicit unapplied ridge policy; the ",
      "view applies it to its self-blocks."
    ), call. = FALSE)
  }
  rows <- lapply(seq_len(nrow(x$block_index)), function(edge) {
    blocks <- .coupling_self_blocks(
      x, edge, tolerance,
      regularization$lambda_left,
      regularization$lambda_right
    )
    left_inverse <- .decomposition_inverse_sqrt(
      blocks$left_regularized, tolerance
    )
    right_inverse <- .decomposition_inverse_sqrt(
      blocks$right_regularized, tolerance
    )
    values <- svd(left_inverse %*% blocks$cross %*% right_inverse,
      nu = 0L, nv = 0L)$d
    if (any(values > 1 + 10 * tolerance)) {
      stop("Regularized covariance blocks yield canonical values above one.",
        call. = FALSE)
    }
    data.frame(
      edge_id = blocks$edge_id,
      mode = seq_along(values),
      canonical_correlation = pmin(1, pmax(0, values)),
      stringsAsFactors = FALSE
    )
  })
  .new_coupling_result(
    "canonical_coupling",
    do.call(rbind, rows),
    x,
    normalization_axis = "experimental_samples",
    summary_axis = "canonical_modes",
    stage_order = c(x$plan$stages$order, "joint_covariance_validation",
      "ridge_regularization", "experimental_sample_whitening", "svd"),
    regularization = regularization,
    units = "canonical_correlation",
    terminology = "multivariate_covariance_coupling_spectrum"
  )
}

.geometry_alignment <- function(x, tolerance = 1e-10) {
  .require_nondegenerate_variation(x)
  if (!isTRUE(x$capabilities$joint_covariance)) {
    stop("Geometry alignment requires a coherent positive joint covariance.",
      call. = FALSE)
  }
  rows <- lapply(seq_len(nrow(x$block_index)), function(edge) {
    blocks <- .coupling_self_blocks(x, edge, tolerance)
    denominator <- sqrt(sum(blocks$left_self^2) * sum(blocks$right_self^2))
    if (!is.finite(denominator) || denominator <= tolerance) {
      stop("Geometry alignment requires nonzero self-geometries.",
        call. = FALSE)
    }
    value <- sum(blocks$cross^2) / denominator
    if (value > 1 + 10 * tolerance) {
      stop("Covariance blocks do not define bounded geometry alignment.",
        call. = FALSE)
    }
    data.frame(edge_id = blocks$edge_id,
      geometry_alignment = min(1, max(0, value)), stringsAsFactors = FALSE)
  })
  .new_coupling_result(
    "geometry_alignment",
    do.call(rbind, rows),
    x,
    normalization_axis = "form_entries",
    summary_axis = "form_entries",
    stage_order = c(x$plan$stages$order, "joint_covariance_validation",
      "form_entry_frobenius_normalization", "geometry_alignment"),
    units = "linear_cka",
    terminology = paste0(
      "static_geometry_alignment_not_dynamic_informational_connectivity"
    )
  )
}

.gaussian_covariance_model <- function(provenance = list()) {
  .validate_effect_provenance(provenance, "Gaussian-model provenance")
  semantic <- list(
    schema_version = 1L,
    family = "joint_gaussian_covariance",
    fixed = TRUE,
    provenance = provenance
  )
  structure(c(semantic[-1L], list(
    signature = paste0("sha256:", digest::digest(
      semantic, algo = "sha256", serialize = TRUE
    ))
  )), class = "effect_gaussian_covariance_model")
}

.validate_gaussian_covariance_model <- function(x) {
  expected <- c("family", "fixed", "provenance", "signature")
  if (!inherits(x, "effect_gaussian_covariance_model") || !is.list(x) ||
      !identical(names(x), expected) ||
      !identical(x$family, "joint_gaussian_covariance") ||
      !identical(x$fixed, TRUE)) {
    stop("Gaussian covariance-model declaration is missing or invalid.",
      call. = FALSE)
  }
  rebuilt <- .gaussian_covariance_model(x$provenance)
  if (!identical(x, rebuilt)) {
    stop("Gaussian covariance-model identity is inconsistent.",
      call. = FALSE)
  }
  rebuilt
}

.gaussian_information <- function(x, regularization, model,
                                  units = c("nats", "bits"),
                                  tolerance = 1e-10) {
  units <- match.arg(units)
  model <- .validate_gaussian_covariance_model(model)
  canonical <- .canonical_coupling(x, regularization, tolerance)
  split_values <- split(
    canonical$values$canonical_correlation,
    canonical$values$edge_id
  )
  rows <- lapply(names(split_values), function(edge_id) {
    rho <- split_values[[edge_id]]
    if (any(rho >= 1 - tolerance)) {
      stop("Gaussian information is infinite at canonical correlation one.",
        call. = FALSE)
    }
    value <- -0.5 * sum(log1p(-(rho^2)))
    if (units == "bits") value <- value / log(2)
    data.frame(edge_id = edge_id, information = value,
      units = units, stringsAsFactors = FALSE)
  })
  .new_coupling_result(
    "gaussian_mutual_information",
    do.call(rbind, rows),
    x,
    normalization_axis = "experimental_samples",
    summary_axis = "canonical_modes",
    stage_order = c(canonical$stage_order,
      "joint_gaussian_information_transform"),
    regularization = regularization,
    units = units,
    terminology = paste0("Gaussian_mutual_information_model:", model$signature)
  )
}

.connectivity_view <- function(
    x, view = c("correlation", "canonical", "geometry_alignment",
                "gaussian_information"),
    regularization = NULL, model = NULL, units = c("nats", "bits"),
    tolerance = 1e-10) {
  view <- match.arg(view)
  switch(view,
    correlation = .pearson_coupling(x, tolerance),
    canonical = {
      if (is.null(regularization)) {
        stop("Canonical connectivity requires an explicit regularization policy.",
          call. = FALSE)
      }
      .canonical_coupling(x, regularization, tolerance)
    },
    geometry_alignment = .geometry_alignment(x, tolerance),
    gaussian_information = {
      if (is.null(regularization) || is.null(model)) {
        stop(paste0(
          "Gaussian-information connectivity requires explicit ",
          "regularization and a Gaussian model declaration."
        ), call. = FALSE)
      }
      .gaussian_information(
        x, regularization, model, match.arg(units), tolerance
      )
    }
  )
}
