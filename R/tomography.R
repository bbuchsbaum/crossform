# Dual-frame reconstruction of global neural evidence operators ------------

.tomography_signature <- function(kind, fields) {
  .sha256_signature(c(
    list(schema_version = 1L, kind = kind), fields
  ))
}

.tomography_stack_frame <- function(frame) {
  frame <- .validate_measurement_frame(frame)
  widths <- vapply(frame$legs, function(leg) nrow(leg$operator), integer(1))
  ends <- cumsum(widths)
  starts <- c(1L, utils::head(ends, -1L) + 1L)
  ranges <- Map(seq.int, starts, ends)
  names(ranges) <- frame$node_ids
  list(
    operator = do.call(rbind, lapply(frame$legs, `[[`, "operator")),
    ranges = ranges,
    dimensions = c(
      measurement = sum(widths),
      neural = frame$source_domain$n_features
    ),
    signature = frame$signature
  )
}

.tomography_resource_plan <- function(form, left_frame, right_frame,
                                      workspace_bytes = NULL) {
  .validate_measurement_form(form, probe = FALSE)
  left_frame <- .validate_measurement_frame(left_frame)
  right_frame <- .validate_measurement_frame(right_frame)
  if (!identical(left_frame$signature, form$left_frame$signature) ||
      !identical(right_frame$signature, form$right_frame$signature)) {
    stop("Tomography frames do not match the measurement-form bases.",
      call. = FALSE)
  }
  if (!is.null(workspace_bytes) &&
      (!is.numeric(workspace_bytes) || length(workspace_bytes) != 1L ||
       is.na(workspace_bytes) || !is.finite(workspace_bytes) ||
       workspace_bytes <= 0)) {
    stop("Tomography workspace must be NULL or one positive finite byte budget.",
      call. = FALSE)
  }
  left_rows <- sum(vapply(left_frame$legs, function(leg) {
    nrow(leg$operator)
  }, integer(1)))
  right_rows <- sum(vapply(right_frame$legs, function(leg) {
    nrow(leg$operator)
  }, integer(1)))
  p_left <- left_frame$source_domain$n_features
  p_right <- right_frame$source_domain$n_features
  dimensions <- as.double(c(left_rows, right_rows, p_left, p_right))
  names(dimensions) <- c(
    "left_measurement", "right_measurement", "left_neural", "right_neural"
  )
  products <- c(
    measured_form = dimensions[["left_measurement"]] *
      dimensions[["right_measurement"]],
    neural_operator = dimensions[["left_neural"]] *
      dimensions[["right_neural"]],
    frame_operators = dimensions[["left_measurement"]] *
      dimensions[["left_neural"]] +
      dimensions[["right_measurement"]] * dimensions[["right_neural"]],
    dual_operators = dimensions[["left_neural"]] *
      dimensions[["left_measurement"]] +
      dimensions[["right_neural"]] * dimensions[["right_measurement"]],
    projections = dimensions[["left_neural"]]^2 +
      dimensions[["right_neural"]]^2
  )
  if (any(!is.finite(products)) ||
      any(dimensions > .Machine$integer.max) ||
      products[["measured_form"]] > .Machine$integer.max ||
      products[["neural_operator"]] > .Machine$integer.max) {
    stop("Tomography dimensions exceed safe dense-matrix limits.",
      call. = FALSE)
  }
  component_bytes <- products * 8
  planned_workspace_bytes <- sum(component_bytes)
  fields <- list(
    dimensions = dimensions,
    component_bytes = component_bytes,
    planned_workspace_bytes = planned_workspace_bytes,
    budget_bytes = if (is.null(workspace_bytes)) NULL else
      as.numeric(workspace_bytes),
    fits_budget = if (is.null(workspace_bytes)) NA else
      planned_workspace_bytes <= workspace_bytes,
    dense_global_operator = TRUE,
    source_plan = form$plan$signature
  )
  structure(c(fields, list(
    signature = .tomography_signature("resource_plan", fields)
  )), class = "effect_tomography_resource_plan")
}

.validate_tomography_resource_plan <- function(x) {
  expected <- c("dimensions", "component_bytes", "planned_workspace_bytes",
    "budget_bytes", "fits_budget", "dense_global_operator", "source_plan",
    "signature")
  if (!inherits(x, "effect_tomography_resource_plan") || !is.list(x) ||
      !identical(names(x), expected) || !is.numeric(x$dimensions) ||
      !is.numeric(x$component_bytes) ||
      !is.numeric(x$planned_workspace_bytes) ||
      length(x$planned_workspace_bytes) != 1L ||
      !identical(x$dense_global_operator, TRUE)) {
    stop("Tomography resource plan is missing or noncanonical.",
      call. = FALSE)
  }
  fields <- x[setdiff(names(x), "signature")]
  if (!identical(x$signature,
      .tomography_signature("resource_plan", fields))) {
    stop("Tomography resource-plan identity is inconsistent.",
      call. = FALSE)
  }
  invisible(x)
}

.tomography_require_budget <- function(plan) {
  .validate_tomography_resource_plan(plan)
  if (identical(plan$fits_budget, FALSE)) {
    stop(sprintf(
      "Dense tomography requires %.0f bytes, exceeding the %.0f-byte budget.",
      plan$planned_workspace_bytes, plan$budget_bytes
    ), call. = FALSE)
  }
  invisible(plan)
}

.assemble_measurement_blocks <- function(form, left_frame, right_frame,
                                         workspace_bytes = NULL) {
  .validate_measurement_form(form, probe = FALSE)
  if (!isTRUE(form$capabilities$complete_edge_set) ||
      !identical(form$edge_completeness, "frame_complete")) {
    .capability_refusal(sprintf(paste0(
      "Reconstructing the global neural operator requires a block for every ",
      "ordered node pair, and this form carries %s (edge completeness: %s). ",
      "A diagonal-only or requested-edge map has discarded the between-node ",
      "directions, so no lossless reconstruction exists to return."
    ), .msg_count(nrow(form$block_index), "edge block"),
      form$edge_completeness),
      capability = "complete_edge_set",
      namespace = "tomography",
      reasons = "edge_set_is_not_frame_complete",
      remedies = paste0(
        "Build the form over every directed node pair, for example with ",
        "`edge_frame()` on `expand.grid(from = nodes$node_ids, ",
        "to = nodes$node_ids)`."
      )
    )
  }
  resource <- .tomography_resource_plan(
    form, left_frame, right_frame, workspace_bytes
  )
  .tomography_require_budget(resource)
  left <- .tomography_stack_frame(left_frame)
  right <- .tomography_stack_frame(right_frame)
  value <- matrix(0, left$dimensions[["measurement"]],
    right$dimensions[["measurement"]])
  for (left_node in names(left$ranges)) {
    for (right_node in names(right$ranges)) {
      position <- which(form$block_index$left == left_node &
        form$block_index$right == right_node)
      if (length(position) != 1L) {
        stop("Complete tomography requires one block for every ordered node pair.",
          call. = FALSE)
      }
      block <- .measurement_block(form, position)
      expected <- c(length(left$ranges[[left_node]]),
        length(right$ranges[[right_node]]))
      if (!identical(dim(block), expected)) {
        stop("Tomography block dimensions disagree with the frame axes.",
          call. = FALSE)
      }
      value[left$ranges[[left_node]], right$ranges[[right_node]]] <- block
    }
  }
  fields <- list(
    value = value,
    left_ranges = left$ranges,
    right_ranges = right$ranges,
    left_frame = left$signature,
    right_frame = right$signature,
    source_plan = form$plan$signature,
    source_receipt = form$receipt$signature,
    resource_plan = resource,
    completeness = "frame_complete"
  )
  structure(c(fields, list(
    signature = .tomography_signature("measured_block_form", fields)
  )), class = "effect_measured_block_form")
}

.validate_measured_block_form <- function(x) {
  expected <- c("value", "left_ranges", "right_ranges", "left_frame",
    "right_frame", "source_plan", "source_receipt", "resource_plan",
    "completeness", "signature")
  if (!inherits(x, "effect_measured_block_form") || !is.list(x) ||
      !identical(names(x), expected) || !is.matrix(x$value) ||
      !is.numeric(x$value) || any(!is.finite(x$value)) ||
      !identical(x$completeness, "frame_complete")) {
    stop("Measured block form is missing or noncanonical.", call. = FALSE)
  }
  .validate_tomography_resource_plan(x$resource_plan)
  fields <- x[setdiff(names(x), "signature")]
  if (!identical(x$signature,
      .tomography_signature("measured_block_form", fields))) {
    stop("Measured block-form identity is inconsistent.", call. = FALSE)
  }
  invisible(x)
}

.tomography_frame_diagnostics <- function(operator, tolerance) {
  decomposition <- svd(operator, nu = 0L, nv = 0L)
  threshold <- if (!length(decomposition$d)) 0 else
    max(decomposition$d) * tolerance
  retained <- decomposition$d > threshold
  rank <- as.integer(sum(retained))
  p <- ncol(operator)
  retained_condition <- if (!rank) Inf else
    max(decomposition$d[retained]) / min(decomposition$d[retained])
  list(
    measurement_dimension = nrow(operator),
    neural_dimension = p,
    rank = rank,
    full_column_rank = rank == p,
    singular_values = decomposition$d,
    relative_cutoff = tolerance,
    absolute_cutoff = threshold,
    condition_number = if (rank == p) retained_condition else Inf,
    retained_condition_number = retained_condition,
    parseval_residual = max(abs(crossprod(operator) - diag(p)))
  )
}

.tomography_pseudoinverse <- function(operator, tolerance) {
  decomposition <- svd(operator,
    nu = min(dim(operator)), nv = min(dim(operator)))
  threshold <- if (!length(decomposition$d)) 0 else
    max(decomposition$d) * tolerance
  retained <- which(decomposition$d > threshold)
  if (!length(retained)) {
    return(matrix(0, ncol(operator), nrow(operator)))
  }
  decomposition$v[, retained, drop = FALSE] %*%
    diag(1 / decomposition$d[retained], length(retained)) %*%
    t(decomposition$u[, retained, drop = FALSE])
}

.tomography_reject <- function(message, diagnostics) {
  condition <- structure(list(
    message = message,
    call = NULL,
    diagnostics = diagnostics
  ), class = c("effect_tomography_rejection", "error", "condition"))
  stop(condition)
}

.new_tomography_result <- function(operator, method, status, lossless,
                                   certified, left_projection,
                                   right_projection, diagnostics,
                                   measured, reference_signature = NULL) {
  fields <- list(
    operator = unname(operator),
    method = method,
    status = status,
    lossless = lossless,
    certified = certified,
    left_projection = unname(left_projection),
    right_projection = unname(right_projection),
    diagnostics = diagnostics,
    measured_form = measured$signature,
    source_plan = measured$source_plan,
    source_receipt = measured$source_receipt,
    reference_signature = reference_signature
  )
  structure(c(fields, list(
    signature = .tomography_signature("reconstruction", fields)
  )), class = "effect_tomography_result")
}

.validate_tomography_result <- function(x) {
  expected <- c("operator", "method", "status", "lossless", "certified",
    "left_projection", "right_projection", "diagnostics", "measured_form",
    "source_plan", "source_receipt", "reference_signature", "signature")
  if (!inherits(x, "effect_tomography_result") || !is.list(x) ||
      !identical(names(x), expected) || !is.matrix(x$operator) ||
      any(!is.finite(x$operator)) || !is.logical(x$lossless) ||
      length(x$lossless) != 1L || !is.logical(x$certified) ||
      length(x$certified) != 1L ||
      !x$method %in% c("parseval", "canonical_dual",
        "projected_pseudoinverse") ||
      !x$status %in% c("exact_algebraic_reconstruction",
        "numerically_certified_exact_reconstruction",
        "projected_reconstruction",
        "numerically_certified_projected_reconstruction")) {
    stop("Tomography result is missing or noncanonical.", call. = FALSE)
  }
  fields <- x[setdiff(names(x), "signature")]
  if (!identical(x$signature,
      .tomography_signature("reconstruction", fields))) {
    stop("Tomography-result identity is inconsistent.", call. = FALSE)
  }
  invisible(x)
}

.reconstruct_neural_evidence <- function(
    form, left_frame, right_frame = left_frame,
    method = c("auto", "parseval", "canonical_dual",
               "projected_pseudoinverse"),
    tolerance = 1e-10, max_condition = 1e8,
    allow_projection = TRUE, workspace_bytes = NULL,
    reference_operator = NULL) {
  method <- match.arg(method)
  if (!is.numeric(tolerance) || length(tolerance) != 1L ||
      is.na(tolerance) || !is.finite(tolerance) || tolerance <= 0 ||
      !is.numeric(max_condition) || length(max_condition) != 1L ||
      is.na(max_condition) || !is.finite(max_condition) ||
      max_condition < 1 || !is.logical(allow_projection) ||
      length(allow_projection) != 1L || is.na(allow_projection)) {
    stop("Tomography tolerances, condition limit, or projection policy are invalid.",
      call. = FALSE)
  }
  .validate_measurement_form(form, probe = FALSE)
  left_frame <- .validate_measurement_frame(left_frame)
  right_frame <- .validate_measurement_frame(right_frame)
  resource <- .tomography_resource_plan(
    form, left_frame, right_frame, workspace_bytes
  )
  .tomography_require_budget(resource)
  if (!isTRUE(form$capabilities$complete_edge_set)) {
    .capability_refusal(sprintf(paste0(
      "Reconstructing the global neural operator requires a block for every ",
      "ordered node pair, and this form carries %s (edge completeness: %s). ",
      "A diagonal-only or requested-edge map has discarded the between-node ",
      "directions, so no lossless reconstruction exists to return."
    ), .msg_count(nrow(form$block_index), "edge block"),
      form$edge_completeness),
      capability = "complete_edge_set",
      namespace = "tomography",
      reasons = "edge_set_is_not_frame_complete",
      remedies = paste0(
        "Build the form over every directed node pair, for example with ",
        "`edge_frame()` on `expand.grid(from = nodes$node_ids, ",
        "to = nodes$node_ids)`."
      )
    )
  }
  left <- .tomography_stack_frame(left_frame)
  right <- .tomography_stack_frame(right_frame)
  left_diagnostics <- .tomography_frame_diagnostics(
    left$operator, tolerance
  )
  right_diagnostics <- .tomography_frame_diagnostics(
    right$operator, tolerance
  )
  frame_diagnostics <- list(
    left = left_diagnostics,
    right = right_diagnostics,
    edge_completeness = form$edge_completeness,
    regularization = list(
      kind = "truncated_svd",
      relative_cutoff = tolerance,
      ridge = 0
    ),
    max_condition = max_condition,
    resource_plan = resource$signature
  )
  full_rank <- left_diagnostics$full_column_rank &&
    right_diagnostics$full_column_rank
  used_method <- method
  if (method == "auto") {
    parseval <- left_diagnostics$parseval_residual <= tolerance &&
      right_diagnostics$parseval_residual <= tolerance
    used_method <- if (full_rank && parseval) "parseval" else
      if (full_rank) "canonical_dual" else "projected_pseudoinverse"
  }
  if (used_method == "parseval" &&
      (!full_rank || left_diagnostics$parseval_residual > tolerance ||
       right_diagnostics$parseval_residual > tolerance)) {
    .tomography_reject(
      "Parseval reconstruction requires two full-rank Parseval frames.",
      frame_diagnostics
    )
  }
  if (used_method == "canonical_dual" && !full_rank) {
    .tomography_reject(
      "Canonical-dual lossless reconstruction requires full column rank.",
      frame_diagnostics
    )
  }
  if (used_method == "projected_pseudoinverse" && full_rank) {
    .tomography_reject(
      "Projected reconstruction is reserved for rank-deficient frames.",
      frame_diagnostics
    )
  }
  if (used_method == "projected_pseudoinverse" && !allow_projection) {
    .tomography_reject(
      "Rank-deficient tomography was rejected because projection is disabled.",
      frame_diagnostics
    )
  }
  conditions <- if (full_rank) {
    c(left_diagnostics$condition_number,
      right_diagnostics$condition_number)
  } else {
    c(left_diagnostics$retained_condition_number,
      right_diagnostics$retained_condition_number)
  }
  if (any(conditions > max_condition)) {
    .tomography_reject(
      "Tomographic reconstruction is too ill-conditioned for the declared limit.",
      frame_diagnostics
    )
  }
  measured <- .assemble_measurement_blocks(
    form, left_frame, right_frame, workspace_bytes
  )
  if (used_method == "parseval") {
    left_inverse <- t(left$operator)
    right_inverse <- t(right$operator)
  } else {
    left_inverse <- .tomography_pseudoinverse(left$operator, tolerance)
    right_inverse <- .tomography_pseudoinverse(right$operator, tolerance)
  }
  operator <- left_inverse %*% measured$value %*% t(right_inverse)
  left_projection <- left_inverse %*% left$operator
  right_projection <- right_inverse %*% right$operator
  lossless <- full_rank && used_method != "projected_pseudoinverse"
  reference_signature <- NULL
  residual <- relative_residual <- NULL
  certified <- !is.null(reference_operator)
  if (certified) {
    if (!is.matrix(reference_operator) || !is.numeric(reference_operator) ||
        !identical(dim(reference_operator), dim(operator)) ||
        any(!is.finite(reference_operator))) {
      stop("Tomography reference operator has invalid axes or values.",
        call. = FALSE)
    }
    expected <- left_projection %*% reference_operator %*% right_projection
    residual <- sqrt(sum((operator - expected)^2))
    scale <- max(sqrt(sum(expected^2)), tolerance)
    relative_residual <- residual / scale
    reference_signature <- .sha256_signature(unname(reference_operator))
    if (!is.finite(relative_residual) || relative_residual > tolerance) {
      diagnostics <- c(frame_diagnostics, list(
        reconstruction_residual = residual,
        relative_reconstruction_residual = relative_residual
      ))
      .tomography_reject(
        "Tomographic reconstruction failed its numerical reference check.",
        diagnostics
      )
    }
  }
  diagnostics <- c(frame_diagnostics, list(
    reconstruction_residual = residual,
    relative_reconstruction_residual = relative_residual,
    left_projection_residual = max(abs(
      left_projection %*% left_projection - left_projection
    )),
    right_projection_residual = max(abs(
      right_projection %*% right_projection - right_projection
    ))
  ))
  status <- if (lossless) {
    if (certified) "numerically_certified_exact_reconstruction" else
      "exact_algebraic_reconstruction"
  } else {
    if (certified) "numerically_certified_projected_reconstruction" else
      "projected_reconstruction"
  }
  result <- .new_tomography_result(
    operator, used_method, status, lossless, certified,
    left_projection, right_projection, diagnostics, measured,
    reference_signature
  )
  .validate_tomography_result(result)
  result
}
