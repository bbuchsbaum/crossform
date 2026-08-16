# Crossed lifting of node measurement decompositions -----------------------

.coherent_configuration_operator <- function(weight, feature_ids, frame_signature,
                                              node, domain) {
  support <- which(weight > 0)
  mass <- sum(weight)
  root_weight <- sqrt(weight[support])
  coherent_direction <- root_weight / sqrt(mass)
  candidate <- cbind(coherent_direction, diag(length(support)))
  basis <- qr.Q(qr(candidate), complete = TRUE)
  if (sum(basis[, 1L] * coherent_direction) < 0) {
    basis[, 1L] <- -basis[, 1L]
  }
  local_operator <- t(basis) %*% diag(root_weight, nrow = length(support))
  operator <- matrix(0, length(support), length(weight))
  operator[, support] <- local_operator
  coordinates <- c("common_mode",
    if (length(support) > 1L) {
      paste0("configuration_mode_", seq_len(length(support) - 1L))
    } else {
      character()
    })
  output <- .measurement_axis(
    coordinates,
    id = paste0(frame_signature, ":", node, ":coherent-configuration"),
    basis_id = paste0("coherent-plus-configuration-subspace:",
      domain$signature),
    units = "sqrt_spatial_weight",
    provenance = list(
      adapter = "coherent_configuration_orthogonal_factor",
      frame = frame_signature,
      node = node
    )
  )
  ranges <- list(
    coherent = 1L,
    configuration = if (length(support) > 1L) {
      as.integer(2:length(support))
    } else {
      integer()
    }
  )
  decomposition <- .measurement_decomposition(
    paste0(frame_signature, ":", node, ":coherent-configuration"),
    output,
    c("coherent", "configuration"),
    ranges,
    provenance = list(
      coherent_orientation = "positive_weighted_common_mode",
      configuration_identity = "weighted_orthogonal_complement"
    )
  )
  decomposition <- .measurement_decomposition_with_orientation(
    decomposition,
    c(coherent = "oriented", configuration = "subspace_only")
  )
  .measurement_leg(
    operator, domain, output,
    support = feature_ids[support],
    estimation = "fixed",
    provenance = list(
      construction = "coherent_configuration_orthogonal_factor",
      frame = frame_signature,
      node = node
    ),
    decomposition = decomposition
  )
}

.measurement_frame_from_additive_decomposition <- function(frame) {
  .validate_frame_for_compile(frame)
  if (!identical(frame$representation, "additive_diagonal")) {
    stop("Coherent/configuration adaptation requires an additive diagonal frame.",
      call. = FALSE)
  }
  weights <- .canonical_additive_weights(frame$weights)
  domain <- frame$domain
  frame_signature <- .additive_frame_signature(frame)
  node_ids <- if (!is.null(frame$index) && is.data.frame(frame$index) &&
      "measurement" %in% names(frame$index)) {
    as.character(frame$index$measurement)
  } else {
    paste0("measurement_", seq_len(nrow(weights)))
  }
  legs <- Map(function(index, node) {
    weight <- .additive_weight_row(weights, index)
    .coherent_configuration_operator(
      weight, domain$feature_ids, frame_signature, node, domain
    )
  }, seq_len(nrow(weights)), node_ids)
  names(legs) <- node_ids
  coverage <- if (inherits(weights, "Matrix")) {
    Matrix::colSums(weights)
  } else {
    colSums(weights)
  }
  injectivity <- if (all(coverage > 0)) {
    "positive_diagonal_coverage"
  } else {
    "not_guaranteed"
  }
  .measurement_frame(legs, injectivity)
}

.measurement_component_ranges <- function(frame_manifest, node) {
  frame_manifest <- .validate_measurement_frame_manifest(frame_manifest)
  if (!node %in% frame_manifest$node_ids) {
    stop("Measurement-decomposition node is absent from its frame.",
      call. = FALSE)
  }
  decomposition <- frame_manifest$legs[[node]]$decomposition
  if (is.null(decomposition)) {
    dimension <- length(frame_manifest$legs[[node]]$output_space$coordinates)
    return(list(
      id = paste0(frame_manifest$legs[[node]]$signature, ":whole"),
      components = "whole",
      ranges = list(whole = seq_len(dimension)),
      orientation = c(whole = "oriented"),
      signature = frame_manifest$legs[[node]]$signature
    ))
  }
  decomposition <- .validate_measurement_decomposition(decomposition)
  list(
    id = decomposition$id,
    components = decomposition$components,
    ranges = decomposition$ranges,
    orientation = decomposition$orientation,
    signature = decomposition$signature
  )
}

.lift_measurement_decomposition <- function(x, edge) {
  .validate_measurement_form(x, probe = FALSE)
  position <- .measurement_edge_position(x$block_index, edge)
  row <- x$block_index[position, , drop = FALSE]
  left <- .measurement_component_ranges(x$left_frame, row$left[[1L]])
  right <- .measurement_component_ranges(x$right_frame, row$right[[1L]])
  value <- .measurement_block(x, position)
  component_index <- expand.grid(
    left_component = left$components,
    right_component = right$components,
    stringsAsFactors = FALSE
  )
  component_index$left_orientation <- left$orientation[
    component_index$left_component
  ]
  component_index$right_orientation <- right$orientation[
    component_index$right_component
  ]
  component_index$raw_interpretation <-
    component_index$left_orientation == "oriented" &
    component_index$right_orientation == "oriented"
  blocks <- stats::setNames(lapply(seq_len(nrow(component_index)), function(i) {
    left_range <- left$ranges[[component_index$left_component[[i]]]]
    right_range <- right$ranges[[component_index$right_component[[i]]]]
    value[left_range, right_range, drop = FALSE]
  }), paste(component_index$left_component,
    component_index$right_component, sep = "::"))
  semantic <- list(
    schema_version = 1L,
    form = x$contract_signature,
    edge_id = row$edge_id[[1L]],
    left_decomposition = left$signature,
    right_decomposition = right$signature,
    index = unclass(component_index)
  )
  structure(list(
    parent_form = x$contract_signature,
    edge_id = row$edge_id[[1L]],
    left_node = row$left[[1L]],
    right_node = row$right[[1L]],
    left_decomposition = left,
    right_decomposition = right,
    index = component_index,
    blocks = blocks,
    original_dimension = dim(value),
    signature = .sha256_signature(semantic)
  ), class = "effect_measurement_decomposition_view")
}

.validate_measurement_decomposition_view <- function(x) {
  expected <- c("parent_form", "edge_id", "left_node", "right_node",
    "left_decomposition", "right_decomposition", "index", "blocks",
    "original_dimension", "signature")
  if (!inherits(x, "effect_measurement_decomposition_view") || !is.list(x) ||
      !identical(names(x), expected) || !is.data.frame(x$index) ||
      !is.list(x$blocks) || nrow(x$index) != length(x$blocks) ||
      !identical(names(x$blocks), paste(
        x$index$left_component, x$index$right_component, sep = "::"
      ))) {
    stop("Measurement-decomposition view is missing or noncanonical.",
      call. = FALSE)
  }
  for (i in seq_len(nrow(x$index))) {
    expected_dimension <- c(
      length(x$left_decomposition$ranges[[x$index$left_component[[i]]]]),
      length(x$right_decomposition$ranges[[x$index$right_component[[i]]]])
    )
    if (!is.matrix(x$blocks[[i]]) || !is.numeric(x$blocks[[i]]) ||
        !identical(dim(x$blocks[[i]]), expected_dimension) ||
        any(!is.finite(x$blocks[[i]]))) {
      stop("A lifted component block has invalid dimensions or values.",
        call. = FALSE)
    }
  }
  semantic <- list(
    schema_version = 1L,
    form = x$parent_form,
    edge_id = x$edge_id,
    left_decomposition = x$left_decomposition$signature,
    right_decomposition = x$right_decomposition$signature,
    index = unclass(x$index)
  )
  expected_signature <- .sha256_signature(semantic)
  if (!identical(x$signature, expected_signature)) {
    stop("Measurement-decomposition view signature is inconsistent.",
      call. = FALSE)
  }
  invisible(x)
}

.recompose_measurement_decomposition <- function(x) {
  .validate_measurement_decomposition_view(x)
  value <- matrix(0, x$original_dimension[[1L]], x$original_dimension[[2L]])
  for (i in seq_len(nrow(x$index))) {
    left_range <- x$left_decomposition$ranges[[
      x$index$left_component[[i]]
    ]]
    right_range <- x$right_decomposition$ranges[[
      x$index$right_component[[i]]
    ]]
    value[left_range, right_range] <- x$blocks[[i]]
  }
  value
}

.measurement_component_block <- function(x, left_component, right_component,
                                         interpretation = c("raw", "coordinate")) {
  .validate_measurement_decomposition_view(x)
  interpretation <- match.arg(interpretation)
  position <- which(
    x$index$left_component == left_component &
      x$index$right_component == right_component
  )
  if (length(position) != 1L) {
    stop("Requested crossed decomposition component does not exist.",
      call. = FALSE)
  }
  if (interpretation == "raw" &&
      !isTRUE(x$index$raw_interpretation[[position]])) {
    stop(paste0(
      "Raw component entries require scientifically oriented bases; this ",
      "component is identified only as a subspace."
    ), call. = FALSE)
  }
  x$blocks[[position]]
}

.decomposition_inverse_sqrt <- function(value, tolerance, ridge = 0) {
  value <- (value + t(value)) / 2
  if (ridge > 0) value <- value + diag(ridge, nrow(value))
  decomposition <- eigen(value, symmetric = TRUE)
  if (min(decomposition$values) < -tolerance) {
    stop("Invariant canonical summaries require positive self-blocks.",
      call. = FALSE)
  }
  retained <- decomposition$values > tolerance
  scale <- numeric(length(decomposition$values))
  scale[retained] <- 1 / sqrt(decomposition$values[retained])
  decomposition$vectors %*% diag(scale, length(scale)) %*%
    t(decomposition$vectors)
}

.measurement_invariant_summary <- function(cross, left_self = NULL,
                                           right_self = NULL,
                                           tolerance = 1e-10, ridge = 0) {
  if (!is.matrix(cross) || !is.numeric(cross) || any(!is.finite(cross))) {
    stop("Invariant summaries require one finite cross block.",
      call. = FALSE)
  }
  singular_values <- svd(cross, nu = 0L, nv = 0L)$d
  canonical_values <- subspace_angles <- numeric()
  geometry_alignment <- NA_real_
  if (!is.null(left_self) || !is.null(right_self)) {
    if (is.null(left_self) || is.null(right_self) ||
        !is.matrix(left_self) || !is.matrix(right_self) ||
        !identical(dim(left_self), c(nrow(cross), nrow(cross))) ||
        !identical(dim(right_self), c(ncol(cross), ncol(cross))) ||
        any(!is.finite(left_self)) || any(!is.finite(right_self))) {
      stop("Canonical summaries require compatible finite self-blocks.",
        call. = FALSE)
    }
    normalized <- .decomposition_inverse_sqrt(
      left_self, tolerance, ridge
    ) %*% cross %*% .decomposition_inverse_sqrt(
      right_self, tolerance, ridge
    )
    raw_canonical <- svd(normalized, nu = 0L, nv = 0L)$d
    if (any(raw_canonical > 1 + 10 * tolerance)) {
      stop(paste0(
        "Canonical values exceed one beyond tolerance; the self-blocks are ",
        "inconsistent with the cross block and will not be silently clipped."
      ), call. = FALSE)
    }
    canonical_values <- pmin(1, pmax(0, raw_canonical))
    subspace_angles <- acos(canonical_values)
    denominator <- sqrt(sum(left_self^2) * sum(right_self^2))
    geometry_alignment <- if (denominator == 0) NA_real_ else
      sum(cross^2) / denominator
  }
  list(
    frobenius_norm = sqrt(sum(cross^2)),
    singular_values = singular_values,
    canonical_values = canonical_values,
    subspace_angles = subspace_angles,
    geometry_alignment = geometry_alignment,
    tolerance = tolerance,
    ridge = ridge,
    basis_invariance = "orthogonal_left_right"
  )
}
