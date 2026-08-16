# Identified neural measurement coordinates --------------------------------

.measurement_identifier <- function(x, label) {
  if (!is.character(x) || length(x) != 1L || is.na(x) || !nzchar(x)) {
    stop(sprintf("%s must be one nonempty identifier.", label), call. = FALSE)
  }
  invisible(x)
}

.measurement_axis <- function(coordinates, id, basis_id = id,
                              units = "arbitrary", provenance = list()) {
  if (!is.character(coordinates) || length(coordinates) < 1L ||
      anyNA(coordinates) || any(!nzchar(coordinates)) ||
      anyDuplicated(coordinates)) {
    stop("Measurement coordinates must be unique nonempty identifiers.",
      call. = FALSE)
  }
  .measurement_identifier(id, "Measurement-axis `id`")
  .measurement_identifier(basis_id, "Measurement-axis `basis_id`")
  if (!is.character(units) || !length(units) %in% c(1L, length(coordinates)) ||
      anyNA(units) || any(!nzchar(units))) {
    stop("Measurement-axis units must contain one unit or one per coordinate.",
      call. = FALSE)
  }
  if (length(units) == 1L) units <- rep(units, length(coordinates))
  .validate_effect_provenance(provenance, "measurement-axis provenance")
  semantic <- list(
    schema_version = 1L,
    id = id,
    coordinates = unname(coordinates),
    basis_id = basis_id,
    units = unname(units),
    provenance = provenance
  )
  structure(c(semantic[-1L], list(
    signature = paste0("sha256:", digest::digest(
      semantic, algo = "sha256", serialize = TRUE, serializeVersion = 2L
    ))
  )), class = "effect_measurement_axis")
}

.validate_measurement_axis <- function(x) {
  expected <- c("id", "coordinates", "basis_id", "units", "provenance",
    "signature")
  if (!inherits(x, "effect_measurement_axis") || !is.list(x) ||
      !identical(names(x), expected)) {
    stop("Measurement-axis fields are missing or noncanonical.",
      call. = FALSE)
  }
  rebuilt <- .measurement_axis(
    x$coordinates, x$id, x$basis_id, x$units, x$provenance
  )
  if (!identical(x, rebuilt)) {
    stop("Measurement-axis identity is inconsistent.", call. = FALSE)
  }
  rebuilt
}

.measurement_leg <- function(operator, source_domain, output_space,
                             support = NULL,
                             estimation = c("fixed", "learned_frozen"),
                             provenance = list(), decomposition = NULL) {
  if (!is.matrix(operator) || !is.numeric(operator) ||
      any(dim(operator) < 1L) || any(!is.finite(operator))) {
    stop("A measurement leg must be a finite nonempty numeric matrix.",
      call. = FALSE)
  }
  operator <- unname(operator)
  storage.mode(operator) <- "double"
  source_domain <- .domain_reference(source_domain)
  output_space <- .validate_measurement_axis(output_space)
  estimation <- match.arg(estimation)
  .validate_effect_provenance(provenance, "measurement-leg provenance")
  if (ncol(operator) != source_domain$n_features ||
      nrow(operator) != length(output_space$coordinates)) {
    stop("Measurement-leg dimensions do not match their identified spaces.",
      call. = FALSE)
  }
  nonzero <- colSums(abs(operator)) > 0
  canonical_support <- source_domain$feature_ids[nonzero]
  if (is.null(support)) support <- canonical_support
  if (length(support) < 1L || anyNA(support) || anyDuplicated(support) ||
      !identical(support, canonical_support)) {
    stop("Measurement-leg support must exactly identify its nonzero source columns.",
      call. = FALSE)
  }
  if (estimation == "learned_frozen") {
    signature <- provenance$training_signature
    if (!identical(provenance$frozen, TRUE) ||
        !is.character(signature) || length(signature) != 1L ||
        is.na(signature) ||
        !grepl("^sha256:[[:xdigit:]]{64}$", signature)) {
      stop(paste0(
        "Learned measurement legs require frozen training provenance and ",
        "a training signature."
      ), call. = FALSE)
    }
  }
  if (!is.null(decomposition)) {
    decomposition <- .validate_measurement_decomposition(decomposition)
    if (!identical(decomposition$output_space, output_space)) {
      stop("Measurement decomposition and leg output spaces differ.",
        call. = FALSE)
    }
  }
  semantic <- list(
    schema_version = 1L,
    kind = "linear",
    source_domain = source_domain,
    output_space = output_space,
    operator = operator,
    support = support,
    estimation = estimation,
    provenance = provenance,
    decomposition = decomposition
  )
  structure(c(semantic[-1L], list(
    signature = paste0("sha256:", digest::digest(
      semantic, algo = "sha256", serialize = TRUE, serializeVersion = 2L
    ))
  )), class = "effect_measurement_leg")
}

.validate_measurement_leg <- function(x, source_domain = NULL) {
  expected <- c("kind", "source_domain", "output_space", "operator",
    "support", "estimation", "provenance", "decomposition", "signature")
  if (!inherits(x, "effect_measurement_leg") || !is.list(x) ||
      !identical(names(x), expected) || !identical(x$kind, "linear")) {
    stop("Measurement-leg fields are missing or noncanonical.",
      call. = FALSE)
  }
  rebuilt <- .measurement_leg(
    x$operator, x$source_domain, x$output_space, x$support,
    x$estimation, x$provenance, x$decomposition
  )
  if (!identical(x, rebuilt)) {
    stop("Measurement-leg identity is inconsistent.", call. = FALSE)
  }
  if (!is.null(source_domain) &&
      !.same_domain_reference(x$source_domain,
        .domain_reference(source_domain))) {
    stop("Measurement leg belongs to a different neural source space.",
      call. = FALSE)
  }
  rebuilt
}

.measurement_leg_metric <- function(x) {
  x <- .validate_measurement_leg(x)
  crossprod(x$operator)
}

.measurement_decomposition <- function(id, output_space, components,
                                       ranges, provenance = list()) {
  .measurement_identifier(id, "Measurement-decomposition `id`")
  output_space <- .validate_measurement_axis(output_space)
  if (!is.character(components) || length(components) < 1L ||
      anyNA(components) || any(!nzchar(components)) ||
      anyDuplicated(components) || !is.list(ranges) ||
      length(ranges) != length(components) ||
      is.null(names(ranges)) || !identical(names(ranges), components) ||
      !all(vapply(ranges, function(x) {
        is.integer(x) && !anyNA(x) &&
          (length(x) == 0L || identical(x, seq.int(min(x), max(x))))
      }, logical(1)))) {
    stop("Measurement decomposition requires named contiguous component ranges.",
      call. = FALSE)
  }
  covered <- unlist(ranges, use.names = FALSE)
  if (!identical(as.integer(covered), seq_along(output_space$coordinates))) {
    stop("Measurement-decomposition ranges must partition the output axis in order.",
      call. = FALSE)
  }
  .validate_effect_provenance(provenance,
    "measurement-decomposition provenance")
  semantic <- list(
    schema_version = 1L,
    id = id,
    output_space = output_space,
    components = components,
    ranges = ranges,
    orientation = stats::setNames(rep("oriented", length(components)),
      components),
    provenance = provenance
  )
  structure(c(semantic[-1L], list(
    signature = paste0("sha256:", digest::digest(
      semantic, algo = "sha256", serialize = TRUE, serializeVersion = 2L
    ))
  )), class = "effect_measurement_decomposition")
}

.validate_measurement_decomposition <- function(x) {
  expected <- c("id", "output_space", "components", "ranges", "orientation",
    "provenance", "signature")
  if (!inherits(x, "effect_measurement_decomposition") || !is.list(x) ||
      !identical(names(x), expected)) {
    stop("Measurement-decomposition fields are missing or noncanonical.",
      call. = FALSE)
  }
  rebuilt <- .measurement_decomposition(
    x$id, x$output_space, x$components, x$ranges, x$provenance
  )
  if (!is.character(x$orientation) ||
      !identical(names(x$orientation), x$components) ||
      any(!x$orientation %in% c("oriented", "subspace_only"))) {
    stop("Measurement-decomposition orientation metadata is invalid.",
      call. = FALSE)
  }
  rebuilt$orientation <- x$orientation
  semantic <- list(
    schema_version = 1L,
    id = rebuilt$id,
    output_space = rebuilt$output_space,
    components = rebuilt$components,
    ranges = rebuilt$ranges,
    orientation = rebuilt$orientation,
    provenance = rebuilt$provenance
  )
  rebuilt$signature <- paste0("sha256:", digest::digest(
    semantic, algo = "sha256", serialize = TRUE, serializeVersion = 2L
  ))
  if (!identical(x, rebuilt)) {
    stop("Measurement-decomposition identity is inconsistent.",
      call. = FALSE)
  }
  rebuilt
}

.measurement_decomposition_with_orientation <- function(
    decomposition, orientation) {
  decomposition <- .validate_measurement_decomposition(decomposition)
  if (!is.character(orientation) ||
      !identical(names(orientation), decomposition$components) ||
      anyNA(orientation) ||
      any(!orientation %in% c("oriented", "subspace_only"))) {
    stop("Every decomposition component requires an orientation status.",
      call. = FALSE)
  }
  decomposition$orientation <- orientation
  semantic <- list(
    schema_version = 1L,
    id = decomposition$id,
    output_space = decomposition$output_space,
    components = decomposition$components,
    ranges = decomposition$ranges,
    orientation = decomposition$orientation,
    provenance = decomposition$provenance
  )
  decomposition$signature <- paste0("sha256:", digest::digest(
    semantic, algo = "sha256", serialize = TRUE, serializeVersion = 2L
  ))
  .validate_measurement_decomposition(decomposition)
}

.direct_sum_measurement_leg <- function(components, id,
                                        basis_id = paste0(id, ":direct-sum"),
                                        provenance = list()) {
  if (!is.list(components) || length(components) < 1L ||
      is.null(names(components)) || anyNA(names(components)) ||
      any(!nzchar(names(components))) || anyDuplicated(names(components))) {
    stop("Direct sums require uniquely named measurement-leg components.",
      call. = FALSE)
  }
  components <- lapply(components, .validate_measurement_leg)
  source <- components[[1L]]$source_domain
  if (!all(vapply(components, function(x) {
    .same_domain_reference(x$source_domain, source)
  }, logical(1)))) {
    stop("Direct-sum components must share one exact neural source space.",
      call. = FALSE)
  }
  coordinates <- unlist(Map(function(name, leg) {
    paste0(name, "::", leg$output_space$coordinates)
  }, names(components), components), use.names = FALSE)
  units <- unlist(lapply(components, function(x) x$output_space$units),
    use.names = FALSE)
  component_signatures <- vapply(components, `[[`, character(1), "signature")
  axis <- .measurement_axis(
    coordinates, id = id, basis_id = basis_id, units = units,
    provenance = list(
      construction = "direct_sum",
      components = unname(component_signatures),
      user = provenance
    )
  )
  widths <- vapply(components, function(x) nrow(x$operator), integer(1))
  ends <- cumsum(widths)
  starts <- c(1L, utils::head(ends, -1L) + 1L)
  ranges <- Map(seq.int, starts, ends)
  names(ranges) <- names(components)
  decomposition <- .measurement_decomposition(
    id, axis, names(components), ranges,
    provenance = list(component_signatures = unname(component_signatures))
  )
  .measurement_leg(
    do.call(rbind, lapply(components, `[[`, "operator")),
    source, axis,
    support = source$feature_ids[
      colSums(abs(do.call(rbind, lapply(components, `[[`, "operator")))) > 0
    ],
    estimation = if (all(vapply(components, function(x) {
      identical(x$estimation, "fixed")
    }, logical(1)))) "fixed" else "learned_frozen",
    provenance = if (all(vapply(components, function(x) {
      identical(x$estimation, "fixed")
    }, logical(1)))) {
      list(construction = "direct_sum", components = unname(component_signatures))
    } else {
      training <- paste0("sha256:", digest::digest(
        component_signatures, algo = "sha256", serialize = TRUE, serializeVersion = 2L
      ))
      list(construction = "direct_sum", components = unname(component_signatures),
        frozen = TRUE, training_signature = training)
    },
    decomposition = decomposition
  )
}

.measurement_frame <- function(legs,
                               injectivity = c("not_guaranteed",
                                                "positive_diagonal_coverage")) {
  injectivity <- match.arg(injectivity)
  if (!is.list(legs) || length(legs) < 1L || is.null(names(legs)) ||
      anyNA(names(legs)) || any(!nzchar(names(legs))) ||
      anyDuplicated(names(legs))) {
    stop("Measurement frames require uniquely named ordered legs.",
      call. = FALSE)
  }
  legs <- lapply(legs, .validate_measurement_leg)
  source <- legs[[1L]]$source_domain
  if (!all(vapply(legs, function(x) {
    .same_domain_reference(x$source_domain, source)
  }, logical(1)))) {
    stop("Every measurement-frame leg must share one exact source space.",
      call. = FALSE)
  }
  operator <- Reduce(`+`, lapply(legs, .measurement_leg_metric))
  operator <- unname(operator)
  coverage_count <- Reduce(`+`, lapply(legs, function(x) {
    as.integer(source$feature_ids %in% x$support)
  }))
  coverage <- list(
    feature_ids = source$feature_ids,
    count = as.integer(coverage_count),
    complete = all(coverage_count > 0L)
  )
  if (injectivity == "positive_diagonal_coverage") {
    off_diagonal <- operator - diag(diag(operator), nrow(operator))
    if (max(abs(off_diagonal)) > 1e-12 || any(diag(operator) <= 0)) {
      stop(paste0(
        "Positive-diagonal injectivity requires a strictly positive diagonal ",
        "frame operator."
      ), call. = FALSE)
    }
    injectivity_record <- list(
      guaranteed = TRUE,
      construction = injectivity,
      lower_bound = min(diag(operator))
    )
  } else {
    injectivity_record <- list(
      guaranteed = FALSE,
      construction = injectivity,
      lower_bound = NULL
    )
  }
  operator_signature <- paste0("sha256:", digest::digest(
    operator, algo = "sha256", serialize = TRUE, serializeVersion = 2L
  ))
  decomposition_identity <- vapply(legs, function(x) {
    if (is.null(x$decomposition)) NA_character_ else x$decomposition$signature
  }, character(1))
  dual <- list(
    kind = "canonical_dual",
    eligible = injectivity_record$guaranteed,
    materialized = FALSE,
    frame_operator = operator_signature,
    stability = "not_certified",
    solver_policy = NULL
  )
  semantic <- list(
    schema_version = 1L,
    source_domain = source,
    node_ids = names(legs),
    leg_signatures = unname(vapply(legs, `[[`, character(1), "signature")),
    frame_operator = operator,
    coverage = coverage,
    decomposition_identity = unname(decomposition_identity),
    injectivity = injectivity_record,
    dual = dual
  )
  structure(c(semantic[c("source_domain", "node_ids")], list(
    legs = legs
  ), semantic[c("frame_operator", "coverage", "decomposition_identity",
    "injectivity", "dual")], list(
    signature = paste0("sha256:", digest::digest(
      semantic, algo = "sha256", serialize = TRUE, serializeVersion = 2L
    ))
  )), class = "effect_measurement_frame")
}

.validate_measurement_frame <- function(x) {
  expected <- c("source_domain", "node_ids", "legs", "frame_operator",
    "coverage", "decomposition_identity", "injectivity", "dual", "signature")
  if (!inherits(x, "effect_measurement_frame") || !is.list(x) ||
      !identical(names(x), expected) ||
      !identical(names(x$legs), x$node_ids)) {
    stop("Measurement-frame fields are missing or noncanonical.",
      call. = FALSE)
  }
  mode <- x$injectivity$construction
  rebuilt <- .measurement_frame(x$legs, injectivity = mode)
  if (!identical(x, rebuilt)) {
    stop("Measurement-frame identity is inconsistent.", call. = FALSE)
  }
  rebuilt
}

.measurement_edges <- function(left, right, left_frame, right_frame = left_frame,
                               weight = NULL) {
  left_frame <- .validate_measurement_frame(left_frame)
  right_frame <- .validate_measurement_frame(right_frame)
  if (!is.character(left) || !is.character(right) || length(left) < 1L ||
      length(left) != length(right) || anyNA(left) || anyNA(right) ||
      any(!nzchar(left)) || any(!nzchar(right)) ||
      any(!left %in% left_frame$node_ids) ||
      any(!right %in% right_frame$node_ids) ||
      anyDuplicated(paste(left, right, sep = "\r"))) {
    stop("Measurement edges must be unique explicit ordered frame-node pairs.",
      call. = FALSE)
  }
  weighted <- !is.null(weight)
  if (is.null(weight)) weight <- rep(1, length(left))
  if (!is.numeric(weight) || length(weight) != length(left) || anyNA(weight) ||
      any(!is.finite(weight))) {
    stop("Measurement-edge weights must be finite when supplied.",
      call. = FALSE)
  }
  table <- data.frame(
    left = unname(left), right = unname(right), weight = unname(weight),
    stringsAsFactors = FALSE
  )
  semantic <- list(
    schema_version = 1L,
    left_frame = left_frame$signature,
    right_frame = right_frame$signature,
    weighted = weighted,
    edges = unclass(table)
  )
  structure(list(
    left_frame = left_frame$signature,
    right_frame = right_frame$signature,
    weighted = weighted,
    edges = table,
    signature = paste0("sha256:", digest::digest(
      semantic, algo = "sha256", serialize = TRUE, serializeVersion = 2L
    ))
  ), class = "effect_measurement_edges")
}

.validate_measurement_edges <- function(x, left_frame, right_frame = left_frame) {
  expected <- c("left_frame", "right_frame", "weighted", "edges", "signature")
  if (!inherits(x, "effect_measurement_edges") || !is.list(x) ||
      !identical(names(x), expected) || !is.logical(x$weighted) ||
      length(x$weighted) != 1L || is.na(x$weighted) ||
      !is.data.frame(x$edges) ||
      !identical(names(x$edges), c("left", "right", "weight"))) {
    stop("Measurement-edge fields are missing or noncanonical.",
      call. = FALSE)
  }
  rebuilt <- .measurement_edges(
    x$edges$left, x$edges$right, left_frame, right_frame,
    if (x$weighted) x$edges$weight else NULL
  )
  if (!identical(x, rebuilt)) {
    stop("Measurement-edge identity is inconsistent.", call. = FALSE)
  }
  rebuilt
}

.reverse_measurement_edges <- function(x, left_frame, right_frame = left_frame) {
  x <- .validate_measurement_edges(x, left_frame, right_frame)
  .measurement_edges(
    x$edges$right, x$edges$left, right_frame, left_frame,
    if (x$weighted) x$edges$weight else NULL
  )
}

.canonical_additive_weights <- function(weights) {
  if (inherits(weights, "Matrix")) {
    value <- methods::as(weights, "CsparseMatrix")
    value <- methods::as(value, "generalMatrix")
    return(Matrix::drop0(value))
  }
  if (!is.matrix(weights) || !is.numeric(weights)) {
    stop("Additive weights must be a numeric matrix or sparse Matrix.",
      call. = FALSE)
  }
  unname(weights)
}

.canonical_additive_weight_identity <- function(weights) {
  weights <- .canonical_additive_weights(weights)
  if (!inherits(weights, "Matrix")) {
    weights <- Matrix::Matrix(weights, sparse = TRUE)
    weights <- .canonical_additive_weights(weights)
  }
  list(
    storage = "canonical_general_csparse",
    i = weights@i,
    p = weights@p,
    x = weights@x,
    dim = weights@Dim
  )
}

.additive_weight_row <- function(weights, row) {
  weights <- .canonical_additive_weights(weights)
  if (!is.numeric(row) || length(row) != 1L || is.na(row) ||
      !is.finite(row) || row %% 1 != 0 || row < 1L || row > nrow(weights)) {
    stop("An additive weight row must identify one measurement.",
      call. = FALSE)
  }
  value <- weights[as.integer(row), , drop = FALSE]
  if (!inherits(value, "Matrix")) {
    return(unname(value[1L, , drop = TRUE]))
  }
  entries <- Matrix::summary(value)
  dense <- numeric(ncol(weights))
  if (nrow(entries)) dense[entries$j] <- entries$x
  dense
}

.additive_frame_signature <- function(frame) {
  .validate_frame_for_compile(frame)
  semantic <- list(
    representation = frame$representation,
    weights = .canonical_additive_weight_identity(frame$weights),
    normalization = frame$normalization,
    domain = frame$domain,
    index = frame$index,
    specification = frame$specification
  )
  paste0("sha256:", digest::digest(semantic, algo = "sha256", serialize = TRUE, serializeVersion = 2L))
}

.measurement_frame_from_additive <- function(
    frame, component = c("total", "coherent")) {
  component <- match.arg(component)
  .validate_frame_for_compile(frame)
  if (!identical(frame$representation, "additive_diagonal")) {
    stop("Only additive diagonal frames have a canonical square-root adapter.",
      call. = FALSE)
  }
  weights <- .canonical_additive_weights(frame$weights)
  domain <- frame$domain
  frame_signature <- .additive_frame_signature(frame)
  node_ids <- if (!is.null(frame$index) &&
      is.data.frame(frame$index) && "measurement" %in% names(frame$index)) {
    as.character(frame$index$measurement)
  } else {
    paste0("measurement_", seq_len(nrow(weights)))
  }
  if (anyNA(node_ids) || any(!nzchar(node_ids)) || anyDuplicated(node_ids)) {
    stop("Additive-frame measurement identities must be unique and nonempty.",
      call. = FALSE)
  }
  legs <- lapply(seq_len(nrow(weights)), function(index) {
    weight <- .additive_weight_row(weights, index)
    support_index <- which(weight > 0)
    node <- node_ids[[index]]
    if (component == "total") {
      operator <- matrix(0, length(support_index), ncol(weights))
      operator[cbind(seq_along(support_index), support_index)] <-
        sqrt(weight[support_index])
      output <- .measurement_axis(
        paste0("feature::", domain$feature_ids[support_index]),
        id = paste0(frame_signature, ":", node, ":sqrt"),
        basis_id = paste0("coordinate-preserving:", domain$signature),
        units = "sqrt_spatial_weight",
        provenance = list(
          adapter = "additive_square_root",
          frame = frame_signature,
          node = node
        )
      )
    } else {
      operator <- matrix(weight / sqrt(sum(weight)), nrow = 1L)
      output <- .measurement_axis(
        "common_mode",
        id = paste0(frame_signature, ":", node, ":coherent"),
        basis_id = paste0("weighted-common-mode:", domain$signature),
        units = "weighted_common_mode",
        provenance = list(
          adapter = "additive_coherent",
          frame = frame_signature,
          node = node
        )
      )
    }
    .measurement_leg(
      operator, domain, output,
      support = domain$feature_ids[support_index],
      estimation = "fixed",
      provenance = list(
        construction = if (component == "total") {
          "coordinate_preserving_square_root"
        } else {
          "coherent_common_mode"
        },
        frame = frame_signature,
        node = node
      )
    )
  })
  names(legs) <- node_ids
  coverage <- if (inherits(weights, "Matrix")) {
    Matrix::colSums(weights)
  } else {
    colSums(weights)
  }
  injectivity <- if (component == "total" && all(coverage > 0)) {
    "positive_diagonal_coverage"
  } else {
    "not_guaranteed"
  }
  .measurement_frame(legs, injectivity)
}
