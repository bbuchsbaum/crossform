# Complete geometry and query-only result contracts -------------------------
#
# `.effect_form_codec_format()`, the codec's recorded storage-format tag, is a
# leaf in R/primitives.R: the memory store here and the file store in
# R/storage.R both stamp it into their manifests, and a result file has no
# business owning a constant the executor's storage layer depends on.

.memory_geometry_store <- function(value, codec = "symmetric_packed") {
  if (!.is_finite_matrix(value)) {
    .input_error("Geometry components must be finite numeric matrices.")
  }
  format <- .effect_form_codec_format(codec)
  structure(
    list(
      dim = dim(value),
      representation = "memory",
      manifest = list(schema_version = 1L, complete = TRUE, dim = dim(value),
        format = format),
      read = function(rows = NULL) {
        if (is.null(rows)) value else value[rows, , drop = FALSE]
      }
    ),
    class = "effect_geometry_store"
  )
}

.block_geometry_store <- function(dim, read, codec = "symmetric_packed") {
  if (!.is_finite_numeric(dim) || length(dim) != 2L || any(dim < 0) ||
      any(dim %% 1 != 0)) {
    .input_error("A block store requires two nonnegative integer dimensions.")
  }
  if (!is.function(read)) {
    .input_error("A block store requires a `read` function.")
  }
  format <- .effect_form_codec_format(codec)
  structure(
    list(
      dim = as.integer(dim),
      representation = "block_backed",
      manifest = list(schema_version = 1L, complete = TRUE,
        dim = as.integer(dim), format = format),
      read = read
    ),
    class = "effect_geometry_store"
  )
}

.as_geometry_store <- function(x, codec = "symmetric_packed") {
  if (inherits(x, "effect_geometry_store")) x else .memory_geometry_store(x, codec)
}

.read_geometry_store <- function(store, rows = NULL) {
  if (!inherits(store, "effect_geometry_store")) {
    .input_error("Invalid geometry store.")
  }
  value <- store$read(rows)
  expected_rows <- if (is.null(rows)) store$dim[[1L]] else length(rows)
  if (!.is_finite_matrix(value) ||
      !identical(dim(value), c(as.integer(expected_rows), store$dim[[2L]]))) {
    .input_error(
      "Geometry store returned a block with invalid shape or values."
    )
  }
  value
}

.effect_result_signature <- function(result_capability, left_space, right_space,
                                     logical_shape, capabilities = NULL,
                                     codec = NULL, component = NULL,
                                     query = NULL, scientific_plan_id) {
  semantic <- list(
    schema_version = 1L,
    result_capability = result_capability,
    left_space = left_space$signature,
    right_space = right_space$signature,
    logical_shape = logical_shape,
    capabilities = capabilities,
    codec = codec,
    component = component,
    query = query,
    scientific_plan_id = scientific_plan_id
  )
  .sha256_signature(semantic)
}

# Universal complete effect-form result. This internal constructor establishes
# result/storage semantics before any two-relation production lowering exists.
effect_form <- function(total, left_space, right_space, receipt, index = NULL,
                        metadata = list(),
                        codec = c("rectangular", "symmetric_packed"),
                        symmetric = FALSE, guaranteed_psd = FALSE,
                        coherent = NULL, marginals = NULL) {
  codec <- match.arg(codec)
  left_space <- .as_effect_space(left_space)
  right_space <- .as_effect_space(right_space)
  self_form <- .same_effect_space(left_space, right_space)
  capabilities <- .effect_form_capabilities(
    self_form = self_form,
    symmetric = symmetric,
    guaranteed_psd = guaranteed_psd,
    coherent = !is.null(coherent)
  )
  total <- .as_geometry_store(total, codec)
  coherent <- if (is.null(coherent)) NULL else .as_geometry_store(coherent, codec)
  if (is.null(index)) index <- seq_len(total$dim[[1L]])
  if (!is.list(metadata)) .input_error("`metadata` must be a list.")
  .validate_execution_receipt(receipt)
  if (!is.null(metadata$scientific_plan_id) &&
      !identical(metadata$scientific_plan_id, receipt$scientific_plan_id)) {
    .contract_error(paste0(
      "Result metadata and execution receipt identify different scientific ",
      "plans."
    ))
  }
  metadata$scientific_plan_id <- receipt$scientific_plan_id
  storage <- unique(c(total$representation,
    if (is.null(coherent)) character() else coherent$representation))
  logical_shape <- as.integer(c(
    length(left_space$coordinates), length(right_space$coordinates)
  ))
  contract_signature <- .effect_result_signature(
    "complete_form", left_space, right_space, logical_shape,
    capabilities = capabilities, codec = codec,
    scientific_plan_id = receipt$scientific_plan_id
  )

  value <- structure(
    list(
      total = total,
      coherent = coherent,
      marginals = marginals,
      left_space = left_space,
      right_space = right_space,
      logical_shape = logical_shape,
      capabilities = capabilities,
      effect_space = if (self_form) left_space else NULL,
      effects = if (self_form) left_space$coordinates else NULL,
      index = index,
      metadata = metadata,
      receipt = receipt,
      contract_signature = contract_signature,
      # Two fields, two frozen vocabularies. `$result_capability` is a
      # capability id ("complete_form" / "query_only") and is hashed into
      # `$contract_signature` via `.effect_result_signature()`.
      # `$completeness` is the human-facing word the print methods show, and
      # its counterpart on a measurement form reads "complete" rather than
      # "full". The words cannot be reconciled without a hash change:
      # `.measurement_contract_signature()` (R/measurement-result.R) and
      # `.evidence_materialization()` (R/evidence-task.R) both hash a
      # `completeness` string, so renaming any of them would invalidate every
      # recorded measurement-form and evidence-task identity. Branch on
      # `$result_capability`, which is the same word everywhere.
      result_capability = "complete_form",
      completeness = "full",
      codec = codec,
      storage = storage
    ),
    class = "effect_form"
  )
  .validate_effect_form(value)
  value
}

.effect_form_capabilities <- function(self_form, symmetric, guaranteed_psd,
                                      coherent = FALSE) {
  flags <- c(self_form, symmetric, guaranteed_psd, coherent)
  if (!is.logical(flags) || length(flags) != 4L || anyNA(flags)) {
    .input_error("Effect-form capabilities must be logical guarantees.")
  }
  if (symmetric && !self_form) {
    .input_error("A symmetric effect form must be a self form.")
  }
  if (guaranteed_psd && !symmetric) {
    .input_error("A guaranteed-PSD effect form must be symmetric.")
  }
  list(
    self_form = self_form,
    symmetric = symmetric,
    guaranteed_psd = guaranteed_psd,
    total = TRUE,
    coherent = coherent,
    configuration = coherent
  )
}

.validate_geometry_store <- function(store, label, codec, probe = TRUE) {
  if (!inherits(store, "effect_geometry_store") ||
      !.is_finite_numeric(store$dim) || length(store$dim) != 2L ||
      anyNA(store$dim) || any(store$dim < 1L) || any(store$dim %% 1 != 0) ||
      !is.function(store$read)) {
    .input_error(sprintf("`%s` is not a valid complete geometry store.", label))
  }
  manifest <- store$manifest
  if (!is.list(manifest) || !identical(manifest$schema_version, 1L) ||
      !isTRUE(manifest$complete) || !identical(manifest$dim, store$dim) ||
      !identical(manifest$format, .effect_form_codec_format(codec))) {
    .input_error(
      sprintf("`%s` has an incomplete or inconsistent store manifest.", label)
    )
  }
  if (isTRUE(probe)) {
    probes <- unique(c(1L, store$dim[[1L]]))
    tryCatch(
      .read_geometry_store(store, probes),
      error = function(error) .input_error(sprintf("`%s` reader cannot supply claimed geometry: %s",
        label, conditionMessage(error)))
    )
  }
  invisible(store)
}

.validate_effect_form <- function(x, probe = TRUE) {
  expected <- c(
    "total", "coherent", "marginals", "left_space", "right_space",
    "logical_shape", "capabilities", "effect_space", "effects", "index",
    "metadata", "receipt", "contract_signature", "result_capability",
    "completeness", "codec", "storage"
  )
  if (!.sealed_fields(x, "effect_form", expected) ||
      !identical(x$result_capability, "complete_form") ||
      !identical(x$completeness, "full")) {
    .input_error("`x` must be a canonical complete effect_form.")
  }
  left_space <- .validate_effect_space(x$left_space)
  right_space <- .validate_effect_space(x$right_space)
  self_form <- .same_effect_space(left_space, right_space)
  if (!identical(x$logical_shape, as.integer(c(
      length(left_space$coordinates), length(right_space$coordinates)
    )))) {
    .contract_error("Effect-form logical shape is inconsistent with its axes.")
  }
  expected_capabilities <- .effect_form_capabilities(
    self_form,
    isTRUE(x$capabilities$symmetric),
    isTRUE(x$capabilities$guaranteed_psd),
    !is.null(x$coherent)
  )
  if (!identical(x$capabilities, expected_capabilities)) {
    .contract_error(
      "Effect-form capabilities are missing, forged, or inconsistent."
    )
  }
  codec <- match.arg(x$codec, c("rectangular", "symmetric_packed"))
  if (codec == "symmetric_packed" && !isTRUE(x$capabilities$symmetric)) {
    .input_error("The symmetric-packed codec requires a symmetric capability.")
  }
  .validate_geometry_store(x$total, "total", codec, probe = probe)
  if (!is.null(x$coherent)) {
    .validate_geometry_store(x$coherent, "coherent", codec, probe = probe)
    if (!identical(x$total$dim, x$coherent$dim)) {
      .input_error("`total` and `coherent` must have identical dimensions.")
    }
  }
  expected_width <- if (codec == "rectangular") {
    prod(x$logical_shape)
  } else {
    q <- x$logical_shape[[1L]]
    if (x$logical_shape[[2L]] != q) NA_real_ else q * (q + 1L) / 2L
  }
  if (!is.finite(expected_width) || x$total$dim[[2L]] != expected_width) {
    message <- if (codec == "symmetric_packed") {
      "Packed geometry width is not triangular or does not match its effect space."
    } else {
      "Rectangular form width does not match its left-by-right logical shape."
    }
    .input_error(message)
  }
  if (length(x$index) != x$total$dim[[1L]] || anyNA(x$index) ||
      anyDuplicated(x$index)) {
    .input_error("`index` must uniquely identify every measurement row.")
  }
  if (!is.list(x$metadata)) .input_error("`metadata` must be a list.")
  .validate_execution_receipt(x$receipt)
  if (!identical(x$metadata$scientific_plan_id,
      x$receipt$scientific_plan_id)) {
    .contract_error(paste0(
      "Result metadata and execution receipt identify different scientific ",
      "plans."
    ))
  }
  expected_signature <- .effect_result_signature(
    "complete_form", left_space, right_space, x$logical_shape,
    capabilities = x$capabilities, codec = codec,
    scientific_plan_id = x$receipt$scientific_plan_id
  )
  .check_signature(
    x$contract_signature, expected_signature,
    "Effect-form contract signature is inconsistent with its claims."
  )
  if (self_form) {
    if (!identical(x$effect_space, left_space) ||
        !identical(x$effects, left_space$coordinates)) {
      .contract_error(
        "Geometry coordinate labels are inconsistent with its effect space."
      )
    }
  } else if (!is.null(x$effect_space) || !is.null(x$effects)) {
    .input_error(
      "Rectangular forms cannot claim one compatibility effect space."
    )
  }
  if (!is.null(x$marginals) && !is.list(x$marginals)) {
    .input_error("Effect-form marginals must be NULL or a list.")
  }
  expected_storage <- unique(c(x$total$representation,
    if (is.null(x$coherent)) character() else x$coherent$representation))
  if (!identical(x$storage, expected_storage)) {
    .contract_error(
      "Geometry storage metadata is inconsistent with its component stores."
    )
  }
  invisible(x)
}

#' Construct a semantically complete effect geometry
#'
#' An `effect_geometry` is the symmetric packed self-form specialization of a
#' complete `effect_form`. Its storage may be in memory or block-backed, while
#' configuration remains exactly total minus coherent.
#'
#' @param total,coherent Packed geometry matrices, with measurements in rows,
#'   or internal geometry stores having the same dimensions.
#' @param marginals Pairing-appropriate signed marginals. Undirected pairings
#'   contain `endpoint`; directed pairings contain `left` and `right`.
#' @param effects An `effect_space()` whose dimension must match the triangular
#'   packed-geometry width and marginal columns. Unique names are accepted as
#'   shorthand for an unspecified-basis space.
#' @param index Optional measurement index with one entry per geometry row.
#' @param receipt The `execution_receipt()` proving how the result was made.
#' @param metadata Optional compact semantic metadata.
#' @return A complete `effect_geometry` and `effect_form`.
#' @keywords internal
effect_geometry <- function(total, coherent, marginals, effects, receipt, index = NULL,
                            metadata = list()) {
  effects <- .as_effect_space(effects)
  value <- effect_form(
    total = total,
    coherent = coherent,
    marginals = marginals,
    left_space = effects,
    right_space = effects,
    receipt = receipt,
    index = index,
    metadata = metadata,
    codec = "symmetric_packed",
    symmetric = TRUE,
    guaranteed_psd = FALSE
  )
  class(value) <- c("effect_geometry", "effect_form")
  .validate_effect_geometry(value)
  value
}

.validate_effect_geometry <- function(x, probe = TRUE) {
  .check_class(
    x, "effect_geometry", "x", what = "a canonical complete effect_geometry"
  )
  .validate_effect_form(x, probe = probe)
  if (!identical(x$codec, "symmetric_packed") ||
      !isTRUE(x$capabilities$self_form) ||
      !isTRUE(x$capabilities$symmetric) || is.null(x$coherent)) {
    .input_error(paste0(
      "An effect_geometry must be a symmetric packed self form with coherent ",
      "data."
    ))
  }
  marginals <- x$marginals
  if (!inherits(marginals, "effect_marginals") || length(marginals) < 1L) {
    .input_error(
      "`marginals` must be a nonempty pairing-appropriate marginal object."
    )
  }
  semantics <- attr(marginals, "semantics")
  expected_names <- switch(semantics,
    undirected_endpoint = "endpoint",
    directed_roles = c("left", "right"),
    .input_error("Marginal semantics are missing or invalid.")
  )
  if (!identical(names(marginals), expected_names)) {
    .contract_error(
      "Marginal members do not match their declared endpoint semantics."
    )
  }
  if (!all(vapply(marginals, function(value) {
    is.matrix(value) && is.numeric(value) && identical(dim(value),
      c(x$total$dim[[1L]], length(x$effects))) && all(is.finite(value)) &&
      identical(colnames(value), x$effects)
  }, logical(1)))) {
    .input_error(
      "Every marginal must match measurement rows and named effect columns."
    )
  }
  invisible(x)
}

#' Construct an explicit query-only effect view
#'
#' An `effect_view` contains derived values but not the geometry from which they
#' arose. It is never returned where a complete `effect_geometry` is promised.
#'
#' @param values A finite numeric measurement-by-view matrix.
#' @param query The compiled query matrix or a descriptive query value.
#' @param component Geometry component used to create the view.
#' @param index Optional measurement index.
#' @param receipt The `execution_receipt()` proving how the view was made.
#' @param metadata Optional compact semantic metadata.
#' @param effects Compatibility self-space binding. Prefer `left_space` and
#'   `right_space` for rectangular views.
#' @param left_space,right_space Optional ordered effect-space bindings.
#' @return A query-only `effect_view`.
#' @keywords internal
effect_view <- function(values, query, component, receipt, index = NULL,
                        metadata = list(), effects = NULL,
                        left_space = NULL, right_space = NULL) {
  .check_matrix(values, "values", what = "a finite numeric matrix")
  component <- match.arg(component, c("total", "coherent", "configuration"))
  if (is.null(index)) index <- seq_len(nrow(values))
  if (length(index) != nrow(values)) {
    .input_error("`index` must have one entry per measurement.")
  }
  if (!is.list(metadata)) {
    .input_error("`metadata` must be a list.")
  }
  .validate_execution_receipt(receipt)
  if (!is.null(effects)) {
    effects <- .validate_effect_space(effects)
    if (is.null(left_space)) left_space <- effects
    if (is.null(right_space)) right_space <- effects
  }
  if (inherits(query, "effect_pair_query")) {
    .validate_query_for_compile(query)
    if (is.null(left_space)) left_space <- query$left_space
    if (is.null(right_space)) right_space <- query$right_space
    if (!.same_effect_space(left_space, query$left_space) ||
        !.same_effect_space(right_space, query$right_space)) {
      .contract_error("Effect-view axes are incompatible with its pair query.")
    }
  } else if (inherits(query, "effect_query") &&
      identical(query$kind, "bilinear") && !is.null(query$effect_space)) {
    if (is.null(left_space)) left_space <- query$effect_space
    if (is.null(right_space)) right_space <- query$effect_space
    if (!.same_effect_space(left_space, query$effect_space) ||
        !.same_effect_space(right_space, query$effect_space)) {
      .contract_error(
        "Effect-view axes are incompatible with its bilinear query."
      )
    }
  }
  if (is.null(left_space) || is.null(right_space)) {
    packed_width <- if (inherits(query, "effect_query")) {
      nrow(query$operator) * (nrow(query$operator) + 1L) / 2L
    } else if (is.matrix(query)) {
      nrow(query)
    } else {
      NA_real_
    }
    effect_count <- (sqrt(8 * packed_width + 1) - 1) / 2
    if (!is.finite(effect_count) || effect_count < 1L || effect_count %% 1 != 0) {
      .input_error(paste0(
        "`effects` is required when the query does not identify a packed ",
        "effect space."
      ))
    }
    inferred <- effect_space(paste0("effect", seq_len(effect_count)))
    left_space <- inferred
    right_space <- inferred
  }
  left_space <- .validate_effect_space(left_space)
  right_space <- .validate_effect_space(right_space)
  self_form <- .same_effect_space(left_space, right_space)
  if (!is.null(metadata$scientific_plan_id) &&
      !identical(metadata$scientific_plan_id, receipt$scientific_plan_id)) {
    .contract_error(
      "Effect-view metadata and receipt identify different scientific plans."
    )
  }
  metadata$scientific_plan_id <- receipt$scientific_plan_id
  logical_shape <- as.integer(c(
    length(left_space$coordinates), length(right_space$coordinates)
  ))
  contract_signature <- .effect_result_signature(
    "query_only", left_space, right_space, logical_shape,
    component = component, query = query,
    scientific_plan_id = receipt$scientific_plan_id
  )

  value <- structure(
    list(
      values = values,
      query = query,
      left_space = left_space,
      right_space = right_space,
      logical_shape = logical_shape,
      effect_space = if (self_form) left_space else NULL,
      component = component,
      index = index,
      metadata = metadata,
      receipt = receipt,
      contract_signature = contract_signature,
      result_capability = "query_only",
      completeness = "query_only"
    ),
    class = "effect_view"
  )
  .validate_effect_view(value)
  value
}

.validate_effect_view <- function(x) {
  expected <- c(
    "values", "query", "left_space", "right_space", "logical_shape",
    "effect_space", "component", "index", "metadata", "receipt",
    "contract_signature", "result_capability", "completeness"
  )
  if (!inherits(x, "effect_view") || inherits(x, "effect_form") ||
      !identical(names(x), expected) ||
      !identical(x$result_capability, "query_only") ||
      !identical(x$completeness, "query_only")) {
    .input_error("`x` must be a canonical query-only effect_view.")
  }
  left_space <- .validate_effect_space(x$left_space)
  right_space <- .validate_effect_space(x$right_space)
  if (!identical(x$logical_shape, as.integer(c(
      length(left_space$coordinates), length(right_space$coordinates)
    )))) {
    .contract_error("Effect-view logical shape is inconsistent with its axes.")
  }
  self_form <- .same_effect_space(left_space, right_space)
  if ((self_form && !identical(x$effect_space, left_space)) ||
      (!self_form && !is.null(x$effect_space))) {
    .contract_error(
      "Effect-view compatibility space is inconsistent with its axes."
    )
  }
  if (!.is_finite_matrix(x$values) || length(x$index) != nrow(x$values) ||
      anyNA(x$index)) {
    .input_error("Effect-view values or measurement index are invalid.")
  }
  .validate_execution_receipt(x$receipt)
  if (!is.list(x$metadata) || !identical(x$metadata$scientific_plan_id,
      x$receipt$scientific_plan_id)) {
    .contract_error(
      "Effect-view metadata and receipt identify different scientific plans."
    )
  }
  expected_signature <- .effect_result_signature(
    "query_only", left_space, right_space, x$logical_shape,
    component = x$component, query = x$query,
    scientific_plan_id = x$receipt$scientific_plan_id
  )
  .check_signature(
    x$contract_signature, expected_signature,
    "Effect-view contract signature is inconsistent with its claims."
  )
  invisible(x)
}

#' Read one component of a complete geometry
#'
#' `geometry_component()` returns the packed geometry rows themselves, for
#' the cases where a view is not enough. Reach for it only after
#' [materialize_geometry()]; a query-only result from [evaluate_geometry()]
#' has no stored component to read.
#'
#' @param x A complete `effect_form` (including an `effect_geometry`).
#' @param component One of `total`, `coherent`, or `configuration`.
#'   `configuration` is computed exactly as `total - coherent`.
#' @param rows Optional measurement rows to read. Block-backed stores read
#'   only the requested rows.
#' @return A numeric matrix with one row per measurement and one column per
#'   packed geometry coordinate (`svec` order for symmetric self forms: the
#'   lower triangle by column, off-diagonal entries scaled by `sqrt(2)`).
#' @seealso [query_geometry()] to apply a linear query instead of reading
#'   packed coordinates, and [rdm()] or [contrast_energy()] for the named
#'   scientific views.
#' @family geometry plans and views
#' @examples
#' domain <- abstract_domain(4, id = "component-example")
#' relation <- relation(
#'   list(run1 = rbind(a = c(1, 0, 2, 1), b = c(0, 1, 1, 0)),
#'        run2 = rbind(a = c(1.1, 0.1, 1.9, 0.8), b = c(0.1, 0.9, 1.2, 0.2))),
#'   domain = domain
#' )
#' geometry <- materialize_geometry(plan_geometry(
#'   relation, compile_frame(regions(c("v1", "v1", "it", "it")), domain),
#'   cross_partitions(relation, independence = "independent")
#' ))
#'
#' # One row per region; three packed coordinates for a 2-effect self form,
#' # namely G[1,1], G[2,1], and G[2,2].
#' geometry_component(geometry, "total")
#'
#' # The coherent/configuration split is an exact partition, not a fit.
#' all.equal(
#'   geometry_component(geometry, "configuration"),
#'   geometry_component(geometry, "total") -
#'     geometry_component(geometry, "coherent")
#' )
#'
#' # Read a single measurement without touching the rest.
#' geometry_component(geometry, "total", rows = 2)
#' @export
geometry_component <- function(x, component = "total", rows = NULL) {
  if (!inherits(x, "effect_form") || !is.list(x) ||
      !inherits(x$total, "effect_geometry_store")) {
    .input_error(paste0(
      "`x` must be a complete effect_geometry or effect_form, not a ",
      "query-only view."
    ))
  }
  component <- match.arg(component, c("total", "coherent", "configuration"))
  if (!is.null(rows) &&
      (!is.numeric(rows) || anyNA(rows) || any(rows %% 1 != 0) ||
       any(rows < 1) || any(rows > x$total$dim[[1L]]))) {
    .input_error("`rows` contains invalid measurement indices.")
  }
  .validate_effect_form(x)
  .require_effect_form_component(x, component)
  .geometry_component_validated(x, component, rows)
}

.require_effect_form_component <- function(x, component) {
  if (!isTRUE(x$capabilities[[component]])) {
    .input_error(
      sprintf("This effect form does not carry the `%s` component.", component)
    )
  }
  invisible(x)
}

.geometry_component_validated <- function(x, component, rows = NULL) {
  total <- if (component != "coherent") .read_geometry_store(x$total, rows) else NULL
  coherent <- if (component != "total") .read_geometry_store(x$coherent, rows) else NULL

  switch(component,
    total = total,
    coherent = coherent,
    configuration = total - coherent
  )
}

#' Apply a linear query to a complete geometry
#'
#' `query_geometry()` projects an already materialized geometry through a
#' fixed linear query. It answers the same question as
#' [evaluate_geometry()] and carries the same view identity; use this form
#' when the complete geometry already exists and several queries will be read
#' from it.
#'
#' @param x A complete `effect_form` (including an `effect_geometry`).
#' @param query An axis-bound `pair_query()`, a compatible
#'   `bilinear_query()`, or a finite physical-coordinate-by-view matrix.
#' @param component Geometry component to query.
#' @param row_block Positive number of measurement rows read at once. This
#'   bounds packed-geometry memory for block-backed stores.
#' @return An `effect_view`: `$values` has one row per measurement and one
#'   column per query column, alongside `$query`, `$component`, `$index`, and
#'   a `$receipt` recording that this view was projected from the parent
#'   estimand. It is a view, not another geometry.
#' @seealso [evaluate_geometry()] for the query-first route that never
#'   materializes geometry, and [geometry_component()] for the packed rows.
#' @family geometry plans and views
#' @examples
#' domain <- abstract_domain(4, id = "query-example")
#' relation <- relation(
#'   list(run1 = rbind(a = c(1, 0, 2, 1), b = c(0, 1, 1, 0)),
#'        run2 = rbind(a = c(1.1, 0.1, 1.9, 0.8), b = c(0.1, 0.9, 1.2, 0.2))),
#'   domain = domain
#' )
#' geometry <- materialize_geometry(plan_geometry(
#'   relation, compile_frame(regions(c("v1", "v1", "it", "it")), domain),
#'   cross_partitions(relation, independence = "independent")
#' ))
#'
#' # The squared cross-generalized distance between the two effects, one row
#' # per region.
#' contrast <- bilinear_query(tcrossprod(c(1, -1)))
#' distance <- query_geometry(geometry, contrast)
#' distance
#' as.data.frame(distance)
#'
#' # The same query reads the two orthogonal modes, which sum back exactly.
#' coherent <- query_geometry(geometry, contrast, component = "coherent")
#' configuration <- query_geometry(
#'   geometry, contrast, component = "configuration"
#' )
#' all.equal(distance$values, coherent$values + configuration$values)
#' @export
query_geometry <- function(x, query, component = "total", row_block = 1024L) {
  validated <- .validate_geometry_query(x, query, component, row_block)
  query <- validated$query
  semantic_query <- validated$semantic_query
  component <- validated$component
  row_block <- validated$row_block

  structured <- .is_pair_difference_query(query)
  output_width <- if (structured) .query_output_width(query) else ncol(query)
  values <- matrix(0, x$total$dim[[1L]], output_width)
  colnames(values) <- .query_output_labels(query)
  packed_positions <- if (structured && identical(x$codec, "symmetric_packed")) {
    # Packed svec layout: columns of the lower triangle in order, diagonal
    # entries unscaled and off-diagonal entries carrying sqrt(2).
    q <- x$logical_shape[[1L]]
    offset <- function(column) (column - 1) * (2 * q - column + 2) / 2
    list(
      ii = offset(query$pair_left) + 1,
      jj = offset(query$pair_right) + 1,
      ij = offset(query$pair_left) + (query$pair_right - query$pair_left + 1),
      cross_scale = sqrt(2)
    )
  } else if (structured) {
    # Rectangular column-major layout of a symmetric self form.
    q <- x$logical_shape[[1L]]
    list(
      ii = (query$pair_left - 1) * q + query$pair_left,
      jj = (query$pair_right - 1) * q + query$pair_right,
      ij = (query$pair_left - 1) * q + query$pair_right,
      cross_scale = 2
    )
  } else {
    NULL
  }
  for (start in .tile_starts(x$total$dim[[1L]], row_block)) {
    rows <- start:min(start + row_block - 1L, x$total$dim[[1L]])
    packed <- .geometry_component_validated(x, component, rows)
    values[rows, ] <- if (structured) {
      distances <- packed[, packed_positions$ii, drop = FALSE] +
        packed[, packed_positions$jj, drop = FALSE] -
        packed_positions$cross_scale *
          packed[, packed_positions$ij, drop = FALSE]
      if (is.null(query$coefficients)) {
        distances
      } else {
        distances %*% t(query$coefficients)
      }
    } else {
      packed %*% query
    }
  }
  # Projection executes a view of the parent estimand; its identity must be
  # route-stable with the fused query-first execution of the same view.
  view_receipt <- .projection_receipt(
    x$receipt,
    .geometry_view_scientific_id(
      x$receipt$scientific_plan_id, component, query
    )
  )
  metadata <- x$metadata
  metadata$scientific_plan_id <- NULL
  metadata$projected_from <- x$receipt$scientific_plan_id
  effect_view(
    values = values,
    query = if (is.null(semantic_query)) query else semantic_query,
    component = component,
    receipt = view_receipt,
    index = x$index,
    metadata = metadata,
    left_space = x$left_space,
    right_space = x$right_space
  )
}

.validate_geometry_query <- function(x, query, component, row_block) {
  if (!inherits(x, "effect_form") || !identical(x$completeness, "full") ||
      !inherits(x$total, "effect_geometry_store") ||
      !.is_finite_numeric(x$total$dim) || length(x$total$dim) != 2L ||
      any(x$total$dim < 1L)) {
    .input_error("`x` must be a structurally complete effect_form.")
  }
  .validate_effect_form(x, probe = FALSE)
  component <- match.arg(component, c("total", "coherent", "configuration"))
  row_block <- .validate_tile_size(row_block, "row_block")
  .require_effect_form_component(x, component)

  semantic_query <- NULL
  if (.is_pair_difference_query(query)) {
    if (!isTRUE(x$capabilities$self_form) ||
        !isTRUE(x$capabilities$symmetric)) {
      .input_error("A pair-difference query requires a symmetric self form.")
    }
    if (!identical(query$effects, x$effects)) {
      .contract_error("The query and geometry effect spaces are incompatible.")
    }
    .validate_effect_form(x)
    return(list(query = query, semantic_query = query,
      component = component, row_block = row_block))
  }
  if (inherits(query, "effect_query")) {
    .validate_query_for_compile(query)
    semantic_query <- query
    if (identical(query$kind, "pair")) {
      if (!.same_effect_space(query$left_space, x$left_space) ||
          !.same_effect_space(query$right_space, x$right_space)) {
        .contract_error("The query and form axis identities are incompatible.")
      }
    } else if (identical(query$kind, "bilinear") && isTRUE(query$fixed)) {
      if (!isTRUE(x$capabilities$self_form) ||
          !isTRUE(x$capabilities$symmetric)) {
        .input_error("A bilinear_query requires a symmetric self form.")
      }
      if (nrow(query$operator) != x$logical_shape[[1L]]) {
        .contract_error(
          "The query operator dimension must equal the experimental dimension."
        )
      }
      if (!is.null(query$effect_space) &&
          !.same_effect_space(query$effect_space, x$left_space)) {
        .contract_error(
          "The query and geometry effect spaces are incompatible."
        )
      }
    } else {
      .input_error(
        "Form operator queries must be fixed pair or bilinear queries."
      )
    }
    query <- .physical_form_query(query$operator, x)
  }
  if (!.is_finite_matrix(query) || nrow(query) < 1L || ncol(query) < 1L) {
    .input_error(
      "`query` must be a finite, nonempty physical-coordinate matrix."
    )
  }
  if (nrow(query) != x$total$dim[[2L]]) {
    .input_error("The query input dimension must equal the stored form width.")
  }
  if (is.null(colnames(query))) colnames(query) <- paste0("view", seq_len(ncol(query)))
  .validate_effect_form(x)
  list(query = query, semantic_query = semantic_query,
    component = component, row_block = row_block)
}

.physical_form_query <- function(operator, x) {
  values <- if (identical(x$codec, "rectangular")) {
    as.vector(operator)
  } else {
    .svec_symmetric(0.5 * (operator + t(operator)))
  }
  matrix(values, ncol = 1L, dimnames = list(NULL, "view1"))
}

# The contrast-view record ---------------------------------------------------
#
# One weight vector reduces a self-form geometry to four aligned per-measurement
# series -- signed marginals, coherent, configuration, total -- plus the
# coherence fraction and the mask saying where that fraction is defined at all.
# It is a record, not a view computation: it lived in R/views.R, which forced
# the executor to call up into a view file on the query-first contrast path.
# The derivation that produces its inputs stays in `contrast_energy()`.
.new_effect_contrast_view <- function(total, coherent, marginals, weights,
                                      index, receipt) {
  total <- drop(total)
  coherent <- drop(coherent)
  configuration <- total - coherent
  signed <- lapply(marginals, function(value) drop(value %*% weights))
  if (identical(attr(marginals, "semantics"), "undirected_endpoint")) {
    signed <- signed$endpoint
  }
  fraction_valid <- is.finite(total) & total > 0 & coherent >= 0 &
    configuration >= 0
  coherence_fraction <- rep(NA_real_, length(total))
  coherence_fraction[fraction_valid] <- coherent[fraction_valid] /
    total[fraction_valid]

  structure(
    list(
      signed = signed,
      coherent = coherent,
      configuration = configuration,
      total = total,
      coherence_fraction = coherence_fraction,
      coherence_fraction_valid = fraction_valid,
      weights = weights,
      index = index,
      receipt = receipt
    ),
    class = "effect_contrast_view"
  )
}
