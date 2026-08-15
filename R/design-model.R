# Declared first-moment design models --------------------------------------

.normalize_partition_matrices <- function(value, name) {
  if (is.matrix(value)) value <- list(value)
  if (!is.list(value) || length(value) < 1L) {
    stop(sprintf("`%s` must provide at least one matrix.", name), call. = FALSE)
  }
  partitions <- names(value)
  if (is.null(partitions) || any(!nzchar(partitions))) {
    stop(sprintf("Every `%s` matrix must have a partition name.", name),
      call. = FALSE)
  }
  partitions <- .validate_partition_names(partitions, length(value), name)
  names(value) <- partitions
  for (partition in partitions) {
    matrix <- value[[partition]]
    if (!is.matrix(matrix) || !is.numeric(matrix) || any(dim(matrix) < 1L) ||
        any(!is.finite(matrix))) {
      stop(sprintf("Partition `%s` `%s` is not a finite nonempty matrix.",
        partition, name), call. = FALSE)
    }
  }
  value
}

.normalize_design_rows <- function(designs, row_ids) {
  partitions <- names(designs)
  if (is.null(row_ids)) {
    row_ids <- lapply(designs, rownames)
  } else if (!is.list(row_ids) && length(partitions) == 1L) {
    row_ids <- list(row_ids)
  }
  if (!is.list(row_ids) || length(row_ids) != length(partitions)) {
    stop("`row_ids` must provide one ordered observation axis per design.",
      call. = FALSE)
  }
  names(row_ids) <- partitions
  for (partition in partitions) {
    ids <- row_ids[[partition]]
    if (is.null(ids) || length(ids) != nrow(designs[[partition]]) ||
        anyNA(ids) || anyDuplicated(ids)) {
      stop(sprintf(
        "Partition `%s` requires unique row identifiers for every design row.",
        partition
      ), call. = FALSE)
    }
    if (!is.null(rownames(designs[[partition]])) &&
        !identical(rownames(designs[[partition]]), as.character(ids))) {
      stop(sprintf(
        "Partition `%s` design row names disagree with `row_ids`.", partition
      ), call. = FALSE)
    }
    rownames(designs[[partition]]) <- as.character(ids)
    row_ids[[partition]] <- unname(ids)
  }
  list(designs = designs, row_ids = row_ids)
}

.normalize_solver_routes <- function(value, partitions) {
  if (!is.character(value) || anyNA(value) ||
      !length(value) %in% c(1L, length(partitions)) ||
      any(!value %in% c("auto", "qr", "svd"))) {
    stop("`solver` must provide `auto`, `qr`, or `svd` per partition.",
      call. = FALSE)
  }
  if (length(value) == 1L) value <- rep(value, length(partitions))
  stats::setNames(unname(value), partitions)
}

.design_compiler_record <- function(protocol, protocol_version, package,
                                    package_version) {
  fields <- list(
    protocol = protocol,
    protocol_version = protocol_version,
    package = package,
    package_version = package_version
  )
  fields <- lapply(names(fields), function(name) {
    .validate_nonempty_id(fields[[name]], name)
  })
  names(fields) <- c("protocol", "protocol_version", "package", "package_version")
  fields
}

#' Declare a semantic design model with compiled routes
#'
#' `design_model_id` covers the semantic mean-model request. Concrete design
#' matrices, coding maps, row order, compiler build, and solver route are kept
#' in `compilation_route_id` and later design receipts.
#'
#' @param specification A portable semantic model declaration.
#' @param conditions The bound [condition_space()].
#' @param designs Named observation-by-coefficient matrices, one per partition.
#' @param parameterizations Named [coefficient_parameterization()] values for
#'   the designs.
#' @param row_ids Ordered observation identifiers per design. Defaults to
#'   design row names; positional designs are refused.
#' @param solver Numerical route per partition.
#' @param protocol,protocol_version Semantic compiler protocol identity.
#' @param package,package_version Compiler implementation receipt fields.
#' @param provenance Portable semantic model provenance.
#' @return An `effect_design_model`.
#' @export
design_model <- function(
    specification, conditions, designs, parameterizations, row_ids = NULL,
    solver = "auto", protocol = "effectagram-semantic-design",
    protocol_version = "1", package = "effectagram",
    package_version = "0.0.0.9000", provenance = list()) {
  if (!is.list(specification)) {
    stop("`specification` must be a portable semantic declaration list.",
      call. = FALSE)
  }
  .validate_effect_provenance(specification, "design specification")
  .validate_effect_provenance(provenance, "design-model provenance")
  conditions <- .validate_condition_space(conditions)
  designs <- .normalize_partition_matrices(designs, "designs")
  partitions <- names(designs)
  rows <- .normalize_design_rows(designs, row_ids)
  designs <- rows$designs
  row_ids <- rows$row_ids

  if (inherits(parameterizations, "effect_coefficient_parameterization")) {
    parameterizations <- rep(list(parameterizations), length(partitions))
  }
  if (!is.list(parameterizations) ||
      length(parameterizations) != length(partitions)) {
    stop("`parameterizations` must provide one coding per design partition.",
      call. = FALSE)
  }
  names(parameterizations) <- partitions
  parameterizations <- lapply(parameterizations,
    .validate_coefficient_parameterization)
  for (partition in partitions) {
    parameterization <- parameterizations[[partition]]
    if (!identical(parameterization$condition_space$signature,
        conditions$signature)) {
      stop(sprintf("Partition `%s` uses a different condition space.",
        partition), call. = FALSE)
    }
    coefficients <- colnames(designs[[partition]])
    if (is.null(coefficients) ||
        !identical(coefficients, parameterization$coefficients)) {
      stop(sprintf(
        "Partition `%s` design columns must equal its coefficient axis.",
        partition
      ), call. = FALSE)
    }
  }
  solver <- .normalize_solver_routes(solver, partitions)
  compiler <- .design_compiler_record(
    protocol, protocol_version, package, package_version
  )
  semantic <- list(
    schema_version = 1L,
    protocol = compiler[c("protocol", "protocol_version")],
    condition_space = conditions,
    specification = specification,
    provenance = provenance
  )
  design_model_id <- .semantic_digest("design-model-sha256:", semantic)
  route <- list(
    schema_version = 1L,
    design_model_id = design_model_id,
    partitions = partitions,
    designs = designs,
    parameterization_ids = vapply(
      parameterizations, `[[`, character(1), "parameterization_id"
    ),
    row_ids = row_ids,
    compiler = compiler,
    solver = solver
  )
  structure(list(
    condition_space = conditions,
    specification = specification,
    partitions = partitions,
    designs = designs,
    parameterizations = parameterizations,
    row_ids = row_ids,
    compiler = compiler,
    solver = solver,
    provenance = provenance,
    design_model_id = design_model_id,
    compilation_route_id = .semantic_digest(
      "design-compilation-route-sha256:", route
    ),
    capabilities = list(
      symbolic_model = TRUE,
      coding_invariant = TRUE,
      row_lineage = TRUE
    )
  ), class = "effect_design_model")
}

#' Bind already compiled raw design matrices
#'
#' This degenerate route binds design values, coefficient order, row identity,
#' and solver into model identity. It makes no semantic coding-invariance claim.
#'
#' @param designs Named observation-by-coefficient design matrices.
#' @param row_ids Ordered observation identifiers per design.
#' @param solver Numerical route per partition.
#' @param provenance Portable provenance for the external construction.
#' @return An `effect_raw_design_model`, also an `effect_design_model`.
#' @export
raw_design_model <- function(designs, row_ids = NULL, solver = "auto",
                             provenance = list()) {
  .validate_effect_provenance(provenance, "raw-design provenance")
  designs <- .normalize_partition_matrices(designs, "designs")
  partitions <- names(designs)
  rows <- .normalize_design_rows(designs, row_ids)
  designs <- rows$designs
  row_ids <- rows$row_ids
  for (partition in partitions) {
    if (is.null(colnames(designs[[partition]]))) {
      stop(sprintf("Partition `%s` raw design requires coefficient names.",
        partition), call. = FALSE)
    }
    .validate_effect_names(
      colnames(designs[[partition]]), ncol(designs[[partition]])
    )
  }
  solver <- .normalize_solver_routes(solver, partitions)
  compiler <- .design_compiler_record(
    "effectagram-raw-design", "1", "effectagram", "0.0.0.9000"
  )
  semantic <- list(
    schema_version = 1L,
    kind = "raw-X",
    partitions = partitions,
    designs = designs,
    row_ids = row_ids,
    provenance = provenance
  )
  design_model_id <- .semantic_digest("raw-design-model-sha256:", semantic)
  route <- list(
    schema_version = 1L,
    design_model_id = design_model_id,
    compiler = compiler,
    solver = solver
  )
  structure(list(
    condition_space = NULL,
    specification = list(kind = "raw-X"),
    partitions = partitions,
    designs = designs,
    parameterizations = NULL,
    row_ids = row_ids,
    compiler = compiler,
    solver = solver,
    provenance = provenance,
    design_model_id = design_model_id,
    compilation_route_id = .semantic_digest(
      "design-compilation-route-sha256:", route
    ),
    capabilities = list(
      symbolic_model = FALSE,
      coding_invariant = FALSE,
      row_lineage = TRUE
    )
  ), class = c("effect_raw_design_model", "effect_design_model"))
}

.validate_design_model <- function(value) {
  expected <- c(
    "condition_space", "specification", "partitions", "designs",
    "parameterizations", "row_ids", "compiler", "solver", "provenance",
    "design_model_id", "compilation_route_id", "capabilities"
  )
  if (!inherits(value, "effect_design_model") || !is.list(value) ||
      !identical(names(value), expected)) {
    stop("Design-model fields are missing or noncanonical.", call. = FALSE)
  }
  raw <- inherits(value, "effect_raw_design_model")
  rebuilt <- if (raw) {
    raw_design_model(
      value$designs,
      row_ids = value$row_ids,
      solver = value$solver,
      provenance = value$provenance
    )
  } else {
    design_model(
      value$specification,
      value$condition_space,
      value$designs,
      value$parameterizations,
      row_ids = value$row_ids,
      solver = value$solver,
      protocol = value$compiler$protocol,
      protocol_version = value$compiler$protocol_version,
      package = value$compiler$package,
      package_version = value$compiler$package_version,
      provenance = value$provenance
    )
  }
  if (!identical(value, rebuilt)) {
    stop("Design-model metadata or identity is inconsistent.", call. = FALSE)
  }
  rebuilt
}
