# Declared first-moment design models --------------------------------------

.normalize_partition_matrices <- function(value, name) {
  if (is.matrix(value)) value <- list(value)
  if (!is.list(value) || length(value) < 1L) {
    .input_error(sprintf("`%s` must provide at least one matrix.", name))
  }
  partitions <- names(value)
  if (is.null(partitions) || any(!nzchar(partitions))) {
    .input_error(sprintf("Every `%s` matrix must have a partition name.", name))
  }
  partitions <- .validate_partition_names(partitions, length(value), name)
  names(value) <- partitions
  for (partition in partitions) {
    matrix <- value[[partition]]
    if (!.is_finite_matrix(matrix) || any(dim(matrix) < 1L)) {
      .input_error(sprintf(
        "Partition `%s` `%s` is not a finite nonempty matrix.",
        partition, name))
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
    .input_error(
      "`row_ids` must provide one ordered observation axis per design."
    )
  }
  names(row_ids) <- partitions
  for (partition in partitions) {
    ids <- row_ids[[partition]]
    if (is.null(ids) || length(ids) != nrow(designs[[partition]]) ||
        anyNA(ids) || anyDuplicated(ids)) {
      .input_error(sprintf(
        "Partition `%s` requires unique row identifiers for every design row.",
        partition
      ))
    }
    if (!is.null(rownames(designs[[partition]])) &&
        !identical(rownames(designs[[partition]]), as.character(ids))) {
      .contract_error(sprintf(
        "Partition `%s` design row names disagree with `row_ids`.", partition
      ))
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
    .input_error("`solver` must provide `auto`, `qr`, or `svd` per partition.")
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

# The public compiler was renamed before release. The semantic protocol token
# predates that rename and already participates in plan hashes, so both spellings
# lower to the original canonical token. Compiler package and route fields keep
# the current `crossform` name and therefore still distinguish the new route.
.canonical_design_protocol <- function(protocol, protocol_version) {
  aliases <- c("effectagram-semantic-design", "crossform-semantic-design")
  if (identical(protocol_version, "1") && protocol %in% aliases) {
    return("effectagram-semantic-design")
  }
  protocol
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
#' @return An `effect_design_model`: a list with the bound `$condition_space`,
#'   the `$specification`, `$partitions`, the compiled `$designs`, their
#'   `$parameterizations`, `$row_ids`, the `$compiler` record, the `$solver`
#'   route, `$provenance`, a `$design_model_id` covering the semantic request
#'   only, a `$compilation_route_id` covering the concrete compilation, and
#'   `$capabilities` with `symbolic_model`, `coding_invariant`, `row_lineage`.
#' @family studies and effect maps
#' @seealso [coefficient_parameterization()] for the coding it carries,
#'   [raw_design_model()] for the route without semantic coding, and
#'   [plan_relation()], which binds this model to a study.
#' @examples
#' conditions <- condition_space(c("face", "body"), basis_id = "cond-mean:v1")
#' design <- cbind(
#'   face = c(1, 0, 1, 0, 1, 0), body = c(0, 1, 0, 1, 0, 1),
#'   drift = seq(-1, 1, length.out = 6L)
#' )
#' rownames(design) <- paste0("scan-", 1:6)
#' map <- cbind(diag(2), drift = 0)
#' dimnames(map) <- list(conditions$coordinates, colnames(design))
#' coding <- coefficient_parameterization(
#'   map, conditions, coding_id = "cell-means-plus-drift"
#' )
#' specification <- list(target = "condition means", nuisance = "linear drift")
#'
#' model <- design_model(
#'   specification, conditions,
#'   designs = list(`run-1` = design), parameterizations = list(`run-1` = coding)
#' )
#' model$capabilities$coding_invariant
#'
#' # Semantic identity ignores the numerical route: switching solvers changes
#' # the compilation receipt, not what was scientifically requested.
#' route <- design_model(
#'   specification, conditions,
#'   designs = list(`run-1` = design), parameterizations = list(`run-1` = coding),
#'   solver = "svd"
#' )
#' c(same_request = identical(model$design_model_id, route$design_model_id),
#'   same_route = identical(
#'     model$compilation_route_id, route$compilation_route_id
#'   ))
#' @export
design_model <- function(
    specification, conditions, designs, parameterizations, row_ids = NULL,
    solver = "auto", protocol = "crossform-semantic-design",
    protocol_version = "1", package = "crossform",
    package_version = "0.0.0.9000", provenance = list()) {
  if (!is.list(specification)) {
    .input_error(
      "`specification` must be a portable semantic declaration list."
    )
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
    .input_error(
      "`parameterizations` must provide one coding per design partition."
    )
  }
  names(parameterizations) <- partitions
  parameterizations <- lapply(parameterizations,
    .validate_coefficient_parameterization)
  for (partition in partitions) {
    parameterization <- parameterizations[[partition]]
    if (!identical(parameterization$condition_space$signature,
        conditions$signature)) {
      .contract_error(sprintf(
        "Partition `%s` uses a different condition space.",
        partition))
    }
    coefficients <- colnames(designs[[partition]])
    if (is.null(coefficients) ||
        !identical(coefficients, parameterization$coefficients)) {
      .input_error(sprintf(
        "Partition `%s` design columns must equal its coefficient axis.",
        partition
      ))
    }
  }
  solver <- .normalize_solver_routes(solver, partitions)
  compiler <- .design_compiler_record(
    protocol, protocol_version, package, package_version
  )
  semantic <- list(
    schema_version = 1L,
    protocol = list(
      protocol = .canonical_design_protocol(
        compiler$protocol, compiler$protocol_version
      ),
      protocol_version = compiler$protocol_version
    ),
    condition_space = conditions,
    specification = specification,
    provenance = provenance
  )
  design_model_id <- .sha256_signature(semantic, "design-model-sha256:")
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
    compilation_route_id = .sha256_signature(route, "design-compilation-route-sha256:"),
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
#' @return An `effect_raw_design_model` (also an `effect_design_model`) with
#'   the same fields as [design_model()], except that `$condition_space` and
#'   `$parameterizations` are `NULL`, `$design_model_id` includes the design
#'   values themselves, and `$capabilities` reports `symbolic_model` and
#'   `coding_invariant` as `FALSE`.
#' @family studies and effect maps
#' @seealso [design_model()] for the semantic route, and [raw_effect_map()],
#'   the effect map this design must be paired with in [plan_relation()].
#' @examples
#' design <- cbind(b1 = c(1, 0, 1, 0), b2 = c(0, 1, 0, 1))
#' rownames(design) <- paste0("scan-", 1:4)
#' model <- raw_design_model(list(`run-1` = design))
#'
#' # The honest cost of skipping the condition space: the numeric design is
#' # part of the identity, so a re-coded design is a different request.
#' model$capabilities$coding_invariant
#' model$capabilities$row_lineage
#'
#' # A raw design must be paired with a raw target on the same column axis.
#' target <- rbind(`b1-b2` = c(1, -1))
#' colnames(target) <- colnames(design)
#' raw_effect_map(target)
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
      .input_error(sprintf(
        "Partition `%s` raw design requires coefficient names.",
        partition))
    }
    .validate_effect_names(
      colnames(designs[[partition]]), ncol(designs[[partition]])
    )
  }
  solver <- .normalize_solver_routes(solver, partitions)
  compiler <- .design_compiler_record(
    "crossform-raw-design", "1", "crossform", "0.0.0.9000"
  )
  semantic <- list(
    schema_version = 1L,
    kind = "raw-X",
    partitions = partitions,
    designs = designs,
    row_ids = row_ids,
    provenance = provenance
  )
  design_model_id <- .sha256_signature(semantic, "raw-design-model-sha256:")
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
    compilation_route_id = .sha256_signature(route, "design-compilation-route-sha256:"),
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
  if (!.sealed_fields(value, "effect_design_model", expected)) {
    .input_error("Design-model fields are missing or noncanonical.")
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
    .contract_error("Design-model metadata or identity is inconsistent.")
  }
  rebuilt
}
