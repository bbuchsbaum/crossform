# Observation ancestry and actual fitting support (layer 2 values) ----------
#
# Identifiers declare observational origins, never equality of numeric arrays.
# A manifest is a caller/ingestion assertion about independent sampling. It
# cannot certify hidden preprocessing or make renamed observations independent.
# This file knows no relations or plans; their collector sits above it.

.origin_ids <- function(x, where, empty = FALSE) {
  if (empty && identical(x, character())) return(x)
  if (!length(x) || !.is_strings(x, unique = TRUE)) {
    .input_error(sprintf("`%s` must contain unique nonempty observation-origin identifiers.", where))
  }
  sort(unname(x), method = "radix")
}

.origin_manifest <- function(value) {
  if (inherits(value, "effect_origin_manifest")) return(.validate_origin_manifest(value))
  if (anyDuplicated(names(value))) .input_error("Origin manifest field names cannot repeat.")
  if (!is.list(value) || !setequal(names(value),
      c("id", "partitions", "independence", "assumption", "dependencies"))) {
    # Dependencies are optional in the declaration, explicit in the value.
    if (is.list(value) && setequal(names(value),
        c("id", "partitions", "independence", "assumption"))) {
      value$dependencies <- list()
    } else {
      .input_error(paste0("An observation-origin manifest needs `id`, `partitions`, ",
        "`independence`, `assumption`, and optional `dependencies`."))
    }
  }
  id <- .validate_nonempty_id(value$id, "origin manifest id")
  assumption <- .validate_nonempty_id(value$assumption, "origin independence assumption")
  if (!identical(value$independence, "independent") &&
      !identical(value$independence, "undeclared")) {
    .input_error("Origin `independence` must be `independent` or `undeclared`.")
  }
  parts <- value$partitions
  if (!is.list(parts) || !length(parts) || !.is_strings(names(parts), unique = TRUE)) {
    .input_error("Origin `partitions` must be a named list of observation-origin sets.")
  }
  parts <- parts[sort(names(parts), method = "radix")]
  parts <- lapply(parts, .origin_ids, where = "partition origins")
  deps <- value$dependencies
  if (!is.list(deps) || (length(deps) &&
      (!.is_strings(names(deps), unique = TRUE) || !all(names(deps) %in% names(parts))))) {
    .input_error("Origin `dependencies` must name declared partitions and their upstream origins.")
  }
  deps <- stats::setNames(lapply(names(parts), function(p) {
    if (is.null(deps[[p]])) character() else .origin_ids(deps[[p]], "upstream origins", empty = TRUE)
  }), names(parts))
  origins <- sort(unique(c(unlist(parts, use.names = FALSE),
    unlist(deps, use.names = FALSE))), method = "radix")
  scope <- list(id = id, origins = origins, independence = value$independence,
    assumption = assumption)
  result <- c(scope, list(partitions = parts, dependencies = deps,
    scope_signature = .sha256_signature(scope, "origin-scope-sha256:")))
  structure(c(result, list(signature = .sha256_signature(result, "origin-manifest-sha256:"))),
    class = "effect_origin_manifest")
}

.validate_origin_manifest <- function(value) {
  fields <- c("id", "origins", "independence", "assumption", "partitions",
    "dependencies", "scope_signature", "signature")
  if (!.sealed_fields(value, "effect_origin_manifest", fields)) {
    .contract_error("Observation-origin manifest fields are missing or noncanonical.")
  }
  declaration <- unclass(value[c("id", "partitions", "independence", "assumption", "dependencies")])
  rebuilt <- .origin_manifest(declaration)
  if (!identical(value, rebuilt)) .contract_error("Observation-origin manifest identity is inconsistent.")
  value
}

.source_origin <- function(manifest, partition) {
  if (!partition %in% names(manifest$partitions)) {
    .input_error(sprintf("The origin manifest has no declaration for partition `%s`.", partition))
  }
  value <- list(scope_signature = manifest$scope_signature,
    observations = manifest$partitions[[partition]], dependencies = manifest$dependencies[[partition]])
  c(value, list(signature = .sha256_signature(value, "source-origin-sha256:")))
}

.validate_source_origin <- function(origin, manifest) {
  if (!is.list(origin) || !identical(names(origin),
      c("scope_signature", "observations", "dependencies", "signature")) ||
      !identical(origin$scope_signature, manifest$scope_signature) ||
      !identical(origin$observations, .origin_ids(origin$observations, "source origins")) ||
      !identical(origin$dependencies, .origin_ids(origin$dependencies, "source dependencies", empty = TRUE)) ||
      !all(c(origin$observations, origin$dependencies) %in% manifest$origins)) {
    .contract_error("Source observation ancestry does not match its declared manifest.")
  }
  .check_signature(origin$signature,
    .sha256_signature(origin[1:3], "source-origin-sha256:"),
    "Source observation ancestry identity is inconsistent.")
  origin
}

.training_support_fields <- c("known", "scope_signature", "assumption", "independence",
  "endpoints", "products", "pairing", "observations", "upstream", "signature")

.training_observations <- function(endpoints, upstream) {
  ids <- c(unlist(lapply(endpoints, function(x) {
    c(x$origin$observations, x$origin$dependencies)
  }), use.names = FALSE), unlist(lapply(upstream, `[[`, "observations"), use.names = FALSE))
  sort(unique(as.character(ids)), method = "radix")
}

.new_geometry_training_support <- function(manifest, endpoints, products, pairing_metadata,
                                            upstream = list()) {
  if (!is.list(upstream)) .input_error("Upstream fitting dependencies must be a list of support records.")
  upstream <- lapply(upstream, .validate_geometry_training_support)
  known <- !is.null(manifest) && all(vapply(upstream, `[[`, FALSE, "known"))
  scope <- if (is.null(manifest)) NULL else manifest$scope_signature
  if (known && any(vapply(upstream, function(x) !identical(x$scope_signature, scope), FALSE))) {
    .capability_refusal("Fitting dependencies come from different observation-origin scopes.",
      capability = "independent_geometry_prediction", namespace = "predictive_geometry",
      reasons = "origin_scope_mismatch", remedies = "Use one declared parent observation manifest.")
  }
  actual_origins <- .training_observations(endpoints, upstream)
  value <- list(known = known, scope_signature = scope,
    assumption = if (is.null(manifest)) "Observation origins are undeclared." else manifest$assumption,
    independence = if (is.null(manifest)) "undeclared" else manifest$independence,
    endpoints = endpoints, products = products, pairing = pairing_metadata,
    observations = actual_origins, upstream = upstream)
  structure(c(value, list(signature = .sha256_signature(value, "training-support-sha256:"))),
    class = "effect_geometry_training_support")
}

.validate_geometry_training_support <- function(value) {
  if (!.sealed_fields(value, "effect_geometry_training_support", .training_support_fields)) {
    .contract_error("Geometry training-support fields are missing or noncanonical.")
  }
  if (!is.list(value$endpoints) || !.is_strings(names(value$endpoints), unique = TRUE) ||
      !is.data.frame(value$products) || !identical(names(value$products), c("left", "right", "weight")) ||
      !all(c(value$products$left, value$products$right) %in% names(value$endpoints)) ||
      any(!is.finite(value$products$weight)) || any(value$products$weight <= 0) ||
      abs(sum(value$products$weight) - 1) > 1e-12 ||
      !.is_flag(value$known) || !is.list(value$upstream)) {
    .contract_error("Geometry training-support metadata are inconsistent.")
  }
  lapply(value$upstream, .validate_geometry_training_support)
  known_origins <- vapply(value$endpoints, function(x) !is.null(x$origin), FALSE)
  expected_known <- !is.null(value$scope_signature) && all(known_origins) &&
    all(vapply(value$upstream, `[[`, FALSE, "known"))
  if (!identical(value$known, expected_known)) {
    .contract_error("Training-support origin status contradicts its dependencies.")
  }
  for (endpoint in value$endpoints) {
    if (!is.list(endpoint) || !identical(names(endpoint),
        c("origin", "source_revision", "extractor_signature")) ||
        !.strong_sha256(endpoint$source_revision) ||
        !.strong_sha256(sub("^extractor-", "", endpoint$extractor_signature))) {
      .contract_error("Training-support source or extractor identities are invalid.")
    }
    if (!is.null(endpoint$origin)) {
      .validate_source_origin(endpoint$origin, list(scope_signature = value$scope_signature,
        origins = value$observations))
    }
  }
  actual_origins <- .training_observations(value$endpoints, value$upstream)
  if (!identical(value$observations, actual_origins)) {
    .contract_error("Geometry training-support origins do not match its actual dependencies.")
  }
  .check_signature(value$signature,
    .sha256_signature(unclass(value[setdiff(.training_support_fields, "signature")]),
      "training-support-sha256:"), "Geometry training-support identity is inconsistent.")
  value
}

# Both axes of independence are checked: each neural product, then every
# evaluation observation against the entire fit/selection dependency set.
.geometry_support_independence <- function(training, evaluation) {
  training <- .validate_geometry_training_support(training)
  evaluation <- .validate_geometry_training_support(evaluation)
  refuse <- function(reason, message) {
    .capability_refusal(message, capability = "independent_geometry_prediction",
      namespace = "predictive_geometry", reasons = reason,
      remedies = "Use disjoint declared observation origins for all fitting and evaluation dependencies.")
  }
  if (!training$known || !evaluation$known) {
    refuse("unknown_observation_origins", "Independent scoring requires declared observation ancestry.")
  }
  if (!identical(training$scope_signature, evaluation$scope_signature)) {
    refuse("origin_scope_mismatch", "Independent scoring requires one common declared observation-origin manifest.")
  }
  for (support in list(training, evaluation)) {
    if (!identical(support$independence, "independent") ||
        !identical(support$pairing$independence, "independent")) {
      refuse("independence_undeclared", "Both observation sampling and neural pairing independence must be declared.")
    }
    for (i in seq_len(nrow(support$products))) {
      a <- support$endpoints[[support$products$left[[i]]]]$origin
      b <- support$endpoints[[support$products$right[[i]]]]$origin
      if (length(intersect(c(a$observations, a$dependencies), c(b$observations, b$dependencies)))) {
        refuse("overlapping_neural_origins", "A declared neural cross-product shares actual observation origins.")
      }
    }
  }
  if (length(intersect(training$observations, evaluation$observations))) {
    refuse("overlapping_training_evaluation_origins", "Evaluation observations overlap fitting or selection dependencies.")
  }
  invisible(TRUE)
}
