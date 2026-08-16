# Fixed factorized measurement bridges -------------------------------------

#' Name a common measurement space
#'
#' @param n_measurements Positive common-coordinate count.
#' @param id Stable nonempty identity.
#' @param provenance Portable fixed-space provenance.
#' @return An immutable `effect_measurement_space`.
#' @export
measurement_space <- function(n_measurements, id, provenance = list()) {
  if (!is.numeric(n_measurements) || length(n_measurements) != 1L ||
      is.na(n_measurements) || !is.finite(n_measurements) ||
      n_measurements < 1L || n_measurements %% 1 != 0 ||
      n_measurements > .Machine$integer.max) {
    stop("`n_measurements` must be one positive integer.", call. = FALSE)
  }
  n_measurements <- as.integer(n_measurements)
  .domain_id(id)
  .validate_effect_provenance(provenance, "measurement-space provenance")
  semantic <- list(
    schema_version = 1L,
    id = id,
    n_measurements = n_measurements,
    provenance = provenance
  )
  structure(c(semantic[-1L], list(
    signature = paste0("sha256:", digest::digest(
      semantic, algo = "sha256", serialize = TRUE, serializeVersion = 2L
    ))
  )), class = "effect_measurement_space")
}

.validate_measurement_space <- function(x) {
  expected <- c("id", "n_measurements", "provenance", "signature")
  if (!inherits(x, "effect_measurement_space") || !is.list(x) ||
      !identical(names(x), expected)) {
    stop("Measurement-space fields are missing or noncanonical.",
      call. = FALSE)
  }
  rebuilt <- measurement_space(x$n_measurements, x$id, x$provenance)
  if (!identical(x, rebuilt)) {
    stop("Measurement-space identity or provenance is inconsistent.",
      call. = FALSE)
  }
  rebuilt
}

#' Construct a fixed factorized bridge between neural spaces
#'
#' The induced cross-space operator is `t(left_leg) %*% right_leg`. Only the
#' factorized fixed legs are public; an arbitrary dense cross operator is not.
#'
#' @param left_leg,right_leg Finite common-coordinate-by-neural-feature
#'   matrices.
#' @param left_domain,right_domain Exact neural domain values or references.
#' @param common_space A named `measurement_space()` shared by both legs.
#' @param provenance Portable provenance for the fixed legs.
#' @return An immutable `effect_measurement_bridge`.
#' @export
measurement_bridge <- function(left_leg, right_leg, left_domain, right_domain,
                               common_space, provenance = list()) {
  if (!is.matrix(left_leg) || !is.numeric(left_leg) ||
      !is.matrix(right_leg) || !is.numeric(right_leg) ||
      any(dim(left_leg) < 1L) || any(dim(right_leg) < 1L) ||
      any(!is.finite(left_leg)) || any(!is.finite(right_leg))) {
    stop("Bridge legs must be finite nonempty numeric matrices.",
      call. = FALSE)
  }
  left_leg <- unname(left_leg)
  right_leg <- unname(right_leg)
  storage.mode(left_leg) <- "double"
  storage.mode(right_leg) <- "double"
  left_domain <- .domain_reference(left_domain)
  right_domain <- .domain_reference(right_domain)
  common_space <- .validate_measurement_space(common_space)
  .validate_effect_provenance(provenance, "measurement-bridge provenance")
  if (nrow(left_leg) != common_space$n_measurements ||
      nrow(right_leg) != common_space$n_measurements) {
    stop("Bridge legs must share the declared common measurement dimension.",
      call. = FALSE)
  }
  if (ncol(left_leg) != left_domain$n_features ||
      ncol(right_leg) != right_domain$n_features) {
    stop("Bridge leg widths must match their exact neural domains.",
      call. = FALSE)
  }
  semantic <- list(
    schema_version = 1L,
    kind = "factorized",
    left_domain = left_domain,
    right_domain = right_domain,
    common_space = common_space,
    left_leg = left_leg,
    right_leg = right_leg,
    provenance = provenance
  )
  structure(c(semantic[-1L], list(
    signature = paste0("sha256:", digest::digest(
      semantic, algo = "sha256", serialize = TRUE, serializeVersion = 2L
    ))
  )), class = "effect_measurement_bridge")
}

.validate_factorized_measurement_bridge <- function(bridge,
                                                    left_domain = NULL,
                                                    right_domain = NULL) {
  expected <- c("kind", "left_domain", "right_domain", "common_space",
    "left_leg", "right_leg", "provenance", "signature")
  if (!inherits(bridge, "effect_measurement_bridge") || !is.list(bridge) ||
      !identical(names(bridge), expected) ||
      !identical(bridge$kind, "factorized")) {
    stop("Factorized measurement-bridge fields are missing or noncanonical.",
      call. = FALSE)
  }
  rebuilt <- measurement_bridge(
    bridge$left_leg, bridge$right_leg,
    bridge$left_domain, bridge$right_domain,
    bridge$common_space, bridge$provenance
  )
  if (!identical(bridge, rebuilt)) {
    stop("Factorized measurement bridge identity is inconsistent.",
      call. = FALSE)
  }
  if (!is.null(left_domain) &&
      !.same_domain_reference(bridge$left_domain,
        .domain_reference(left_domain))) {
    stop("Bridge left neural-space identity is incompatible.", call. = FALSE)
  }
  if (!is.null(right_domain) &&
      !.same_domain_reference(bridge$right_domain,
        .domain_reference(right_domain))) {
    stop("Bridge right neural-space identity is incompatible.", call. = FALSE)
  }
  rebuilt
}

.validate_measurement_bridge <- function(bridge, left = NULL, right = NULL) {
  if (!inherits(bridge, "effect_measurement_bridge")) {
    stop("`bridge` must be a crossform measurement bridge.", call. = FALSE)
  }
  if (identical(bridge$kind, "identity")) {
    if (is.null(left) || is.null(right)) {
      stop("Identity bridge validation requires both relation sides.",
        call. = FALSE)
    }
    return(.validate_identity_measurement_bridge(bridge, left, right))
  }
  left_domain <- if (is.null(left)) NULL else
    if (inherits(left, "effect_relation")) left$domain else left
  right_domain <- if (is.null(right)) NULL else
    if (inherits(right, "effect_relation")) right$domain else right
  .validate_factorized_measurement_bridge(
    bridge, left_domain, right_domain
  )
}

#' Reverse a factorized measurement bridge
#'
#' @param bridge A factorized `measurement_bridge()`.
#' @return The exact bridge with source domains and legs exchanged.
#' @export
reverse_bridge <- function(bridge) {
  bridge <- .validate_factorized_measurement_bridge(bridge)
  measurement_bridge(
    bridge$right_leg, bridge$left_leg,
    bridge$right_domain, bridge$left_domain,
    bridge$common_space,
    provenance = bridge$provenance
  )
}

.apply_measurement_bridge <- function(left_relations, right_relations, bridge) {
  bridge <- .validate_factorized_measurement_bridge(bridge)
  validate_family <- function(relations, width, side) {
    if (!is.list(relations) || length(relations) < 1L ||
        is.null(names(relations))) {
      stop(sprintf("Bridged %s relations must be a named family.", side),
        call. = FALSE)
    }
    for (value in relations) {
      if (!is.matrix(value) || !is.numeric(value) || ncol(value) != width ||
          any(!is.finite(value))) {
        stop(sprintf("Bridged %s relations have invalid neural width.", side),
          call. = FALSE)
      }
    }
  }
  validate_family(left_relations, ncol(bridge$left_leg), "left")
  validate_family(right_relations, ncol(bridge$right_leg), "right")
  list(
    left = lapply(left_relations, function(value) value %*% t(bridge$left_leg)),
    right = lapply(right_relations, function(value) value %*% t(bridge$right_leg)),
    common_space = bridge$common_space,
    bridge_signature = bridge$signature
  )
}

.bridge_joint_gram <- function(left, right, bridge) {
  bridge <- .validate_factorized_measurement_bridge(bridge)
  if (!is.matrix(left) || !is.matrix(right) ||
      ncol(left) != ncol(bridge$left_leg) ||
      ncol(right) != ncol(bridge$right_leg)) {
    stop("Joint bridge inputs do not match the bridge legs.", call. = FALSE)
  }
  measured <- rbind(left %*% t(bridge$left_leg),
    right %*% t(bridge$right_leg))
  measured %*% t(measured)
}
