# Lazy experimental-neural relations --------------------------------------

.matrix_response_source <- function(value) {
  if (!is.matrix(value) || !is.numeric(value) || any(dim(value) < 1L) ||
      any(!is.finite(value))) {
    stop("Matrix response sources must be finite and nonempty.", call. = FALSE)
  }
  revision <- paste0(
    "sha256:", digest::digest(value, algo = "sha256", serialize = TRUE)
  )
  structure(
    list(
      dim = as.integer(dim(value)),
      read = function(features) value[, features, drop = FALSE],
      kind = "matrix",
      descriptor = .source_descriptor(
        kind = "matrix",
        dim = dim(value),
        access = "coordinator",
        stable_revision = revision,
        spec = list()
      )
    ),
    class = "effect_response_source"
  )
}

.function_response_source <- function(read, dim) {
  if (!is.function(read)) stop("Function sources must be functions.", call. = FALSE)
  if (!is.numeric(dim) || length(dim) != 2L || any(!is.finite(dim)) ||
      any(dim < 1L) || any(dim %% 1 != 0)) {
    stop("Function source dimensions must be two positive integers.", call. = FALSE)
  }
  structure(
    list(dim = as.integer(dim), read = read, kind = "function", descriptor = NULL),
    class = "effect_response_source"
  )
}

.descriptor_response_source <- function(descriptor) {
  descriptor <- .validate_source_descriptor(descriptor)
  read <- if (identical(descriptor$kind, "file_matrix")) {
    function(features) {
      .with_source_descriptor(descriptor, function(handle) handle$read(features))
    }
  } else {
    function(features) {
      stop("Shared response sources require an explicit executor adapter.",
        call. = FALSE)
    }
  }
  structure(
    list(
      dim = descriptor$dim,
      read = read,
      kind = descriptor$kind,
      descriptor = descriptor
    ),
    class = "effect_response_source"
  )
}

#' Construct a lazy experimental-neural relation
#'
#' Each source is read only by neural feature block. Extractors map its
#' observation rows into one common named experimental space. Omitting
#' `extract` declares that sources already contain effect-by-feature matrices.
#'
#' @param sources A matrix, or a named list of matrix, function, or
#'   `effect_source_descriptor` response sources.
#' @param extract NULL, one `effect_extractor`, or one extractor per partition.
#' @param effects Effect names for already estimated sources.
#' @param source_dims Required dimensions for function sources, as one
#'   two-element vector per partition.
#' @param partitions Optional partition names when not supplied by `sources`.
#' @param domain Optional `effect_domain`; its identity is recorded without
#'   resampling or changing neural features.
#' @param domain_id Stable neural-domain identity.
#' @param capabilities Optional source-capability values, one per partition.
#' @param provenance Optional provenance metadata.
#' @return An `effect_relation`.
#' @export
relation <- function(sources, extract = NULL, effects = NULL,
                     source_dims = NULL, partitions = NULL,
                     domain = NULL, domain_id = "abstract", capabilities = NULL,
                     provenance = list()) {
  if (is.matrix(sources)) sources <- list(sources)
  if (!is.list(sources) || length(sources) < 1L) {
    stop("`sources` must be a matrix or nonempty list.", call. = FALSE)
  }
  if (is.null(partitions)) partitions <- names(sources)
  if (is.null(partitions) || any(!nzchar(partitions))) {
    partitions <- paste0("partition", seq_along(sources))
  }
  if (!is.character(partitions) || length(partitions) != length(sources) ||
      anyNA(partitions) || any(!nzchar(partitions)) || anyDuplicated(partitions)) {
    stop("Partitions must have unique nonempty names.", call. = FALSE)
  }
  names(sources) <- partitions
  if (!is.null(domain)) {
    .validate_domain(domain)
    if (!identical(domain_id, "abstract") && !identical(domain_id, domain$id)) {
      stop("`domain` and `domain_id` identify different neural domains.",
        call. = FALSE)
    }
    domain_id <- domain$id
  }
  if (!is.character(domain_id) || length(domain_id) != 1L ||
      is.na(domain_id) || !nzchar(domain_id)) {
    stop("`domain_id` must be one nonempty identifier.", call. = FALSE)
  }
  if (!is.list(provenance)) stop("`provenance` must be a list.", call. = FALSE)

  if (is.null(source_dims)) source_dims <- vector("list", length(sources))
  if (!is.list(source_dims) || length(source_dims) != length(sources)) {
    stop("`source_dims` must provide one entry per source.", call. = FALSE)
  }
  compiled_sources <- Map(function(source, dim) {
    if (is.matrix(source)) {
      .matrix_response_source(source)
    } else if (inherits(source, "effect_source_descriptor")) {
      .descriptor_response_source(source)
    } else {
      .function_response_source(source, dim)
    }
  }, sources, source_dims)
  names(compiled_sources) <- partitions

  if (is.null(extract)) {
    q <- compiled_sources[[1L]]$dim[[1L]]
    inferred_effects <- if (is.null(effects) && is.matrix(sources[[1L]])) {
      rownames(sources[[1L]])
    } else {
      effects
    }
    effects <- .validate_effect_names(inferred_effects, q)
    extractors <- lapply(compiled_sources, function(source) {
      if (source$dim[[1L]] != q) {
        stop("Precomputed effect sources must share one effect dimension.",
          call. = FALSE)
      }
      effect_extractor(diag(q), effects, estimator = "identity")
    })
  } else {
    extractors <- if (inherits(extract, "effect_extractor")) {
      rep(list(extract), length(sources))
    } else {
      extract
    }
    if (!is.list(extractors) || length(extractors) != length(sources)) {
      stop("`extract` must supply one extractor per partition.", call. = FALSE)
    }
    extractors <- lapply(extractors, .validate_effect_extractor)
    effects <- extractors[[1L]]$effects
  }
  names(extractors) <- partitions
  for (partition in partitions) {
    if (extractors[[partition]]$n_observations !=
        compiled_sources[[partition]]$dim[[1L]]) {
      stop("Extractor observations must match their response source.", call. = FALSE)
    }
    if (!identical(extractors[[partition]]$effects, effects)) {
      stop("All partition extractors must share one named effect space.",
        call. = FALSE)
    }
  }
  feature_counts <- vapply(compiled_sources, function(x) x$dim[[2L]], integer(1))
  if (length(unique(feature_counts)) != 1L) {
    stop("All relation sources must share one neural feature dimension.",
      call. = FALSE)
  }
  if (is.null(capabilities) && all(vapply(compiled_sources, function(source) {
    !is.null(source$descriptor)
  }, logical(1)))) {
    capabilities <- lapply(compiled_sources, function(source) {
      descriptor <- source$descriptor
      source_capabilities(
        block_read = TRUE,
        reopenable = descriptor$access %in% c("reopenable", "shared"),
        thread_safe = FALSE,
        stable_revision = descriptor$stable_revision
      )
    })
  }
  if (!is.null(capabilities)) {
    if (inherits(capabilities, "effect_source_capabilities")) {
      capabilities <- rep(list(capabilities), length(sources))
    }
    if (!is.list(capabilities) || length(capabilities) != length(sources)) {
      stop("`capabilities` must supply one value per partition.", call. = FALSE)
    }
    capabilities <- lapply(capabilities, .validate_source_capabilities)
    names(capabilities) <- partitions
    for (partition in partitions) {
      descriptor <- compiled_sources[[partition]]$descriptor
      capability <- capabilities[[partition]]
      if (isTRUE(capability$reopenable) &&
          (is.null(descriptor) || identical(descriptor$access, "coordinator"))) {
        stop("Reopenable source capabilities require a reopenable descriptor.",
          call. = FALSE)
      }
      if (!is.null(descriptor) &&
          !identical(tolower(capability$stable_revision),
            tolower(descriptor$stable_revision))) {
        stop("Source capability and descriptor revisions must agree.",
          call. = FALSE)
      }
    }
  }

  structure(
    list(
      sources = compiled_sources,
      extractors = extractors,
      effects = effects,
      partitions = partitions,
      n_features = feature_counts[[1L]],
      domain_id = domain_id,
      capabilities = capabilities,
      provenance = provenance
    ),
    class = "effect_relation"
  )
}

#' Read one experimental-neural relation block
#'
#' @param x An `effect_relation`.
#' @param partition One partition name or index.
#' @param features Unique neural feature indices.
#' @return A finite effect-by-feature matrix.
#' @export
relation_block <- function(x, partition, features) {
  .validate_relation(x)
  if (is.numeric(partition) && length(partition) == 1L && !is.na(partition) &&
      partition %% 1 == 0 && partition >= 1L && partition <= length(x$partitions)) {
    partition <- x$partitions[[as.integer(partition)]]
  }
  if (!is.character(partition) || length(partition) != 1L ||
      is.na(partition) || !partition %in% x$partitions) {
    stop("`partition` must identify one relation partition.", call. = FALSE)
  }
  if (!is.numeric(features) || length(features) < 1L || anyNA(features) ||
      any(!is.finite(features)) || any(features %% 1 != 0) ||
      any(features < 1L) || any(features > x$n_features) ||
      anyDuplicated(features)) {
    stop("`features` must be unique valid neural feature indices.", call. = FALSE)
  }
  features <- as.integer(features)
  source <- x$sources[[partition]]
  response <- source$read(features)
  if (!is.matrix(response) || !is.numeric(response) ||
      !identical(dim(response), c(source$dim[[1L]], length(features))) ||
      any(!is.finite(response))) {
    stop("Response source returned an invalid observation-by-feature block.",
      call. = FALSE)
  }
  value <- x$extractors[[partition]]$map %*% response
  if (any(!is.finite(value))) {
    stop("Effect extraction produced non-finite relation values.", call. = FALSE)
  }
  dimnames(value) <- list(x$effects, NULL)
  value
}

.validate_relation <- function(x) {
  expected <- c("sources", "extractors", "effects", "partitions", "n_features",
    "domain_id", "capabilities", "provenance")
  if (!inherits(x, "effect_relation") || !is.list(x) ||
      !identical(names(x), expected) || !is.list(x$sources) ||
      !is.list(x$extractors) || length(x$sources) < 1L ||
      length(x$sources) != length(x$extractors) ||
      !identical(names(x$sources), x$partitions) ||
      !identical(names(x$extractors), x$partitions)) {
    stop("Relation fields are missing or noncanonical.", call. = FALSE)
  }
  .validate_effect_names(x$effects, length(x$effects))
  if (!is.numeric(x$n_features) || length(x$n_features) != 1L ||
      is.na(x$n_features) || !is.finite(x$n_features) ||
      x$n_features < 1L || x$n_features %% 1 != 0) {
    stop("Relation feature metadata is invalid.", call. = FALSE)
  }
  for (partition in x$partitions) {
    source <- x$sources[[partition]]
    if (!inherits(source, "effect_response_source") ||
        !is.numeric(source$dim) || length(source$dim) != 2L ||
        source$dim[[2L]] != x$n_features || !is.function(source$read) ||
        !(is.null(source$descriptor) ||
          inherits(source$descriptor, "effect_source_descriptor"))) {
      stop("Relation response-source metadata is invalid.", call. = FALSE)
    }
    if (!is.null(source$descriptor)) {
      descriptor <- .validate_source_descriptor(source$descriptor)
      if (!identical(descriptor$dim, source$dim)) {
        stop("Relation source descriptor dimensions are inconsistent.",
          call. = FALSE)
      }
    }
    extractor <- .validate_effect_extractor(x$extractors[[partition]])
    if (extractor$n_observations != source$dim[[1L]] ||
        !identical(extractor$effects, x$effects)) {
      stop("Relation extractor metadata is inconsistent.", call. = FALSE)
    }
  }
  if (!is.null(x$capabilities)) {
    if (!is.list(x$capabilities) ||
        length(x$capabilities) != length(x$partitions) ||
        !identical(names(x$capabilities), x$partitions)) {
      stop("Relation source capabilities are inconsistent.", call. = FALSE)
    }
    lapply(x$capabilities, .validate_source_capabilities)
  }
  invisible(x)
}

.relation_source_descriptors <- function(x, require_reopenable = FALSE) {
  .validate_relation(x)
  descriptors <- lapply(x$sources, `[[`, "descriptor")
  names(descriptors) <- x$partitions
  if (require_reopenable && any(vapply(descriptors, function(descriptor) {
    is.null(descriptor) || identical(descriptor$access, "coordinator")
  }, logical(1)))) {
    stop("Worker-read execution requires reopenable or shared source descriptors.",
      call. = FALSE)
  }
  descriptors
}
