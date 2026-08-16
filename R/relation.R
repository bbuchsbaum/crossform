# Lazy experimental-neural relations --------------------------------------

.matrix_response_source <- function(value, label = NULL) {
  where <- if (is.null(label)) "" else sprintf("Partition `%s`: ", label)
  if (!is.matrix(value) || !is.numeric(value)) {
    stop(sprintf("%sa response source must be a numeric matrix; received %s.",
      where, .msg_value(value)), call. = FALSE)
  }
  if (any(dim(value) < 1L)) {
    stop(sprintf(
      "%sa response source must be nonempty; received a %d x %d matrix.",
      where, nrow(value), ncol(value)), call. = FALSE)
  }
  if (any(!is.finite(value))) {
    bad <- !is.finite(value)
    rows <- which(apply(bad, 1L, any))
    stop(sprintf(paste0(
      "%s%d of %d values %s non-finite (NA, NaN, or Inf), in %s %s. ",
      "crossform never drops or imputes them silently: repair or exclude ",
      "those features before building the relation."
    ), where, sum(bad), length(value),
      if (sum(bad) == 1L) "is" else "are",
      if (length(rows) == 1L) "row" else "rows",
      .msg_positions(apply(bad, 1L, any), rownames(value))), call. = FALSE)
  }
  revision <- .sha256_signature(value)
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
#' A relation made from precomputed effects supports point evidence but carries
#' no residual error channel. Beta matrices alone cannot recover within-
#' participant analytic uncertainty. If later work needs a learned neural
#' metric or [rdm_sampling_covariance()], start from raw responses with
#' [lm_relation_fit()]. Subject-level resampling answers a population question;
#' it does not recreate discarded first-level residual information.
#'
#' @param sources A matrix, or a named list of matrix, function, or
#'   `effect_source_descriptor` response sources.
#' @param extract NULL, one `effect_extractor`, or one extractor per partition.
#' @param effects An `effect_space()` for already estimated sources, or unique
#'   names used as shorthand for an unspecified-basis effect space.
#' @param source_dims Required dimensions for function sources, as one
#'   two-element vector per partition.
#' @param partitions Optional partition names when not supplied by `sources`.
#' @param domain Optional `effect_domain`; its identity is recorded without
#'   resampling or changing neural features.
#' @param domain_id Stable neural-domain identity.
#' @param capabilities Optional source-capability values, one per partition.
#' @param provenance Optional provenance metadata.
#' @return An `effect_relation`: a list with the compiled `$sources`, one
#'   `$extractors` entry per partition, the shared `$effect_space` and its
#'   `$effects` labels, `$partitions`, `$n_features`, the `$domain` reference
#'   and `$domain_id`, optional source `$capabilities`, and `$provenance`. No
#'   neural values have been read.
#' @section Structure:
#' A relation declares what will be read, so its elements are identities and
#' shapes rather than values.
#'
#' - `$effects`: the effect labels, in the order every partition estimates
#'   them and the order unnamed contrast weights are read in.
#' - `$effect_space`: the shared `effect_space()` those labels coordinate,
#'   carrying the basis identity, units, and signature every partition
#'   agrees on.
#' - `$partitions`: the partition names, in the order they were supplied.
#'   These are what [cross_partitions()] and [pairing()] name.
#' - `$n_features`: the number of neural features every partition spans.
#' - `$domain`: the neural domain the source columns are bound to, carrying
#'   its `id`, `n_features`, and `feature_ids`; `$domain_id` repeats that
#'   `id`.
#' - `$capabilities`: one [source_capabilities()] per partition when the
#'   sources declared them, otherwise `NULL`.
#' - `$provenance`: the metadata supplied at construction, unchanged.
#'
#' `$sources` and `$extractors` are the compiled read path [relation_block()]
#' uses; they and any element not listed here are internal and may change.
#' @family relation planning and fitting
#' @seealso [lm_relation_fit()] when raw responses are available and a residual
#'   channel is needed, [relation_block()] to read one block, and
#'   [plan_geometry()] for the next step.
#' @examples
#' # Precomputed condition betas: one effect-by-feature matrix per run.
#' set.seed(1)
#' domain <- abstract_domain(5L, id = "relation-example")
#' betas <- lapply(c("run-1", "run-2"), function(partition) {
#'   matrix(rnorm(15), 3L, 5L,
#'     dimnames = list(c("face", "body", "tool"), NULL))
#' })
#' names(betas) <- c("run-1", "run-2")
#'
#' point <- relation(betas, domain = domain)
#' point$effects
#' point$partitions
#'
#' # Betas alone carry no residual information, so anything needing the error
#' # channel refuses here rather than inventing a standard error.
#' catch_refusal(residual_df(point, "run-1"))$capability
#' @export
relation <- function(sources, extract = NULL, effects = NULL,
                     source_dims = NULL, partitions = NULL,
                     domain = NULL, domain_id = "abstract", capabilities = NULL,
                     provenance = list()) {
  if (missing(sources)) {
    stop(paste0(
      "`sources` is required: pass one effect-by-feature matrix per ",
      "partition, named by effect, as a named list such as ",
      "`relation(list(run1 = betas1, run2 = betas2), domain = domain)`."
    ), call. = FALSE)
  }
  if (is.matrix(sources)) sources <- list(sources)
  if (!is.list(sources) || length(sources) < 1L) {
    stop(sprintf(paste0(
      "`sources` must be one effect-by-feature matrix or a nonempty named ",
      "list with one entry per partition; received %s."
    ), .msg_value(sources)), call. = FALSE)
  }
  if (is.null(partitions)) partitions <- names(sources)
  if (is.null(partitions) || any(!nzchar(partitions))) {
    partitions <- paste0("partition", seq_along(sources))
  }
  if (!is.character(partitions) || length(partitions) != length(sources) ||
      anyNA(partitions) || any(!nzchar(partitions)) || anyDuplicated(partitions)) {
    stop(sprintf(paste0(
      "Partitions must have unique nonempty names; received %s for %s%s."
    ), .msg_count(length(partitions), "name"),
      .msg_count(length(sources), "source"),
      if (anyDuplicated(partitions)) {
        sprintf(", with %s repeated",
          .msg_names(unique(partitions[duplicated(partitions)])))
      } else {
        ""
      }), call. = FALSE)
  }
  names(sources) <- partitions
  domain_reference <- NULL
  if (!is.null(domain)) {
    domain_reference <- .domain_reference(domain)
    if (!identical(domain_id, "abstract") &&
        !identical(domain_id, domain_reference$id)) {
      stop("`domain` and `domain_id` identify different neural domains.",
        call. = FALSE)
    }
    domain_id <- domain_reference$id
  }
  if (!is.character(domain_id) || length(domain_id) != 1L ||
      is.na(domain_id) || !nzchar(domain_id)) {
    stop("`domain_id` must be one nonempty identifier.", call. = FALSE)
  }
  .validate_effect_provenance(provenance, "relation provenance")

  if (is.null(source_dims)) source_dims <- vector("list", length(sources))
  if (!is.list(source_dims) || length(source_dims) != length(sources)) {
    stop("`source_dims` must provide one entry per source.", call. = FALSE)
  }
  declared_effect_space <- !is.null(effects)
  if (is.null(extract)) {
    first_dimension <- if (is.matrix(sources[[1L]])) {
      nrow(sources[[1L]])
    } else if (inherits(sources[[1L]], "effect_source_descriptor")) {
      .validate_source_descriptor(sources[[1L]])$dim[[1L]]
    } else {
      source_dims[[1L]][[1L]]
    }
    if (is.null(effects)) {
      if (!all(vapply(sources, is.matrix, logical(1)))) {
        stop("Non-matrix precomputed sources require an explicit effect space.",
          call. = FALSE)
      }
      names_first <- .matrix_effect_names(sources[[1L]], required = TRUE,
        label = partitions[[1L]])
      effects <- effect_space(names_first)
    } else {
      effects <- .as_effect_space(effects, first_dimension)
    }
    for (partition in partitions) {
      source <- sources[[partition]]
      source_dimension <- if (is.matrix(source)) nrow(source) else if (
        inherits(source, "effect_source_descriptor")) {
        .validate_source_descriptor(source)$dim[[1L]]
      } else {
        source_dims[[match(partition, partitions)]][[1L]]
      }
      if (!identical(as.integer(source_dimension),
          as.integer(length(effects$coordinates)))) {
        stop(sprintf(paste0(
          "Partition `%s` has %s but the effect space declares %s (%s). ",
          "Every partition must estimate the same effects."
        ), partition, .msg_count(source_dimension, "row"),
          .msg_count(length(effects$coordinates), "effect"),
          .msg_names(effects$coordinates)), call. = FALSE)
      }
      if (is.matrix(source)) {
        source_names <- .matrix_effect_names(source,
          required = !declared_effect_space, label = partition)
        if (!is.null(source_names)) {
          missing <- setdiff(effects$coordinates, source_names)
          extra <- setdiff(source_names, effects$coordinates)
          if (length(missing) || length(extra)) {
            detail <- c(
              if (length(missing)) sprintf("%s absent", .msg_names(missing)),
              if (length(extra)) sprintf("%s unexpected", .msg_names(extra))
            )
            stop(sprintf(paste0(
              "Partition `%s` does not carry the declared effects (%s): %s."
            ), partition, .msg_names(effects$coordinates),
              paste(detail, collapse = "; ")), call. = FALSE)
          }
          source <- source[match(effects$coordinates, source_names), , drop = FALSE]
        } else {
          rownames(source) <- effects$coordinates
        }
        sources[[partition]] <- source
      }
    }
  }

  compiled_sources <- Map(function(source, dim, label) {
    if (is.matrix(source)) {
      .matrix_response_source(source, label)
    } else if (inherits(source, "effect_source_descriptor")) {
      .descriptor_response_source(source)
    } else {
      .function_response_source(source, dim)
    }
  }, sources, source_dims, partitions)
  names(compiled_sources) <- partitions

  if (is.null(extract)) {
    q <- compiled_sources[[1L]]$dim[[1L]]
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
    extractor_space <- extractors[[1L]]$effect_space
    if (!is.null(effects)) {
      declared <- .as_effect_space(effects, length(extractor_space$coordinates))
      if (!.same_effect_space(declared, extractor_space)) {
        stop("The declared and extractor effect spaces are incompatible.",
          call. = FALSE)
      }
    }
    effects <- extractor_space
  }
  names(extractors) <- partitions
  for (partition in partitions) {
    if (extractors[[partition]]$n_observations !=
        compiled_sources[[partition]]$dim[[1L]]) {
      stop("Extractor observations must match their response source.", call. = FALSE)
    }
    if (!.same_effect_space(extractors[[partition]]$effect_space, effects)) {
      stop("All partition extractors must share one identical effect space.",
        call. = FALSE)
    }
  }
  feature_counts <- vapply(compiled_sources, function(x) x$dim[[2L]], integer(1))
  if (length(unique(feature_counts)) != 1L) {
    stop(sprintf(paste0(
      "All relation sources must span the same neural features, but their ",
      "column counts differ: %s. Every partition must be estimated on one ",
      "common domain."
    ), paste(sprintf("`%s` has %d", partitions, feature_counts),
      collapse = ", ")), call. = FALSE)
  }
  if (is.null(domain_reference)) {
    domain_reference <- .positional_domain_reference(feature_counts[[1L]],
      domain_id)
  } else if (!identical(domain_reference$n_features, feature_counts[[1L]])) {
    stop(sprintf(paste0(
      "The relation sources have %s but domain `%s` declares %s. The source ",
      "columns must be exactly the domain's features, in domain order."
    ), .msg_count(feature_counts[[1L]], "column"), domain_reference$id,
      .msg_count(domain_reference$n_features, "feature")), call. = FALSE)
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
      effect_space = effects,
      effects = effects$coordinates,
      partitions = partitions,
      n_features = feature_counts[[1L]],
      domain = domain_reference,
      domain_id = domain_reference$id,
      capabilities = capabilities,
      provenance = provenance
    ),
    class = "effect_relation"
  )
}

#' Read one experimental-neural relation block
#'
#' This is where neural values are actually read. Only the requested feature
#' columns are pulled from the source and mapped through that partition's
#' extractor, which is what keeps out-of-memory sources bounded.
#'
#' @param x An `effect_relation` or `effect_relation_fit`.
#' @param partition One partition name or index.
#' @param features Unique neural feature indices.
#' @return A finite effect-by-feature matrix with one row per effect
#'   coordinate, in effect-space order, and one column per requested feature.
#' @family relation planning and fitting
#' @seealso [relation()] and [lm_relation_fit()] for the objects it reads, and
#'   [residual_block()] for the matching residual read.
#' @examples
#' set.seed(1)
#' domain <- abstract_domain(5L, id = "relation-block-example")
#' betas <- matrix(rnorm(15), 3L, 5L,
#'   dimnames = list(c("face", "body", "tool"), NULL))
#' point <- relation(list(`run-1` = betas), domain = domain)
#'
#' # Ask for two features only; the source is read for those columns alone.
#' round(relation_block(point, "run-1", c(1L, 4L)), 3)
#'
#' # Partitions may be named or given by index, and rows always come back in
#' # effect-space order regardless of how the source was stored.
#' rownames(relation_block(point, 1L, seq_len(5L)))
#' @export
relation_block <- function(x, partition, features) {
  if (inherits(x, "effect_relation_fit")) {
    .validate_relation_fit(x, deep = FALSE)
    x <- x$relation
  }
  .relation_block_with_reader(x, partition, features,
    function(partition, features) x$sources[[partition]]$read(features))
}

.relation_block_with_reader <- function(x, partition, features, read_response,
                                        validate = TRUE) {
  if (isTRUE(validate)) .validate_relation(x, deep = FALSE)
  if (!is.function(read_response)) {
    stop("Relation response reader must be a function.", call. = FALSE)
  }
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
  response <- read_response(partition, features)
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
  dimnames(value) <- list(x$effect_space$coordinates, NULL)
  value
}

.validate_relation <- function(x, deep = TRUE) {
  if (!inherits(x, "effect_relation")) {
    stop(sprintf(paste0(
      "Expected an `effect_relation` from `relation()` or the `$relation` of ",
      "an `lm_relation_fit()`; received %s."
    ), .msg_value(x)), call. = FALSE)
  }
  expected <- c("sources", "extractors", "effect_space", "effects",
    "partitions", "n_features", "domain", "domain_id", "capabilities",
    "provenance")
  if (!inherits(x, "effect_relation") || !is.list(x) ||
      !identical(names(x), expected) || !is.list(x$sources) ||
      !is.list(x$extractors) || length(x$sources) < 1L ||
      length(x$sources) != length(x$extractors) ||
      !identical(names(x$sources), x$partitions) ||
      !identical(names(x$extractors), x$partitions)) {
    stop("Relation fields are missing or noncanonical.", call. = FALSE)
  }
  effect_space <- if (isTRUE(deep)) {
    .validate_effect_space(x$effect_space)
  } else {
    .validate_effect_space_shallow(x$effect_space)
  }
  if (!identical(x$effects, effect_space$coordinates)) {
    stop("Relation coordinate labels are inconsistent with its effect space.",
      call. = FALSE)
  }
  if (!is.numeric(x$n_features) || length(x$n_features) != 1L ||
      is.na(x$n_features) || !is.finite(x$n_features) ||
      x$n_features < 1L || x$n_features %% 1 != 0) {
    stop("Relation feature metadata is invalid.", call. = FALSE)
  }
  domain <- if (isTRUE(deep)) {
    .validate_domain_reference(x$domain)
  } else {
    .validate_domain_reference_shallow(x$domain)
  }
  if (!identical(domain$n_features, as.integer(x$n_features)) ||
      !identical(domain$id, x$domain_id)) {
    stop("Relation metadata is inconsistent with its exact neural domain.",
      call. = FALSE)
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
      descriptor <- if (isTRUE(deep)) {
        .validate_source_descriptor(source$descriptor)
      } else {
        source$descriptor
      }
      if (!inherits(descriptor, "effect_source_descriptor") ||
          !is.list(descriptor) ||
          !identical(names(descriptor),
            c("kind", "dim", "access", "stable_revision", "spec"))) {
        stop("Relation source descriptor metadata are invalid.",
          call. = FALSE)
      }
      if (!identical(descriptor$dim, source$dim)) {
        stop("Relation source descriptor dimensions are inconsistent.",
          call. = FALSE)
      }
    }
    extractor <- if (isTRUE(deep)) {
      .validate_effect_extractor(x$extractors[[partition]])
    } else {
      .validate_effect_extractor_shallow(x$extractors[[partition]])
    }
    if (extractor$n_observations != source$dim[[1L]] ||
        (!identical(extractor$effect_space$signature, effect_space$signature) ||
         !identical(extractor$effect_space, effect_space))) {
      stop("Relation extractor metadata is inconsistent.", call. = FALSE)
    }
  }
  if (!is.null(x$capabilities)) {
    if (!is.list(x$capabilities) ||
        length(x$capabilities) != length(x$partitions) ||
        !identical(names(x$capabilities), x$partitions)) {
      stop("Relation source capabilities are inconsistent.", call. = FALSE)
    }
    if (isTRUE(deep)) {
      lapply(x$capabilities, .validate_source_capabilities)
    } else if (!all(vapply(x$capabilities, function(value) {
      inherits(value, "effect_source_capabilities") && is.list(value) &&
        identical(names(value), c("block_read", "reopenable", "thread_safe",
          "stable_revision")) && .strong_sha256(value$stable_revision)
    }, logical(1)))) {
      stop("Relation source capabilities are inconsistent.", call. = FALSE)
    }
  }
  if (isTRUE(deep)) {
    .validate_effect_provenance(x$provenance, "relation provenance")
  } else if (!is.list(x$provenance)) {
    stop("Relation provenance must be a list.", call. = FALSE)
  }
  invisible(x)
}

.validate_effect_space_shallow <- function(x) {
  expected <- c("coordinates", "basis_id", "units", "scale", "provenance",
    "signature")
  if (!inherits(x, "effect_space") || !is.list(x) ||
      !identical(names(x), expected) || !is.character(x$coordinates) ||
      length(x$coordinates) < 1L || anyNA(x$coordinates) ||
      any(!nzchar(x$coordinates)) || anyDuplicated(x$coordinates) ||
      !.strong_sha256(x$signature)) {
    stop("Effect-space fields are missing or noncanonical.", call. = FALSE)
  }
  x
}

.validate_domain_reference_shallow <- function(x) {
  expected <- c("id", "n_features", "feature_ids", "coordinate_units",
    "geometry_signature", "signature")
  if (!inherits(x, "effect_domain_reference") || !is.list(x) ||
      !identical(names(x), expected) || !is.character(x$id) ||
      length(x$id) != 1L || is.na(x$id) || !nzchar(x$id) ||
      !is.integer(x$n_features) || length(x$n_features) != 1L ||
      x$n_features < 1L || length(x$feature_ids) != x$n_features ||
      !.strong_sha256(x$signature)) {
    stop("Domain-reference fields are missing or noncanonical.",
      call. = FALSE)
  }
  x
}

.validate_effect_extractor_shallow <- function(x) {
  expected <- c("map", "effect_space", "effects", "n_observations",
    "estimator", "diagnostics")
  if (!inherits(x, "effect_extractor") || !is.list(x) ||
      !identical(names(x), expected) || !is.matrix(x$map) ||
      !is.numeric(x$map) || any(dim(x$map) < 1L) ||
      !is.integer(x$n_observations) || length(x$n_observations) != 1L ||
      x$n_observations != ncol(x$map)) {
    stop("Extractor fields are missing or noncanonical.", call. = FALSE)
  }
  effect_space <- .validate_effect_space_shallow(x$effect_space)
  if (!identical(x$effects, effect_space$coordinates) ||
      nrow(x$map) != length(effect_space$coordinates)) {
    stop("Extractor coordinate metadata are inconsistent.", call. = FALSE)
  }
  x
}

.relation_family_identity <- function(x) {
  .validate_relation(x)
  capabilities <- .relation_source_capabilities(x)
  semantic <- list(
    schema_version = 1L,
    source_revisions = vapply(capabilities, `[[`, character(1),
      "stable_revision"),
    extractors = lapply(x$extractors, function(value) {
      list(
        map = value$map,
        effect_space = value$effect_space,
        estimator = value$estimator,
        diagnostics = value$diagnostics
      )
    }),
    effect_space = x$effect_space,
    partitions = x$partitions,
    n_features = x$n_features,
    domain = x$domain,
    provenance = x$provenance
  )
  .sha256_signature(semantic)
}

.identity_measurement_bridge <- function(left, right) {
  .validate_relation(left)
  .validate_relation(right)
  if (!.same_domain_reference(left$domain, right$domain) ||
      !identical(as.integer(left$n_features), as.integer(right$n_features))) {
    stop(paste0(
      "Distinct neural spaces require an explicit compatible measurement ",
      "bridge before source admission."
    ), call. = FALSE)
  }
  semantic <- list(
    schema_version = 1L,
    kind = "identity",
    left_domain = left$domain,
    right_domain = right$domain,
    common_space = left$domain,
    n_features = as.integer(left$n_features),
    provenance = list(rule = "exact_shared_neural_feature_identity")
  )
  signature <- .sha256_signature(semantic)
  structure(
    c(semantic[-1L], list(signature = signature)),
    class = "effect_measurement_bridge"
  )
}

.validate_identity_measurement_bridge <- function(bridge, left, right) {
  expected <- c(
    "kind", "left_domain", "right_domain", "common_space", "n_features",
    "provenance", "signature"
  )
  if (!inherits(bridge, "effect_measurement_bridge") || !is.list(bridge) ||
      !identical(names(bridge), expected) || !identical(bridge$kind, "identity")) {
    stop("Identity measurement-bridge fields are missing or noncanonical.",
      call. = FALSE)
  }
  rebuilt <- .identity_measurement_bridge(left, right)
  if (!identical(bridge, rebuilt)) {
    stop("Identity measurement bridge is inconsistent with its relation sides.",
      call. = FALSE)
  }
  invisible(bridge)
}

.matrix_effect_names <- function(x, required = FALSE, label = NULL) {
  where <- if (is.null(label)) "A precomputed partition" else
    sprintf("Partition `%s`", label)
  value <- rownames(x)
  if (is.null(value)) {
    if (required) {
      stop(sprintf(paste0(
        "%s has no row names, so its %s cannot be identified. Name the rows ",
        "of every beta matrix, or declare the shared space with ",
        "`effects = effect_space(...)`."
      ), where, .msg_count(nrow(x), "effect")), call. = FALSE)
    }
    return(NULL)
  }
  if (length(value) != nrow(x) || anyNA(value) || any(!nzchar(value)) ||
      anyDuplicated(value)) {
    stop(sprintf(paste0(
      "%s has row names that are incomplete or repeated%s; each row must ",
      "name one distinct effect."
    ), where, if (anyDuplicated(value)) {
      sprintf(" (%s repeated)",
        .msg_names(unique(value[duplicated(value)])))
    } else {
      ""
    }), call. = FALSE)
  }
  value
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
