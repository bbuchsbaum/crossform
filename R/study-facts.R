# First-moment factual records ---------------------------------------------

.validate_partition_names <- function(value, expected = NULL,
                                      name = "partitions") {
  if (!is.character(value) || length(value) < 1L || anyNA(value) ||
      any(!nzchar(value)) || anyDuplicated(value)) {
    stop(sprintf("`%s` must contain unique nonempty identifiers.", name),
      call. = FALSE)
  }
  if (!is.null(expected) && length(value) != expected) {
    stop(sprintf("`%s` has the wrong length.", name), call. = FALSE)
  }
  unname(value)
}

.canonical_fact_table <- function(value, name) {
  if (!is.data.frame(value) || nrow(value) < 1L || ncol(value) < 1L ||
      is.null(names(value)) || anyNA(names(value)) || any(!nzchar(names(value))) ||
      anyDuplicated(names(value))) {
    stop(sprintf("`%s` must be a nonempty data frame with unique column names.",
      name), call. = FALSE)
  }
  bad <- vapply(value, function(column) {
    is.list(column) || is.matrix(column) || is.data.frame(column) ||
      anyNA(column) || (is.numeric(column) && any(!is.finite(column)))
  }, logical(1))
  if (any(bad)) {
    stop(sprintf(
      "`%s` columns must be complete finite atomic or factor values; invalid: %s.",
      name, paste(names(value)[bad], collapse = ", ")
    ), call. = FALSE)
  }
  rownames(value) <- NULL
  value
}

.fact_table_schema <- function(value) {
  lapply(value, function(column) list(
    class = class(column),
    storage = typeof(column),
    levels = if (is.factor(column)) levels(column) else NULL
  ))
}

#' Define one partition's observation axis
#'
#' @param observation_id Unique ordered observation identifiers.
#' @param partition One partition identifier.
#' @param time Optional strictly increasing finite observation times.
#' @param units One physical time unit when `time` is supplied.
#' @param provenance Portable acquisition provenance.
#' @return An `effect_observation_index`.
#' @export
observation_index <- function(observation_id, partition, time = NULL,
                              units = NULL, provenance = list()) {
  if (length(observation_id) < 1L || anyNA(observation_id) ||
      anyDuplicated(observation_id)) {
    stop("`observation_id` must contain unique complete identifiers.",
      call. = FALSE)
  }
  if (!(is.character(observation_id) || is.numeric(observation_id) ||
      is.integer(observation_id))) {
    stop("`observation_id` must be character or numeric.", call. = FALSE)
  }
  partition <- .validate_nonempty_id(partition, "partition")
  if (is.null(time)) {
    if (!is.null(units)) {
      stop("`units` requires an observation `time` axis.", call. = FALSE)
    }
    timing <- FALSE
  } else {
    if (!is.numeric(time) || length(time) != length(observation_id) ||
        anyNA(time) || any(!is.finite(time)) || any(diff(time) <= 0)) {
      stop("`time` must be finite, strictly increasing, and match observations.",
        call. = FALSE)
    }
    units <- .validate_nonempty_id(units, "units")
    time <- as.numeric(time)
    timing <- TRUE
  }
  .validate_effect_provenance(provenance, "observation-index provenance")
  semantic <- list(
    schema_version = 1L,
    observation_id = unname(observation_id),
    partition = partition,
    time = time,
    units = units,
    timing = timing,
    provenance = provenance
  )
  structure(c(semantic[-1L], list(
    signature = .semantic_digest("observation-index-sha256:", semantic)
  )), class = "effect_observation_index")
}

.validate_observation_index <- function(value) {
  expected <- c(
    "observation_id", "partition", "time", "units", "timing",
    "provenance", "signature"
  )
  if (!inherits(value, "effect_observation_index") || !is.list(value) ||
      !identical(names(value), expected)) {
    stop("Observation-index fields are missing or noncanonical.", call. = FALSE)
  }
  rebuilt <- observation_index(
    value$observation_id, value$partition, value$time, value$units,
    value$provenance
  )
  if (!identical(value, rebuilt)) {
    stop("Observation-index metadata or signature is inconsistent.",
      call. = FALSE)
  }
  rebuilt
}

.observation_source <- function(source, dim) {
  if (is.matrix(source)) return(.matrix_response_source(source))
  if (inherits(source, "effect_source_descriptor")) {
    return(.descriptor_response_source(source))
  }
  .function_response_source(source, dim)
}

#' Bind raw observation sources to indexes and a neural domain
#'
#' The constructor validates only source metadata and axes. Neural values are
#' read later by [estimate_relation()].
#'
#' @param sources A matrix or named list of matrix, function, or
#'   `effect_source_descriptor` sources.
#' @param index One [observation_index()] or a named list, one per partition.
#' @param domain The exact neural domain.
#' @param source_dims Dimensions required for function sources.
#' @param partitions Optional partition order.
#' @param capabilities Optional [source_capabilities()] per partition. Function
#'   sources require explicit capabilities.
#' @param provenance Portable observation provenance.
#' @return An `effect_observations` fact object.
#' @export
observations <- function(sources, index, domain, source_dims = NULL,
                         partitions = NULL, capabilities = NULL,
                         provenance = list()) {
  if (is.matrix(sources) || inherits(sources, "effect_source_descriptor") ||
      is.function(sources)) {
    sources <- list(sources)
  }
  if (!is.list(sources) || length(sources) < 1L) {
    stop("`sources` must provide at least one observation source.",
      call. = FALSE)
  }
  if (inherits(index, "effect_observation_index")) index <- list(index)
  if (!is.list(index) || length(index) != length(sources)) {
    stop("`index` must provide one observation index per source.",
      call. = FALSE)
  }
  index <- lapply(index, .validate_observation_index)
  if (is.null(partitions)) {
    partitions <- names(sources)
    if (is.null(partitions) || any(!nzchar(partitions))) {
      partitions <- vapply(index, `[[`, character(1), "partition")
    }
  }
  partitions <- .validate_partition_names(
    partitions, length(sources), "partitions"
  )
  names(sources) <- names(index) <- partitions
  for (partition in partitions) {
    if (!identical(index[[partition]]$partition, partition)) {
      stop(sprintf(
        "Observation index `%s` declares partition `%s`.",
        partition, index[[partition]]$partition
      ), call. = FALSE)
    }
  }
  domain <- .domain_reference(domain)
  .validate_effect_provenance(provenance, "observations provenance")

  if (is.null(source_dims)) source_dims <- vector("list", length(sources))
  if (!is.list(source_dims) || length(source_dims) != length(sources)) {
    stop("`source_dims` must provide one entry per source.", call. = FALSE)
  }
  compiled <- Map(.observation_source, sources, source_dims)
  names(compiled) <- partitions
  for (partition in partitions) {
    expected <- c(
      length(index[[partition]]$observation_id), domain$n_features
    )
    if (!identical(compiled[[partition]]$dim, as.integer(expected))) {
      stop(sprintf(
        paste0(
          "Partition `%s` source dimensions are %d x %d; its observation ",
          "index and neural domain require %d x %d."
        ),
        partition,
        compiled[[partition]]$dim[[1L]], compiled[[partition]]$dim[[2L]],
        expected[[1L]], expected[[2L]]
      ), call. = FALSE)
    }
  }

  if (is.null(capabilities)) {
    if (any(vapply(compiled, function(value) is.null(value$descriptor),
      logical(1)))) {
      stop("Function observation sources require explicit `capabilities`.",
        call. = FALSE)
    }
    capabilities <- lapply(compiled, function(value) {
      descriptor <- value$descriptor
      source_capabilities(
        block_read = TRUE,
        reopenable = descriptor$access %in% c("reopenable", "shared"),
        thread_safe = FALSE,
        stable_revision = descriptor$stable_revision
      )
    })
  } else if (inherits(capabilities, "effect_source_capabilities")) {
    capabilities <- rep(list(capabilities), length(partitions))
  }
  if (!is.list(capabilities) || length(capabilities) != length(partitions)) {
    stop("`capabilities` must provide one value per observation source.",
      call. = FALSE)
  }
  capabilities <- lapply(capabilities, .validate_source_capabilities)
  names(capabilities) <- partitions
  for (partition in partitions) {
    descriptor <- compiled[[partition]]$descriptor
    capability <- capabilities[[partition]]
    if (!is.null(descriptor) && !identical(
        tolower(descriptor$stable_revision),
        tolower(capability$stable_revision)
      )) {
      stop(sprintf(
        "Partition `%s` source descriptor and capability revisions differ.",
        partition
      ), call. = FALSE)
    }
    if (isTRUE(capability$reopenable) &&
        (is.null(descriptor) || identical(descriptor$access, "coordinator"))) {
      stop("Reopenable capabilities require a reopenable descriptor.",
        call. = FALSE)
    }
  }

  source_facts <- lapply(partitions, function(partition) {
    source <- compiled[[partition]]
    list(
      kind = source$kind,
      dim = source$dim,
      descriptor = source$descriptor,
      capabilities = capabilities[[partition]]
    )
  })
  names(source_facts) <- partitions
  semantic <- list(
    schema_version = 1L,
    partitions = partitions,
    indexes = index,
    domain = domain,
    sources = source_facts,
    provenance = provenance
  )
  structure(list(
    sources = compiled,
    indexes = index,
    partitions = partitions,
    domain = domain,
    n_features = domain$n_features,
    capabilities = capabilities,
    provenance = provenance,
    observations_id = .semantic_digest("observations-sha256:", semantic)
  ), class = "effect_observations")
}

.validate_observations <- function(value, deep = TRUE) {
  expected <- c(
    "sources", "indexes", "partitions", "domain", "n_features",
    "capabilities", "provenance", "observations_id"
  )
  if (!inherits(value, "effect_observations") || !is.list(value) ||
      !identical(names(value), expected) ||
      !identical(names(value$sources), value$partitions) ||
      !identical(names(value$indexes), value$partitions) ||
      !identical(names(value$capabilities), value$partitions)) {
    stop("Observation fields are missing or noncanonical.", call. = FALSE)
  }
  .validate_partition_names(value$partitions)
  domain <- if (isTRUE(deep)) {
    .validate_domain_reference(value$domain)
  } else {
    .validate_domain_reference_shallow(value$domain)
  }
  if (!identical(value$n_features, domain$n_features)) {
    stop("Observation neural-domain metadata are inconsistent.", call. = FALSE)
  }
  source_facts <- lapply(value$partitions, function(partition) {
    index <- .validate_observation_index(value$indexes[[partition]])
    capability <- .validate_source_capabilities(value$capabilities[[partition]])
    source <- value$sources[[partition]]
    if (!inherits(source, "effect_response_source") || !is.list(source) ||
        !is.numeric(source$dim) || length(source$dim) != 2L ||
        !identical(source$dim, as.integer(c(
          length(index$observation_id), domain$n_features
        ))) || !is.function(source$read)) {
      stop("Observation source metadata are inconsistent.", call. = FALSE)
    }
    descriptor <- source$descriptor
    if (!is.null(descriptor)) {
      descriptor <- .validate_source_descriptor(descriptor)
      if (!identical(descriptor$dim, source$dim) || !identical(
          tolower(descriptor$stable_revision),
          tolower(capability$stable_revision)
        )) {
        stop("Observation source revision metadata are inconsistent.",
          call. = FALSE)
      }
    }
    list(
      kind = source$kind,
      dim = source$dim,
      descriptor = descriptor,
      capabilities = capability
    )
  })
  names(source_facts) <- value$partitions
  .validate_effect_provenance(value$provenance, "observations provenance")
  semantic <- list(
    schema_version = 1L,
    partitions = value$partitions,
    indexes = value$indexes,
    domain = domain,
    sources = source_facts,
    provenance = value$provenance
  )
  if (!identical(value$observations_id,
      .semantic_digest("observations-sha256:", semantic))) {
    stop("Observation identity is inconsistent.", call. = FALSE)
  }
  value
}

#' Declare a typed event record
#'
#' @param data A nonempty event data frame.
#' @param partition Column naming the observation partition.
#' @param event_id Column containing event identifiers.
#' @param onset,duration Optional timing columns. Supply both or neither.
#' @param units Physical time unit when timing columns are supplied.
#' @param provenance Portable event provenance.
#' @return An `effect_events` fact object. Column roles remain model-relative.
#' @export
observation_events <- function(data, partition = "partition",
                               event_id = "event_id", onset = "onset",
                               duration = "duration", units = "seconds",
                               provenance = list()) {
  data <- .canonical_fact_table(data, "data")
  identifiers <- c(partition, event_id)
  if (!all(identifiers %in% names(data))) {
    stop("Event data must contain the partition and event-id columns.",
      call. = FALSE)
  }
  partition <- .validate_nonempty_id(partition, "partition")
  event_id <- .validate_nonempty_id(event_id, "event_id")
  if (xor(is.null(onset), is.null(duration))) {
    stop("Supply both event timing columns or neither.", call. = FALSE)
  }
  timing <- !is.null(onset)
  if (timing) {
    onset <- .validate_nonempty_id(onset, "onset")
    duration <- .validate_nonempty_id(duration, "duration")
    if (!all(c(onset, duration) %in% names(data)) ||
        !is.numeric(data[[onset]]) || !is.numeric(data[[duration]]) ||
        any(data[[duration]] < 0)) {
      stop("Event onset and duration columns must be finite numeric values with nonnegative duration.",
        call. = FALSE)
    }
    units <- .validate_nonempty_id(units, "units")
  } else {
    units <- NULL
  }
  keys <- paste(data[[partition]], data[[event_id]], sep = "\r")
  if (anyDuplicated(keys)) {
    stop("Event identifiers must be unique within partition.", call. = FALSE)
  }
  .validate_effect_provenance(provenance, "event provenance")
  semantic <- list(
    schema_version = 1L,
    data = data,
    schema = .fact_table_schema(data),
    partition_column = partition,
    event_id_column = event_id,
    onset_column = onset,
    duration_column = duration,
    units = units,
    timing = timing,
    provenance = provenance
  )
  structure(c(semantic[-1L], list(
    events_id = .semantic_digest("events-sha256:", semantic)
  )), class = "effect_events")
}

.validate_events <- function(value) {
  expected <- c(
    "data", "schema", "partition_column", "event_id_column", "onset_column",
    "duration_column", "units", "timing", "provenance", "events_id"
  )
  if (!inherits(value, "effect_events") || !is.list(value) ||
      !identical(names(value), expected)) {
    stop("Event fields are missing or noncanonical.", call. = FALSE)
  }
  rebuilt <- observation_events(
    value$data,
    partition = value$partition_column,
    event_id = value$event_id_column,
    onset = value$onset_column,
    duration = value$duration_column,
    units = value$units,
    provenance = value$provenance
  )
  if (!identical(value, rebuilt)) {
    stop("Event metadata or identity is inconsistent.", call. = FALSE)
  }
  rebuilt
}

#' Declare observation-level confounds and censor facts
#'
#' @param data A nonempty data frame with one row per observation.
#' @param partition,observation_id Columns binding rows to observation indexes.
#' @param censor Optional logical censor column. `TRUE` means retained.
#' @param provenance Portable confound provenance.
#' @return An `effect_observation_confounds` fact object.
#' @export
observation_confounds <- function(
    data, partition = "partition", observation_id = "observation_id",
    censor = NULL, provenance = list()) {
  data <- .canonical_fact_table(data, "data")
  if (!all(c(partition, observation_id) %in% names(data))) {
    stop("Confounds must contain partition and observation-id columns.",
      call. = FALSE)
  }
  partition <- .validate_nonempty_id(partition, "partition")
  observation_id <- .validate_nonempty_id(observation_id, "observation_id")
  if (!is.null(censor)) {
    censor <- .validate_nonempty_id(censor, "censor")
    if (!censor %in% names(data) || !is.logical(data[[censor]])) {
      stop("`censor` must identify a complete logical column.", call. = FALSE)
    }
  }
  keys <- paste(data[[partition]], data[[observation_id]], sep = "\r")
  if (anyDuplicated(keys)) {
    stop("Confound rows must be unique by partition and observation.",
      call. = FALSE)
  }
  .validate_effect_provenance(provenance, "confound provenance")
  semantic <- list(
    schema_version = 1L,
    data = data,
    schema = .fact_table_schema(data),
    partition_column = partition,
    observation_id_column = observation_id,
    censor_column = censor,
    provenance = provenance
  )
  structure(c(semantic[-1L], list(
    confounds_id = .semantic_digest("observation-confounds-sha256:", semantic)
  )), class = "effect_observation_confounds")
}

.validate_observation_confounds <- function(value) {
  expected <- c(
    "data", "schema", "partition_column", "observation_id_column",
    "censor_column", "provenance", "confounds_id"
  )
  if (!inherits(value, "effect_observation_confounds") || !is.list(value) ||
      !identical(names(value), expected)) {
    stop("Observation-confound fields are missing or noncanonical.",
      call. = FALSE)
  }
  rebuilt <- observation_confounds(
    value$data,
    partition = value$partition_column,
    observation_id = value$observation_id_column,
    censor = value$censor_column,
    provenance = value$provenance
  )
  if (!identical(value, rebuilt)) {
    stop("Observation-confound metadata or identity is inconsistent.",
      call. = FALSE)
  }
  rebuilt
}

#' Declare nested partition axes
#'
#' Columns are ordered from the leaf partition outward. Each child level must
#' map to exactly one parent level.
#'
#' @param data A data frame with one row per leaf partition and one column per
#'   nested axis.
#' @param leaf Name of the leaf partition column; defaults to the first column.
#' @param provenance Portable hierarchy provenance.
#' @return An `effect_partition_hierarchy`.
#' @export
partition_hierarchy <- function(data, leaf = names(data)[[1L]],
                                provenance = list()) {
  data <- .canonical_fact_table(data, "data")
  if (!is.character(leaf) || length(leaf) != 1L || is.na(leaf) ||
      !leaf %in% names(data) || !identical(leaf, names(data)[[1L]])) {
    stop("`leaf` must name the first hierarchy column.", call. = FALSE)
  }
  for (column in names(data)) {
    data[[column]] <- as.character(data[[column]])
    if (any(!nzchar(data[[column]]))) {
      stop("Hierarchy levels must be nonempty identifiers.", call. = FALSE)
    }
  }
  if (anyDuplicated(data[[leaf]])) {
    stop("The leaf hierarchy axis must identify each partition once.",
      call. = FALSE)
  }
  parent_maps <- list()
  if (ncol(data) > 1L) {
    for (index in seq_len(ncol(data) - 1L)) {
      child <- names(data)[[index]]
      parent <- names(data)[[index + 1L]]
      split_parent <- split(data[[parent]], data[[child]])
      if (any(lengths(lapply(split_parent, unique)) != 1L)) {
        stop(sprintf(
          "Hierarchy axis `%s` does not nest uniquely within `%s`.",
          child, parent
        ), call. = FALSE)
      }
      parent_maps[[child]] <- vapply(split_parent, unique, character(1))
    }
  }
  levels <- lapply(data, function(column) unique(unname(column)))
  .validate_effect_provenance(provenance, "partition-hierarchy provenance")
  semantic <- list(
    schema_version = 1L,
    data = data,
    leaf = leaf,
    axes = names(data),
    levels = levels,
    parent_maps = parent_maps,
    provenance = provenance
  )
  structure(c(semantic[-1L], list(
    signature = .semantic_digest("partition-hierarchy-sha256:", semantic)
  )), class = "effect_partition_hierarchy")
}

.validate_partition_hierarchy <- function(value) {
  expected <- c(
    "data", "leaf", "axes", "levels", "parent_maps", "provenance",
    "signature"
  )
  if (!inherits(value, "effect_partition_hierarchy") || !is.list(value) ||
      !identical(names(value), expected)) {
    stop("Partition-hierarchy fields are missing or noncanonical.",
      call. = FALSE)
  }
  rebuilt <- partition_hierarchy(
    value$data, leaf = value$leaf, provenance = value$provenance
  )
  if (!identical(value, rebuilt)) {
    stop("Partition-hierarchy metadata or signature is inconsistent.",
      call. = FALSE)
  }
  rebuilt
}
