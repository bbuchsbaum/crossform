# First-moment factual records ---------------------------------------------

.validate_partition_names <- function(value, expected = NULL,
                                      name = "partitions") {
  if (!.is_strings(value, unique = TRUE) || length(value) < 1L) {
    .input_error(sprintf(
      "`%s` must contain unique nonempty identifiers%s.", name,
      if (is.character(value) && anyDuplicated(value)) {
        sprintf("; %s appears more than once",
          .msg_names(unique(value[duplicated(value)])))
      } else {
        sprintf("; received %s", .msg_value(value))
      }))
  }
  if (!is.null(expected) && length(value) != expected) {
    .input_error(sprintf("`%s` must supply %s; received %s.", name,
      .msg_count(expected, "name"), .msg_count(length(value), "name")))
  }
  unname(value)
}

.canonical_fact_table <- function(value, name) {
  if (!is.data.frame(value) || nrow(value) < 1L || ncol(value) < 1L ||
      is.null(names(value)) || anyNA(names(value)) || any(!nzchar(names(value))) ||
      anyDuplicated(names(value))) {
    .input_error(sprintf(
      "`%s` must be a nonempty data frame with unique column names; received %s.",
      name, .msg_value(value)))
  }
  bad <- vapply(value, function(column) {
    is.list(column) || is.matrix(column) || is.data.frame(column) ||
      anyNA(column) || (is.numeric(column) && any(!is.finite(column)))
  }, logical(1))
  if (any(bad)) {
    .input_error(sprintf(
      "`%s` columns must be complete finite atomic or factor values; invalid: %s.",
      name, paste(names(value)[bad], collapse = ", ")
    ))
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
#' Records what was acquired, in what order, and on which clock, for a single
#' partition such as a run. Build one index per response matrix before calling
#' [observations()].
#'
#' @param observation_id Unique ordered observation identifiers.
#' @param partition One partition identifier.
#' @param time Optional strictly increasing finite observation times.
#' @param units One physical time unit when `time` is supplied.
#' @param provenance Portable acquisition provenance.
#' @return An `effect_observation_index`: a list with `$observation_id`,
#'   `$partition`, `$time` and `$units` (both `NULL` when untimed), the
#'   `$timing` flag, `$provenance`, and a `$signature`.
#' @family typed observation facts
#' @seealso [observations()] to bind indexes to response sources, and
#'   [observation_events()] whose clock must use the same units.
#' @examples
#' # One acquisition axis: six scans on a 2 s clock.
#' index <- observation_index(
#'   seq_len(6L), partition = "run-1",
#'   time = seq(0, by = 2, length.out = 6L), units = "seconds"
#' )
#' index$partition
#' index$timing
#'
#' # The clock must be strictly increasing, so a repeated acquisition time is
#' # rejected here rather than silently reordered later.
#' try(observation_index(1:3, "run-1", time = c(0, 2, 2), units = "seconds"))
#' @export
observation_index <- function(observation_id, partition, time = NULL,
                              units = NULL, provenance = list()) {
  if (length(observation_id) < 1L || anyNA(observation_id) ||
      anyDuplicated(observation_id)) {
    .input_error("`observation_id` must contain unique complete identifiers.")
  }
  if (!(is.character(observation_id) || is.numeric(observation_id) ||
      is.integer(observation_id))) {
    .input_error("`observation_id` must be character or numeric.")
  }
  partition <- .validate_nonempty_id(partition, "partition")
  if (is.null(time)) {
    if (!is.null(units)) {
      .input_error("`units` requires an observation `time` axis.")
    }
    timing <- FALSE
  } else {
    if (!.is_finite_numeric(time) || length(time) != length(observation_id) ||
        anyNA(time) || any(diff(time) <= 0)) {
      .input_error(
        "`time` must be finite, strictly increasing, and match observations."
      )
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
    signature = .sha256_signature(semantic, "observation-index-sha256:")
  )), class = "effect_observation_index")
}

.validate_observation_index <- function(value) {
  expected <- c(
    "observation_id", "partition", "time", "units", "timing",
    "provenance", "signature"
  )
  if (!.sealed_fields(value, "effect_observation_index", expected)) {
    .input_error("Observation-index fields are missing or noncanonical.")
  }
  rebuilt <- observation_index(
    value$observation_id, value$partition, value$time, value$units,
    value$provenance
  )
  if (!identical(value, rebuilt)) {
    .contract_error("Observation-index metadata or signature is inconsistent.")
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
#' @return An `effect_observations` fact object: a list with the compiled
#'   `$sources`, `$indexes`, `$partitions`, the `$domain` reference,
#'   `$n_features`, per-partition `$capabilities`, `$provenance`, and an
#'   `$observations_id` covering all of them.
#' @family typed observation facts
#' @seealso [observation_index()] for the axes, [study()] to bind events and
#'   confounds to these observations, and [file_matrix_source()] for
#'   out-of-memory sources.
#' @examples
#' set.seed(1)
#' domain <- abstract_domain(4L, id = "observations-example")
#' indexes <- list(
#'   `run-1` = observation_index(
#'     1:6, "run-1", time = seq(0, 10, by = 2), units = "seconds"
#'   ),
#'   `run-2` = observation_index(
#'     1:6, "run-2", time = seq(0, 10, by = 2), units = "seconds"
#'   )
#' )
#' responses <- lapply(indexes, function(index) matrix(rnorm(24), 6L, 4L))
#'
#' record <- observations(responses, indexes, domain)
#' record$partitions
#' record$n_features
#'
#' # No neural values were read: only shapes, axes, and source revisions were
#' # checked, so a mis-shaped run is caught before any fit is attempted.
#' short <- responses
#' short$`run-2` <- short$`run-2`[1:5, , drop = FALSE]
#' try(observations(short, indexes, domain))
#' @export
observations <- function(sources, index, domain, source_dims = NULL,
                         partitions = NULL, capabilities = NULL,
                         provenance = list()) {
  if (is.matrix(sources) || inherits(sources, "effect_source_descriptor") ||
      is.function(sources)) {
    sources <- list(sources)
  }
  if (missing(index) || missing(domain)) {
    .input_error(paste0(
      "`index` and `domain` are both required: `observations()` binds one ",
      "`observation_index()` per source to the neural domain those sources ",
      "were sampled on."
    ))
  }
  if (!is.list(sources) || length(sources) < 1L) {
    .input_error(sprintf(paste0(
      "`sources` must provide at least one observation-by-feature source, ",
      "as a matrix, function, descriptor, or a list of them; received %s."
    ), .msg_value(sources)))
  }
  if (inherits(index, "effect_observation_index")) index <- list(index)
  if (!is.list(index) || length(index) != length(sources)) {
    .input_error(sprintf(paste0(
      "`index` must provide one `observation_index()` per source: received ",
      "%s for %s."
    ), if (is.list(index)) .msg_count(length(index), "index") else
      .msg_value(index),
      .msg_count(length(sources), "source")))
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
      .contract_error(sprintf(
        "Observation index `%s` declares partition `%s`.",
        partition, index[[partition]]$partition
      ))
    }
  }
  domain <- .domain_reference(domain)
  .validate_effect_provenance(provenance, "observations provenance")

  if (is.null(source_dims)) source_dims <- vector("list", length(sources))
  if (!is.list(source_dims) || length(source_dims) != length(sources)) {
    .input_error("`source_dims` must provide one entry per source.")
  }
  compiled <- Map(.observation_source, sources, source_dims)
  names(compiled) <- partitions
  for (partition in partitions) {
    expected <- c(
      length(index[[partition]]$observation_id), domain$n_features
    )
    if (!identical(compiled[[partition]]$dim, as.integer(expected))) {
      .input_error(sprintf(
        paste0(
          "Partition `%s` source dimensions are %d x %d; its observation ",
          "index and neural domain require %d x %d."
        ),
        partition,
        compiled[[partition]]$dim[[1L]], compiled[[partition]]$dim[[2L]],
        expected[[1L]], expected[[2L]]
      ))
    }
  }

  if (is.null(capabilities)) {
    if (any(vapply(compiled, function(value) is.null(value$descriptor),
      logical(1)))) {
      .input_error(
        "Function observation sources require explicit `capabilities`."
      )
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
    .input_error(
      "`capabilities` must provide one value per observation source."
    )
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
      .contract_error(sprintf(
        "Partition `%s` source descriptor and capability revisions differ.",
        partition
      ))
    }
    if (isTRUE(capability$reopenable) &&
        (is.null(descriptor) || identical(descriptor$access, "coordinator"))) {
      .input_error("Reopenable capabilities require a reopenable descriptor.")
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
    observations_id = .sha256_signature(semantic, "observations-sha256:")
  ), class = "effect_observations")
}

.validate_observations <- function(value, deep = TRUE) {
  expected <- c(
    "sources", "indexes", "partitions", "domain", "n_features",
    "capabilities", "provenance", "observations_id"
  )
  if (!.sealed_fields(value, "effect_observations", expected) ||
      !identical(names(value$sources), value$partitions) ||
      !identical(names(value$indexes), value$partitions) ||
      !identical(names(value$capabilities), value$partitions)) {
    .input_error("Observation fields are missing or noncanonical.")
  }
  .validate_partition_names(value$partitions)
  domain <- if (isTRUE(deep)) {
    .validate_domain_reference(value$domain)
  } else {
    .validate_domain_reference_shallow(value$domain)
  }
  if (!identical(value$n_features, domain$n_features)) {
    .contract_error("Observation neural-domain metadata are inconsistent.")
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
      .input_error("Observation source metadata are inconsistent.")
    }
    descriptor <- source$descriptor
    if (!is.null(descriptor)) {
      descriptor <- .validate_source_descriptor(descriptor)
      if (!identical(descriptor$dim, source$dim) || !identical(
          tolower(descriptor$stable_revision),
          tolower(capability$stable_revision)
        )) {
        .input_error("Observation source revision metadata are inconsistent.")
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
      .sha256_signature(semantic, "observations-sha256:"))) {
    .contract_error("Observation identity is inconsistent.")
  }
  value
}

#' Declare a typed event record
#'
#' Records the experimental events as facts, keeping every supplied column. It
#' assigns no model roles: whether a column is a target, a nuisance term, or
#' unused is decided later by the design model.
#'
#' @param data A nonempty event data frame.
#' @param partition Column naming the observation partition.
#' @param event_id Column containing event identifiers.
#' @param onset,duration Optional timing columns. Supply both or neither.
#' @param units Physical time unit when timing columns are supplied.
#' @param provenance Portable event provenance.
#' @return An `effect_events` fact object: a list with the canonical `$data`
#'   and its `$schema`, the `$partition_column`, `$event_id_column`,
#'   `$onset_column` and `$duration_column` role names, `$units`, the `$timing`
#'   flag, `$provenance`, and an `$events_id`. Column roles remain
#'   model-relative.
#' @family typed observation facts
#' @seealso [study()] to bind events to observations, [bids_events()] to read
#'   them from BIDS TSV files, and [observation_confounds()] for the
#'   observation-level table.
#' @examples
#' events <- data.frame(
#'   partition = rep(c("run-1", "run-2"), each = 2L),
#'   event_id = paste0("e", 1:4),
#'   onset = c(0, 6, 0, 6),
#'   duration = 0.5,
#'   condition = c("face", "body", "face", "body")
#' )
#' record <- observation_events(events)
#' record$timing
#' record$units
#'
#' # `condition` is preserved but carries no model role yet; the design model,
#' # not this record, decides what is a target or a nuisance term.
#' names(record$data)
#'
#' # Event ids must be unique within a partition.
#' duplicated_ids <- events
#' duplicated_ids$event_id <- c("e1", "e1", "e3", "e4")
#' try(observation_events(duplicated_ids))
#' @export
observation_events <- function(data, partition = "partition",
                               event_id = "event_id", onset = "onset",
                               duration = "duration", units = "seconds",
                               provenance = list()) {
  if (missing(data)) {
    .input_error(paste0(
      "`data` is required: pass an event table with one row per event and, ",
      "at minimum, the partition and event-id columns."
    ))
  }
  data <- .canonical_fact_table(data, "data")
  identifiers <- c(partition, event_id)
  if (!all(identifiers %in% names(data))) {
    .input_error(sprintf(paste0(
      "`data` is missing the %s column%s. `observation_events()` reads the ",
      "partition and event id from columns you name; `data` has %s."
    ), .msg_names(setdiff(identifiers, names(data))),
      if (length(setdiff(identifiers, names(data))) == 1L) "" else "s",
      .msg_names(names(data))))
  }
  partition <- .validate_nonempty_id(partition, "partition")
  event_id <- .validate_nonempty_id(event_id, "event_id")
  if (xor(is.null(onset), is.null(duration))) {
    .input_error("Supply both event timing columns or neither.")
  }
  timing <- !is.null(onset)
  if (timing) {
    onset <- .validate_nonempty_id(onset, "onset")
    duration <- .validate_nonempty_id(duration, "duration")
    if (!all(c(onset, duration) %in% names(data)) ||
        !is.numeric(data[[onset]]) || !is.numeric(data[[duration]]) ||
        any(data[[duration]] < 0)) {
      .input_error(paste0(
        "Event onset and duration columns must be finite numeric values with ",
        "nonnegative duration."
      ))
    }
    units <- .validate_nonempty_id(units, "units")
  } else {
    units <- NULL
  }
  keys <- paste(data[[partition]], data[[event_id]], sep = "\r")
  if (anyDuplicated(keys)) {
    .input_error("Event identifiers must be unique within partition.")
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
    events_id = .sha256_signature(semantic, "events-sha256:")
  )), class = "effect_events")
}

.validate_events <- function(value) {
  expected <- c(
    "data", "schema", "partition_column", "event_id_column", "onset_column",
    "duration_column", "units", "timing", "provenance", "events_id"
  )
  if (!.sealed_fields(value, "effect_events", expected)) {
    .input_error("Event fields are missing or noncanonical.")
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
    .contract_error("Event metadata or identity is inconsistent.")
  }
  rebuilt
}

#' Declare observation-level confounds and censor facts
#'
#' Records one row per observation. Censoring is never inferred from motion or
#' outlier columns: to exclude observations you must name an explicit logical
#' retain column.
#'
#' @param data A nonempty data frame with one row per observation.
#' @param partition,observation_id Columns binding rows to observation indexes.
#' @param censor Optional logical censor column. `TRUE` means retained.
#' @param provenance Portable confound provenance.
#' @return An `effect_observation_confounds` fact object: a list with the
#'   canonical `$data` and its `$schema`, the `$partition_column`,
#'   `$observation_id_column` and `$censor_column` role names, `$provenance`,
#'   and a `$confounds_id`.
#' @family typed observation facts
#' @seealso [study()], which joins these rows to the observation axis and
#'   records the resulting row lineage, and [bids_confounds()] for
#'   fMRIPrep-style TSV input.
#' @examples
#' confounds <- data.frame(
#'   partition = "run-1",
#'   observation_id = 1:6,
#'   framewise_displacement = c(0.1, 0.2, 0.9, 0.1, 0.1, 0.3),
#'   retained = c(TRUE, TRUE, FALSE, TRUE, TRUE, TRUE)
#' )
#'
#' # Naming the censor column is what makes the exclusion a declared fact.
#' record <- observation_confounds(confounds, censor = "retained")
#' record$censor_column
#' sum(record$data$retained)
#'
#' # Without `censor`, the motion column is kept but nothing is excluded: no
#' # censor policy is inferred from it.
#' is.null(observation_confounds(confounds)$censor_column)
#' @export
observation_confounds <- function(
    data, partition = "partition", observation_id = "observation_id",
    censor = NULL, provenance = list()) {
  if (missing(data)) {
    .input_error(paste0(
      "`data` is required: pass a confound table with one row per ",
      "observation and, at minimum, the partition and observation-id columns."
    ))
  }
  data <- .canonical_fact_table(data, "data")
  if (!all(c(partition, observation_id) %in% names(data))) {
    .input_error(sprintf(paste0(
      "`data` is missing the %s column%s. `observation_confounds()` reads the ",
      "partition and observation id from columns you name; `data` has %s."
    ), .msg_names(setdiff(c(partition, observation_id), names(data))),
      if (length(setdiff(c(partition, observation_id), names(data))) == 1L) {
        ""
      } else {
        "s"
      }, .msg_names(names(data))))
  }
  partition <- .validate_nonempty_id(partition, "partition")
  observation_id <- .validate_nonempty_id(observation_id, "observation_id")
  if (!is.null(censor)) {
    censor <- .validate_nonempty_id(censor, "censor")
    if (!censor %in% names(data) || !is.logical(data[[censor]])) {
      .input_error("`censor` must identify a complete logical column.")
    }
  }
  keys <- paste(data[[partition]], data[[observation_id]], sep = "\r")
  if (anyDuplicated(keys)) {
    .input_error("Confound rows must be unique by partition and observation.")
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
    confounds_id = .sha256_signature(semantic, "observation-confounds-sha256:")
  )), class = "effect_observation_confounds")
}

.validate_observation_confounds <- function(value) {
  expected <- c(
    "data", "schema", "partition_column", "observation_id_column",
    "censor_column", "provenance", "confounds_id"
  )
  if (!.sealed_fields(value, "effect_observation_confounds", expected)) {
    .input_error("Observation-confound fields are missing or noncanonical.")
  }
  rebuilt <- observation_confounds(
    value$data,
    partition = value$partition_column,
    observation_id = value$observation_id_column,
    censor = value$censor_column,
    provenance = value$provenance
  )
  if (!identical(value, rebuilt)) {
    .contract_error(
      "Observation-confound metadata or identity is inconsistent."
    )
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
#' @return An `effect_partition_hierarchy`: a list with the canonical `$data`,
#'   the `$leaf` column name, the ordered `$axes`, the distinct `$levels` per
#'   axis, the child-to-parent `$parent_maps`, `$provenance`, and a
#'   `$signature`.
#' @family typed observation facts
#' @seealso [study()], which binds the hierarchy to observations, and
#'   [study_axis()] to select one axis from it.
#' @examples
#' # Runs nested in sessions nested in one subject, leaf column first.
#' hierarchy <- partition_hierarchy(data.frame(
#'   partition = c("run-1", "run-2", "run-3", "run-4"),
#'   session = c("ses-1", "ses-1", "ses-2", "ses-2"),
#'   subject = "sub-01"
#' ))
#' hierarchy$axes
#' hierarchy$parent_maps$partition
#'
#' # Nesting must be exact: a session that belongs to two subjects is refused
#' # rather than quietly flattened.
#' crossed <- data.frame(
#'   partition = c("run-1", "run-2"),
#'   session = c("ses-1", "ses-1"),
#'   subject = c("sub-01", "sub-02")
#' )
#' try(partition_hierarchy(crossed))
#' @export
partition_hierarchy <- function(data, leaf = names(data)[[1L]],
                                provenance = list()) {
  data <- .canonical_fact_table(data, "data")
  if (!.is_string(leaf, allow_empty = TRUE) || !leaf %in% names(data) ||
      !identical(leaf, names(data)[[1L]])) {
    .input_error("`leaf` must name the first hierarchy column.")
  }
  for (column in names(data)) {
    data[[column]] <- as.character(data[[column]])
    if (any(!nzchar(data[[column]]))) {
      .input_error("Hierarchy levels must be nonempty identifiers.")
    }
  }
  if (anyDuplicated(data[[leaf]])) {
    .input_error("The leaf hierarchy axis must identify each partition once.")
  }
  parent_maps <- list()
  if (ncol(data) > 1L) {
    for (index in seq_len(ncol(data) - 1L)) {
      child <- names(data)[[index]]
      parent <- names(data)[[index + 1L]]
      split_parent <- split(data[[parent]], data[[child]])
      if (any(lengths(lapply(split_parent, unique)) != 1L)) {
        .input_error(sprintf(
          "Hierarchy axis `%s` does not nest uniquely within `%s`.",
          child, parent
        ))
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
    signature = .sha256_signature(semantic, "partition-hierarchy-sha256:")
  )), class = "effect_partition_hierarchy")
}

.validate_partition_hierarchy <- function(value) {
  expected <- c(
    "data", "leaf", "axes", "levels", "parent_maps", "provenance",
    "signature"
  )
  if (!.sealed_fields(value, "effect_partition_hierarchy", expected)) {
    .input_error("Partition-hierarchy fields are missing or noncanonical.")
  }
  rebuilt <- partition_hierarchy(
    value$data, leaf = value$leaf, provenance = value$provenance
  )
  if (!identical(value, rebuilt)) {
    .contract_error(
      "Partition-hierarchy metadata or signature is inconsistent."
    )
  }
  rebuilt
}
