# Study binding -------------------------------------------------------------

.study_refusal <- function(message, capability, reasons, remedies) {
  .capability_refusal(
    message,
    capability = capability,
    namespace = "study",
    reasons = reasons,
    remedies = remedies
  )
}

.default_partition_hierarchy <- function(partitions) {
  partition_hierarchy(data.frame(
    partition = partitions,
    stringsAsFactors = FALSE
  ))
}

.event_partitions <- function(events) {
  if (is.null(events)) return(character())
  as.character(events$data[[events$partition_column]])
}

.bind_confounds <- function(observations, confounds) {
  lineage <- stats::setNames(vector("list", length(observations$partitions)),
    observations$partitions)
  for (partition in observations$partitions) {
    index <- observations$indexes[[partition]]
    ids <- index$observation_id
    if (is.null(confounds)) {
      lineage[[partition]] <- data.frame(
        source_row = seq_along(ids),
        observation_id = ids,
        confound_row = rep(NA_integer_, length(ids)),
        retained = rep(TRUE, length(ids)),
        check.names = FALSE
      )
      next
    }
    table <- confounds$data
    selected <- as.character(table[[confounds$partition_column]]) == partition
    rows <- which(selected)
    candidate <- table[[confounds$observation_id_column]][selected]
    if (!identical(class(candidate), class(ids)) ||
        !identical(typeof(candidate), typeof(ids))) {
      .study_refusal(
        sprintf("Partition `%s` uses incompatible observation-id types.",
          partition),
        capability = "aligned_observations",
        reasons = c(
          paste0("Observation index type: ", paste(class(ids), collapse = "/"),
            "."),
          paste0("Confound id type: ", paste(class(candidate), collapse = "/"),
            ".")
        ),
        remedies = paste0(
          "Encode both observation-id columns with the same type before ",
          "constructing `study()`."
        )
      )
    }
    matched <- match(ids, candidate)
    extra <- which(!candidate %in% ids)
    if (anyNA(matched) || length(extra)) {
      missing_ids <- ids[is.na(matched)]
      extra_ids <- candidate[extra]
      reasons <- c(
        if (length(missing_ids)) paste0(
          "Missing confound rows: ", paste(missing_ids, collapse = ", "), "."
        ),
        if (length(extra_ids)) paste0(
          "Extra confound rows: ", paste(extra_ids, collapse = ", "), "."
        )
      )
      .study_refusal(
        sprintf("Partition `%s` confounds do not cover its observation axis.",
          partition),
        capability = "aligned_observations",
        reasons = reasons,
        remedies = paste0(
          "Provide exactly one confound row for every observation in ",
          "partition `", partition, "`."
        )
      )
    }
    ordered_rows <- rows[matched]
    retained <- if (is.null(confounds$censor_column)) {
      rep(TRUE, length(ids))
    } else {
      table[[confounds$censor_column]][ordered_rows]
    }
    lineage[[partition]] <- data.frame(
      source_row = seq_along(ids),
      observation_id = ids,
      confound_row = as.integer(ordered_rows),
      retained = retained,
      check.names = FALSE
    )
  }
  lineage
}

.bind_event_clocks <- function(observations, events, tolerance) {
  coverage <- stats::setNames(vector("list", length(observations$partitions)),
    observations$partitions)
  if (is.null(events) || !isTRUE(events$timing)) {
    for (partition in observations$partitions) {
      coverage[[partition]] <- list(
        timing_resolved = FALSE,
        reason = if (is.null(events)) "no_event_record" else "untimed_events"
      )
    }
    return(list(resolved = FALSE, coverage = coverage))
  }

  event_partition <- as.character(events$data[[events$partition_column]])
  for (partition in observations$partitions) {
    index <- observations$indexes[[partition]]
    selected <- event_partition == partition
    if (!any(selected)) {
      coverage[[partition]] <- list(
        timing_resolved = isTRUE(index$timing),
        event_count = 0L,
        units = if (isTRUE(index$timing)) index$units else NULL
      )
      next
    }
    if (!isTRUE(index$timing)) {
      .study_refusal(
        sprintf("Partition `%s` has timed events but no observation clock.",
          partition),
        capability = "timing_resolved",
        reasons = paste0(
          "Events are measured in `", events$units,
          "`; the observation index has no time axis."
        ),
        remedies = paste0(
          "Supply observation times in `observation_index()` or use a ",
          "compiler for genuinely untimed reduced observations."
        )
      )
    }
    if (!identical(index$units, events$units)) {
      .study_refusal(
        sprintf("Partition `%s` event and observation clocks use different units.",
          partition),
        capability = "timing_resolved",
        reasons = c(
          paste0("Observation unit: `", index$units, "`."),
          paste0("Event unit: `", events$units, "`.")
        ),
        remedies = "Convert both clocks to one physical unit before binding."
      )
    }
    onset <- events$data[[events$onset_column]][selected]
    duration <- events$data[[events$duration_column]][selected]
    ending <- onset + duration
    lower <- min(index$time)
    upper <- max(index$time)
    outside <- onset < lower - tolerance | ending > upper + tolerance
    if (any(outside)) {
      event_ids <- events$data[[events$event_id_column]][selected][outside]
      .study_refusal(
        sprintf("Partition `%s` events extend beyond observation coverage.",
          partition),
        capability = "timing_resolved",
        reasons = c(
          sprintf("Observation coverage is [%s, %s] %s.",
            format(lower), format(upper), index$units),
          paste0("Out-of-coverage events: ", paste(event_ids, collapse = ", "),
            ".")
        ),
        remedies = paste0(
          "Correct event timing or declare an observation window covering ",
          "the complete event intervals."
        )
      )
    }
    coverage[[partition]] <- list(
      timing_resolved = TRUE,
      event_count = as.integer(sum(selected)),
      units = index$units,
      observation_range = c(lower, upper),
      event_range = c(min(onset), max(ending))
    )
  }
  list(
    resolved = all(vapply(coverage, function(value) {
      isTRUE(value$timing_resolved)
    }, logical(1))),
    coverage = coverage
  )
}

#' Bind observations, events, confounds, clocks, and partition axes
#'
#' `study()` records factual correspondence only. It does not assign model
#' roles, choose a sampling unit, infer independence, or select a
#' `generalizes_over` axis.
#'
#' @param observations An [observations()] object.
#' @param events Optional [observation_events()] record.
#' @param confounds Optional [observation_confounds()] record.
#' @param hierarchy Optional [partition_hierarchy()]. A leaf-only hierarchy is
#'   constructed from observation partitions when omitted.
#' @param clock_tolerance Nonnegative finite tolerance for clock coverage.
#' @param provenance Portable binding provenance.
#' @return An `effect_study`: a list with the validated `$observations`,
#'   `$events`, `$confounds` and `$hierarchy`, the ordered `$partitions`, a
#'   per-partition `$lineage` data frame recording which observations censoring
#'   retained, `$clock_coverage`, `$clock_tolerance`, `$capabilities`,
#'   `$provenance`, and a `$study_id`.
#' @family studies and effect maps
#' @seealso [observations()], [observation_events()],
#'   [observation_confounds()], [partition_hierarchy()] for the inputs;
#'   [study_capabilities()] to inspect what was established; and
#'   [plan_relation()] for the next step.
#' @examples
#' set.seed(1)
#' domain <- abstract_domain(4L, id = "study-example")
#' indexes <- list(
#'   `run-1` = observation_index(
#'     1:6, "run-1", time = seq(0, 10, by = 2), units = "seconds"
#'   ),
#'   `run-2` = observation_index(
#'     1:6, "run-2", time = seq(0, 10, by = 2), units = "seconds"
#'   )
#' )
#' record <- observations(
#'   lapply(indexes, function(index) matrix(rnorm(24), 6L, 4L)), indexes, domain
#' )
#' events <- observation_events(data.frame(
#'   partition = rep(names(indexes), each = 2L),
#'   event_id = paste0("e", 1:4), onset = c(0, 6, 0, 6), duration = 0.5,
#'   condition = c("face", "body", "face", "body")
#' ))
#'
#' # Binding checks correspondence only; it assigns no model roles.
#' facts <- study(record, events)
#' facts
#' study_capabilities(facts)
#'
#' # An event outside the acquired window is refused, not silently clipped.
#' late <- observation_events(data.frame(
#'   partition = "run-1", event_id = "e5", onset = 40, duration = 1
#' ))
#' catch_refusal(study(record, late))$capability
#' @export
study <- function(observations, events = NULL, confounds = NULL,
                  hierarchy = NULL, clock_tolerance = 0,
                  provenance = list()) {
  if (missing(observations)) {
    stop(paste0(
      "`observations` is required: pass the `observations()` record whose ",
      "partitions the events, confounds, and hierarchy are checked against."
    ), call. = FALSE)
  }
  observations <- .validate_observations(observations)
  if (!is.null(events)) events <- .validate_events(events)
  if (!is.null(confounds)) {
    confounds <- .validate_observation_confounds(confounds)
  }
  if (is.null(hierarchy)) {
    hierarchy <- .default_partition_hierarchy(observations$partitions)
  } else {
    hierarchy <- .validate_partition_hierarchy(hierarchy)
  }
  if (!is.numeric(clock_tolerance) || length(clock_tolerance) != 1L ||
      is.na(clock_tolerance) || !is.finite(clock_tolerance) ||
      clock_tolerance < 0) {
    stop("`clock_tolerance` must be one nonnegative finite number.",
      call. = FALSE)
  }
  .validate_effect_provenance(provenance, "study provenance")

  leaf <- hierarchy$data[[hierarchy$leaf]]
  if (!identical(leaf, observations$partitions)) {
    .study_refusal(
      "The hierarchy leaf axis does not equal the ordered observation partitions.",
      capability = "aligned_observations",
      reasons = c(
        paste0("Observation partitions: ",
          paste(observations$partitions, collapse = ", "), "."),
        paste0("Hierarchy leaves: ", paste(leaf, collapse = ", "), ".")
      ),
      remedies = paste0(
        "Provide one hierarchy row per observation partition in the same ",
        "declared order."
      )
    )
  }
  known <- observations$partitions
  unknown_events <- setdiff(unique(.event_partitions(events)), known)
  if (length(unknown_events)) {
    .study_refusal(
      "The event record names partitions absent from observations.",
      capability = "aligned_observations",
      reasons = paste0("Unknown event partitions: ",
        paste(unknown_events, collapse = ", "), "."),
      remedies = "Correct event partition labels or supply matching observations."
    )
  }
  if (!is.null(confounds)) {
    confound_partitions <- as.character(
      confounds$data[[confounds$partition_column]]
    )
    unknown_confounds <- setdiff(unique(confound_partitions), known)
    if (length(unknown_confounds)) {
      .study_refusal(
        "The confound record names partitions absent from observations.",
        capability = "aligned_observations",
        reasons = paste0("Unknown confound partitions: ",
          paste(unknown_confounds, collapse = ", "), "."),
        remedies = paste0(
          "Correct confound partition labels or supply matching observations."
        )
      )
    }
  }

  lineage <- .bind_confounds(observations, confounds)
  clocks <- .bind_event_clocks(observations, events, clock_tolerance)
  capabilities <- list(
    aligned_observations = TRUE,
    timing_resolved = clocks$resolved,
    partition_hierarchy = TRUE,
    stable_source_revision = all(vapply(
      observations$capabilities,
      function(value) .strong_sha256(value$stable_revision),
      logical(1)
    ))
  )
  semantic <- list(
    schema_version = 1L,
    observations_id = observations$observations_id,
    events_id = if (is.null(events)) NULL else events$events_id,
    confounds_id = if (is.null(confounds)) NULL else confounds$confounds_id,
    hierarchy_signature = hierarchy$signature,
    lineage = lineage,
    clock_coverage = clocks$coverage,
    clock_tolerance = as.numeric(clock_tolerance),
    provenance = provenance
  )
  structure(list(
    observations = observations,
    events = events,
    confounds = confounds,
    hierarchy = hierarchy,
    partitions = observations$partitions,
    lineage = lineage,
    clock_coverage = clocks$coverage,
    clock_tolerance = as.numeric(clock_tolerance),
    capabilities = capabilities,
    provenance = provenance,
    study_id = .sha256_signature(semantic, "study-sha256:")
  ), class = "effect_study")
}

.validate_study <- function(value, deep = TRUE) {
  expected <- c(
    "observations", "events", "confounds", "hierarchy", "partitions",
    "lineage", "clock_coverage", "clock_tolerance", "capabilities",
    "provenance", "study_id"
  )
  if (!inherits(value, "effect_study") || !is.list(value) ||
      !identical(names(value), expected)) {
    stop("Study fields are missing or noncanonical.", call. = FALSE)
  }
  observations <- .validate_observations(value$observations, deep = deep)
  events <- if (is.null(value$events)) NULL else .validate_events(value$events)
  confounds <- if (is.null(value$confounds)) {
    NULL
  } else {
    .validate_observation_confounds(value$confounds)
  }
  hierarchy <- .validate_partition_hierarchy(value$hierarchy)
  if (!identical(value$partitions, observations$partitions) ||
      !identical(hierarchy$data[[hierarchy$leaf]], value$partitions)) {
    stop("Study partition metadata are inconsistent.", call. = FALSE)
  }
  expected_lineage <- .bind_confounds(observations, confounds)
  expected_clocks <- .bind_event_clocks(
    observations, events, value$clock_tolerance
  )
  expected_capabilities <- list(
    aligned_observations = TRUE,
    timing_resolved = expected_clocks$resolved,
    partition_hierarchy = TRUE,
    stable_source_revision = all(vapply(
      observations$capabilities,
      function(item) .strong_sha256(item$stable_revision),
      logical(1)
    ))
  )
  if (!identical(value$lineage, expected_lineage) ||
      !identical(value$clock_coverage, expected_clocks$coverage) ||
      !identical(value$capabilities, expected_capabilities)) {
    stop("Study lineage or capabilities are inconsistent.", call. = FALSE)
  }
  .validate_effect_provenance(value$provenance, "study provenance")
  semantic <- list(
    schema_version = 1L,
    observations_id = observations$observations_id,
    events_id = if (is.null(events)) NULL else events$events_id,
    confounds_id = if (is.null(confounds)) NULL else confounds$confounds_id,
    hierarchy_signature = hierarchy$signature,
    lineage = expected_lineage,
    clock_coverage = expected_clocks$coverage,
    clock_tolerance = value$clock_tolerance,
    provenance = value$provenance
  )
  if (!identical(value$study_id, .sha256_signature(semantic, "study-sha256:"))) {
    stop("Study identity is inconsistent.", call. = FALSE)
  }
  value
}

#' Select a typed study axis
#'
#' This accessor supplies vocabulary only. The returned object does not select
#' a generalization target or assert independence.
#'
#' @param x An [study()].
#' @param name One hierarchy axis name.
#' @return An `effect_study_axis`: a list with the axis `$name`, its distinct
#'   `$levels`, the enclosing `$parent` axis (`NULL` at the outermost axis),
#'   and a `$signature` tied to the study hierarchy.
#' @family studies and effect maps
#' @seealso [partition_hierarchy()] which declares the available axes, and
#'   [cross_partitions()], where a generalization axis is actually chosen.
#' @examples
#' set.seed(1)
#' domain <- abstract_domain(4L, id = "study-axis-example")
#' partitions <- c("run-1", "run-2")
#' indexes <- lapply(partitions, function(partition) {
#'   observation_index(1:6, partition)
#' })
#' names(indexes) <- partitions
#' record <- observations(
#'   lapply(indexes, function(index) matrix(rnorm(24), 6L, 4L)), indexes, domain
#' )
#' facts <- study(record, hierarchy = partition_hierarchy(data.frame(
#'   partition = partitions, session = c("ses-1", "ses-2")
#' )))
#'
#' # Vocabulary only: naming an axis does not assert independence across it.
#' axis <- study_axis(facts, "session")
#' axis
#' axis$levels
#'
#' # An axis that was never declared is refused rather than invented.
#' catch_refusal(study_axis(facts, "subject"))$capability
#' @export
study_axis <- function(x, name) {
  x <- .validate_study(x, deep = FALSE)
  name <- .validate_nonempty_id(name, "name")
  if (!name %in% x$hierarchy$axes) {
    .study_refusal(
      sprintf("Study axis `%s` is not declared.", name),
      capability = "declared_generalization_axis",
      reasons = paste0("Available axes: ",
        paste(x$hierarchy$axes, collapse = ", "), "."),
      remedies = "Choose a declared axis or correct the partition hierarchy."
    )
  }
  position <- match(name, x$hierarchy$axes)
  parent <- if (position < length(x$hierarchy$axes)) {
    x$hierarchy$axes[[position + 1L]]
  } else {
    NULL
  }
  semantic <- list(
    schema_version = 1L,
    hierarchy_signature = x$hierarchy$signature,
    name = name,
    levels = x$hierarchy$levels[[name]],
    parent = parent
  )
  structure(c(semantic[-1L], list(
    signature = .sha256_signature(semantic, "study-axis-sha256:")
  )), class = "effect_study_axis")
}

#' Inspect factual study capabilities
#'
#' Reports what the binding actually established, so downstream calls can
#' require a guarantee instead of assuming it.
#'
#' @param x An [study()].
#' @return A one-row data frame with logical columns `aligned_observations`,
#'   `timing_resolved`, `partition_hierarchy`, and `stable_source_revision`.
#' @family studies and effect maps
#' @seealso [study()], and [compiler_conformance()] for the equivalent report
#'   on a compiled relation plan.
#' @examples
#' set.seed(1)
#' domain <- abstract_domain(4L, id = "study-capabilities-example")
#' index <- observation_index(
#'   1:6, "run-1", time = seq(0, 10, by = 2), units = "seconds"
#' )
#' record <- observations(
#'   list(`run-1` = matrix(rnorm(24), 6L, 4L)), list(`run-1` = index), domain
#' )
#'
#' # With no event record there is no timing to resolve, and the report says so
#' # rather than defaulting to TRUE.
#' study_capabilities(study(record))
#'
#' # Binding timed events on the same clock earns `timing_resolved`.
#' events <- observation_events(data.frame(
#'   partition = "run-1", event_id = c("e1", "e2"), onset = c(0, 6),
#'   duration = 0.5
#' ))
#' study_capabilities(study(record, events))$timing_resolved
#' @export
study_capabilities <- function(x) {
  x <- .validate_study(x, deep = FALSE)
  as.data.frame(x$capabilities, check.names = FALSE)
}

#' @export
format.effect_study <- function(x, ...) {
  x <- .validate_study(x, deep = FALSE)
  sprintf(
    "study<%d partitions; timing %s>",
    length(x$partitions),
    if (isTRUE(x$capabilities$timing_resolved)) "resolved" else "unresolved"
  )
}

#' @export
print.effect_study <- function(x, ...) {
  cat(format(x), "\n", sep = "")
  invisible(x)
}

#' @export
format.effect_study_axis <- function(x, ...) {
  sprintf("study_axis<%s; %d levels>", x$name, length(x$levels))
}

#' @export
print.effect_study_axis <- function(x, ...) {
  cat(format(x), "\n", sep = "")
  invisible(x)
}
