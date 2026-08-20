# BIDS file adapters -------------------------------------------------------

.bids_partition_files <- function(files, partitions, what) {
  if (!.is_strings(files) || length(files) < 1L) {
    .input_error(sprintf("`%s` must contain one or more file paths.", what))
  }
  if (is.null(partitions)) partitions <- names(files)
  # The adapter's own argument check, spelled out rather than borrowed from
  # the study-facts validator it used to call: what an adapter needs here is
  # an ordinary guard on a character vector, and an external adapter writes
  # exactly these three lines. The partition *column* it eventually builds is
  # checked again, for real, by `observation_events()`.
  if (!.is_strings(partitions, unique = TRUE) || length(partitions) < 1L) {
    .input_error(sprintf(
      "`partitions` must contain unique nonempty identifiers%s.",
      if (is.character(partitions) && anyDuplicated(partitions)) {
        sprintf("; %s appears more than once",
          .msg_names(unique(partitions[duplicated(partitions)])))
      } else {
        sprintf("; received %s", .msg_value(partitions))
      }))
  }
  if (length(partitions) != length(files)) {
    .input_error(sprintf("`partitions` must supply %s; received %s.",
      .msg_count(length(files), "name"),
      .msg_count(length(partitions), "name")))
  }
  partitions <- unname(partitions)
  if (!is.null(names(files)) && any(nzchar(names(files))) &&
      !identical(names(files), partitions)) {
    .input_error("Named BIDS files must follow the declared partition order.")
  }
  files <- vapply(files, normalizePath, character(1), mustWork = TRUE)
  names(files) <- partitions
  files
}

.read_bids_table <- function(path) {
  utils::read.delim(
    path, header = TRUE, sep = "\t", quote = "", comment.char = "",
    check.names = FALSE, stringsAsFactors = FALSE, na.strings = character()
  )
}

.bids_file_provenance <- function(files) {
  lapply(files, function(path) list(
    file = basename(path),
    revision = .sha256_file(path)
  ))
}

#' Import BIDS task events as typed event facts
#'
#' The adapter preserves arbitrary BIDS columns and adds private partition and
#' event-key columns required by the generic [observation_events()] contract.
#' Partition identity is explicit rather than inferred from filenames.
#'
#' @param files Character event-TSV paths, one per partition.
#' @param partitions Explicit ordered partition identifiers. Named `files` may
#'   supply these names.
#' @param units Physical onset and duration units.
#' @return An [observation_events()] fact object whose `$data` holds every
#'   original BIDS column plus the private `.bids_partition` and
#'   `.bids_event_id` key columns, and whose `$provenance` records each file's
#'   basename and content revision.
#' @seealso [observation_events()] for the generic contract,
#'   [bids_confounds()] for the confound side, and [bids_study()] to bind both
#'   into a study.
#' @examples
#' # Stand in for one run's events.tsv.
#' path <- tempfile(fileext = ".tsv")
#' utils::write.table(
#'   data.frame(
#'     onset = c(0, 6), duration = 0.5, trial_type = c("face", "body")
#'   ),
#'   path, sep = "\t", row.names = FALSE, quote = FALSE
#' )
#'
#' # Partition identity is declared by the caller, never parsed from the
#' # filename, and every original column survives the import.
#' record <- crossform:::bids_events(c(`run-1` = path))
#' record$timing
#' names(record$data)
#' record$data$.bids_event_id
#' unlink(path)
#' @keywords internal
bids_events <- function(files, partitions = names(files), units = "seconds") {
  files <- .bids_partition_files(files, partitions, "files")
  tables <- lapply(names(files), function(partition) {
    value <- .read_bids_table(files[[partition]])
    if (!all(c("onset", "duration") %in% names(value)) ||
        !is.numeric(value$onset) || !is.numeric(value$duration)) {
      .capability_refusal(
        sprintf("BIDS events for `%s` lack resolved numeric timing.", partition),
        capability = "timing_resolved",
        namespace = "bids_adapter",
        reasons = "The adapter requires numeric `onset` and `duration` columns.",
        remedies = "Resolve missing or nonnumeric event timing before import."
      )
    }
    if (any(c(".bids_partition", ".bids_event_id") %in% names(value))) {
      .input_error(
        "BIDS tables may not use crossform's private adapter columns."
      )
    }
    value$.bids_partition <- partition
    value$.bids_event_id <- sprintf("event-%06d", seq_len(nrow(value)))
    value
  })
  data <- do.call(rbind, tables)
  rownames(data) <- NULL
  observation_events(
    data,
    partition = ".bids_partition",
    event_id = ".bids_event_id",
    onset = "onset",
    duration = "duration",
    units = units,
    provenance = list(
      adapter = "BIDS task events",
      files = .bids_file_provenance(files)
    )
  )
}

#' Import fMRIPrep-style confounds as observation facts
#'
#' All columns are preserved without assigning them target or nuisance roles.
#' If `censor` is supplied it must name a complete logical retain column; no
#' censor policy is inferred from motion or outlier columns.
#'
#' @param files Character confound-TSV paths, one per partition.
#' @param partitions Explicit ordered partition identifiers.
#' @param observation_ids Optional named list of observation identifiers. Row
#'   numbers are used when omitted.
#' @param censor Optional logical retain-column name already present in every
#'   input table.
#' @return An [observation_confounds()] fact object whose `$data` holds every
#'   original confound column plus the private `.bids_partition` and
#'   `.bids_observation_id` key columns, with `$censor_column` set only when
#'   `censor` was supplied.
#' @seealso [observation_confounds()] for the generic contract,
#'   [bids_events()] for the event side, and [bids_study()] to bind both.
#' @examples
#' # Stand in for one run's fMRIPrep confounds table.
#' path <- tempfile(fileext = ".tsv")
#' utils::write.table(
#'   data.frame(
#'     framewise_displacement = c(0, 0.1, 0.9, 0.2),
#'     retained = c(TRUE, TRUE, FALSE, TRUE)
#'   ),
#'   path, sep = "\t", row.names = FALSE, quote = FALSE
#' )
#'
#' # Naming the retain column is required: no censor policy is inferred from
#' # the motion columns, however large they are.
#' record <- crossform:::bids_confounds(c(`run-1` = path), censor = "retained")
#' record$censor_column
#' sum(record$data$retained)
#'
#' # Naming a column that is absent or not logical refuses explicitly.
#' catch_refusal(
#'   crossform:::bids_confounds(
#'     c(`run-1` = path), censor = "framewise_displacement"
#'   )
#' )$capability
#' unlink(path)
#' @keywords internal
bids_confounds <- function(files, partitions = names(files),
                           observation_ids = NULL, censor = NULL) {
  files <- .bids_partition_files(files, partitions, "files")
  if (!is.null(observation_ids) &&
      (!is.list(observation_ids) ||
       !setequal(names(observation_ids), names(files)))) {
    .input_error(
      "`observation_ids` must be a named list covering every partition."
    )
  }
  tables <- lapply(names(files), function(partition) {
    value <- .read_bids_table(files[[partition]])
    if (any(c(".bids_partition", ".bids_observation_id") %in% names(value))) {
      .input_error(
        "Confound tables may not use crossform's private adapter columns."
      )
    }
    ids <- if (is.null(observation_ids)) {
      seq_len(nrow(value))
    } else {
      observation_ids[[partition]]
    }
    if (length(ids) != nrow(value)) {
      .capability_refusal(
        sprintf("Confound rows for `%s` do not cover its observation ids.",
          partition),
        capability = "aligned_observations",
        namespace = "bids_adapter",
        reasons = sprintf("Found %d rows for %d declared observations.",
          nrow(value), length(ids)),
        remedies = "Supply exactly one confound row per acquired volume."
      )
    }
    if (!is.null(censor) &&
        (!censor %in% names(value) || !is.logical(value[[censor]]))) {
      .capability_refusal(
        sprintf("Confound censor policy is unresolved for `%s`.", partition),
        capability = "censoring_declared",
        namespace = "bids_adapter",
        reasons = "The named censor column is absent or is not logical.",
        remedies = "Create an explicit logical retain column before import."
      )
    }
    value$.bids_partition <- partition
    value$.bids_observation_id <- ids
    value
  })
  data <- do.call(rbind, tables)
  rownames(data) <- NULL
  observation_confounds(
    data,
    partition = ".bids_partition",
    observation_id = ".bids_observation_id",
    censor = censor,
    provenance = list(
      adapter = "fMRIPrep-style confounds",
      files = .bids_file_provenance(files)
    )
  )
}

#' Bind BIDS files into a generic study
#'
#' BIDS is an adapter boundary, not the crossform object model: the result is
#' the same [study()] that can be built from any event and confound source.
#'
#' @param observations A declared [observations()] record.
#' @param event_files Event-TSV paths in observation partition order.
#' @param confound_files Optional confound-TSV paths in the same order.
#' @param partitions Explicit ordered partition identifiers.
#' @param observation_ids Optional confound-row identifiers. Defaults to the
#'   identifiers already declared by `observations`.
#' @param censor Optional explicit logical retain-column name.
#' @param hierarchy Optional [partition_hierarchy()].
#' @param units Physical event-time unit.
#' @return A generic [study()] object, identical in structure to one built from
#'   any other event and confound source, with `$provenance$adapter` recording
#'   the BIDS route.
#' @family typed observation facts
#' @seealso [study()] for the object returned, [bids_events()] and
#'   [bids_confounds()] for the individual file adapters, and
#'   [plan_relation()] for the next step.
#' @examples
#' set.seed(1)
#' domain <- abstract_domain(3L, id = "bids-study-example")
#' index <- observation_index(
#'   1:4, "run-1", time = seq(0, by = 2, length.out = 4L), units = "seconds"
#' )
#' record <- observations(
#'   list(`run-1` = matrix(rnorm(12), 4L, 3L)), list(`run-1` = index), domain
#' )
#'
#' path <- tempfile(fileext = ".tsv")
#' utils::write.table(
#'   data.frame(
#'     onset = c(0, 4), duration = 0.5, trial_type = c("face", "body")
#'   ),
#'   path, sep = "\t", row.names = FALSE, quote = FALSE
#' )
#'
#' # BIDS stays at the file boundary: what comes back is the ordinary study
#' # object, with the event clock checked against the acquisition clock.
#' facts <- bids_study(record, event_files = c(`run-1` = path))
#' class(facts)
#' study_capabilities(facts)$timing_resolved
#' unlink(path)
#' @export
bids_study <- function(observations, event_files, confound_files = NULL,
                       partitions = observations$partitions,
                       observation_ids = NULL, censor = NULL,
                       hierarchy = NULL, units = "seconds") {
  # Only the class is checked here. The record itself is validated by
  # `study()` below, which refuses a malformed `observations` in exactly the
  # words it always did -- one frame later, and without the adapter having to
  # reach for the validator that call already runs on intake.
  .check_class(observations, "effect_observations", "observations",
    from = "observations()")
  if (!identical(partitions, observations$partitions)) {
    .input_error(
      "BIDS `partitions` must equal the observation partition axis in order."
    )
  }
  event_record <- bids_events(event_files, partitions, units)
  confound_record <- NULL
  if (!is.null(confound_files)) {
    if (is.null(observation_ids)) {
      observation_ids <- lapply(observations$indexes, `[[`, "observation_id")
    }
    confound_record <- bids_confounds(
      confound_files, partitions, observation_ids, censor
    )
  }
  study(
    observations,
    events = event_record,
    confounds = confound_record,
    hierarchy = hierarchy,
    # Frozen because study provenance participates in downstream plan identity.
    provenance = list(adapter = "crossform::bids_study")
  )
}
