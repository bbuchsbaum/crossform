# BIDS file adapters -------------------------------------------------------

.bids_partition_files <- function(files, partitions, what) {
  if (!is.character(files) || length(files) < 1L || anyNA(files) ||
      any(!nzchar(files))) {
    stop(sprintf("`%s` must contain one or more file paths.", what),
      call. = FALSE)
  }
  if (is.null(partitions)) partitions <- names(files)
  partitions <- .validate_partition_names(partitions, length(files),
    "partitions")
  if (!is.null(names(files)) && any(nzchar(names(files))) &&
      !identical(names(files), partitions)) {
    stop("Named BIDS files must follow the declared partition order.",
      call. = FALSE)
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
    revision = .file_sha256(path)
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
#' @return An [observation_events()] fact object.
#' @export
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
      stop("BIDS tables may not use crossform's private adapter columns.",
        call. = FALSE)
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
#' @return An [observation_confounds()] fact object.
#' @export
bids_confounds <- function(files, partitions = names(files),
                           observation_ids = NULL, censor = NULL) {
  files <- .bids_partition_files(files, partitions, "files")
  if (!is.null(observation_ids) &&
      (!is.list(observation_ids) ||
       !setequal(names(observation_ids), names(files)))) {
    stop("`observation_ids` must be a named list covering every partition.",
      call. = FALSE)
  }
  tables <- lapply(names(files), function(partition) {
    value <- .read_bids_table(files[[partition]])
    if (any(c(".bids_partition", ".bids_observation_id") %in% names(value))) {
      stop("Confound tables may not use crossform's private adapter columns.",
        call. = FALSE)
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
#' @return A generic [study()] object.
#' @export
bids_study <- function(observations, event_files, confound_files = NULL,
                       partitions = observations$partitions,
                       observation_ids = NULL, censor = NULL,
                       hierarchy = NULL, units = "seconds") {
  observations <- .validate_observations(observations)
  if (!identical(partitions, observations$partitions)) {
    stop("BIDS `partitions` must equal the observation partition axis in order.",
      call. = FALSE)
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
    provenance = list(adapter = "effectagram::bids_study")
  )
}
