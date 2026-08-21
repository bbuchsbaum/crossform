# Opening a relation's sources ------------------------------------------------
#
# Layer 2 (values). A source session is not source vocabulary: it is what a
# *relation* does with its sources, opening one handle per distinct descriptor,
# reading feature blocks through it, and closing every handle exactly once.
# The distinction is the whole point of the value order. `R/source.R` defines
# what a source is and how one descriptor is opened; this file spends that
# vocabulary on a relation, and so it may call both `source.R` and
# `relation.R`. While these two hundred lines sat in `source.R` the one call
# they make to `.validate_relation()` ran the other way -- the file that
# defines a source reached up to the file that defines a relation over
# sources -- and that single edge was enough to hold `source.R` inside the
# nine-file value tangle.

.open_relation_source_session <- function(relation,
                                          open_descriptor = .open_source_descriptor,
                                          shared_opener = NULL,
                                          validate = TRUE) {
  if (isTRUE(validate)) .validate_relation(relation)
  if (!is.function(open_descriptor)) {
    .input_error("`open_descriptor` must be a function.")
  }
  handles <- list()
  partition_keys <- stats::setNames(vector("list", length(relation$partitions)),
    relation$partitions)
  facts <- new.env(parent = emptyenv())
  facts$closed <- FALSE
  facts$close_attempts <- 0L
  facts$read_count <- stats::setNames(integer(length(relation$partitions)),
    relation$partitions)
  facts$bytes_read <- stats::setNames(numeric(length(relation$partitions)),
    relation$partitions)

  close_opened <- function() {
    if (!facts$closed) {
      facts$close_attempts <- facts$close_attempts + 1L
      failures <- character()
      for (handle in rev(handles)) {
        failure <- tryCatch({
          .close_source_handle(handle)
          NULL
        }, error = function(error) conditionMessage(error))
        if (!is.null(failure)) failures <- c(failures, failure)
      }
      facts$closed <- TRUE
      if (length(failures)) {
        .input_error(paste("Source-session cleanup failed:",
          paste(failures, collapse = "; ")))
      }
    }
    invisible(NULL)
  }

  for (partition in relation$partitions) {
    descriptor <- relation$sources[[partition]]$descriptor
    if (is.null(descriptor) || identical(descriptor$access, "coordinator")) {
      partition_keys[partition] <- list(NULL)
      next
    }
    key <- .source_descriptor_key(descriptor)
    partition_keys[[partition]] <- key
    if (is.null(handles[[key]])) {
      handle <- tryCatch(
        open_descriptor(descriptor,
          expected_revision = relation$capabilities[[partition]]$stable_revision,
          shared_opener = shared_opener),
        error = function(error) {
          cleanup_error <- tryCatch({
            close_opened()
            NULL
          }, error = function(cleanup) conditionMessage(cleanup))
          if (!is.null(cleanup_error)) {
            .input_error(
              paste(conditionMessage(error), cleanup_error, sep = "; ")
            )
          }
          stop(error)
        }
      )
      handles[[key]] <- handle
    }
  }

  read <- function(partition, features) {
    if (facts$closed) .input_error("Source session is closed.")
    if (!partition %in% relation$partitions) {
      .input_error("Source session partition is invalid.")
    }
    features <- .validate_source_features(features, relation$n_features)
    key <- partition_keys[[partition]]
    value <- if (is.null(key)) {
      relation$sources[[partition]]$read(features)
    } else {
      handles[[key]]$read(features)
    }
    facts$read_count[[partition]] <- facts$read_count[[partition]] + 1L
    facts$bytes_read[[partition]] <- facts$bytes_read[[partition]] +
      prod(as.double(dim(value))) * 8
    value
  }

  summary <- function() {
    descriptors <- lapply(relation$sources, `[[`, "descriptor")
    access <- vapply(descriptors, function(descriptor) {
      if (is.null(descriptor)) "opaque_coordinator" else descriptor$access
    }, character(1))
    list(
      access_mode = stats::setNames(access, relation$partitions),
      distinct_owned_handles = length(handles),
      read_count = facts$read_count,
      bytes_read = facts$bytes_read,
      closed = facts$closed,
      close_attempts = facts$close_attempts
    )
  }

  structure(list(read = read, close = close_opened, summary = summary),
    class = "effect_source_session")
}

.close_source_session <- function(session) {
  if (!inherits(session, "effect_source_session") || !is.function(session$close)) {
    .input_error("`session` must be a crossform source session.")
  }
  session$close()
}

# One two-sided source session, shared by the compiled effect task below and
# by the raw measurement task in R/measurement-kernel.R.
#
# The two callers differ in exactly three things, which are the three
# parameters here: the class the session carries, the noun in its side error,
# and whether the relations are revalidated on the way in. Everything else --
# the shared-relation shortcut, the guarantee that a failure opening the right
# side closes the left one, the idempotent two-sided close, the read
# dispatcher -- is the same session, and was written out twice.
#
# The task itself arrives already validated. This is a layer-2 file and both
# task types are plan-layer values: `.validate_evidence_task()` runs in the
# plan validators, and the compiled task is a derived projection of a
# validated evidence task, so re-checking it here would be a values file
# auditing a plan.
.open_two_sided_source_session <- function(task, class, side_noun,
                                           open_descriptor, shared_opener,
                                           validate = TRUE) {
  open_side <- function(relation) {
    .open_relation_source_session(relation,
      open_descriptor = open_descriptor, shared_opener = shared_opener,
      validate = validate)
  }
  reject_side <- function() {
    .input_error(sprintf(
      "%s source side must be `left` or `right`.", side_noun
    ))
  }
  if (isTRUE(task$same_relation)) {
    session <- open_side(task$left_relation)
    return(structure(
      list(
        read = function(side, partition, features) {
          if (!side %in% c("left", "right")) reject_side()
          session$read(partition, features)
        },
        close = session$close,
        summary = session$summary
      ),
      class = class
    ))
  }
  left <- open_side(task$left_relation)
  right <- tryCatch(
    open_side(task$right_relation),
    error = function(error) {
      left$close()
      stop(error)
    }
  )
  closed <- FALSE
  close_both <- function() {
    if (!closed) {
      on.exit(right$close(), add = TRUE)
      left$close()
      closed <<- TRUE
    }
    invisible(NULL)
  }
  structure(
    list(
      read = function(side, partition, features) {
        if (identical(side, "left")) {
          left$read(partition, features)
        } else if (identical(side, "right")) {
          right$read(partition, features)
        } else {
          reject_side()
        }
      },
      close = close_both,
      summary = function() list(left = left$summary(), right = right$summary())
    ),
    class = class
  )
}

.open_effect_task_source_session <- function(task,
                                             open_descriptor = .open_source_descriptor,
                                             shared_opener = NULL,
                                             validate = TRUE) {
  .open_two_sided_source_session(
    task, "effect_task_source_session", "Task",
    open_descriptor = open_descriptor, shared_opener = shared_opener,
    validate = validate
  )
}
