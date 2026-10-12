# Consistent voice for user-facing validation messages ----------------------
#
# A refusal or a shape error is only useful when the reader can act on it, so
# every user-facing message in the public entry points states what was
# received, what was expected, and — when one exists — the call that fixes it.
# These helpers keep that voice uniform without adding a dependency: declarative
# sentences, backticks around argument and function names, no exclamation
# marks, and identifiers quoted so an empty or whitespace label is visible.

# Quote a set of identifiers for a message, truncating long sets so an error
# raised on a thousand-effect space stays readable.
.msg_names <- function(x, max_shown = 8L) {
  x <- as.character(x)
  if (!length(x)) return("none")
  shown <- x[seq_len(min(length(x), max_shown))]
  out <- paste0("`", shown, "`", collapse = ", ")
  if (length(x) > max_shown) {
    out <- paste0(out, ", and ", length(x) - max_shown, " more")
  }
  out
}

# Show a content hash as its algorithm plus first `chars` hex digits. A whole
# digest is never printed: the point of quoting one in an error is to let the
# reader see that two identities differ, and twelve digits does that.
.msg_signature <- function(x, chars = 12L) {
  # Spelled out rather than delegated to `.is_string()`: this file is called
  # by `check.R`, so calling back into it would make layer 1 mutually
  # recursive for the sake of one clause.
  if (is.null(x) || !is.character(x) || length(x) != 1L || is.na(x) ||
      !nzchar(x)) {
    return("none")
  }
  if (grepl("^[A-Za-z0-9_-]+:[[:xdigit:]]{16,}$", x)) {
    algorithm <- sub(":.*$", "", x)
    digits <- sub("^[A-Za-z0-9_-]+:", "", x)
    if (nchar(digits) > chars) {
      return(paste0(algorithm, ":", substr(digits, 1L, chars), "..."))
    }
    return(x)
  }
  if (nchar(x) > chars * 2L) paste0(substr(x, 1L, chars * 2L), "...") else x
}

# "1 effect" / "3 effects", so messages do not read "1 effects".
.msg_count <- function(n, singular, plural = paste0(singular, "s")) {
  n <- as.integer(n)
  paste0(n, " ", if (identical(n, 1L)) singular else plural)
}

# A short, honest description of an unexpected value. Class comes first
# because that is what the reader needs in order to find the wrong variable.
.msg_value <- function(x) {
  if (is.null(x)) return("NULL")
  if (is.matrix(x)) {
    return(sprintf("a %d x %d %s matrix", nrow(x), ncol(x), typeof(x)))
  }
  if (is.array(x)) {
    return(sprintf("a %s array with dimensions %s", typeof(x),
      paste(dim(x), collapse = " x ")))
  }
  if (is.data.frame(x)) {
    return(sprintf("a data frame with %s and %s",
      .msg_count(nrow(x), "row"), .msg_count(ncol(x), "column")))
  }
  if (is.function(x)) return("a function")
  primary <- class(x)[[1L]]
  if (identical(primary, "numeric")) primary <- "numeric"
  if (primary %in% c("numeric", "integer", "character", "logical", "double",
      "complex", "list", "factor")) {
    article <- if (primary %in% c("integer")) "an" else "a"
    return(sprintf("%s %s vector of length %d", article, primary, length(x)))
  }
  sprintf("an object of class `%s`", paste(class(x), collapse = "/"))
}

# The same job as `.msg_value()`, but for the scalar case a guard usually
# fails on: when one value was wanted and one value arrived, naming its class
# and length says nothing the reader did not already know, so show the value
# itself. `.msg_value()` remains the description for everything else, and is
# left untouched because existing messages quote it verbatim.
.msg_received <- function(x) {
  if (is.null(x)) return("NULL")
  if (!is.atomic(x) || length(x) != 1L) return(.msg_value(x))
  if (is.na(x)) return(sprintf("NA (%s)", .msg_value(x)))
  if (is.character(x)) {
    if (!nzchar(x)) return("an empty string")
    shown <- if (nchar(x) > 40L) paste0(substr(x, 1L, 40L), "...") else x
    return(sprintf("\"%s\"", shown))
  }
  if (is.logical(x)) return(if (x) "TRUE" else "FALSE")
  if (is.numeric(x)) return(format(x, scientific = FALSE))
  .msg_value(x)
}

# Measurement selectors: identifiers versus positions -------------------------
#
# A frame's measurements have two names. Their *position* is the row of
# `$weights` (and of every view's values), 1..n. Their *identifier* is
# `$index$measurement`, which is what every view prints in its `measurement`
# column: domain feature ids for voxelwise and searchlight frames (full-volume
# indices on a masked `volume_domain()`), region labels, or `"family::node"`.
# On an abstract domain or a full mask the two coincide; on a masked volume
# they do not, and a bare `117` copied from a print would silently select the
# 117th row instead. So a bare number is read as a position only where it
# cannot be confused with a printed identifier, and is refused otherwise; the
# two wrappers below say which one is meant.

#' Select measurements by identifier or by position
#'
#' Wrap the measurements you pass to `at =` in [rdm_sampling_covariance()],
#' or to `measurement =` and `highlight =` in [plot_views], to say whether
#' the numbers are *identifiers* or *row positions*.
#'
#' Every compiled frame names its measurements twice. The **identifier** is
#' the frame's `$index$measurement`, carried into each view as `$index` and
#' printed in the view's `measurement` column: domain feature ids for
#' [voxelwise()] and [searchlights()] (the full-volume voxel indices of a
#' masked [volume_domain()]), region labels for [regions()], and
#' `"family::node"` for a [frame_family()]. The **position** is the row of
#' the frame's `$weights` and of the view's values, `1..n`: what
#' `which.max(view$total)` or `order(view$total)` returns.
#'
#' On an [abstract_domain()] with default ids, or a volume without holes, the
#' two coincide and a bare number works either way. When the identifiers are
#' whole numbers that differ from the positions -- any masked volume -- a bare
#' number is ambiguous: `117` copied from a printed view is voxel 117, but as
#' a position it is some other voxel. crossform does not guess. It refuses a
#' bare number there, and asks for one of:
#'
#' - `measurement_id(117)`, or the character identifier `"117"`: the
#'   measurement printed as `117`.
#' - `measurement_position(which.max(view$total))`: the row position.
#'
#' Identifiers that are not whole numbers (region names, `"family::node"`)
#' cannot be mistaken for positions, so bare numbers stay positions there.
#'
#' @param x For `measurement_id()`, one or more measurement identifiers,
#'   numeric or character, exactly as a view's `measurement` column prints
#'   them (equivalently, values of the frame's `$index$measurement`). For
#'   `measurement_position()`, one or more whole-number row positions.
#' @return A classed selector, `crossform_measurement_id` or
#'   `crossform_measurement_position`, holding `x`. It is resolved against
#'   the frame or view it is passed with; an identifier the frame does not
#'   carry, or a position outside `1..n`, is an input error.
#' @family sampling uncertainty
#' @examples
#' mask <- array(TRUE, c(4L, 4L, 2L))
#' mask[1:2, , 1] <- FALSE
#' domain <- volume_domain(mask, spacing = c(3, 3, 3), id = "masked")
#' frame <- compile_frame(searchlights(3.1), domain)
#'
#' # Identifiers are full-volume voxel indices; positions are rows 1..n.
#' head(frame$index$measurement)
#'
#' # The same measurement, named both ways.
#' measurement_id(11)
#' measurement_position(match(11, frame$index$measurement))
#' @export
measurement_id <- function(x) {
  if (!(is.numeric(x) || is.character(x)) || is.matrix(x) || !length(x) ||
      anyNA(x) || (is.character(x) && any(!nzchar(x)))) {
    .input_error(sprintf(paste0(
      "`measurement_id()` takes one or more measurement identifiers as a ",
      "numeric or character vector without missing values; received %s."
    ), .msg_value(x)),
      arg = "x", received = .msg_value(x),
      expected = "measurement identifiers, as a view prints them")
  }
  structure(list(values = unclass(x)), class = "crossform_measurement_id")
}

#' @rdname measurement_id
#' @export
measurement_position <- function(x) {
  if (!is.numeric(x) || is.matrix(x) || !length(x) || anyNA(x) ||
      any(!is.finite(x)) || any(x %% 1 != 0)) {
    .input_error(sprintf(paste0(
      "`measurement_position()` takes one or more whole-number row ",
      "positions; received %s."
    ), .msg_value(x)),
      arg = "x", received = .msg_value(x),
      expected = "whole-number row positions")
  }
  structure(list(values = as.numeric(unclass(x))),
    class = "crossform_measurement_position")
}

#' @export
print.crossform_measurement_id <- function(x, ...) {
  cat("<measurement_id> ", paste(format(x$values, scientific = FALSE,
    trim = TRUE), collapse = ", "), "\n", sep = "")
  invisible(x)
}

#' @export
print.crossform_measurement_position <- function(x, ...) {
  cat("<measurement_position> ", paste(format(x$values, scientific = FALSE,
    trim = TRUE), collapse = ", "), "\n", sep = "")
  invisible(x)
}

# The identifiers a frame prints for its measurements, or NULL when it carries
# none (positions are then the only names).
.msg_frame_measurement_ids <- function(frame) {
  if (is.null(frame)) return(NULL)
  index <- frame$index
  if (is.data.frame(index) && "measurement" %in% names(index)) {
    return(index$measurement)
  }
  if (is.atomic(index) && length(index)) return(index)
  NULL
}

# One comparable text key per identifier, so `117`, `117L` and `"117"` name
# the same measurement and large voxel indices never print as `1e+05`.
.msg_id_key <- function(x) {
  if (is.factor(x)) x <- as.character(x)
  if (is.numeric(x)) {
    whole <- is.finite(x) & x %% 1 == 0
    out <- as.character(x)
    out[whole] <- sprintf("%.0f", x[whole])
    return(out)
  }
  trimws(as.character(x))
}

# Whether a bare number could be read as a printed identifier of a different
# measurement: the identifiers are all whole numbers and are not 1..n in order.
.msg_ids_shadow_positions <- function(ids, n) {
  if (is.null(ids) || length(ids) != n) return(FALSE)
  text <- .msg_id_key(ids)
  if (!all(grepl("^-?[0-9]+$", text))) return(FALSE)
  !identical(text, as.character(seq_len(n)))
}

# Resolve a measurement selector to row positions. Character values and
# `measurement_id()` are identifiers; `measurement_position()` is positions;
# a bare number is a position unless the identifiers shadow positions, where
# it is refused as ambiguous. `ids = NULL` keeps the historical positional
# contract for internal callers that hold no identifiers.
.msg_resolve_measurements <- function(value, ids, measurements,
                                      argument = "at", subject = "plan") {
  measurements <- as.integer(measurements)
  if (inherits(value, "crossform_measurement_position")) {
    return(.msg_measurement_positions(value$values, measurements, argument,
      subject))
  }
  if (inherits(value, "crossform_measurement_id") || is.character(value)) {
    requested <- if (is.character(value)) value else value$values
    if (!length(requested) || anyNA(requested)) {
      .input_error(sprintf(paste0(
        "`%s` must name one or more measurement identifiers; received %s."
      ), argument, .msg_value(requested)),
        arg = argument, received = .msg_value(requested),
        expected = "measurement identifiers")
    }
    if (is.null(ids)) ids <- seq_len(measurements)
    # A repeat is named here, by identifier, so the caller sees the value they
    # typed rather than the row position it resolved to.
    if (anyDuplicated(.msg_id_key(requested))) {
      repeated <- .msg_id_key(requested)
      .input_error(sprintf(paste0(
        "`%s` must name each measurement once; %s repeated."
      ), argument, .msg_names(unique(repeated[duplicated(repeated)]))),
        arg = argument, received = "repeated measurement identifiers",
        expected = "unique measurement identifiers")
    }
    positions <- match(.msg_id_key(requested), .msg_id_key(ids))
    if (anyNA(positions)) {
      .input_error(sprintf(paste0(
        "`%s` names %s that the %s does not carry: %s. Identifiers are the ",
        "values in the view's `measurement` column (the frame's ",
        "`$index$measurement`), for example %s."
      ), argument,
        .msg_count(sum(is.na(positions)), "measurement identifier"),
        subject, .msg_names(as.character(requested[is.na(positions)])),
        .msg_names(as.character(utils::head(ids, 3L)))),
        arg = argument,
        received = .msg_names(as.character(requested[is.na(positions)])),
        expected = "identifiers from the view's `measurement` column")
    }
    return(as.integer(positions))
  }
  if (is.numeric(value) && .msg_ids_shadow_positions(ids, measurements)) {
    .msg_ambiguous_measurement(value, ids, measurements, argument, subject)
  }
  .msg_measurement_positions(value, measurements, argument, subject)
}

# The refusal for a bare number where identifiers shadow positions. It names
# both spellings with the caller's own value, so either can be pasted back.
.msg_ambiguous_measurement <- function(value, ids, measurements,
                                       argument = "at", subject = "plan") {
  shown <- if (length(value) == 1L) {
    format(value, scientific = FALSE)
  } else {
    sprintf("c(%s%s)", paste(format(utils::head(value, 4L),
      scientific = FALSE, trim = TRUE), collapse = ", "),
      if (length(value) > 4L) ", ..." else "")
  }
  .input_error(sprintf(paste0(
    "`%s` = %s is ambiguous: this %s's measurements are identified by ",
    "whole numbers (%s, ...) that are not their row positions 1..%d, so a ",
    "bare number could mean either. Say which: `%s = measurement_id(%s)` ",
    "(or the character id) selects the measurement printed with that ",
    "identifier in a view's `measurement` column; `%s = ",
    "measurement_position(%s)` selects that row, as returned by ",
    "`which.max()` or `order()` on a view's values."
  ), argument, shown, subject,
    paste(utils::head(.msg_id_key(ids), 3L), collapse = ", "),
    as.integer(measurements), argument, shown, argument, shown),
    arg = argument, received = sprintf("the bare number %s", shown),
    expected = "`measurement_id()` or `measurement_position()`")
}

.msg_measurement_positions <- function(value, measurements, argument,
                                       subject) {
  if (length(value) == 1L) {
    .msg_measurement_index(value, measurements, argument, subject)
  } else {
    .msg_measurement_indices(value, measurements, argument, subject)
  }
}

# One measurement position, checked against the number of measurements the
# caller's plan or frame actually has. The message names the argument, the
# value received, and the range, because "out of bounds" alone leaves the
# reader guessing which of the two numbers is wrong.
.msg_measurement_index <- function(value, measurements, argument = "at",
                                   subject = "plan") {
  measurements <- as.integer(measurements)
  usable <- is.numeric(value) && length(value) == 1L && !is.na(value) &&
    is.finite(value) && value %% 1 == 0
  if (!usable) {
    .input_error(sprintf(paste0(
      "`%s` must be one measurement index in 1..%d; received %s."
    ), argument, measurements, .msg_value(value)))
  }
  if (value < 1L || value > measurements) {
    .input_error(sprintf(paste0(
      "`%s` = %s is outside the %s's 1..%d measurements."
    ), argument, format(value, scientific = FALSE), subject, measurements))
  }
  as.integer(value)
}

# One or more measurement positions. A length-1 vector is the historical
# single-node contract; longer vectors name an explicit batch.
.msg_measurement_indices <- function(value, measurements, argument = "at",
                                     subject = "plan", ids = NULL) {
  # Identifier-aware callers pass the frame's identifiers; selectors and
  # character ids are resolved there, and a bare number is refused where it
  # could be mistaken for a printed identifier.
  if (!is.null(ids) || is.character(value) ||
      inherits(value, c("crossform_measurement_id",
        "crossform_measurement_position"))) {
    return(.msg_resolve_measurements(value, ids, measurements, argument,
      subject))
  }
  # A single index keeps the scalar message, which names the offending value.
  if (length(value) == 1L) {
    return(.msg_measurement_index(value, measurements, argument, subject))
  }
  measurements <- as.integer(measurements)
  usable <- is.numeric(value) && length(value) >= 1L &&
    all(is.finite(value)) && all(value %% 1 == 0)
  if (!usable) {
    .input_error(sprintf(paste0(
      "`%s` must be one or more measurement indices in 1..%d; received %s."
    ), argument, measurements, .msg_value(value)))
  }
  value <- as.integer(value)
  bad <- value < 1L | value > measurements
  if (any(bad)) {
    .input_error(sprintf(paste0(
      "`%s` contains indices outside the %s's 1..%d measurements: %s."
    ), argument, subject, measurements, .msg_positions(bad)))
  }
  value
}

# Which positions of a numeric vector are not usable, for a message that can
# point at the offending entries rather than restating the rule.
.msg_positions <- function(index, labels = NULL, max_shown = 6L) {
  index <- which(index)
  if (!length(index)) return("none")
  if (!is.null(labels) && length(labels) >= max(index)) {
    return(.msg_names(labels[index], max_shown))
  }
  shown <- index[seq_len(min(length(index), max_shown))]
  out <- paste(shown, collapse = ", ")
  if (length(index) > max_shown) {
    out <- paste0(out, ", and ", length(index) - max_shown, " more")
  }
  out
}
