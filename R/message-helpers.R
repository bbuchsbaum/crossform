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
    return(sprintf("a %s vector of length %d", primary, length(x)))
  }
  sprintf("an object of class `%s`", paste(class(x), collapse = "/"))
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
    stop(sprintf(paste0(
      "`%s` must be one measurement index in 1..%d; received %s."
    ), argument, measurements, .msg_value(value)), call. = FALSE)
  }
  if (value < 1L || value > measurements) {
    stop(sprintf(paste0(
      "`%s` = %s is outside the %s's 1..%d measurements."
    ), argument, format(value, scientific = FALSE), subject, measurements),
      call. = FALSE)
  }
  as.integer(value)
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
