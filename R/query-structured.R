# Structured rank-one pair-difference queries --------------------------------
#
# Every RDM edge is the rank-one operator (e_i - e_j)(e_i - e_j)^T, and a
# fixed linear RSA readout is a fixed map from pair space to a few output
# views. Materializing those queries as a dense packed
# `q(q+1)/2`-by-`q(q-1)/2` matrix costs O(q^4) memory before any data are
# read, and forces a dense q-by-q operator matmul per view inside every
# kernel. The structured representation stores only the pair index vectors
# plus the optional pair-space coefficient map; kernels evaluate one edge
# with two row differences and a Hadamard product:
#   tr(u u^T L K R^T) = (u^T L) K (u^T R)^T with u = e_i - e_j.

.pair_difference_query <- function(effects, pairs = NULL,
                                   coefficients = NULL, labels = NULL) {
  q <- length(effects)
  if (!is.character(effects) || q < 2L || anyNA(effects) ||
      any(!nzchar(effects)) || anyDuplicated(effects)) {
    stop("Pair-difference queries require at least two named effects.",
      call. = FALSE)
  }
  if (is.null(pairs)) {
    pairs <- t(utils::combn(seq_len(q), 2L))
  }
  if (is.character(pairs)) {
    resolved <- matrix(match(pairs, effects), nrow(pairs), ncol(pairs))
    if (anyNA(resolved)) {
      stop("Pair names must identify declared experimental effects.",
        call. = FALSE)
    }
    pairs <- resolved
  }
  if (!is.matrix(pairs) || ncol(pairs) != 2L || nrow(pairs) < 1L ||
      !is.numeric(pairs) || anyNA(pairs) || any(pairs %% 1 != 0) ||
      any(pairs < 1L) || any(pairs > q)) {
    stop("`pairs` must be a two-column matrix of effect indices or names.",
      call. = FALSE)
  }
  storage.mode(pairs) <- "integer"
  if (any(pairs[, 1L] == pairs[, 2L])) {
    stop("Pair-difference queries require two distinct effects per pair.",
      call. = FALSE)
  }
  swapped <- pairs[, 1L] > pairs[, 2L]
  if (any(swapped)) {
    pairs[swapped, ] <- pairs[swapped, c(2L, 1L), drop = FALSE]
  }
  if (anyDuplicated(paste(pairs[, 1L], pairs[, 2L], sep = "\r"))) {
    stop("`pairs` contains duplicate effect pairs.", call. = FALSE)
  }
  if (is.null(labels)) {
    labels <- paste(effects[pairs[, 1L]], "-", effects[pairs[, 2L]])
  }
  if (!is.character(labels) || length(labels) != nrow(pairs) ||
      anyNA(labels)) {
    stop("Pair labels must name every pair.", call. = FALSE)
  }
  if (!is.null(coefficients)) {
    if (!is.matrix(coefficients) || !is.numeric(coefficients) ||
        ncol(coefficients) != nrow(pairs) || nrow(coefficients) < 1L ||
        any(!is.finite(coefficients))) {
      stop(paste0(
        "`coefficients` must be a finite output-by-pair matrix aligned ",
        "with `pairs`."
      ), call. = FALSE)
    }
  }
  structure(list(
    kind = "pair_differences",
    effects = effects,
    pair_left = pairs[, 1L],
    pair_right = pairs[, 2L],
    coefficients = coefficients,
    pair_labels = labels
  ), class = "effect_pair_difference_query")
}

.is_pair_difference_query <- function(query) {
  inherits(query, "effect_pair_difference_query")
}

.query_output_width <- function(query) {
  if (is.null(query)) {
    stop("A query is required to report an output width.", call. = FALSE)
  }
  if (.is_pair_difference_query(query)) {
    if (is.null(query$coefficients)) {
      length(query$pair_left)
    } else {
      nrow(query$coefficients)
    }
  } else {
    ncol(query)
  }
}

.query_output_labels <- function(query) {
  if (.is_pair_difference_query(query)) {
    if (is.null(query$coefficients)) {
      query$pair_labels
    } else {
      labels <- rownames(query$coefficients)
      if (is.null(labels)) {
        labels <- paste0("view", seq_len(nrow(query$coefficients)))
      }
      labels
    }
  } else {
    labels <- colnames(query)
    if (is.null(labels)) labels <- paste0("view", seq_len(ncol(query)))
    labels
  }
}

# The estimand-bearing content of a query, independent of its physical
# representation, used by route-stable view identities.
.query_identity_semantic <- function(query) {
  if (.is_pair_difference_query(query)) {
    list(
      kind = "pair_differences",
      pair_left = query$pair_left,
      pair_right = query$pair_right,
      coefficients = if (is.null(query$coefficients)) {
        NULL
      } else {
        unname(query$coefficients)
      }
    )
  } else {
    unname(query)
  }
}

.query_payload_bytes <- function(query) {
  if (is.null(query)) {
    0
  } else if (.is_pair_difference_query(query)) {
    8 * (2 * length(query$pair_left) +
      if (is.null(query$coefficients)) 0 else length(query$coefficients))
  } else {
    8 * length(query)
  }
}

# Per-edge structured evaluation for effect-by-feature relation blocks:
# returns the pair-by-feature difference product (u^T L) * (u^T R).
.pair_difference_edge_products <- function(query, left, right,
                                           pair_indices = NULL) {
  pair_indices <- if (is.null(pair_indices)) {
    seq_along(query$pair_left)
  } else {
    as.integer(pair_indices)
  }
  dl <- left[query$pair_left[pair_indices], , drop = FALSE] -
    left[query$pair_right[pair_indices], , drop = FALSE]
  dr <- right[query$pair_left[pair_indices], , drop = FALSE] -
    right[query$pair_right[pair_indices], , drop = FALSE]
  dl * dr
}

# Validate a structured query against a self-form task's effect axes.
.validate_pair_difference_for_task <- function(query, left_effects,
                                               right_effects, same_relation) {
  if (!isTRUE(same_relation) || !identical(left_effects, right_effects) ||
      !identical(query$effects, left_effects)) {
    stop(paste0(
      "Pair-difference queries are self-form queries: their effects must ",
      "match one shared relation effect axis."
    ), call. = FALSE)
  }
  invisible(query)
}
