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
  if (!.is_strings(effects, unique = TRUE) || q < 2L) {
    .input_error(sprintf(paste0(
      "Pair-difference views need at least two distinct named effects; the ",
      "relation declares %s (%s). A distance is a comparison between two ",
      "effects."
    ), .msg_count(q, "effect"), .msg_names(effects)))
  }
  if (is.null(pairs)) {
    pairs <- t(utils::combn(seq_len(q), 2L))
  }
  if (is.character(pairs) || is.factor(pairs)) {
    if (!is.matrix(pairs)) {
      .input_error(sprintf(paste0(
        "`pairs` must be a two-column matrix naming the effect pairs to ",
        "report, for example `cbind(\"%s\", \"%s\")`; received %s."
      ), effects[[1L]], effects[[2L]], .msg_value(pairs)))
    }
    pairs <- matrix(as.character(pairs), nrow(pairs), ncol(pairs))
    resolved <- matrix(match(pairs, effects), nrow(pairs), ncol(pairs))
    if (anyNA(resolved)) {
      .input_error(sprintf(paste0(
        "`pairs` names %s, which %s a declared experimental effect. The ",
        "relation declares %s."
      ), .msg_names(unique(pairs[is.na(resolved)])),
        if (length(unique(pairs[is.na(resolved)])) == 1L) "is not" else
          "are not",
        .msg_names(effects)))
    }
    pairs <- resolved
  }
  if (!is.matrix(pairs) || ncol(pairs) != 2L || nrow(pairs) < 1L ||
      !is.numeric(pairs) || anyNA(pairs) || any(pairs %% 1 != 0) ||
      any(pairs < 1L) || any(pairs > q)) {
    .input_error(sprintf(paste0(
      "`pairs` must be a two-column matrix of effect names or of indices in ",
      "1:%d; received %s."
    ), q, .msg_value(pairs)))
  }
  storage.mode(pairs) <- "integer"
  if (any(pairs[, 1L] == pairs[, 2L])) {
    .input_error(sprintf(paste0(
      "`pairs` asks for the distance from %s to itself, which is zero by ",
      "construction; every pair must name two distinct effects."
    ), .msg_names(unique(effects[pairs[pairs[, 1L] == pairs[, 2L], 1L]]))))
  }
  swapped <- pairs[, 1L] > pairs[, 2L]
  if (any(swapped)) {
    pairs[swapped, ] <- pairs[swapped, c(2L, 1L), drop = FALSE]
  }
  if (anyDuplicated(paste(pairs[, 1L], pairs[, 2L], sep = "\r"))) {
    .input_error("`pairs` contains duplicate effect pairs.")
  }
  if (is.null(labels)) {
    labels <- paste(effects[pairs[, 1L]], "-", effects[pairs[, 2L]])
  }
  if (!is.character(labels) || length(labels) != nrow(pairs) ||
      anyNA(labels)) {
    .input_error("Pair labels must name every pair.")
  }
  if (!is.null(coefficients)) {
    if (!.is_finite_matrix(coefficients) ||
        ncol(coefficients) != nrow(pairs) || nrow(coefficients) < 1L) {
      .input_error(paste0(
        "`coefficients` must be a finite output-by-pair matrix aligned ",
        "with `pairs`."
      ))
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
    .input_error("A query is required to report an output width.")
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
    .input_error(paste0(
      "Pair-difference queries are self-form queries: their effects must ",
      "match one shared relation effect axis."
    ))
  }
  invisible(query)
}
