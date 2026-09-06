# Model-coordinate geometry: the reader (layer 5) ----------------------------
#
# `design/model-coordinate-geometry-contract.md` sections 4, 5 and 6. Sits
# beside `latent.R` and calls down or sideways only: `geometry-plan.R` and
# `geometry-entry.R` to re-plan and execute the lowered relation (the way
# `population-heterogeneity.R` re-plans through `plan_geometry()`),
# `relation.R` and `model-basis.R` (layer 2) for the lowering itself,
# `latent.R` (layer 5, sideways) for the rank-limited PSD projection and its
# two-part moved-mass accounting, and `format-results.R` for printing.
#
# Two executions serve every reading: one materializes the lowered form, one
# reads two trace queries from the original plan. The view takes the ORIGINAL
# geometry plan and the model basis, lowers the plan's relation through the basis
# (`B~ = Q'B`, composing the basis map with each partition's extractor),
# materializes the lowered form `S_x = Q' G_x Q` per measurement, and reads
# the two trace queries `tr(H G_x H)` and `tr((H - P) G_x)` from the original
# plan as one packed query bank. From those it reports:
#
#   * the TRACE DECOMPOSITION (section 5), on the signed estimation layer:
#     centered total = model-addressable + model-orthogonal, exact at every
#     measurement, with addressable = tr(S_x). The orthogonal term is read
#     from the original plan and is not available from the lowered form
#     alone; the object says which plan it came from.
#
#   * one STRUCTURE FIT (section 4), on the latent layer: isotropic and
#     diagonal are nonnegative quadratic programs in the model-only Gram
#     `M = (T'T)^2` with `g_j = R[, j]' S_x R[, j]` read from the compressed
#     form; block is exact block coordinate descent with each block's update
#     the rank-limited PSD projection; shared IS `latent_geometry()`'s
#     rank-budgeted truncation, called through `.latent_rank_psd_form()` and
#     not implemented twice. The returned object is `C` in model coordinates
#     (never a factor `W`), the residual `||S_x - A||_F^2` in model
#     coordinates, and the fitted energy `tr(A)`.
#
# What a fit tests is on the record beside it (section 4's table), because
# the four structures answer four different questions and only the isotropic
# one tests the model geometries as supplied. A shared fit on a basis that
# spans every centered direction tests nothing about the models (section 6,
# Theorem 3) and is refused through `.model_basis_refuse_saturated()`; above
# a declared warning fraction it is admitted and labelled `near_saturated`.
# No fraction, cumulative curve or "variance explained" is formed on the
# signed `S_x`: the fitted energy is a latent quantity, printed with the
# latent reading line, and never as a share of a signed denominator.

.model_geometry_structures <- c("isotropic", "diagonal", "block", "shared")

.model_geometry_hypotheses <- c(
  isotropic = "the truncated model geometries with nonnegative weights",
  diagonal = "the model axes with unequal, possibly sparse weights",
  block = "an independent learned metric per model",
  shared = paste0("a learned rank-limited readout of the joint model span, ",
    "not the supplied geometries")
)

# With one model the block structure has one block, which is the whole span:
# the block fit is then the shared fit under another name, and it is
# described and gated as the shared fit (Theorem 3) rather than admitted as
# "an independent metric per model".
.model_geometry_hypothesis <- function(structure, n_models) {
  if (identical(structure, "block") && n_models == 1L) {
    return(paste0("a learned rank-limited readout of the single model's span ",
      "(with one model the block fit is the shared fit), not the supplied ",
      "geometry"))
  }
  unname(.model_geometry_hypotheses[[structure]])
}

.model_geometry_learns_rotation <- function(structure, n_models) {
  identical(structure, "shared") ||
    (identical(structure, "block") && n_models == 1L)
}

# Nonnegative least squares in quadratic form: minimize 0.5 c'Mc - g'c over
# c >= 0, for a PSD M. Lawson and Hanson's active-set method, written out
# because the problems here are tiny (one variable per model coordinate) and
# a package dependency for a forty-line routine is the wrong trade. At the
# solution the KKT conditions hold: c >= 0, the gradient M c - g >= 0, and
# their elementwise product is zero; the tests check exactly that.
.nnls_quadratic <- function(M, g, tolerance = 1e-12, max_iterations = NULL) {
  m <- length(g)
  if (is.null(max_iterations)) max_iterations <- 50L * max(1L, m)
  scale <- max(1, max(abs(g)), max(abs(M)))
  c <- numeric(m)
  passive <- logical(m)
  solve_passive <- function(index) {
    z <- numeric(m)
    if (!length(index)) return(z)
    block <- M[index, index, drop = FALSE]
    z[index] <- tryCatch(solve(block, g[index]), error = function(e) {
      sv <- svd(block)
      keep <- sv$d > 1e-12 * max(sv$d)
      sv$v[, keep, drop = FALSE] %*% ((t(sv$u[, keep, drop = FALSE]) %*%
        g[index]) / sv$d[keep])
    })
    z
  }
  iterations <- 0L
  repeat {
    gradient <- as.numeric(M %*% c) - g
    candidates <- which(!passive & gradient < -tolerance * scale)
    if (!length(candidates) || iterations >= max_iterations) break
    j <- candidates[which.min(gradient[candidates])]
    passive[j] <- TRUE
    repeat {
      iterations <- iterations + 1L
      index <- which(passive)
      z <- solve_passive(index)
      # Positivity is tested exactly, as Lawson and Hanson do; `tolerance`
      # is the optimality threshold only, so a coarse tolerance cannot make
      # the entering variable look blocked and stall the loop.
      if (all(z[index] > 0)) {
        c <- z
        break
      }
      blocking <- index[z[index] <= 0]
      alpha <- min(c[blocking] / (c[blocking] - z[blocking]))
      c <- c + alpha * (z - c)
      passive[c <= 0] <- FALSE
      c[!passive] <- 0
      if (!any(passive) || iterations >= max_iterations) break
    }
  }
  c[c < 0] <- 0
  list(solution = c, iterations = iterations,
    converged = iterations < max_iterations)
}

# The model-only pieces every fit shares, computed once per basis: the
# stacked model Gram in model coordinates, its elementwise square, each
# model's column set, and for the block structure an orthonormal basis of
# each model's span in Q-coordinates with the square `R_i` that carries block
# coordinates back to that model's own coordinates.
.model_geometry_model_side <- function(basis) {
  if (ncol(basis$R) > nrow(basis$R)) {
    .capability_refusal(paste0(
      "This reader reports coefficients in model coordinates, which are ",
      "not uniquely identified when model spans overlap. The model basis ",
      "remains valid for invariant form prediction."),
      capability = "identified_model_coordinates", namespace = "model_coordinate",
      reasons = c("model_span_overlap", sprintf("span_overlap_rank = %d.",
        ncol(basis$R) - nrow(basis$R))),
      remedies = c("Reduce the ranks for a coefficient-based model_geometry() reading.",
        "Drop or merge redundant models for that coefficient question.",
        "Use fit_geometry() for invariant form prediction with overlapping models."))
  }
  R <- unname(basis$R)
  columns <- basis$columns
  gram <- crossprod(R)
  blocks <- lapply(columns, function(cols) {
    sv <- svd(R[, cols, drop = FALSE])
    keep <- sv$d > basis$tolerance * sv$d[[1L]]
    Q_i <- sv$u[, keep, drop = FALSE]
    list(columns = cols, Q = Q_i, R = t(Q_i) %*% R[, cols, drop = FALSE])
  })
  list(R = R, columns = columns, gram = gram, square = gram^2, blocks = blocks,
    R_inverse = solve(R))
}

.model_geometry_block_diagonal <- function(parts, columns, m) {
  C <- matrix(0, m, m)
  for (i in seq_along(parts)) {
    cols <- columns[[i]]
    C[cols, cols] <- parts[[i]]
  }
  C
}

# Exact block coordinate descent (section 4.2, oracle section O8). Each
# block update is the rank-limited PSD projection of that block's residual,
# which is the exact minimizer over the block with the others held fixed, so
# the objective never increases; the sweep stops when no block moves.
.model_geometry_block_descent <- function(S, model, ranks, tolerance,
                                          max_sweeps) {
  blocks <- model$blocks
  q <- nrow(S)
  A <- lapply(blocks, function(block) matrix(0, ncol(block$Q), ncol(block$Q)))
  lift <- function(i) blocks[[i]]$Q %*% A[[i]] %*% t(blocks[[i]]$Q)
  scale <- max(1, max(abs(S)))
  converged <- FALSE
  sweeps <- 0L
  while (sweeps < max_sweeps) {
    sweeps <- sweeps + 1L
    moved <- 0
    for (i in seq_along(blocks)) {
      others <- matrix(0, q, q)
      for (j in seq_along(blocks)[-i]) others <- others + lift(j)
      residual <- t(blocks[[i]]$Q) %*% (S - others) %*% blocks[[i]]$Q
      updated <- .latent_rank_psd_form(residual, ranks[[i]])$form
      moved <- max(moved, max(abs(updated - A[[i]])))
      A[[i]] <- updated
    }
    if (moved <= tolerance * scale || length(blocks) == 1L) {
      # One block has nothing to alternate with: its first exact update is
      # the minimizer.
      converged <- TRUE
      break
    }
  }
  fitted <- matrix(0, q, q)
  for (i in seq_along(blocks)) fitted <- fitted + lift(i)
  parts <- lapply(seq_along(blocks), function(i) {
    R_i_inverse <- solve(blocks[[i]]$R)
    R_i_inverse %*% A[[i]] %*% t(R_i_inverse)
  })
  list(A = fitted, C = .model_geometry_block_diagonal(parts, model$columns,
    ncol(model$R)), iterations = sweeps, converged = converged)
}

# One measurement, one structure. Returns the fitted form `A` in
# Q-coordinates, the model-side metric `C`, and the fit accounting.
.model_geometry_fit <- function(S, structure, model, ranks, tolerance,
                                max_sweeps) {
  R <- model$R
  m <- ncol(R)
  clipped <- NA_real_
  truncated <- NA_real_
  iterations <- NA_integer_
  converged <- TRUE
  weights <- NULL
  if (structure %in% c("isotropic", "diagonal")) {
    g <- vapply(seq_len(m), function(j) {
      drop(t(R[, j]) %*% S %*% R[, j])
    }, numeric(1))
    if (structure == "diagonal") {
      fit <- .nnls_quadratic(model$square, g, tolerance)
      weights <- fit$solution
      C <- diag(weights, m)
    } else {
      k <- length(model$columns)
      g_iso <- vapply(model$columns, function(cols) sum(g[cols]), numeric(1))
      M_iso <- matrix(0, k, k)
      for (i in seq_len(k)) for (l in seq_len(k)) {
        M_iso[i, l] <- sum(model$square[model$columns[[i]], model$columns[[l]]])
      }
      fit <- .nnls_quadratic(M_iso, g_iso, tolerance)
      weights <- fit$solution
      C <- diag(rep(weights, times = lengths(model$columns)), m)
    }
    iterations <- fit$iterations
    converged <- fit$converged
    A <- R %*% C %*% t(R)
  } else if (structure == "block") {
    fit <- .model_geometry_block_descent(S, model, ranks, tolerance, max_sweeps)
    A <- fit$A
    C <- fit$C
    iterations <- fit$iterations
    converged <- fit$converged
  } else {
    projected <- .latent_rank_psd_form(S, ranks)
    A <- projected$form
    C <- model$R_inverse %*% A %*% t(model$R_inverse)
    clipped <- projected$clipped_negative_mass
    truncated <- projected$truncated_positive_mass
    iterations <- 1L
  }
  C <- (C + t(C)) / 2
  list(
    A = A, C = C, weights = weights,
    residual = sum((S - A)^2),
    fitted_energy = sum(diag(A)),
    iterations = as.integer(iterations),
    converged = isTRUE(converged),
    clipped_negative_mass = clipped,
    truncated_positive_mass = truncated
  )
}

# The rank budgets a structure takes: none for isotropic and diagonal, one
# per model for block, one for shared. Requested budgets are clamped to the
# dimension they act on and the clamp is recorded, as in the basis.
.model_geometry_ranks <- function(rank, basis, structure) {
  sizes <- lengths(basis$columns)
  q <- basis$dimension
  if (structure %in% c("isotropic", "diagonal")) {
    if (!is.null(rank)) {
      .input_error(sprintf(paste0(
        "`rank` applies to the block and shared structures; the %s ",
        "structure learns no rotation and has no rank budget."
      ), structure), arg = "rank", received = .msg_value(rank),
        expected = "NULL")
    }
    return(list(requested = NULL, effective = NULL, clamped = FALSE,
      total = as.integer(q)))
  }
  if (structure == "shared") {
    requested <- if (is.null(rank)) q else {
      .check_count(rank, "rank", what = "NULL or one positive whole number")
    }
    effective <- min(requested, q)
    return(list(requested = as.integer(requested),
      effective = as.integer(effective), clamped = requested > effective,
      total = as.integer(effective)))
  }
  requested <- if (is.null(rank)) sizes else {
    if (!is.numeric(rank) || !length(rank) %in% c(1L, length(sizes)) ||
        anyNA(rank) || any(!is.finite(rank)) || any(rank < 1) ||
        any(rank %% 1 != 0)) {
      .input_error(sprintf(paste0(
        "`rank` must be NULL, one positive whole number, or one per model ",
        "(%s); received %s."
      ), .msg_names(names(sizes)), .msg_value(rank)),
        arg = "rank", received = .msg_value(rank),
        expected = "NULL, one positive whole number, or one per model")
    }
    if (!is.null(names(rank))) {
      if (!setequal(names(rank), names(sizes)) || anyDuplicated(names(rank))) {
        .input_error(sprintf("`rank` names (%s) must match the model names (%s).",
          .msg_names(names(rank)), .msg_names(names(sizes))),
          arg = "rank", received = .msg_names(names(rank)),
          expected = .msg_names(names(sizes)))
      }
      rank <- rank[names(sizes)]
    }
    if (length(rank) == 1L) rank <- rep(rank, length(sizes))
    stats::setNames(as.integer(unname(rank)), names(sizes))
  }
  effective <- pmin(requested, sizes)
  list(requested = requested, effective = stats::setNames(as.integer(effective),
    names(sizes)), clamped = any(requested > effective),
    total = as.integer(sum(effective)))
}

# Degrees of freedom of the structure (section 6(a)): the free parameters of
# the fitted `C`, beside those of the unconstrained rank-s PSD form on n - 1
# centered directions, so a reader can see how much the model family
# constrains. A rank-s PSD matrix on d dimensions has d s - s (s - 1) / 2
# free parameters.
.model_geometry_df <- function(structure, basis, ranks) {
  psd_df <- function(d, s) d * s - s * (s - 1) / 2
  sizes <- lengths(basis$columns)
  q <- basis$dimension
  n <- basis$n_conditions
  fit <- switch(structure,
    isotropic = length(sizes),
    diagonal = sum(sizes),
    block = sum(psd_df(sizes, ranks$effective)),
    shared = psd_df(q, ranks$effective)
  )
  s_free <- min(ranks$total, n - 1L)
  list(fit = as.numeric(fit), free = as.numeric(psd_df(n - 1L, s_free)))
}

# The edge-disjoint cross-fit of the learned readout (contract section 7).
# For each evaluation edge (a, b) of the plan's pairing, the rank-s
# eigenbasis V of the compressed form is learned from the edges among the
# OTHER partitions and the energy tr(V' S^(a,b) V) is read on the edge
# itself; the edges are then reduced with the pairing's own weights. Both
# the evaluation edges and the training edges are the DECLARED pairing's
# edges: a pairing states exactly which partition products a plan may form,
# so the readout is never trained on a product the pairing did not declare,
# and a pairing in which some edge has no disjoint edge (a star) refuses. V
# is a function of data independent of the products it weights, so the
# energy is signed, on the estimation layer, and zero in expectation under
# pure noise (oracle section O7) when the pairing's partitions are
# independent, where the plug-in latent energy is not. The unit is the
# partition edge, the same unit `metric_training_policy()` uses for a
# learned neural metric, and the four-partition floor is the one
# `heterogeneity()` enforces: an evaluation edge needs two partitions and its
# training edges need two more. This version defines the cross-fit for the
# shared structure only; the isotropic, diagonal and block plug-in energies
# are fitted quantities too and are not unbiased under the null.
.model_geometry_cross_fit <- function(x, lowered, lowered_plan, rank, component,
                                      policy) {
  refuse <- function(detail, remedies) {
    .capability_refusal(paste0(
      "The cross-fitted model energy learns the readout on partition edges ",
      "disjoint from the edge it is evaluated on, and the declared pairing ",
      "does not provide them: ", detail, "."
    ),
      capability = "cross_fitted_model_energy",
      namespace = "model_coordinate",
      reasons = c("insufficient_disjoint_edges", detail),
      remedies = c(remedies,
        paste0("Read the plug-in rank-limited fit (`training = NULL`), which is ",
          "a latent quantity and not an unbiased estimate.")))
  }
  independence <- attr(x$pairing, "independence", exact = TRUE)
  if (identical(independence, "undeclared")) independence <- NULL
  generalizes_over <- attr(x$pairing, "generalizes_over", exact = TRUE)
  # The evaluation edges are the pairing's own cross-partition edges, each
  # orientation of one unordered edge folded into it with its weight added.
  declared <- as.data.frame(x$pairing)
  declared <- declared[declared$left != declared$right, , drop = FALSE]
  if (!nrow(declared)) {
    refuse("no cross-partition evaluation edges: the pairing declares only self products",
      "Pair distinct partitions with `cross_partitions()` or `pairing()`.")
  }
  key <- vapply(seq_len(nrow(declared)), function(i) {
    paste(sort(c(declared$left[[i]], declared$right[[i]])), collapse = "\r")
  }, character(1))
  # Edges keep the pairing's declared order (first appearance of each
  # unordered edge), so `by_edge` and `edges` read in the order the user
  # wrote the pairing.
  edges <- do.call(rbind, lapply(unique(key), function(k) {
    group <- declared[key == k, , drop = FALSE]
    data.frame(left = group$left[[1L]], right = group$right[[1L]],
      weight = sum(group$weight), stringsAsFactors = FALSE)
  }))
  rownames(edges) <- NULL
  partitions <- unique(c(edges$left, edges$right))
  if (length(partitions) < 4L) {
    refuse(sprintf("partitions: %d, floor: 4", length(partitions)),
      "Pair at least four independent partitions.")
  }
  weights <- edges$weight / sum(edges$weight)
  form_for <- function(over) {
    plan <- plan_geometry(lowered, x$frame, over, compute = x$compute,
      metric = x$metric_schedule$metric)
    if (!identical(plan$metric_schedule$signature,
        lowered_plan$metric_schedule$signature) ||
        !identical(plan$codec, lowered_plan$codec) ||
        !identical(plan$packed_width, lowered_plan$packed_width) ||
        !identical(plan$measurements, lowered_plan$measurements)) {
      .invariant_error(paste0(
        "Restricting the lowered plan to a partition edge changed it beyond ",
        "its pairing, so the training and evaluation forms would not describe ",
        "the same measurements."
      ))
    }
    geometry_component(materialize_geometry(plan), component)
  }
  q <- length(lowered$effects)
  by_edge <- matrix(NA_real_, lowered_plan$measurements, nrow(edges))
  trained <- character(nrow(edges))
  training_edges <- integer(nrow(edges))
  for (e in seq_len(nrow(edges))) {
    evaluation <- c(edges$left[[e]], edges$right[[e]])
    disjoint <- !(edges$left %in% evaluation) & !(edges$right %in% evaluation)
    if (!any(disjoint)) {
      refuse(sprintf("edge %s|%s has no disjoint training edge in the pairing",
        evaluation[[1L]], evaluation[[2L]]),
        "Declare a pairing in which every edge has an edge disjoint from it.")
    }
    train <- edges[disjoint, , drop = FALSE]
    trained[[e]] <- paste(sort(unique(c(train$left, train$right))),
      collapse = ",")
    training_edges[[e]] <- nrow(train)
    train_packed <- form_for(pairing(train$left, train$right,
      weight = train$weight, directed = FALSE, independence = independence,
      generalizes_over = generalizes_over))
    eval_packed <- form_for(pairing(evaluation[[1L]], evaluation[[2L]],
      directed = FALSE, independence = independence,
      generalizes_over = generalizes_over))
    for (i in seq_len(nrow(by_edge))) {
      S_train <- .unsvec_symmetric(train_packed[i, ], q)
      S_eval <- .unsvec_symmetric(eval_packed[i, ], q)
      V <- eigen(S_train, symmetric = TRUE)$vectors[, seq_len(rank), drop = FALSE]
      by_edge[i, e] <- sum(diag(t(V) %*% S_eval %*% V))
    }
  }
  colnames(by_edge) <- paste(edges$left, edges$right, sep = "|")
  list(
    energy = as.numeric(by_edge %*% weights),
    by_edge = by_edge,
    edges = data.frame(left = edges$left, right = edges$right,
      weight = weights, trained_on = trained, training_edges = training_edges,
      stringsAsFactors = FALSE),
    rank = as.integer(rank),
    policy = policy$kind,
    policy_signature = policy$signature,
    independence = if (is.null(independence)) "undeclared" else independence,
    layer = "estimation"
  )
}

.model_geometry_source <- function(x, basis) {
  x <- .check_class(x, "effect_geometry_plan", "x", from = "plan_geometry()")
  # A self-form plan; the rectangular refusal is the shared one.
  .self_geometry_source(x, "A model geometry")
  schedule <- x$metric_schedule
  if (!schedule$kind %in% c("implicit_identity_before_frame",
      "fixed_metric_before_frame")) {
    .capability_refusal(sprintf(paste0(
      "A model geometry lowers the plan's relation and re-plans it under the ",
      "same metric, and this plan's metric schedule (`%s`) cannot be ",
      "re-planned: a whitened or learned metric freezes coordinates or ",
      "statistics on the original effect axis, and the compiler does not ",
      "materialize a complete lowered form under a learned local metric."
    ), schedule$kind),
      capability = "model_coordinate_lowering",
      namespace = "model_coordinate",
      reasons = c("metric_schedule_not_lowerable",
        paste0("metric_schedule_kind:", schedule$kind)),
      remedies = paste0("Plan the geometry on the implicit identity metric or ",
        "one declared fixed `neural_metric()` with `composition = \"native\"`."))
  }
  relation <- x$task$left_relation
  if (!identical(relation$effect_space$coordinates, basis$conditions)) {
    .input_error(sprintf(paste0(
      "The model basis is declared over conditions (%s) that are not the ",
      "plan's effects in the plan's order (%s). Build the basis with ",
      "`conditions = <the relation's effect_space>`."
    ), .msg_names(basis$conditions), .msg_names(relation$effect_space$coordinates)),
      arg = "basis", received = .msg_names(basis$conditions),
      expected = .msg_names(relation$effect_space$coordinates))
  }
  if (.is_model_coordinate_space(relation$effect_space)) {
    .input_error(paste0(
      "`x` is already a lowered plan: pass the plan on the original relation, ",
      "so the model-orthogonal term can be read from it."
    ), arg = "x", received = "a plan on a model-coordinate effect space",
      expected = "a plan on the original condition effects")
  }
  x
}

#' Read a geometry in model coordinates: the trace split and a structure fit
#'
#' A model geometry reads one crossvalidated geometry through a
#' [model_basis()]. It lowers the plan's relation into model coordinates
#' (`B~ = Q'B`, composing the basis with each partition's extractor),
#' materializes the complete compressed form `S_x = Q' G_x Q` at every
#' measurement, and reads two fixed trace queries from the original plan.
#' From those it reports two things that live on different layers and are
#' kept apart.
#'
#' On the **signed estimation layer**, the trace decomposition
#' \eqn{\operatorname{tr}(H G_x H) = \operatorname{tr}(P G_x) +
#' \operatorname{tr}((H - P) G_x)}{tr(H G H) = tr(P G) + tr((H - P) G)}:
#' the centered total energy of the geometry equals its model-addressable
#' energy, which is `tr(S_x)`, plus its model-orthogonal energy. Both terms
#' are fixed bilinear queries with vanishing marginals, so both are unbiased
#' and signed, and the identity is exact at every measurement. It mirrors
#' `coherent + configuration = total` along the representational axis
#' instead of the spatial one.
#'
#' On the **latent layer**, one structure fit of `S_x` to the model family:
#' the nonnegative model-side metric `C` in the parameterization
#' `G_x ~ T C T'`, with `T = Q R` the stacked model factor. The four
#' structures are readings of one form and answer four different questions:
#'
#' | structure | `C` | tests |
#' |---|---|---|
#' | `isotropic` | one nonnegative weight per model | the truncated model geometries with nonnegative weights |
#' | `diagonal` | one nonnegative weight per model coordinate | the model axes with unequal, possibly sparse weights |
#' | `block` | one PSD rank-limited block per model | an independent learned metric per model (with one model, the shared fit) |
#' | `shared` | one PSD rank-limited matrix over the joint span | a learned readout of the joint model span, not the supplied geometries |
#'
#' Only the isotropic structure tests the model geometries as supplied. The
#' shared fit is [latent_geometry()]'s rank-budgeted truncation of `S_x`
#' (`[S_x]_{+,s}`, the closed-form Frobenius minimizer over PSD forms of
#' rank at most `s`), not a second implementation of it. Every fit returns
#' `C` and never a factor `W`, because `W` is identified only up to rotation.
#' The fitted energy `tr(A)` and the residual `||S_x - A||_F^2` are latent
#' quantities and are never reported as fractions of a signed denominator.
#'
#' @param x An `effect_geometry_plan` from [plan_geometry()] on the
#'   **original** relation (condition effects, not model coordinates), a
#'   complete self form under the implicit identity metric or one fixed
#'   `neural_metric()`. The view re-plans the lowered relation with the same
#'   frame, pairing, metric and compute policy and executes both; there is no
#'   separate compute override, because the plan already carries one.
#' @param basis A [model_basis()] declared over the plan's effects in the
#'   plan's order (build it with `conditions = relation$effect_space`).
#' @param structure Which model-side metric to fit; see the table above.
#' @param rank The rank budget: `NULL` for none (block: each model's own
#'   dimension; shared: the basis dimension), one positive whole number
#'   (shared: the budget; block: shared by every model), or one per model for
#'   block, named by model. Budgets are clamped to the dimension they act on
#'   and the clamp is recorded. Not accepted for isotropic or diagonal.
#' @param component Geometry component to read and fit.
#' @param warning_fraction Span fraction at or above which a shared fit is
#'   admitted but labelled `near_saturated`. At a span fraction of one the
#'   shared fit is refused (see Refusals).
#' @param tolerance Positive relative tolerance for the nonnegative and block
#'   solvers' convergence.
#' @param max_sweeps Positive cap on block coordinate descent sweeps per
#'   measurement.
#' @param training `NULL`, or [metric_training_policy()] with kind
#'   `"exclude_evaluation"`, which additionally cross-fits the learned
#'   readout of the shared structure: for each evaluation edge of the plan's
#'   pairing the rank-`s` eigenbasis of the compressed form is learned from
#'   the pairing's edges disjoint from it (never from a product the pairing
#'   did not declare), the energy is read on the edge itself, and the edges
#'   are reduced with the pairing's weights. The result, `$cross_fit`, is
#'   signed and on the estimation layer: zero in expectation under pure noise
#'   when the pairing's partitions are independent, where the plug-in fitted
#'   energy is not. It needs a pairing over at least four partitions in which
#'   every edge has a disjoint edge, and refuses otherwise. Defined for the
#'   shared structure only in this version. A tie in the training spectrum
#'   straddling the rank cut would make the readout platform dependent; for
#'   continuous data that is a measure-zero event.
#' @return An `effect_model_geometry`.
#' @section Structure:
#' One row per spatial measurement, with the estimation-layer split and the
#' latent-layer fit kept in separately named fields.
#'
#' - `$addressable`, `$orthogonal`, `$centered_total`: the signed trace
#'   decomposition per measurement; `$identity_error` is the largest
#'   departure of `addressable + orthogonal` from `centered_total`.
#' - `$fit`: a data frame per measurement with `residual` (Frobenius, in
#'   model coordinates), `fitted_energy` (`tr(A)`), `iterations`,
#'   `converged` (block descent needs `max_sweeps` of at least two to
#'   observe a quiet sweep unless the basis has one model), and for the
#'   shared structure the two-part moved mass `clipped_negative_mass` and
#'   `truncated_positive_mass`.
#' - `$coefficients`: the fitted `C`, a model-coordinate by model-coordinate
#'   by measurement array, coordinates named `<model>.<k>` as in the basis's
#'   `$R`. For `isotropic` and `diagonal`, `$weights` additionally holds the
#'   nonnegative weights per model or per coordinate.
#' - `$cross_fit` (with `training`): `energy` per measurement, the
#'   per-edge energies `by_edge`, the `edges` with their weights, the
#'   partitions each readout was trained on and how many pairing edges
#'   trained it, the `rank`, the policy, and the pairing's declared
#'   `independence`. It is listed under `$layer$estimation`; the plug-in
#'   `$fit$fitted_energy` stays on the latent layer, and the two are never
#'   the same number. The execution `$receipt` carries none of this; its
#'   schema is sealed.
#' - `$structure`, `$hypothesis`, `$tests_supplied_geometries`: which fit
#'   was made and what it tests.
#' - `$rank`: the requested and effective budgets and whether any was
#'   clamped.
#' - `$df`: `fit`, the free parameters of the structure, beside `free`, those
#'   of the unconstrained rank-limited PSD form on `n - 1` centered
#'   directions; `$model_sizes` are the per-model dimensions they were
#'   computed from.
#' - `$span_fraction`, `$near_saturated`: the basis's span fraction and the
#'   warning label.
#' - `$basis_signature`, `$source_scientific_plan_id`,
#'   `$lowered_scientific_plan_id`: the identities the reading was built
#'   from; `$receipt` is derived from the lowered execution.
#' - `$reading`: the latent reading line, which applies to `$fit`,
#'   `$coefficients` and `$weights` and not to the trace split.
#'
#' Any element not listed here is internal and may change.
#' @section Refusals:
#' Each is an `effect_capability_refusal` (see [catch_refusal()]) in
#' namespace `"model_coordinate"`. A plan on a whitened or learned metric
#' schedule refuses capability `"model_coordinate_lowering"` with reason
#' `metric_schedule_not_lowerable`; a shared fit on a basis whose span
#' fraction is one refuses capability `"model_geometry_test"` with reason
#' `model_span_saturated`, because a full-span basis reproduces the best
#' rank-limited PSD approximation of the centered geometry whatever models
#' produced it; a block fit over a single model is the shared fit under
#' another name and is gated the same way; a cross-fit on fewer than four
#' partitions refuses capability `"cross_fitted_model_energy"` with reason
#' `insufficient_disjoint_edges`.
#' @seealso [model_basis()] for the basis; [latent_geometry()] for the
#'   rank-budgeted projection the shared fit delegates to; [geometry_spectrum()]
#'   for the signed spectrum of a lowered geometry; [rsa()] for the fixed
#'   linear reading of model RDMs on the original effects.
#' @family geometry plans and views
#' @examples
#' example <- example_fmri_effects()
#' relation <- example$fit$relation
#' plan <- plan_geometry(relation, example$frame,
#'   cross_partitions(relation, independence = "independent",
#'     generalizes_over = "run"))
#'
#' # The README category model as a basis over the plan's four conditions:
#' # its Gram has rank one, so the basis spans one of three centered
#' # directions.
#' basis <- model_basis(list(category = example$model_rdm),
#'   distance = "squared_euclidean", conditions = relation$effect_space)
#' basis$span_fraction
#'
#' # The trace split at every searchlight, and the isotropic fit: one
#' # nonnegative weight on the category geometry.
#' reading <- model_geometry(plan, basis)
#' reading
#' head(as.data.frame(reading))
#' reading$identity_error
#'
#' # The shared fit learns a readout of the model span; here the span is one
#' # dimensional so it coincides with the isotropic fit's energy.
#' shared <- model_geometry(plan, basis, structure = "shared")
#' all.equal(shared$fit$fitted_energy, reading$fit$fitted_energy)
#' @export
model_geometry <- function(x, basis,
                           structure = c("isotropic", "diagonal", "block",
                             "shared"),
                           rank = NULL,
                           component = c("total", "coherent", "configuration"),
                           warning_fraction = 0.9, tolerance = 1e-10,
                           max_sweeps = 1000L, training = NULL) {
  structure <- match.arg(structure)
  component <- match.arg(component)
  policy <- NULL
  if (!is.null(training)) {
    policy <- .validate_metric_training_policy(training)
    if (!identical(structure, "shared")) {
      .input_error(sprintf(paste0(
        "`training` cross-fits the learned readout of the shared structure ",
        "only in this version. The %s structure's plug-in energy is a fitted ",
        "quantity too and is not unbiased under the null; read it on the ",
        "latent layer, or pass `structure = \"shared\"` for the cross-fit."
      ), structure), arg = "training", received = .msg_value(training),
        expected = "NULL, or `structure = \"shared\"`")
    }
    if (!identical(policy$kind, "exclude_evaluation")) {
      .input_error(paste0(
        "The cross-fitted model energy admits only ",
        "`metric_training_policy(\"exclude_evaluation\")`: the readout is ",
        "learned on partition edges disjoint from the evaluation edge, and no ",
        "residual-reuse policy has a meaning on the model side."
      ), arg = "training", received = policy$kind,
        expected = "\"exclude_evaluation\"")
    }
  }
  .check_number(tolerance, "tolerance", positive = TRUE)
  if (tolerance >= 1) {
    .input_error("`tolerance` is relative and must be below one.",
      arg = "tolerance", received = .msg_value(tolerance),
      expected = "a positive number below one")
  }
  max_sweeps <- .check_count(max_sweeps, "max_sweeps")
  .check_number(warning_fraction, "warning_fraction", positive = TRUE)
  if (warning_fraction > 1) {
    .input_error("`warning_fraction` is a span fraction and must be at most one.",
      arg = "warning_fraction", received = .msg_value(warning_fraction),
      expected = "a number in (0, 1]")
  }
  basis <- .validate_model_basis(basis)
  model <- .model_geometry_model_side(basis)
  x <- .model_geometry_source(x, basis)
  ranks <- .model_geometry_ranks(rank, basis, structure)
  if (.model_geometry_learns_rotation(structure, length(basis$columns))) {
    .model_basis_refuse_saturated(basis, structure)
  }
  near_saturated <- basis$span_fraction >= warning_fraction

  # The lowered plan: the same frame, pairing, metric and compute policy on
  # the relation lowered through the basis.
  relation <- x$task$left_relation
  lowered <- .relation_lowered(relation, basis)
  lowered_plan <- plan_geometry(lowered, x$frame, x$pairing,
    compute = x$compute, metric = x$metric_schedule$metric)
  if (!identical(lowered_plan$measurements, x$measurements) ||
      !identical(lowered_plan$metric_schedule$kind, x$metric_schedule$kind)) {
    .invariant_error(paste0(
      "Lowering the relation changed the plan beyond its effect axis, so the ",
      "compressed form and the trace queries would not describe the same ",
      "measurements."
    ))
  }
  compressed <- materialize_geometry(lowered_plan)
  packed <- geometry_component(compressed, component)
  index <- compressed$index

  # The two trace queries on the original plan, as one packed bank: the
  # centering operator and the model-orthogonal projector. Their Euclidean
  # inner product with the packed form is the trace against the operator,
  # because both sides carry the codec's sqrt(2) off-diagonals.
  n <- basis$n_conditions
  H <- diag(n) - 1 / n
  P <- basis$Q %*% t(basis$Q)
  bank <- cbind(centered = .svec_symmetric(H), orthogonal = .svec_symmetric(H - P))
  traces <- evaluate_geometry(x, query = bank, component = component)
  centered_total <- as.numeric(traces$values[, 1L])
  orthogonal <- as.numeric(traces$values[, 2L])
  if (!identical(NROW(traces$index), length(index))) {
    .invariant_error("The trace queries and the compressed form disagree on the measurements.")
  }

  q <- basis$dimension
  fit_ranks <- ranks$effective
  measurements <- nrow(packed)
  m <- ncol(model$R)
  coefficients <- array(NA_real_, c(m, m, measurements),
    dimnames = list(colnames(basis$R), colnames(basis$R), NULL))
  addressable <- numeric(measurements)
  fit <- data.frame(
    residual = numeric(measurements),
    fitted_energy = numeric(measurements),
    iterations = integer(measurements),
    converged = logical(measurements),
    clipped_negative_mass = rep(NA_real_, measurements),
    truncated_positive_mass = rep(NA_real_, measurements)
  )
  weights <- NULL
  if (structure %in% c("isotropic", "diagonal")) {
    weight_names <- if (structure == "isotropic") names(basis$columns) else {
      colnames(basis$R)
    }
    weights <- matrix(NA_real_, measurements, length(weight_names),
      dimnames = list(NULL, weight_names))
  }
  for (i in seq_len(measurements)) {
    S <- .unsvec_symmetric(packed[i, ], q)
    addressable[[i]] <- sum(diag(S))
    one <- .model_geometry_fit(S, structure, model, fit_ranks, tolerance,
      max_sweeps)
    coefficients[, , i] <- one$C
    fit$residual[[i]] <- one$residual
    fit$fitted_energy[[i]] <- one$fitted_energy
    fit$iterations[[i]] <- one$iterations
    fit$converged[[i]] <- one$converged
    fit$clipped_negative_mass[[i]] <- one$clipped_negative_mass
    fit$truncated_positive_mass[[i]] <- one$truncated_positive_mass
    if (!is.null(weights)) weights[i, ] <- one$weights
  }
  if (structure != "shared") {
    fit$clipped_negative_mass <- NULL
    fit$truncated_positive_mass <- NULL
  }
  scale <- max(1, max(abs(centered_total)))
  identity_error <- max(abs(addressable + orthogonal - centered_total))
  if (identity_error > 1e-10 * scale) {
    .invariant_error(sprintf(paste0(
      "The trace decomposition does not close: max |addressable + orthogonal ",
      "- centered total| = %.3g."
    ), identity_error))
  }

  cross_fit <- if (is.null(policy)) NULL else {
    .model_geometry_cross_fit(x, lowered, lowered_plan, ranks$effective,
      component, policy)
  }
  scientific_plan_id <- .sha256_signature(list(
    schema_version = 1L,
    role = "model_geometry",
    source = x$scientific_plan_id,
    lowered = lowered_plan$scientific_plan_id,
    basis = basis$signature,
    structure = structure,
    rank = ranks$effective,
    component = component,
    training = if (is.null(policy)) NULL else policy$signature
  ), "model-geometry-sha256:")
  receipt <- .projection_receipt(compressed$receipt, scientific_plan_id)
  .new_effect_model_geometry(
    structure = structure, ranks = ranks, component = component,
    addressable = addressable, orthogonal = orthogonal,
    centered_total = centered_total, identity_error = identity_error,
    fit = fit, coefficients = coefficients, weights = weights,
    training = if (is.null(policy)) NULL else policy$kind,
    cross_fit = cross_fit,
    df = .model_geometry_df(structure, basis, ranks), basis = basis,
    near_saturated = near_saturated, warning_fraction = warning_fraction,
    index = index, receipt = receipt,
    source_scientific_plan_id = x$scientific_plan_id,
    lowered_scientific_plan_id = lowered_plan$scientific_plan_id
  )
}

.model_geometry_fields <- c(
  "structure", "hypothesis", "tests_supplied_geometries", "rank", "component",
  "addressable", "orthogonal", "centered_total", "identity_error", "fit",
  "coefficients", "weights", "training", "cross_fit", "df", "model_sizes",
  "span_fraction", "near_saturated",
  "warning_fraction", "basis_signature", "basis_models", "conditions",
  "index", "receipt", "source_scientific_plan_id",
  "lowered_scientific_plan_id", "orthogonal_source", "layer", "reading",
  "contract_signature"
)

# The signature covers every value the reader reports, not only the claims
# about them, so an edited energy, weight, coefficient or degree of freedom
# is a signature mismatch even where no other field contradicts it.
.model_geometry_contract_signature <- function(x) {
  .sha256_signature(list(
    schema_version = 1L,
    result_capability = "model_coordinate_geometry",
    structure = x$structure,
    rank = x$rank,
    component = x$component,
    measurements = length(x$addressable),
    basis_signature = x$basis_signature,
    training = x$training,
    scientific_plan_id = x$receipt$scientific_plan_id,
    values = .sha256_signature(list(
      addressable = x$addressable, orthogonal = x$orthogonal,
      centered_total = x$centered_total, identity_error = x$identity_error,
      fit = x$fit, coefficients = x$coefficients, weights = x$weights,
      df = x$df, model_sizes = x$model_sizes, span_fraction = x$span_fraction,
      near_saturated = x$near_saturated, warning_fraction = x$warning_fraction,
      cross_fit = x$cross_fit, index = x$index
    ))
  ))
}

.new_effect_model_geometry <- function(structure, ranks, component,
                                       addressable, orthogonal,
                                       centered_total, identity_error, fit,
                                       coefficients, weights, training,
                                       cross_fit, df, basis,
                                       near_saturated, warning_fraction,
                                       index, receipt,
                                       source_scientific_plan_id,
                                       lowered_scientific_plan_id) {
  value <- structure(list(
    structure = structure,
    hypothesis = .model_geometry_hypothesis(structure, length(basis$columns)),
    tests_supplied_geometries = identical(structure, "isotropic"),
    rank = ranks[c("requested", "effective", "clamped")],
    component = component,
    addressable = addressable,
    orthogonal = orthogonal,
    centered_total = centered_total,
    identity_error = identity_error,
    fit = fit,
    coefficients = coefficients,
    weights = weights,
    training = training,
    cross_fit = cross_fit,
    df = df,
    model_sizes = stats::setNames(as.integer(lengths(basis$columns)),
      names(basis$columns)),
    span_fraction = basis$span_fraction,
    near_saturated = near_saturated,
    warning_fraction = warning_fraction,
    basis_signature = basis$signature,
    basis_models = basis$models$model,
    conditions = basis$conditions,
    index = index,
    receipt = receipt,
    source_scientific_plan_id = source_scientific_plan_id,
    lowered_scientific_plan_id = lowered_scientific_plan_id,
    orthogonal_source = "original plan, operator H - P",
    layer = list(
      estimation = c("addressable", "orthogonal", "centered_total",
        if (!is.null(cross_fit)) "cross_fit"),
      latent = c("fit", "coefficients", "weights")
    ),
    reading = .latent_reading_line
  ), class = "effect_model_geometry")
  value$contract_signature <- .model_geometry_contract_signature(value)
  .validate_effect_model_geometry(value)
  value
}

.validate_effect_model_geometry <- function(x) {
  if (!.sealed_fields(x, "effect_model_geometry", .model_geometry_fields)) {
    .input_error("`x` must be a canonical `effect_model_geometry`.")
  }
  n_models <- length(x$basis_models)
  if (!x$structure %in% .model_geometry_structures ||
      !identical(x$hypothesis,
        .model_geometry_hypothesis(x$structure, n_models)) ||
      !identical(x$tests_supplied_geometries, identical(x$structure, "isotropic")) ||
      !identical(x$reading, .latent_reading_line) ||
      !identical(x$orthogonal_source, "original plan, operator H - P") ||
      !identical(x$layer, list(
        estimation = c("addressable", "orthogonal", "centered_total",
          if (!is.null(x$cross_fit)) "cross_fit"),
        latent = c("fit", "coefficients", "weights")))) {
    .contract_error("A model geometry misstates what its structure tests.")
  }
  if (!is.integer(x$model_sizes) || !identical(names(x$model_sizes), x$basis_models) ||
      !identical(sum(x$model_sizes), dim(x$coefficients)[[1L]]) ||
      !.is_number(x$span_fraction) || !.is_number(x$warning_fraction) ||
      !identical(x$near_saturated, x$span_fraction >= x$warning_fraction)) {
    .contract_error("A model geometry's saturation label is not the one its span fraction gives.")
  }
  q <- dim(x$coefficients)[[1L]]
  n <- length(x$conditions)
  psd_df <- function(d, s) d * s - s * (s - 1) / 2
  derived_df <- switch(x$structure,
    isotropic = list(fit = as.numeric(n_models), total = q),
    diagonal = list(fit = as.numeric(sum(x$model_sizes)), total = q),
    block = list(fit = as.numeric(sum(psd_df(x$model_sizes, x$rank$effective))),
      total = sum(x$rank$effective)),
    shared = list(fit = as.numeric(psd_df(q, x$rank$effective)),
      total = x$rank$effective)
  )
  derived_df$free <- as.numeric(psd_df(n - 1L, min(derived_df$total, n - 1L)))
  if (!identical(x$df, derived_df[c("fit", "free")])) {
    .contract_error("A model geometry's degrees of freedom are not its structure's.")
  }
  if (x$structure %in% c("isotropic", "diagonal")) {
    if (!is.matrix(x$weights) || nrow(x$weights) != length(x$addressable)) {
      .input_error("A nonnegative structure fit must carry its weights.")
    }
    expanded <- if (x$structure == "isotropic") {
      t(apply(x$weights, 1, function(w) rep(w, times = x$model_sizes)))
    } else {
      x$weights
    }
    diagonals <- t(vapply(seq_len(length(x$addressable)), function(i) {
      diag(matrix(x$coefficients[, , i], q, q))
    }, numeric(q)))
    if (max(abs(unname(expanded) - unname(diagonals))) > 1e-10 * max(1, max(abs(expanded)))) {
      .contract_error("A model geometry's weights are not the diagonal of its coefficients.")
    }
  } else if (!is.null(x$weights)) {
    .input_error("Only the isotropic and diagonal structures carry weights.")
  }
  measurements <- length(x$addressable)
  if (!.is_finite_numeric(x$addressable) || !.is_finite_numeric(x$orthogonal) ||
      !.is_finite_numeric(x$centered_total) ||
      length(x$orthogonal) != measurements ||
      length(x$centered_total) != measurements ||
      NROW(x$index) != measurements || !is.data.frame(x$fit) ||
      nrow(x$fit) != measurements || !is.array(x$coefficients) ||
      !identical(dim(x$coefficients)[[3L]], measurements)) {
    .input_error("Model-geometry fields have inconsistent shapes.")
  }
  scale <- max(1, max(abs(x$centered_total)))
  closure <- max(abs(x$addressable + x$orthogonal - x$centered_total))
  if (closure > 1e-10 * scale || !isTRUE(all.equal(closure, x$identity_error))) {
    .contract_error(
      "A model geometry's trace decomposition does not close as it claims."
    )
  }
  if (any(x$fit$fitted_energy < -1e-10 * scale) || any(x$fit$residual < 0)) {
    .contract_error("A model geometry's fitted energy or residual is negative.")
  }
  if (!is.null(x$weights) && any(x$weights < 0)) {
    .contract_error("A nonnegative structure fit carries a negative weight.")
  }
  if (xor(is.null(x$training), is.null(x$cross_fit)) ||
      (!is.null(x$training) && !identical(x$training, "exclude_evaluation"))) {
    .contract_error("A model geometry's training policy and cross-fit disagree.")
  }
  if (!is.null(x$cross_fit)) {
    cf <- x$cross_fit
    if (!is.list(cf) || !.is_finite_numeric(cf$energy) ||
        length(cf$energy) != measurements || !is.matrix(cf$by_edge) ||
        nrow(cf$by_edge) != measurements || !is.data.frame(cf$edges) ||
        nrow(cf$edges) != ncol(cf$by_edge) ||
        !isTRUE(all.equal(sum(cf$edges$weight), 1, tolerance = 1e-12)) ||
        !identical(cf$policy, "exclude_evaluation") ||
        !identical(cf$layer, "estimation") ||
        !"cross_fit" %in% x$layer$estimation) {
      .input_error("Model-geometry cross-fit fields are inconsistent.")
    }
    reduced <- as.numeric(cf$by_edge %*% cf$edges$weight)
    if (max(abs(reduced - cf$energy)) > 1e-10 * max(1, max(abs(cf$energy)))) {
      .contract_error(
        "A cross-fitted energy is not the weighted reduction of its edges."
      )
    }
    overlapping <- vapply(seq_len(nrow(cf$edges)), function(e) {
      trained <- strsplit(cf$edges$trained_on[[e]], ",", fixed = TRUE)[[1L]]
      any(c(cf$edges$left[[e]], cf$edges$right[[e]]) %in% trained) ||
        length(trained) < 2L
    }, logical(1))
    if (!identical(cf$rank, x$rank$effective) ||
        any(cf$edges$left == cf$edges$right) || any(overlapping) ||
        !is.integer(cf$edges$training_edges) || any(cf$edges$training_edges < 1L) ||
        !.is_string(cf$independence)) {
      .contract_error(
        "A cross-fit's edge bookkeeping contradicts the disjoint-edge discipline."
      )
    }
  }
  .validate_execution_receipt(x$receipt)
  .check_signature(x$contract_signature, .model_geometry_contract_signature(x),
    "Model-geometry contract signature is inconsistent with its claims.")
  invisible(x)
}

#' @export
as.data.frame.effect_model_geometry <- function(x, row.names = NULL,
                                                optional = FALSE, ...) {
  .validate_effect_model_geometry(x)
  columns <- cbind(
    addressable = x$addressable, orthogonal = x$orthogonal,
    centered_total = x$centered_total,
    fit_residual = x$fit$residual, fitted_energy = x$fit$fitted_energy,
    converged = x$fit$converged
  )
  if (!is.null(x$cross_fit)) {
    columns <- cbind(columns, cross_fit_energy = x$cross_fit$energy)
  }
  if (!is.null(x$weights)) {
    weights <- x$weights
    colnames(weights) <- paste0("weight_", colnames(weights))
    columns <- cbind(columns, weights)
  }
  .bind_result_values(x$index, columns)
}

#' @export
format.effect_model_geometry <- function(x, ...) {
  .format_counted_result("effect_model_geometry", length(x$addressable),
    sprintf("%s fit over %s", x$structure,
      .msg_count(length(x$basis_models), "model")))
}

.model_geometry_rank_line <- function(x) {
  if (is.null(x$rank$effective)) return(NULL)
  effective <- x$rank$effective
  text <- if (length(effective) > 1L) {
    paste(sprintf("%s %d", names(effective), effective), collapse = ", ")
  } else {
    as.character(effective)
  }
  paste0(text, if (isTRUE(x$rank$clamped)) " (clamped)" else "")
}

.model_geometry_signed_line <- function(values) {
  sprintf("median %s (range %s to %s)", .pf_num(stats::median(values), 4L),
    .pf_num(min(values), 4L), .pf_num(max(values), 4L))
}

#' @export
print.effect_model_geometry <- function(x, ...) {
  .validate_effect_model_geometry(x)
  unconverged <- sum(!x$fit$converged)
  .format_result_preview(x, "effect_model_geometry", fields = list(
    component = x$component,
    structure = .pf_wrap(c(x$structure, paste0("tests ", x$hypothesis))),
    basis = sprintf("%s (%s); span fraction %s%s",
      .msg_count(length(x$basis_models), "model"),
      paste(x$basis_models, collapse = ", "), .pf_num(x$span_fraction, 3L),
      if (isTRUE(x$near_saturated)) ", near saturated" else ""),
    rank = .model_geometry_rank_line(x),
    addressable = .pf_wrap(c(paste0("signed tr(S_x), ",
      .model_geometry_signed_line(x$addressable)),
      paste0("orthogonal ", .model_geometry_signed_line(x$orthogonal)),
      sprintf("split closes to %s", .pf_num(x$identity_error, 2L)))),
    cross_fit = if (is.null(x$cross_fit)) NULL else .pf_wrap(c(
      paste0("signed rank-", x$cross_fit$rank, " energy, ",
        .model_geometry_signed_line(x$cross_fit$energy)),
      sprintf("readout learned on edges disjoint from each of %s",
        .msg_count(nrow(x$cross_fit$edges), "evaluation edge")),
      "exclude_evaluation; estimation layer, not the plug-in fit"
    )),
    fit = .pf_wrap(c(
      paste0("fitted energy ", .model_geometry_signed_line(x$fit$fitted_energy)),
      paste0("residual ", .model_geometry_signed_line(x$fit$residual)),
      sprintf("df %s of %s free", .pf_num(x$df$fit, 3L), .pf_num(x$df$free, 3L)),
      if (unconverged) sprintf("%d not converged", unconverged)
    )),
    reading = .pf_wrap(c(.latent_reading_line,
      "applies to fit, coefficients and weights"))
  ), ...)
  cat(sprintf("  %-14s%s\n", "next:", "as.data.frame(x), x$coefficients, x$fit"))
  invisible(x)
}
