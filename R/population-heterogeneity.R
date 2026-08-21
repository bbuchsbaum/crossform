# Population heterogeneity ----------------------------------------------------
#
# Layer 5 (results / views). `heterogeneity()` is the reader that turns a
# planned or estimated population form into the `N`-by-`N` subject Gram of
# `design/population-form-contract.md` (`population-form-v1`) sections 5 and 6,
# its spectrum, its subject loadings, and --- at nodes the caller names --- the
# geometry of its modes.
#
# Its edges all point down or sideways: `population-driver.R` (layer 4) for the
# streaming tile, the group index, the packed coordinate table and the standard
# population receipt; `population-plan.R` and `geometry-plan.R` (layer 3) for
# the restricted-pairing re-plan the cross-fitted estimator needs; `latent.R`
# (layer 5, sideways) for the PSD projection and its moved-mass discipline; and
# `print-methods.R` / `format-results.R` (layer 5, sideways) for the printing
# primitives. Nothing in either printing file refers back here, which is what
# keeps the file call graph acyclic.
#
# What this file computes, and why in this shape:
#
#   * Per group node `u`, `K_u = D_u D_u^T`, where `D_u` is the `N`-by-`p`
#     matrix of subject deviations from the group fit at `u` in `svec`
#     coordinates. `K = sum_u K_u` is the global Gram, and the sum is the whole
#     reason it streams: section 5 measures node-additivity at `3.55e-15`, so
#     the accumulator is `N`-by-`N` while the `N`-by-`m`-by-`p` stack is never
#     built.
#   * `D` is the residual of the group model, not merely the deviation from the
#     mean. Under the default `~ 1` the two coincide, and under any other
#     design the residual is the one that answers "heterogeneity beyond what
#     the covariates already explain". The divisor is the residual degrees of
#     freedom, which is `N - 1` for the intercept-only case the contract's
#     oracle uses.
#   * The `sqrt(2)` off-diagonal packing of `symmetric_packed` is load-bearing
#     and not cosmetic. Euclidean inner products of packed rows are Frobenius
#     inner products of forms only with it (section 5, measured `3.11e-15`
#     against `2.96e+00` for a naive `upper.tri` packing), so a Gram built on
#     the naive codec is not a Gram of geometries at all. `crossform`'s codec
#     is already correct; `tests/testthat/test-population-heterogeneity.R`
#     pins it here so a future storage change cannot break this layer silently.
#   * The plug-in Gram is biased. Within-subject sampling noise lands entirely
#     in the heterogeneity term, measured at **+62.7 %** of the true
#     between-subject trace (section 6.2). The cross-fitted Gram of section 6.4
#     removes it exactly, by taking every inner product across two independent
#     partition halves --- including the diagonal, which is what makes it
#     unbiased for the *noiseless* Gram rather than merely for its off-diagonal
#     part.
#
# The sink is not in the Gram. It is a row of a population result and it is
# always in budget units, while the group nodes are in the plan's declared
# semantics; an inner product mixing the two would be adding a mass to a
# density. It is excluded, the exclusion is recorded, and the sink budget stays
# readable on the receipt where a differentially failing transport is supposed
# to be visible.

# The two estimators, and the two coordinate axes a Gram can be accumulated
# along. `cross_fit` is first because it is the one section 6.4 makes
# normative; `plug_in` exists because a caller may only have one partition
# split, and because the contract admits it for directions.
.heterogeneity_estimators <- c("cross_fit", "plug_in")
.heterogeneity_spaces <- c("packed_form", "query_bank")

# Section 6.2's measured plug-in inflation, quoted rather than recomputed. It
# is a property of the estimator and of the design that produced the oracle's
# fixture, so it travels on the receipt as a *reference* and never as this
# run's own bias, which nothing here can know.
.heterogeneity_plug_in_reference <- list(
  measured_inflation = 0.627,
  quantity = "trace of the heterogeneity term",
  source = paste0(
    "design/population-form-contract.md section 6.2; ",
    "design/oracles/population-geometry-split.R section P6.b"
  )
)

.heterogeneity_estimator_phrase <- function(estimator) {
  switch(estimator,
    cross_fit = "cross-partition inner products (section 6.4)",
    plug_in = "same-partition inner products (biased, section 6.2)",
    estimator
  )
}

# The estimator's own refusals -------------------------------------------------

# Section 6.4's construction needs two *independent* halves, and independence
# here is a statement about partitions: two halves that share a run share its
# noise, and an inner product across them is no longer cross-fitted. So a half
# is a set of partitions, the halves are disjoint, and each one needs at least
# two partitions of its own --- a half with one partition has no cross-
# generalized estimate inside it, and `cross_partitions()` would refuse it.
# Four partitions is therefore the floor for the cross-fitted route.
.heterogeneity_partition_refusal <- function(reasons, remedies) {
  .capability_refusal(paste0(
    "The cross-fitted subject Gram needs two independent halves of each ",
    "participant's partitions, and each half needs at least two partitions ",
    "so that it carries a cross-generalized estimate of its own. Splitting a ",
    "participant's partitions any finer would put the same run on both sides ",
    "of an inner product, which is exactly the within-subject noise the ",
    "cross-fit exists to remove -- `population-form-v1` section 6.4 measures ",
    "the plug-in inflation it would reintroduce at +62.7%."
  ),
    capability = "cross_fitted_subject_gram",
    namespace = "population_heterogeneity",
    reasons = reasons,
    remedies = c(remedies,
      paste0(
        "Read the plug-in Gram instead with `estimator = \"plug_in\"`. Its ",
        "loading directions are admissible and named as such; its spectrum ",
        "is inflated by within-subject noise and carries no `n_eff`."
      )))
}

# `plan_geometry()` derives the whole plan from the relation, the frame and the
# pairing, so a restricted pairing is a new plan over the same source. What must
# not change is everything else the plan fixes: a re-plan that lands on a
# different metric schedule, codec or packed width is measuring something else,
# and the halves would not be comparable. The check is written against the
# schedule's own signature rather than against a list of fields, so a future
# schedule kind cannot slip through by adding one.
.heterogeneity_replan <- function(subject, partitions, label) {
  schedule <- subject$metric_schedule
  if (!schedule$kind %in%
      c("implicit_identity_before_frame", "fixed_metric_before_frame")) {
    .heterogeneity_partition_refusal(
      reasons = c(
        paste0("metric_schedule_not_repartitionable:", label),
        paste0("metric_schedule_kind:", schedule$kind)
      ),
      remedies = paste0(
        "A whitened or learned metric schedule freezes coordinates or ",
        "statistics from the whole partition set, so a half of it is not the ",
        "same measurement. Re-plan the participants on ",
        "`composition = \"native\"` with a fixed metric if the cross-fitted ",
        "Gram is the estimand you want."
      ))
  }
  # The half inherits the declared semantics of the pairing it restricts: a
  # half that silently declared independence the original did not would put a
  # claim on the estimand that nobody made. `pairing()` stores an undeclared
  # axis as `"undeclared"` and takes it back as `NULL`.
  independence <- attr(subject$pairing, "independence", exact = TRUE)
  if (identical(independence, "undeclared")) independence <- NULL
  over <- cross_partitions(partitions, independence = independence,
    generalizes_over = attr(subject$pairing, "generalizes_over", exact = TRUE))
  restricted <- plan_geometry(subject$task$left_relation, subject$frame, over,
    compute = subject$compute, metric = schedule$metric)
  if (!identical(restricted$metric_schedule$signature, schedule$signature) ||
      !identical(restricted$codec, subject$codec) ||
      !identical(restricted$packed_width, subject$packed_width) ||
      !identical(restricted$measurements, subject$measurements)) {
    .contract_error(sprintf(paste0(
      "Restricting participant `%s` to a partition half changed the plan ",
      "beyond its pairing: the two halves would measure different things and ",
      "their inner product would not be a cross-fit. This is a crossform ",
      "bug, not a data problem."
    ), label))
  }
  restricted
}

# The split. `NULL` interleaves --- odd-numbered partitions to the first half,
# even to the second --- because runs are usually acquired in order and a
# contiguous split confounds the two halves with anything that drifts across a
# session. A caller who knows better names the halves.
.heterogeneity_partition_split <- function(plan, partitions) {
  labels <- names(plan$subjects)
  available <- lapply(plan$subjects,
    function(subject) subject$task$left_relation$partitions)
  if (is.null(partitions)) {
    split <- lapply(available, function(value) list(
      a = value[seq.int(1L, length(value), by = 2L)],
      b = value[seq.int(2L, length(value), by = 2L)]
    ))
    kind <- "interleaved_partitions"
  } else {
    if (!is.list(partitions) || length(partitions) != 2L ||
        !all(vapply(partitions, .is_strings, logical(1), unique = TRUE))) {
      .input_error(paste0(
        "`partitions` names the two halves the cross-fit reads across: pass ",
        "a list of two character vectors of partition identifiers, or NULL ",
        "to interleave each participant's own partitions."
      ), arg = "partitions", received = .msg_value(partitions),
        expected = "a list of two character vectors of partition identifiers")
    }
    halves <- list(a = as.character(partitions[[1L]]),
      b = as.character(partitions[[2L]]))
    shared <- intersect(halves$a, halves$b)
    if (length(shared)) {
      .heterogeneity_partition_refusal(
        reasons = paste0("partition_in_both_halves:", shared),
        remedies = paste0(
          "Give each partition to exactly one half; ",
          .msg_names(shared), " appears in both."
        ))
    }
    missing_from <- unlist(lapply(labels, function(label) {
      absent <- setdiff(c(halves$a, halves$b), available[[label]])
      if (length(absent)) paste0(label, ":", absent) else NULL
    }), use.names = FALSE)
    if (length(missing_from)) {
      .heterogeneity_partition_refusal(
        reasons = paste0("partition_absent_for_subject:", missing_from),
        remedies = paste0(
          "Name only partitions every participant carries, or pass NULL and ",
          "let each participant's own partitions be interleaved."
        ))
    }
    split <- stats::setNames(rep(list(halves), length(labels)), labels)
    kind <- "declared"
  }
  short <- unlist(lapply(labels, function(label) {
    sizes <- lengths(split[[label]])
    if (any(sizes < 2L)) {
      paste0(label, ":", sum(sizes), "_partitions_split_",
        paste(sizes, collapse = "_"))
    } else {
      NULL
    }
  }), use.names = FALSE)
  if (length(short)) {
    .heterogeneity_partition_refusal(
      reasons = c("half_carries_fewer_than_two_partitions",
        paste0("subject_split:", short)),
      remedies = paste0(
        "Plan the participants over at least four partitions, so that each ",
        "half holds a cross-generalized estimate."
      ))
  }
  list(split = split, kind = kind)
}

# One half's population plan: the same participants, the same transports, the
# same group model, over a restricted pairing. It is built with the public
# planner rather than by editing the sealed plan, so every admission the group
# estimand rests on is re-checked and the half gets its own scientific identity.
.heterogeneity_half_plan <- function(plan, split, side) {
  labels <- names(plan$subjects)
  subjects <- stats::setNames(lapply(labels, function(label) {
    .heterogeneity_replan(plan$subjects[[label]], split[[label]][[side]], label)
  }), labels)
  half <- plan_population(subjects, plan$transport,
    model = plan$model$formula, data = plan$data,
    normalization = plan$normalization, compute = plan$compute,
    allow_nonconservative = plan$allow_nonconservative)
  if (!identical(half$model$matrix, plan$model$matrix)) {
    .contract_error(paste0(
      "The two partition halves must be fitted under one group design; ",
      "re-planning produced a different model matrix. This is a crossform ",
      "bug, not a data problem."
    ))
  }
  half
}

# A functional transport that was itself learned from the data has seen some
# partitions, and if those overlap a half then that half's inner product is not
# fully held out. crossform cannot repair it, and it is not an error --- the
# transport is a declared input --- so it is measured and recorded, and the
# print says so.
.heterogeneity_transport_overlap <- function(plan, split) {
  labels <- names(plan$subjects)
  overlap <- unlist(lapply(labels, function(label) {
    declared <- plan$transport[[label]]$provenance$cross_fit
    if (is.null(declared)) return(NULL)
    shared <- intersect(as.character(declared),
      unlist(split[[label]], use.names = FALSE))
    if (length(shared)) paste0(label, ":", shared) else NULL
  }), use.names = FALSE)
  if (is.null(overlap)) character() else overlap
}

# The Gram ---------------------------------------------------------------------

# Deviations from the group fit, with unresolved columns held out rather than
# solved through. A column carrying a non-finite response is a measurement:
# density semantics returns `NA` at a group node reached by no native mass, and
# it must not become a zero deviation, which would read as perfect agreement
# between participants at a node nobody measured.
.heterogeneity_residual <- function(model, block) {
  finite <- apply(is.finite(block), 2L, all)
  residuals <- matrix(0, nrow(block), ncol(block))
  if (any(finite)) {
    residuals[, finite] <- qr.resid(model$qr, block[, finite, drop = FALSE])
  }
  list(residuals = residuals, finite = finite)
}

# The streaming loop, shared by both passes and by both estimators. `plans` is
# one plan for the plug-in route and two --- the partition halves --- for the
# cross-fitted one; `visit` sees the kept columns of each plan's tile and
# accumulates whatever the pass is accumulating. Nothing of size
# `N x m x p` exists at any point: the largest array in flight is the group
# stack for one tile, `N x (m + 1) x tile`, one per plan.
.heterogeneity_stream <- function(plans, component, width, coordinate_labels,
                                  nodes, keep, starts, tile, visit) {
  states <- lapply(plans, function(plan)
    .population_tile_state(plan, width, coordinate_labels))
  for (start in starts) {
    columns <- start:min(start + tile - 1L, width)
    take <- rep(keep, times = length(columns))
    blocks <- vector("list", length(plans))
    for (side in seq_along(plans)) {
      streamed <- .population_stream_tile(plans[[side]], component, columns,
        width, coordinate_labels, nodes, states[[side]])
      states[[side]] <- streamed$state
      blocks[[side]] <- streamed$stack[, take, drop = FALSE]
      rm(streamed)
    }
    visit(columns, blocks)
    rm(blocks)
  }
  states
}

# Pass one: `K = sum_u K_u`, accumulated tile by tile.
#
# For one plan this is `R R^T`, the plug-in Gram of the residual forms. For two
# it is `(R_A R_B^T + R_B R_A^T) / 2`, section 6.4's `Gamma-hat` with the group
# residual maker applied on the subject axis. The residual maker is a fixed
# symmetric idempotent operator on that axis, so it passes through the
# expectation: `E[M Gamma-hat M] = M Gamma_true M`, and the cross-fit's
# unbiasedness survives the group fit unchanged.
.heterogeneity_gram <- function(plans, model, component, width,
                                coordinate_labels, nodes, keep, starts, tile) {
  labels <- names(plans[[1L]]$subjects)
  gram <- matrix(0, length(labels), length(labels),
    dimnames = list(labels, labels))
  unresolved <- 0L
  states <- .heterogeneity_stream(plans, component, width, coordinate_labels,
    nodes, keep, starts, tile, function(columns, blocks) {
      residual <- lapply(blocks, .heterogeneity_residual, model = model)
      finite <- Reduce(`&`, lapply(residual, `[[`, "finite"))
      unresolved <<- unresolved + sum(!finite)
      kept <- lapply(residual, function(value)
        value$residuals[, finite, drop = FALSE])
      gram <<- gram + if (length(kept) == 1L) {
        tcrossprod(kept[[1L]])
      } else {
        cross <- tcrossprod(kept[[1L]], kept[[2L]])
        (cross + t(cross)) / 2
      }
    })
  # Symmetric by construction; symmetrized anyway so that `eigen(symmetric =)`
  # is reading the matrix this function claims to have built and not a
  # floating-point neighbour of it.
  list(gram = (gram + t(gram)) / 2, unresolved = unresolved, states = states)
}

# Pass two: the mode geometry, at the nodes the caller named.
#
# Section 5's mode recovery is `v = C^T u`: the coordinate-space mode is the
# loading-weighted combination of the subject deviations, so the `D x D`
# covariance is never formed. The full mode spans every retained group node;
# only the selected nodes' slices are kept, and the norm is accumulated over
# all of them so that the slices returned are slices of one unit-norm mode and
# their relative sizes mean something.
#
# Under the cross-fitted estimator the deviations of the two halves are
# averaged before the contraction. That is a plug-in reconstruction of a
# *direction*, which section 6.3 admits explicitly and requires be named; the
# eigenvalues beside it still come from the cross-fitted Gram.
.heterogeneity_modes <- function(plans, model, component, width,
                                 coordinate_labels, nodes, keep, starts, tile,
                                 loadings, selected, node_labels,
                                 coordinate_names) {
  n_mode <- ncol(loadings)
  n_kept <- sum(keep)
  forms <- array(0, c(length(selected), n_mode, width),
    dimnames = stats::setNames(
      list(node_labels[selected], colnames(loadings), coordinate_names),
      c("node", "mode", "coordinate")
    ))
  # Which of the selected `(node, coordinate)` cells the Gram was allowed to
  # see. A cell held out of the accumulation contributes nothing to the mode
  # and must come back absent, not zero: zero is a form, and it would read as
  # "this mode has no geometry here" at a node nobody measured.
  resolved <- matrix(FALSE, length(selected), width)
  norms <- numeric(n_mode)
  .heterogeneity_stream(plans, component, width, coordinate_labels, nodes,
    keep, starts, tile, function(columns, blocks) {
      residual <- lapply(blocks, .heterogeneity_residual, model = model)
      finite <- Reduce(`&`, lapply(residual, `[[`, "finite"))
      deviations <- Reduce(`+`, lapply(residual, `[[`, "residuals")) /
        length(residual)
      deviations[, !finite] <- 0
      # `(n_kept * length(columns))` by `n_mode`, the mode's coordinates on
      # this tile.
      weighted <- crossprod(deviations, loadings)
      norms <<- norms + colSums(weighted^2)
      # The tile stacked node fastest, so row `(c - 1) * n_kept + u` is node
      # `u` of the tile's `c`th coordinate.
      pick <- as.integer(outer(selected,
        (seq_along(columns) - 1L) * n_kept, `+`))
      resolved[, columns] <<- matrix(finite[pick], length(selected),
        length(columns))
      for (mode in seq_len(n_mode)) {
        forms[, mode, columns] <<- matrix(weighted[pick, mode],
          length(selected), length(columns))
      }
    })
  .heterogeneity_scale_modes(forms, resolved, sqrt(norms))
}

# One place where a reconstructed mode is turned into the unit-norm direction
# that is reported, and where the two absences are written in. The norm is
# taken over every retained group node, not only the selected ones, so a
# selected node's slice carries its share of one mode rather than being
# renormalized into a mode of its own.
.heterogeneity_scale_modes <- function(forms, resolved, scale) {
  for (mode in seq_along(scale)) {
    values <- forms[, mode, ]
    # Held-out cells first: they contributed zero to the norm above, so
    # blanking them here leaves the unit-norm accounting exactly as it was.
    values[!resolved] <- NA_real_
    forms[, mode, ] <- if (scale[[mode]] > 0) {
      values / scale[[mode]]
    } else {
      # A mode with no mass at the retained nodes has no direction to report.
      NA_real_
    }
  }
  forms
}

# Eigenvector signs are arbitrary, and an arbitrary sign that changes between
# runs makes two identical analyses look different. The convention is fixed
# here and applied before the mode forms are reconstructed, so a loading and
# the geometry it weights always carry the same sign.
.heterogeneity_orient <- function(vectors) {
  for (mode in seq_len(ncol(vectors))) {
    column <- vectors[, mode]
    if (column[[which.max(abs(column))]] < 0) vectors[, mode] <- -column
  }
  vectors
}

# The latent layer -------------------------------------------------------------

# Section 6.5, which is `conservative-geometry-v1` section 6 applied one level
# up: the heterogeneity spectrum is signed, and every functional that treats it
# as a nonnegative partition -- a variance-explained fraction, a cumulative
# curve, an effective number of modes -- is defined only on a named projection
# of it. `latent_geometry()`'s own operators and its moved-mass accounting are
# reused rather than re-derived, so the two layers cannot drift on what a PSD
# projection costs.
.heterogeneity_latent <- function(spectrum, method) {
  projected <- .latent_psd_projection(matrix(spectrum, nrow = 1L))
  list(
    method = method,
    operator = unname(.latent_projection_operators[[method]]),
    spectrum = as.numeric(projected$spectrum[1L, ]),
    cumulative = as.numeric(projected$cumulative[1L, ]),
    n_eff = projected$n_eff[[1L]],
    moved_mass = projected$moved_mass[[1L]],
    moved_share = projected$moved_share[[1L]],
    negative_modes = .heterogeneity_negative_modes(spectrum)
  )
}

# How many modes are negative *as a statement about the matrix* rather than
# about the eigensolver. The plug-in Gram is `R R^T` and is positive
# semidefinite by construction, so its trailing eigenvalues are numerical
# zeros that land on either side of zero; counting one of those as negative
# would report an indefiniteness the estimator cannot have. The cross-fitted
# Gram is genuinely indefinite (section 6.4) and its negative eigenvalues sit
# far outside this band. The band is used only to report a *count* and never
# to clip a value: the projection in `.latent_psd_projection()` stays exact,
# as `conservative-geometry-v1` section 6 requires.
.heterogeneity_negative_modes <- function(spectrum) {
  sum(spectrum < -1e-10 * max(abs(spectrum), 0))
}

# Section 6.3 and section 6.4 are normative about which Gram may carry a
# derived number: loading *directions* may be read from either, provided the
# source is named, but any variance-explained figure, cumulative curve or
# `n_eff` must come from the cross-fitted Gram. The plug-in spectrum is
# reported --- withholding the eigenvalues of a matrix the record also carries
# would be theatre --- and it is reported beside the estimator's name, the
# measured inflation, and this refusal, which travels with the record exactly
# as the transported-covariance refusal travels with a population result.
.heterogeneity_latent_refusal <- function() {
  list(
    capability = "plug_in_spectrum_functionals",
    namespace = "population_heterogeneity",
    reasons = c(
      "within_subject_noise_booked_as_heterogeneity",
      "plug_in_trace_inflated_62_7_percent_on_the_contract_fixture"
    ),
    remedies = c(
      paste0(
        "Read `$loadings`, which section 6.3 admits from either Gram: the ",
        "plug-in spectrum is inflated roughly isotropically, so the leading ",
        "directions survive even where the eigenvalues do not."
      ),
      paste0(
        "Re-read the population with `estimator = \"cross_fit\"` on the ",
        "plan, which computes every inner product across two independent ",
        "partition halves and is unbiased for the noiseless subject Gram."
      )
    ),
    # Section 6.4's own list is "any spectrum, variance-explained figure or
    # `n_eff`", and the spectrum is the one item of the three crossform
    # reports anyway. The departure is stated here and in the Rd rather than
    # left for a reader to notice by comparing the two lists: the record
    # carries the Gram, the validator re-derives the spectrum from it, and
    # withholding a number two lines of arithmetic away from a field the
    # object publishes would be theatre rather than discipline. What is
    # withheld is everything the inflation actually distorts.
    message = paste0(
      "An effective number of modes, a variance-explained fraction and a ",
      "cumulative curve are refused on a plug-in heterogeneity spectrum. The ",
      "plug-in estimator books every participant's own within-subject ",
      "sampling noise as between-participant heterogeneity -- measured at ",
      "+62.7% of the true between-subject trace on the contract's own ",
      "fixture, which is a reference and not this run's own bias -- so the ",
      "fractions would be ",
      "fractions of an inflated whole, and the effective count would count ",
      "noise modes. The spectrum itself is reported, named as plug-in."
    )
  )
}

# Node selection ---------------------------------------------------------------

# Which group nodes get a reconstructed mode form. `NULL` is none, and it is
# the default because the reconstruction costs a second streaming pass over
# every participant; the Gram, the spectrum and the loadings come out of the
# first pass alone.
.heterogeneity_nodes <- function(nodes, node_labels) {
  if (is.null(nodes)) return(integer())
  if (is.numeric(nodes)) {
    positions <- suppressWarnings(as.integer(nodes))
    if (anyNA(positions) || any(positions != nodes) ||
        any(positions < 1L) || any(positions > length(node_labels))) {
      .input_error(sprintf(paste0(
        "`nodes` given as positions must index the %s the population ",
        "estimates; the sink is not among them."
      ), .msg_count(length(node_labels), "group node")),
        arg = "nodes", received = .msg_value(nodes),
        expected = sprintf("whole numbers in 1:%d", length(node_labels)))
    }
    return(sort(unique(positions)))
  }
  if (!.is_strings(nodes)) {
    .input_error(paste0(
      "`nodes` selects the group nodes whose mode geometry is reconstructed: ",
      "pass their identifiers, their positions, or NULL for none."
    ), arg = "nodes", received = .msg_value(nodes),
      expected = "group node identifiers, positions, or NULL")
  }
  absent <- setdiff(nodes, node_labels)
  if (length(absent)) {
    .capability_refusal(sprintf(paste0(
      "`nodes` names %s the population does not carry. The sink is not ",
      "selectable: it is unmapped territory in budget units, not a location, ",
      "and it carries no geometry to reconstruct."
    ), .msg_count(length(absent), "group node")),
      capability = "selected_group_nodes",
      namespace = "population_heterogeneity",
      reasons = paste0("group_node_absent:", absent),
      remedies = paste0(
        "Choose from the estimated group nodes: ", .msg_names(node_labels), "."
      ))
  }
  sort(unique(match(nodes, node_labels)))
}

# Identity ---------------------------------------------------------------------

# The record's own scientific identity. The estimator is in it because the
# cross-fitted and plug-in Grams are two estimands and not two renderings of
# one; the partition split is in it for the same reason section 6.4 requires it
# recorded; the node selection is *not*, because it decides what is displayed
# rather than what was estimated.
.heterogeneity_scientific_id <- function(parent, estimator, space, component,
                                         split, projection) {
  .sha256_signature(list(
    schema_version = 1L,
    contract_version = "population-form-v1",
    role = "population_heterogeneity",
    parent = parent,
    estimator = estimator,
    space = space,
    component = component,
    split = split,
    projection = projection
  ), "population-sha256:")
}

# The record -------------------------------------------------------------------
#
# FIELD CONTRACT --- `effect_population_heterogeneity`
#
#   $estimator     "cross_fit" or "plug_in"; read this before any spectrum
#   $space         "packed_form" (the streamed plan route, `svec` coordinates)
#                  or "query_bank" (the query-first route, one coordinate per
#                  declared contrast)
#   $gram          N x N, the residual subject Gram, divided by $residual_df
#   $spectrum      its eigenvalues, decreasing, named mode1..modeN; SIGNED, and
#                  routinely indefinite under `cross_fit` (section 6.2)
#   $loadings      N x N eigenvectors, subjects by mode, sign-oriented
#   $mode_forms    node x mode x coordinate, or NULL when no `nodes` were
#                  selected; every `[node, mode, ]` slice is one packed form
#                  under `$space == "packed_form"`
#   $nodes         the selected node identifiers, or NULL
#   $latent        the PSD-projected descriptive layer (n_eff, cumulative,
#                  moved mass), or NULL under `plug_in`, where section 6.4
#                  refuses it and `$receipt$latent_refusal` says so
#   $index         the retained group nodes; the sink is never among them
#   $effects       the q effect names, or NULL in query space
#   $coordinates   the axis the Gram was accumulated along
#   $receipt       $estimator, $gram, $cross_fit, $projection, $streaming, plus
#                  the standard population receipt blocks
.population_heterogeneity_fields <- c(
  "estimator", "space", "gram", "spectrum", "loadings", "mode_forms", "nodes",
  "latent", "subjects", "residual_df", "index", "effects", "coordinates",
  "component", "ledger", "semantics", "normalization", "receipt",
  "scientific_plan_id"
)

# The Gram's decomposition, taken once and shared by the record, the mode
# reconstruction and the receipt. Eigenvector signs are fixed here, so a
# loading and the geometry it weights can never disagree about which end of a
# mode is which.
.heterogeneity_decompose <- function(gram) {
  if (!all(is.finite(gram))) {
    .contract_error(paste0(
      "The subject Gram carries non-finite entries. Unresolved group nodes ",
      "are held out of the accumulation rather than summed through, so this ",
      "is a crossform bug and not a property of the data."
    ))
  }
  labels <- rownames(gram)
  modes <- paste0("mode", seq_along(labels))
  decomposition <- eigen(gram, symmetric = TRUE)
  loadings <- .heterogeneity_orient(decomposition$vectors)
  dimnames(loadings) <- list(labels, modes)
  list(spectrum = stats::setNames(decomposition$values, modes),
    loadings = loadings)
}

.new_population_heterogeneity <- function(estimator, space, gram, decomposition,
                                          latent, mode_forms, selected, index,
                                          effects, coordinates, component,
                                          semantics, normalization,
                                          residual_df, receipt, identity) {
  labels <- rownames(gram)
  spectrum <- decomposition$spectrum
  loadings <- decomposition$loadings
  value <- structure(list(
    estimator = estimator,
    space = space,
    gram = gram,
    spectrum = spectrum,
    loadings = loadings,
    mode_forms = mode_forms,
    nodes = if (length(selected)) as.character(index$node)[selected],
    latent = latent,
    subjects = labels,
    residual_df = as.integer(residual_df),
    index = index,
    effects = effects,
    coordinates = coordinates,
    component = component,
    ledger = unname(.population_ledger_names[[component]]),
    semantics = semantics,
    normalization = normalization,
    receipt = receipt,
    scientific_plan_id = identity
  ), class = "effect_population_heterogeneity")
  .validate_population_heterogeneity(value)
  value
}

.validate_population_heterogeneity <- function(x) {
  if (!inherits(x, "effect_population_heterogeneity")) {
    .input_error(sprintf(paste0(
      "Expected an `effect_population_heterogeneity` from `heterogeneity()`; ",
      "received %s."
    ), .msg_value(x)), arg = "x", received = .msg_value(x),
      expected = "an `effect_population_heterogeneity` from `heterogeneity()`")
  }
  if (!.sealed_fields(x, "effect_population_heterogeneity",
        .population_heterogeneity_fields) ||
      !.is_string(x$estimator) || !x$estimator %in% .heterogeneity_estimators ||
      !.is_string(x$space) || !x$space %in% .heterogeneity_spaces ||
      !is.matrix(x$gram) || !is.numeric(x$gram) ||
      !.is_strings(x$subjects, unique = TRUE) ||
      !is.integer(x$residual_df) || !is.data.frame(x$index) ||
      !is.data.frame(x$coordinates) ||
      !.is_string(x$component) || !.is_string(x$ledger) ||
      !.is_string(x$semantics) || !.is_string(x$normalization) ||
      !is.list(x$receipt) ||
      !.strong_sha256(sub("^population-", "", x$scientific_plan_id))) {
    .input_error(
      "Population-heterogeneity fields are missing or noncanonical."
    )
  }
  n <- length(x$subjects)
  if (!identical(dim(x$gram), c(n, n)) ||
      !identical(dimnames(x$gram), list(x$subjects, x$subjects)) ||
      length(x$spectrum) != n || !identical(dim(x$loadings), c(n, n)) ||
      !identical(rownames(x$loadings), x$subjects)) {
    .contract_error(paste0(
      "A heterogeneity record's Gram, spectrum and loadings must agree with ",
      "its participant list."
    ))
  }
  # The spectrum is a re-derivable function of the Gram, so it is re-derived
  # rather than trusted: a plausible spectrum beside a Gram that does not
  # produce it is exactly the kind of number this record exists to make
  # impossible. The tolerance is the eigensolver's, scaled by the Gram.
  expected <- sort(eigen(x$gram, symmetric = TRUE, only.values = TRUE)$values,
    decreasing = TRUE)
  scale <- max(1, max(abs(expected)))
  if (max(abs(unname(x$spectrum) - expected)) > 1e-8 * scale) {
    .contract_error(paste0(
      "A heterogeneity record's spectrum is not the spectrum of the Gram it ",
      "carries."
    ))
  }
  if (any(x$index$sink)) {
    .contract_error(paste0(
      "The sink is unmapped territory in budget units, not a location, and ",
      "it is never a row of a heterogeneity index."
    ))
  }
  if (identical(x$estimator, "plug_in") && !is.null(x$latent)) {
    .contract_error(paste0(
      "`population-form-v1` section 6.4 refuses a variance-explained figure ",
      "or an effective mode count from the plug-in Gram; this record carries ",
      "a latent layer built on one."
    ))
  }
  if (!is.null(x$mode_forms)) {
    if (!is.array(x$mode_forms) || length(dim(x$mode_forms)) != 3L ||
        !identical(names(dimnames(x$mode_forms)),
          c("node", "mode", "coordinate")) ||
        dim(x$mode_forms)[[1L]] != length(x$nodes) ||
        dim(x$mode_forms)[[3L]] != nrow(x$coordinates)) {
      .contract_error(paste0(
        "A heterogeneity record's mode forms do not agree with their ",
        "selected group nodes and coordinate table."
      ))
    }
  } else if (!is.null(x$nodes)) {
    .contract_error(paste0(
      "A heterogeneity record naming selected group nodes must carry the ",
      "mode forms reconstructed at them."
    ))
  }
  invisible(x)
}

# The receipt ------------------------------------------------------------------
#
# The standard population blocks, from the same builder the executors use, plus
# the four things a reader of *this* record would have to know: which estimator
# formed the Gram and how it was normalized, which partitions formed the two
# halves, what a PSD projection cost, and what the streaming bound actually was.
.heterogeneity_receipt <- function(base, estimator, space, gram, split,
                                   overlap, half_ids, latent, streaming,
                                   unresolved, nodes, divisor) {
  # The two fields this record redefines are dropped from the standard blocks
  # rather than shadowed: `c()` keeps both copies of a duplicated name and `$`
  # returns the first, so leaving them in would publish the executor's
  # streaming bound and its unresolved count under this run's names.
  c(base[setdiff(names(base), c("streaming", "unresolved_columns"))], list(
    estimator = estimator,
    estimator_reads = .heterogeneity_estimator_phrase(estimator),
    gram = list(
      space = space,
      # Section 6.1's naming hazard, written out: the exact identity takes a
      # `1/N` divisor and the estimator takes the residual degrees of freedom.
      # They are two objects, and the divisor is named rather than left to be
      # inferred from context.
      divisor = "residual_df",
      divisor_value = as.double(divisor),
      subjects = as.integer(nrow(gram)),
      retained_nodes = as.integer(nodes),
      sink_excluded = TRUE,
      sink_reason = "sink_is_budget_units_at_no_location",
      # Section 6.3 admits a plug-in *direction* provided the source is
      # named, and under the cross-fitted estimator the mode geometry is
      # reconstructed from the average of the two halves' deviations. The
      # eigenvalues beside it still come from the cross-fitted Gram, and this
      # field is what stops a reader inferring the second from the first.
      mode_reconstruction = if (identical(estimator, "cross_fit")) {
        "plug_in_average_of_partition_halves"
      } else {
        "plug_in_single_read"
      },
      packing = if (identical(space, "packed_form")) {
        "symmetric_packed_sqrt2_off_diagonal"
      } else {
        "declared_query_bank"
      },
      node_additive = TRUE,
      # Which read the standard blocks below describe. A cross-fitted run reads
      # every participant twice, and one set of per-participant receipts cannot
      # describe both; the budget certificate is the exception, and it is the
      # worst of the two halves rather than either one.
      receipts_from = if (identical(estimator, "cross_fit")) {
        "partition_half_a"
      } else {
        "single_read"
      },
      budget_certificate_scope = if (identical(estimator, "cross_fit")) {
        "worst_of_both_halves"
      } else {
        "single_read"
      }
    ),
    cross_fit = if (identical(estimator, "cross_fit")) {
      list(
        split = split$kind,
        partitions = lapply(split$split, function(halves)
          list(a = halves$a, b = halves$b)),
        half_plan_ids = half_ids,
        transport_partition_overlap = overlap
      )
    },
    plug_in_reference = if (identical(estimator, "plug_in")) {
      .heterogeneity_plug_in_reference
    },
    latent_refusal = if (identical(estimator, "plug_in")) {
      .heterogeneity_latent_refusal()
    },
    projection = latent,
    streaming = streaming,
    unresolved_columns = as.integer(unresolved)
  ))
}

# The verb ---------------------------------------------------------------------

#' Decompose a population form into consensus and heterogeneity
#'
#' `heterogeneity()` builds the subject Gram of
#' `design/population-form-contract.md` (`population-form-v1`) sections 5 and
#' 6: the \eqn{N \times N}{N x N} matrix of inner products between
#' participants' deviations from the group fit, its spectrum, the subject
#' loadings on its modes, and --- at group nodes the caller names --- the
#' geometry of those modes.
#'
#' Geometry-space covariance across `N` participants has rank at most
#' \eqn{N-1}{N-1}, so the whole spectrum lives in an \eqn{N \times N}{N x N}
#' matrix and the \eqn{D \times D}{D x D} covariance never has to be formed.
#' The Gram is additive over group nodes, so it accumulates while streaming and
#' the \eqn{N \times m \times p}{N x m x p} stack is never built either.
#'
#' @section What the deviations are:
#' At group node `u` the record holds `K_u = D_u D_u^T`, where `D_u` is the
#' \eqn{N \times p}{N x p} matrix of participant deviations *from the fitted
#' value under the group model*, in the packed `svec` coordinates of
#' `effect-form-v1` section 3. Under the default `~ 1` design that is each
#' participant minus the group mean; under any other design it is what is left
#' after the covariates. The global Gram is `sum_u K_u`, divided by the residual
#' degrees of freedom.
#'
#' The `sqrt(2)` on the packed off-diagonals is load-bearing: only with it is
#' the Euclidean inner product of two packed rows the Frobenius inner product
#' of the two forms. A Gram built on a naive `upper.tri` packing is off by an
#' \eqn{O(1)}{O(1)} amount and is not a Gram of geometries at all.
#'
#' The sink is excluded. It is unmapped territory reported in budget units
#' under either semantics, so an inner product including it would add a mass to
#' a density; `$receipt$gram$sink_excluded` records the exclusion and the sink
#' budget stays readable on `$receipt$sink_budget`.
#'
#' @section The two estimators:
#' The plug-in Gram books every participant's own within-subject sampling noise
#' as between-participant heterogeneity. Measured over 2 000 Monte Carlo
#' replications of the contract's fixture, it inflates the heterogeneity trace
#' by **+62.7 %** (section 6.2). The cross-fitted Gram takes every inner
#' product across two independent halves of each participant's partitions,
#' \eqn{\widehat\Gamma_{ii'} = \tfrac12(\langle Z_i^{A},Z_{i'}^{B}\rangle +
#' \langle Z_i^{B},Z_{i'}^{A}\rangle)}{Gamma_ii' = (<Z_i^A, Z_i'^B> + <Z_i^B,
#' Z_i'^A>)/2}, which is unbiased for the Gram of the noiseless participant
#' vectors --- on the diagonal too, because the two halves are independent.
#'
#' The cross-fitted Gram is **indefinite by construction** (measured in 100 %
#' of the contract's replications), and a single draw of its trace is not an
#' estimate of the between-participant trace: measured standard deviation
#' `3.58` against a mean of `10.81`. Both facts are reported and never
#' repaired.
#'
#' `estimator = "cross_fit"` re-plans each participant over two disjoint halves
#' of their partitions, so it needs at least four partitions per participant and
#' costs a second streaming read. `estimator = "plug_in"` reads one pass and is
#' admissible for loading *directions*, which section 6.3 measures surviving
#' the inflation; its `$latent` layer is `NULL`, and the refusal that says why
#' travels on `$receipt$latent_refusal` as a record rather than an error.
#'
#' **One departure from section 6.4, stated rather than buried.** Its literal
#' text refuses "any spectrum, variance-explained figure or `n_eff`" from the
#' plug-in Gram; `crossform` reports the plug-in **spectrum** and refuses the
#' other two. The Gram itself is a field of the record and the validator
#' re-derives its eigenvalues from it, so withholding them would be theatre
#' rather than discipline --- while a fraction of an inflated whole, or an
#' effective count over noise modes, is a number that misleads whatever it is
#' labelled. The spectrum is reported beside the estimator's name, the
#' measured inflation and the refusal, never alone.
#'
#' Under `"cross_fit"` the mode *geometry* is reconstructed from the average of
#' the two halves' deviations --- a plug-in reconstruction of a direction,
#' which section 6.3 admits provided the source is named. It is named on
#' `$receipt$gram$mode_reconstruction` and on the printed `modes` line. The
#' eigenvalues beside it still come from the cross-fitted Gram.
#'
#' @section The latent layer:
#' The heterogeneity spectrum is signed. Every functional that treats it as a
#' nonnegative partition --- a fraction, a cumulative curve, an effective mode
#' count --- is defined only on a named projection of it, exactly as
#' [latent_geometry()] requires one level down (section 6.5). `$latent` carries
#' the projected spectrum, its cumulative curve, its `n_eff`, and the mass the
#' projection moved; `$spectrum` stays signed and carries none of them.
#'
#' @section Refusals:
#' Each is an `effect_capability_refusal` (see [catch_refusal()]).
#'
#' * `cross_fitted_subject_gram` --- fewer than four partitions, a declared
#'   split that shares a partition between halves or names one a participant
#'   does not carry, or a metric schedule that cannot be re-planned onto a
#'   half.
#' * `identified_heterogeneity` --- a saturated group design, where no residual
#'   degree of freedom is left for participants to differ in.
#' * `query_space_cross_fit` --- `estimator = "cross_fit"` on an estimated
#'   result, which carries one partition split and cannot supply a second.
#' * `complete_form_normalization` --- a plan declaring `"unit_budget"`, as in
#'   [materialize_population()].
#' * `selected_group_nodes` --- a `nodes` identifier the population does not
#'   carry.
#'
#' One further refusal is *recorded* rather than raised, because the record it
#' belongs to is still worth having: `plug_in_spectrum_functionals` sits on
#' `$receipt$latent_refusal` of every plug-in record and carries the reasons
#' and remedies for the absent `$latent` layer.
#'
#' @param x An `effect_population_plan` from [plan_population()], read through
#'   the streamed complete-form route, or an `effect_population_result` from
#'   [estimate_population()], read in the space of its own query bank.
#' @param ... Passed to the method; every method refuses an unmatched argument
#'   rather than absorbing it.
#' @param estimator `"cross_fit"` (the default on a plan) or `"plug_in"`. See
#'   the section above; the choice is part of the record's identity.
#' @param component Which component of each participant's conservative geometry
#'   to carry: `"total"` (the default), `"coherent"` or `"configuration"`, as
#'   in [materialize_population()].
#' @param nodes Group nodes whose mode geometry is reconstructed, by identifier
#'   or position. `NULL`, the default, reconstructs none: the Gram, the
#'   spectrum and the loadings come out of the first pass, and the mode forms
#'   cost a second one.
#' @param modes Number of leading modes to reconstruct at those nodes. `NULL`
#'   takes `min(3, residual df)`.
#' @param partitions Optional list of two character vectors naming the
#'   partition halves the cross-fit reads across. `NULL` interleaves each
#'   participant's own partitions, which balances the halves against anything
#'   drifting across a session.
#' @param projection The named PSD projection the latent layer is built with,
#'   from [latent_geometry()]'s closed set. It is validated on either method
#'   and enters the record's identity, but it does no work on a result: the
#'   query-space route is plug-in, so that record has no latent layer to
#'   project.
#' @param coordinate_tile Optional positive number of packed coordinates held
#'   in flight at once, as in [materialize_population()]. The tile bounds
#'   memory and is not part of the estimand.
#' @return An `effect_population_heterogeneity`: a sealed record carrying the
#'   `$gram`, its signed `$spectrum` and subject `$loadings`, the `$mode_forms`
#'   reconstructed at the selected `$nodes`, the `$latent` projection, and a
#'   `$receipt` recording the estimator, the partition split, the streaming
#'   bound and every participant's read. `as.data.frame()` returns the subject
#'   loadings in long form.
#'
#'   `$mode_forms` is `NA` at a `(node, coordinate)` cell the Gram could not
#'   see --- a group node reached by no native mass under density semantics ---
#'   because zero there is a form and would read as agreement rather than as
#'   absence. `$receipt$unresolved_columns` counts them.
#' @references `design/population-form-contract.md` (`population-form-v1`),
#'   sections 5 and 6.
#' @family population transports
#' @seealso [estimate_population()] and [materialize_population()] for the
#'   consensus half of the same decomposition, and [latent_geometry()] for the
#'   single-participant form of the projection discipline.
#' @examples
#' # Four participants over four runs, so each cross-fit half holds two.
#' effects <- effect_space(c("face", "house"), basis_id = "heterogeneity:v1")
#' subject <- function(id, n, gain) {
#'   domain <- abstract_domain(n, coordinates = cbind(x = seq_len(n) - 1),
#'     feature_ids = paste0("f", seq_len(n)), id = id)
#'   run <- function(seed) {
#'     set.seed(seed)
#'     matrix(gain * stats::rnorm(2 * n), 2, n,
#'       dimnames = list(c("face", "house"), NULL))
#'   }
#'   rel <- relation(
#'     list(run1 = run(1), run2 = run(2), run3 = run(3), run4 = run(4)),
#'     effects = effects, domain = domain
#'   )
#'   plan_geometry(rel, compile_frame(voxelwise(), domain),
#'     cross_partitions(rel))
#' }
#' carrier <- function(n) anatomical_transport(
#'   native_coords = cbind(seq_len(n) - 1), group_coords = cbind(c(0, 3)),
#'   semantics = "budget", radius = 2
#' )
#' sizes <- c(s01 = 5L, s02 = 6L, s03 = 7L, s04 = 8L)
#' gains <- c(s01 = 1, s02 = 1.4, s03 = 0.7, s04 = 1.1)
#' subjects <- stats::setNames(lapply(names(sizes), function(id)
#'   subject(id, sizes[[id]], gains[[id]])), names(sizes))
#' plan <- plan_population(subjects, lapply(sizes, carrier))
#'
#' # The cross-fitted subject Gram: one N-by-N matrix, streamed. Naming a
#' # group node reconstructs that node's mode geometry in a second pass.
#' split <- heterogeneity(plan, nodes = "group1", modes = 1L)
#' split
#'
#' # Signed and routinely indefinite; the nonnegative functionals live on the
#' # named projection beside it.
#' split$spectrum
#' split$latent$n_eff
#' split$latent$moved_mass
#'
#' # Subject loadings on the leading mode, and that mode's geometry at the
#' # named group node, unpacked from its `svec` coordinates.
#' split$loadings[, "mode1"]
#' crossform:::.unsvec_symmetric(split$mode_forms["group1", "mode1", ], 2L)
#' split$receipt$gram$mode_reconstruction
#'
#' # The plug-in Gram costs one pass instead of two and is admissible for
#' # directions; it carries no `n_eff`, and says why.
#' plug <- heterogeneity(plan, estimator = "plug_in")
#' is.null(plug$latent)
#' plug$receipt$latent_refusal$capability
#' @export
heterogeneity <- function(x, ...) UseMethod("heterogeneity")

#' @rdname heterogeneity
#' @export
heterogeneity.effect_population_plan <- function(x,
                                                 estimator = c("cross_fit",
                                                   "plug_in"),
                                                 component = c("total",
                                                   "coherent",
                                                   "configuration"),
                                                 nodes = NULL, modes = NULL,
                                                 partitions = NULL,
                                                 projection = "psd_projection",
                                                 coordinate_tile = NULL, ...) {
  .check_no_extra_arguments("heterogeneity", ...)
  .validate_population_plan(x, deep = FALSE)
  estimator <- match.arg(estimator)
  component <- match.arg(component)
  projection <- .latent_projection_method(projection)
  # The complete-form route reads a form and not a query, so it inherits
  # `materialize_population()`'s normalization refusal for the same reason: a
  # per-coordinate divisor is a weight varying along the coordinate axis.
  .population_admit_complete_form(x)

  effects <- x$subjects[[1L]]$task$left_relation$effect_space$coordinates
  coordinates <- .population_packed_index(effects)
  width <- nrow(coordinates)
  coordinate_labels <- coordinates$coordinate
  tile <- .population_coordinate_tile(x, width, coordinate_tile)
  starts <- .tile_starts(width, tile)

  full_index <- .population_result_index(x)
  keep <- !full_index$sink
  index <- full_index[keep, , drop = FALSE]
  rownames(index) <- NULL
  node_labels <- as.character(index$node)
  residual_df <- .heterogeneity_residual_df(x)
  selected <- .heterogeneity_nodes(nodes, node_labels)

  split <- list(kind = NA_character_, split = NULL)
  overlap <- character()
  half_ids <- NULL
  plans <- list(x)
  if (identical(estimator, "cross_fit")) {
    split <- .heterogeneity_partition_split(x, partitions)
    overlap <- .heterogeneity_transport_overlap(x, split$split)
    plans <- list(
      .heterogeneity_half_plan(x, split$split, "a"),
      .heterogeneity_half_plan(x, split$split, "b")
    )
    half_ids <- c(a = plans[[1L]]$scientific_plan_id,
      b = plans[[2L]]$scientific_plan_id)
  }

  accumulated <- .heterogeneity_gram(plans, x$model, component, width,
    coordinate_labels, nrow(full_index), keep, starts, tile)
  gram <- accumulated$gram / residual_df
  decomposition <- .heterogeneity_decompose(gram)
  latent <- if (identical(estimator, "cross_fit")) {
    .heterogeneity_latent(unname(decomposition$spectrum), projection)
  }

  mode_forms <- NULL
  if (length(selected)) {
    n_mode <- .heterogeneity_mode_count(modes, residual_df)
    mode_forms <- .heterogeneity_modes(plans, x$model, component, width,
      coordinate_labels, nrow(full_index), keep, starts, tile,
      decomposition$loadings[, seq_len(n_mode), drop = FALSE],
      selected, node_labels, coordinate_labels)
  }

  state <- accumulated$states[[1L]]
  base <- .population_receipt(x,
    .population_readout("complete_form", coordinate_labels), component, 0,
    state$receipts, state$native_totals, state$sink_budget,
    matrix(TRUE, length(x$subjects), width,
      dimnames = list(names(x$subjects), coordinate_labels)),
    do.call(pmax, lapply(accumulated$states, `[[`, "deviations")),
    accumulated$unresolved)
  passes <- length(starts) * length(plans) * (1L + (length(selected) > 0L))
  identity <- .heterogeneity_scientific_id(x$scientific_plan_id, estimator,
    "packed_form", component,
    if (identical(estimator, "cross_fit")) split$split, projection)
  .new_population_heterogeneity(estimator, "packed_form", gram, decomposition,
    latent, mode_forms, selected, index, effects, coordinates, component,
    x$semantics, x$normalization, residual_df,
    .heterogeneity_receipt(base, estimator, "packed_form", gram, split,
      overlap, half_ids, latent,
      list(
        coordinate_tile = as.integer(tile),
        packed_width = as.integer(width),
        passes_per_subject = as.integer(passes),
        # `length(plans)` stacks are live at once: the cross-fitted route
        # holds both partition halves' tiles in order to take their cross
        # product, and a bound quoting one of them understates its own peak
        # by a factor of two.
        group_stack_doubles = as.double(length(plans)) *
          length(x$subjects) * nrow(full_index) * tile,
        gram_doubles = as.double(length(x$subjects))^2,
        refused_dense_doubles = as.double(length(x$subjects)) *
          nrow(index) * width
      ),
      accumulated$unresolved, nrow(index), residual_df),
    identity)
}

#' @rdname heterogeneity
#' @export
heterogeneity.effect_population_result <- function(x,
                                                   estimator = c("plug_in",
                                                     "cross_fit"),
                                                   nodes = NULL, modes = NULL,
                                                   projection =
                                                     "psd_projection", ...) {
  .check_no_extra_arguments("heterogeneity", ...)
  .validate_population_result(x)
  estimator <- match.arg(estimator)
  projection <- .latent_projection_method(projection)
  if (identical(estimator, "cross_fit")) {
    .capability_refusal(paste0(
      "An estimated population result holds one reading of each ",
      "participant, so there is no second, independent partition half to ",
      "take an inner product against. The cross-fitted Gram of ",
      "`population-form-v1` section 6.4 is a property of the *plan*: it ",
      "re-reads every participant over two disjoint halves of their ",
      "partitions, which a finished result can no longer do."
    ),
      capability = "query_space_cross_fit",
      namespace = "population_heterogeneity",
      reasons = c("result_carries_one_partition_split",
        "cross_fit_requires_two_independent_reads"),
      remedies = c(
        paste0(
          "Call `heterogeneity()` on the `plan_population()` estimand ",
          "instead, which streams both halves."
        ),
        paste0(
          "Keep the query-space reading with `estimator = \"plug_in\"`, ",
          "whose loading directions section 6.3 admits."
        )
      ))
  }
  if (!identical(x$basis, "query_bank")) {
    .capability_refusal(paste0(
      "A complete-form population result carries no participant-level ",
      "values: they are the `N x (m + 1) x p` array the streamed route exists ",
      "in order not to build, so there are no deviations here to take a Gram ",
      "of."
    ),
      capability = "heterogeneity_readout",
      namespace = "population_heterogeneity",
      reasons = c("complete_form_result_has_no_subject_values",
        paste0("basis:", x$basis)),
      remedies = paste0(
        "Call `heterogeneity()` on the `plan_population()` estimand, which ",
        "streams the participant deviations without materializing them."
      ))
  }
  keep <- !x$index$sink
  index <- x$index[keep, , drop = FALSE]
  rownames(index) <- NULL
  node_labels <- as.character(index$node)
  selected <- .heterogeneity_nodes(nodes, node_labels)
  labels <- dimnames(x$values)[[3L]]
  residual_df <- as.integer(x$residual_df)
  if (residual_df < 1L) .heterogeneity_saturated_refusal(length(labels))

  # `$residuals` is already the deviation from the group fit at every node and
  # query, so the query-space Gram needs no refit: it is one contraction of an
  # array the result carries. The array is `node x query x subject` with the
  # node axis fastest, so the transpose stacks subjects on rows and
  # node-within-query on columns.
  block <- t(matrix(x$residuals[keep, , , drop = FALSE],
    sum(keep) * dim(x$residuals)[[2L]], length(labels),
    dimnames = list(NULL, labels)))
  finite <- apply(is.finite(block), 2L, all)
  block[, !finite] <- 0
  gram <- tcrossprod(block) / residual_df
  gram <- (gram + t(gram)) / 2
  decomposition <- .heterogeneity_decompose(gram)

  queries <- rownames(x$queries)
  coordinates <- data.frame(coordinate = queries, stringsAsFactors = FALSE)
  mode_forms <- NULL
  if (length(selected)) {
    n_mode <- .heterogeneity_mode_count(modes, residual_df)
    loadings <- decomposition$loadings[, seq_len(n_mode), drop = FALSE]
    weighted <- crossprod(block, loadings)
    mode_forms <- array(NA_real_, c(length(selected), n_mode, length(queries)),
      dimnames = stats::setNames(
        list(node_labels[selected], paste0("mode", seq_len(n_mode)), queries),
        c("node", "mode", "coordinate")
      ))
    pick <- as.integer(outer(selected,
      (seq_along(queries) - 1L) * sum(keep), `+`))
    for (mode in seq_len(n_mode)) {
      mode_forms[, mode, ] <- matrix(weighted[pick, mode], length(selected),
        length(queries))
    }
    mode_forms <- .heterogeneity_scale_modes(mode_forms,
      matrix(finite[pick], length(selected), length(queries)),
      sqrt(colSums(weighted^2)))
  }

  identity <- .heterogeneity_scientific_id(x$scientific_plan_id, estimator,
    "query_bank", x$component, NULL, projection)
  .new_population_heterogeneity(estimator, "query_bank", gram, decomposition,
    NULL, mode_forms, selected, index, NULL, coordinates, x$component,
    x$semantics, x$normalization, residual_df,
    .heterogeneity_receipt(x$receipt, estimator, "query_bank", gram,
      list(kind = NA_character_, split = NULL), character(), NULL, NULL,
      list(
        coordinate_tile = NA_integer_,
        packed_width = as.integer(length(queries)),
        passes_per_subject = 0L,
        group_stack_doubles = as.double(length(labels)) * nrow(x$index) *
          length(queries),
        gram_doubles = as.double(length(labels))^2,
        refused_dense_doubles = as.double(length(labels)) * nrow(index) *
          length(queries)
      ),
      sum(!finite), nrow(index), residual_df),
    identity)
}

# A saturated group design leaves no residual degree of freedom, so there is no
# deviation for participants to differ in: the fit passes through every
# participant exactly and the Gram is identically zero. That is a refusal and
# not a zero, because a zero here would read as "these participants agree".
.heterogeneity_saturated_refusal <- function(n_subject) {
  .capability_refusal(sprintf(paste0(
    "The group design is saturated: %s with a design of the same rank leaves ",
    "no residual degree of freedom, so the fitted value passes through every ",
    "participant and there is no deviation to decompose. crossform refuses ",
    "rather than reporting a Gram of zeros, which would read as agreement ",
    "between participants rather than as an absence of evidence about it."
  ), .msg_count(n_subject, "participant")),
    capability = "identified_heterogeneity",
    namespace = "population_heterogeneity",
    reasons = "group_design_leaves_no_residual_degrees_of_freedom",
    remedies = paste0(
      "Fit a smaller group model, or plan the population with more ",
      "participants than the design has columns."
    ))
}

.heterogeneity_residual_df <- function(plan) {
  residual_df <- length(plan$subjects) - plan$model$rank
  if (residual_df < 1L) .heterogeneity_saturated_refusal(length(plan$subjects))
  residual_df
}

.heterogeneity_mode_count <- function(modes, residual_df) {
  # The Gram's rank is bounded by the residual degrees of freedom: the residual
  # maker annihilates the design's column space on the subject axis, so modes
  # beyond it are numerical zeros and reconstructing one would return noise
  # normalized to unit norm.
  if (is.null(modes)) return(min(3L, residual_df))
  modes <- .check_count(modes, "modes", max = residual_df,
    what = paste0("a number of leading modes no larger than the residual ",
      "degrees of freedom (", residual_df, ")"))
  modes
}

# Presentation -----------------------------------------------------------------

#' @export
as.data.frame.effect_population_heterogeneity <- function(x,
                                                          row.names = NULL,
                                                          optional = FALSE,
                                                          ...) {
  .validate_population_heterogeneity(x)
  loadings <- x$loadings
  data.frame(
    subject = rep(rownames(loadings), times = ncol(loadings)),
    mode = rep(colnames(loadings), each = nrow(loadings)),
    eigenvalue = rep(unname(x$spectrum), each = nrow(loadings)),
    loading = as.numeric(loadings),
    stringsAsFactors = FALSE
  )
}

#' @export
format.effect_population_heterogeneity <- function(x, ...) {
  .pf_inline("effect_population_heterogeneity", x$estimator,
    .msg_count(length(x$subjects), "subject"),
    paste0(nrow(x$index), " group nodes"),
    paste0(x$space, " Gram"))
}

#' @export
print.effect_population_heterogeneity <- function(x, ...) {
  .validate_population_heterogeneity(x)
  leading <- seq_len(min(3L, length(x$spectrum)))
  .pf_emit("effect_population_heterogeneity", list(
    estimator = paste0(x$estimator, " -- ",
      .heterogeneity_estimator_phrase(x$estimator)),
    gram = paste0(length(x$subjects), " x ", length(x$subjects), " over ",
      nrow(x$index), " group nodes, ", nrow(x$coordinates), " ",
      if (identical(x$space, "packed_form")) "packed coordinates" else
        "queries", "; residual df ", x$residual_df),
    ledger = x$ledger,
    subjects = .pf_set(x$subjects, max = 4L),
    spectrum = paste0(.pf_num(unname(x$spectrum)[leading]),
      if (length(x$spectrum) > length(leading)) {
        paste0(" (+", length(x$spectrum) - length(leading), " more)")
      },
      ", ", .heterogeneity_negative_modes(unname(x$spectrum)),
      " negative"),
    n_eff = if (is.null(x$latent)) {
      "refused on the plug-in spectrum"
    } else {
      paste0(.pf_num(x$latent$n_eff), " modes after ", x$latent$method,
        ", moved mass ", .pf_num(x$latent$moved_mass), " (",
        .pf_num(100 * x$latent$moved_share), "%)")
    },
    `cross-fit` = if (identical(x$estimator, "cross_fit")) {
      splits <- x$receipt$cross_fit$partitions
      halves <- splits[[1L]]
      paste0(sub("_partitions$", "", x$receipt$cross_fit$split), ": [",
        .pf_set(halves$a, max = 3L), "] x [", .pf_set(halves$b, max = 3L), "]",
        # Only when it matters: the split is per participant, and the first
        # one stands for all of them exactly when they agree.
        if (length(unique(splits)) > 1L) " (first of several)")
    },
    modes = if (is.null(x$mode_forms)) {
      "not reconstructed (pass `nodes =`)"
    } else {
      # Under `cross_fit` the geometry is a plug-in direction taken from the
      # two halves' average while the eigenvalues beside it are cross-fitted,
      # and section 6.3 admits that only where a reader will see it said.
      paste0(dim(x$mode_forms)[[2L]], " at ", .pf_set(x$nodes, max = 3L),
        if (identical(x$estimator, "cross_fit")) {
          " (direction: half average)"
        })
    },
    unresolved = if (x$receipt$unresolved_columns) {
      paste0(x$receipt$unresolved_columns,
        " node-coordinate cells held out of the Gram")
    },
    estimand = .pf_sig(x$scientific_plan_id)
  ))
  if (identical(x$estimator, "plug_in")) {
    .pf_note(paste0(
      "plug-in: within-subject sampling noise is booked as heterogeneity. ",
      "On the contract's own fixture that inflated the trace by +62.7% ",
      "(population-form-v1 section 6.2); this run's inflation is unknown. ",
      "Directions are admissible, the spectrum is not, and there is no n_eff."
    ))
  } else {
    .pf_note(paste0(
      "cross-fitted: indefinite by construction, and one draw of its trace ",
      "is not an estimate of the between-subject trace ",
      "(population-form-v1 section 6.4). Reported as-is."
    ))
  }
  if (length(x$receipt$cross_fit$transport_partition_overlap)) {
    .pf_note(paste0(
      "a functional transport declares cross-fit partitions that overlap ",
      "these halves, so those reads are not fully held out: ",
      .pf_set(x$receipt$cross_fit$transport_partition_overlap, max = 3L)
    ))
  }
  cat(sprintf("  %-14s%s\n", "next:",
    "as.data.frame(x), x$loadings, x$gram"))
  invisible(x)
}
