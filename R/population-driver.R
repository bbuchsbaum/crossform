# Population execution -------------------------------------------------------
#
# Layer 4 (execution, and the records execution produces). `plan_population()`
# names the group estimand; this file is the only thing that computes a number
# from one. It is the population counterpart of `execution-driver.R`, and it
# calls sideways into that file's runner rather than up into
# `evaluate_geometry()`: the public entry points sit in layer 5, above the
# records an executor emits, and `crossnobis-driver.R` takes the same route for
# the same reason.
#
# `design/population-form-contract.md` (`population-form-v1`) is normative:
#
#   * the evaluation order is query -> transport -> OLS, and it is *forced*
#     rather than chosen: with heterogeneous native frames there is no common
#     node axis, so fit-then-transport is not even definable (section 3,
#     "scope of the tensor argument"). Query and transport act on different
#     tensor factors, so the answer does not depend on which of the two runs
#     first --- but which one *can* run first does (section 3);
#   * the transported total, sink included, is the native total exactly, for
#     signed ledgers. That is a certificate asserted per subject at fit time
#     with a tolerance scaled by the ledger's L1 norm, not by its signed total
#     (section 2);
#   * the budget normalization is a closed set of three that differ in argmax,
#     so it is applied here exactly as section 4.1 writes it and never as a
#     display option;
#   * transported components carry the names of section 8.1
#     (`transported_total`, `native_coherent_ledger`,
#     `native_configuration_ledger`, `sink_budget`), and the bare names
#     `coherent` and `configuration` are forbidden on a transported result.
#
# Why query-first is not an optimization. Complete geometry at a native node is
# `p = q(q+1)/2` packed coordinates; a bank of `K` queries is `K` numbers. The
# naive route materializes `n_i x p` per subject and transports `(m+1) x p`;
# this one holds `n_i x K` and transports `(m+1) x K`, and by section 3 the two
# agree exactly. For a four-condition study read through four contrasts that is
# ten coordinates against four, and the ratio grows like `q^2 / K`.
#
# Measured against the naive route written out in
# `tests/testthat/test-population-executor.R` (`pex_dense_oracle()`), on
# `N = 8` participants, `m = 500` group nodes, `K = 4` queries, native frames
# of 560-900 nodes, as the growth in R's vector high-water mark over the run:
#
#                      query-first   naive
#   q  =  4 (p = 10):    24.0 MiB   24.4 MiB
#   q  =  8 (p = 36):    27.4 MiB   39.1 MiB
#   q  = 12 (p = 78):    30.9 MiB   49.9 MiB
#
# The gap is the packed width, and it is bounded below by the fact that both
# routes run the same kernel: at four conditions the saving is 1.5 % and at
# twelve it is 38 %. The retained arrays are the cleaner statement of the same
# thing --- the transported group stack is `N(m+1)K` doubles here against
# `N(m+1)p` there (125 KiB vs 2.4 MiB at `q = 12`), and the largest native
# block held at once is `n_i K` against `2 n_i p`, the factor two being the
# separate `total` and `coherent` stores a materialized geometry keeps.

# The three admitted components and the transported ledger name each one takes.
# `population-form-v1` section 8.1 is the whole reason this table exists: the
# transported total *is* the group node's own total under a feature-additive
# metric, while the transported coherent part is a ledger of native-node
# coherence carried to a group location and is a different object from any
# group-node common mode. Naming them apart is what stops a reader inferring
# the second from the first.
.population_ledger_names <- c(
  total = "transported_total",
  coherent = "native_coherent_ledger",
  configuration = "native_configuration_ledger"
)

# Section 11: the budget-preservation tolerance is `1e-12 * sum_x |c_x|`. The
# L1 norm and not the total, because the total is a signed sum that may sit
# arbitrarily close to zero while its summands are large, and a
# relative-to-total tolerance would be unbounded exactly where the estimate is
# weakest.
.population_budget_tolerance <- 1e-12

# The query bank ---------------------------------------------------------------
#
# One bank parser for the whole package. `.sampling_query_bank()` is D8's, and
# reusing it is what makes `estimate_population(plan, queries)` and
# `sampling_covariance(x, queries = )` agree on what a bank *is* --- a numeric
# contrast vector, a K-by-q matrix of contrast rows, or a list of them, aligned
# against the relation's declared effect order and named once. A second parser
# here would be a second dialect.
.population_query_bank <- function(plan, queries) {
  if (missing(queries)) {
    .input_error(paste0(
      "`queries` is required: pass the bank of contrasts the group estimate ",
      "is about, for example `queries = rbind(faces = c(1, -1))`. A ",
      "population fit is a fit *of* a query --- there is no default readout ",
      "of a transported geometry, and guessing one would name an estimand ",
      "the plan never sealed."
    ),
      arg = "queries", received = "no argument",
      expected = "a contrast vector, a K-by-q contrast matrix, or a list")
  }
  effects <- plan$subjects[[1L]]$task$left_relation$effect_space$coordinates
  bank <- .sampling_query_bank(queries, effects)
  width <- length(effects) * (length(effects) + 1L) / 2L
  packed <- matrix(vapply(seq_len(nrow(bank$matrix)), function(k) {
    .svec_symmetric(tcrossprod(bank$matrix[k, ]))
  }, numeric(width)), nrow = width,
    dimnames = list(NULL, bank$labels))
  # The packed operator is exactly what `contrast_energy()` lowers a contrast
  # to (`bilinear_query(tcrossprod(weights))`), so a one-query bank and
  # `contrast_energy()` read the same physical query and not merely the same
  # arithmetic.
  list(matrix = bank$matrix, labels = bank$labels, packed = packed,
    effects = effects)
}

# Normalization ----------------------------------------------------------------
#
# Section 4.1's group functionals, written out rather than templated: `none` is
# the mean subject ledger and `unit_budget` is the *unweighted* mean of
# *unit-normalized* ledgers, and no single weighted-mean template covers both
# (the contract's own section 4.1 records an earlier draft getting this wrong).
# Both are per-subject scalar rescalings of the response, so both commute with
# transport (section 4.4, measured at 1.1e-16) --- which is why they are
# applied here, to the transported ledger, and give the same answer as
# rescaling the native one.
#
# The scalar is per subject *and per query*: `T_i` is the native total of one
# queried ledger, and a bank of K queries has K of them.
.population_scale <- function(plan, native, budget_floor) {
  totals <- colSums(native)
  if (identical(plan$normalization, "none")) {
    return(list(totals = totals, scale = rep(1, length(totals)),
      admitted = rep(TRUE, length(totals))))
  }
  # Section 4.3: `unit_budget`'s divisor is a signed estimate, so it can sit
  # near zero or be negative, and a subject whose total is not bounded away
  # from zero is marked NA rather than emitting a divergent share. The
  # criterion is scale-free --- `|T_i| >= floor * sum_x |c_ix|` --- because a
  # signed total may be a near-cancellation of large summands, and that is
  # exactly the case an absolute floor cannot see.
  norm <- colSums(abs(native))
  admitted <- abs(totals) > budget_floor * norm & totals != 0
  scale <- ifelse(admitted, 1 / totals, NA_real_)
  list(totals = totals, scale = scale, admitted = admitted)
}

# Section 4.1 admits `unit_budget` because it conserves: `sum_j g_j = 1`
# exactly. That is a statement about a *budget*, and density conserves nothing
# at all (section 1.3, measured: the density columns of the contract's own
# fixture satisfy no conservation law). Dividing a density profile by a native
# budget total therefore produces a number whose documented property is false,
# which is the one thing a receipt may not do.
.population_admit_normalization <- function(plan) {
  if (!identical(plan$normalization, "unit_budget") ||
      identical(plan$semantics, "budget")) {
    return(invisible(NULL))
  }
  .capability_refusal(paste0(
    "`normalization = \"unit_budget\"` divides each participant's ledger by ",
    "its own native budget total, and it is admitted because the result ",
    "conserves: every participant contributes budget one, so the group ",
    "ledger sums to one exactly. Density semantics transports budget per unit ",
    "transported territory, which satisfies no conservation law, so the ",
    "quotient conserves nothing and the share it looks like is a share of ",
    "nothing. crossform will not print a number under a name whose stated ",
    "property it does not have."
  ),
    capability = "budget_normalization_semantics",
    namespace = "population_execution",
    reasons = c(
      "unit_budget_requires_budget_semantics",
      paste0("declared_semantics:", plan$semantics)
    ),
    remedies = c(
      paste0(
        "Plan the population with `normalization = \"none\"`, which is the ",
        "mean subject density profile and makes no conservation claim."
      ),
      paste0(
        "Declare `semantics = \"budget\"` on every transport if the ",
        "attribution *share* is the estimand you want."
      )
    ))
}

# Budget preservation ----------------------------------------------------------
#
# Section 2's certificate, asserted per subject at fit time on the raw
# transported ledger --- before any normalization, because a per-subject scalar
# rescales both sides and would make the check vacuous.
#
# The argument is `P 1 = 1` and never `c >= 0`, so this holds for the signed
# crossvalidated ledgers a conservative geometry actually produces, and it
# holds for every component: the coherent and configuration ledgers are
# ledgers too. It is checked rather than assumed because deleting the sink
# column from a well-formed `P` loses 48.1 % of the native total on the
# contract's own fixture and leaves no trace anywhere else.
.population_budget_certificate <- function(label, native, carried) {
  native_total <- colSums(native)
  carried_total <- colSums(carried)
  norm <- pmax(colSums(abs(native)), .Machine$double.eps)
  deviation <- abs(carried_total - native_total) / norm
  worst <- which.max(deviation)
  if (deviation[[worst]] > .population_budget_tolerance) {
    .contract_error(sprintf(paste0(
      "Budget preservation failed for participant `%s` on query `%s`: the ",
      "transported total including the sink is %s against a native total of ",
      "%s, a relative deviation of %s. A row-stochastic transport preserves ",
      "the signed total exactly, so a failure here means the operator is ",
      "ill-formed rather than imprecise, and crossform refuses rather than ",
      "renormalizing it."
    ), label, colnames(native)[[worst]], format(carried_total[[worst]]),
      format(native_total[[worst]]), format(deviation[[worst]])))
  }
  max(deviation)
}

# The group index --------------------------------------------------------------
#
# The plan's `m` group nodes plus the sink, which is a row of the result and
# not a footnote: section 3.3 fits it and section 1.1 requires it materialized
# even when it is identically zero, because unmapped territory has to be a
# number a reader can see. `units` is what stops a density result reading as if
# its sink were a density: the sink is reported in budget units under either
# semantics, since a density of unmapped territory has no referent
# (section 1.3).
.population_result_index <- function(plan) {
  group <- plan$group_index
  # A group index is an open table --- a caller may carry any column on it ---
  # so the two the result adds have to be checked rather than assumed free.
  # Overwriting a caller's own `sink` or `units` column would replace declared
  # metadata with the result's bookkeeping and leave nothing saying so;
  # `location_transport()` refuses a `coord1` clash for the same reason.
  clash <- intersect(c("sink", "units"), names(group))
  if (length(clash)) {
    .capability_refusal(sprintf(paste0(
      "The group index already carries %s, which is where a population ",
      "result records whether a row is the sink and which units its values ",
      "are in. crossform will not overwrite a declared column with its own ",
      "bookkeeping."
    ), paste0("`", clash, "`", collapse = " and ")),
      capability = "reserved_group_index_columns",
      namespace = "population_execution",
      reasons = paste0("group_index_column_reserved:", clash),
      remedies = paste0(
        "Rename ", paste0("`", clash, "`", collapse = " and "),
        " in the `group_index` every transport was built against, then ",
        "re-plan the population."
      ))
  }
  sink <- group[1L, , drop = FALSE]
  for (column in names(sink)) {
    sink[[column]] <- group[[column]][NA_integer_]
  }
  sink$node <- .transport_sink_label
  index <- rbind(group, sink)
  index$sink <- c(rep(FALSE, nrow(group)), TRUE)
  index$units <- c(rep(plan$semantics, nrow(group)), "budget")
  rownames(index) <- NULL
  index
}

# Uncertainty ------------------------------------------------------------------
#
# What the transported covariance would need, and why it is refused.
#
# The transported value at group node `j` is `b_j = sum_x P_xj z_x`, so
#
#   Cov(b_j) = sum_x sum_x' P_xj P_x'j Cov(z_x, z_x')
#
# and every term with `x != x'` is a *cross-node* sampling covariance. D8
# refuses exactly that object under capability `cross_node_sampling_covariance`
# (reasons `spatial_covariance_model_unavailable` and
# `conservation_gives_no_uncertainty`), and `population-form-v1` section 4.5
# records that the refusal stands. A transport with more than one native row
# per group column --- which is every transport that does anything --- needs
# those terms, so the transported covariance cannot be assembled from what
# exists. What *can* be carried honestly is the per-native-node marginal block
# D8 does supply, untransported and labelled as such.
#
# The one case where the cross terms vanish is a group column fed by exactly
# one native row, where `Cov(b_j) = P_xj^2 Cov(z_x)`. E8 implements exactly
# that case and nothing wider (`.population_within_uncertainty()` below): a
# real transport is a mixture of the two, so the result marks which group
# columns carry a covariance and which carry an absence, per participant,
# rather than presenting one number for both. The general refusal is unchanged
# and still travels with the result.
.population_transport_covariance_refusal <- function() {
  list(
    capability = "transported_sampling_covariance",
    namespace = "population_execution",
    reasons = c(
      "cross_node_sampling_covariance_unavailable",
      "transported_value_mixes_native_nodes"
    ),
    remedies = c(
      paste0(
        "Read `$uncertainty$native` as what it is: one K-by-K marginal ",
        "covariance per native node per participant, before transport."
      ),
      paste0(
        "The refusal lifts for a group column fed by exactly one native row, ",
        "where the cross-node terms carry weight zero; a general transport ",
        "has no such column."
      )
    ),
    message = paste0(
      "The sampling covariance of a transported value is refused. A group ",
      "node's value is a weighted sum over many native nodes, so its ",
      "covariance needs the covariance *between* native nodes, and crossform ",
      "has no route that produces one: the analytic sampling law is local to ",
      "one measurement, and a conservative frame supplies no substitute, ",
      "because conservation is a law about estimates and not about their ",
      "uncertainty. The per-native-node blocks are carried untransported."
    )
  )
}

# The between-subject layer's ingredients, recorded once at fit time.
#
# `population_uncertainty()` needs two things the result does not otherwise
# carry: the unscaled covariance `(X'X)^-1` of the group design, and the names
# of the columns it is indexed by. Everything else the between-subject layer
# needs --- the residuals, the residual df, the coefficients --- is already on
# the result, so the layer is computed when it is read rather than stored as
# three more `node x query x term` arrays beside the coefficients.
#
# The design does not vary along the node or the query axis (that is why one
# factorization serves the whole fit), so one `p`-by-`p` matrix serves every
# group node and every query. It is stored unpivoted and named, because
# `qr.coef()` returns coefficients in the design's original column order and a
# standard error indexed differently from its estimate is a bug waiting to be
# printed.
.population_between_uncertainty <- function(model) {
  rank <- as.integer(model$rank)
  columns <- model$columns
  unscaled <- matrix(NA_real_, length(columns), length(columns),
    dimnames = list(columns, columns))
  # `plan_population()` refuses a rank-deficient design, so the pivot is a
  # permutation of every column and the inverse exists; the guard is here so
  # that a future admission of aliasing degrades to `NA` rather than to a
  # silently wrong block.
  if (rank == length(columns)) {
    pivot <- model$pivot[seq_len(rank)]
    triangle <- qr.R(model$qr)[seq_len(rank), seq_len(rank), drop = FALSE]
    unscaled[pivot, pivot] <- chol2inv(triangle)
  }
  list(
    scope = "between_subject",
    basis = "group_ols_residual",
    estimator = "subject_residual_covariance_about_the_group_fit",
    terms = columns,
    unscaled_covariance = unscaled,
    # Not a hedge and not conditional on anything: the ratio of a transported
    # ledger to a between-subject standard error is a t *statistic*, and
    # calling its tail probability a p-value would import a normal-errors,
    # correctly-specified-model, homogeneous-transport assumption that no part
    # of this package has checked on real data.
    calibration = "uncalibrated",
    calibration_evidence = paste0(
      "benchmarks/run-population-null-coverage.R measures the coverage of ",
      "the nominal 95% interval under a correctly specified group model; ",
      "see benchmarks/POPULATION-NULL-COVERAGE.md."
    )
  )
}

# The within-subject layer, admitted only where it is exact.
#
# For participant `i` the transported value at group column `u` is a fixed
# linear functional `w' z_i` of that participant's native query values, with
# `w` the column of the transport operator (divided by the transported row
# mass under `"density"`, which is data-free). Its variance is
# `sum_x sum_x' w_x w_x' Cov(z_ix, z_ix')`, and D8 supplies only the diagonal
# blocks `Cov(z_ix)`.
#
# Where the column has exactly one nonzero coefficient the cross terms carry
# weight zero and `Var(w'z) = w_x^2 Var(z_ix)` is **exact**: no independence
# is assumed, and none is claimed. Where it has two or more, dropping the
# cross terms would assume the native nodes are uncorrelated --- false for
# overlapping searchlight supports under spatially correlated noise, and
# biased in a known direction (positive spatial correlation makes the diagonal
# sum an *under*-estimate), which is a worse failure than an absence. So the
# layer is admitted per participant and per group column, and refused
# elsewhere with the reason recorded rather than defaulted.
#
# Two whole-layer gates come before the per-column one:
#
#   `same_data_ratio_normalization` --- under `"unit_budget"` the response is
#     the ledger divided by a total estimated from the same data, so its
#     variance needs the variance of that divisor and its covariance with the
#     numerator. Section 4.3 records that the per-participant standard error
#     does not exist; a delta-method ratio variance would be inventing one.
#
#   `native_node_labels_unaligned` --- the covariance batch is named by the
#     nodes `rdm_sampling_covariance(at = )` read, and the transport is named
#     by its own `native_index`. Binding two differently-named node sets by
#     position would attach node `x`'s error bar to node `y`'s coefficient,
#     which is the same error `aligned_subject_uncertainty` refuses one level
#     up. This one is *recorded*, not raised: the untransported blocks are
#     still honest and still shipped, and only the carve-out is unavailable.
.population_within_uncertainty <- function(plan, bank, native, index) {
  nodes <- as.character(index$node)
  labels <- names(plan$subjects)
  record <- function(status, reasons, remedies) {
    list(
      scope = "transported_single_source_column",
      basis = "query_bank",
      status = status,
      exactness = "exact_where_admitted",
      assumption = "none: the cross-node terms carry weight zero",
      admitted = NULL,
      refusal = list(
        capability = "transported_sampling_covariance",
        namespace = "population_execution",
        reasons = reasons,
        remedies = remedies
      )
    )
  }
  if (!identical(plan$normalization, "none")) {
    return(record("refused", "same_data_ratio_normalization", paste0(
      "Estimate the population under `normalization = \"none\"`. A ",
      "`unit_budget` share is a ratio of two quantities read from the same ",
      "data, and section 4.3 records that the standard error of its divisor ",
      "does not exist yet."
    )))
  }
  queries <- bank$labels
  variance <- array(NA_real_, c(length(nodes), length(queries), length(labels)),
    dimnames = stats::setNames(list(nodes, queries, labels),
      c("node", "query", "subject")))
  admitted <- matrix(FALSE, length(nodes), length(labels),
    dimnames = list(nodes, labels))
  coefficient <- matrix(NA_real_, length(nodes), length(labels),
    dimnames = list(nodes, labels))
  source_node <- matrix(NA_character_, length(nodes), length(labels),
    dimnames = list(nodes, labels))
  unaligned <- character()
  for (label in labels) {
    carrier <- plan$transport[[label]]
    blocks <- native[[label]]
    native_nodes <- as.character(carrier$native_index$node)
    if (!setequal(names(blocks), native_nodes)) {
      unaligned <- c(unaligned, label)
      next
    }
    # Column support is read off the sparse structure rather than by
    # densifying: the operator is `n_i` by `m + 1` and `n_i` is a whole native
    # frame, so `as.matrix()` here would allocate the one array the population
    # layer exists not to build. `drop0()` first, because a stored zero is not
    # a source and would make a single-source column look like two.
    operator <- Matrix::drop0(methods::as(carrier$matrix, "CsparseMatrix"))
    # Under `"density"` the group rows are divided by transported mass, which
    # is data-free, so the carried value is still a fixed linear functional of
    # the native values and the carve-out is unchanged. The sink is left alone
    # (section 1.3): a density of unmapped territory has no referent.
    divisor <- rep(1, ncol(operator))
    if (identical(carrier$semantics, "density")) {
      mass <- as.numeric(Matrix::crossprod(carrier$matrix, carrier$row_mass))
      group <- seq_len(ncol(operator) - 1L)
      divisor[group] <- ifelse(mass[group] > 0, mass[group], NA_real_)
    }
    starts <- operator@p
    for (column in seq_along(nodes)) {
      if (starts[[column + 1L]] - starts[[column]] != 1L) next
      position <- starts[[column]] + 1L
      weight <- operator@x[[position]] / divisor[[column]]
      if (!is.finite(weight) || weight == 0) next
      source <- native_nodes[[operator@i[[position]] + 1L]]
      diagonal <- sampling_covariance(blocks[[source]], "diagonal")
      admitted[column, label] <- TRUE
      coefficient[column, label] <- weight
      source_node[column, label] <- source
      variance[column, , label] <- weight^2 * as.numeric(diagonal[queries])
    }
  }
  if (length(unaligned)) {
    return(record("refused", paste0("native_node_labels_unaligned:",
      paste(unaligned, collapse = ",")), paste0(
      "Name the transport's `native_index` with the node identifiers ",
      "`rdm_sampling_covariance(at = )` reported --- for a voxelwise frame ",
      "those are the domain's feature identifiers --- so the two node sets ",
      "bind by name rather than by position."
    )))
  }
  value <- record(
    if (any(admitted)) "partial" else "refused",
    "transported_value_mixes_native_nodes",
    paste0(
      "A group column fed by two or more native rows needs the covariance ",
      "*between* those rows, which D8 refuses (capability ",
      "`cross_node_sampling_covariance`). Use a transport that assigns at ",
      "most one native node to each group node if the transported variance ",
      "is what you need."
    )
  )
  value$admitted <- admitted
  value$coefficient <- coefficient
  value$source_node <- source_node
  value$variance <- variance
  value$admitted_columns <- as.integer(sum(admitted))
  value$columns <- as.integer(length(admitted))
  value
}

.population_uncertainty <- function(plan, bank, index, uncertainty) {
  layers <- list(
    between = .population_between_uncertainty(plan$model),
    within = NULL,
    # Section 7 keeps the two layers apart, and so does this record: a
    # within-subject sampling variance and a between-subject residual variance
    # answer different questions, and a sum of the two would answer neither
    # without a variance-components model nobody here has fitted.
    separation = paste0(
      "between_subject and within_subject are reported separately and are ",
      "never pooled"
    )
  )
  if (is.null(uncertainty)) return(layers)
  labels <- names(plan$subjects)
  if (!is.list(uncertainty) || inherits(uncertainty, "effect_sampling_covariance") ||
      inherits(uncertainty, "effect_sampling_covariance_batch")) {
    .input_error(paste0(
      "`uncertainty` must be a named list of `rdm_sampling_covariance()` ",
      "results, one per participant. A single covariance belongs to one ",
      "participant's own frame, so pairing it with the population by ",
      "position would attach one person's error bars to another's estimate."
    ), arg = "uncertainty", received = .msg_value(uncertainty),
      expected = "a named list of sampling covariances, one per participant")
  }
  given <- names(uncertainty)
  if (!.is_strings(given, unique = TRUE) || !setequal(given, labels)) {
    .capability_refusal(paste0(
      "Every participant needs exactly one sampling covariance and every ",
      "covariance needs a participant. The names are the binding, exactly as ",
      "they are for the transports."
    ),
      capability = "aligned_subject_uncertainty",
      namespace = "population_execution",
      reasons = c(
        paste0("subject_without_uncertainty:",
          setdiff(labels, as.character(given))),
        paste0("uncertainty_without_subject:",
          setdiff(as.character(given), labels))
      ),
      remedies = paste0(
        "Name `uncertainty` with the participant identifiers of the plan: ",
        .msg_names(labels), "."
      ))
  }
  # D8's own route, per subject. It refuses an uncentred contrast --- the
  # crossvalidated distance basis carries no information about one --- and that
  # refusal propagates rather than being caught, because a bank the covariance
  # cannot express is a bank the caller has to know about.
  native <- lapply(uncertainty[labels], function(value) {
    sampling_covariance(value, queries = bank$matrix)
  })
  names(native) <- labels
  layers$within <- .population_within_uncertainty(plan, bank, native, index)
  c(list(
    native = native,
    basis = "query_bank",
    scope = "native_node_marginal",
    transported = .population_transport_covariance_refusal()
  ), layers)
}

# The readout axis -------------------------------------------------------------
#
# Both population executors fit the same model on the same group nodes; what
# separates them is the axis they read the transported geometry along. A
# `"query_bank"` readout is `K` declared contrasts; a `"complete_form"` readout
# is the `p = q(q+1)/2` packed `svec` coordinates of the whole form. One record
# describes either, so the receipt, the identity and the validator branch on a
# field rather than on which verb was called.
.population_result_bases <- c("query_bank", "complete_form")

.population_readout <- function(basis, labels, queries = NULL) {
  list(basis = basis, labels = labels, queries = queries)
}

# Identity ---------------------------------------------------------------------
#
# The result's own scientific identity: its parent estimand, the axis it was
# read along, the component it names, and the normalization criterion that
# decided which participants contributed a share. `budget_floor` is in it
# because it decides which subjects are marked NA, and a different set of
# contributing participants is a different population.
#
# `basis` is in the digest because the two readouts name two estimands, not two
# renderings of one: a query bank fixes `K` declared contrasts and a complete
# form fixes none, and the second is not recoverable from the first.
.population_result_id <- function(plan, readout, component, budget_floor) {
  .sha256_signature(list(
    schema_version = 1L,
    contract_version = "population-form-v1",
    role = "population_result",
    basis = readout$basis,
    parent = plan$scientific_plan_id,
    component = component,
    ledger = unname(.population_ledger_names[[component]]),
    queries = unname(readout$queries),
    query_labels = readout$labels,
    normalization = plan$normalization,
    budget_floor = if (identical(plan$normalization, "unit_budget")) {
      as.numeric(budget_floor)
    } else {
      NULL
    }
  ), "population-sha256:")
}

# The verb ---------------------------------------------------------------------

#' Estimate a planned population form through a bank of queries
#'
#' `estimate_population()` is the authorized group-level execution verb: it
#' takes the estimand [plan_population()] sealed, reads each participant's
#' geometry through one bank of `K` contrasts, carries the queried ledgers onto
#' the shared group nodes with that participant's declared transport, and fits
#' the group model once at every group node and query.
#'
#' @section The order, and why it is not a choice:
#' The evaluation order is **query, then transport, then fit**. Query and
#' transport act on different axes of the same array --- the query mixes
#' experimental coordinates, the transport mixes spatial nodes --- so they
#' commute exactly, and under a subject-constant weight operator the fit
#' commutes with both. What is *not* free is which may run first: participants
#' have different native frames, so there is no common node axis to fit at
#' before transporting, and fit-then-transport is not definable at all. The
#' order is recorded on the plan rather than selected here.
#'
#' Reading the query first is what makes the route cheap. Complete geometry at
#' a native node is \eqn{q(q+1)/2}{q(q+1)/2} packed coordinates; a bank of
#' `K` queries is `K` numbers, and complete geometry is never allocated.
#'
#' @section What the result is a ledger of:
#' `component` selects which part of each participant's conservative geometry
#' is carried, and the transported object takes its own name
#' (`population-form-v1` section 8.1) because it is not the group-node quantity
#' the bare name would suggest.
#'
#' | `component` | `$ledger` | what it is |
#' |---|---|---|
#' | `"total"` | `transported_total` | a genuine group-node total for that participant |
#' | `"coherent"` | `native_coherent_ledger` | native-node coherent evidence carried to a group node |
#' | `"configuration"` | `native_configuration_ledger` | native-node configuration evidence carried to a group node |
#'
#' `native_coherent_ledger + native_configuration_ledger = transported_total`
#' holds exactly; what fails is reading the coherent ledger as a group-node
#' common mode. Transport carries *nodes*, so there is no group frame and no
#' group-node weight vector for such a mode to be defined against.
#'
#' The sink is a row of the result, labelled `<sink>`, not a discarded
#' remainder: it is fitted like any other column because covariate-linked sink
#' mass is how a differentially failing transport becomes visible. It is never
#' a value at a location, and it is always reported in budget units.
#'
#' @section Budget preservation:
#' Under `"budget"` semantics each participant's transported total, sink
#' included, equals its native total exactly. The certificate is asserted per
#' participant and per query at fit time, with a tolerance of `1e-12` times the
#' ledger's \eqn{L^1}{L1} norm --- not its total, which is a signed sum that
#' may sit near zero while its summands are large. A participant that fails it
#' has an ill-formed transport, and the run refuses rather than renormalizing.
#'
#' @section Normalization:
#' `"none"` fits the ledgers as they are, so participants with more evidence
#' weigh more. `"unit_budget"` divides each participant's ledger by its own
#' native total first, so each contributes budget one; the divisor is a
#' *signed* estimate, so a participant whose total is not bounded away from
#' zero by `budget_floor` is marked `NA` for that query rather than emitting a
#' divergent share. The default `budget_floor = 0` admits every nonzero total,
#' which is the weakest honest criterion and is recorded as such in the
#' receipt: the constant is an open maintainer decision, and it needs a
#' per-participant standard error that does not exist yet. Both normalizations
#' are computed against a total read from the same data as the ledger, which
#' the receipt declares (`budget_estimate = "same_data_ratio"`).
#'
#' @section Uncertainty:
#' `$uncertainty` always carries `$between`, the ingredients of the
#' between-subject layer: the group design's unscaled covariance
#' \eqn{(X'X)^{-1}}{(X'X)^-1} and the terms it is indexed by. That is what
#' [population_uncertainty()] needs and the result does not otherwise hold;
#' the standard errors themselves are computed when they are read rather than
#' stored as three more `node`-by-`query`-by-`term` arrays.
#'
#' Supplying `uncertainty` adds two more blocks. `$native` is D8's own
#' per-native-node covariance per participant, carried **untransported**, with
#' `$transported` recording the refusal: a group node's value mixes native
#' nodes, so its covariance needs the covariance *between* them, and no route
#' produces one. `$within` is the part of that refusal E8 lifts --- a group
#' column fed by exactly one native row, where the cross-node terms carry
#' weight zero --- admitted per participant and per column and absent
#' elsewhere. No independence assumption is made anywhere.
#'
#' The two layers are reported separately and are never pooled; there is no
#' field holding their sum. See [population_uncertainty()].
#'
#' @section Refusals:
#' Each is an `effect_capability_refusal` (see [catch_refusal()]).
#'
#' * `budget_normalization_semantics` --- `"unit_budget"` on a `"density"`
#'   plan. Density conserves nothing, so the share would be a share of
#'   nothing.
#' * `aligned_subject_uncertainty` --- `uncertainty` names that do not match
#'   the participants.
#' * `distance_basis_query_bank` --- an uncentred contrast in the bank when
#'   `uncertainty` is supplied. Point estimates admit uncentred contrasts;
#'   their sampling law does not.
#'
#' @param plan An `effect_population_plan` from [plan_population()].
#' @param queries The bank of contrasts to read: a numeric contrast vector, a
#'   `K`-by-`q` matrix with one contrast per row, or a list of contrast
#'   vectors, over the participants' shared experimental effects. Named rows
#'   name the queries; unnamed ones are numbered. This is the same bank
#'   [sampling_covariance()] takes.
#' @param component Which component of each participant's conservative
#'   geometry to carry: `"total"` (the default), `"coherent"`, or
#'   `"configuration"`. See the table above for the name the transported
#'   object takes.
#' @param uncertainty Optional named list of [rdm_sampling_covariance()]
#'   results, one per participant, lowered onto the query bank and carried
#'   **untransported**. The covariance of a transported value is refused
#'   except where a group column is fed by exactly one native row, and both
#'   the refusal and the admitted carve-out travel with the result. The
#'   between-subject layer needs none of this and is always available.
#' @param budget_floor Nonnegative relative floor for the `"unit_budget"`
#'   divisor: a participant contributes a share for query `k` only when
#'   \eqn{|T_i| > \mathrm{floor} \cdot \sum_x |c_{ix}|}{|T_i| > floor * sum_x
#'   |c_ix|}. Ignored under `"none"`.
#' @return An `effect_population_result`: a sealed record carrying
#'   `$coefficients` (a `node`-by-`query`-by-`term` array), `$values`,
#'   `$fitted` and `$residuals` (`node`-by-`query`-by-`subject`), the
#'   `$index` of group nodes plus the sink, `$queries`, `$ledger`,
#'   `$uncertainty` (see the section above; read it through
#'   [population_uncertainty()]), and a `$receipt` recording every
#'   participant's read, its
#'   transport signature, the budget certificate, and the normalization.
#'   `as.data.frame()` returns the coefficient table in long form.
#' @references `design/population-form-contract.md` (`population-form-v1`),
#'   sections 2, 3, 4 and 8.
#' @family population transports
#' @seealso [plan_population()] for the estimand this executes,
#'   [transport_values()] for the carrying step, [population_uncertainty()]
#'   for the error bars on what it produced, and [contrast_energy()] for
#'   the single-participant reading of the same query.
#' @examples
#' # Four participants on different native frames, two shared group nodes.
#' effects <- effect_space(c("face", "house"), basis_id = "population:v1")
#' subject <- function(id, n, gain) {
#'   domain <- abstract_domain(n, coordinates = cbind(x = seq_len(n) - 1),
#'     feature_ids = paste0("f", seq_len(n)), id = id)
#'   values <- function(divisor) matrix(
#'     gain * seq_len(2 * n) / (n * divisor), 2, n,
#'     dimnames = list(c("face", "house"), NULL)
#'   )
#'   rel <- relation(list(run1 = values(1), run2 = values(2)),
#'     effects = effects, domain = domain)
#'   plan_geometry(rel, compile_frame(voxelwise(), domain),
#'     cross_partitions(rel))
#' }
#' carrier <- function(n) anatomical_transport(
#'   native_coords = cbind(seq_len(n) - 1), group_coords = cbind(c(0, 3)),
#'   semantics = "budget"
#' )
#' sizes <- c(s01 = 5L, s02 = 6L, s03 = 7L, s04 = 8L)
#' gains <- c(s01 = 1, s02 = 1.4, s03 = 0.7, s04 = 1.1)
#' subjects <- stats::setNames(lapply(names(sizes), function(id)
#'   subject(id, sizes[[id]], gains[[id]])), names(sizes))
#' plan <- plan_population(subjects, lapply(sizes, carrier))
#'
#' # One bank of two contrasts, read at every group node in one pass.
#' fit <- estimate_population(plan, rbind(
#'   `face-house` = c(1, -1), `face+house` = c(1, 1)
#' ))
#' fit
#'
#' # The group mean at each node and query, with the sink as its own row.
#' fit$coefficients[, , "(Intercept)"]
#'
#' # Budget preservation: each participant's transported total, sink
#' # included, is its native total.
#' fit$receipt$budget$max_relative_deviation < 1e-12
#'
#' # The transported object is named for what it is.
#' fit$ledger
#' @export
estimate_population <- function(plan, queries,
                                component = c("total", "coherent",
                                  "configuration"),
                                uncertainty = NULL, budget_floor = 0) {
  .validate_population_plan(plan, deep = FALSE)
  component <- match.arg(component)
  budget_floor <- .check_number(budget_floor, "budget_floor",
    nonnegative = TRUE,
    what = "one nonnegative relative floor on the unit-budget divisor")
  .population_admit_normalization(plan)
  bank <- .population_query_bank(plan, queries)

  labels <- names(plan$subjects)
  n_subject <- length(labels)
  n_query <- length(bank$labels)
  index <- .population_result_index(plan)
  nodes <- nrow(index)
  node_labels <- as.character(index$node)

  # The stack is `N` by `(m + 1) * K`, and it is the only object here whose
  # size grows with the group node count. Nothing native survives the loop:
  # each participant's `n_i` by `K` view is contracted to `(m + 1)` by `K` and
  # dropped before the next participant is opened, so peak memory is set by
  # the group side rather than by the largest native frame.
  stack <- matrix(NA_real_, n_subject, nodes * n_query,
    dimnames = list(labels, NULL))
  receipts <- vector("list", n_subject)
  names(receipts) <- labels
  native_totals <- matrix(NA_real_, n_subject, n_query,
    dimnames = list(labels, bank$labels))
  sink_budget <- native_totals
  admitted <- matrix(TRUE, n_subject, n_query,
    dimnames = list(labels, bank$labels))
  deviations <- numeric(n_subject)

  for (position in seq_len(n_subject)) {
    label <- labels[[position]]
    carrier <- plan$transport[[label]]
    # Sideways into the executor's own runner, exactly as `crossnobis()` does:
    # `evaluate_geometry()` is the layer-5 entry point for the same call, and
    # a driver must not reach up into the views it feeds.
    view <- .run_geometry_compiler(
      plan$subjects[[label]], query = bank$packed, component = component
    )
    native <- view$values
    colnames(native) <- bank$labels
    receipts[[label]] <- view$receipt

    carried <- transport_values(carrier, native)
    if (identical(plan$semantics, "budget")) {
      deviations[[position]] <-
        .population_budget_certificate(label, native, carried)
    }
    native_totals[position, ] <- colSums(native)
    sink_budget[position, ] <- as.numeric(
      Matrix::crossprod(carrier$matrix[, ncol(carrier$matrix), drop = FALSE],
        native)
    )

    scaling <- .population_scale(plan, native, budget_floor)
    admitted[position, ] <- scaling$admitted
    carried <- carried * rep(scaling$scale, each = nodes)
    stack[position, ] <- as.numeric(carried)
    rm(view, native, carried)
  }

  fit <- .population_ols(plan$model, stack)
  shape <- function(values, third, third_name) {
    array(t(values), dim = c(nodes, n_query, length(third)),
      dimnames = stats::setNames(
        list(node_labels, bank$labels, third),
        c("node", "query", third_name)
      ))
  }

  readout <- .population_readout("query_bank", bank$labels, bank$matrix)
  value <- structure(list(
    basis = "query_bank",
    coefficients = shape(fit$coefficients, plan$model$columns, "term"),
    values = shape(stack, labels, "subject"),
    fitted = shape(fit$fitted, labels, "subject"),
    residuals = shape(fit$residuals, labels, "subject"),
    residual_df = as.integer(n_subject - plan$model$rank),
    index = index,
    queries = bank$matrix,
    component = component,
    ledger = unname(.population_ledger_names[[component]]),
    semantics = plan$semantics,
    normalization = plan$normalization,
    uncertainty = .population_uncertainty(plan, bank, index, uncertainty),
    receipt = .population_receipt(
      plan, readout, component, budget_floor, receipts, native_totals,
      sink_budget, admitted, deviations, fit$unresolved
    ),
    scientific_plan_id = .population_result_id(
      plan, readout, component, budget_floor
    )
  ), class = "effect_population_result")
  .validate_population_result(value)
  value
}

# The complete form ------------------------------------------------------------
#
# `estimate_population()` reads `K` declared contrasts; this route reads the
# whole `q`-by-`q` form at every group node. The two share the fit, the
# transport, the budget certificate and the receipt, and differ in exactly one
# thing: the axis the transported geometry is read along.
#
# Why the identity is the right query. Reading a complete packed form *is*
# reading it through the `p`-by-`p` identity, and a coordinate tile is a column
# block of that identity. So this route runs the same fused kernel
# `estimate_population()` runs a contrast bank through, handed a different
# matrix -- the two executors cannot drift, because there is only one execution
# path underneath them. It also means the native `n_i`-by-`p` packed block is
# never allocated: the kernel contracts against the tile it was given.
#
# What that costs, stated plainly: `ceiling(p / tile)` passes over each
# participant's source instead of one. The alternative -- materialize each
# participant to a block store once and read coordinate tiles back off disk --
# buys the passes back and pays `n_i * p` doubles of durable storage per
# participant plus a temporary-file lifetime. This route pays passes and keeps
# the peak; `design/population-form-contract.md` section 5's whole point is
# that the population layer must stream, and streaming is the thing that has to
# hold when `q` is large.

# The packed `svec` coordinate index, as a table rather than as bare labels.
# `symmetric_packed` is the Frobenius-consistent codec of `effect-form-v1`
# section 3 -- lower triangle by column, off-diagonal entries carrying `sqrt(2)`
# -- and section 5 of the population contract records that the `sqrt(2)` is
# load-bearing rather than cosmetic: a naive `upper.tri` packing gets the inner
# product wrong by an `O(1)` amount, not by a tolerance. Carrying the scale
# beside every coordinate is what lets a reader of the result confirm which
# codec produced it without reverse-engineering the widths.
.population_packed_index <- function(effects) {
  q <- length(effects)
  rows <- unlist(lapply(seq_len(q), function(column) column:q),
    use.names = FALSE)
  columns <- rep(seq_len(q), times = q:1)
  data.frame(
    coordinate = paste0(effects[rows], ":", effects[columns]),
    row = as.integer(rows),
    column = as.integer(columns),
    scale = ifelse(rows == columns, 1, sqrt(2)),
    stringsAsFactors = FALSE
  )
}

# One column block of the `p`-by-`p` identity: the physical query whose fused
# execution returns exactly those packed coordinates.
.population_coordinate_query <- function(width, columns, labels) {
  probe <- matrix(0, width, length(columns),
    dimnames = list(NULL, labels[columns]))
  probe[cbind(columns, seq_along(columns))] <- 1
  probe
}

# The tile, from the population plan's own compute policy unless the caller
# names one.
#
# With no declared workspace the default is the compiler's own
# (`.select_compiler_memory_plan()`: `min(64, output_width)`), so a population
# run tiles the coordinate axis the same way a single-participant one does.
# With a declared workspace the tile is the largest that keeps the two arrays
# linear in it -- the `N`-by-`(m + 1) * tile` group stack and the widest
# `n_i`-by-`tile` native block -- inside the budget.
#
# It clamps at one rather than at zero. A budget too small for a single
# coordinate is a real refusal, but it is not this function's to raise: each
# participant's own plan carries its own `compute_policy()`, and the geometry
# compiler already refuses when the conservative memory plan exceeds it. The
# population policy governs the group side, so the honest thing here is to hand
# back the smallest tile that does any work and let the subject-level budget
# speak for the subject-level arrays.
.population_coordinate_tile <- function(plan, width, coordinate_tile) {
  if (!is.null(coordinate_tile)) {
    return(min(.validate_tile_size(coordinate_tile, "coordinate_tile"), width))
  }
  budget <- plan$compute$workspace_bytes
  if (is.null(budget)) return(min(64L, width))
  nodes <- nrow(plan$group_index) + 1L
  native <- max(vapply(plan$subjects, function(subject)
    as.double(subject$measurements), numeric(1)))
  per_coordinate <- 8 * (length(plan$subjects) * nodes + native)
  as.integer(max(1L, min(width, floor(budget / per_coordinate))))
}

# Section 4.1's normalizations are group functionals of a *queried* ledger:
# `T_i = sum_x c_ix` is the native total of one contrast, and section 4.1 is
# explicit that the scalar is per participant and per query. A complete form has
# no query, so the only per-coordinate divisor available is the coordinate's own
# native total -- and that is a weight varying along the coordinate axis, which
# section 3.2 measures breaking query-fit commutation by 25.5 % of the largest
# coefficient. The assembled form would then not be queryable: reading it with a
# contrast would disagree with `estimate_population()` on the same contrast, and
# the object's one advertised property would be false.
.population_admit_complete_form <- function(plan) {
  if (identical(plan$normalization, "none")) return(invisible(NULL))
  .capability_refusal(sprintf(paste0(
    "`normalization = \"%s\"` divides each participant's ledger by the native ",
    "total *of the query being read*, and a complete form is not read through ",
    "a query. The only divisor a form could use is per packed coordinate, ",
    "which is a weight varying along the coordinate axis -- exactly the case ",
    "`population-form-v1` section 3.2 measures breaking the query-fit ",
    "commutation. The assembled form would no longer answer a contrast the ",
    "way `estimate_population()` answers it, so crossform refuses rather than ",
    "shipping a form whose defining property does not hold."
  ), plan$normalization),
    capability = "complete_form_normalization",
    namespace = "population_execution",
    reasons = c(
      "normalization_varies_along_the_coordinate_axis",
      paste0("declared_normalization:", plan$normalization)
    ),
    remedies = c(
      paste0(
        "Plan the population with `normalization = \"none\"`, the mean ",
        "participant ledger in native evidence units, which is a per-",
        "participant scalar of one and commutes with every query."
      ),
      paste0(
        "Keep the normalization and read the contrasts you care about with ",
        "`estimate_population(plan, queries)`, where the divisor is the ",
        "native total of that query and section 4.1 defines it."
      )
    ))
}

# One coordinate tile of the transported group stack ---------------------------
#
# The per-tile subject loop, extracted so that everything which streams the
# complete form reads it through *one* path. `materialize_population()` fits
# the tile; `heterogeneity()` (`R/population-heterogeneity.R`,
# `population-form-v1` sections 5 and 6) residualizes it and accumulates an
# `N`-by-`N` Gram from it. Two loops over the same source would be two chances
# to disagree about what a transported coordinate is, and the Gram's entire
# claim is that it is the second moment of the same numbers the coefficients
# are fitted from.
#
# `state` is what survives a tile: one execution receipt per participant, the
# native totals and sink budget per packed coordinate, the running worst budget
# deviation, and the native row counts the streaming bound quotes. It is
# returned rather than mutated, so the caller owns the accumulation and the
# helper stays a function of its arguments.
.population_tile_state <- function(plan, width, coordinate_labels) {
  labels <- names(plan$subjects)
  n_subject <- length(labels)
  totals <- matrix(NA_real_, n_subject, width,
    dimnames = list(labels, coordinate_labels))
  list(
    labels = labels,
    receipts = stats::setNames(vector("list", n_subject), labels),
    native_totals = totals,
    sink_budget = totals,
    deviations = numeric(n_subject),
    native_rows = integer(n_subject),
    # `TRUE` until the first tile has been streamed: one receipt per
    # participant, taken from that tile. Later tiles read the same source under
    # the same plan and would record the same provenance `ceiling(p / tile)`
    # times over.
    first = TRUE
  )
}

.population_stream_tile <- function(plan, component, columns, width,
                                    coordinate_labels, nodes, state) {
  probe <- .population_coordinate_query(width, columns, coordinate_labels)
  # `N` by `(m + 1) * length(columns)`, rebuilt per tile and dropped with it.
  stack <- matrix(NA_real_, length(state$labels), nodes * length(columns),
    dimnames = list(state$labels, NULL))
  for (position in seq_along(state$labels)) {
    label <- state$labels[[position]]
    carrier <- plan$transport[[label]]
    view <- .run_geometry_compiler(
      plan$subjects[[label]], query = probe, component = component
    )
    native <- view$values
    colnames(native) <- coordinate_labels[columns]
    if (state$first) {
      state$receipts[[label]] <- view$receipt
      state$native_rows[[position]] <- nrow(native)
    }
    carried <- transport_values(carrier, native)
    if (identical(plan$semantics, "budget")) {
      state$deviations[[position]] <- max(state$deviations[[position]],
        .population_budget_certificate(label, native, carried))
    }
    state$native_totals[position, columns] <- colSums(native)
    # The sink row is already in budget units under either semantics:
    # `transport_values()` divides the `m` group rows by transported mass and
    # leaves the sink alone, because a density of unmapped territory has no
    # referent (section 1.3).
    state$sink_budget[position, columns] <- carried[nodes, ]
    stack[position, ] <- as.numeric(carried)
    rm(view, native, carried)
  }
  state$first <- FALSE
  list(stack = stack, state = state)
}

#' Materialize a planned population form at every group node
#'
#' `materialize_population()` is the complete-geometry counterpart of
#' [estimate_population()]. Where that verb reads each participant's geometry
#' through a declared bank of `K` contrasts, this one carries the **whole**
#' \eqn{q \times q}{q x q} form to the group nodes and fits the group model at
#' every packed coordinate, so the result is a coefficient *form* per group node
#' and model term rather than a coefficient per query.
#'
#' It stands to [estimate_population()] exactly as [materialize_geometry()]
#' stands to [evaluate_geometry()]: two verbs, because complete materialization
#' is a different memory regime and a different estimand, and because a caller
#' who forgets a query bank should get [estimate_population()]'s refusal rather
#' than a silently enormous complete form.
#'
#' @section How it streams, and what that bounds:
#' Reading a complete packed form is reading it through the identity, and a
#' coordinate tile is a column block of the identity. For each tile the route
#' opens each participant in turn, contracts that participant's native nodes
#' against the tile (never allocating the full \eqn{n_i \times p}{n_i x p}
#' packed block), transports the tile's rows with that participant's operator,
#' stacks the transported tile across participants, and solves the plan's
#' prefactored QR once for the tile. Coefficient tiles come out and are written
#' into the packed coefficient forms.
#'
#' Nothing of size \eqn{N \times (m+1) \times p}{N x (m+1) x p} is ever
#' allocated. The arrays in flight are the group stack,
#' \eqn{N \times (m+1) \times \mathrm{tile}}{N x (m+1) x tile}, and one native
#' block, \eqn{\max_i n_i \times \mathrm{tile}}{max_i n_i x tile}; the durable
#' output is \eqn{(m+1) \times k \times p}{(m+1) x k x p} for `k` model terms,
#' which is smaller than the refused array whenever the design is not saturated.
#' The receipt records all four numbers under `$streaming`. The cost is
#' \eqn{\lceil p/\mathrm{tile}\rceil}{ceiling(p/tile)} passes over each
#' participant's source instead of one.
#'
#' @section What the result does not carry:
#' There are no `$values`, `$fitted` or `$residuals` arrays. Those are indexed
#' by participant, group node and packed coordinate --- they *are* the
#' \eqn{N \times (m+1) \times p}{N x (m+1) x p} object this route exists in
#' order not to build. Read the participant-level transported values through
#' [estimate_population()] with the contrasts you care about, where the
#' subject axis costs `K` and not `p`.
#'
#' @section Normalization:
#' Only `normalization = "none"` is admitted; see the refusal below. Budget and
#' density semantics are both admitted, and the budget certificate of
#' `population-form-v1` section 2 is asserted per participant and per packed
#' coordinate, exactly as [estimate_population()] asserts it per query.
#'
#' @section Refusals:
#' Each is an `effect_capability_refusal` (see [catch_refusal()]).
#'
#' * `complete_form_normalization` --- a plan declaring `"unit_budget"` or
#'   `"precision_weighted"`. Their divisor is the native total of a query, and
#'   a form is not read through one; a per-coordinate divisor would break the
#'   query-fit commutation the form's queryability rests on.
#' * `reserved_group_index_columns` --- a group index already carrying `sink`
#'   or `units`, as in [estimate_population()].
#'
#' @param plan An `effect_population_plan` from [plan_population()].
#' @param component Which component of each participant's conservative geometry
#'   to carry: `"total"` (the default), `"coherent"`, or `"configuration"`. The
#'   transported object takes the `population-form-v1` section 8.1 name for it,
#'   as in [estimate_population()].
#' @param coordinate_tile Optional positive number of packed coordinates held
#'   in flight at once. `NULL` takes the tile from the plan's
#'   [compute_policy()]. The tile bounds memory and is not part of the
#'   estimand: it does not enter `$scientific_plan_id`, and two tile sizes
#'   agree to [numerical_contract()]`$atol`. Not bit-identically --- the tile
#'   is a blocking choice, and the contract's `bitwise_across_blocking` is
#'   `FALSE` because it changes the column count of one fused contraction.
#' @return An `effect_population_result` with `$basis` `"complete_form"`,
#'   carrying `$coefficient_forms` (a `node`-by-`term`-by-`coordinate` array
#'   whose last axis is the packed `svec` form, so every `[node, term, ]` slice
#'   is one packed \eqn{q \times q}{q x q} coefficient form), `$effects`,
#'   the `$coordinates` table describing the packed codec, the `$index` of
#'   group nodes plus the sink, `$ledger`, and a `$receipt` recording every
#'   participant's read, its transport signature, the budget certificate and
#'   the streaming bound. `as.data.frame()` returns the coefficient table in
#'   long form.
#' @references `design/population-form-contract.md` (`population-form-v1`),
#'   sections 2, 3, 4, 5 and 8.
#' @family population transports
#' @seealso [estimate_population()] for the query-first path over the same
#'   plan, [plan_population()] for the estimand both execute, and
#'   [materialize_geometry()] for the single-participant analogue.
#' @examples
#' # Four participants on different native frames, two shared group nodes.
#' effects <- effect_space(c("face", "house"), basis_id = "population:v1")
#' subject <- function(id, n, gain) {
#'   domain <- abstract_domain(n, coordinates = cbind(x = seq_len(n) - 1),
#'     feature_ids = paste0("f", seq_len(n)), id = id)
#'   values <- function(divisor) matrix(
#'     gain * seq_len(2 * n) / (n * divisor), 2, n,
#'     dimnames = list(c("face", "house"), NULL)
#'   )
#'   rel <- relation(list(run1 = values(1), run2 = values(2)),
#'     effects = effects, domain = domain)
#'   plan_geometry(rel, compile_frame(voxelwise(), domain),
#'     cross_partitions(rel))
#' }
#' carrier <- function(n) anatomical_transport(
#'   native_coords = cbind(seq_len(n) - 1), group_coords = cbind(c(0, 3)),
#'   semantics = "budget"
#' )
#' sizes <- c(s01 = 5L, s02 = 6L, s03 = 7L, s04 = 8L)
#' gains <- c(s01 = 1, s02 = 1.4, s03 = 0.7, s04 = 1.1)
#' subjects <- stats::setNames(lapply(names(sizes), function(id)
#'   subject(id, sizes[[id]], gains[[id]])), names(sizes))
#' plan <- plan_population(subjects, lapply(sizes, carrier))
#'
#' form <- materialize_population(plan)
#' form
#'
#' # The group-mean coefficient form at the first group node: three packed
#' # coordinates, the off-diagonal one carrying the codec's root two.
#' form$coordinates
#' form$coefficient_forms["group1", "(Intercept)", ]
#'
#' # Querying the assembled form with a contrast is the same number the
#' # query-first executor returns for that contrast.
#' contrast <- c(1, -1)
#' packed <- crossform:::.svec_symmetric(tcrossprod(contrast))
#' assembled <- form$coefficient_forms[, "(Intercept)", ] %*% packed
#' queried <- estimate_population(plan, rbind(`face-house` = contrast))
#' max(abs(assembled - queried$coefficients[, "face-house", "(Intercept)"]))
#'
#' # The tile bounds memory and nothing else: one coordinate at a time gives
#' # the same forms, and the same estimand identity.
#' tiled <- materialize_population(plan, coordinate_tile = 1L)
#' max(abs(tiled$coefficient_forms - form$coefficient_forms))
#' identical(tiled$scientific_plan_id, form$scientific_plan_id)
#' tiled$receipt$streaming
#' @export
materialize_population <- function(plan,
                                   component = c("total", "coherent",
                                     "configuration"),
                                   coordinate_tile = NULL) {
  .validate_population_plan(plan, deep = FALSE)
  component <- match.arg(component)
  .population_admit_complete_form(plan)

  effects <- plan$subjects[[1L]]$task$left_relation$effect_space$coordinates
  coordinates <- .population_packed_index(effects)
  width <- nrow(coordinates)
  coordinate_labels <- coordinates$coordinate
  tile <- .population_coordinate_tile(plan, width, coordinate_tile)

  labels <- names(plan$subjects)
  n_subject <- length(labels)
  index <- .population_result_index(plan)
  nodes <- nrow(index)
  terms <- plan$model$columns

  # The one durable array. It is `(m + 1) * k * p` doubles against the
  # `N * (m + 1) * p` the naive route would hold, and `k <= N` because
  # `plan_population()` refuses a rank-deficient design.
  forms <- array(NA_real_, c(nodes, length(terms), width),
    dimnames = stats::setNames(
      list(as.character(index$node), terms, coordinate_labels),
      c("node", "term", "coordinate")
    ))

  unresolved <- 0L
  starts <- .tile_starts(width, tile)
  state <- .population_tile_state(plan, width, coordinate_labels)

  for (start in starts) {
    columns <- start:min(start + tile - 1L, width)
    streamed <- .population_stream_tile(plan, component, columns, width,
      coordinate_labels, nodes, state)
    state <- streamed$state
    fit <- .population_ols(plan$model, streamed$stack, coefficients_only = TRUE)
    unresolved <- unresolved + fit$unresolved
    # `as.numeric(carried)` stacked the tile column-major with the node index
    # fastest, so the fit's response columns unfold as `node` by `coordinate`;
    # the permutation puts the packed axis last, where every slice along it is
    # one form.
    forms[, , columns] <- aperm(
      array(t(fit$coefficients), c(nodes, length(columns), length(terms))),
      c(1L, 3L, 2L)
    )
    rm(streamed, fit)
  }
  receipts <- state$receipts
  native_totals <- state$native_totals
  sink_budget <- state$sink_budget
  deviations <- state$deviations
  native_rows <- state$native_rows

  readout <- .population_readout("complete_form", coordinate_labels)
  value <- structure(list(
    basis = "complete_form",
    coefficient_forms = forms,
    residual_df = as.integer(n_subject - plan$model$rank),
    index = index,
    effects = effects,
    coordinates = coordinates,
    component = component,
    ledger = unname(.population_ledger_names[[component]]),
    semantics = plan$semantics,
    normalization = plan$normalization,
    uncertainty = NULL,
    receipt = .population_receipt(
      plan, readout, component, 0, receipts, native_totals, sink_budget,
      matrix(TRUE, n_subject, width,
        dimnames = list(labels, coordinate_labels)),
      deviations, unresolved,
      streaming = list(
        coordinate_tile = as.integer(tile),
        packed_width = as.integer(width),
        passes_per_subject = length(starts),
        group_stack_doubles = as.double(n_subject) * nodes * tile,
        native_block_doubles = as.double(max(native_rows)) * tile,
        output_doubles = as.double(nodes) * length(terms) * width,
        refused_dense_doubles = as.double(n_subject) * nodes * width
      )
    ),
    scientific_plan_id = .population_result_id(plan, readout, component, 0)
  ), class = "effect_population_result")
  .validate_population_result(value)
  value
}

# The group fit ----------------------------------------------------------------
#
# One backsolve for the whole result. The plan stored the design's pivoted QR
# because the design varies along neither the node axis nor the coordinate
# axis, so the same factorization serves every group node and every query; the
# per-(node, query) loop section 3 describes is the column loop inside LAPACK
# rather than an R one, and `Theta` comes back for all `(m + 1) * K` responses
# at once.
#
# A column carrying a non-finite response is solved around rather than through.
# Two things produce one, and both are measurements: density semantics returns
# `NA` at a group node reached by no native mass (never `0`, because zero is a
# measurement and absence is not), and `unit_budget` marks a participant whose
# signed native total is not bounded away from zero. Dropping the participant
# instead would silently change which population the estimate is about.
#
# `coefficients_only` exists for the streamed complete-form route and for
# nothing else. There the fitted values and residuals are `N` by
# `(m + 1) * tile` each, and retaining them would triple the one array whose
# size the tile is there to bound; the route documents their absence from the
# result rather than paying for them and dropping them.
.population_ols <- function(model, stack, coefficients_only = FALSE) {
  finite <- apply(is.finite(stack), 2L, all)
  coefficients <- matrix(NA_real_, ncol(model$matrix), ncol(stack),
    dimnames = list(model$columns, NULL))
  fitted <- NULL
  residuals <- NULL
  if (!coefficients_only) {
    fitted <- matrix(NA_real_, nrow(stack), ncol(stack))
    residuals <- fitted
  }
  if (any(finite)) {
    responses <- stack[, finite, drop = FALSE]
    coefficients[, finite] <- qr.coef(model$qr, responses)
    if (!coefficients_only) {
      fitted[, finite] <- qr.fitted(model$qr, responses)
      residuals[, finite] <- qr.resid(model$qr, responses)
    }
  }
  list(coefficients = coefficients, fitted = fitted, residuals = residuals,
    unresolved = as.integer(sum(!finite)))
}

# The receipt ------------------------------------------------------------------
#
# What a reader would have to know to say what these numbers mean, and what
# they would have to check to disbelieve them: who was read and under which
# plan identity, how each participant was carried and under what provenance,
# whether the budget closed, which normalization made incommensurable native
# totals commensurable and on what criterion, and --- section 8.1's required
# print line --- which native frame family the ledger belongs to.
.population_receipt <- function(plan, readout, component, budget_floor,
                                receipts, native_totals, sink_budget,
                                admitted, deviations, unresolved,
                                streaming = NULL) {
  # The per-query numbers stay as `subject`-by-`query` matrices rather than
  # becoming list columns of the audit table: a list column is not portable,
  # `location_transport()` refuses one in an index for exactly that reason, and
  # a matrix reads the same way whether the bank holds one query or twenty.
  subjects <- plan$subject_index
  subjects$budget_deviation <- if (identical(plan$semantics, "budget")) {
    deviations
  } else {
    rep(NA_real_, nrow(subjects))
  }
  subjects$contributes <- rowSums(admitted) == ncol(admitted)
  list(
    contract_version = "population-form-v1",
    population_plan_id = plan$scientific_plan_id,
    population_signature = plan$signature,
    subjects = subjects,
    native_total = native_totals,
    sink_budget = sink_budget,
    subject_receipts = receipts,
    component = component,
    ledger = unname(.population_ledger_names[[component]]),
    frame = .population_frame_identity(plan),
    semantics = plan$semantics,
    evaluation_order = plan$fit$evaluation_order,
    fit = plan$fit,
    basis = readout$basis,
    # The axis the transported geometry was read along, named the same way for
    # either executor. `queries` stays beside it and stays a *bank* of declared
    # contrasts: a complete form declares none, and reusing the word for the
    # packed coordinates would let a reader believe a contrast was chosen.
    readout = list(basis = readout$basis, labels = readout$labels,
      width = length(readout$labels)),
    queries = if (identical(readout$basis, "query_bank")) readout$labels,
    # Present only on the streamed complete-form route: the coordinate tile the
    # run was bounded by, and the arrays that bound implies. This is the
    # auditable form of the acceptance claim -- the dense stack the route
    # refuses to build is quoted beside the one it did build.
    streaming = streaming,
    normalization = list(
      mode = plan$normalization,
      # Section 4.3(a): the divisor is read from the same data as the ledger,
      # so it is declared descriptively rather than presented as an
      # independent estimate.
      budget_estimate = if (identical(plan$normalization, "unit_budget")) {
        "same_data_ratio"
      } else {
        NA_character_
      },
      budget_floor = if (identical(plan$normalization, "unit_budget")) {
        as.numeric(budget_floor)
      } else {
        NA_real_
      },
      floor_criterion = if (identical(plan$normalization, "unit_budget")) {
        "abs(native_total) > budget_floor * sum(abs(native_ledger))"
      } else {
        NA_character_
      },
      # One logical per participant and query: `TRUE` where that
      # participant's ledger was admitted a share. Reported as the admissions
      # rather than the exclusions so the empty case reads as "everyone
      # contributed" instead of as an absent field.
      admitted = admitted
    ),
    budget = list(
      asserted = identical(plan$semantics, "budget"),
      tolerance = .population_budget_tolerance,
      scale = "relative_to_ledger_l1_norm",
      max_relative_deviation = if (identical(plan$semantics, "budget")) {
        max(deviations)
      } else {
        NA_real_
      }
    ),
    unresolved_columns = unresolved
  )
}

# Section 8.1's required print line needs the native frame family the ledger
# belongs to. The transports carry per-row frame metadata when the frame was
# built by `frame_family()`; when they do not, the plan's own audit table still
# names the declared normalization, and saying "undeclared" is the honest
# answer rather than inventing one.
.population_frame_identity <- function(plan) {
  families <- unlist(lapply(plan$transport, function(carrier) {
    if ("family" %in% names(carrier$native_index)) {
      unique(as.character(carrier$native_index$family))
    } else {
      NULL
    }
  }), use.names = FALSE)
  scales <- unlist(lapply(plan$transport, function(carrier) {
    if ("scale" %in% names(carrier$native_index)) {
      unique(as.numeric(carrier$native_index$scale))
    } else {
      NULL
    }
  }), use.names = FALSE)
  list(
    family = if (length(families)) sort(unique(families)) else NA_character_,
    scale = if (length(scales)) sort(unique(scales)) else NA_real_,
    normalization = sort(unique(plan$subject_index$declared_normalization)),
    row_metadata = all(plan$subject_index$row_metadata)
  )
}

# The record -------------------------------------------------------------------
#
# FIELD CONTRACT --- `effect_population_result`
#
# `$basis` is the discriminator and must be read before any other field. It is
# one of `.population_result_bases`, and the two bases carry different arrays
# because they read the transported geometry along different axes. Everything
# below is stable; a view built on these names will not have to branch on which
# verb produced the result beyond branching on `$basis`.
#
#   basis = "query_bank"   --- `estimate_population(plan, queries)`
#     $coefficients            node x query x term
#     $values                  node x query x subject  (transported, scaled)
#     $fitted, $residuals      node x query x subject
#     $queries                 K x q contrast matrix, rows named by query
#
#   basis = "complete_form" --- `materialize_population(plan)`
#     $coefficient_forms       node x term x coordinate; the PACKED AXIS IS
#                              LAST, so `[node, term, ]` is one `svec` form
#                              ready for `.unsvec_symmetric(., q)`
#     $effects                 the q effect coordinate names
#     $coordinates             data frame: coordinate, row, column, scale
#                              (`scale` is 1 on the diagonal and sqrt(2) off it,
#                              the Frobenius-consistent codec of
#                              `effect-form-v1` section 3)
#     -- deliberately absent: $values, $fitted, $residuals. Indexed by
#     participant, node and packed coordinate they *are* the N x (m+1) x p
#     array the streamed route exists in order not to build.
#
#   both bases
#     $residual_df             N - rank(X)
#     $index                   group nodes plus the sink row, with `sink` and
#                              `units` columns; row order matches the node axis
#     $component, $ledger      section 8.1 component and its transported name
#     $semantics               "budget" or "density"
#     $normalization           the plan's section 4.1 mode
#     $uncertainty             D8 passthrough or NULL
#     $receipt                 $basis, $readout, $subjects, $native_total,
#                              $sink_budget, $budget, $normalization, $frame,
#                              $fit, $streaming (complete_form only)
#     $scientific_plan_id      "population-sha256:..."

.population_result_fields <- function(basis) {
  switch(basis,
    query_bank = c("basis", "coefficients", "values", "fitted", "residuals",
      "residual_df", "index", "queries", "component", "ledger", "semantics",
      "normalization", "uncertainty", "receipt", "scientific_plan_id"),
    complete_form = c("basis", "coefficient_forms", "residual_df", "index",
      "effects", "coordinates", "component", "ledger", "semantics",
      "normalization", "uncertainty", "receipt", "scientific_plan_id"),
    character()
  )
}

# The axis-shape check, one branch per basis. Both end at the same question:
# do the arrays agree with the group node index, the readout axis and the
# participant list the receipt names.
.validate_population_result_shape <- function(x) {
  if (identical(x$basis, "complete_form")) {
    if (!is.array(x$coefficient_forms) ||
        length(dim(x$coefficient_forms)) != 3L ||
        !is.data.frame(x$coordinates) || !.is_strings(x$effects) ||
        !identical(names(dimnames(x$coefficient_forms)),
          c("node", "term", "coordinate"))) {
      .input_error("Population-result fields are missing or noncanonical.")
    }
    q <- length(x$effects)
    shape <- dim(x$coefficient_forms)
    if (shape[[1L]] != nrow(x$index) ||
        shape[[2L]] > nrow(x$receipt$subjects) ||
        shape[[3L]] != nrow(x$coordinates) ||
        nrow(x$coordinates) != q * (q + 1L) / 2L ||
        !identical(names(x$coordinates),
          c("coordinate", "row", "column", "scale"))) {
      .contract_error(paste0(
        "Population-result coefficient forms do not agree with their group ",
        "node index, packed coordinate table and participant list."
      ))
    }
    return(invisible(x))
  }
  if (!is.array(x$coefficients) || length(dim(x$coefficients)) != 3L ||
      !is.array(x$values) || length(dim(x$values)) != 3L ||
      !is.matrix(x$queries)) {
    .input_error("Population-result fields are missing or noncanonical.")
  }
  shape <- dim(x$coefficients)
  if (shape[[1L]] != nrow(x$index) || shape[[2L]] != nrow(x$queries) ||
      !identical(dim(x$values), dim(x$fitted)) ||
      !identical(dim(x$values), dim(x$residuals)) ||
      dim(x$values)[[1L]] != shape[[1L]] ||
      dim(x$values)[[2L]] != shape[[2L]] ||
      dim(x$values)[[3L]] != nrow(x$receipt$subjects)) {
    .contract_error(paste0(
      "Population-result arrays do not agree with their group node index, ",
      "query bank and participant list."
    ))
  }
  invisible(x)
}

.validate_population_result <- function(x) {
  if (!inherits(x, "effect_population_result")) {
    .input_error(sprintf(paste0(
      "Expected an `effect_population_result` from `estimate_population()` ",
      "or `materialize_population()`; received %s."
    ), .msg_value(x)), arg = "x", received = .msg_value(x),
      expected = paste0("an `effect_population_result` from ",
        "`estimate_population()` or `materialize_population()`"))
  }
  if (!is.list(x) || !.is_string(x$basis) ||
      !x$basis %in% .population_result_bases ||
      !.sealed_fields(x, "effect_population_result",
        .population_result_fields(x$basis)) ||
      !is.data.frame(x$index) ||
      !.is_string(x$component) || !.is_string(x$ledger) ||
      !.is_string(x$semantics) || !.is_string(x$normalization) ||
      !is.integer(x$residual_df) || !is.list(x$receipt) ||
      !identical(x$receipt$basis, x$basis) ||
      !.strong_sha256(sub("^population-", "", x$scientific_plan_id))) {
    .input_error("Population-result fields are missing or noncanonical.")
  }
  .validate_population_result_shape(x)
  if (!identical(unname(.population_ledger_names[[x$component]]), x$ledger)) {
    .contract_error(paste0(
      "A transported component must carry the ledger name ",
      "`population-form-v1` section 8.1 assigns it; this result names a ",
      "different one."
    ))
  }
  if (sum(x$index$sink) != 1L ||
      !identical(as.character(x$index$node[x$index$sink]),
        .transport_sink_label)) {
    .contract_error(paste0(
      "A population result carries exactly one sink row, materialized even ",
      "when it is empty."
    ))
  }
  invisible(x)
}
