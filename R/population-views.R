# Views on an estimated population form ---------------------------------------
#
# Layer 5 (results / views). `population-driver.R` computes the group numbers;
# this file reads them. Nothing here executes anything: every value it returns
# is a fixed linear combination of coefficients `estimate_population()` already
# produced, which is why the whole file is arithmetic over `$coefficients` and
# a set of refusals for the combinations that arithmetic cannot honestly reach.
#
# The one fact that makes these verbs possible is `population-form-v1`
# section 3. Query, transport and the group fit act on three different tensor
# factors, so under the OLS default they commute: the fitted coefficient of a
# combined query is the same combination of the fitted coefficients of the
# bank's queries. A population view is therefore a *selection or recombination
# in query space*, exact rather than approximate, and it needs no second pass
# over anybody's data.
#
# Three consequences shape the file.
#
#  1. **The bank is a basis, and a basis has a span.** A result estimated
#     through `K` contrasts holds the group coefficients of `K` packed
#     operators `svec(c_k c_k^T)`. A view is derivable exactly when its own
#     packed operator lies in the span of those `K`, and the test is a least
#     squares solve with a residual check --- not a set membership test on the
#     contrast vectors, which would refuse `2c` when `c` is in the bank and
#     would miss that the full pairwise bank reaches every RDM and every RSA
#     regression. Outside the span the verb refuses and names the remedy, and
#     it never re-runs anything. What "outside" means is a relative residual
#     threshold, stated at `.population_span_tolerance` below with what it does
#     and does not buy.
#
#     The bank of all `q(q-1)/2` pairwise difference contrasts spans exactly
#     the quadratic forms of zero-sum contrasts --- so `rdm()` and `rsa()` are
#     always derivable from it, and no uncentred contrast ever is. That is the
#     same fact that makes a centered Gram recoverable from an RDM, read one
#     level up.
#
#  2. **The response is not linear along the query axis under every
#     normalization.** Under `normalization = "unit_budget"` each query column
#     carries its *own* divisor --- `T_ik` is the native total of query `k`'s
#     ledger, and a bank of `K` queries has `K` of them
#     (`population-form-v1` section 4.1, and `.population_scale()` in the
#     driver). Under `"none"` nothing is rescaled and every combination is
#     exact; the transport's `semantics` does not enter this at all, because a
#     density divisor `P^T mu` is a property of the transport and not of the
#     query. Under `"unit_budget"` a mixture of columns carries no denominator
#     and is refused. A *multiple* of one estimated query is the exception and
#     is exact --- `s c` has ledger `s c` over total `s T` --- so it is
#     admitted with the weight forced to one, which is the whole of
#     `.population_require_query_linearity()`.
#
#  3. **The sink is a row of every view.** Section 1.1 requires it
#     materialized, section 3.3 fits it and forbids reporting it as a value at
#     a location, and section 2's budget law is what makes the aggregate add
#     up only when it is present. It is carried, labelled `<sink>`, and
#     `contribution()` gives it its own group rather than letting it join a
#     territory.
#
# Labelling is `population-form-v1` section 8.1 throughout: a transported
# component is a `native_coherent_ledger` or a `native_configuration_ledger`,
# and no value, column name or printed label of a population view carries the
# bare name `coherent` or `configuration`. Two refusal messages do quote
# `component = "coherent"`, and that is the scope of the rule rather than an
# exception to it: section 8.1 forbids the bare name *as a name for a
# transported quantity*, and `component` there is the estimation-time argument
# a caller has to type to get a different run.
#
# The two bases, and why only one of them can refuse. `$basis` is the
# discriminator of the field contract in `population-driver.R`, and every verb
# below reads the result through `.population_query_basis()` and
# `.population_basis_slice()` rather than through `$coefficients` or
# `$coefficient_forms` directly.
#
#   `"query_bank"`     --- `estimate_population(plan, queries)`. The basis is
#                          the `K` packed operators the bank named, and a view
#                          outside their span is refused.
#   `"complete_form"`  --- `materialize_population(plan)`. The basis is the
#                          whole packed coordinate system, so `alpha` *is* the
#                          packed target, every symmetric query is in the span
#                          by construction, and no span refusal can fire. The
#                          `p`-by-`p` identity is never built: the solver
#                          short-circuits, which is what keeps the route
#                          usable at the `q` the streamed executor exists for.
#
# So the same four verbs read either result, and the only difference a caller
# sees is that a complete form never asks them to re-estimate.

# The span test's tolerance, relative to the target operator's own norm. It is
# looser than the `1e-12` of section 11 on purpose: those are exact identities
# asserted on a computed number, and this is a *rank* question answered by a
# least squares solve, whose residual for a target genuinely in the span is
# machine epsilon times the basis condition number.
#
# It is a threshold and not a proof, and the honest statement of what it buys
# is this: a target whose relative residual is at or below `1e-9` is reported
# as the least squares fit in the span, which for a target genuinely in the
# span is the target and for a target `1e-10` away from it is the projection.
# Measured on a pairwise bank over three effects, an uncentred `(1, 1, 0)`
# misses by `9.4e-1` and a near-zero-sum `(1, -1, 1e-8)` by `5.8e-9` --- the
# separation is decisive for the queries anyone asks and is not decisive
# against a query deliberately built to sit just outside. Nothing downstream
# treats it as a proof of exactness: the route is recorded, and a caller who
# needs an exact query rather than a near one estimates it.
.population_span_tolerance <- 1e-9

# The selection tolerance: when one entry of an `alpha` column counts as the
# only one, and when a weight counts as one. A view whose packed target *is* a
# basis column takes an exact unit weight rather than a least squares
# approximation to it, so that a `unit_budget` selection reads its column
# bit-exactly rather than through `1 - 3e-16`. The support test is relative to
# the column's own largest weight; the unit test is absolute, because `1` is a
# constant and not a scale.
.population_selection_tolerance <- 1e-12

# The basis ---------------------------------------------------------------

# What the result actually holds a group coefficient for, in packed
# coordinates. For a query bank that is one column per query,
# `svec(c_k c_k^T)`, built by the same expression `.population_query_bank()`
# uses in the driver so that a view and the run it reads agree on what a query
# *is*. For a complete form it is the packed coordinate system itself, and
# `$packed` is `NULL` because the identity is never materialized.
.population_query_basis <- function(x) {
  if (identical(x$basis, "complete_form")) {
    return(list(
      kind = "complete_form",
      packed = NULL,
      labels = as.character(x$coordinates$coordinate),
      effects = x$effects,
      width = nrow(x$coordinates)
    ))
  }
  effects <- colnames(x$queries)
  width <- length(effects) * (length(effects) + 1L) / 2L
  packed <- matrix(vapply(seq_len(nrow(x$queries)), function(k) {
    .svec_symmetric(tcrossprod(x$queries[k, ]))
  }, numeric(width)), nrow = width,
    dimnames = list(NULL, rownames(x$queries)))
  list(
    kind = "query_bank",
    packed = packed,
    labels = rownames(x$queries),
    effects = effects,
    # The number of basis columns, not the packed width: what a reader needs
    # to know is how many estimated things a view is a combination of.
    width = ncol(packed)
  )
}

# Which single basis column a view column is reached through, or `NA` when it
# needs more than one. The threshold is relative to the column's own largest
# weight, so the answer does not change when a view is scaled.
.population_column_support <- function(weights) {
  scale <- max(abs(weights))
  if (!is.finite(scale) || scale <= 0) return(NA_integer_)
  nonzero <- which(abs(weights) > .population_selection_tolerance * scale)
  if (length(nonzero) == 1L) nonzero else NA_integer_
}

# The three ways a query-bank view reaches its query, in increasing order of
# what they assume. A *selection* reads one estimated column as it stands. A
# *scaled selection* reads one estimated column at a different scale --- the
# same functional, a different multiple of it. A *recombination* mixes
# columns, and is exact only where the response is linear along the query
# axis (`.population_require_query_linearity()`).
.population_basis_route <- function(basis, alpha) {
  if (identical(basis$kind, "complete_form")) return("coordinate_form")
  if (is.null(alpha)) return("whole_bank")
  support <- vapply(seq_len(ncol(alpha)), function(column) {
    .population_column_support(alpha[, column])
  }, integer(1))
  if (anyNA(support)) return("recombination")
  weights <- vapply(seq_along(support), function(column) {
    alpha[support[[column]], column]
  }, numeric(1))
  if (all(abs(weights - 1) <= .population_selection_tolerance)) {
    "selection"
  } else {
    "scaled_selection"
  }
}

# One node-by-basis matrix: the group coefficients of `term` along whichever
# axis the result was read on. This and `.population_query_basis()` are the
# only two places that touch the basis-specific arrays of the field contract.
.population_basis_slice <- function(x, term) {
  slice <- if (identical(x$basis, "complete_form")) {
    array(x$coefficient_forms[, term, ], dim(x$coefficient_forms)[c(1L, 3L)])
  } else {
    array(x$coefficients[, , term], dim(x$coefficients)[1:2])
  }
  # The basis labels and the readout axis of the coefficient array are the
  # same axis in the same order, and every combination below is positional, so
  # the correspondence is asserted rather than assumed. The driver guarantees
  # it today; a reshaped result would otherwise combine the right numbers
  # under the wrong names.
  if (!identical(dimnames(slice)[[2L]], NULL) ||
      ncol(slice) != .population_readout_width(x)) {
    .contract_error(paste0(
      "A population result's readout axis does not agree with the basis its ",
      "views are read through."
    ))
  }
  dimnames(slice) <- list(as.character(x$index$node), NULL)
  slice
}

# The width of the axis a result was read along, from the field contract's own
# discriminator.
.population_readout_width <- function(x) {
  if (identical(x$basis, "complete_form")) return(nrow(x$coordinates))
  nrow(x$queries)
}

# One packed symmetric operator per column, from one contrast vector per
# column. `svec` is the package's Frobenius-consistent codec, so the inner
# product of a packed operator with a packed form is the bilinear query it
# names -- the `sqrt(2)` on the off-diagonals is load-bearing here exactly as
# `population-form-v1` section 5 says it is.
.population_packed_targets <- function(contrasts, labels) {
  width <- ncol(contrasts) * (ncol(contrasts) + 1L) / 2L
  matrix(vapply(seq_len(nrow(contrasts)), function(k) {
    .svec_symmetric(tcrossprod(contrasts[k, ]))
  }, numeric(width)), nrow = width, dimnames = list(NULL, labels))
}

# Solve `packed %*% alpha = targets`, one column of `alpha` per view column,
# and refuse the columns that miss. An exact basis column is matched before
# the solve so that a selection is a selection and not a least squares
# approximation to one (see `.population_selection_tolerance`).
.population_span_solve <- function(x, basis, targets, what) {
  if (identical(basis$kind, "complete_form")) {
    # In the packed coordinate basis every symmetric operator *is* its own
    # coefficient vector: the residual is identically zero, no query can be
    # outside the span, and building the `p`-by-`p` identity to discover that
    # would cost `p^2` doubles at the one `q` where `p` is large.
    dimnames(targets) <- list(basis$labels, colnames(targets))
    return(targets)
  }
  packed <- basis$packed
  scale <- pmax(sqrt(colSums(targets^2)), .Machine$double.eps)
  alpha <- matrix(0, ncol(packed), ncol(targets),
    dimnames = list(basis$labels, colnames(targets)))
  pending <- logical(ncol(targets))
  for (column in seq_len(ncol(targets))) {
    distance <- sqrt(colSums((packed - targets[, column])^2)) / scale[[column]]
    hit <- which(distance <= .population_selection_tolerance)
    if (length(hit) == 1L) {
      alpha[hit[[1L]], column] <- 1
    } else {
      pending[[column]] <- TRUE
    }
  }
  if (any(pending)) {
    solved <- qr.coef(qr(packed), targets[, pending, drop = FALSE])
    solved[is.na(solved)] <- 0
    alpha[, pending] <- solved
  }
  residual <- sqrt(colSums((targets - packed %*% alpha)^2)) / scale
  outside <- colnames(targets)[residual > .population_span_tolerance]
  if (length(outside)) {
    .population_span_refusal(x, basis, what, outside, residual)
  }
  alpha
}

.population_span_refusal <- function(x, basis, what, outside, residual) {
  .capability_refusal(sprintf(paste0(
    "This population was estimated through a bank of %s (%s), and %s of %s ",
    "is not a combination of them (%s): %s. A population view is a selection ",
    "or recombination in query space -- exact, because the query, the ",
    "transport and the group fit commute -- and there is nothing in the ",
    "result to recombine into a query the bank never named. crossform will ",
    "not re-run anybody's geometry behind a reader verb, and it will not ",
    "project the query onto the span and report the projection under the ",
    "name of the query."
  ), .msg_count(length(basis$labels), "query", "queries"),
    .msg_names(basis$labels),
    .msg_count(length(outside), "column"), what, .msg_names(outside),
    sprintf("worst relative residual %s", format(signif(max(residual), 3L)))),
    capability = "population_query_in_bank",
    namespace = "population_views",
    reasons = c(
      "query_outside_estimated_bank",
      paste0("query_not_derivable:", outside)
    ),
    remedies = c(
      paste0(
        "Re-estimate the population with the contrast in the bank: ",
        "`estimate_population(plan, queries = rbind(<the bank>, ",
        "<the contrast>))`. The bank is the estimand, so widening it is a ",
        "new run rather than a new reading of this one."
      ),
      paste0(
        "A bank of all pairwise difference contrasts spans exactly the ",
        "quadratic forms of zero-sum contrasts, so it reaches every RDM and ",
        "every RSA regression over these effects and no uncentred contrast ",
        "at all. If uncentred contrasts are the question, estimate them."
      )
    ))
}

# Where a combination along the basis axis is exact, and what to do with the
# one case where it is exact for a different reason than the solver assumed.
#
# `"none"` rescales nothing, so the response is linear in the query operator
# and every `alpha` is exact. `"unit_budget"` divides each participant's
# ledger by *that query's* own native total (`population-form-v1` section 4.1,
# `.population_scale()` in the driver), so the divisor varies along the query
# axis and a mixture of columns is not the estimate of the mixed query.
#
# A *scalar multiple* of one estimated query is the exception, and it is exact
# rather than merely admissible. The requested query `H = s H_k` has native
# ledger `s c_k` and native total `s T_ik`, so its unit-budget response is
# `s b / (s T) = b / T` --- exactly column `k`, for every nonzero `s`, and
# with the same budget-floor admission because `|sT| > f * sum|s c|` iff
# `|T| > f * sum|c|`. So the weight is forced to one rather than the `s`
# the span solve returned: `s` would be right under `"none"` and is wrong
# here, and reporting it would be reporting `s` times a share under the name
# of a share. That is the whole content of the refusal below, applied instead
# of refusing where the answer is known.
#
# The gate is an allowlist. `precision_weighted` is refused by
# `plan_population()` today and so cannot reach here, but if that gate ever
# lifts its per-participant weights would rescale the response again, and a
# denylist would admit the recombination silently.
.population_require_query_linearity <- function(x, alpha, what) {
  if (identical(x$normalization, "none")) return(alpha)
  if (!identical(x$normalization, "unit_budget")) {
    .capability_refusal(sprintf(paste0(
      "A population view combines estimated columns, which is exact only ",
      "where the response is linear along the basis axis. This population ",
      "was fitted under `normalization = \"%s\"`, and no view knows what ",
      "that rescaling does to a combination."
    ), x$normalization),
      capability = "population_view_query_linearity",
      namespace = "population_views",
      reasons = paste0("unhandled_normalization:", x$normalization),
      remedies = paste0(
        "Plan the population with `normalization = \"none\"`, where the ",
        "response is linear along the basis axis and every combination is ",
        "exact."
      ))
  }
  support <- vapply(seq_len(ncol(alpha)), function(column) {
    .population_column_support(alpha[, column])
  }, integer(1))
  if (anyNA(support)) {
    combined <- colnames(alpha)[is.na(support)]
    .capability_refusal(sprintf(paste0(
      "This population was fitted under `normalization = \"unit_budget\"`, ",
      "which divides each participant's ledger by *that query's* own native ",
      "total. The divisor therefore varies along the query axis, the ",
      "response is not linear in the query operator, and %s of %s (%s) ",
      "would mix columns carrying different divisors --- a number with no ",
      "denominator, presented as a share. A query that is a multiple of one ",
      "estimated query is exact and is admitted; a mixture of several is ",
      "not."
    ), .msg_count(length(combined), "column"), what, .msg_names(combined)),
      capability = "population_view_query_linearity",
      namespace = "population_views",
      reasons = c(
        "unit_budget_divisor_is_query_specific",
        paste0("query_is_a_recombination:", combined)
      ),
      remedies = c(
        paste0(
          "Plan the population with `normalization = \"none\"`, the mean ",
          "subject ledger in native evidence units, where the response is ",
          "linear along the query axis and every combination is exact."
        ),
        paste0(
          "Re-estimate with the contrast itself in the bank, so the view ",
          "selects one estimated query rather than mixing several."
        )
      ))
  }
  selected <- alpha * 0
  selected[cbind(support, seq_along(support))] <- 1
  selected
}

# The term ----------------------------------------------------------------

# Which group coefficient the view reads. `$coefficients` is
# node-by-query-by-term, and a view is a view of one term: the group mean
# under the default `~ 1`, a group difference or a slope under anything else.
# There is no default beyond the single-term case, because picking one of
# several silently would put an unnamed estimand on the page.
.population_view_term <- function(x, term) {
  terms <- if (identical(x$basis, "complete_form")) {
    dimnames(x$coefficient_forms)[[2L]]
  } else {
    dimnames(x$coefficients)[[3L]]
  }
  if (is.null(term)) {
    if (length(terms) == 1L) return(terms[[1L]])
    .input_error(sprintf(paste0(
      "`term` is required: the group model fits %s (%s), and a population ",
      "view reads one of them at every group node. Which coefficient is ",
      "part of what the number means, so it is named rather than defaulted."
    ), .msg_count(length(terms), "term"), .msg_names(terms)),
      arg = "term", received = "no argument",
      expected = .msg_names(terms))
  }
  if (.is_string(term)) {
    position <- match(term, terms)
    if (is.na(position)) {
      .input_error(sprintf(
        "`term = \"%s\"` names no column of the group model; it fits %s.",
        term, .msg_names(terms)),
        arg = "term", received = sprintf("\"%s\"", term),
        expected = .msg_names(terms))
    }
    return(terms[[position]])
  }
  if (is.numeric(term) && length(term) == 1L && !is.na(term) &&
      term %% 1 == 0 && term >= 1 && term <= length(terms)) {
    return(terms[[as.integer(term)]])
  }
  .input_error(sprintf(paste0(
    "`term` must name one column of the group model (%s) or give its ",
    "position; received %s."
  ), .msg_names(terms), .msg_value(term)),
    arg = "term", received = .msg_value(term),
    expected = .msg_names(terms))
}

# Combine the estimated query columns at one term. Exact zeros are dropped
# rather than multiplied in: an unresolved node-query cell is `NA` --- density
# semantics at a group node reached by no native mass, or a `unit_budget`
# participant below the floor --- and multiplying it by zero would spread that
# absence across every view column that does not use it.
.population_combine <- function(x, term, alpha) {
  slice <- .population_basis_slice(x, term)
  values <- vapply(seq_len(ncol(alpha)), function(column) {
    used <- which(alpha[, column] != 0)
    if (!length(used)) return(rep(0, nrow(slice)))
    drop(slice[, used, drop = FALSE] %*% alpha[used, column])
  }, numeric(nrow(slice)))
  matrix(values, nrow(slice), ncol(alpha),
    dimnames = list(as.character(x$index$node), colnames(alpha)))
}

# The record --------------------------------------------------------------

.population_view_fields <- c(
  "values", "view", "term", "ledger", "native_ledger", "semantics",
  "normalization", "index", "columns", "query", "receipt",
  "scientific_plan_id", "metadata"
)

# What a view's identity is a hash *of*: the packed operators it reads and the
# names it reads them under. Not the arguments that produced them. A contrast
# energy is quadratic in its contrast, so `c(1, -1, 0)` and `c(-1, 1, 0)` are
# one view and would be two identities under a hash of the weights; the same
# holds for an RSA design that differs only in which predictors were called
# nuisance. The operator is what the result is queried with, so it is what the
# identity is over.
.population_view_descriptor <- function(targets) {
  list(operators = unname(targets), columns = colnames(targets))
}

# A view's own scientific identity: the population plan it descends from, the
# result that estimated it, the reader verb, the ledger it is a ledger of, the
# term it reads, and the operators it reads with. It mirrors
# `.geometry_view_scientific_id()` --- two ways of reaching the same view hash
# identically, and two different views never do.
.population_view_scientific_id <- function(x, view, term, descriptor) {
  .sha256_signature(list(
    schema_version = 1L,
    contract_version = "population-form-v1",
    role = "population_view",
    population_plan_id = x$receipt$population_plan_id,
    parent = x$scientific_plan_id,
    view = view,
    ledger = x$ledger,
    term = term,
    descriptor = descriptor
  ), "population-sha256:")
}

# Everything a reader would have to know to say what a population view means:
# what the transport did (semantics, provenance, cross-fit), which native
# frame family the ledger belongs to and therefore which "coherent" this is a
# ledger of (section 8.1's required print line), whether the budget closed,
# which normalization made incommensurable participants commensurable, and how
# the view's own query was reached from the estimated basis
# (`.population_basis_route()`).
.population_view_receipt <- function(x, view, term, basis, alpha, identity) {
  list(
    contract_version = "population-form-v1",
    view = view,
    population_plan_id = x$receipt$population_plan_id,
    population_result_id = x$scientific_plan_id,
    ledger = x$ledger,
    native_ledger = !identical(x$ledger, "transported_total"),
    term = term,
    transport = list(
      semantics = x$semantics,
      provenance = sort(unique(as.character(x$receipt$subjects$provenance))),
      cross_fit = sort(unique(stats::na.omit(
        as.character(x$receipt$subjects$cross_fit)
      )))
    ),
    frame = x$receipt$frame,
    subjects = x$receipt$subjects,
    normalization = x$receipt$normalization,
    budget = x$receipt$budget,
    fit = x$receipt$fit,
    evaluation_order = x$receipt$evaluation_order,
    # How the view's own query was reached from the estimated basis, and
    # `$route` is the field to read first. The coefficients are recorded for a
    # query bank, where `K` is small by construction and the numbers are the
    # whole audit. They are `NULL` in two cases, both of which `$route` names:
    # `"coordinate_form"`, where `alpha` *is* the packed query already carried
    # by `$query` and a `p`-by-column copy of it would put the largest array
    # in the run inside a receipt; and `"whole_bank"`, where the view is the
    # basis itself and the coefficients would be an identity.
    basis = list(
      kind = basis$kind,
      width = basis$width,
      queries = if (identical(basis$kind, "query_bank")) basis$labels,
      coefficients = if (identical(basis$kind, "query_bank")) alpha,
      route = .population_basis_route(basis, alpha),
      tolerance = .population_span_tolerance
    ),
    scientific_plan_id = identity
  )
}

.new_population_view <- function(x, view, values, columns, term, query,
                                 basis, alpha, descriptor, index = NULL,
                                 metadata = list()) {
  if (is.null(index)) index <- x$index
  identity <- .population_view_scientific_id(x, view, term, descriptor)
  value <- structure(list(
    values = values,
    view = view,
    term = term,
    ledger = x$ledger,
    native_ledger = !identical(x$ledger, "transported_total"),
    semantics = x$semantics,
    normalization = x$normalization,
    index = index,
    columns = columns,
    query = query,
    receipt = .population_view_receipt(x, view, term, basis, alpha, identity),
    scientific_plan_id = identity,
    metadata = metadata
  ), class = "effect_population_view")
  .validate_population_view(value)
  value
}

.validate_population_view <- function(x) {
  if (!inherits(x, "effect_population_view")) {
    .input_error(sprintf(paste0(
      "Expected an `effect_population_view` from a reader verb on an ",
      "`effect_population_result`; received %s."
    ), .msg_value(x)), arg = "x", received = .msg_value(x),
      expected = "an `effect_population_view`")
  }
  if (!.sealed_fields(x, "effect_population_view", .population_view_fields) ||
      !is.matrix(x$values) || !is.numeric(x$values) ||
      !.is_string(x$view) || !.is_string(x$term) || !.is_string(x$ledger) ||
      !.is_flag(x$native_ledger) || !.is_string(x$semantics) ||
      !.is_string(x$normalization) || !is.data.frame(x$index) ||
      !is.data.frame(x$columns) || !is.list(x$receipt) ||
      !.strong_sha256(sub("^population-", "", x$scientific_plan_id))) {
    .input_error("Population-view fields are missing or noncanonical.")
  }
  if (nrow(x$values) != nrow(x$index) ||
      ncol(x$values) != nrow(x$columns) ||
      !identical(colnames(x$values), as.character(x$columns$column))) {
    .contract_error(paste0(
      "Population-view values do not agree with their group node index and ",
      "column table."
    ))
  }
  if (!x$ledger %in% .population_ledger_names ||
      !identical(x$native_ledger,
        !identical(x$ledger, "transported_total"))) {
    .contract_error(paste0(
      "A population view must carry the transported ledger name ",
      "`population-form-v1` section 8.1 assigns the component it reads, and ",
      "the flag saying whether that name is a native-frame ledger."
    ))
  }
  if (sum(x$index$sink) != 1L ||
      !identical(as.character(x$index$node[x$index$sink]),
        .transport_sink_label)) {
    .contract_error(paste0(
      "A population view carries exactly one sink row, materialized even ",
      "when it is empty."
    ))
  }
  invisible(x)
}

# The verbs ---------------------------------------------------------------

#' Reader verbs on an estimated population form
#'
#' `contrast_energy()`, `rdm()`, `rsa()` and `contribution()` are S3 generics.
#' Applied to an [estimate_population()] result they read the *group*
#' coefficients that run already produced: a population view is a fixed linear
#' combination of estimated query columns, not a second execution, and no
#' participant's geometry is opened again.
#'
#' @section Why a combination is exact:
#' The query, the transport and the group fit act on the experimental,
#' spatial and participant axes respectively, so under the OLS default they
#' commute (`population-form-v1` section 3). The fitted coefficient of a
#' combined query is therefore the same combination of the fitted
#' coefficients of the bank's queries, exactly rather than approximately.
#'
#' @section What the estimated basis reaches:
#' A result from [estimate_population()] (`$basis` `"query_bank"`) holds the
#' group coefficient of `K` packed operators
#' \eqn{\mathrm{svec}(c_k c_k^\top)}{svec(c_k c_k^T)}. A view is derivable
#' exactly when its own packed operator lies in their span, which is a wider
#' set than the bank's contrast vectors: `2c` is derivable when `c` is in the
#' bank, and the bank of all \eqn{q(q-1)/2}{q(q-1)/2} pairwise difference
#' contrasts spans **exactly** the quadratic forms of zero-sum contrasts. So a
#' full pairwise bank reaches every [rdm()] and every [rsa()] regression over
#' those effects, and reaches no uncentred contrast at all. (Other banks reach
#' other things: three uncentred contrasts can reach a fourth. The span is the
#' criterion; the pairwise case is the one worth stating because it is the one
#' `rdm()` and `rsa()` rest on.) Outside the span the verb refuses and names
#' the re-estimation remedy rather than reporting the nearest thing it can
#' reach under the query's name.
#'
#' A result from [materialize_population()] (`$basis` `"complete_form"`) holds
#' the assembled packed coefficient form, whose basis is the whole coordinate
#' system. Every symmetric query is in its span by construction, so these
#' verbs query the form directly and no span refusal can fire.
#'
#' @section `unit_budget` admits one estimated query at a time:
#' Under `normalization = "unit_budget"` each query column is divided by *that
#' query's* own native total (`population-form-v1` section 4.1), so the
#' divisor varies along the query axis and the response is not linear in the
#' query operator. A view that reads **one** estimated query is exact and is
#' admitted, including a multiple of one: `s H_k` has ledger `s c_k` over
#' total `s T_{ik}`, so its share is column `k` for every nonzero `s`, and the
#' weight is set to one rather than to the \eqn{s}{s} a span solve returns. A
#' view that **mixes** estimated queries is refused: the mixture carries no
#' denominator, and a share of nothing is not a share.
#'
#' @section The sink:
#' Every view carries the sink as a row of `$index`, labelled `<sink>`, marked
#' by `$index$sink` and reported in budget units. It is never a value at a
#' location; it is there because [contribution()]'s ledger adds up only when
#' unmapped territory is a number in the table.
#'
#' @section Labelling:
#' A transported component is a `native_coherent_ledger` or a
#' `native_configuration_ledger`, never a bare `coherent` or `configuration`
#' (`population-form-v1` section 8.1). `$ledger` carries the name, `$receipt`
#' carries it beside the native frame family the ledger belongs to and the
#' transport that carried it, and no output column or printed label of a
#' population view uses the bare name.
#' `native_coherent_ledger + native_configuration_ledger = transported_total`
#' holds exactly; what fails is reading the coherent ledger as a group-node
#' common mode.
#'
#' @section Refusals:
#' Each is an `effect_capability_refusal` in namespace `"population_views"`
#' (see [catch_refusal()]).
#'
#' * `population_query_in_bank` --- the view's packed operator is outside the
#'   span of the estimated bank. Remedy: re-estimate with the contrast in the
#'   bank.
#' * `population_view_query_linearity` --- a view mixing estimated queries on
#'   a `unit_budget` result, or any normalization no view knows how to combine
#'   under.
#' * `population_component_fixed_at_estimation` --- `component`, which is an
#'   argument of [estimate_population()] and not of a view.
#' * `budget_semantics` --- [contribution()] on a `"density"` population.
#'   Density satisfies no conservation law, so a territory sum of densities is
#'   a share of nothing.
#' * `conservative_frame` --- [contribution()] on a population planned with
#'   `allow_nonconservative = TRUE` over a frame that is not column
#'   normalized.
#' * `sink_is_not_a_territory` --- a `by` that labels a group node `<sink>`.
#'   The sink is appended as its own row automatically.
#' * `nondestructive_decomposition` and `guaranteed_psd` --- `remove_univariate`
#'   and `normalize`, refused for the same reasons the per-participant views
#'   refuse them.
#'
#' @param x An `effect_population_result` from [estimate_population()] or
#'   [materialize_population()], or -- for [contribution()] -- an
#'   `effect_population_view` from one of the other three verbs.
#' @param term Which column of the group model to read, by name or position.
#'   Required when the group model fits more than one; the default `~ 1` fits
#'   only `(Intercept)`, the group mean, and needs no argument.
#' @param weights One finite contrast weight per experimental effect, named or
#'   in the effect order of `$queries`.
#' @param pairs Optional two-column matrix of effect pairs, by name or index.
#'   The default reports every unordered pair.
#' @param models,nuisance Named model and nuisance dissimilarity matrices over
#'   the experimental effects, as [rsa()] takes them.
#' @param intercept Whether the RSA design carries an intercept column.
#' @param component Must be omitted. The geometry component is fixed when the
#'   population is estimated -- it is `component =` on
#'   [estimate_population()] -- because a transported ledger of one component
#'   holds no information about another.
#' @param normalize Must be omitted; refused exactly as [rdm()] refuses it.
#' @param remove_univariate Must be omitted or `FALSE`; refused exactly as
#'   [contrast_energy()] refuses it.
#' @param by The grouping of group nodes, as one label per group node (the
#'   sink excluded, since it is always its own group) or one column name of
#'   `using` or of the view's `$index`.
#' @param using Optional per-group-node metadata table naming the territories
#'   a view cannot carry itself. It is **joined on its `node` column**, which
#'   must hold the group node identifiers of `$index`, so row order does not
#'   matter and a table without that column is refused rather than bound by
#'   position. Meaningful only when `by` names one of its columns.
#' @param ... Unused; present so the methods match their generics.
#' @return An `effect_population_view`: `$values`, one row per group node plus
#'   the sink and one column per view column; `$columns` naming those columns;
#'   `$index`, the group node table with the sink marked and its units; the
#'   `$ledger` name, `$term`, `$semantics`, `$normalization`; a `$receipt`
#'   carrying the transport, the native frame family, the budget certificate,
#'   the normalization and the basis coefficients that reached this view; and
#'   a `$scientific_plan_id` derived from the population plan and the view's
#'   own parameters. `as.data.frame()` returns the long table.
#' @references `design/population-form-contract.md` (`population-form-v1`),
#'   sections 2, 3, 4 and 8.
#' @seealso [estimate_population()] and [materialize_population()] for the
#'   runs these read, and [contrast_energy()], [rdm()], [rsa()] and
#'   [contribution()] for the per-participant verbs of the same names.
#' @family population transports
#' @name population_views
#' @examples
#' # Three participants on different native frames, two group nodes, and a
#' # bank spanning all three pairwise contrasts.
#' effects <- effect_space(c("face", "house", "tool"), basis_id = "popview:v1")
#' subject <- function(id, n) {
#'   domain <- abstract_domain(n, coordinates = cbind(x = seq_len(n) - 1),
#'     feature_ids = paste0("f", seq_len(n)), id = id)
#'   values <- function(divisor) matrix(
#'     seq_len(3 * n) / (n * divisor), 3, n,
#'     dimnames = list(c("face", "house", "tool"), NULL)
#'   )
#'   rel <- relation(list(run1 = values(1), run2 = values(1.4)),
#'     effects = effects, domain = domain)
#'   plan_geometry(rel, compile_frame(voxelwise(), domain),
#'     cross_partitions(rel))
#' }
#' carrier <- function(n) anatomical_transport(
#'   native_coords = cbind(seq_len(n) - 1), group_coords = cbind(c(0, 4)),
#'   semantics = "budget"
#' )
#' sizes <- c(s01 = 5L, s02 = 6L, s03 = 7L)
#' subjects <- stats::setNames(
#'   lapply(names(sizes), function(id) subject(id, sizes[[id]])), names(sizes)
#' )
#' fit <- estimate_population(
#'   plan_population(subjects, lapply(sizes, carrier)),
#'   rbind(`face-house` = c(1, -1, 0), `face-tool` = c(1, 0, -1),
#'     `house-tool` = c(0, 1, -1))
#' )
#'
#' # One contrast, selected out of the estimated bank.
#' contrast_energy(fit, c(face = 1, house = -1, tool = 0))
#'
#' # The three distances, and an RSA regression on them: both are inside the
#' # span of a full pairwise bank, so neither needs a second run.
#' rdm(fit)$values
#' animacy <- matrix(c(0, 0, 1, 0, 0, 1, 1, 1, 0), 3, 3,
#'   dimnames = list(c("face", "house", "tool"), c("face", "house", "tool")))
#' as.data.frame(rsa(fit, models = list(animacy = animacy)))
#'
#' # The ledger adds up over the group nodes and the sink.
#' ledger <- contribution(fit, by = c("anterior", "posterior"))
#' ledger
#'
#' # An uncentred contrast is outside the span of a zero-sum bank, and is
#' # refused with the re-estimation remedy rather than projected onto it.
#' refusal <- catch_refusal(contrast_energy(fit, c(1, 1, 0)))
#' refusal$capability
NULL

#' @rdname population_views
#' @export
contrast_energy.effect_population_result <- function(x, weights,
                                                     remove_univariate = FALSE,
                                                     term = NULL, ...) {
  .check_no_extra_arguments("contrast_energy", ...)
  .validate_population_result(x)
  if (missing(weights)) {
    .input_error(paste0(
      "`weights` is required: pass one finite weight per experimental ",
      "effect, for example `contrast_energy(fit, c(face = 1, house = -1))`. ",
      "The contrast must be reachable from the bank the population was ",
      "estimated through."
    ),
      arg = "weights", received = "no argument",
      expected = "one finite weight per experimental effect")
  }
  if (!isFALSE(remove_univariate)) {
    .population_refuse_remove_univariate()
  }
  basis <- .population_query_basis(x)
  weights <- .align_contrast(weights, basis$effects)
  term <- .population_view_term(x, term)
  targets <- .population_packed_targets(
    matrix(weights, nrow = 1L), "energy"
  )
  alpha <- .population_span_solve(x, basis, targets, "this contrast energy")
  alpha <- .population_require_query_linearity(
    x, alpha, "this contrast energy"
  )
  .new_population_view(
    x, "contrast_energy", .population_combine(x, term, alpha),
    data.frame(column = "energy", role = "contrast", stringsAsFactors = FALSE),
    term, weights, basis, alpha,
    descriptor = .population_view_descriptor(targets)
  )
}

#' @rdname population_views
#' @export
rdm.effect_population_result <- function(x, component = NULL, pairs = NULL,
                                         normalize = NULL, term = NULL, ...) {
  .check_no_extra_arguments("rdm", ...)
  .validate_population_result(x)
  .population_refuse_component(component, "rdm")
  if (!is.null(normalize)) {
    .population_refuse_normalize()
  }
  basis <- .population_query_basis(x)
  term <- .population_view_term(x, term)
  query <- .pair_difference_query(basis$effects, pairs = pairs)
  contrasts <- matrix(0, length(query$pair_left), length(basis$effects))
  contrasts[cbind(seq_along(query$pair_left), query$pair_left)] <- 1
  contrasts[cbind(seq_along(query$pair_right), query$pair_right)] <- -1
  targets <- .population_packed_targets(contrasts, query$pair_labels)
  alpha <- .population_span_solve(x, basis, targets, "these distances")
  alpha <- .population_require_query_linearity(x, alpha, "these distances")
  .new_population_view(
    x, "rdm", .population_combine(x, term, alpha),
    data.frame(
      column = query$pair_labels, role = "distance",
      left = basis$effects[query$pair_left],
      right = basis$effects[query$pair_right],
      stringsAsFactors = FALSE
    ),
    term, query, basis, alpha,
    descriptor = .population_view_descriptor(targets)
  )
}

#' @rdname population_views
#' @export
rsa.effect_population_result <- function(x, models, nuisance = NULL,
                                         intercept = TRUE, component = NULL,
                                         term = NULL, ...) {
  .check_no_extra_arguments("rsa", ...)
  .validate_population_result(x)
  .population_refuse_component(component, "rsa")
  if (missing(models)) {
    .input_error(paste0(
      "`models` is required: pass one dissimilarity matrix over the ",
      "experimental effects, or a named list of them, for example ",
      "`rsa(fit, models = list(category = m))`."
    ),
      arg = "models", received = "no argument",
      expected = "one dissimilarity matrix, or a named list of them")
  }
  if (!.is_flag(intercept)) {
    .input_error(sprintf("`intercept` must be TRUE or FALSE; received %s.",
      .msg_value(intercept)),
      arg = "intercept", received = .msg_value(intercept),
      expected = "TRUE or FALSE")
  }
  basis <- .population_query_basis(x)
  term <- .population_view_term(x, term)
  models <- .validate_rdm_models(models, basis$effects, "models")
  nuisance <- .validate_rdm_models(nuisance, basis$effects, "nuisance")
  if (any(names(models) %in% names(nuisance))) {
    .input_error(sprintf(paste0(
      "Model and nuisance names must be distinct; %s appears in both. Each ",
      "name becomes one coefficient column."
    ), .msg_names(intersect(names(models), names(nuisance)))),
      arg = "nuisance",
      received = .msg_names(intersect(names(models), names(nuisance))),
      expected = "names distinct from `models`")
  }
  predictors <- c(models, nuisance)
  design <- do.call(cbind, lapply(predictors, .rdm_vector))
  colnames(design) <- names(predictors)
  roles <- c(rep("model", length(models)), rep("nuisance", length(nuisance)))
  if (intercept) {
    design <- cbind(`(Intercept)` = 1, design)
    roles <- c("intercept", roles)
  }
  qr_design <- qr(design, LAPACK = FALSE)
  if (qr_design$rank != ncol(design)) {
    .input_error(.rsa_rank_deficiency_message(design, qr_design, intercept),
      arg = "models", received = "a rank-deficient design",
      expected = "model and nuisance RDMs of full column rank")
  }
  coefficient_map <- .thin_qr_coefficient_map(qr_design)
  rownames(coefficient_map) <- colnames(design)
  # The RSA coefficient is a fixed linear map of the pair distances, and each
  # pair distance is a fixed packed operator, so an RSA term is one packed
  # operator too --- the whole regression is a change of basis in query space
  # and never a second read of anybody's geometry.
  query <- .pair_difference_query(basis$effects, coefficients = coefficient_map)
  contrasts <- matrix(0, length(query$pair_left), length(basis$effects))
  contrasts[cbind(seq_along(query$pair_left), query$pair_left)] <- 1
  contrasts[cbind(seq_along(query$pair_right), query$pair_right)] <- -1
  targets <- .population_packed_targets(contrasts, query$pair_labels) %*%
    t(coefficient_map)
  colnames(targets) <- colnames(design)
  alpha <- .population_span_solve(x, basis, targets, "these RSA coefficients")
  alpha <- .population_require_query_linearity(
    x, alpha, "these RSA coefficients"
  )
  .new_population_view(
    x, "rsa", .population_combine(x, term, alpha),
    data.frame(column = colnames(design), role = roles,
      stringsAsFactors = FALSE),
    term, query, basis, alpha,
    descriptor = .population_view_descriptor(targets)
  )
}

.population_refuse_component <- function(component, verb) {
  if (is.null(component)) return(invisible(NULL))
  .capability_refusal(sprintf(paste0(
    "`%s()` on a population result takes no `component`: which part of each ",
    "participant's conservative geometry was carried is fixed when the ",
    "population is estimated (`component =` on `estimate_population()`), and ",
    "a transported ledger of one component holds no information about ",
    "another. Reading a second component is a second run."
  ), verb),
    capability = "population_component_fixed_at_estimation",
    namespace = "population_views",
    reasons = "component_is_an_estimation_time_choice",
    remedies = paste0(
      "Call `estimate_population(plan, queries, component = \"coherent\")` ",
      "(or `\"configuration\"`) and read the view off that result."
    ))
}

.population_refuse_normalize <- function() {
  .capability_refusal(paste0(
    "`rdm()` reports signed squared distances and will not apply ",
    "correlation-style normalization, at the group level no less than at ",
    "the participant level: crossvalidated diagonal estimates can be zero ",
    "or negative, so dividing by them is not conventional `1 - Pearson` ",
    "distance. A transported ledger inherits the sign, and a group fit of ",
    "signed values is signed too."
  ),
    capability = "guaranteed_psd",
    namespace = "population_views",
    reasons = "signed_cross_generalized_diagonals",
    remedies = "Use the signed squared-distance RDM.")
}

.population_refuse_remove_univariate <- function() {
  .capability_refusal(paste0(
    "`contrast_energy()` does not remove univariate signal, and a population ",
    "result could not honour the request in any case: the component that was ",
    "carried is fixed at estimation time, and the transported ledger holds ",
    "no separable common spatial mode to delete. Estimate the population on ",
    "the component you mean instead."
  ),
    capability = "nondestructive_decomposition",
    namespace = "population_views",
    reasons = "univariate_removal_changes_estimand",
    remedies = paste0(
      "Choose the component with `estimate_population(plan, queries, ",
      "component = \"configuration\")`, whose transported ledger is the ",
      "orthogonal remainder."
    ))
}

# `contribution()` over group nodes ---------------------------------------
#
# D4's aggregation, one level up. The machinery is literally D4's --- the
# metadata join, the grouping resolver, the empty/`NA` group refusal, the
# column sums, the provenance record --- because the operation is the same
# operation: a conservative frame partitions one fixed budget, transport
# carries that partition onto the group nodes, and adding shares up over a
# territory is reading the ledger.
#
# Two things are new at the group level.
#
#  1. **The sink is always its own group.** Section 3.3 forbids reporting it
#     as a value at a location, and section 2's law is what makes the
#     aggregate exact: `sum_u Theta_{j,u}` over the group nodes *and the sink*
#     is the coefficient of the fit on the participants' native totals. Merge
#     the sink into a territory and the territory stops being a place; drop it
#     and the ledger stops adding up. So it is carried, alone, and a grouping
#     that tries to name it is refused.
#  2. **Density conserves nothing.** Section 1.3 measures it: the density
#     columns of the contract's own fixture satisfy no conservation law at
#     all. A territory sum of densities is a share of nothing, so
#     `contribution()` refuses a `"density"` population outright rather than
#     printing a number whose documented property is false --- the same
#     judgement `.population_admit_normalization()` makes in the driver.

.population_require_budget <- function(x) {
  if (identical(x$semantics, "budget")) return(invisible(NULL))
  .capability_refusal(paste0(
    "`contribution()` adds shares of one fixed budget up over a territory, ",
    "and this population was transported under `semantics = \"density\"`. ",
    "Density is transported budget per unit transported territory; it ",
    "satisfies no conservation law (`population-form-v1` section 1.3, ",
    "measured), so a territory sum of densities is not a share of anything ",
    "and the group nodes' sums do not add back to any total."
  ),
    capability = "budget_semantics",
    namespace = "population_views",
    reasons = c(
      "density_conserves_no_budget",
      "declared_semantics:density"
    ),
    remedies = paste0(
      "Build every participant's transport with `semantics = \"budget\"` ",
      "and re-plan the population; the sink then accounts for the territory ",
      "the density normalization was hiding."
    ))
}

.population_require_conservative <- function(x) {
  declared <- as.character(x$receipt$frame$normalization)
  declared <- declared[!is.na(declared)]
  if (length(declared) && all(declared == "conservative")) {
    return(invisible("conservative"))
  }
  .capability_refusal(sprintf(paste0(
    "`contribution()` needs conservative participant frames, and this ",
    "population was planned over %s. A locally normalized frame is a ",
    "detection map: every node reports the mean evidence density inside its ",
    "own support, overlapping neighbourhoods double-count the features they ",
    "share, and transporting those values onto group nodes carries the ",
    "double counting with them. Only a column-normalized conservative frame ",
    "partitions one fixed global budget, which is what makes a territory sum ",
    "a ledger entry."
  ), if (!length(declared)) {
    "frames whose normalization it does not record"
  } else {
    sprintf("frames declaring %s", .msg_names(declared))
  }),
    capability = "conservative_frame",
    namespace = "population_views",
    reasons = if (!length(declared)) {
      "frame_normalization_not_recorded"
    } else {
      paste0("frame_is_not_column_normalized:", setdiff(declared,
        "conservative"))
    },
    remedies = paste0(
      "Re-plan the population over frames built with `normalization = ",
      "\"conservative\"` and without `allow_nonconservative = TRUE`, then ",
      "re-estimate."
    ))
}

# The whole estimated basis at one term, as a view. `contribution()` on a
# population *result* aggregates this: one column per estimated query, or one
# per packed coordinate on a complete form, which is the widest ledger the
# result can honestly show. Section 2's law holds columnwise on either axis --
# the argument is `P 1 = 1` and never touches which functional is carried --
# so the aggregate adds up in both cases.
#
# The values come straight off the basis slice rather than through an identity
# `alpha`: on a complete form that identity is `p`-by-`p`.
.population_ledger_view <- function(x, term) {
  basis <- .population_query_basis(x)
  values <- .population_basis_slice(x, term)
  colnames(values) <- basis$labels
  form <- identical(basis$kind, "complete_form")
  .new_population_view(
    x, "ledger", values,
    data.frame(column = basis$labels,
      role = if (form) "coordinate" else "query",
      stringsAsFactors = FALSE),
    term, if (form) x$coordinates else x$queries, basis, NULL,
    descriptor = list(
      basis = basis$kind,
      labels = basis$labels,
      queries = if (form) NULL else unname(x$queries)
    )
  )
}

.population_group_label <- function(expression) {
  if (is.symbol(expression)) as.character(expression) else "group"
}

# `.contribution_scientific_id()`'s content, under the population namespace
# prefix and with the membership added. Two things differ from D4's version,
# and both are deliberate.
#
#  * The prefix is not reusable: `population-sha256:` is what marks a value as
#    living under `population-form-v1`, and every population record's
#    validator pins it.
#  * The *assignment* is hashed, not only the group keys. `by = c("a","a","b")`
#    and `by = c("a","b","b")` name the same two territories over the same
#    nodes and are two different aggregates; an identity that could not tell
#    them apart would let one be presented under the other's id.
.population_contribution_id <- function(parent, label, keys, members) {
  .sha256_signature(list(
    schema_version = 1L,
    contract_version = "population-form-v1",
    role = "population_contribution",
    parent = parent,
    aggregated_by = label,
    groups = as.character(keys),
    members = lapply(members, as.character)
  ), "population-sha256:")
}

# D4's `.contribution_using_table()` joins on a `measurement` column, because
# a per-measurement view's identity column is called `measurement`. A group
# node's is called `node`, and an atlas table naming group nodes is the whole
# point of `using` here -- so `node` is renamed to the column the join looks
# for, and a table carrying neither is refused rather than joined by row
# order. Row-order joins are how a permuted atlas becomes a wrong aggregate
# with no error anywhere, which is the same reason `plan_population()` refuses
# to bind covariates by position.
.population_using_keys <- function(using, by) {
  if (is.null(using)) return(NULL)
  if (!.is_string(by)) {
    .input_error(paste0(
      "`using` is only meaningful when `by` names one of its columns. `by` ",
      "was given as a grouping vector, which already says which territory ",
      "every group node belongs to."
    ),
      arg = "using", received = "a table alongside a grouping vector",
      expected = "`by` as a column name, or `using = NULL`")
  }
  if (!is.data.frame(using)) return(using)
  if (!"node" %in% names(using)) {
    .input_error(sprintf(paste0(
      "`using` must name the group nodes it describes, in a `node` column ",
      "carrying the identifiers of the population's own group node index. ",
      "Received a table of %s. crossform will not join per-node metadata by ",
      "row order: a permuted atlas would produce a wrong territory ledger ",
      "with nothing anywhere to see. (A frame family's `$index` is keyed by ",
      "`measurement` and indexes *native* nodes, so it is not the table this ",
      "argument takes.)"
    ), .msg_names(names(using))),
      arg = "using", received = .msg_names(names(using)),
      expected = "a `node` column naming each group node")
  }
  if ("measurement" %in% names(using)) {
    .input_error(paste0(
      "`using` carries both `node` and `measurement`, and crossform will not ",
      "guess which one names the group nodes. Keep `node` and rename or drop ",
      "the other."
    ),
      arg = "using", received = "both `node` and `measurement`",
      expected = "a `node` column and no `measurement` column")
  }
  names(using)[names(using) == "node"] <- "measurement"
  using
}

.population_contribution <- function(x, by, using, label) {
  .validate_population_view(x)
  .population_require_budget(x)
  normalization <- .population_require_conservative(x)
  index <- x$index
  sink <- which(index$sink)
  nodes <- which(!index$sink)
  ids <- as.character(index$node[nodes])
  table <- .contribution_using_table(.population_using_keys(using, by), ids)
  resolved <- .contribution_grouping(
    by, label, table, list(index = index[nodes, , drop = FALSE]), length(nodes)
  )
  if (any(as.character(resolved$values) %in% .transport_sink_label)) {
    .capability_refusal(sprintf(paste0(
      "The grouping `%s` labels a group node `%s`, which is the sink's own ",
      "label. Unmapped native territory is carried as its own row and is ",
      "never merged into a territory: it is a mass-accounting column, not a ",
      "value at a location (`population-form-v1` section 3.3), and giving it ",
      "a place would put budget that reached nowhere inside somewhere."
    ), resolved$label, .transport_sink_label),
      capability = "sink_is_not_a_territory",
      namespace = "population_views",
      reasons = "grouping_names_the_sink",
      remedies = paste0(
        "Rename that group; the sink is appended as its own row of the ",
        "aggregate automatically, so the ledger still adds up."
      ))
  }
  grouped <- .contribution_groups(resolved$values, resolved$label)
  rows <- c(
    split(nodes, grouped$groups),
    stats::setNames(list(sink), .transport_sink_label)
  )
  keys <- c(as.character(grouped$keys), .transport_sink_label)
  aggregate_index <- data.frame(
    node = keys,
    sink = c(rep(FALSE, length(rows) - 1L), TRUE),
    units = "budget",
    n_nodes = as.integer(lengths(rows, use.names = FALSE)),
    stringsAsFactors = FALSE, check.names = FALSE
  )
  metadata <- x$metadata
  metadata$aggregated_from <- x$scientific_plan_id
  # `budget_exact` holds for every component, because section 2's argument is
  # `P 1 = 1` and never touches which ledger is carried. What the native
  # ledgers give up is comparability, not exactness: a coherent budget is a
  # share of *that frame's* coherent mass, so it is recorded frame-relative
  # exactly as `contribution()` records it at the participant level.
  # The sink is a row of the aggregate but not a measurement of it: section
  # 3.3 forbids reading it as a value at a location, so the provenance counts
  # the group nodes that carry one.
  metadata$aggregation <- .contribution_provenance(
    resolved$label, keys, rows[-length(rows)], normalization,
    budget_exact = "values",
    frame_relative_components = if (x$native_ledger) x$ledger else character(),
    masked = character()
  )
  metadata$aggregation$sink_group <- .transport_sink_label
  metadata$aggregation$source_view <- x$view
  identity <- .population_contribution_id(
    x$scientific_plan_id, resolved$label, keys,
    lapply(rows, function(subset) index$node[subset])
  )
  value <- structure(list(
    values = .contribution_group_column_sums(x$values, rows),
    view = "contribution",
    term = x$term,
    ledger = x$ledger,
    native_ledger = x$native_ledger,
    semantics = x$semantics,
    normalization = x$normalization,
    index = aggregate_index,
    columns = x$columns,
    query = x$query,
    receipt = utils::modifyList(x$receipt, list(
      view = "contribution",
      aggregated_from = x$scientific_plan_id,
      scientific_plan_id = identity
    )),
    scientific_plan_id = identity,
    metadata = metadata
  ), class = "effect_population_view")
  .validate_population_view(value)
  value
}

#' @rdname population_views
#' @export
contribution.effect_population_result <- function(x, by, using = NULL,
                                                  term = NULL, ...) {
  .check_no_extra_arguments("contribution", ...)
  .validate_population_result(x)
  label <- .population_group_label(substitute(by))
  if (missing(by)) {
    .population_missing_grouping()
  }
  .population_contribution(
    .population_ledger_view(x, .population_view_term(x, term)),
    by, using, label
  )
}

#' @rdname population_views
#' @export
contribution.effect_population_view <- function(x, by, using = NULL, ...) {
  .check_no_extra_arguments("contribution", ...)
  label <- .population_group_label(substitute(by))
  if (missing(by)) {
    .population_missing_grouping()
  }
  .population_contribution(x, by, using, label)
}

.population_missing_grouping <- function() {
  .input_error(paste0(
    "`by` is required: name the territory to add the transported ledger up ",
    "over, either as one label per group node (`contribution(fit, by = ",
    "network)`) or as a column of the group node index (`contribution(fit, ",
    "by = \"network\", using = atlas)`). The sink is appended as its own row ",
    "and is never part of a territory."
  ),
    arg = "by", received = "no argument",
    expected = "one group label per group node, or a column name")
}
