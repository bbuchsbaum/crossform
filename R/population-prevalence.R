# Population prevalence --------------------------------------------------------
#
# Layer 5 (results / views). `population_prevalence()` reads an
# `estimate_population()` result and reports two fractions over participants:
# how many of them carry a transported ledger value above a threshold at each
# group node and query, and how many of them point the same way as the rest of
# the group across the whole query bank at each node.
#
# Its edges all point down or sideways: `population-driver.R` (layer 4) for the
# result validator, the group index and the section 8.1 ledger names;
# `check.R`, `conditions.R`, `message-helpers.R` and `primitives.R` (layer 1)
# for the guards, the refusals and the packed codec the query-bank Gram is
# computed in; and `print-methods.R` (layer 5, sideways) for the printing
# primitives and for section 8.1's required frame and transport lines. It also
# reads `.latent_reading_line` from `latent.R` (layer 5, sideways) -- a
# constant, so the file graph records no edge for it, which is exactly why the
# coupling is written down here. Like `population-heterogeneity.R` it defines
# its own `print`, `format` and `as.data.frame` so that neither printing file
# has to refer back here and the file call graph stays acyclic.
#
# WHY THIS IS A LATENT-LAYER RECORD, THOUGH IT PROJECTS NOTHING.
#
# `design/conservative-geometry-contract.md` section 6 confines fractions,
# cumulative curves and effective counts to a declared projection of the signed
# estimates, and `design/population-form-contract.md` section 6.5 applies that
# one level up. A prevalence fraction takes no eigenvalue truncation, so it is
# not a PSD projection -- but it is the same *kind* of operation and it fails
# in the same direction. `value > threshold` is a per-participant sign clamp:
# it discards the magnitude of a crossvalidated estimate and keeps the sign,
# which is the part of it that carries the noise. Three consequences, and each
# is on the record rather than in a footnote:
#
#  1. THE NULL IS ONE HALF, NOT ZERO. A group node and query at which nothing
#     reproduces gives every participant an independent coin flip, so the
#     fraction concentrates at `0.5`. A reader who compares a prevalence
#     against zero is reading a noise floor as a signal. `$reference` carries
#     the number and the print states it.
#  2. NO ERROR BAR FOLLOWS FROM IT. A count of participants is not an estimate
#     of a population quantity here: the participants are not a sample from a
#     declared superpopulation of transports, and the binomial variance of a
#     thresholded crossvalidated estimate is not the sampling variance of
#     anything the plan named. This file computes no standard error, no
#     interval and no p-value, and `population_uncertainty()` -- the separate
#     inferential layer, with its own measured calibration record -- is
#     untouched by constructing one of these.
#  3. THE THRESHOLD IS EXACT AND ABSOLUTE. The comparison is strictly greater,
#     with no relative tolerance, which is the guard every per-node fraction in
#     the package uses and the open contract decision
#     (`conservative-geometry-v1` section 11.4 gap G3) that D7 also inherits
#     rather than settling privately. A cell that is zero only to within
#     rounding is counted by its rounding.
#
# WHAT IS NOT COMPUTABLE HERE, AND WHY IT IS A RECORD RATHER THAN AN ERROR.
#
# The natural form-space reading -- the fraction of participants whose
# transported *form* has positive Frobenius inner product with the group mean
# form at a node -- needs per-participant transported forms, and no population
# result ships them. `estimate_population()` retains each participant's value
# only in the `K` coordinates the declared query bank names;
# `materialize_population()` retains no participant axis at all, because
# bounding the `N`-by-`(m+1)`-by-`p` array is the whole point of the streamed
# route. So the refusal is a property of the shipped objects and not of the
# arguments, and it travels on every record at
# `$receipt$form_prevalence_refusal` exactly as `heterogeneity()`'s
# `$receipt$latent_refusal` does: a reader learns that the form-space reading
# exists and why this record does not carry it.

.population_prevalence_measures <- c("sign", "alignment")

# What the readout inner product is, and when it is a Frobenius inner product.
#
# Participant `i`'s value at query `k` is `<h_k, y_i>` with `h_k` the packed
# `svec(w_k w_k^T)` the driver lowers the contrast to and `y_i` the packed
# transported form. Stacking the bank as `H` (`K` by `p`), the profile is
# `v_i = H y_i` and `<v_i, v_j> = y_i^T H^T H y_j`. That is the Frobenius inner
# product of the two forms *restricted to the bank's span* exactly when the
# packed queries are orthonormal, i.e. when `H H^T = I_K`; otherwise it is a
# Euclidean inner product in coordinates the bank chose, which is a perfectly
# well defined thing to count signs of and is not a statement about geometry.
# The departure is measured and reported rather than assumed either way: a
# two-condition contrast bank is not orthonormal (`||w w^T||_F = ||w||^2`), and
# a reader who is told "alignment" without being told which inner product will
# supply the geometric one from memory.
.population_prevalence_readout_gram <- function(weights) {
  packed <- vapply(seq_len(nrow(weights)), function(k) {
    .svec_symmetric(tcrossprod(weights[k, ]))
  }, numeric(ncol(weights) * (ncol(weights) + 1L) / 2L))
  gram <- crossprod(matrix(packed, ncol = nrow(weights)))
  list(gram = gram,
    deviation = max(abs(gram - diag(nrow(weights)))))
}

# Which queries to report. The same shape as `population_uncertainty()`'s
# `term`, because a caller who has learned one selector should not have to
# learn a second, and the same refusal text: a selector that silently returned
# nothing would produce an empty record rather than a question.
.population_prevalence_query <- function(x, query) {
  labels <- dimnames(x$values)[[2L]]
  if (is.null(query)) return(labels)
  if (is.numeric(query)) {
    positions <- suppressWarnings(as.integer(query))
    if (anyNA(positions) || any(positions != query) || any(positions < 1L) ||
        any(positions > length(labels))) {
      .input_error(sprintf(paste0(
        "`query` given as positions must index the %d queries the population ",
        "was read through."
      ), length(labels)), arg = "query", received = .msg_value(query),
        expected = sprintf("whole numbers in 1:%d", length(labels)))
    }
    return(labels[sort(unique(positions))])
  }
  if (!.is_strings(query) || !all(query %in% labels)) {
    .input_error(sprintf(paste0(
      "`query` must name queries from the bank the population was read ",
      "through: %s."
    ), .msg_names(labels)), arg = "query", received = .msg_value(query),
      expected = paste0("one or more of ", .msg_names(labels)))
  }
  labels[sort(unique(match(query, labels)))]
}

# The refusals -----------------------------------------------------------------

.population_prevalence_refuse_basis <- function() {
  .capability_refusal(paste0(
    "A prevalence is a count over participants, and the streamed ",
    "complete-form route retains no participant axis: it carries the group ",
    "model's coefficient forms and nothing indexed by participant, because ",
    "bounding the array indexed by participant, group node and packed ",
    "coordinate is the reason that route exists. The coefficients it does ",
    "carry cannot say how many participants stood behind them."
  ),
    capability = "population_participant_values",
    namespace = "population_prevalence",
    reasons = "streamed_route_retains_no_participant_values",
    remedies = paste0(
      "Re-read the same population through `estimate_population(plan, ",
      "queries)` with the contrasts you want prevalence for. That route ",
      "keeps `$values`, indexed by group node, query and participant, which ",
      "is the array a prevalence is counted over."
    ))
}

# `term` exists on this verb only to be refused, on the precedent of
# `remove_univariate`, `normalize` and `component` on a population view: a
# caller arriving from `population_uncertainty(x, term = )` will reach for it,
# and under a silent `...` the reach would return a number answering a
# different question.
.population_prevalence_refuse_term <- function() {
  .capability_refusal(paste0(
    "A prevalence has no group-model term axis. `$values` holds one ",
    "transported ledger per participant, not one per column of the group ",
    "model: the terms are properties of the fit across participants, and no ",
    "participant carries a separate value for each of them. A count of ",
    "participants indexed by term would have to invent the decomposition it ",
    "was indexing."
  ),
    capability = "participant_term_decomposition",
    namespace = "population_prevalence",
    reasons = "participant_values_carry_no_term_axis",
    remedies = c(
      paste0(
        "Drop `term`. `query = ` selects the axis participant values do ",
        "have, and the record reports every group node and query on it."
      ),
      paste0(
        "For the term-indexed group coefficients and their error bars, read ",
        "`population_uncertainty(x, term = )`. That is the inferential ",
        "layer and it is deliberately a different object: a prevalence ",
        "carries no standard error, interval or p-value."
      )
    ))
}

# Section 4.3's divisor is a *signed* native total, so `unit_budget` rescales a
# participant whose total is negative by a negative number and flips the sign
# of every value it carries. Each participant's share of its own budget is
# still well defined -- that is what the normalization means -- but the sign of
# one participant's share and the sign of another's are then not the same
# question, and a count that mixes them is a count of nothing. This is refused
# rather than marked, because unlike a coverage floor there is no reading of
# the number that survives the mixture.
.population_prevalence_refuse_orientation <- function(flipped) {
  .capability_refusal(sprintf(paste0(
    "A `unit_budget` population divides each participant's ledger by that ",
    "participant's own *signed* native total (`population-form-v1` section ",
    "4.3), and %s has a negative total. Their transported values are ",
    "therefore sign-flipped relative to the others, so counting how many ",
    "participants are positive counts two different questions at once: ",
    "\"above zero\" means \"above the group\" for one participant and ",
    "\"below\" for another."
  ), .msg_count(length(flipped), "participant")),
    capability = "comparable_ledger_orientation",
    namespace = "population_prevalence",
    reasons = c("unit_budget_divisor_negative",
      paste0("participant:", flipped)),
    remedies = c(
      paste0(
        "Plan the population with `normalization = \"none\"`, which is the ",
        "mean subject ledger in native evidence units and rescales no ",
        "participant, then re-read it."
      ),
      paste0(
        "The signed totals are on `x$receipt$native_total`, one row per ",
        "participant and one column per query, so which participants are ",
        "affected and by how much is readable without re-running anything."
      )
    ))
}

# The form-space reading, refused as a record. It is built here rather than
# inline so the record, the print and the documentation quote one object.
.population_prevalence_form_refusal <- function(basis) {
  list(
    capability = "participant_form_prevalence",
    namespace = "population_prevalence",
    reasons = c(
      "participant_transported_forms_not_retained",
      paste0("basis:", basis)
    ),
    remedies = c(
      paste0(
        "Read `$alignment`, which is the same question asked in the ",
        "coordinates the query bank does retain: a participant's whole ",
        "profile across the bank against the rest of the group's. It is a ",
        "Frobenius inner product of the transported forms restricted to the ",
        "bank's span exactly when the packed queries are orthonormal, and ",
        "`$alignment$readout_gram_deviation` says how far this bank is from ",
        "that."
      ),
      paste0(
        "For a form-space reading with an estimator behind it, use ",
        "`heterogeneity(plan, nodes = )`. Its `$loadings` give every ",
        "participant a signed weight on each mode of the subject Gram, in ",
        "packed form coordinates. It answers about deviations from the group ",
        "fit rather than about agreement with the group mean, and it has to ",
        "stream every participant twice to do it."
      )
    ),
    message = paste0(
      "The fraction of participants whose transported *form* has positive ",
      "Frobenius inner product with the group mean form is refused: no ",
      "population result carries per-participant transported forms. The ",
      "query-bank route retains each participant's value only in the ",
      "coordinates the declared bank names, and the streamed complete-form ",
      "route retains no participant axis at all."
    )
  )
}

# The arithmetic ---------------------------------------------------------------

# Sign prevalence, per group node and query.
#
# The denominator is the number of participants with a *finite* value at the
# cell, counted per cell rather than fixed at `N`: density semantics returns
# `NA` at a group node reached by no native mass and `unit_budget` returns
# `NA` for a participant whose divisor was not admitted, and dividing by `N`
# through either would report a fraction of a population that was not there.
# A cell with no finite value at all is `NA`, not `0`, on
# `.coherence_fraction()`'s discipline: a fraction that was not earned is
# withheld.
#
# An exactly zero value is counted in the denominator and not in the numerator.
# Under budget semantics that is the honest reading -- a group node no native
# mass reached received budget zero, and zero is a measurement -- but it does
# mean a node with poor coverage reports a low prevalence for a reason that is
# not about the effect. `$coverage$contributing` is the number that separates
# the two, and the print carries it beside the fraction rather than below it.
.population_prevalence_sign <- function(values, threshold) {
  finite <- is.finite(values)
  resolved <- apply(finite, c(1L, 2L), sum)
  above <- finite & values > threshold
  count <- apply(above, c(1L, 2L), sum)
  fraction <- count / resolved
  fraction[resolved == 0L] <- NA_real_
  contributing <- apply(finite & values != 0, c(1L, 2L), sum)
  storage.mode(resolved) <- "integer"
  storage.mode(count) <- "integer"
  storage.mode(contributing) <- "integer"
  list(fraction = fraction, count = count, resolved = resolved,
    contributing = contributing)
}

# Alignment prevalence, per group node.
#
# THE REFERENCE EXCLUDES THE PARTICIPANT BEING SCORED. With the plain group
# mean the inner product `<v_i, mean(v)>` contains the self term
# `||v_i||^2 / n`, which is positive whatever the participant does, so the
# fraction would be inflated by construction and would read exactly `1` at
# `n = 1`. The leave-one-out mean makes every product a *cross-participant*
# one, which is the same device `population-form-v1` section 6.2 uses to kill
# the consensus-squared term and section 7.1 uses to define `V^C`. Two
# participants are therefore the floor, and one is `NA` with the reason
# recorded rather than a `1`.
#
# A participant is scored only where its whole profile across the selected
# queries is finite: an inner product taken over the subset of coordinates that
# happened to resolve is an inner product in a different space for every
# participant.
.population_prevalence_alignment <- function(values) {
  nodes <- dim(values)[[1L]]
  fraction <- rep(NA_real_, nodes)
  count <- rep(NA_integer_, nodes)
  resolved <- integer(nodes)
  for (u in seq_len(nodes)) {
    profile <- matrix(values[u, , ], dim(values)[[2L]], dim(values)[[3L]])
    ok <- colSums(!is.finite(profile)) == 0L
    n <- sum(ok)
    resolved[[u]] <- n
    if (n < 2L) next
    kept <- profile[, ok, drop = FALSE]
    total <- rowSums(kept)
    # `<v_i, (S - v_i) / (n - 1)>`, written as one pass over the columns.
    products <- (colSums(kept * total) - colSums(kept^2)) / (n - 1L)
    count[[u]] <- sum(products > 0)
    fraction[[u]] <- count[[u]] / n
  }
  list(fraction = fraction, count = count, resolved = resolved)
}

# Section 7.5's coverage number, and section 14.3's open floor.
#
# What is counted here is what the shipped result can support: participants
# whose transported value at a node is finite and not exactly zero. Section
# 7.5's `group_node_subject_coverage` is a property of the transport operators
# themselves -- participants contributing nonzero *mass* -- and a population
# result does not carry them, so this is a proxy read off the values and is
# named as one. It agrees with the definition except where a participant's own
# native ledger vanishes at every native row the node draws from.
#
# The floor is an argument with no default, because section 14.3 records the
# threshold as an open maintainer decision: the contract requires the number to
# be reported and the marking mechanism to exist, and legislating a constant
# here would be this file deciding a study-design question on the reader's
# behalf.
.population_prevalence_coverage <- function(sign, nodes, floor) {
  minimum <- if (length(sign$contributing)) {
    apply(sign$contributing, 1L, min)
  } else {
    integer(0)
  }
  below <- if (is.null(floor)) character(0) else nodes[minimum < floor]
  list(
    definition = "participants_with_a_finite_nonzero_transported_value",
    proxy_for = "population-form-v1_section_7.5_group_node_subject_coverage",
    contributing = sign$contributing,
    minimum = as.integer(minimum),
    floor = if (is.null(floor)) NA_integer_ else as.integer(floor),
    below_floor = below
  )
}

# The record -------------------------------------------------------------------
#
# FIELD CONTRACT --- `effect_population_prevalence`
#
#   $layer      "latent_descriptive"; read this before any fraction
#   $reading    the one sentence every latent-layer record in the package
#               carries, shared with `latent_geometry()` so the two cannot
#               drift
#   $measures   the fractions this record carries, by name
#   $reference  0.5 -- the value a pure-noise cell reports, not 0
#   $sign       node x query: $fraction, $count, $resolved, $contributing
#   $alignment  per node: $fraction, $count, $resolved, plus the reference and
#               the readout inner product this bank defines
#   $coverage   the section 7.5 proxy, its minimum, the declared floor and the
#               nodes below it
#   $threshold  the comparison value, in ledger units, applied strictly
#   $index      the retained group nodes; the sink is never among them
#   $receipt    the standard population blocks, plus $prevalence and the
#               $form_prevalence_refusal that travels with every record
.population_prevalence_fields <- c(
  "layer", "reading", "measures", "reference", "sign", "alignment",
  "coverage", "threshold", "queries", "query_labels", "subjects", "index",
  "component", "ledger", "semantics", "normalization", "receipt",
  "scientific_plan_id"
)

# The parent is the *result's* identity and not the plan's, because a plan read
# through two components is two results and their prevalences are two different
# objects; a plan id would collide them under one signature. The threshold is in
# because it names a different estimand; the query selection is, because the
# alignment inner product is taken over exactly the selected coordinates; the
# coverage floor is not, because it decides what is *marked* rather than what
# was counted, which is the same line `heterogeneity()` draws at its node
# selection. The parent travels on the receipt so the validator can re-derive
# the identity rather than merely check its shape.
.population_prevalence_id <- function(parent, queries, threshold) {
  .sha256_signature(list(
    schema_version = 1L,
    contract_version = "population-form-v1",
    role = "population_prevalence",
    parent = parent,
    queries = queries,
    threshold = as.numeric(threshold),
    comparison = "strictly_greater",
    layer = "latent_descriptive",
    measures = .population_prevalence_measures
  ), "population-sha256:")
}

.population_prevalence_receipt <- function(base, threshold, queries, readout,
                                           coverage, subjects, nodes, parent) {
  # The standard blocks are carried unchanged, `$queries` included: on a
  # population receipt that field names the bank the *run* was read through,
  # and overwriting it with the subset this record counted over would change
  # what the field means rather than what it holds. The selection goes in the
  # prevalence block beside the threshold it belongs with.
  c(base, list(
    prevalence = list(
      layer = "latent_descriptive",
      parent = parent,
      queries = queries,
      measures = .population_prevalence_measures,
      comparison = "value > threshold",
      threshold = as.numeric(threshold),
      threshold_scale = "ledger_units",
      threshold_tolerance = "exact",
      # `conservative-geometry-v1` section 11.4 gap G3 is the open decision
      # about whether the package's nonnegativity guards take a relative
      # tolerance. This layer inherits it rather than making a private one,
      # the same way `latent_geometry()` does.
      threshold_tolerance_decision = "conservative-geometry-v1_11.4_G3",
      null_reference = 0.5,
      subjects = as.integer(subjects),
      retained_nodes = as.integer(nodes),
      sink_excluded = TRUE,
      sink_reason = "sink_is_budget_units_at_no_location",
      alignment_reference = "leave_one_out_participant_mean",
      alignment_inner_product = "query_readout_euclidean",
      alignment_frobenius_equivalent = readout$deviation <= 1e-10,
      readout_gram_deviation = readout$deviation,
      coverage_floor = coverage$floor,
      inference = "none_derivable"
    ),
    form_prevalence_refusal = .population_prevalence_form_refusal(base$basis)
  ))
}

.new_population_prevalence <- function(sign, alignment, coverage, readout,
                                       threshold, queries, weights, subjects,
                                       index, x, receipt, identity) {
  value <- structure(list(
    layer = "latent_descriptive",
    reading = .latent_reading_line,
    measures = .population_prevalence_measures,
    reference = 0.5,
    sign = sign,
    alignment = c(alignment, list(
      reference = "leave_one_out_participant_mean",
      inner_product = "query_readout_euclidean",
      readout_gram_deviation = readout$deviation,
      frobenius_equivalent = readout$deviation <= 1e-10
    )),
    coverage = coverage,
    threshold = as.numeric(threshold),
    queries = weights,
    query_labels = queries,
    subjects = subjects,
    index = index,
    component = x$component,
    ledger = x$ledger,
    semantics = x$semantics,
    normalization = x$normalization,
    receipt = receipt,
    scientific_plan_id = identity
  ), class = "effect_population_prevalence")
  .validate_population_prevalence(value)
  value
}

.validate_population_prevalence <- function(x) {
  if (!inherits(x, "effect_population_prevalence")) {
    .input_error(sprintf(paste0(
      "Expected an `effect_population_prevalence` from ",
      "`population_prevalence()`; received %s."
    ), .msg_value(x)), arg = "x", received = .msg_value(x),
      expected = paste0("an `effect_population_prevalence` from ",
        "`population_prevalence()`"))
  }
  if (!.sealed_fields(x, "effect_population_prevalence",
        .population_prevalence_fields) ||
      !identical(x$layer, "latent_descriptive") ||
      !identical(x$reading, .latent_reading_line) ||
      !identical(x$measures, .population_prevalence_measures) ||
      !identical(x$reference, 0.5) ||
      !is.list(x$sign) || !is.list(x$alignment) || !is.list(x$coverage) ||
      !.is_number(x$threshold) || !.is_strings(x$query_labels) ||
      !.is_strings(x$subjects, unique = TRUE) || !is.data.frame(x$index) ||
      !.is_string(x$component) || !.is_string(x$ledger) ||
      !.is_string(x$semantics) || !.is_string(x$normalization) ||
      !is.list(x$receipt) ||
      !.strong_sha256(sub("^population-", "", x$scientific_plan_id))) {
    .input_error("Population-prevalence fields are missing or noncanonical.")
  }
  nodes <- nrow(x$index)
  shape <- c(nodes, length(x$query_labels))
  if (!identical(dim(x$sign$fraction), shape) ||
      !identical(dim(x$sign$count), shape) ||
      !identical(dim(x$sign$resolved), shape) ||
      !identical(dim(x$coverage$contributing), shape) ||
      length(x$alignment$fraction) != nodes ||
      length(x$coverage$minimum) != nodes) {
    .contract_error(paste0(
      "A prevalence record's fractions must be indexed by the group nodes it ",
      "retained and the queries it counted over."
    ))
  }
  if (any(x$index$sink)) {
    .contract_error(paste0(
      "The sink is unmapped territory in budget units, not a location, and ",
      "it is never a row of a prevalence index."
    ))
  }
  # A fraction outside [0, 1], or one whose denominator is larger than the
  # participant list, is not a prevalence of anything. The counts are
  # re-divided rather than range-checked: `$fraction` is a function of
  # `$count` and `$resolved`, both of which the record publishes, so a
  # plausible fraction beside counts that do not produce it is detectable and
  # must be detected.
  for (block in list(x$sign, x$alignment)) {
    finite <- is.finite(block$fraction)
    if (any(block$count[finite] < 0L) ||
        any(block$count[finite] > block$resolved[finite]) ||
        any(block$resolved > length(x$subjects))) {
      .contract_error(paste0(
        "A prevalence count must lie between zero and the number of ",
        "participants resolved at the cell it counts."
      ))
    }
    derived <- block$count[finite] / block$resolved[finite]
    if (length(derived) &&
        max(abs(derived - block$fraction[finite])) > 1e-12) {
      .contract_error(paste0(
        "A prevalence fraction is not the one its own count and denominator ",
        "produce."
      ))
    }
    if (any(is.finite(block$fraction) &
        (block$fraction < 0 | block$fraction > 1))) {
      .contract_error("A prevalence fraction leaves [0, 1].")
    }
  }
  # The separation the whole layer rests on, asserted rather than documented:
  # nothing on this record may look like an error bar. A field named `se`,
  # `t`, `lower`, `upper` or `p_value` appearing anywhere in it would be a
  # reader's licence to test a count of participants against a null it does
  # not have.
  inferential <- c("se", "t", "lower", "upper", "p_value", "p", "statistic",
    "level", "confidence")
  present <- intersect(inferential,
    unlist(lapply(x[c("sign", "alignment", "coverage")], names)))
  if (length(present)) {
    .contract_error(sprintf(paste0(
      "A prevalence record carries %s, which names an inferential quantity. ",
      "This layer computes none: `population_uncertainty()` is the separate ",
      "layer that does."
    ), .msg_names(present)))
  }
  .check_signature(
    x$scientific_plan_id,
    .population_prevalence_id(x$receipt$prevalence$parent, x$query_labels,
      x$threshold),
    "Population-prevalence identity is inconsistent with its claims."
  )
  invisible(x)
}

# The verb ---------------------------------------------------------------------

#' Count how many participants a transported group value stands on
#'
#' `population_prevalence()` reports two fractions over the participants of an
#' [estimate_population()] result. **Sign prevalence** is, at each group node
#' and query, the fraction of participants whose transported ledger value
#' exceeds `threshold`. **Alignment prevalence** is, at each group node, the
#' fraction of participants whose whole profile across the query bank points
#' the same way as the rest of the group's.
#'
#' Both are descriptive, both live on the latent layer of
#' `design/conservative-geometry-contract.md` section 6 as
#' `design/population-form-contract.md` section 6.5 applies it, and neither
#' carries or admits an error bar. The record says so in `$layer`, in
#' `$reading`, and on every printed line.
#'
#' @section Why this is a latent-layer record:
#' The layer projects nothing --- there is no eigenvalue truncation here and no
#' moved mass to account for --- but `value > threshold` is the same kind of
#' operation and fails in the same direction. It is a per-participant sign
#' clamp: it discards the magnitude of a crossvalidated estimate and keeps the
#' sign, which is the part that carries the noise. **A group node and query at
#' which nothing reproduces reports a fraction near `0.5`, not near `0`**,
#' because every participant contributes an independent coin flip. `$reference`
#' carries that number so the comparison is not made against zero by habit.
#'
#' The threshold is applied strictly (`>`) and absolutely, in the ledger's own
#' units, with no relative tolerance. That is the guard every per-node fraction
#' in the package uses, and whether these guards should take a relative
#' tolerance is an open contract decision
#' (`conservative-geometry-contract.md` section 11.4, gap G3) that this layer
#' inherits rather than settling privately. A cell that is zero only to within
#' rounding is counted by its rounding.
#'
#' @section There is no inference here, and none is derivable:
#' A count of participants is not an estimate of a population quantity on this
#' layer. The participants are not a sample from a declared superpopulation of
#' transports, and the binomial variance of a thresholded crossvalidated
#' estimate is not the sampling variance of anything the plan sealed. This
#' function computes no standard error, no interval and no p-value, the record
#' refuses to carry a field named like one, and constructing a prevalence
#' leaves the result's own `$uncertainty` untouched.
#' [population_uncertainty()] is the separate inferential layer, and it carries
#' its own measured calibration record.
#'
#' @section What alignment is an inner product in:
#' Participant `i`'s profile at a group node is `v_i = H y_i`, with `H` the
#' packed queries the driver lowers the bank to and `y_i` the participant's
#' packed transported form. The scored quantity is
#' \eqn{\langle v_i, \bar v_{-i}\rangle}{<v_i, mean of the others>} against the
#' **leave-one-out** mean of the other participants. The reference excludes the
#' participant being scored on purpose: the plain group mean puts the positive
#' self term \eqn{\lVert v_i\rVert^2/n}{||v_i||^2 / n} inside every product, so
#' the fraction would be inflated by construction and would read `1` at
#' `n = 1`. Excluding it makes every product a cross-participant one, which is
#' the same device `population-form-v1` sections 6.2 and 7.1 use.
#'
#' That inner product is the Frobenius inner product of the transported forms
#' *restricted to the bank's span* exactly when the packed queries are
#' orthonormal, and a Euclidean inner product in coordinates the bank chose
#' otherwise. `$alignment$readout_gram_deviation` measures the departure and
#' `$alignment$frobenius_equivalent` states the verdict; a two-condition
#' contrast bank is not orthonormal.
#'
#' @section The sink, and coverage:
#' The sink is excluded, on the same grounds [heterogeneity()] excludes it: it
#' is unmapped territory reported in budget units at no location, so "the
#' fraction of participants whose unmapped mass is positive" is a statement
#' about transport failure and not about a place. The exclusion is recorded at
#' `$receipt$prevalence$sink_excluded` and the per-participant sink budget
#' stays readable at `$receipt$sink_budget`.
#'
#' A prevalence is uninterpretable without the number of participants behind
#' it. `$coverage$contributing` counts, per node and query, the participants
#' whose transported value is finite and not exactly zero --- a proxy, read off
#' the shipped values, for `population-form-v1` section 7.5's
#' `group_node_subject_coverage`, which is a property of the transport
#' operators a result does not carry. `coverage_floor` marks the nodes below a
#' declared floor; it has no default because section 14.3 records the threshold
#' itself as an open maintainer decision.
#'
#' @section Refusals:
#' Each is an `effect_capability_refusal` in namespace
#' `"population_prevalence"` (see [catch_refusal()]).
#'
#' * `population_participant_values` --- a `materialize_population()` result.
#'   The streamed complete-form route retains no participant axis, by design.
#' * `participant_term_decomposition` --- a non-`NULL` `term`. Participants
#'   carry one transported ledger, not one per group-model column.
#' * `comparable_ledger_orientation` --- a `"unit_budget"` population in which
#'   some participant's signed native total is negative, so its values are
#'   sign-flipped relative to the others and a sign count would mix two
#'   questions.
#'
#' One refusal is carried rather than raised. The form-space reading --- the
#' fraction of participants whose transported *form* has positive Frobenius
#' inner product with the group mean form --- needs per-participant transported
#' forms, and no population result ships them: the query-bank route retains
#' each participant's value only in the coordinates the declared bank names,
#' and the streamed route retains no participant axis. Capability
#' `participant_form_prevalence` travels on
#' `$receipt$form_prevalence_refusal` with its reasons and remedies.
#'
#' @param x An `effect_population_result` from [estimate_population()].
#' @param query Which queries to count over, by name or position. `NULL` (the
#'   default) counts over every query in the bank. The alignment inner product
#'   is taken over exactly the selected coordinates.
#' @param threshold The value a participant's transported ledger must strictly
#'   exceed to be counted, in the ledger's own units. `0` is the default and
#'   counts participants with a positive contribution.
#' @param coverage_floor Optional whole number. Group nodes whose smallest
#'   contributing-participant count falls below it are listed in
#'   `$coverage$below_floor`. `NULL` (the default) reports the counts and marks
#'   nothing.
#' @param term Refused. Present so that a caller arriving from
#'   [population_uncertainty()] gets a reason rather than silence.
#' @return An `effect_population_prevalence`: `$sign` holding `$fraction`,
#'   `$count`, `$resolved` and `$contributing` (each a `node`-by-`query`
#'   matrix); `$alignment` holding `$fraction`, `$count`, `$resolved` and the
#'   inner product it was taken in, one entry per group node; `$coverage`;
#'   `$reference`, the fraction a pure-noise cell reports; `$layer` and
#'   `$reading`, which mark the record as descriptive; `$queries` and
#'   `$query_labels`, the selected rows of the bank it counted over; and the
#'   `$index`, `$ledger`, `$semantics`, `$normalization` and `$receipt` of the
#'   result it read. `as.data.frame(x, measure = )` returns one measure in long
#'   form.
#' @references `design/population-form-contract.md` (`population-form-v1`),
#'   sections 6.5, 7.5 and 8.1;
#'   `design/conservative-geometry-contract.md` (`conservative-geometry-v1`),
#'   section 6.
#' @family population transports
#' @seealso [latent_geometry()] for the projection layer this one shares its
#'   discipline with, [population_uncertainty()] for the separate inferential
#'   layer, and [heterogeneity()] for the form-space reading of how
#'   participants differ.
#' @examples
#' # Six participants on different native frames, carrying one planted effect
#' # every run of every participant shares (face over house) and a phase-shifted
#' # pattern that reproduces in none of them. Two group nodes: one every
#' # participant reaches and one only the larger native frames do.
#' effects <- effect_space(c("face", "house", "tool"), basis_id = "popprev:v1")
#' subject <- function(id, n, phase) {
#'   domain <- abstract_domain(n, coordinates = cbind(x = seq_len(n) - 1),
#'     feature_ids = paste0("f", seq_len(n)), id = id)
#'   planted <- outer(c(face = 0.8, house = -0.8, tool = 0), rep(1, n))
#'   values <- function(shift) matrix(cos(seq_len(3 * n) + shift), 3, n,
#'     dimnames = list(c("face", "house", "tool"), NULL)) + planted
#'   rel <- relation(list(run1 = values(0), run2 = values(phase)),
#'     effects = effects, domain = domain)
#'   plan_geometry(rel, compile_frame(voxelwise(), domain),
#'     cross_partitions(rel))
#' }
#' carrier <- function(n) anatomical_transport(
#'   native_coords = cbind(seq_len(n) - 1), group_coords = cbind(c(0, 6)),
#'   semantics = "budget", radius = 2
#' )
#' sizes <- c(s01 = 4L, s02 = 5L, s03 = 6L, s04 = 7L, s05 = 8L, s06 = 9L)
#' phases <- c(s01 = 0.2, s02 = 1.1, s03 = 2.4, s04 = 0.6, s05 = 3, s06 = 1.7)
#' subjects <- stats::setNames(lapply(names(sizes), function(id)
#'   subject(id, sizes[[id]], phases[[id]])), names(sizes))
#' plan <- plan_population(subjects, lapply(sizes, carrier))
#' fit <- estimate_population(plan,
#'   rbind(`face-house` = c(1, -1, 0), `face-tool` = c(1, 0, -1)))
#'
#' shared <- population_prevalence(fit)
#' shared
#'
#' # The planted contrast is carried by every participant; the one that
#' # reproduces in none of them sits at the 0.5 the print warns about.
#' shared$sign$fraction
#'
#' # And how many participants were there to be counted at each node.
#' shared$coverage$contributing
#'
#' # A term axis does not exist on participant values, and asking says why.
#' catch_refusal(population_prevalence(fit, term = "(Intercept)"))$capability
#' @export
population_prevalence <- function(x, query = NULL, threshold = 0,
                                  coverage_floor = NULL, term = NULL) {
  .validate_population_result(x)
  if (!is.null(term)) .population_prevalence_refuse_term()
  if (!identical(x$basis, "query_bank")) {
    .population_prevalence_refuse_basis()
  }
  threshold <- .check_number(threshold, "threshold",
    what = "one finite value in the ledger's own units")
  if (!is.null(coverage_floor)) {
    coverage_floor <- .check_count(coverage_floor, "coverage_floor",
      max = dim(x$values)[[3L]],
      what = paste0("one whole number of participants no larger than the ",
        dim(x$values)[[3L]], " the population carries"))
  }
  queries <- .population_prevalence_query(x, query)

  # Section 4.3's signed divisor, checked before anything is counted: a
  # sign-flipped participant makes the count itself meaningless, so this is a
  # gate and not a mark.
  if (identical(x$normalization, "unit_budget")) {
    totals <- x$receipt$native_total[, queries, drop = FALSE]
    admitted <- x$receipt$normalization$admitted[, queries, drop = FALSE]
    flipped <- rownames(totals)[rowSums(admitted & totals < 0) > 0L]
    if (length(flipped)) .population_prevalence_refuse_orientation(flipped)
  }

  keep <- !x$index$sink
  index <- x$index[keep, , drop = FALSE]
  rownames(index) <- NULL
  values <- x$values[keep, queries, , drop = FALSE]
  node_labels <- as.character(index$node)
  subjects <- dimnames(x$values)[[3L]]

  sign <- .population_prevalence_sign(values, threshold)
  alignment <- .population_prevalence_alignment(values)
  for (field in c("fraction", "count", "resolved", "contributing")) {
    dimnames(sign[[field]]) <- list(node = node_labels, query = queries)
  }
  for (field in c("fraction", "count", "resolved")) {
    names(alignment[[field]]) <- node_labels
  }
  weights <- x$queries[queries, , drop = FALSE]
  readout <- .population_prevalence_readout_gram(weights)
  coverage <- .population_prevalence_coverage(sign, node_labels,
    coverage_floor)
  names(coverage$minimum) <- node_labels

  identity <- .population_prevalence_id(x$scientific_plan_id, queries,
    threshold)
  receipt <- .population_prevalence_receipt(x$receipt, threshold, queries,
    readout, coverage, length(subjects), nrow(index), x$scientific_plan_id)
  .new_population_prevalence(sign, alignment, coverage, readout, threshold,
    queries, weights, subjects, index, x, receipt, identity)
}

# Presentation -----------------------------------------------------------------

#' @export
as.data.frame.effect_population_prevalence <- function(
    x, row.names = NULL, optional = FALSE, ...,
    measure = c("sign", "alignment")) {
  .validate_population_prevalence(x)
  measure <- match.arg(measure)
  # One measure at a time, on `population_uncertainty()`'s precedent and for
  # the same structural reason: sign prevalence is indexed by node and query
  # and alignment prevalence has no query axis at all, so a single table
  # holding both would have to repeat one of them down the other's rows.
  if (identical(measure, "alignment")) {
    return(data.frame(
      node = names(x$alignment$fraction),
      prevalence = unname(x$alignment$fraction),
      participants = unname(x$alignment$count),
      resolved = unname(x$alignment$resolved),
      stringsAsFactors = FALSE
    ))
  }
  nodes <- rownames(x$sign$fraction)
  data.frame(
    node = rep(nodes, times = length(x$query_labels)),
    query = rep(x$query_labels, each = length(nodes)),
    prevalence = as.numeric(x$sign$fraction),
    participants = as.integer(x$sign$count),
    resolved = as.integer(x$sign$resolved),
    contributing = as.integer(x$sign$contributing),
    stringsAsFactors = FALSE
  )
}

.pf_prevalence_spread <- function(fraction) {
  finite <- fraction[is.finite(fraction)]
  if (!length(finite)) {
    return("undefined; no cell resolved a participant")
  }
  masked <- length(fraction) - length(finite)
  paste0(
    sprintf("median %s (range %s to %s)", .pf_num(stats::median(finite), 3L),
      .pf_num(min(finite), 3L), .pf_num(max(finite), 3L)),
    if (masked) sprintf(", %s masked", .msg_count(masked, "cell")) else ""
  )
}

.pf_prevalence_coverage <- function(x) {
  coverage <- x$coverage
  floor <- coverage$floor
  if (!length(coverage$minimum)) {
    return("no group node survived the sink exclusion")
  }
  paste0("minimum ", min(coverage$minimum), " of ",
    .msg_count(length(x$subjects), "participant"), " contributing",
    if (is.na(floor)) {
      "; no floor declared"
    } else if (length(coverage$below_floor)) {
      paste0("; ", .msg_count(length(coverage$below_floor), "node"),
        " below floor ", floor, " (", .pf_set(coverage$below_floor, max = 3L),
        ")")
    } else {
      paste0("; every node meets floor ", floor)
    })
}

.pf_prevalence_readout <- function(x) {
  paste0(.msg_count(length(x$query_labels), "query", "queries"), "; ",
    if (isTRUE(x$alignment$frobenius_equivalent)) {
      "packed queries orthonormal, Frobenius on the span"
    } else {
      paste0("bank Gram off identity by ",
        .pf_num(x$alignment$readout_gram_deviation, 3L), ", not Frobenius")
    })
}

#' @export
format.effect_population_prevalence <- function(x, ...) {
  .pf_inline("effect_population_prevalence", x$ledger,
    .msg_count(length(x$subjects), "participant"),
    paste0(nrow(x$index), " group nodes"),
    paste0("threshold ", .pf_num(x$threshold, 3L)))
}

#' @export
print.effect_population_prevalence <- function(x, ...) {
  .validate_population_prevalence(x)
  .pf_emit("effect_population_prevalence", list(
    layer = "latent descriptive (fractions over participants)",
    ledger = .pf_population_ledger(x),
    participants = paste0(length(x$subjects), " (",
      .pf_set(x$subjects, max = 4L), ")"),
    `group nodes` = paste0(nrow(x$index), " (sink excluded)"),
    queries = .pf_set(x$query_labels, max = 4L),
    threshold = paste0("value > ", .pf_num(x$threshold, 3L),
      " (ledger units, exact)"),
    sign = .pf_prevalence_spread(x$sign$fraction),
    alignment = paste0(.pf_prevalence_spread(x$alignment$fraction),
      ", leave-one-out reference"),
    readout = .pf_prevalence_readout(x),
    coverage = .pf_prevalence_coverage(x),
    frame = .pf_population_frame(x),
    transport = .pf_population_carrier(x),
    estimand = .pf_sig(x$scientific_plan_id),
    reading = x$reading
  ))
  .pf_note(paste0(
    "a cell at which nothing reproduces reports a fraction near ", x$reference,
    ", not near 0: thresholding a signed crossvalidated estimate keeps the ",
    "sign and discards the magnitude, and the sign is the noisy part."
  ))
  .pf_note(paste0(
    "descriptive only. No standard error, interval or p-value is attached to ",
    "a count of participants, and none follows from one; ",
    "population_uncertainty() is the separate inferential layer."
  ))
  if (!identical(x$component, "total")) {
    .pf_note(paste0(
      "native_coherent_ledger + native_configuration_ledger = ",
      "transported_total holds exactly; this is native-node ", x$component,
      " evidence carried to a group location, not a group-node common mode."
    ))
  }
  cat(sprintf("  %-14s%s\n", "next:",
    "as.data.frame(x, measure), x$coverage, x$sign$resolved"))
  invisible(x)
}
