# Population plans -----------------------------------------------------------
#
# Layer 3 (plans). A population plan is the group-level *estimand*: N
# participants' compiled conservative geometry, one declared transport each
# onto one shared group node set, and one group model over the subject axis.
# It reads no data, opens no source, and fits nothing. What it does is name
# what a group number will mean, seal that name into a content-addressed
# identity, and refuse every combination that would make the name false.
#
# `design/population-form-contract.md` (`population-form-v1`) is normative
# throughout:
#
#   * the transport is part of the estimand, not a preprocessing step, so it
#     enters the scientific identity rather than only the physical signature
#     (§1.5);
#   * OLS is the default group fit because subject-constant weights are the
#     commuting class --- per-node weights break transport-then-fit and
#     per-coordinate weights break query-then-fit, both by ~25 % on the
#     contract's own fixtures (§3.1--§3.3);
#   * with heterogeneous native frames there is no common node axis, so
#     fit-then-transport is not even definable and the evaluation order is
#     forced rather than chosen (§3, "scope of the tensor argument");
#   * the budget normalization is a closed set of three which differ in
#     argmax, so it is a plan-identity field with `none` the only default
#     (§4).
#
# Nothing here computes a group number; that is the executor's job. The file
# calls sideways to `geometry-plan.R` for the subject plan validator and
# downward to `transport.R`, `frame.R`, `compute-policy.R` and the guards.

.population_normalizations <- c("none", "unit_budget", "precision_weighted")

# Every per-row frame metadata column `population-form-v1` §9.1 requires of a
# transported native index. Recorded per subject rather than refused on:
# §9.1's readiness table has eleven rows and two of them (`composition_and_root`,
# `latent_projection_named`) are supplied by nothing in the tree today, so a
# constructor that refused the whole table would refuse every input that can
# currently exist. What the plan can honestly do is carry the readiness of the
# rows that *are* supplied into the printed record.
.population_row_metadata <- c("family", "scale", "center", "alpha")

# `precision_weighted` requires `pi_i = 1 / Var(T_i-hat)`: the variance of a
# *conserved budget*, which `population-form-v1` §4.5 shows needs the full
# cross-node sampling covariance and not a sum of per-node margins. That is
# the object `conservative-geometry-v1` gap G8 records as missing and ticket
# D8 owes.
#
# The gate is decided statically against the installed tree rather than
# probed at runtime, so that a plan's admissible normalizations do not depend
# on which session sealed it. Checked at this file's authoring commit: D8's
# query-bank route (`sampling_covariance(x, queries =)`) has landed, and it is
# *not* the missing object. Its own contract says so --- each covariance is
# local to one measurement, and `scope = "cross_measurement"` refuses with
# capability `cross_node_sampling_covariance`, reasons
# `spatial_covariance_model_unavailable` and `conservation_gives_no_uncertainty`
# (`.require_local_sampling_scope()`). Per-row variances therefore do not add
# to the variance of a conserved budget. The route that would be needed is
# exactly the one D8 refuses, so the gate stands rather than being lifted by
# the arrival of a `queries =` argument that answers a different question.
#
# D8's own second remedy points here: "supply cross-node precision externally;
# the population layer declares that as its own input rather than deriving it
# here." That is §4.5's interim route, and it is a `precision =` argument this
# constructor does not yet have.
#
# Lifting it is `population-form-v1` §14.1's open maintainer decision and
# needs two things together: a per-subject cross-node covariance (an explicit
# spatial model, being the cross-support residual second moment and the frame
# overlap it induces), and a `precision =` argument here carrying the
# externally supplied, provenance-bearing vector §4.5 permits in the interim.
.population_precision_refusal <- function() {
  .capability_refusal(paste0(
    "`normalization = \"precision_weighted\"` needs a per-subject precision ",
    "1 / Var(T_i), and the variance of a conserved budget requires the full ",
    "cross-node sampling covariance rather than a sum of per-node margins. ",
    "crossform has no route that produces one -- `sampling_covariance()` ",
    "reads a query bank at one measurement and refuses ",
    "`scope = \"cross_measurement\"` for want of a spatial model -- and it ",
    "refuses to synthesize one. Conservation is a law about estimates, not ",
    "about their uncertainty, so a conserved budget buys no error bars and ",
    "per-measurement blocks must not be assembled into a joint covariance: a ",
    "precision built that way would silently reweight subjects by the wrong ",
    "quantity, and the contract measures the three normalizations disagreeing ",
    "by up to 94 % and moving the argmax."
  ),
    capability = "precision_weighted_normalization",
    namespace = "population_plans",
    reasons = c(
      "per_subject_budget_variance_unavailable",
      "cross_node_sampling_covariance_route_absent"
    ),
    remedies = c(
      paste0(
        "Use `normalization = \"none\"` for the mean subject ledger in ",
        "native evidence units, or `\"unit_budget\"` for the mean ",
        "attribution share."
      ),
      paste0(
        "The gate lifts when a cross-node sampling covariance exists (gap G8; ",
        "`design/population-form-contract.md` sections 4.5 and 14.1). D8's ",
        "query bank is not it: it refuses the cross-measurement scope that ",
        "route would need."
      )
    ))
}

# Subjects, canonically ordered ------------------------------------------------
#
# The list is sorted by subject name with a locale-independent radix order, so
# that two sessions that listed the same participants in different orders seal
# the same estimand. Sorting is what lets the identity quote the subject ids
# directly: because the order is canonical, hashing the *named* vector both
# makes the identity order-invariant and pins which transport belongs to which
# participant, which a bare sorted set of ids would not.
.population_subjects <- function(subjects) {
  if (!is.list(subjects) || inherits(subjects, "effect_geometry_plan")) {
    .input_error(paste0(
      "`subjects` must be a named list of compiled geometry plans, one per ",
      "participant."
    ), arg = "subjects", received = .msg_value(subjects),
      expected = "a named list of `effect_geometry_plan` values")
  }
  labels <- names(subjects)
  if (length(subjects) < 2L) {
    .input_error(paste0(
      "A population plan needs at least two participants; received ",
      .msg_count(length(subjects), "subject"), ". Every population quantity ",
      "the contract defines -- the cross-participant consensus term, the ",
      "subject Gram, the group fit -- is a contrast across participants, and ",
      "there is none to take with one."
    ), arg = "subjects", received = .msg_count(length(subjects), "subject"),
      expected = "at least two subjects")
  }
  if (!.is_strings(labels, unique = TRUE) ||
      length(labels) != length(subjects)) {
    .input_error(paste0(
      "`subjects` must be named with one nonempty, unique participant ",
      "identifier each. The names are what bind a transport and a design row ",
      "to a participant; positions alone cannot survive a group result."
    ), arg = "subjects", expected = "unique nonempty names")
  }
  wrong <- labels[!vapply(subjects, inherits, logical(1),
    "effect_geometry_plan")]
  if (length(wrong)) {
    .input_error(sprintf(paste0(
      "Every entry of `subjects` must be an `effect_geometry_plan` from ",
      "`plan_geometry()`; %s %s not."
    ), .msg_names(wrong), if (length(wrong) == 1L) "is" else "are"),
      arg = "subjects",
      expected = "compiled geometry plans from `plan_geometry()`")
  }
  subjects[order(labels, method = "radix")]
}

# The three properties of a subject plan the population layer needs before it
# can call the transported number a group quantity. Each is collected across
# every subject and refused as one condition naming all of them, because a
# study that miscompiled one plan usually miscompiled several and reporting
# only the first turns one fix into six rounds.
.population_admit_subject_plans <- function(subjects) {
  labels <- names(subjects)
  for (plan in subjects) .validate_geometry_plan(plan, deep = FALSE)

  rectangular <- labels[vapply(subjects, function(plan)
    !identical(plan$codec, "symmetric_packed"), logical(1))]
  if (length(rectangular)) {
    .capability_refusal(paste0(
      "The population layer stores transported geometry in the ",
      "`symmetric_packed` codec, whose root-two off-diagonals are what make ",
      "the packed inner product equal the Frobenius inner product -- the ",
      "subject-Gram trick and every population spectrum read off it are ",
      "statements about that isometry. A rectangular cross-axis plan has no ",
      "symmetric form to pack, so there is nothing here for a group node to ",
      "carry."
    ),
      capability = "symmetric_packed_population",
      namespace = "population_plans",
      reasons = paste0("rectangular_subject_plan:", rectangular),
      remedies = paste0(
        "Supply self-form geometry plans (`plan_geometry(x, at, over)` with ",
        "no `right =`) for every participant."
      ))
  }

  learned <- labels[vapply(subjects, function(plan)
    identical(plan$metric_schedule$kind, "learned_local_before_frame"),
    logical(1))]
  if (length(learned)) {
    .capability_refusal(paste0(
      "Population plans do not yet admit subject plans on a learned neural ",
      "metric. A locally learned metric is estimated inside each ",
      "participant's own frame, so the transported ledgers of two subjects ",
      "are budgets under two different, data-dependent metrics and their ",
      "difference at a group node is not an estimate of anything. Carrying ",
      "it honestly needs a declared cross-participant metric composition, ",
      "which no contract fixes today."
    ),
      capability = "fixed_metric_population",
      namespace = "population_plans",
      reasons = paste0("learned_metric_subject_plan:", learned),
      remedies = paste0(
        "Plan each participant's geometry on the implicit identity metric or ",
        "on one declared fixed `neural_metric()`, and transport that."
      ))
  }

  spaces <- vapply(subjects, function(plan)
    plan$task$left_relation$effect_space$signature, character(1))
  reference <- spaces[[1L]]
  disagree <- labels[spaces != reference]
  if (length(disagree)) {
    .capability_refusal(paste0(
      "Every participant's geometry must be measured over one common ",
      "experimental space. Transport carries *nodes*, never coordinates: it ",
      "cannot align two participants whose forms are indexed by different ",
      "effects, different orders, different bases, or different units, and a ",
      "group form assembled from them would add entries that name different ",
      "quantities."
    ),
      capability = "common_experimental_space",
      namespace = "population_plans",
      reasons = c(
        paste0("reference_effect_space:", reference),
        paste0("disagreeing_subject:", disagree, ":", spaces[disagree])
      ),
      remedies = paste0(
        "Declare one `effect_space()` -- the same coordinates in the same ",
        "order, basis, units and scale -- and build every participant's ",
        "relation against it."
      ))
  }
  invisible(subjects)
}

# `population-form-v1` consumes a *conservative* result: the node-sum law of
# `conservative-geometry-v1` section 2 is what makes a native ledger a budget,
# and budget preservation under transport (section 2 here) is a statement about
# that budget. A locally normalized frame double-counts overlapping
# neighbourhoods, so its node values are not contributions to any whole-brain
# quantity and transporting them produces a number with no conserved total
# behind it.
#
# `declared_normalization` rather than `normalization` is the right field to
# read: a frame that has had a diagonal metric folded into its weights reports
# `normalization = "none"` while still being conservative under that metric,
# and `frame_conservation()` already resolves the two.
.population_admit_conservation <- function(subjects, conservation,
                                           allow_nonconservative) {
  labels <- names(subjects)
  declared <- vapply(conservation, function(report)
    report$declared_normalization, character(1))
  conserved <- vapply(conservation, function(report)
    isTRUE(report$conserved), logical(1))
  not_declared <- labels[declared != "conservative"]
  not_certified <- labels[declared == "conservative" & !conserved]
  if (!length(not_declared) && !length(not_certified)) {
    return(invisible(NULL))
  }
  if (isTRUE(allow_nonconservative)) {
    return(invisible(NULL))
  }
  .capability_refusal(paste0(
    "A population plan carries conservative subject geometry: the native ",
    "ledger has to be a budget before transporting it can preserve one. A ",
    "frame that is not conservative has no conserved total behind its node ",
    "values -- under `\"local\"` normalization overlapping neighbourhoods ",
    "double-count shared features -- so the transported group number is not ",
    "a share of anything, and the sink column stops being an accounting of ",
    "what went missing."
  ),
    capability = "conservative_subject_geometry",
    namespace = "population_plans",
    reasons = c(
      if (length(not_declared)) {
        paste0("normalization_not_conservative:", not_declared, ":",
          declared[not_declared])
      },
      if (length(not_certified)) {
        paste0("conservation_certificate_failed:", not_certified)
      }
    ),
    remedies = c(
      paste0(
        "Compile each participant's frame with ",
        "`normalization = \"conservative\"` -- the default for ",
        "`voxelwise()`, `searchlights()` and `regions()` -- or fold a ",
        "diagonal metric into a conservative frame."
      ),
      paste0(
        "Pass `allow_nonconservative = TRUE` to plan it anyway. The flag is ",
        "recorded on the plan and enters its scientific identity, because a ",
        "non-conservative population estimand is a different estimand and ",
        "not a relaxed setting on the same one."
      )
    ))
}

# Transports ------------------------------------------------------------------

.population_transports <- function(transport, subjects) {
  labels <- names(subjects)
  if (!is.list(transport) || inherits(transport, "effect_location_transport")) {
    .input_error(paste0(
      "`transport` must be a named list of location transports, one per ",
      "participant."
    ), arg = "transport", received = .msg_value(transport),
      expected = "a named list of `effect_location_transport` values")
  }
  given <- names(transport)
  if (!.is_strings(given, unique = TRUE) || length(given) != length(transport)) {
    .input_error(paste0(
      "`transport` must be named with one nonempty, unique participant ",
      "identifier each, matching the names of `subjects`."
    ), arg = "transport", expected = "unique nonempty names")
  }
  wrong <- given[!vapply(transport, inherits, logical(1),
    "effect_location_transport")]
  if (length(wrong)) {
    .input_error(sprintf(paste0(
      "Every entry of `transport` must be an `effect_location_transport`; ",
      "%s %s not."
    ), .msg_names(wrong), if (length(wrong) == 1L) "is" else "are"),
      arg = "transport",
      expected = "transports from `location_transport()` or a named rule")
  }
  unmatched <- setdiff(labels, given)
  extra <- setdiff(given, labels)
  if (length(unmatched) || length(extra)) {
    .capability_refusal(paste0(
      "Every participant needs exactly one transport and every transport ",
      "needs a participant. The names are the binding: a transport is an ",
      "operator on one participant's own native frame, so pairing it with ",
      "another participant's geometry by position would carry evidence ",
      "between people without anything in the result saying so."
    ),
      capability = "aligned_subject_transports",
      namespace = "population_plans",
      reasons = c(
        if (length(unmatched)) paste0("subject_without_transport:", unmatched),
        if (length(extra)) paste0("transport_without_subject:", extra)
      ),
      remedies = paste0(
        "Name `subjects` and `transport` with the same participant ",
        "identifiers: ", .msg_names(labels), "."
      ))
  }
  transport <- transport[labels]
  for (value in transport) .validate_location_transport(value)

  native <- vapply(transport, function(value) nrow(value$matrix), integer(1))
  measurements <- vapply(subjects, function(plan) as.integer(plan$measurements),
    integer(1))
  mismatched <- labels[native != measurements]
  if (length(mismatched)) {
    .capability_refusal(paste0(
      "A transport's native rows are the participant's own conservative ",
      "frame nodes, in that frame's order, so their count has to equal the ",
      "number of measurements the geometry plan produces. A transport of the ",
      "wrong height is an operator on some other frame, and applying it ",
      "would relabel evidence rather than carry it."
    ),
      capability = "aligned_subject_transports",
      namespace = "population_plans",
      reasons = paste0("native_size_mismatch:", mismatched, ":",
        native[mismatched], "_rows_for_", measurements[mismatched],
        "_measurements"),
      remedies = paste0(
        "Build each transport against the frame its participant's plan was ",
        "compiled at, for example from `plan$frame$index$measurement`."
      ))
  }
  transport
}

# All transports must land on one group node set, with the same semantics.
# The group index carries the labels and coordinates every coverage and
# displacement diagnostic is defined against, so two subjects reaching
# different node sets are not two views of one population; and budget and
# density are two different estimands built from the same operator, so a
# population mixing them would average a signed sum against a ratio.
.population_group_nodes <- function(transport) {
  labels <- names(transport)
  reference <- transport[[1L]]$group_index
  disagree <- labels[vapply(transport, function(value)
    !identical(value$group_index, reference), logical(1))]
  mixed <- vapply(transport, function(value) value$semantics, character(1))
  disagreeing_semantics <- labels[mixed != mixed[[1L]]]
  if (length(disagree) || length(disagreeing_semantics)) {
    .capability_refusal(paste0(
      "Every transport must land on one shared group node set under one ",
      "semantics. The group index is what a coverage count, a displacement ",
      "and a group estimate are all indexed by; two participants reaching ",
      "different nodes have nowhere to meet, and budget and density are two ",
      "different estimands built from the same operator, so a plan mixing ",
      "them would combine a signed sum with a ratio."
    ),
      capability = "shared_group_nodes",
      namespace = "population_plans",
      reasons = c(
        paste0("reference_subject:", labels[[1L]]),
        if (length(disagree)) paste0("group_index_disagrees:", disagree),
        if (length(disagreeing_semantics)) {
          paste0("semantics_disagrees:", disagreeing_semantics, ":",
            mixed[disagreeing_semantics])
        }
      ),
      remedies = paste0(
        "Build every participant's transport against one `group_index` -- ",
        "the same node identifiers and coordinates in the same order -- and ",
        "declare the same `semantics` for all of them."
      ))
  }
  reference
}

# The group model --------------------------------------------------------------

# One row per participant, and the row has to say which participant it is.
# Trusting row order here is the same mistake as pairing transports by
# position: it would bind a covariate to the wrong person and leave nothing in
# the result that could show it.
.population_model_data <- function(data, labels) {
  if (is.null(data)) {
    return(data.frame(row.names = labels))
  }
  if (!is.data.frame(data)) {
    .input_error(paste0(
      "`data` must be a data frame with one row per participant, or NULL for ",
      "an intercept-only group model."
    ), arg = "data", received = .msg_value(data),
      expected = "a data frame with one row per subject")
  }
  data <- as.data.frame(data, stringsAsFactors = FALSE)
  automatic <- identical(rownames(data), as.character(seq_len(nrow(data))))
  key <- if ("subject" %in% names(data)) {
    as.character(data$subject)
  } else if (!automatic) {
    rownames(data)
  } else {
    NULL
  }
  if (is.null(key) || anyNA(key) || anyDuplicated(key) ||
      !setequal(key, labels)) {
    .capability_refusal(paste0(
      "The group design needs one row per participant, and the row has to ",
      "name the participant it belongs to. crossform will not bind ",
      "covariates by row order: a silently misaligned design row attaches a ",
      "covariate to the wrong person and leaves no trace of it in the group ",
      "estimate."
    ),
      capability = "aligned_subject_rows",
      namespace = "population_plans",
      reasons = c(
        paste0("subjects:", .msg_names(labels)),
        paste0("data_rows:", if (is.null(key)) "unlabelled" else
          .msg_names(key))
      ),
      remedies = paste0(
        "Give `data` a `subject` column, or row names, holding exactly the ",
        "names of `subjects` -- one row each, in any order."
      ))
  }
  data <- data[match(labels, key), , drop = FALSE]
  rownames(data) <- labels
  data
}

# The model matrix is built once, here, and its pivoted QR is stored with it.
# The QR is the group fit's whole numerical content under OLS: the same
# factorization serves every group node and every packed coordinate, because
# the design does not vary along either axis. Storing it on the plan is what
# makes that a fact about the estimand rather than a caching decision made
# inside a loop -- and the pivot is what lets a rank failure name the
# offending columns instead of reporting a number.
.population_model <- function(model, data, labels) {
  if (!inherits(model, "formula")) {
    .input_error(paste0(
      "`model` must be a one-sided group-level formula, evaluated against ",
      "one row per participant; the default `~ 1` is the group mean."
    ), arg = "model", received = .msg_value(model),
      expected = "a one-sided formula such as `~ 1` or `~ group + age`")
  }
  data <- .population_model_data(data, labels)
  terms <- stats::terms(model, data = data)
  if (!identical(as.integer(attr(terms, "response")), 0L)) {
    .capability_refusal(paste0(
      "`model` must be one-sided. The response of a population fit is the ",
      "transported geometry itself -- one column per group node and packed ",
      "coordinate -- so a left-hand side would name a second, conflicting ",
      "response."
    ),
      capability = "identified_group_model",
      namespace = "population_plans",
      reasons = "model_formula_has_response",
      remedies = "Drop the left-hand side, for example `~ group + age`.")
  }
  variables <- all.vars(terms)
  absent <- setdiff(variables, names(data))
  if (length(absent)) {
    .capability_refusal(sprintf(paste0(
      "The group model names %s that `data` does not carry. Group covariates ",
      "are read from the per-participant table only, never from the calling ",
      "environment: a plan whose design depends on a variable in a user's ",
      "workspace could not be re-sealed from what it stores."
    ), .msg_count(length(absent), "variable")),
      capability = "identified_group_model",
      namespace = "population_plans",
      reasons = paste0("model_variable_absent:", absent),
      remedies = paste0(
        "Add ", .msg_names(absent), " to `data`, one value per ",
        "participant."
      ))
  }
  incomplete <- variables[vapply(variables, function(name)
    anyNA(data[[name]]), logical(1))]
  if (length(incomplete)) {
    .capability_refusal(paste0(
      "Group covariates must be complete. crossform will not drop a ",
      "participant to make a design fit: dropping one silently changes which ",
      "population the group estimate is about, and the change would not ",
      "appear anywhere in the result."
    ),
      capability = "identified_group_model",
      namespace = "population_plans",
      reasons = paste0("model_variable_incomplete:", incomplete),
      remedies = paste0(
        "Supply every participant's value, or drop the participant from ",
        "`subjects`, `transport` and `data` together so the plan records the ",
        "smaller population."
      ))
  }
  frame <- stats::model.frame(terms, data = data, na.action = stats::na.fail)
  design <- stats::model.matrix(terms, frame)
  if (!all(is.finite(design))) {
    .capability_refusal(paste0(
      "The group model matrix holds a non-finite entry, so the design is not ",
      "a linear map on the subject axis and no fit against it is defined."
    ),
      capability = "identified_group_model",
      namespace = "population_plans",
      reasons = "model_matrix_not_finite",
      remedies = paste0(
        "Check the covariate transforms in `model` -- a `log()` of a ",
        "non-positive value is the usual cause."
      ))
  }
  factorization <- qr(design)
  rank <- as.integer(factorization$rank)
  columns <- ncol(design)
  if (rank < columns) {
    aliased <- colnames(design)[factorization$pivot[(rank + 1L):columns]]
    .capability_refusal(sprintf(paste0(
      "The group design is rank deficient: %d columns span a %d-dimensional ",
      "space, so its coefficients are not identified and infinitely many ",
      "answers fit the same data equally well. crossform refuses rather than ",
      "dropping columns for you, because which column is dropped decides what ",
      "the remaining coefficients mean."
    ), columns, rank),
      capability = "identified_group_model",
      namespace = "population_plans",
      reasons = c(
        "model_rank_deficient",
        if (nrow(design) < columns) {
          paste0("subjects_fewer_than_columns:", nrow(design), "_of_", columns)
        },
        paste0("aliased_column:", aliased)
      ),
      remedies = paste0(
        "Remove ", .msg_names(aliased), " from `model`, or recode the ",
        "covariate it duplicates. A per-participant factor is the usual ",
        "cause: it saturates the subject axis and leaves nothing for the ",
        "group fit to estimate."
      ))
  }
  text <- gsub("[[:space:]]+", " ",
    trimws(paste(deparse(model, width.cutoff = 500L), collapse = " ")))
  # The stored formula's environment is reset to `baseenv()`: a sealed plan
  # must not keep a user's workspace alive through a promise it never
  # evaluates again, and `effect_space()` refuses environments in provenance
  # for the same reason.
  stored <- model
  environment(stored) <- baseenv()
  used <- names(frame)
  list(
    formula = stored,
    formula_text = text,
    terms = attr(terms, "term.labels"),
    intercept = as.integer(attr(terms, "intercept")),
    variables = used,
    columns = colnames(design),
    assign = as.integer(attr(design, "assign")),
    matrix = design,
    qr = factorization,
    rank = rank,
    pivot = as.integer(factorization$pivot),
    data = data,
    # One digest per model-frame column, keyed by name and taken in a
    # locale-independent name order, so that the identity follows the
    # covariate *content* and not the order the columns happened to sit in.
    data_digest = vapply(
      frame[order(used, method = "radix")],
      function(column) .sha256_signature(unname(column)), character(1)
    )
  )
}

# Identity ---------------------------------------------------------------------

# The scientific identity of a population estimand. Every ingredient is
# something a reader would have to know to say what the group number means:
# which participants' geometry (their plan ids), how each was carried (the
# transport signatures, which already cover the operator's own contents, its
# semantics, its declared row mass and its provenance), what was regressed on
# what (the canonicalized formula and a content digest of the covariates), how
# the incommensurable per-subject budgets were made commensurable, and which
# fit and evaluation order produced it.
#
# The subject and transport vectors are *named*, and the plan has already
# sorted its subjects into a canonical order, so the identity is invariant to
# the order the participants were listed in while still pinning which
# transport belongs to which participant.
.population_plan_scientific_id <- function(subjects, transport, model,
                                           normalization, fit,
                                           allow_nonconservative) {
  .sha256_signature(list(
    schema_version = 1L,
    contract_version = "population-form-v1",
    role = "population_plan",
    subjects = vapply(subjects, function(plan) plan$scientific_plan_id,
      character(1)),
    transports = vapply(transport, function(value) value$signature,
      character(1)),
    model = list(
      formula = model$formula_text,
      terms = model$terms,
      intercept = model$intercept,
      columns = model$columns,
      data = model$data_digest
    ),
    normalization = normalization,
    fit = fit,
    allow_nonconservative = allow_nonconservative
  ), "population-sha256:")
}

# The physical signature, kept separate from the estimand exactly as
# `.geometry_plan_signature()` keeps them apart: a compute policy changes how
# the same group number is produced, never what it is. The audit tables and
# the covariate table ride here rather than in the estimand -- they are
# derived from, or wider than, the sealed ingredients above -- so that
# tampering with a printed record still breaks a signature.
.population_plan_signature <- function(scientific_plan_id, compute,
                                       subject_index, group_index, data) {
  .sha256_signature(list(
    schema_version = 1L,
    scientific_plan_id = scientific_plan_id,
    compute = unclass(compute),
    subject_index = subject_index,
    group_index = group_index,
    data = data
  ))
}

# One row per participant: what was carried, from what frame, under what
# provenance, and how much of the participant's territory never arrived.
# `sink_territory` is a property of the operator alone and needs no data,
# which is what makes it printable from a plan.
.population_subject_index <- function(subjects, transport, conservation) {
  labels <- names(subjects)
  rows <- lapply(labels, function(label) {
    plan <- subjects[[label]]
    carrier <- transport[[label]]
    report <- conservation[[label]]
    sink <- .transport_sink_territory(carrier)
    cross_fit <- carrier$provenance$cross_fit
    data.frame(
      subject = label,
      measurements = as.integer(plan$measurements),
      declared_normalization = report$declared_normalization,
      metric_folded = isTRUE(report$metric_folded),
      conserved = isTRUE(report$conserved),
      max_deviation = as.numeric(report$max_deviation),
      semantics = carrier$semantics,
      provenance = carrier$provenance$method,
      cross_fit = if (is.null(cross_fit)) NA_character_ else
        paste(cross_fit, collapse = "; "),
      sink_territory = as.numeric(sink$share),
      all_sink_rows = as.integer(sink$all_sink),
      row_metadata = all(.population_row_metadata %in%
        names(carrier$native_index)),
      plan_id = plan$scientific_plan_id,
      transport_signature = carrier$signature,
      stringsAsFactors = FALSE
    )
  })
  index <- do.call(rbind, rows)
  rownames(index) <- NULL
  index
}

# Constructor ------------------------------------------------------------------

#' Plan a population form over transported conservative geometry
#'
#' `plan_population()` names a group-level estimand: each participant's
#' compiled conservative geometry, one declared [location_transport()] each
#' onto one shared set of group nodes, and one group-level model over the
#' subject axis. It reads no data and fits nothing --- it is a sealed
#' declaration of what a group number will mean, and a set of refusals for
#' the combinations that would make that meaning false.
#'
#' @section What the plan fixes:
#' Three choices decide what the group number *is*, and all three enter the
#' plan's scientific identity rather than being applied silently.
#'
#' * **The transport.** The group form is not "the population geometry"; it
#'   is the population geometry *as resolved by* `P`. Two transports give two
#'   different, equally valid estimands.
#' * **The normalization.** Participants have incommensurable native budgets,
#'   so combining them requires a declared choice: `"none"` is the mean
#'   subject ledger in native evidence units (participants with more evidence
#'   weigh more), and `"unit_budget"` is the mean attribution *share*, in
#'   which every participant contributes budget one regardless of how much
#'   evidence they have. On the contract's own fixture the two disagree by up
#'   to 69 % and pick different group nodes as the maximum.
#' * **The fit and its order.** The group fit is OLS. Subject-constant
#'   weights are the commuting class: query, transport and fit are then linear
#'   maps on distinct axes and every evaluation order gives the same answer.
#'   With heterogeneous native frames there is no common node axis at all, so
#'   transport must precede the fit; the order is recorded rather than chosen.
#'
#' @section Refusals:
#' Each is an `effect_capability_refusal` carrying the missing capability, all
#' unmet reasons, and a remedy (see [catch_refusal()]).
#'
#' * `common_experimental_space` --- participants' geometry indexed by
#'   different effects, orders, bases, units or scales.
#' * `symmetric_packed_population` --- a rectangular cross-axis subject plan.
#' * `fixed_metric_population` --- a subject plan on a learned neural metric.
#' * `conservative_subject_geometry` --- a subject frame that is not
#'   conservative, either by declaration or after a diagonal metric fold,
#'   unless `allow_nonconservative = TRUE`.
#' * `aligned_subject_transports` --- transport and subject names that do not
#'   match, or a transport whose native rows do not count the participant's
#'   frame measurements.
#' * `shared_group_nodes` --- transports landing on different group node sets,
#'   or declaring different `semantics`.
#' * `aligned_subject_rows` --- a `data` table whose rows do not name their
#'   participants.
#' * `identified_group_model` --- a two-sided formula, an absent or incomplete
#'   covariate, a non-finite design, or a rank-deficient design (the refusal
#'   names the aliased columns).
#' * `precision_weighted_normalization` --- see below.
#'
#' @section The `precision_weighted` gate:
#' `"precision_weighted"` weighs each participant by
#' \eqn{\pi_i = 1/\mathrm{Var}(\hat T_i)}, the precision of a *conserved
#' budget*. That variance needs the full cross-node sampling covariance rather
#' than a sum of per-node margins, and crossform has no route that produces
#' one: [sampling_covariance()] reads a query bank at one measurement and
#' refuses `scope = "cross_measurement"` for want of a spatial model. The mode
#' is therefore refused rather than approximated --- overlapping conservative
#' rows are strongly correlated, so a precision assembled from per-node
#' margins would reweight participants by the wrong quantity, and the contract
#' measures the three normalizations disagreeing by up to 94 %. The mode stays
#' in the closed set, and in this argument's default, because the set *is* the
#' plan identity; what is gated is admitting it.
#'
#' @param subjects A named list of compiled `effect_geometry_plan` values from
#'   [plan_geometry()], one per participant, over conservative frames on a
#'   common experimental space. Names are participant identifiers and bind
#'   `transport` and `data`. The list is stored in a canonical name order, so
#'   the plan's identity does not depend on the order it was given in.
#' @param transport A named list of [location_transport()] values, one per
#'   participant with matching names, all landing on one shared group node set
#'   under one `semantics`.
#' @param model A one-sided group-level formula evaluated against one row per
#'   participant. The default `~ 1` is the group mean.
#' @param data A data frame with one row per participant, identified by a
#'   `subject` column or by row names holding the participant identifiers.
#'   `NULL` (the default) supplies the empty table an intercept-only model
#'   needs.
#' @param normalization How incommensurable per-participant budgets are made
#'   commensurable: `"none"` (the default), `"unit_budget"`, or
#'   `"precision_weighted"` (gated, see above).
#' @param compute A [compute_policy()]. It enters the plan's physical
#'   signature, never its scientific identity.
#' @param allow_nonconservative Plan on non-conservative subject frames
#'   anyway. Recorded on the plan and part of its scientific identity, because
#'   the result is a different estimand rather than a relaxed setting on the
#'   same one.
#' @return An `effect_population_plan`: a sealed list carrying `$subjects`,
#'   `$transport`, the shared `$group_index`, the transport `$semantics`, the
#'   `$model` record (formula, canonical text, term labels, model matrix, its
#'   pivoted `$qr`, `$rank` and `$pivot`, and per-covariate content digests),
#'   `$data`, `$normalization`, the `$fit` marker, `$allow_nonconservative`,
#'   a per-participant `$subject_index` audit table, `$compute`, the
#'   `$scientific_plan_id` naming the estimand, and the `$signature` covering
#'   how it will be executed.
#' @references `design/population-form-contract.md` (`population-form-v1`),
#'   sections 1--4 and 9.
#' @family population transports
#' @seealso [location_transport()] and [anatomical_transport()] for the
#'   operator, [transport_values()] for the contraction it performs, and
#'   [plan_geometry()] for the per-participant input.
#' @examples
#' # Two participants, different native frames, one shared pair of group nodes.
#' effects <- effect_space(c("face", "house"), basis_id = "conditions:v1")
#' subject_plan <- function(id, n) {
#'   domain <- abstract_domain(n, coordinates = cbind(x = seq_len(n) - 1),
#'     feature_ids = paste0("f", seq_len(n)), id = id)
#'   values <- function() matrix(seq_len(2 * n) / n, 2, n,
#'     dimnames = list(c("face", "house"), NULL))
#'   rel <- relation(list(run1 = values(), run2 = values() / 2),
#'     effects = effects, domain = domain)
#'   plan_geometry(rel, compile_frame(voxelwise(), domain),
#'     cross_partitions(rel))
#' }
#' carrier <- function(n) anatomical_transport(
#'   native_coords = cbind(seq_len(n) - 1), group_coords = cbind(c(0, 3)),
#'   semantics = "budget"
#' )
#'
#' plan <- plan_population(
#'   subjects = list(s01 = subject_plan("s01", 5L), s02 = subject_plan("s02", 6L)),
#'   transport = list(s01 = carrier(5L), s02 = carrier(6L))
#' )
#' plan
#' plan$subject_index[, c("subject", "declared_normalization", "sink_territory")]
#'
#' # A group covariate: rows name the participant they belong to.
#' covariates <- data.frame(age = c(24, 31), row.names = c("s01", "s02"))
#' with_age <- plan_population(
#'   subjects = list(s01 = subject_plan("s01", 5L), s02 = subject_plan("s02", 6L)),
#'   transport = list(s01 = carrier(5L), s02 = carrier(6L)),
#'   model = ~ age, data = covariates
#' )
#' with_age$model$columns
#'
#' # The transport is part of the estimand, so changing it moves the identity.
#' identical(plan$scientific_plan_id, with_age$scientific_plan_id)
#' @export
plan_population <- function(subjects, transport, model = ~ 1, data = NULL,
                            normalization = c("none", "unit_budget",
                              "precision_weighted"),
                            compute = compute_policy(),
                            allow_nonconservative = FALSE) {
  normalization <- match.arg(normalization)
  .check_flag(allow_nonconservative, "allow_nonconservative")
  compute <- .validate_compute_policy(compute)

  subjects <- .population_subjects(subjects)
  .population_admit_subject_plans(subjects)
  conservation <- lapply(subjects, function(plan) frame_conservation(plan$frame))
  .population_admit_conservation(subjects, conservation, allow_nonconservative)

  transport <- .population_transports(transport, subjects)
  group_index <- .population_group_nodes(transport)
  semantics <- transport[[1L]]$semantics

  if (identical(normalization, "precision_weighted")) {
    .population_precision_refusal()
  }

  model <- .population_model(model, data, names(subjects))

  # `population-form-v1` section 3.3: OLS is the default group fit because
  # subject-constant weights are the commuting class, and with heterogeneous
  # native frames there is no common node axis, so the transport must precede
  # the fit. Both facts are recorded, and both enter the estimand, so that a
  # later weighted variant -- which gives up query or aggregation commutation
  # -- cannot inherit this plan's identity.
  fit <- list(
    kind = "ols",
    weights = "subject_constant",
    commuting = TRUE,
    evaluation_order = "transport_then_fit"
  )
  subject_index <- .population_subject_index(subjects, transport, conservation)
  scientific_plan_id <- .population_plan_scientific_id(
    subjects, transport, model, normalization, fit, allow_nonconservative
  )
  value <- structure(list(
    subjects = subjects,
    transport = transport,
    group_index = group_index,
    semantics = semantics,
    model = model[setdiff(names(model), "data")],
    data = model$data,
    normalization = normalization,
    fit = fit,
    allow_nonconservative = allow_nonconservative,
    subject_index = subject_index,
    compute = compute,
    scientific_plan_id = scientific_plan_id,
    signature = .population_plan_signature(
      scientific_plan_id, compute, subject_index, group_index, model$data
    )
  ), class = "effect_population_plan")
  .validate_population_plan(value, deep = FALSE)
  value
}

.population_plan_fields <- c(
  "subjects", "transport", "group_index", "semantics", "model", "data",
  "normalization", "fit", "allow_nonconservative", "subject_index",
  "compute", "scientific_plan_id", "signature"
)

.population_model_fields <- c(
  "formula", "formula_text", "terms", "intercept", "variables", "columns",
  "assign", "matrix", "qr", "rank", "pivot", "data_digest"
)

.validate_population_plan <- function(x, deep = TRUE) {
  if (!inherits(x, "effect_population_plan")) {
    .input_error(sprintf(paste0(
      "Expected an `effect_population_plan` from `plan_population()`; ",
      "received %s."
    ), .msg_value(x)), arg = "x", received = .msg_value(x),
      expected = "an `effect_population_plan` from `plan_population()`")
  }
  if (.validated_before(x, "population_plan", deep)) return(invisible(x))
  model <- x$model
  if (!.sealed_fields(x, "effect_population_plan", .population_plan_fields) ||
      !is.list(x$subjects) || !is.list(x$transport) ||
      !identical(names(x$subjects), names(x$transport)) ||
      !is.data.frame(x$group_index) || !is.data.frame(x$data) ||
      !is.data.frame(x$subject_index) ||
      nrow(x$subject_index) != length(x$subjects) ||
      nrow(x$data) != length(x$subjects) ||
      !.is_string(x$semantics) || !x$semantics %in% c("budget", "density") ||
      !.is_string(x$normalization) ||
      !x$normalization %in% .population_normalizations ||
      !.is_flag(x$allow_nonconservative) ||
      !identical(names(model), .population_model_fields) ||
      !inherits(model$formula, "formula") || !.is_string(model$formula_text) ||
      !is.matrix(model$matrix) || !inherits(model$qr, "qr") ||
      !is.integer(model$rank) || !is.integer(model$pivot) ||
      !identical(model$columns, colnames(model$matrix)) ||
      nrow(model$matrix) != length(x$subjects) ||
      !identical(x$fit, list(kind = "ols", weights = "subject_constant",
        commuting = TRUE, evaluation_order = "transport_then_fit")) ||
      !.strong_sha256(sub("^population-", "", x$scientific_plan_id)) ||
      !.strong_sha256(x$signature)) {
    .input_error("Population-plan fields are missing or noncanonical.")
  }
  if (model$rank != ncol(model$matrix)) {
    .contract_error(paste0(
      "A sealed population plan carries an identified group design; this one ",
      "records a rank below its column count."
    ))
  }
  for (plan in x$subjects) .validate_geometry_plan(plan, deep = deep)
  for (value in x$transport) .validate_location_transport(value)
  compute <- .validate_compute_policy(x$compute)
  .check_signature(
    x$scientific_plan_id,
    .population_plan_scientific_id(
      x$subjects, x$transport, model, x$normalization, x$fit,
      x$allow_nonconservative
    ),
    "Population-plan scientific identity is inconsistent."
  )
  .check_signature(
    x$signature,
    .population_plan_signature(
      x$scientific_plan_id, compute, x$subject_index, x$group_index, x$data
    ),
    "Population-plan execution signature is inconsistent."
  )
  .record_validated(x, "population_plan", deep)
  invisible(x)
}
