# Sampling-covariance plan identities ---------------------------------------

.sampling_record <- function(kind, fields) {
  if (!is.character(kind) || length(kind) != 1L || is.na(kind) ||
      !nzchar(kind) || !is.list(fields) || is.null(names(fields)) ||
      any(!nzchar(names(fields))) || anyDuplicated(names(fields))) {
    stop("Sampling-plan record fields are invalid.", call. = FALSE)
  }
  semantic <- c(list(schema_version = 1L, kind = kind), fields)
  structure(c(list(kind = kind), fields, list(
    signature = paste0("sha256:", digest::digest(
      semantic, algo = "sha256", serialize = TRUE, serializeVersion = 2L
    ))
  )), class = "effect_sampling_record")
}

.validate_sampling_record <- function(x, kind = NULL) {
  if (!inherits(x, "effect_sampling_record") || !is.list(x) ||
      length(x) < 2L || !identical(names(x)[[1L]], "kind") ||
      !identical(utils::tail(names(x), 1L), "signature") ||
      !is.character(x$kind) || length(x$kind) != 1L || is.na(x$kind) ||
      !nzchar(x$kind) || !.strong_sha256(x$signature) ||
      (!is.null(kind) && !identical(x$kind, kind))) {
    stop("Sampling-plan record is missing or noncanonical.", call. = FALSE)
  }
  fields <- x[setdiff(names(x), c("kind", "signature"))]
  rebuilt <- .sampling_record(x$kind, unclass(fields))
  if (!identical(x, rebuilt)) {
    stop("Sampling-plan record identity is inconsistent.", call. = FALSE)
  }
  x
}

.sampling_operation <- function(
    kind = c("diagonal", "selected_entries", "apply", "quadratic_form",
             "transport", "materialize"), argument = NULL) {
  kind <- match.arg(kind)
  needs_argument <- kind %in% c(
    "selected_entries", "apply", "quadratic_form", "transport"
  )
  if (needs_argument && is.null(argument)) {
    stop(sprintf("Sampling operation `%s` requires its exact query argument.",
      kind), call. = FALSE)
  }
  if (!needs_argument && !is.null(argument)) {
    stop(sprintf("Sampling operation `%s` does not accept a query argument.",
      kind), call. = FALSE)
  }
  if (!is.null(argument)) {
    if (!(is.atomic(argument) || is.matrix(argument)) ||
        anyNA(argument) || (is.numeric(argument) && any(!is.finite(argument)))) {
      stop("Sampling-operation arguments must be finite atomic values.",
        call. = FALSE)
    }
    # Preserve named query axes: names and dimnames identify evidence inputs
    # and transported outputs and are therefore semantic, not display noise.
  }
  .sampling_record("operation", list(
    operation = kind,
    argument = argument
  ))
}

.sampling_materialization <- function(
    kind = c("query_only", "dense_covariance")) {
  kind <- match.arg(kind)
  .sampling_record("materialization", list(mode = kind))
}

.sampling_target <- function(
    kind = c("point_estimate", "null", "linear_hypothesis"),
    policy = NULL, specification = NULL) {
  kind <- match.arg(kind)
  if (is.null(policy)) {
    policy <- switch(kind,
      point_estimate = "nonnegative_plugin",
      null = "fixed_zero",
      linear_hypothesis = "declared_operator"
    )
  }
  if (!is.character(policy) || length(policy) != 1L || is.na(policy) ||
      !nzchar(policy)) {
    stop("A sampling target requires one named evaluation policy.",
      call. = FALSE)
  }
  if (kind == "linear_hypothesis" && is.null(specification)) {
    stop("A linear-hypothesis target requires an exact specification.",
      call. = FALSE)
  }
  if (kind != "linear_hypothesis" && !is.null(specification)) {
    stop("Only a linear-hypothesis target accepts a specification.",
      call. = FALSE)
  }
  if (!is.null(specification) &&
      (!(is.atomic(specification) || is.matrix(specification)) ||
       anyNA(specification) ||
       (is.numeric(specification) && any(!is.finite(specification))))) {
    stop("Sampling-target specifications must be finite atomic values.",
      call. = FALSE)
  }
  .sampling_record("target", list(
    target = kind,
    policy = policy,
    specification = if (is.null(specification)) NULL else unname(specification)
  ))
}

.sampling_evidence_descriptor <- function(x) {
  if (inherits(x, "effect_geometry_plan")) {
    .validate_geometry_plan(x, deep = FALSE)
    relation <- x$task$left_relation
    metric_schedule <- x$metric_schedule
    metric_status <- if (!is.null(metric_schedule$metric) &&
        isTRUE(metric_schedule$metric$capabilities$learned_frozen)) {
      "learned"
    } else {
      "fixed"
    }
    metric_identity <- metric_schedule$signature
    evidence_kind <- "geometry_plan"
    plan_id <- x$scientific_plan_id
    frame_signature <- .additive_frame_signature(x$frame)
    pairing <- x$pairing
  } else if (inherits(x, "effect_crossnobis_plan")) {
    .validate_crossnobis_plan(x, deep = FALSE)
    relation <- x$task$left_relation
    metric_schedule <- x$metric_schedule
    metric_status <- if (isTRUE(metric_schedule$capabilities$learned) ||
        isTRUE(metric_schedule$capabilities$
          calibration_requires_metric_uncertainty)) {
      "learned"
    } else {
      "fixed"
    }
    metric_identity <- metric_schedule$signature
    evidence_kind <- "crossnobis_plan"
    plan_id <- x$scientific_plan_id
    frame_signature <- .additive_frame_signature(x$frame)
    pairing <- x$pairing
  } else {
    stop(paste0(
      "Evidence sampling currently binds a compiled geometry or crossnobis ",
      "plan, not an uncompiled matrix or relation."
    ), call. = FALSE)
  }
  record <- .sampling_record("evidence", list(
    evidence_kind = evidence_kind,
    plan_id = plan_id,
    task_id = x$task$task_id,
    relation_id = .relation_family_identity(relation),
    same_relation = x$task$same_relation,
    frame_signature = frame_signature,
    pairing = .metric_pairing_identity(pairing),
    metric_identity = metric_identity,
    metric_status = metric_status
  ))
  list(record = record, relation = relation, pairing = pairing)
}

.sampling_error_channel <- function(x, relation) {
  .validate_relation(relation, deep = FALSE)
  relation_id <- .relation_family_identity(relation)
  if (inherits(x, "effect_relation_fit")) {
    .validate_relation_fit(x, deep = FALSE)
    if (!identical(.relation_family_identity(x$relation), relation_id)) {
      stop(paste0(
        "The error channel is not identity-bound to the relation used by ",
        "the evidence plan."
      ), call. = FALSE)
    }
    present <- vapply(x$error_models, Negate(is.null), logical(1))
    model_kinds <- unique(vapply(x$error_models[present], `[[`,
      character(1), "kind"))
    sampling_units <- unique(vapply(x$error_models[present], `[[`,
      character(1), "sampling_unit"))
    record <- .sampling_record("error_channel", list(
      channel = if (all(present)) "relation_fit" else "absent",
      relation_id = relation_id,
      channel_identity = x$signature,
      model = if (all(present) && length(model_kinds) == 1L) {
        model_kinds[[1L]]
      } else if (all(present)) {
        "declared_alternative"
      } else {
        "absent"
      },
      partition_models = vapply(x$error_models, function(model) {
        if (is.null(model)) NA_character_ else model$signature
      }, character(1)),
      sampling_units = sampling_units,
      available = all(present)
    ))
    return(list(record = record, source = x))
  }
  if (!inherits(x, "effect_relation")) {
    stop("An error channel must be a relation fit or an explicitly bare relation.",
      call. = FALSE)
  }
  .validate_relation(x, deep = FALSE)
  if (!identical(.relation_family_identity(x), relation_id)) {
    stop("The bare relation differs from the evidence-plan relation.",
      call. = FALSE)
  }
  record <- .sampling_record("error_channel", list(
    channel = "absent",
    relation_id = relation_id,
    channel_identity = NULL,
    model = "absent",
    partition_models = stats::setNames(
      rep(NA_character_, length(x$partitions)), x$partitions
    ),
    sampling_units = character(),
    available = FALSE
  ))
  list(record = record, source = x)
}

.sampling_equal_error_structure <- function(channel) {
  .validate_sampling_record(channel$record, "error_channel")
  if (!isTRUE(channel$record$available) ||
      !inherits(channel$source, "effect_relation_fit")) {
    return(FALSE)
  }
  models <- channel$source$error_models
  reference <- models[[1L]]
  all(vapply(models, function(model) {
    identical(model$kind, reference$kind) &&
      identical(model$scale_convention, reference$scale_convention) &&
      identical(model$sampling_unit, reference$sampling_unit) &&
      identical(model$residual_df, reference$residual_df) &&
      isTRUE(all.equal(model$effect_covariance,
        reference$effect_covariance, tolerance = 1e-12,
        check.attributes = TRUE))
  }, logical(1)))
}

.sampling_partition_record <- function(pairing, partitions,
                                       model = c("equal", "heterogeneous")) {
  model <- match.arg(model)
  .validate_pairing(pairing)
  partitions <- as.character(partitions)
  if (length(partitions) < 2L || anyNA(partitions) ||
      any(!nzchar(partitions)) || anyDuplicated(partitions)) {
    stop("Sampling partitions must have unique nonempty identities.",
      call. = FALSE)
  }
  expected <- utils::combn(sort(partitions), 2L)
  expected_keys <- paste(expected[1L, ], expected[2L, ], sep = "\r")
  observed_keys <- paste(
    pmin(pairing$left, pairing$right),
    pmax(pairing$left, pairing$right), sep = "\r"
  )
  complete_all_pairs <- !isTRUE(attr(pairing, "directed")) &&
    identical(sort(observed_keys), sort(expected_keys))
  equal_weights <- max(abs(pairing$weight - 1 / nrow(pairing))) <= 1e-14
  .sampling_record("partition", list(
    model = model,
    partitions = partitions,
    count = as.integer(length(partitions)),
    pairing = .metric_pairing_identity(pairing),
    complete_all_pairs = complete_all_pairs,
    equal_weights = equal_weights,
    endpoint_independence = identical(
      attr(pairing, "independence", exact = TRUE), "independent"
    ),
    edge_products_independent = FALSE
  ))
}

.sampling_metric_record <- function(identity, status = c("fixed", "learned"),
                                    uncertainty = NULL,
                                    role = c("metric", "bridge")) {
  status <- match.arg(status)
  role <- match.arg(role)
  if (!is.character(identity) || length(identity) != 1L || is.na(identity) ||
      !nzchar(identity)) {
    stop("A sampling metric or bridge requires one stable identity.",
      call. = FALSE)
  }
  if (status == "fixed") {
    if (is.null(uncertainty)) uncertainty <- "not_applicable"
    if (!identical(uncertainty, "not_applicable")) {
      stop("A fixed metric must declare uncertainty `not_applicable`.",
        call. = FALSE)
    }
  } else {
    if (is.null(uncertainty) ||
        !uncertainty %in% c("ignored", "propagated")) {
      stop("A learned metric must declare uncertainty `ignored` or `propagated`.",
        call. = FALSE)
    }
  }
  .sampling_record("metric", list(
    role = role,
    status = status,
    uncertainty = uncertainty,
    identity = identity
  ))
}

.sampling_plan_semantic <- function(evidence, error_channel, metric,
                                    partition, sampling_axis, target,
                                    spatial_scope, operation,
                                    materialization, capabilities, reasons) {
  list(
    schema_version = 1L,
    contract = "evidence-sampling-v1",
    evidence = evidence$signature,
    error_channel = error_channel$signature,
    metric = metric$signature,
    partition = partition$signature,
    sampling_axis = sampling_axis,
    target = target$signature,
    spatial_scope = spatial_scope,
    operation = operation$signature,
    materialization = materialization$signature,
    capabilities = capabilities,
    reasons = reasons
  )
}

.sampling_unavailability_reasons <- function(
    evidence, error_channel, metric, partition, sampling_axis,
    spatial_scope, equal_error_structure) {
  reasons <- character()
  if (!isTRUE(error_channel$available)) {
    reasons <- c(reasons, "missing_error_channel")
  }
  if (!identical(error_channel$model, "separable_glm") &&
      isTRUE(error_channel$available)) {
    reasons <- c(reasons, "unsupported_error_model")
  }
  if (!isTRUE(evidence$same_relation)) {
    reasons <- c(reasons, "cross_relation_law_not_admitted")
  }
  if (!identical(metric$role, "metric")) {
    reasons <- c(reasons, "bridge_law_not_admitted")
  }
  if (identical(metric$status, "learned") &&
      !identical(metric$uncertainty, "propagated")) {
    reasons <- c(reasons, "learned_metric_uncertainty_unpropagated")
  }
  if (identical(metric$status, "learned")) {
    reasons <- c(reasons, "learned_metric_law_not_admitted")
  }
  # Equality of the per-partition error structures is a property of the error
  # channel. With no channel there is nothing to compare, so reporting it
  # alongside `missing_error_channel` would state a second, independent-looking
  # failure that is really the same one, and would contradict the partition
  # model the capabilities record shows.
  if (!identical(partition$model, "equal") ||
      (isTRUE(error_channel$available) && !isTRUE(equal_error_structure))) {
    reasons <- c(reasons, "heterogeneous_partition_model")
  }
  if (!isTRUE(partition$complete_all_pairs) ||
      !isTRUE(partition$equal_weights)) {
    reasons <- c(reasons, "not_equal_weight_all_pairs")
  }
  if (!isTRUE(partition$endpoint_independence)) {
    reasons <- c(reasons, "endpoint_independence_not_declared")
  }
  if (is.null(sampling_axis) || !is.character(sampling_axis) ||
      length(sampling_axis) != 1L || is.na(sampling_axis) ||
      !nzchar(sampling_axis) ||
      (length(error_channel$sampling_units) > 0L &&
       !identical(error_channel$sampling_units, sampling_axis))) {
    reasons <- c(reasons, "sampling_axis_missing_or_inconsistent")
  }
  if (identical(spatial_scope, "modeled")) {
    reasons <- c(reasons, "spatial_covariance_model_unavailable")
  }
  unique(reasons)
}

.compile_evidence_sampling_plan <- function(
    evidence, error_channel = NULL,
    partition_model = c("equal", "heterogeneous"),
    sampling_axis = NULL,
    target = .sampling_target("point_estimate"),
    spatial_scope = c("local_marginals", "modeled"),
    operation = .sampling_operation("diagonal"),
    materialization = .sampling_materialization("query_only"),
    metric_uncertainty = NULL, metric_role = c("metric", "bridge")) {
  partition_model <- match.arg(partition_model)
  spatial_scope <- match.arg(spatial_scope)
  metric_role <- match.arg(metric_role)
  descriptor <- .sampling_evidence_descriptor(evidence)
  if (is.null(error_channel)) error_channel <- descriptor$relation
  channel <- .sampling_error_channel(error_channel, descriptor$relation)
  operation <- .validate_sampling_record(operation, "operation")
  materialization <- .validate_sampling_record(
    materialization, "materialization"
  )
  target <- .validate_sampling_record(target, "target")
  if (operation$operation == "materialize" &&
      materialization$mode != "dense_covariance") {
    stop("A materialize operation requires dense-covariance materialization.",
      call. = FALSE)
  }
  if (operation$operation != "materialize" &&
      materialization$mode == "dense_covariance") {
    stop("Dense covariance materialization requires a materialize operation.",
      call. = FALSE)
  }
  if (is.null(metric_uncertainty) &&
      identical(descriptor$record$metric_status, "learned")) {
    # Nothing in evidence-sampling-v1 propagates metric-estimation
    # uncertainty; record that honestly so the plan compiles and the
    # admission layer refuses with the full reason set instead of a
    # construction-shape error.
    metric_uncertainty <- "ignored"
  }
  metric <- .sampling_metric_record(
    descriptor$record$metric_identity,
    descriptor$record$metric_status,
    metric_uncertainty,
    metric_role
  )
  partition <- .sampling_partition_record(
    descriptor$pairing, descriptor$relation$partitions, partition_model
  )
  if (is.null(sampling_axis) &&
      length(channel$record$sampling_units) == 1L) {
    sampling_axis <- channel$record$sampling_units[[1L]]
  }
  if (is.null(sampling_axis)) {
    # A pairing that declares `generalizes_over` has already named the axis the
    # estimand generalizes across. Reading it here means a plan with no error
    # channel does not additionally report a missing sampling axis it was in
    # fact given.
    declared_axis <- attr(
      descriptor$pairing, "generalizes_over", exact = TRUE
    )
    if (!is.null(declared_axis)) sampling_axis <- declared_axis
  }
  if (!is.null(sampling_axis) &&
      (!is.character(sampling_axis) || length(sampling_axis) != 1L ||
       is.na(sampling_axis) || !nzchar(sampling_axis))) {
    stop("`sampling_axis` must be NULL or one nonempty identity.",
      call. = FALSE)
  }
  equal_error_structure <- .sampling_equal_error_structure(channel)
  reasons <- .sampling_unavailability_reasons(
    descriptor$record, channel$record, metric, partition, sampling_axis,
    spatial_scope, equal_error_structure
  )
  available <- length(reasons) == 0L
  capabilities <- list(
    sampling_covariance = if (available) "available" else "unavailable",
    error_channel = channel$record$channel,
    error_model = channel$record$model,
    metric_role = metric$role,
    metric_status = metric$status,
    metric_uncertainty = metric$uncertainty,
    partition_model = partition$model,
    sampling_axis = sampling_axis,
    calibration_target = target$target,
    spatial_covariance = if (available) "local_marginals" else "unavailable",
    requested_operation = operation$operation,
    materialization = materialization$mode
  )
  semantic <- .sampling_plan_semantic(
    descriptor$record, channel$record, metric, partition, sampling_axis,
    target, spatial_scope, operation, materialization, capabilities, reasons
  )
  scientific_plan_id <- paste0("sampling-sha256:", digest::digest(
    semantic, algo = "sha256", serialize = TRUE, serializeVersion = 2L
  ))
  value <- structure(list(
    contract = "evidence-sampling-v1",
    evidence_plan = evidence,
    error_source = channel$source,
    evidence = descriptor$record,
    error_channel = channel$record,
    metric = metric,
    partition = partition,
    sampling_axis = sampling_axis,
    target = target,
    spatial_scope = spatial_scope,
    operation = operation,
    materialization = materialization,
    capabilities = capabilities,
    unavailable_reasons = reasons,
    scientific_plan_id = scientific_plan_id,
    signature = paste0("sha256:", digest::digest(list(
      schema_version = 1L,
      scientific_plan_id = scientific_plan_id,
      evidence_object = descriptor$record$signature,
      error_source = channel$record$signature
    ), algo = "sha256", serialize = TRUE, serializeVersion = 2L))
  ), class = "effect_evidence_sampling_plan")
  .validate_evidence_sampling_plan(value, deep = TRUE)
  value
}

.rebind_evidence_sampling_plan <- function(x, operation, materialization) {
  .validate_evidence_sampling_plan(x, deep = FALSE)
  .compile_evidence_sampling_plan(
    x$evidence_plan,
    error_channel = x$error_source,
    partition_model = x$partition$model,
    sampling_axis = x$sampling_axis,
    target = x$target,
    spatial_scope = x$spatial_scope,
    operation = operation,
    materialization = materialization,
    metric_uncertainty = x$metric$uncertainty,
    metric_role = x$metric$role
  )
}

.validate_evidence_sampling_plan <- function(x, deep = TRUE) {
  if (.validated_before(x, "evidence_sampling_plan", deep)) return(invisible(x))
  expected <- c("contract", "evidence_plan", "error_source", "evidence",
    "error_channel", "metric", "partition", "sampling_axis", "target",
    "spatial_scope", "operation", "materialization", "capabilities",
    "unavailable_reasons", "scientific_plan_id", "signature")
  if (!inherits(x, "effect_evidence_sampling_plan") || !is.list(x) ||
      !identical(names(x), expected) ||
      !identical(x$contract, "evidence-sampling-v1") ||
      !x$spatial_scope %in% c("local_marginals", "modeled") ||
      !is.list(x$capabilities) || !is.character(x$unavailable_reasons) ||
      !.strong_sha256(sub("^sampling-", "", x$scientific_plan_id)) ||
      !.strong_sha256(x$signature)) {
    stop("Evidence-sampling-plan fields are missing or noncanonical.",
      call. = FALSE)
  }
  evidence <- .sampling_evidence_descriptor(x$evidence_plan)
  channel <- .sampling_error_channel(x$error_source, evidence$relation)
  .validate_sampling_record(x$evidence, "evidence")
  .validate_sampling_record(x$error_channel, "error_channel")
  metric <- .validate_sampling_record(x$metric, "metric")
  partition <- .validate_sampling_record(x$partition, "partition")
  target <- .validate_sampling_record(x$target, "target")
  operation <- .validate_sampling_record(x$operation, "operation")
  materialization <- .validate_sampling_record(
    x$materialization, "materialization"
  )
  if (!identical(x$evidence, evidence$record) ||
      !identical(x$error_channel, channel$record) ||
      !identical(x$metric$identity, evidence$record$metric_identity) ||
      !identical(x$partition$pairing, evidence$record$pairing)) {
    stop("Evidence-sampling-plan sources or identities are inconsistent.",
      call. = FALSE)
  }
  reasons <- .sampling_unavailability_reasons(
    evidence$record, channel$record, metric, partition, x$sampling_axis,
    x$spatial_scope, .sampling_equal_error_structure(channel)
  )
  available <- length(reasons) == 0L
  expected_capabilities <- list(
    sampling_covariance = if (available) "available" else "unavailable",
    error_channel = channel$record$channel,
    error_model = channel$record$model,
    metric_role = metric$role,
    metric_status = metric$status,
    metric_uncertainty = metric$uncertainty,
    partition_model = partition$model,
    sampling_axis = x$sampling_axis,
    calibration_target = target$target,
    spatial_covariance = if (available) "local_marginals" else "unavailable",
    requested_operation = operation$operation,
    materialization = materialization$mode
  )
  if (!identical(x$unavailable_reasons, reasons) ||
      !identical(x$capabilities, expected_capabilities)) {
    stop("Evidence-sampling capabilities are inconsistent.", call. = FALSE)
  }
  if (isTRUE(deep)) {
    semantic <- .sampling_plan_semantic(
      evidence$record, channel$record, metric, partition, x$sampling_axis,
      target, x$spatial_scope, operation, materialization,
      expected_capabilities, reasons
    )
    expected_id <- paste0("sampling-sha256:", digest::digest(
      semantic, algo = "sha256", serialize = TRUE, serializeVersion = 2L
    ))
    expected_signature <- paste0("sha256:", digest::digest(list(
      schema_version = 1L,
      scientific_plan_id = expected_id,
      evidence_object = evidence$record$signature,
      error_source = channel$record$signature
    ), algo = "sha256", serialize = TRUE, serializeVersion = 2L))
    if (!identical(x$scientific_plan_id, expected_id) ||
        !identical(x$signature, expected_signature)) {
      stop("Evidence-sampling-plan identity is inconsistent.", call. = FALSE)
    }
  }
  .record_validated(x, "evidence_sampling_plan", deep)
  invisible(x)
}

.sampling_refusal_details <- function(reason) {
  switch(reason,
    missing_error_channel = list(
      why = paste0(
        "this evidence plan has only a precomputed relation and no error ",
        "channel. Refit raw observations with `lm_relation_fit()` or supply ",
        "a validated, identity-bound external error channel; beta matrices ",
        "alone cannot recover residual uncertainty"
      ),
      remedy = "Refit raw observations with `lm_relation_fit()`."
    ),
    unsupported_error_model = list(
      why = paste0(
        "the retained error channel is not a separable GLM residual model, ",
        "and evidence-sampling-v1 admits no other error law"
      ),
      remedy = paste0(
        "Fit with `lm_relation_fit()`, which records the separable GLM ",
        "channel this law requires."
      )
    ),
    cross_relation_law_not_admitted = list(
      why = paste0(
        "the left and right relations differ, and the analytic sampling law ",
        "is admitted only for self-forms built from one relation"
      ),
      remedy = "Build the plan from a single relation on both sides."
    ),
    bridge_law_not_admitted = list(
      why = paste0(
        "the neural operator plays a bridge role, and no analytic sampling ",
        "law is admitted for bridge readouts"
      ),
      remedy = "Use a same-space neural metric instead of a bridge."
    ),
    learned_metric_uncertainty_unpropagated = list(
      why = paste0(
        "the neural metric was learned and its estimation uncertainty is ",
        "not propagated"
      ),
      remedy = paste0(
        "Use one common fixed metric, or retain the learned-metric result ",
        "as a signed point estimate without an analytic law."
      )
    ),
    learned_metric_law_not_admitted = list(
      why = paste0(
        "evidence-sampling-v1 admits only a common fixed metric, and this ",
        "plan's neural metric was learned"
      ),
      remedy = paste0(
        "Use a geometry plan with one common fixed metric, or retain the ",
        "learned-metric result as a signed point estimate."
      )
    ),
    heterogeneous_partition_model = list(
      why = paste0(
        "the equal-partition analytic law is not valid for this partition ",
        "model or for partitions with unequal error structure"
      ),
      remedy = paste0(
        "Use partitions sharing one error model, or treat the result as a ",
        "point estimate."
      )
    ),
    not_equal_weight_all_pairs = list(
      why = paste0(
        "the partition pairing is not the complete, equally weighted set of ",
        "unordered partition pairs that the analytic law covers"
      ),
      remedy = paste0(
        "Pair with `cross_partitions(..., independence = \"independent\")` ",
        "over all partitions when that assumption is justified."
      )
    ),
    endpoint_independence_not_declared = list(
      why = paste0(
        "the pairing does not declare its partition endpoints statistically ",
        "independent"
      ),
      remedy = paste0(
        "Declare `cross_partitions(..., independence = \"independent\")` or ",
        "`pairing(..., independence = \"independent\")` only when the ",
        "partition endpoints are truly independent."
      )
    ),
    sampling_axis_missing_or_inconsistent = list(
      why = paste0(
        "no single sampling axis is declared, or the declared axis ",
        "conflicts with the error channel's recorded sampling unit"
      ),
      remedy = paste0(
        "Declare one `sampling_axis` that matches the error channel's ",
        "sampling unit."
      )
    ),
    spatial_covariance_model_unavailable = list(
      why = paste0(
        "cross-location sampling covariance requires an explicit spatial ",
        "model; local marginal variances do not create one"
      ),
      remedy = paste0(
        "Request `spatial_scope = \"local_marginals\"`; no spatial ",
        "covariance model is admitted yet."
      )
    ),
    list(
      why = sprintf("requirement `%s` is unmet", reason),
      remedy = character()
    )
  )
}

.require_sampling_covariance <- function(x) {
  .validate_evidence_sampling_plan(x, deep = FALSE)
  if (identical(x$capabilities$sampling_covariance, "available")) {
    return(invisible(TRUE))
  }
  reasons <- x$unavailable_reasons
  details <- lapply(reasons, .sampling_refusal_details)
  why <- vapply(details, function(detail) detail$why, character(1))
  remedies <- unique(unlist(lapply(details, function(detail) detail$remedy)))
  message <- if (length(why) == 1L) {
    paste0("Sampling covariance is unavailable because ", why, ".")
  } else {
    paste0(
      "Sampling covariance is unavailable; ", length(why),
      " requirements of the admitted analytic law are unmet:\n",
      paste0("* ", why, ".", collapse = "\n")
    )
  }
  .capability_refusal(
    message,
    capability = "sampling_covariance",
    namespace = "evidence_sampling",
    reasons = reasons,
    remedies = remedies
  )
}
