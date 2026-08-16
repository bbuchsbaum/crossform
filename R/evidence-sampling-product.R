# Relation-fit specialization of evidence-sampling-v1 ----------------------

.sampling_distance_contrasts <- function(effects) {
  if (!is.character(effects) || length(effects) < 2L || anyNA(effects) ||
      any(!nzchar(effects)) || anyDuplicated(effects)) {
    stop("Sampling distances require at least two identified effects.",
      call. = FALSE)
  }
  pairs <- utils::combn(seq_along(effects), 2L)
  value <- matrix(0, ncol(pairs), length(effects))
  for (distance in seq_len(ncol(pairs))) {
    value[distance, pairs[1L, distance]] <- 1
    value[distance, pairs[2L, distance]] <- -1
  }
  labels <- paste(effects[pairs[1L, ]], effects[pairs[2L, ]], sep = " - ")
  dimnames(value) <- list(labels, effects)
  list(
    matrix = value,
    labels = labels,
    pairs = data.frame(
      left = effects[pairs[1L, ]],
      right = effects[pairs[2L, ]],
      stringsAsFactors = FALSE
    )
  )
}

.sampling_common_effect_covariance <- function(fit) {
  .validate_relation_fit(fit, deep = FALSE)
  models <- fit$error_models
  reference <- models[[1L]]$effect_covariance
  if (!all(vapply(models, function(model) {
    isTRUE(all.equal(model$effect_covariance, reference,
      tolerance = 1e-12, check.attributes = TRUE))
  }, logical(1)))) {
    stop(paste0(
      "The fixed-metric Eq. 13 specialization requires one equal ",
      "effect-coordinate covariance across partitions."
    ), call. = FALSE)
  }
  reference
}

.sampling_target_signal <- function(plan, node, metric) {
  target <- plan$target
  if (target$target == "null") {
    return(matrix(0, length(plan$error_source$relation$effects),
      nrow(metric$value), dimnames = list(
        plan$error_source$relation$effects, node$support
      )))
  }
  if (target$target != "point_estimate" ||
      !identical(target$policy, "partition_mean_plugin")) {
    stop(paste0(
      "The fixed-metric Eq. 13 product path requires target `null` or ",
      "point-estimate policy `partition_mean_plugin`; no target is chosen ",
      "silently."
    ), call. = FALSE)
  }
  relation <- plan$error_source$relation
  estimates <- lapply(relation$partitions, function(partition) {
    relation_block(relation, partition, node$support_positions)
  })
  Reduce(`+`, estimates) / length(estimates)
}

.sampling_local_metric <- function(evidence_plan, node) {
  schedule <- evidence_plan$metric_schedule
  if (identical(schedule$kind, "implicit_identity_before_frame")) {
    base <- neural_metric(
      diag(length(node$support)), evidence_plan$frame$domain,
      support = node$support,
      inverse = diag(length(node$support)),
      provenance = list(construction = "implicit_identity_sampling_metric")
    )
    return(.compose_frame_metric(
      evidence_plan$frame, base, node$position
    )$metric)
  }
  metric <- .validate_neural_metric(schedule$metric, deep = FALSE)
  positions <- match(node$support, metric$support)
  if (anyNA(positions)) {
    stop("A sampling frame support falls outside its fixed neural metric.",
      call. = FALSE)
  }
  restriction_provenance <- list(
    construction = "fixed_metric_support_restriction",
    parent_metric = metric$signature
  )
  # A learned-frozen parent keeps its frozen training provenance through the
  # support restriction; dropping it here would make `neural_metric()` refuse
  # with a provenance message that contradicts what the user supplied.
  inherited <- intersect(
    c("frozen", "training_signature"), names(metric$provenance)
  )
  restriction_provenance[inherited] <- metric$provenance[inherited]
  base <- neural_metric(
    metric$value[positions, positions, drop = FALSE],
    metric$domain, node$support,
    inverse = NULL,
    estimation = metric$estimation,
    tolerance = metric$tolerance,
    provenance = restriction_provenance
  )
  .compose_frame_metric(evidence_plan$frame, base, node$position)$metric
}

.sampling_whitened_components <- function(signal, metric,
                                          residual_covariance) {
  metric <- .validate_neural_metric(metric, deep = FALSE)
  if (!metric$capabilities$positive_definite) {
    stop("Sampling covariance requires a positive-definite fixed metric.",
      call. = FALSE)
  }
  factor <- chol(metric$value)
  # chol() returns R with K = R'R. Row patterns transform by R' so their
  # Euclidean inner product is exactly b K b'.
  signal_whitened <- signal %*% t(factor)
  residual_whitened <- factor %*% residual_covariance %*% t(factor)
  list(signal = signal_whitened, residual = residual_whitened)
}

.compile_fixed_metric_rdm_sampling <- function(
    evidence, error_channel,
    target = .sampling_target(
      "point_estimate", policy = "partition_mean_plugin"
    ),
    operation = .sampling_operation("diagonal"),
    materialization = .sampling_materialization("query_only")) {
  plan <- .compile_evidence_sampling_plan(
    evidence, error_channel = error_channel,
    partition_model = "equal", target = target,
    spatial_scope = "local_marginals", operation = operation,
    materialization = materialization
  )
  .require_sampling_covariance(plan)
  plan
}

.fixed_metric_rdm_sampling_resources <- function(
    plan, residual_workspace_bytes = 512 * 1024^2,
    strategy = c("node_local", "shared_pair_statistics")) {
  .validate_evidence_sampling_plan(plan, deep = FALSE)
  .require_sampling_covariance(plan)
  strategy <- match.arg(strategy)
  evidence <- plan$evidence_plan
  fit <- plan$error_source
  if (!inherits(evidence, "effect_geometry_plan") ||
      !inherits(fit, "effect_relation_fit")) {
    stop("Sampling resources require one geometry plan and relation fit.",
      call. = FALSE)
  }
  statistics <- if (identical(strategy, "shared_pair_statistics") &&
      !is.null(evidence$frame$support_index)) {
    residual_pair_statistics(
      fit, evidence$frame, workspace_bytes = residual_workspace_bytes
    )
  } else if (identical(strategy, "shared_pair_statistics")) {
    stop(paste0(
      "Shared pair statistics require a frame with an explicit support ",
      "index; use node-local resources for this frame."
    ), call. = FALSE)
  } else {
    NULL
  }
  structure(list(
    plan = plan$scientific_plan_id,
    relation_fit = fit$signature,
    frame = .additive_frame_signature(evidence$frame),
    strategy = strategy,
    residual_statistics = statistics,
    signature = .sha256_signature(list(
      schema_version = 1L,
      plan = plan$scientific_plan_id,
      relation_fit = fit$signature,
      frame = .additive_frame_signature(evidence$frame),
      strategy = strategy,
      residual_statistics = if (is.null(statistics)) NULL else
        statistics$signature
    ))
  ), class = "effect_sampling_resources")
}

.validate_fixed_metric_rdm_sampling_resources <- function(x, plan) {
  expected <- c("plan", "relation_fit", "frame", "strategy",
    "residual_statistics", "signature")
  .validate_evidence_sampling_plan(plan, deep = FALSE)
  if (!inherits(x, "effect_sampling_resources") || !is.list(x) ||
      !identical(names(x), expected) ||
      !identical(x$plan, plan$scientific_plan_id) ||
      !identical(x$relation_fit, plan$error_source$signature) ||
      !identical(x$frame,
        .additive_frame_signature(plan$evidence_plan$frame)) ||
      !x$strategy %in% c("node_local", "shared_pair_statistics") ||
      !.strong_sha256(x$signature)) {
    stop("Sampling-resource fields or plan identity are inconsistent.",
      call. = FALSE)
  }
  if (!is.null(x$residual_statistics)) {
    .validate_residual_pair_statistics(x$residual_statistics, deep = FALSE)
  }
  if (identical(x$strategy, "node_local") &&
      !is.null(x$residual_statistics)) {
    stop("Node-local sampling resources cannot contain frame-wide pair statistics.",
      call. = FALSE)
  }
  if (identical(x$strategy, "shared_pair_statistics") &&
      is.null(x$residual_statistics)) {
    stop("Shared-pair sampling resources are missing their pair statistics.",
      call. = FALSE)
  }
  expected_signature <- .sha256_signature(list(
    schema_version = 1L,
    plan = x$plan,
    relation_fit = x$relation_fit,
    frame = x$frame,
    strategy = x$strategy,
    residual_statistics = if (is.null(x$residual_statistics)) NULL else
      x$residual_statistics$signature
  ))
  if (!identical(x$signature, expected_signature)) {
    stop("Sampling-resource identity is inconsistent.", call. = FALSE)
  }
  x
}

.sampling_node_residual_covariance <- function(plan, node_value, resources) {
  fit <- plan$error_source
  if (!is.null(resources$residual_statistics)) {
    statistics <- resources$residual_statistics
    index <- plan$evidence_plan$frame$support_index
    return(.local_residual_scope_covariance(
      statistics, index, node_value$support_positions,
      fit$relation$partitions
    ))
  }
  residuals <- lapply(fit$relation$partitions, function(partition) {
    crossprod(residual_block(fit, partition, node_value$support_positions))
  })
  total_df <- sum(vapply(fit$relation$partitions, function(partition) {
    residual_df(fit, partition)
  }, integer(1)))
  Reduce(`+`, residuals) / total_df
}

.fixed_metric_rdm_sampling_covariance <- function(
    plan, node = 1L, resources = NULL) {
  .validate_evidence_sampling_plan(plan)
  .require_sampling_covariance(plan)
  evidence <- plan$evidence_plan
  if (!inherits(evidence, "effect_geometry_plan") ||
      plan$evidence$metric_status != "fixed") {
    stop("This specialization requires a fixed-metric geometry plan.",
      call. = FALSE)
  }
  .validate_geometry_plan(evidence, deep = FALSE)
  read_node <- .frame_metric_node_accessor(evidence$frame)
  node_value <- read_node(node)
  if (is.null(resources)) {
    resources <- .fixed_metric_rdm_sampling_resources(plan)
  }
  resources <- .validate_fixed_metric_rdm_sampling_resources(
    resources, plan
  )
  metric <- .sampling_local_metric(evidence, node_value)
  fit <- plan$error_source
  effect_covariance <- .sampling_common_effect_covariance(fit)
  residual_covariance <- .sampling_node_residual_covariance(
    plan, node_value, resources
  )
  total_df <- sum(vapply(fit$relation$partitions, function(partition) {
    residual_df(fit, partition)
  }, integer(1)))
  distances <- .sampling_distance_contrasts(fit$relation$effects)
  signal <- .sampling_target_signal(plan, node_value, metric)
  whitened <- .sampling_whitened_components(
    signal, metric, residual_covariance
  )
  .sampling_covariance_from_components(
    plan,
    contrasts = distances$matrix,
    signal_patterns = whitened$signal,
    effect_covariance = effect_covariance,
    residual_covariance = whitened$residual,
    normalization = 1,
    labels = distances$labels,
    source = list(
      specialization = "fixed_metric_equal_partition_rdm",
      node = as.character(node_value$id),
      support = node_value$support,
      metric = metric$signature,
      effect_covariance = fit$error_models[[1L]]$signature,
      residual_df = as.integer(total_df),
      target = plan$target$signature,
      local_only = TRUE
    )
  )
}

.execute_fixed_metric_rdm_sampling <- function(plan, node = 1L,
                                               max_bytes = 512 * 1024^2,
                                               resources = NULL) {
  covariance <- .fixed_metric_rdm_sampling_covariance(
    plan, node, resources
  )
  .execute_evidence_sampling_plan(plan, covariance, max_bytes)
}

#' Construct exact analytic sampling covariance for crossvalidated distances
#'
#' Builds the local sampling-covariance form from Diedrichsen, Provost, and
#' Zareamoghaddam (2016), Eq. 13, for an equal-weight all-partition-pairs RDM.
#' The admitted specialization requires a geometry plan built from an
#' `lm_relation_fit()` relation, equal partition error structures, independent
#' partition endpoints, and one common fixed neural metric. The result remains
#' factorized and queryable; it does not allocate the complete distance-by-
#' distance covariance.
#'
#' Writing \eqn{\mu_r} for the whitened contrast pattern of distance \eqn{r},
#' \eqn{\Sigma_w} for the residual covariance in the same whitened
#' coordinates, \eqn{\Xi} for the effect-coordinate contrast cross-products,
#' and \eqn{M} for the partition count, the law evaluated here is
#' \deqn{\mathrm{Cov}(d_r, d_s) = \frac{4}{M}\,\Xi_{rs}\,
#'   \mu_r\Sigma_w\mu_s^\top
#'   + \frac{2}{M(M-1)}\,\Xi_{rs}^2\,\mathrm{tr}(\Sigma_w\Sigma_w).}
#' The residual covariance enters the signal term as a metric on the whitened
#' patterns, so the law is correct for a general anisotropic \eqn{\Sigma_w}
#' and not only for the spherical case.
#'
#' @param x An `effect_geometry_plan` using `cross_partitions()` and a fixed
#'   neural metric.
#' @param fit The identity-bound `effect_relation_fit` that supplied
#'   `x$task$left_relation`.
#' @param target Whether the signal-dependent term is evaluated at the
#'   partition-mean plug-in estimate (`"plugin"`) or at the fixed zero null
#'   (`"null"`). These are different calibration policies and are never
#'   chosen implicitly. `"plugin"` substitutes the partition mean of the
#'   *estimates* for the unknown signal, and because
#'   \eqn{E[\hat\mu_r\Sigma_w\hat\mu_s^\top] = \mu_r\Sigma_w\mu_s^\top +
#'   \Xi_{rs}\,\mathrm{tr}(\Sigma_w\Sigma_w)/M}, its signal term is biased
#'   upward by \eqn{4\Xi_{rs}^2\mathrm{tr}(\Sigma_w\Sigma_w)/M^2}, an
#'   inflation that shrinks like \eqn{1/M^2} and is largest when the noise
#'   dominates the true distances. Prefer `"null"` for calibrating a test of
#'   no effect, where it is exact; use `"plugin"` when reporting uncertainty
#'   around an estimated nonzero distance and read it as mildly conservative.
#' @param at One measurement position or identifier in the compiled frame.
#' @param residual_strategy `"node_local"` reads residual blocks only for the
#'   requested measurement. `"shared_pair_statistics"` explicitly compiles
#'   reusable residual pair sufficient statistics for a batch of overlapping
#'   supports.
#' @param residual_workspace_bytes Positive crossform-owned workspace budget
#'   for shared residual pair statistics.
#' @return An `effect_rdm_sampling_covariance`, queryable by
#'   [sampling_covariance()]. It contains within-measurement uncertainty only;
#'   it does not imply covariance between spatial locations.
#' @references Diedrichsen J, Provost S, Zareamoghaddam H (2016),
#'   "On the distribution of cross-validated Mahalanobis distances",
#'   especially Eqs. 10, 13, and 35 and Section 5.1.
#'   \doi{10.48550/arXiv.1607.01371}
#' @seealso [sampling_covariance()] to query the result,
#'   [sampling_capabilities()] to ask whether the law is available before
#'   provoking a refusal, and [rdm()] for the point estimates it describes.
#' @family sampling uncertainty
#' @examples
#' set.seed(23)
#' design <- model.matrix(~ 0 + factor(rep(c("a", "b"), each = 4)))
#' colnames(design) <- c("a", "b")
#' effects <- diag(2)
#' rownames(effects) <- colnames(design)
#' truth <- rbind(a = c(0.4, 0.1, 0), b = c(0, 0.2, 0.5))
#' responses <- setNames(lapply(seq_len(3), function(run) {
#'   design %*% truth + matrix(rnorm(8 * 3, sd = 0.2), 8, 3)
#' }), paste0("run", seq_len(3)))
#' domain <- abstract_domain(3, id = "sampling-example")
#' fit <- lm_relation_fit(
#'   responses, design, effects, effect_names = rownames(truth),
#'   sampling_unit = "trial", domain = domain
#' )
#' metric <- noise_precision(
#'   diag(3), domain, covariance = diag(3),
#'   provenance = list(source = "fixed-example-metric")
#' )
#' plan <- plan_geometry(
#'   fit$relation, compile_frame(whole_brain(), domain),
#'   cross_partitions(fit$relation, independence = "independent"),
#'   metric = metric
#' )
#' uncertainty <- rdm_sampling_covariance(plan, fit, target = "null")
#' sqrt(sampling_covariance(uncertainty))
#' @export
rdm_sampling_covariance <- function(
    x, fit,
    target,
    at = 1L,
    residual_strategy = c("node_local", "shared_pair_statistics"),
    residual_workspace_bytes = 512 * 1024^2) {
  if (missing(fit)) {
    stop(paste0(
      "`fit` is required: pass the `lm_relation_fit()` whose `$relation` ",
      "built `x`, so the analytic law can read its residual error channel."
    ), call. = FALSE)
  }
  if (missing(target)) {
    .capability_refusal(paste0(
      "Choose `target = \"plugin\"` for the partition-mean signal policy ",
      "or `target = \"null\"` for the fixed-zero null; the ",
      "signal-dependent covariance target is not inferred."
    ),
      capability = "calibration_target_declared",
      namespace = "evidence_sampling",
      reasons = "calibration_target_not_declared",
      remedies = "Pass `target = \"plugin\"` or `target = \"null\"`."
    )
  }
  descriptor <- .sampling_evidence_descriptor(x)
  if (identical(descriptor$record$metric_status, "learned")) {
    .capability_refusal(paste0(
      "Analytic RDM sampling covariance is unavailable because this plan's ",
      "neural metric was learned. Version 0.1 does not propagate metric-",
      "estimation uncertainty; use a geometry plan with one common fixed ",
      "metric or retain this result as a signed point estimate."
    ),
      capability = "fixed_metric_sampling_law",
      namespace = "evidence_sampling",
      reasons = "learned_metric_law_not_admitted",
      remedies = paste0(
        "Use one common fixed metric, or keep the signed point estimate ",
        "without an analytic law."
      )
    )
  }
  if (inherits(x, "effect_crossnobis_plan")) {
    .capability_refusal(paste0(
      "Analytic RDM sampling covariance is not routed for crossnobis plans; ",
      "compile the same fixed metric into `plan_geometry(metric = )` to use ",
      "the admitted fixed-metric law."
    ),
      capability = "fixed_metric_sampling_law",
      namespace = "evidence_sampling",
      reasons = "crossnobis_plan_not_routed",
      remedies = paste0(
        "Build a geometry plan with `plan_geometry(metric = )` carrying the ",
        "same fixed metric."
      )
    )
  }
  target <- match.arg(target, c("plugin", "null"))
  residual_strategy <- match.arg(residual_strategy)
  # Checked here, against the plan the caller passed, so the message can name
  # the plan's own measurement count instead of a compiled node position.
  at <- .msg_measurement_index(
    at,
    if (is.null(x$measurements)) nrow(x$frame$weights) else x$measurements,
    argument = "at", subject = "plan"
  )
  target_record <- if (identical(target, "null")) {
    .sampling_target("null")
  } else {
    .sampling_target("point_estimate", policy = "partition_mean_plugin")
  }
  plan <- .compile_fixed_metric_rdm_sampling(
    x, fit, target = target_record
  )
  resources <- .fixed_metric_rdm_sampling_resources(
    plan, residual_workspace_bytes = residual_workspace_bytes,
    strategy = residual_strategy
  )
  value <- .fixed_metric_rdm_sampling_covariance(
    plan, node = at, resources = resources
  )
  class(value) <- c("effect_rdm_sampling_covariance", class(value))
  value
}

#' Ask whether the analytic sampling law is available before provoking it
#'
#' A scientist should be able to ask what the uncertainty channel will grant
#' instead of discovering refusals one at a time. This inspection compiles
#' the evidence-sampling admission for a plan and error source and reports
#' every unmet requirement with its remedy, without computing anything.
#'
#' @param x A compiled `effect_geometry_plan` (or crossnobis plan).
#' @param fit Optional `effect_relation_fit` supplying the error channel;
#'   omitting it probes the plan's bare relation, which has no channel.
#' @return An `effect_sampling_capabilities` list: `available`, the full
#'   `capabilities` record, and a `reasons` data frame with one row per
#'   unmet requirement (`reason`, `why`, `remedy`).
#' @seealso [rdm_sampling_covariance()], the call this inspection describes,
#'   and [catch_refusal()] for branching on a refusal that has already
#'   happened.
#' @family sampling uncertainty
#' @examples
#' # Ask before provoking: the fixture's fit retains a residual channel, so
#' # the analytic law is admitted.
#' example <- example_fmri_effects()
#' plan <- plan_geometry(
#'   example$fit$relation, example$frame,
#'   cross_partitions(example$fit$relation, independence = "independent")
#' )
#' capabilities <- sampling_capabilities(plan, example$fit)
#' capabilities$available
#' capabilities
#'
#' # Probing the bare relation instead reports every unmet requirement with
#' # its remedy, rather than failing one refusal at a time.
#' without_fit <- sampling_capabilities(plan)
#' without_fit$reasons
#' @export
sampling_capabilities <- function(x, fit = NULL) {
  plan <- .compile_evidence_sampling_plan(x, error_channel = fit)
  reasons <- plan$unavailable_reasons
  details <- lapply(reasons, .sampling_refusal_details)
  structure(list(
    available = identical(
      plan$capabilities$sampling_covariance, "available"
    ),
    capabilities = plan$capabilities,
    reasons = data.frame(
      reason = as.character(reasons),
      why = vapply(details, function(detail) detail$why, character(1)),
      remedy = vapply(details, function(detail) {
        if (length(detail$remedy)) detail$remedy else NA_character_
      }, character(1)),
      stringsAsFactors = FALSE
    )
  ), class = "effect_sampling_capabilities")
}

#' @export
print.effect_sampling_capabilities <- function(x, ...) {
  cat("<effect_sampling_capabilities>\n")
  cat("  analytic sampling law:",
    if (x$available) "available" else "unavailable", "\n")
  cat("  metric:", x$capabilities$metric_status,
    "| partitions:", x$capabilities$partition_model,
    "| error channel:", x$capabilities$error_channel, "\n")
  if (nrow(x$reasons)) {
    cat("  unmet requirements:\n")
    for (row in seq_len(nrow(x$reasons))) {
      cat("  *", x$reasons$reason[[row]], "-", x$reasons$why[[row]], "\n")
      if (!is.na(x$reasons$remedy[[row]])) {
        cat("      remedy:", x$reasons$remedy[[row]], "\n")
      }
    }
    if (identical(x$capabilities$error_channel, "absent")) {
      cat("  note: requirements that describe the error channel itself",
        "cannot be\n        evaluated until one exists, and are not listed.\n")
    }
  }
  invisible(x)
}

#' Query an exact factorized sampling-covariance form
#'
#' Applies exact linear covariance operations without materializing the full
#' covariance unless `operation = "materialize"` is requested explicitly.
#' A returned diagonal is a vector of variances, not standard errors or an
#' automatic confidence interval.
#'
#' The values returned inherit the calibration target chosen when `x` was
#' built. Under `target = "plugin"` they carry the documented upward bias of
#' the partition-mean plug-in policy; see the `target` argument of
#' [rdm_sampling_covariance()].
#'
#' @param x An object from [rdm_sampling_covariance()].
#' @param operation One exact covariance operation.
#' @param query Operation-specific argument: a two-column integer matrix for
#'   `selected_entries`, an evidence-coordinate vector for `apply` or
#'   `quadratic_form`, or an output-by-evidence matrix for `transport`.
#' @param max_bytes Positive payload/workspace budget used only for explicit
#'   dense materialization.
#' @return A named variance vector (`"diagonal"`), selected covariance vector,
#'   covariance action, scalar quadratic form, transported covariance, or
#'   explicitly materialized covariance matrix, according to `operation`.
#'   Names come from the distance labels carried by `x`.
#' @seealso [rdm_sampling_covariance()] to build `x`, and
#'   [sampling_capabilities()] to check the law is available first.
#' @family sampling uncertainty
#' @examples
#' # Three conditions in three runs, with a fixed identity noise metric.
#' set.seed(11)
#' conditions <- c("face", "house", "tool")
#' design <- model.matrix(~ 0 + factor(rep(conditions, each = 3)))
#' colnames(design) <- conditions
#' effects <- diag(3)
#' rownames(effects) <- conditions
#' truth <- rbind(
#'   face = c(0.6, 0.1, 0), house = c(0, 0.5, 0.2), tool = c(0.1, 0, 0.6)
#' )
#' responses <- setNames(lapply(1:3, function(run) {
#'   design %*% truth + matrix(rnorm(9 * 3, sd = 0.3), 9, 3)
#' }), paste0("run", 1:3))
#' domain <- abstract_domain(3, id = "sampling-covariance-example")
#' fit <- lm_relation_fit(
#'   responses, design, effects, effect_names = conditions,
#'   sampling_unit = "trial", domain = domain
#' )
#' plan <- plan_geometry(
#'   fit$relation, compile_frame(whole_brain(), domain),
#'   cross_partitions(fit$relation, independence = "independent"),
#'   metric = noise_precision(diag(3), domain, covariance = diag(3))
#' )
#' uncertainty <- rdm_sampling_covariance(plan, fit, target = "null")
#'
#' # The diagonal is a vector of variances, so take a square root yourself if
#' # you want standard errors. No interval is implied.
#' round(sqrt(sampling_covariance(uncertainty)), 4)
#'
#' # Distances that share a condition covary; read one entry without
#' # building the full matrix.
#' uncertainty$labels
#' round(sampling_covariance(
#'   uncertainty, "selected_entries", query = cbind(1, 2)
#' ), 6)
#'
#' # The variance of a fixed contrast of distances, again without the matrix.
#' round(sampling_covariance(
#'   uncertainty, "quadratic_form", query = c(1, -1, 0)
#' ), 6)
#'
#' # Materialization is explicit, never a silent fallback.
#' round(sampling_covariance(uncertainty, "materialize"), 6)
#' @export
sampling_covariance <- function(
    x,
    operation = c("diagonal", "selected_entries", "apply",
      "quadratic_form", "transport", "materialize"),
    query = NULL,
    max_bytes = 512 * 1024^2) {
  if (!inherits(x, "effect_rdm_sampling_covariance")) {
    stop("`x` must come from `rdm_sampling_covariance()`.", call. = FALSE)
  }
  .validate_sampling_covariance(x, deep = TRUE)
  operation <- match.arg(operation)
  materialization <- if (identical(operation, "materialize")) {
    .sampling_materialization("dense_covariance")
  } else {
    .sampling_materialization("query_only")
  }
  operation_record <- .sampling_operation(operation, query)
  plan <- .rebind_evidence_sampling_plan(
    x$plan, operation_record, materialization
  )
  covariance <- x
  covariance$plan <- plan
  covariance$signature <- .sampling_covariance_signature(covariance)
  .execute_evidence_sampling_plan(plan, covariance, max_bytes)
}
