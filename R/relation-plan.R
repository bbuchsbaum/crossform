# First-moment relation plans and execution --------------------------------

.relation_plan_aliases <- function(design, tolerance) {
  decomposition <- svd(design, nu = min(dim(design)), nv = ncol(design))
  keep <- decomposition$d > tolerance * max(decomposition$d)
  if (all(keep)) return(character())
  names <- colnames(design)
  null <- decomposition$v[, !keep, drop = FALSE]
  unique(vapply(seq_len(ncol(null)), function(index) {
    values <- null[, index]
    active <- abs(values) > tolerance * 10 * max(1, max(abs(values)))
    paste(names[active], collapse = ", ")
  }, character(1)))
}

.assert_model_rows <- function(study, model) {
  for (partition in study$partitions) {
    observed <- study$observations$indexes[[partition]]$observation_id
    declared <- model$row_ids[[partition]]
    if (!identical(observed, declared)) {
      .study_refusal(
        sprintf("Partition `%s` design rows do not bind its observation axis.",
          partition),
        capability = "aligned_observations",
        reasons = c(
          paste0("Observation ids: ", paste(observed, collapse = ", "), "."),
          paste0("Design row ids: ", paste(declared, collapse = ", "), ".")
        ),
        remedies = paste0(
          "Compile the design against the study observation ids in their ",
          "declared order."
        )
      )
    }
  }
  invisible(TRUE)
}

.validate_design_receipt <- function(value) {
  expected <- c(
    "contract_version", "partition", "study_id", "design_model_id",
    "effect_map_id", "observation_model_id", "compiler", "design",
    "coefficient_axis", "lowered_target", "effect_space", "row_lineage",
    "censoring", "solver_policy", "solver", "tolerance", "rank",
    "aliases", "residual_df", "observation_whitener", "capabilities",
    "design_receipt_id"
  )
  if (!inherits(value, "effect_design_receipt") || !is.list(value) ||
      !identical(names(value), expected)) {
    stop("Design-receipt fields are missing or noncanonical.", call. = FALSE)
  }
  if (!is.matrix(value$design) || !is.matrix(value$lowered_target) ||
      !is.character(value$coefficient_axis) ||
      !identical(colnames(value$design), value$coefficient_axis) ||
      !identical(colnames(value$lowered_target), value$coefficient_axis) ||
      !inherits(value$effect_space, "effect_space") ||
      !is.data.frame(value$row_lineage) || !is.list(value$capabilities)) {
    stop("Design-receipt values or axes are invalid.", call. = FALSE)
  }
  .validate_effect_space(value$effect_space)
  semantic <- c(list(schema_version = 1L), value[setdiff(
    names(value), "design_receipt_id"
  )])
  if (!identical(value$design_receipt_id,
      .semantic_digest("design-receipt-sha256:", semantic))) {
    stop("Design-receipt identity is inconsistent.", call. = FALSE)
  }
  value
}

.compile_relation_receipts <- function(study, model, effects,
                                       observation_model,
                                       tolerance) {
  partitions <- study$partitions
  full_rows <- stats::setNames(vapply(partitions, function(partition) {
    length(study$observations$indexes[[partition]]$observation_id)
  }, integer(1)), partitions)
  retained_rows <- stats::setNames(lapply(partitions, function(partition) {
    which(study$lineage[[partition]]$retained)
  }), partitions)
  if (any(lengths(retained_rows) < 1L)) {
    missing <- partitions[lengths(retained_rows) < 1L]
    .study_refusal(
      "Censoring removed every observation from a partition.",
      capability = "estimable_effects",
      reasons = paste0("No retained rows in partition: ", missing, "."),
      remedies = "Revise the censor policy or omit the empty partition."
    )
  }
  observation <- .observation_model_for_partitions(
    observation_model, partitions, full_rows, retained_rows
  )
  lowered <- receipts <- stats::setNames(
    vector("list", length(partitions)), partitions
  )
  raw <- inherits(model, "effect_raw_design_model")
  if (raw) effects <- .validate_lowered_effect_map(effects)

  for (partition in partitions) {
    rows <- retained_rows[[partition]]
    design <- model$designs[[partition]][rows, , drop = FALSE]
    lowered[[partition]] <- if (raw) {
      if (!inherits(effects, "effect_raw_map")) {
        stop("Raw design models require `raw_effect_map()`.", call. = FALSE)
      }
      if (!identical(colnames(effects$target), colnames(design))) {
        .capability_refusal(
          sprintf("Partition `%s` raw target and design columns differ.",
            partition),
          capability = "valid_effect_lowering",
          namespace = "relation_plan",
          reasons = c(
            paste0("Design: ", paste(colnames(design), collapse = ", "), "."),
            paste0("Target: ", paste(colnames(effects$target), collapse = ", "),
              ".")
          ),
          remedies = "Supply a raw target on the exact design coefficient axis."
        )
      }
      effects
    } else {
      lower_effect_map(effects, model$parameterizations[[partition]])
    }
    whitener <- .resolve_observation_whitener(
      observation$whiteners[[partition]], NULL, nrow(design),
      legacy_supplied = FALSE
    )
    compiled <- .compile_lm_estimator(
      design,
      lowered[[partition]]$target,
      whitener,
      lowered[[partition]]$effect_space,
      tolerance,
      partition = partition,
      solver = model$solver[[partition]]
    )
    diagnostics <- compiled$extractor$diagnostics
    lineage <- study$lineage[[partition]][rows, , drop = FALSE]
    rownames(lineage) <- NULL
    censoring <- list(
      declared = !is.null(study$confounds) &&
        !is.null(study$confounds$censor_column),
      column = if (is.null(study$confounds)) {
        NULL
      } else {
        study$confounds$censor_column
      },
      input_rows = as.integer(full_rows[[partition]]),
      retained_rows = as.integer(rows),
      excluded_rows = as.integer(setdiff(
        seq_len(full_rows[[partition]]), rows
      ))
    )
    capabilities <- list(
      aligned_observations = TRUE,
      complete_row_lineage = TRUE,
      estimable_effects = TRUE,
      symbolic_effects = !raw,
      coding_invariant = !raw,
      fixed_observation_model =
        observation_model$capabilities$fixed_observation_model,
      learned_observation_model =
        observation_model$capabilities$learned_observation_model,
      residual_channel =
        observation_model$capabilities$fixed_observation_model,
      analytic_effect_covariance =
        observation_model$capabilities$analytic_effect_covariance,
      portable_design_receipt = TRUE
    )
    semantic <- list(
      schema_version = 1L,
      contract_version = "first-moment-relation-v1",
      partition = partition,
      study_id = study$study_id,
      design_model_id = model$design_model_id,
      effect_map_id = lowered[[partition]]$effect_map_id,
      observation_model_id = observation_model$observation_model_id,
      compiler = model$compiler,
      design = design,
      coefficient_axis = colnames(design),
      lowered_target = lowered[[partition]]$target,
      effect_space = lowered[[partition]]$effect_space,
      row_lineage = lineage,
      censoring = censoring,
      solver_policy = model$solver[[partition]],
      solver = diagnostics$solver,
      tolerance = as.numeric(tolerance),
      rank = diagnostics$rank,
      aliases = .relation_plan_aliases(design, tolerance),
      residual_df = compiled$residual_df,
      observation_whitener = compiled$observation_whitener,
      capabilities = capabilities
    )
    receipts[[partition]] <- structure(c(semantic[-1L], list(
      design_receipt_id = .semantic_digest("design-receipt-sha256:", semantic)
    )), class = "effect_design_receipt")
  }
  list(
    receipts = receipts,
    lowered = lowered,
    retained_rows = retained_rows,
    sampling_unit = observation$sampling_unit,
    whiteners = observation$whiteners
  )
}

#' Plan a first-moment experimental--neural relation
#'
#' The plan composes factual study identity, a semantic or raw design model,
#' named effects, and the observation-model assumptions. Design compilation and
#' estimability are validated without reading neural values.
#'
#' @param study A [study()].
#' @param model A [design_model()] or [raw_design_model()].
#' @param effects An [effect_map()] for a semantic model, or [raw_effect_map()]
#'   for a raw design.
#' @param observation_model An [observation_model()].
#' @param tolerance Positive numerical rank and estimability tolerance.
#' @return An `effect_relation_plan`.
#' @export
plan_relation <- function(study, model, effects, observation_model,
                          tolerance = sqrt(.Machine$double.eps)) {
  study <- .validate_study(study)
  model <- .validate_design_model(model)
  observation_model <- .validate_observation_model(observation_model)
  if (!identical(model$partitions, study$partitions)) {
    .study_refusal(
      "Design-model partitions do not equal study partitions.",
      capability = "aligned_observations",
      reasons = c(
        paste0("Study: ", paste(study$partitions, collapse = ", "), "."),
        paste0("Model: ", paste(model$partitions, collapse = ", "), ".")
      ),
      remedies = "Compile one design for every study partition in study order."
    )
  }
  .assert_model_rows(study, model)
  raw <- inherits(model, "effect_raw_design_model")
  if (raw) {
    effects <- .validate_lowered_effect_map(effects)
    if (!inherits(effects, "effect_raw_map")) {
      stop("A raw design requires `raw_effect_map()`.", call. = FALSE)
    }
  } else {
    effects <- .validate_effect_map(effects)
    if (!identical(effects$condition_space$signature,
        model$condition_space$signature)) {
      .capability_refusal(
        "The design model and effect map bind different condition spaces.",
        capability = "valid_effect_lowering",
        namespace = "relation_plan",
        reasons = c(
          paste0("Model: ", model$condition_space$signature, "."),
          paste0("Effects: ", effects$condition_space$signature, ".")
        ),
        remedies = "Declare the model and effects on one condition space."
      )
    }
  }
  if (!is.numeric(tolerance) || length(tolerance) != 1L ||
      is.na(tolerance) || !is.finite(tolerance) || tolerance <= 0) {
    stop("`tolerance` must be one positive finite number.", call. = FALSE)
  }
  effect_map_id <- effects$effect_map_id
  semantic <- list(
    schema_version = 1L,
    contract_version = "first-moment-relation-v1",
    study_id = study$study_id,
    design_model_id = model$design_model_id,
    effect_map_id = effect_map_id,
    observation_model_id = observation_model$observation_model_id
  )
  relation_plan_id <- .semantic_digest("relation-plan-sha256:", semantic)
  compiled <- .compile_relation_receipts(
    study, model, effects, observation_model, tolerance
  )
  capabilities <- list(
    symbolic_model = !raw,
    symbolic_effects = !raw,
    coding_invariant = !raw,
    aligned_observations = TRUE,
    portable_design_receipts = TRUE,
    fixed_observation_model =
      observation_model$capabilities$fixed_observation_model,
    learned_observation_model =
      observation_model$capabilities$learned_observation_model,
    analytic_effect_covariance =
      observation_model$capabilities$analytic_effect_covariance
  )
  structure(list(
    study = study,
    model = model,
    effects = effects,
    observation_model = observation_model,
    partitions = study$partitions,
    lowered_effects = compiled$lowered,
    design_receipts = compiled$receipts,
    retained_rows = compiled$retained_rows,
    sampling_unit = compiled$sampling_unit,
    whiteners = compiled$whiteners,
    tolerance = as.numeric(tolerance),
    capabilities = capabilities,
    relation_plan_id = relation_plan_id
  ), class = "effect_relation_plan")
}

.validate_relation_plan <- function(value) {
  expected <- c(
    "study", "model", "effects", "observation_model", "partitions",
    "lowered_effects", "design_receipts", "retained_rows", "sampling_unit",
    "whiteners", "tolerance", "capabilities", "relation_plan_id"
  )
  if (!inherits(value, "effect_relation_plan") || !is.list(value) ||
      !identical(names(value), expected)) {
    stop("Relation-plan fields are missing or noncanonical.", call. = FALSE)
  }
  rebuilt <- plan_relation(
    value$study, value$model, value$effects, value$observation_model,
    tolerance = value$tolerance
  )
  if (!identical(value, rebuilt)) {
    stop("Relation-plan metadata or identity is inconsistent.", call. = FALSE)
  }
  rebuilt
}

#' Inspect portable design receipts
#'
#' @param x An [plan_relation()] result.
#' @return A named list of `effect_design_receipt` values.
#' @export
relation_plan_receipts <- function(x) {
  x <- .validate_relation_plan(x)
  lapply(x$design_receipts, .validate_design_receipt)
}

.planned_observation_sources <- function(plan) {
  observations <- plan$study$observations
  sources <- lapply(plan$partitions, function(partition) {
    source <- observations$sources[[partition]]
    rows <- plan$retained_rows[[partition]]
    force(source)
    force(rows)
    function(features) source$read(features)[rows, , drop = FALSE]
  })
  names(sources) <- plan$partitions
  dimensions <- lapply(plan$partitions, function(partition) {
    c(length(plan$retained_rows[[partition]]), observations$n_features)
  })
  names(dimensions) <- plan$partitions
  capabilities <- lapply(plan$partitions, function(partition) {
    source <- observations$capabilities[[partition]]
    revision <- .semantic_digest("sha256:", list(
      source_revision = source$stable_revision,
      retained_rows = plan$retained_rows[[partition]],
      design_receipt_id =
        plan$design_receipts[[partition]]$design_receipt_id
    ))
    source_capabilities(
      block_read = source$block_read,
      reopenable = FALSE,
      thread_safe = FALSE,
      stable_revision = revision
    )
  })
  names(capabilities) <- plan$partitions
  list(sources = sources, dimensions = dimensions, capabilities = capabilities)
}

#' Estimate a planned relation family
#'
#' This is the authorized first-moment execution verb. It emits the existing
#' `effect_relation_fit` contract, so every second-moment plan and view remains
#' unchanged.
#'
#' @param x An [plan_relation()] result.
#' @param ... Reserved for future execution policies.
#' @return An `effect_relation_fit` whose identity binds the relation plan,
#'   source revisions, realized row lineage, and design receipts.
#' @export
estimate_relation <- function(x, ...) {
  if (length(list(...))) {
    stop(
      "No additional `estimate_relation()` arguments are currently supported.",
      call. = FALSE
    )
  }
  plan <- .validate_relation_plan(x)
  sources <- .planned_observation_sources(plan)
  designs <- lapply(plan$design_receipts, `[[`, "design")
  targets <- lapply(plan$design_receipts, `[[`, "lowered_target")
  fit_provenance <- list(
    contract_version = "first-moment-relation-v1",
    relation_plan_id = plan$relation_plan_id,
    design_receipt_ids = vapply(
      plan$design_receipts, `[[`, character(1), "design_receipt_id"
    ),
    study_id = plan$study$study_id,
    observation_model_id = plan$observation_model$observation_model_id,
    source_revisions = vapply(
      sources$capabilities, `[[`, character(1), "stable_revision"
    )
  )
  fitted <- lm_relation_fit(
    sources$sources,
    design = designs,
    effects = targets,
    observation_whitener = plan$whiteners,
    effect_names = plan$lowered_effects[[1L]]$effect_space,
    tolerance = plan$tolerance,
    source_dims = sources$dimensions,
    partitions = plan$partitions,
    domain = plan$study$observations$domain,
    capabilities = sources$capabilities,
    sampling_unit = plan$sampling_unit,
    provenance = fit_provenance,
    solver = plan$model$solver
  )
  if (isTRUE(plan$observation_model$capabilities$learned_observation_model)) {
    fitted <- relation_fit(
      fitted$relation,
      error_models = NULL,
      provenance = c(fit_provenance, list(
        learned_observation_model = TRUE,
        training_revision = plan$observation_model$training_revision,
        analytic_error_channel = "withheld"
      ))
    )
  }
  .validate_relation_fit(fitted)
  fitted
}

#' @export
format.effect_relation_plan <- function(x, ...) {
  x <- .validate_relation_plan(x)
  sprintf(
    "relation_plan<%d partitions; %s; uncertainty %s>",
    length(x$partitions),
    if (isTRUE(x$capabilities$symbolic_model)) "semantic" else "raw-X-T",
    if (isTRUE(x$capabilities$analytic_effect_covariance)) {
      "analytic"
    } else {
      "withheld"
    }
  )
}

#' @export
print.effect_relation_plan <- function(x, ...) {
  cat(format(x), "\n", sep = "")
  invisible(x)
}
