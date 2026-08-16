# fmridesign compiler adapter ----------------------------------------------

.fmridesign_semantic_specification <- function(model) {
  frame <- model$sampling_frame
  terms <- lapply(model$terms, function(term) {
    hrf <- attr(term, "hrfspec")
    hrf_function <- if (is.null(hrf)) NULL else hrf$hrf
    list(
      term_tag = attr(term, "term_tag"),
      variable = term$varname,
      conditions = unname(term$condition_levels),
      hrf = if (is.null(hrf)) NULL else list(
        name = attr(hrf_function, "name"),
        basis_count = attr(hrf_function, "nbasis"),
        span = attr(hrf_function, "span"),
        parameters = attr(hrf_function, "params"),
        precision = hrf$precision,
        summate = hrf$summate,
        normalize = hrf$normalize
      )
    )
  })
  list(
    grammar = "fmridesign-event-model",
    terms = terms,
    sampling_frame = list(
      block_lengths = unname(frame$blocklens),
      repetition_times = unname(frame$TR),
      start_times = unname(frame$start_time),
      precision = frame$precision
    ),
    drop_empty = model$model_spec$drop_empty,
    convolution_precision = model$model_spec$precision
  )
}

.fmridesign_block_map <- function(partitions, model, block_map) {
  count <- length(model$sampling_frame$blocklens)
  if (is.null(block_map)) {
    if (count != length(partitions)) {
      .capability_refusal(
        "The fmridesign block axis cannot be bound positionally.",
        capability = "aligned_observations",
        namespace = "relation_compiler",
        reasons = sprintf(
          "The model has %d blocks and the study has %d partitions.",
          count, length(partitions)
        ),
        remedies = "Supply a named `block_map` from study partitions to model blocks."
      )
    }
    block_map <- seq_len(count)
    names(block_map) <- partitions
  }
  if (!is.numeric(block_map) || length(block_map) != length(partitions) ||
      is.null(names(block_map)) || !setequal(names(block_map), partitions) ||
      anyNA(block_map) || any(block_map %% 1 != 0) ||
      any(block_map < 1L | block_map > count) || anyDuplicated(block_map)) {
    stop("`block_map` must map every study partition to one unique model block.",
      call. = FALSE)
  }
  stats::setNames(as.integer(block_map[partitions]), partitions)
}

.fmridesign_assert_event_binding <- function(model, study, block_map) {
  record <- study$events
  if (is.null(record) || !isTRUE(record$timing)) {
    .capability_refusal(
      "The fmridesign adapter requires timed study events.",
      capability = "study_bound_compilation",
      namespace = "relation_compiler",
      reasons = "A compiled event model cannot be identified against absent or untimed facts.",
      remedies = "Bind the timed event table used to construct the fmridesign model into `study()`."
    )
  }
  data <- record$data
  partition_values <- as.character(data[[record$partition_column]])
  for (term_name in names(model$terms)) {
    term <- model$terms[[term_name]]
    if (is.null(term$varname) || !term$varname %in% names(data) ||
        !all(term$subset)) {
      .capability_refusal(
        sprintf("fmridesign term `%s` cannot be identified against study facts.",
          term_name),
        capability = "study_bound_compilation",
        namespace = "relation_compiler",
        reasons = paste0(
          "Automatic conformance currently requires an unsubsetted event term ",
          "whose variable is present in the event table."
        ),
        remedies = "Use an unsubsetted event term or add a certified term-specific adapter."
      )
    }
    for (partition in study$partitions) {
      selected <- partition_values == partition
      model_selected <- term$blockids == block_map[[partition]]
      declared_condition <- paste0(
        term$varname, ".", as.character(data[[term$varname]][selected])
      )
      compiled_condition <- term$condition_levels[term$condition_ids[model_selected]]
      same <- identical(as.numeric(data[[record$onset_column]][selected]),
          as.numeric(term$onsets[model_selected])) &&
        identical(as.numeric(data[[record$duration_column]][selected]),
          as.numeric(term$durations[model_selected])) &&
        identical(declared_condition, unname(compiled_condition))
      if (!same) {
        .capability_refusal(
          sprintf("fmridesign term `%s` disagrees with partition `%s` facts.",
            term_name, partition),
          capability = "study_bound_compilation",
          namespace = "relation_compiler",
          reasons = paste0(
            "Compiled onsets, durations, or condition labels differ from the ",
            "event record bound into the study."
          ),
          remedies = "Recompile the fmridesign model from the study event facts."
        )
      }
    }
  }
  invisible(TRUE)
}

.fmridesign_semantic_map <- function(design, conditions, explicit = NULL) {
  coefficients <- colnames(design)
  if (!is.null(explicit)) {
    if (!is.matrix(explicit) || !is.numeric(explicit) ||
        !identical(rownames(explicit), conditions$coordinates) ||
        !identical(colnames(explicit), coefficients)) {
      stop(paste0(
        "An explicit `semantic_map` must be a numeric condition-by-coefficient ",
        "matrix on the declared axes."
      ), call. = FALSE)
    }
    return(explicit)
  }
  metadata <- as.data.frame(fmridesign::design_meta(design),
    stringsAsFactors = FALSE)
  if (!all(c("name", "condition", "basis_ix", "basis_total") %in%
      names(metadata)) || !identical(metadata$name, coefficients)) {
    .capability_refusal(
      "The fmridesign column metadata are insufficient for symbolic lowering.",
      capability = "valid_effect_lowering",
      namespace = "relation_compiler",
      reasons = "Required name, condition, and basis metadata are absent or misaligned.",
      remedies = "Supply an explicit condition-by-coefficient `semantic_map`."
    )
  }
  map <- matrix(0, length(conditions$coordinates), length(coefficients),
    dimnames = list(conditions$coordinates, coefficients))
  for (coordinate in conditions$coordinates) {
    candidate <- which(metadata$condition == coordinate &
      (is.na(metadata$basis_ix) | metadata$basis_ix == 1L))
    if (length(candidate) != 1L) {
      .capability_refusal(
        sprintf("Condition `%s` has %d primary fmridesign columns.",
          coordinate, length(candidate)),
        capability = "valid_effect_lowering",
        namespace = "relation_compiler",
        reasons = "Automatic lowering requires exactly one primary-basis coefficient per condition.",
        remedies = "Supply an explicit condition-by-coefficient `semantic_map`."
      )
    }
    map[coordinate, candidate] <- 1
  }
  map
}

#' Compile a semantic design model from fmridesign
#'
#' This exact-version adapter binds a compiled `fmridesign::event_model()` to
#' the event facts and observation rows of a study. Scientific condition
#' effects remain declared in condition space; concrete design coding is a
#' compilation receipt.
#'
#' @param model A `fmridesign` event model from the certified package version.
#' @param study A bound [study()] containing the exact timed event facts.
#' @param basis_id,units,scale Scientific coordinate declarations passed to
#'   [condition_space()].
#' @param block_map Optional named integer map from study partitions to
#'   fmridesign blocks.
#' @param semantic_map Optional condition-by-coefficient matrix, or named list
#'   per partition, for models that cannot be lowered from column metadata.
#' @param specification Optional portable semantic declaration. When omitted,
#'   the adapter derives one from the pinned fmridesign model.
#' @param solver Numerical route recorded in the design receipt.
#' @return An [design_model()] suitable for [plan_relation()], whose
#'   `$condition_space` holds the fmridesign condition names, whose `$designs`
#'   hold one compiled matrix per study partition with study observation ids as
#'   row names, and whose `$compiler` records the pinned fmridesign version.
#' @family relation planning and fitting
#' @seealso [design_model()] for the generic contract, [plan_relation()] for
#'   the next step, and [compiler_conformance()] to check the receipt this
#'   adapter produced.
#' @examples
#' # This adapter is certified against exactly one fmridesign version and
#' # refuses any other, so the example runs only under that version.
#' if (requireNamespace("fmridesign", quietly = TRUE) &&
#'     identical(as.character(utils::packageVersion("fmridesign")), "0.6.0")) {
#'   set.seed(1)
#'   domain <- abstract_domain(3L, id = "fmridesign-example")
#'   index <- observation_index(
#'     1:20, "run-1", time = seq(0, by = 2, length.out = 20L), units = "seconds"
#'   )
#'   events <- data.frame(
#'     partition = "run-1", event_id = paste0("e", 1:4),
#'     onset = c(2, 10, 18, 26), duration = 0,
#'     condition = c("face", "body", "face", "body")
#'   )
#'   facts <- study(
#'     observations(list(`run-1` = matrix(rnorm(60), 20L, 3L)),
#'       list(`run-1` = index), domain),
#'     observation_events(events)
#'   )
#'
#'   # The event model must be the one compiled from these exact facts; the
#'   # adapter re-checks onsets, durations, and condition labels against them.
#'   external <- fmridesign::event_model(
#'     onset ~ fmridesign::hrf(condition), data = events,
#'     block = rep(1L, nrow(events)),
#'     sampling_frame = fmridesign::sampling_frame(
#'       blocklens = 20L, TR = 2, start_time = 0
#'     ),
#'     durations = events$duration
#'   )
#'   model <- fmridesign_design_model(
#'     external, facts, basis_id = "canonical-hrf-amplitude",
#'     units = "percent-signal-change"
#'   )
#'   print(model$condition_space$coordinates)
#'   print(model$capabilities$coding_invariant)
#'
#'   # Effects are now named in condition space, so the contrast below does not
#'   # depend on which columns the HRF compiler emitted.
#'   weights <- rbind(`face-body` = c(1, -1))
#'   colnames(weights) <- model$condition_space$coordinates
#'   plan <- plan_relation(
#'     facts, model, effect_map(weights, model$condition_space),
#'     observation_model("ols", sampling_unit = "scan")
#'   )
#'   print(all(as.matrix(compiler_conformance(plan)[-1L])))
#' }
#' @export
fmridesign_design_model <- function(
    model, study, basis_id, units, scale = 1, block_map = NULL,
    semantic_map = NULL, specification = NULL, solver = "auto") {
  version <- .require_adapter_version("fmridesign", "0.6.0")
  study <- .validate_study(study)
  if (!inherits(model, "event_model")) {
    stop("`model` must inherit from `fmridesign::event_model()`.",
      call. = FALSE)
  }
  block_map <- .fmridesign_block_map(study$partitions, model, block_map)
  .fmridesign_assert_event_binding(model, study, block_map)
  coordinates <- unname(fmridesign::conditions(model))
  conditions <- condition_space(
    coordinates, basis_id = basis_id, units = units, scale = scale,
    provenance = list(compiler = "fmridesign", version = version)
  )
  if (is.matrix(semantic_map)) {
    semantic_map <- rep(list(semantic_map), length(study$partitions))
    names(semantic_map) <- study$partitions
  }
  if (!is.null(semantic_map) &&
      (!is.list(semantic_map) || !setequal(names(semantic_map),
        study$partitions))) {
    stop("`semantic_map` must be one matrix or a named matrix per partition.",
      call. = FALSE)
  }
  designs <- parameterizations <- stats::setNames(
    vector("list", length(study$partitions)), study$partitions
  )
  for (partition in study$partitions) {
    compiled <- fmridesign::design_matrix(
      model, blockid = block_map[[partition]]
    )
    explicit <- if (is.null(semantic_map)) NULL else semantic_map[[partition]]
    map <- .fmridesign_semantic_map(compiled, conditions, explicit)
    design <- as.matrix(compiled)
    storage.mode(design) <- "double"
    rownames(design) <- as.character(
      study$observations$indexes[[partition]]$observation_id
    )
    designs[[partition]] <- design
    parameterizations[[partition]] <- coefficient_parameterization(
      map, conditions,
      coding_id = paste0("fmridesign-", version, "-", partition),
      provenance = list(
        compiler = "fmridesign",
        version = version,
        block = block_map[[partition]]
      )
    )
  }
  if (is.null(specification)) {
    specification <- .fmridesign_semantic_specification(model)
  }
  design_model(
    specification = specification,
    conditions = conditions,
    designs = designs,
    parameterizations = parameterizations,
    row_ids = lapply(study$observations$indexes, `[[`, "observation_id"),
    solver = solver,
    protocol = "first-moment-relation-v1/fmridesign",
    protocol_version = "1",
    package = "fmridesign",
    package_version = version,
    provenance = list(
      # Frozen because semantic model provenance participates in plan identity.
      adapter = "crossform::fmridesign_design_model",
      block_map = as.list(block_map)
    )
  )
}
