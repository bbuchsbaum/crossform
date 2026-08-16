first_moment_vertical_fixture <- function(
    coding = c("cell", "treatment"), solver = "qr",
    observation_kind = c("fixed_gls", "learned_frozen_gls", "ols"),
    seed = 2026081517L) {
  coding <- match.arg(coding)
  observation_kind <- match.arg(observation_kind)
  set.seed(seed)
  partitions <- paste0("run-", 1:4)
  counts <- c(28L, 30L, 32L, 34L)
  retained_count <- 24L
  condition_names <- c("face", "body", "house", "tool")
  domain <- abstract_domain(60L, id = "first-moment-vertical-slice:v1")

  indexes <- sources <- designs_cell <- whiteners <- vector(
    "list", length(partitions)
  )
  names(indexes) <- names(sources) <- names(designs_cell) <-
    names(whiteners) <- partitions
  event_tables <- confound_tables <- vector("list", length(partitions))
  truth <- matrix(rnorm(length(condition_names) * domain$n_features, sd = 0.35),
    length(condition_names), domain$n_features,
    dimnames = list(condition_names, domain$feature_ids))
  nuisance_truth <- matrix(rnorm(2L * domain$n_features, sd = 0.08), 2L,
    domain$n_features)

  for (index in seq_along(partitions)) {
    partition <- partitions[[index]]
    count <- counts[[index]]
    observation_id <- seq_len(count)
    time <- seq(0, by = 2, length.out = count)
    condition <- factor(rep(condition_names, length.out = count),
      levels = condition_names)
    indicator <- stats::model.matrix(~ condition - 1)
    colnames(indicator) <- condition_names
    drift <- rep(seq(-1, 1, length.out = retained_count),
      length.out = count)
    motion <- rep(sin(seq(0, 2 * pi, length.out = retained_count)),
      length.out = count)
    design <- cbind(indicator, drift = drift, motion = motion)
    rownames(design) <- as.character(observation_id)
    response <- design[, condition_names, drop = FALSE] %*% truth +
      design[, c("drift", "motion"), drop = FALSE] %*% nuisance_truth +
      matrix(rnorm(count * domain$n_features, sd = 0.45), count,
        domain$n_features)

    retained_whitener <- diag(seq(0.85, 1.15, length.out = retained_count))
    whitener <- diag(count)
    whitener[seq_len(retained_count), seq_len(retained_count)] <-
      retained_whitener
    indexes[[partition]] <- observation_index(
      observation_id, partition, time = time, units = "seconds"
    )
    sources[[partition]] <- response
    designs_cell[[partition]] <- design
    whiteners[[partition]] <- whitener
    event_tables[[index]] <- data.frame(
      partition = partition,
      event_id = paste0(partition, "-scan-", observation_id),
      onset = time,
      duration = 0,
      condition = as.character(condition),
      stringsAsFactors = FALSE
    )
    confound_tables[[index]] <- data.frame(
      partition = partition,
      observation_id = observation_id,
      motion = motion,
      retained = observation_id <= retained_count,
      stringsAsFactors = FALSE
    )
  }

  observation_record <- observations(sources, indexes, domain)
  event_record <- observation_events(do.call(rbind, event_tables))
  confound_record <- observation_confounds(
    do.call(rbind, confound_tables), censor = "retained"
  )
  hierarchy <- partition_hierarchy(data.frame(
    partition = partitions,
    run = partitions,
    session = rep(c("session-1", "session-2"), each = 2L),
    subject = "subject-1",
    stringsAsFactors = FALSE
  ))
  study_value <- study(
    observation_record, event_record, confound_record, hierarchy,
    provenance = list(fixture = "first-moment-vertical-slice:v1")
  )
  conditions <- condition_space(
    condition_names,
    basis_id = "scan-level-condition-mean:v1",
    units = "arbitrary-BOLD"
  )
  cell_map <- cbind(diag(length(condition_names)), drift = 0, motion = 0)
  rownames(cell_map) <- condition_names
  colnames(cell_map)[seq_along(condition_names)] <- condition_names

  treatment_map <- cbind(
    intercept = 1,
    body_minus_face = c(0, 1, 0, 0),
    house_minus_face = c(0, 0, 1, 0),
    tool_minus_face = c(0, 0, 0, 1),
    drift = 0,
    motion = 0
  )
  rownames(treatment_map) <- condition_names
  coding_matrix <- rbind(
    c(1, 0, 0, 0, 0, 0),
    c(1, 1, 0, 0, 0, 0),
    c(1, 0, 1, 0, 0, 0),
    c(1, 0, 0, 1, 0, 0),
    c(0, 0, 0, 0, 1, 0),
    c(0, 0, 0, 0, 0, 1)
  )

  designs <- parameterizations <- vector("list", length(partitions))
  names(designs) <- names(parameterizations) <- partitions
  for (partition in partitions) {
    if (identical(coding, "cell")) {
      designs[[partition]] <- designs_cell[[partition]]
      parameterizations[[partition]] <- coefficient_parameterization(
        cell_map, conditions, coding_id = "cell-means-plus-nuisance"
      )
    } else {
      designs[[partition]] <- designs_cell[[partition]] %*% coding_matrix
      colnames(designs[[partition]]) <- colnames(treatment_map)
      rownames(designs[[partition]]) <- rownames(designs_cell[[partition]])
      parameterizations[[partition]] <- coefficient_parameterization(
        treatment_map, conditions, coding_id = "treatment-plus-nuisance"
      )
    }
  }
  model <- design_model(
    specification = list(
      target = "scan-level condition means",
      nuisance = c("linear drift", "declared motion covariate"),
      censor_policy = "retain first 24 observations per partition"
    ),
    conditions = conditions,
    designs = designs,
    parameterizations = parameterizations,
    row_ids = lapply(indexes, `[[`, "observation_id"),
    solver = solver,
    protocol = "first-moment-vertical-slice",
    protocol_version = "1",
    package_version = "1"
  )
  identity_effects <- diag(length(condition_names))
  dimnames(identity_effects) <- list(condition_names, condition_names)
  effects <- effect_map(identity_effects, conditions)
  independence <- "partitions independently acquired conditional on model"
  observation <- switch(observation_kind,
    ols = observation_model(
      "ols", sampling_unit = "scan", independence = independence
    ),
    fixed_gls = observation_model(
      "fixed_gls", sampling_unit = "scan", whitener = whiteners,
      independence = independence
    ),
    learned_frozen_gls = observation_model(
      "learned_frozen_gls", sampling_unit = "scan", whitener = whiteners,
      independence = independence,
      training_revision = paste0("sha256:", paste(rep("e", 64), collapse = "")),
      training_provenance = list(
        method = "fixture-frozen-whitener",
        training_partition = "independent-training-data"
      )
    )
  )
  plan <- plan_relation(study_value, model, effects, observation)
  frame <- compile_frame(
    regions(rep(paste0("roi-", 1:5), each = 12L), normalization = "local"),
    domain
  )
  over <- cross_partitions(
    partitions,
    independence = "independent",
    generalizes_over = "run"
  )
  list(
    version = "first-moment-vertical-slice:v1",
    partitions = partitions,
    counts = counts,
    retained_count = retained_count,
    sources = sources,
    truth = truth,
    study = study_value,
    conditions = conditions,
    model = model,
    effects = effects,
    observation = observation,
    plan = plan,
    frame = frame,
    over = over
  )
}

first_moment_direct_blocks <- function(fixture, tolerance = 1e-12) {
  stats::setNames(lapply(fixture$plan$partitions, function(partition) {
    receipt <- fixture$plan$design_receipts[[partition]]
    rows <- fixture$plan$retained_rows[[partition]]
    response <- fixture$sources[[partition]][rows, , drop = FALSE]
    whitener <- fixture$plan$whiteners[[partition]]
    extractor <- receipt$lowered_target %*%
      relation_plan_inverse(whitener %*% receipt$design, tolerance) %*%
      whitener
    value <- extractor %*% response
    colnames(value) <- NULL
    value
  }), fixture$plan$partitions)
}

first_moment_direct_forms <- function(blocks, frame, over) {
  pairs <- utils::combn(seq_along(blocks), 2L)
  lapply(seq_len(nrow(frame$weights)), function(node) {
    weight <- as.numeric(frame$weights[node, ])
    Reduce(`+`, lapply(seq_len(ncol(pairs)), function(edge) {
      left <- blocks[[pairs[1L, edge]]]
      right <- blocks[[pairs[2L, edge]]]
      cross <- (left * rep(weight, each = nrow(left))) %*% t(right)
      0.5 * (cross + t(cross))
    })) / ncol(pairs)
  })
}

first_moment_direct_rdm <- function(forms) {
  do.call(rbind, lapply(forms, function(form) {
    pairs <- utils::combn(seq_len(nrow(form)), 2L)
    vapply(seq_len(ncol(pairs)), function(edge) {
      left <- pairs[1L, edge]
      right <- pairs[2L, edge]
      form[left, left] + form[right, right] - 2 * form[left, right]
    }, numeric(1))
  }))
}
