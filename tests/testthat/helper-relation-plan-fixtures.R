relation_plan_inverse <- function(value, tolerance = 1e-12) {
  decomposition <- svd(value, nu = min(dim(value)), nv = ncol(value))
  keep <- decomposition$d > tolerance * max(decomposition$d)
  decomposition$v[, keep, drop = FALSE] %*%
    (t(decomposition$u[, keep, drop = FALSE]) / decomposition$d[keep])
}

relation_plan_fixture <- function(solver = "qr", coding = "cell",
                                  observation_kind = "ols",
                                  use_functions = FALSE) {
  bound <- bound_study_fixture(use_functions)
  conditions <- condition_space(
    c("face", "place", "object"),
    basis_id = "canonical-hrf-amplitude",
    units = "percent-signal-change"
  )
  cell_map <- cbind(
    face = c(1, 0, 0),
    place = c(0, 1, 0),
    object = c(0, 0, 1),
    drift = c(0, 0, 0)
  )
  rownames(cell_map) <- conditions$coordinates
  treatment_map <- cbind(
    intercept = c(1, 1, 1),
    place_minus_face = c(0, 1, 0),
    object_minus_face = c(0, 0, 1),
    drift = c(0, 0, 0)
  )
  rownames(treatment_map) <- conditions$coordinates

  designs <- parameterizations <- list()
  coding_matrix <- rbind(
    c(1, 0, 0, 0),
    c(1, 1, 0, 0),
    c(1, 0, 1, 0),
    c(0, 0, 0, 1)
  )
  for (partition in bound$fixture$partitions) {
    count <- length(bound$fixture$indexes[[partition]]$observation_id)
    condition <- rep(conditions$coordinates, length.out = count)
    semantic <- stats::model.matrix(~ condition - 1)
    colnames(semantic) <- sub("condition", "", colnames(semantic))
    semantic <- semantic[, conditions$coordinates, drop = FALSE]
    cell <- cbind(semantic, drift = seq(-0.91, 1.07, length.out = count))
    rownames(cell) <- as.character(
      bound$fixture$indexes[[partition]]$observation_id
    )
    if (identical(coding, "cell")) {
      designs[[partition]] <- cell
      parameterizations[[partition]] <- coefficient_parameterization(
        cell_map, conditions, coding_id = "cell-means"
      )
    } else {
      treatment <- cell %*% coding_matrix
      colnames(treatment) <- colnames(treatment_map)
      rownames(treatment) <- rownames(cell)
      designs[[partition]] <- treatment
      parameterizations[[partition]] <- coefficient_parameterization(
        treatment_map, conditions, coding_id = "treatment"
      )
    }
  }
  model <- design_model(
    specification = list(
      mean = "canonical condition amplitudes",
      nuisance = "linear drift"
    ),
    conditions = conditions,
    designs = designs,
    parameterizations = parameterizations,
    row_ids = lapply(bound$fixture$indexes, `[[`, "observation_id"),
    solver = solver,
    package_version = if (identical(coding, "cell")) "1.0" else "2.0"
  )
  weights <- rbind(
    place_minus_object = c(0, 1, -1),
    face_vs_others = c(1, -0.5, -0.5)
  )
  colnames(weights) <- conditions$coordinates
  effects <- effect_map(weights, conditions)

  whiteners <- stats::setNames(lapply(bound$fixture$counts, function(count) {
    diag(seq(0.81, 1.19, length.out = count))
  }), bound$fixture$partitions)
  observation <- switch(observation_kind,
    ols = observation_model("ols", sampling_unit = "scan"),
    fixed_gls = observation_model(
      "fixed_gls", sampling_unit = "scan", whitener = whiteners
    ),
    learned_frozen_gls = observation_model(
      "learned_frozen_gls", sampling_unit = "scan", whitener = whiteners,
      training_revision = paste0(
        "sha256:", paste(rep("d", 64), collapse = "")
      ),
      training_provenance = list(method = "AR1", source = "training-split")
    )
  )
  list(
    bound = bound,
    conditions = conditions,
    model = model,
    effects = effects,
    observation = observation,
    plan = plan_relation(
      study(
        bound$observations,
        bound$events,
        bound$confounds,
        bound$hierarchy
      ),
      model,
      effects,
      observation
    )
  )
}
