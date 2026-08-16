# Question-first generated example -----------------------------------------

.with_example_seed <- function(seed, code) {
  if (!is.numeric(seed) || length(seed) != 1L || is.na(seed) ||
      !is.finite(seed) || seed %% 1 != 0 || seed < 0 ||
      seed > .Machine$integer.max) {
    stop("`seed` must be one nonnegative integer.", call. = FALSE)
  }
  global <- .GlobalEnv
  had_seed <- exists(".Random.seed", envir = global, inherits = FALSE)
  if (had_seed) old_seed <- get(".Random.seed", envir = global, inherits = FALSE)
  on.exit({
    if (had_seed) {
      assign(".Random.seed", old_seed, envir = global)
    } else if (exists(".Random.seed", envir = global, inherits = FALSE)) {
      rm(".Random.seed", envir = global)
    }
  }, add = TRUE)
  set.seed(as.integer(seed))
  force(code)
}

.example_count <- function(x, name, minimum) {
  if (!is.numeric(x) || length(x) != 1L || is.na(x) || !is.finite(x) ||
      x %% 1 != 0 || x < minimum || x > .Machine$integer.max) {
    stop(sprintf("`%s` must be one integer at least %d.", name, minimum),
      call. = FALSE)
  }
  as.integer(x)
}

#' Generate a small task-fMRI relation with known spatial truth
#'
#' This is the executable newcomer fixture. Four conditions are observed in
#' independent runs over a small volume. A reproducible animate-versus-
#' inanimate multivariate pattern is planted near the volume centre, and raw
#' trial responses are retained through [lm_relation_fit()] so both point
#' geometry and the admitted fixed-metric uncertainty path are available.
#'
#' The returned object is generated, not empirical data. Use the Haxby exemplar
#' under `exemplars/haxby2001` for the public-data parity workflow.
#'
#' @param seed Nonnegative integer random seed. The caller's random-number state
#'   is restored on exit.
#' @param dimensions Three volume dimensions, each at least three.
#' @param partitions Number of independent runs, at least two.
#' @param trials_per_condition Trials per condition and run, at least two.
#' @param noise_sd Positive residual standard deviation.
#' @param spacing Three positive voxel spacings in millimetres.
#' @param searchlight_radius Positive searchlight radius in millimetres.
#' @return A named list containing `fit`, its full `domain`, a compiled
#'   searchlight `frame`, the `contrast` weights, a category `model_rdm`, and
#'   `truth` metadata identifying the planted features and noiseless patterns.
#' @examples
#' example <- example_fmri_effects()
#' plan <- plan_geometry(
#'   example$fit$relation, example$frame,
#'   cross_partitions(example$fit$relation, independence = "independent")
#' )
#' distances <- rdm(plan)
#' dim(distances$values)
#' @export
example_fmri_effects <- function(
    seed = 20260814L,
    dimensions = c(7L, 7L, 5L),
    partitions = 4L,
    trials_per_condition = 8L,
    noise_sd = 0.6,
    spacing = c(3, 3, 3),
    searchlight_radius = 4) {
  if (!is.numeric(dimensions) || length(dimensions) != 3L ||
      anyNA(dimensions) || any(!is.finite(dimensions)) ||
      any(dimensions %% 1 != 0) || any(dimensions < 3) ||
      any(dimensions > .Machine$integer.max)) {
    stop("`dimensions` must contain three integers of at least three.",
      call. = FALSE)
  }
  dimensions <- as.integer(dimensions)
  partitions <- .example_count(partitions, "partitions", 2L)
  trials_per_condition <- .example_count(
    trials_per_condition, "trials_per_condition", 2L
  )
  if (!is.numeric(noise_sd) || length(noise_sd) != 1L || is.na(noise_sd) ||
      !is.finite(noise_sd) || noise_sd <= 0) {
    stop("`noise_sd` must be one positive finite number.", call. = FALSE)
  }
  if (!is.numeric(spacing) || length(spacing) != 3L || anyNA(spacing) ||
      any(!is.finite(spacing)) || any(spacing <= 0)) {
    stop("`spacing` must contain three positive finite values.", call. = FALSE)
  }
  if (!is.numeric(searchlight_radius) || length(searchlight_radius) != 1L ||
      is.na(searchlight_radius) || !is.finite(searchlight_radius) ||
      searchlight_radius <= 0) {
    stop("`searchlight_radius` must be one positive finite number.",
      call. = FALSE)
  }

  .with_example_seed(seed, {
    conditions <- c("face", "body", "house", "tool")
    mask <- array(TRUE, dimensions)
    domain <- volume_domain(
      # Frozen dataset identity: changing this branding token would change
      # every plan built from the otherwise unchanged generated fixture.
      mask, spacing = spacing, id = "effectagram:generated-fmri-example"
    )
    centre <- colMeans(domain$coordinates)
    displacement <- sweep(domain$coordinates, 2L, centre, `-`)
    distance <- sqrt(rowSums(displacement^2))
    planted <- which(distance <= 1.5 * min(spacing))
    profile <- numeric(domain$n_features)
    local <- displacement[planted, 1L] - 0.7 * displacement[planted, 2L] +
      0.4 * displacement[planted, 3L]
    if (all(local == 0)) local <- seq_along(planted) - mean(seq_along(planted))
    local <- local / sqrt(mean(local^2))
    profile[planted] <- local

    contrast_weights <- stats::setNames(
      c(0.5, 0.5, -0.5, -0.5), conditions
    )
    true_patterns <- outer(2 * contrast_weights, profile)
    dimnames(true_patterns) <- list(conditions, domain$feature_ids)
    condition <- factor(
      rep(conditions, each = trials_per_condition), levels = conditions
    )
    design <- stats::model.matrix(~ 0 + condition)
    colnames(design) <- conditions
    targets <- diag(length(conditions))
    dimnames(targets) <- list(conditions, conditions)
    sources <- stats::setNames(lapply(seq_len(partitions), function(run) {
      design %*% true_patterns + matrix(
        stats::rnorm(nrow(design) * domain$n_features, sd = noise_sd),
        nrow(design), domain$n_features
      )
    }), paste0("run", seq_len(partitions)))
    effects <- effect_space(
      conditions, basis_id = "generated-condition-means:v1",
      units = "arbitrary-BOLD"
    )
    fit <- lm_relation_fit(
      sources, design, targets, effect_names = effects,
      domain = domain, sampling_unit = "trial",
      provenance = list(fixture = "example_fmri_effects", seed = as.integer(seed))
    )
    frame <- compile_frame(
      searchlights(searchlight_radius, normalization = "local"), domain
    )
    signal_measurements <- which(
      Matrix::rowSums(frame$weights[, planted, drop = FALSE]) > 0
    )
    model_rdm <- outer(
      sign(contrast_weights), sign(contrast_weights),
      function(left, right) as.numeric(left != right)
    )
    dimnames(model_rdm) <- list(conditions, conditions)
    list(
      fit = fit,
      domain = domain,
      frame = frame,
      contrast = contrast_weights,
      model_rdm = model_rdm,
      truth = list(
        planted_features = planted,
        planted_feature_ids = domain$feature_ids[planted],
        signal_measurements = signal_measurements,
        contrast_pattern = profile,
        condition_patterns = true_patterns,
        noise_sd = noise_sd,
        seed = as.integer(seed)
      )
    )
  })
}
