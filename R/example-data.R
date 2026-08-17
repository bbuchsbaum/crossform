# Question-first generated example -----------------------------------------

.with_example_seed <- function(seed, code) {
  if (!.is_number(seed) || seed %% 1 != 0 || seed < 0 ||
      seed > .Machine$integer.max) {
    .input_error("`seed` must be one nonnegative integer.")
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
  if (!.is_number(x) || x %% 1 != 0 || x < minimum || x > .Machine$integer.max) {
    .input_error(
      sprintf("`%s` must be one integer at least %d.", name, minimum)
    )
  }
  as.integer(x)
}

#' Generate a small task-fMRI relation with known spatial truth
#'
#' This is the executable newcomer fixture. Four conditions are observed in
#' independent runs over a small volume, and raw trial responses are retained
#' through [lm_relation_fit()] so both point geometry and the admitted
#' fixed-metric uncertainty path are available.
#'
#' Two animate-versus-inanimate blocks are planted, at opposite ends of the
#' longest axis and never overlapping. Both carry the same per-voxel amplitude;
#' only the sign structure differs, so the two halves of the energy
#' decomposition are each demonstrated by one block:
#'
#' - the **pattern** block alternates sign between neighboring voxels within
#'   each axial slice, so its frame-weighted average nearly cancels and its
#'   reproducible energy is almost entirely `configuration`;
#' - the **mean** block shifts every voxel the same way, so its energy is
#'   almost entirely `coherent`.
#'
#' The returned object is generated, not empirical data. Use the Haxby exemplar
#' under `exemplars/haxby2001` for the public-data parity workflow.
#'
#' @param seed Nonnegative integer random seed. The caller's random-number state
#'   is restored on exit.
#' @param dimensions Three volume dimensions, each at least three. The default
#'   is the smallest volume in which the two planted blocks, and the
#'   searchlights that touch them, stay disjoint; smaller volumes still work
#'   but the blocks crowd each other.
#' @param partitions Number of independent runs, at least two.
#' @param trials_per_condition Trials per condition and run, at least two.
#' @param noise_sd Positive residual standard deviation.
#' @param spacing Three positive voxel spacings in millimeters.
#' @param searchlight_radius Positive searchlight radius in millimeters.
#' @param plant Which blocks to plant: `"pattern"`, `"mean"`, or both (the
#'   default). Dropping one leaves its feature and measurement sets empty; the
#'   noise draw is unchanged either way, so the two settings differ only in the
#'   planted signal.
#' @return An `effect_example_effects` list.
#' @section Structure:
#' The returned list has six public elements.
#'
#' - `$fit`: an [lm_relation_fit()] carrying a residual channel, so
#'   [rdm_sampling_covariance()] is admitted. `$fit$relation` is what
#'   [plan_geometry()] takes.
#' - `$domain`: the full [volume_domain()] the fit was made over.
#' - `$frame`: a searchlight [compile_frame()] over that domain, one
#'   measurement per voxel.
#' - `$contrast`: named animate-versus-inanimate weights over the four
#'   conditions.
#' - `$model_rdm`: a condition-by-condition category model for [rsa()].
#' - `$truth`: what was planted, listed below.
#'
#' `$truth` holds feature indices into the domain, measurement indices into the
#' frame, and the generating settings.
#'
#' - `$planted_features`, `$mean_features`: the two disjoint voxel sets, as
#'   positions in the domain. `$planted_feature_ids` gives the same voxels as
#'   domain feature identifiers.
#' - `$pattern_measurements`, `$mean_measurements`: the searchlights that
#'   overlap each block. `$signal_measurements` is their sorted union, which is
#'   what a map-reading example should highlight.
#' - `$contrast_pattern`: the planted per-voxel contrast profile over the whole
#'   domain, zero outside the two blocks.
#' - `$condition_patterns`: the noiseless condition-by-voxel means.
#' - `$noise_sd`, `$seed`: the generating settings.
#' @family geometry plans and views
#' @seealso [plan_geometry()] for the next step, then [contrast_energy()],
#'   [rdm()], or [rsa()]; and [rdm_sampling_covariance()], which the retained
#'   residual channel makes available.
#' @examples
#' example <- example_fmri_effects()
#'
#' # The planted contrast and the two blocks it was planted in are carried
#' # alongside the data, so any result can be checked against ground truth.
#' example$contrast
#' c(pattern = length(example$truth$planted_features),
#'   mean = length(example$truth$mean_features))
#'
#' # Everything a second-moment question needs is already built: a fit with
#' # residuals, a compiled searchlight frame, and the cross-run pairing.
#' plan <- plan_geometry(
#'   example$fit$relation, example$frame,
#'   cross_partitions(example$fit$relation, independence = "independent")
#' )
#' distances <- rdm(plan)
#' dim(distances$values)
#'
#' # Each block reproduces the half of the decomposition it was built from.
#' energy <- contrast_energy(plan, example$contrast)
#' round(c(
#'   pattern_configuration = max(energy$configuration[
#'     example$truth$pattern_measurements]),
#'   mean_coherent = max(energy$coherent[example$truth$mean_measurements]),
#'   elsewhere = max(energy$total[-example$truth$signal_measurements])
#' ), 3)
#' @export
example_fmri_effects <- function(
    seed = 20260814L,
    dimensions = c(8L, 7L, 5L),
    partitions = 4L,
    trials_per_condition = 8L,
    noise_sd = 0.6,
    spacing = c(3, 3, 3),
    searchlight_radius = 4,
    plant = c("pattern", "mean")) {
  if (!.is_finite_numeric(dimensions) || length(dimensions) != 3L ||
      anyNA(dimensions) || any(dimensions %% 1 != 0) || any(dimensions < 3) ||
      any(dimensions > .Machine$integer.max)) {
    .input_error("`dimensions` must contain three integers of at least three.")
  }
  dimensions <- as.integer(dimensions)
  partitions <- .example_count(partitions, "partitions", 2L)
  trials_per_condition <- .example_count(
    trials_per_condition, "trials_per_condition", 2L
  )
  .check_number(noise_sd, "noise_sd", positive = TRUE)
  if (!.is_finite_numeric(spacing) || length(spacing) != 3L ||
      anyNA(spacing) || any(spacing <= 0)) {
    .input_error("`spacing` must contain three positive finite values.")
  }
  .check_number(searchlight_radius, "searchlight_radius", positive = TRUE)
  plant <- match.arg(plant, c("pattern", "mean"), several.ok = TRUE)

  .with_example_seed(seed, {
    conditions <- c("face", "body", "house", "tool")
    mask <- array(TRUE, dimensions)
    domain <- volume_domain(
      # Frozen dataset identity: changing this branding token would change
      # every plan built from the otherwise unchanged generated fixture.
      mask, spacing = spacing, id = "crossform:generated-fmri-example"
    )
    # The two blocks sit one voxel in from either end of the longest axis and
    # are centered on the other two. At the default dimensions that keeps the
    # blocks, and the searchlights touching them, disjoint.
    center <- colMeans(domain$coordinates)
    long_axis <- which.max(dimensions)
    axis_levels <- sort(unique(domain$coordinates[, long_axis]))
    block_center <- function(level) {
      value <- center
      value[[long_axis]] <- level
      value
    }
    in_ball <- function(origin) {
      displacement <- sweep(domain$coordinates, 2L, origin, `-`)
      which(sqrt(rowSums(displacement^2)) <= 1.5 * min(spacing))
    }
    pattern_features <- in_ball(
      block_center(axis_levels[[length(axis_levels) - 1L]])
    )
    mean_features <- setdiff(in_ball(block_center(axis_levels[[2L]])),
      pattern_features)
    if (!"pattern" %in% plant) pattern_features <- integer(0)
    if (!"mean" %in% plant) mean_features <- integer(0)

    # Voxel grid positions, used only for the alternating sign pattern.
    grid <- round(sweep(domain$coordinates, 2L, spacing, `/`))
    profile <- numeric(domain$n_features)
    # Alternating within each axial slice: the six-neighbor searchlight
    # average of this profile is -1/7 of the center value, so almost none of
    # the block's reproducible energy survives into the coherent component.
    profile[pattern_features] <- (-1)^(grid[pattern_features, 1L] +
      grid[pattern_features, 2L])
    # The same amplitude with one sign everywhere: a pure regional mean shift.
    profile[mean_features] <- 1

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
    touching <- function(features) {
      if (!length(features)) {
        return(integer(0))
      }
      which(Matrix::rowSums(frame$weights[, features, drop = FALSE]) > 0)
    }
    pattern_measurements <- touching(pattern_features)
    mean_measurements <- touching(mean_features)
    signal_measurements <- sort(
      union(pattern_measurements, mean_measurements)
    )
    model_rdm <- outer(
      sign(contrast_weights), sign(contrast_weights),
      function(left, right) as.numeric(left != right)
    )
    dimnames(model_rdm) <- list(conditions, conditions)
    structure(list(
      fit = fit,
      domain = domain,
      frame = frame,
      contrast = contrast_weights,
      model_rdm = model_rdm,
      truth = list(
        planted_features = pattern_features,
        planted_feature_ids = domain$feature_ids[pattern_features],
        mean_features = mean_features,
        mean_feature_ids = domain$feature_ids[mean_features],
        pattern_measurements = pattern_measurements,
        mean_measurements = mean_measurements,
        signal_measurements = signal_measurements,
        contrast_pattern = profile,
        condition_patterns = true_patterns,
        noise_sd = noise_sd,
        seed = as.integer(seed)
      )
    ), class = c("effect_example_effects", "list"))
  })
}
