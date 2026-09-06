# Fixed-total coherent/configuration mixture generator.
#
# This file is simulation infrastructure, not package code. It uses only base R
# and returns truth separately from any Crossform estimator output.

.mixture_scalar <- function(value, name, lower = -Inf, upper = Inf,
                            open_lower = FALSE) {
  valid <- is.numeric(value) && length(value) == 1L && is.finite(value) &&
    value <= upper && if (open_lower) value > lower else value >= lower
  if (!valid) {
    interval <- if (open_lower) "(" else "["
    stop(name, " must be one finite number in ", interval, lower, ", ",
         upper, "].", call. = FALSE)
  }
  as.numeric(value)
}

.mixture_direction <- function(value, n_features, name) {
  if (!is.numeric(value) || length(value) != n_features ||
      any(!is.finite(value))) {
    stop(name, " must contain exactly ", n_features,
         " finite numeric values.", call. = FALSE)
  }
  as.numeric(value)
}

.mixture_seeded_configuration <- function(n_features, seed, tolerance) {
  old_exists <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (old_exists) old_seed <- get(".Random.seed", envir = .GlobalEnv)
  on.exit({
    if (old_exists) {
      assign(".Random.seed", old_seed, envir = .GlobalEnv)
    } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  }, add = TRUE)
  set.seed(seed)
  value <- stats::rnorm(n_features)
  value <- value - mean(value)
  if (sqrt(sum(value^2)) <= tolerance) {
    value <- rep(c(-1, 1), length.out = n_features)
    value <- value - mean(value)
  }
  value / sqrt(sum(value^2))
}

fixed_total_mixture <- function(n_features = 16L, total_magnitude = 1,
                                theta = pi / 4, seed = 20260821L,
                                mean_direction = NULL,
                                configuration_direction = NULL,
                                allow_nonorthogonal = FALSE,
                                allow_unnormalized = FALSE,
                                allow_noncanonical_mean = FALSE,
                                tolerance = 1e-10) {
  if (!is.numeric(n_features) || length(n_features) != 1L ||
      !is.finite(n_features) || n_features != as.integer(n_features) ||
      n_features < 2L) {
    stop("n_features must be one whole number of at least 2.", call. = FALSE)
  }
  n_features <- as.integer(n_features)
  total_magnitude <- .mixture_scalar(
    total_magnitude, "total_magnitude", 0, Inf, open_lower = TRUE
  )
  theta <- .mixture_scalar(theta, "theta", 0, pi / 2)
  if (!is.numeric(seed) || length(seed) != 1L || !is.finite(seed) ||
      seed != as.integer(seed)) {
    stop("seed must be one finite whole number.", call. = FALSE)
  }
  seed <- as.integer(seed)
  tolerance <- .mixture_scalar(tolerance, "tolerance", 0, Inf,
                               open_lower = TRUE)

  canonical_mean <- rep(1 / sqrt(n_features), n_features)
  if (is.null(mean_direction)) {
    mean_direction <- canonical_mean
  } else {
    mean_direction <- .mixture_direction(
      mean_direction, n_features, "mean_direction"
    )
  }
  if (is.null(configuration_direction)) {
    configuration_direction <- .mixture_seeded_configuration(
      n_features, seed, tolerance
    )
  } else {
    configuration_direction <- .mixture_direction(
      configuration_direction, n_features, "configuration_direction"
    )
  }

  mean_norm <- sqrt(sum(mean_direction^2))
  configuration_norm <- sqrt(sum(configuration_direction^2))
  inner_product <- sum(mean_direction * configuration_direction)
  canonical_gap <- max(abs(mean_direction - canonical_mean))
  configuration_mean <- mean(configuration_direction)

  if (!allow_unnormalized &&
      (abs(mean_norm - 1) > tolerance ||
       abs(configuration_norm - 1) > tolerance)) {
    stop("mean and configuration directions must have Euclidean norm 1; ",
         "set allow_unnormalized = TRUE only for an explicit negative ",
         "fixture.", call. = FALSE)
  }
  if (!allow_nonorthogonal && abs(inner_product) > tolerance) {
    stop("mean and configuration directions must be orthogonal; set ",
         "allow_nonorthogonal = TRUE only for an explicit negative fixture.",
         call. = FALSE)
  }
  if (!allow_noncanonical_mean &&
      (canonical_gap > tolerance || abs(configuration_mean) > tolerance)) {
    stop("mean_direction must be the normalized constant mode and ",
         "configuration_direction must have zero mean; set ",
         "allow_noncanonical_mean = TRUE only for an explicit negative ",
         "fixture.", call. = FALSE)
  }

  pattern <- sqrt(total_magnitude) * (
    cos(theta) * mean_direction + sin(theta) * configuration_direction
  )
  coherent_projector <- tcrossprod(canonical_mean)
  configuration_projector <- diag(n_features) - coherent_projector
  expected <- c(
    total = sum(pattern^2),
    coherent = drop(crossprod(pattern, coherent_projector %*% pattern)),
    configuration = drop(crossprod(
      pattern, configuration_projector %*% pattern
    ))
  )
  expected <- c(
    expected,
    coherent_share = expected[["coherent"]] / expected[["total"]],
    configuration_share = expected[["configuration"]] / expected[["total"]]
  )

  effects <- rbind(condition_a = pattern / 2, condition_b = -pattern / 2)
  contrast <- c(condition_a = 1, condition_b = -1)
  basis_id <- paste0(
    "fixed-total-mixture-v1:p", n_features, ":seed", seed,
    ":theta", format(theta, digits = 17, scientific = FALSE)
  )
  strict_basis <- !allow_nonorthogonal && !allow_unnormalized &&
    !allow_noncanonical_mean

  structure(list(
    effect_pattern = pattern,
    effect_matrix = effects,
    contrast = contrast,
    truth = list(
      total = unname(expected[["total"]]),
      coherent = unname(expected[["coherent"]]),
      configuration = unname(expected[["configuration"]]),
      coherent_share = unname(expected[["coherent_share"]]),
      configuration_share = unname(expected[["configuration_share"]]),
      component_spectrum = unname(expected[c("coherent", "configuration")]),
      neural_spectrum = c(unname(expected[["total"]]),
                          rep(0, n_features - 1L))
    ),
    metadata = list(
      schema_version = "matched-mixture-v1",
      basis_id = basis_id,
      basis = list(
        mean = mean_direction,
        configuration = configuration_direction,
        gram = crossprod(cbind(mean_direction, configuration_direction)),
        canonical_mean_gap = canonical_gap,
        configuration_mean = configuration_mean
      ),
      normalization = list(
        direction = "euclidean_unit",
        frame = "none",
        coherent_projector = "11T/n_features"
      ),
      total_magnitude = total_magnitude,
      theta = theta,
      theta_degrees = theta * 180 / pi,
      seed = seed,
      expected_spectrum = unname(expected[c(
        "total", "coherent", "configuration",
        "coherent_share", "configuration_share"
      )]),
      strict_basis = strict_basis,
      explicit_relaxations = c(
        nonorthogonal = allow_nonorthogonal,
        unnormalized = allow_unnormalized,
        noncanonical_mean = allow_noncanonical_mean
      )
    )
  ), class = "crossform_matched_mixture")
}

fixed_total_mixture_grid <- function(theta = c(0, pi / 6, pi / 4,
                                                pi / 3, pi / 2), ...) {
  if (!is.numeric(theta) || !length(theta) || any(!is.finite(theta))) {
    stop("theta must contain finite numeric mixture angles.", call. = FALSE)
  }
  lapply(theta, function(value) fixed_total_mixture(theta = value, ...))
}
