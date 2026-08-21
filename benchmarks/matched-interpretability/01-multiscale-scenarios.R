# Conservative multiscale scenarios for the fixed-total mixture generator.
# Source 00-mixture-generator.R before this file.

matched_multiscale_scenarios <- function(
    n_features = 17L,
    radii = c(0.01, 1.01, 2.01, 4.01),
    alpha = rep(1 / length(radii), length(radii)),
    total_magnitude = 4,
    seed = 20260822L,
    id = "matched-multiscale-v1") {
  if (!exists("fixed_total_mixture", mode = "function", inherits = TRUE)) {
    stop("Source 00-mixture-generator.R before building multiscale scenarios.",
         call. = FALSE)
  }
  if (!is.numeric(n_features) || length(n_features) != 1L ||
      n_features != as.integer(n_features) || n_features < 9L) {
    stop("n_features must be one whole number of at least 9.", call. = FALSE)
  }
  n_features <- as.integer(n_features)
  if (!is.numeric(radii) || length(radii) < 3L || any(!is.finite(radii)) ||
      any(radii <= 0) || is.unsorted(radii, strictly = TRUE)) {
    stop("radii must be at least three finite, positive, strictly increasing values.",
         call. = FALSE)
  }
  if (!is.numeric(alpha) || length(alpha) != length(radii) ||
      any(!is.finite(alpha)) || any(alpha <= 0) || abs(sum(alpha) - 1) > 1e-12) {
    stop("alpha must contain one positive weight per radius and sum to one.",
         call. = FALSE)
  }

  coordinates <- cbind(x = seq_len(n_features) - 1, y = 0)
  domain <- abstract_domain(
    n_features, coordinates = coordinates,
    feature_ids = paste0("x", seq_len(n_features)), id = id
  )
  labels <- paste0("radius-", format(radii, trim = TRUE, scientific = FALSE))
  frames <- stats::setNames(lapply(radii, function(radius) {
    compile_frame(searchlights(radius, "conservative"), domain)
  }), labels)
  family_arguments <- c(frames, list(alpha = stats::setNames(alpha, labels)))
  family <- do.call(frame_family, family_arguments)

  alternating <- rep(c(1, -1), length.out = n_features)
  alternating <- alternating - mean(alternating)
  alternating <- alternating / sqrt(sum(alternating^2))
  specifications <- data.frame(
    scenario = c("broad_coherent", "mixed_broad_fine", "fine_configuration"),
    theta = c(0, pi / 4, pi / 2),
    stringsAsFactors = FALSE
  )
  scenarios <- stats::setNames(lapply(seq_len(nrow(specifications)), function(i) {
    fixed_total_mixture(
      n_features = n_features,
      total_magnitude = total_magnitude,
      theta = specifications$theta[[i]],
      seed = seed,
      configuration_direction = alternating
    )
  }), specifications$scenario)

  expected_rows <- list()
  for (scenario in names(scenarios)) {
    pattern <- scenarios[[scenario]]$effect_pattern
    for (scale in seq_along(frames)) {
      sparse_weights <- frames[[scale]]$weights
      weights <- as.matrix(sparse_weights)
      row_mass <- rowSums(weights)
      total <- sum(weights %*% pattern^2)
      coherent <- sum(drop(weights %*% pattern)^2 / row_mass)
      configuration <- total - coherent
      support_size <- rowSums(weights != 0)
      expected_rows[[length(expected_rows) + 1L]] <- data.frame(
        scenario = scenario,
        family = names(frames)[[scale]],
        scale = radii[[scale]],
        alpha = alpha[[scale]],
        total = alpha[[scale]] * total,
        coherent = alpha[[scale]] * coherent,
        configuration = alpha[[scale]] * configuration,
        coherent_share = coherent / total,
        configuration_share = configuration / total,
        n_nodes = nrow(weights),
        min_support = min(support_size),
        max_support = max(support_size),
        boundary_nodes = sum(support_size < max(support_size)),
        overlapping_features = sum(colSums(weights != 0) > 1L),
        sparse_weights = inherits(sparse_weights, "sparseMatrix"),
        stringsAsFactors = FALSE
      )
    }
  }
  expected <- do.call(rbind, expected_rows)
  rownames(expected) <- NULL

  list(
    schema_version = "matched-multiscale-v1",
    domain = domain,
    coordinates = coordinates,
    frames = frames,
    family = family,
    scenarios = scenarios,
    specifications = specifications,
    expected = expected,
    metadata = list(
      n_features = n_features,
      radii = radii,
      alpha = alpha,
      total_magnitude = total_magnitude,
      seed = as.integer(seed),
      configuration_basis = "centered alternating spatial frequency",
      fine_min_configuration_share_at_first_nonpoint = 0.85,
      mixed_transition_threshold = 0.45,
      mixed_expected_first_crossing = radii[[3L]],
      peak_scale_tolerance_steps = 0L,
      frame_normalization = "conservative",
      metric = "identity"
    )
  )
}

matched_multiscale_plan <- function(bundle, scenario,
                                    family = bundle$family,
                                    id_suffix = scenario) {
  if (!scenario %in% names(bundle$scenarios)) {
    stop("Unknown scenario `", scenario, "`.", call. = FALSE)
  }
  fixture <- bundle$scenarios[[scenario]]
  blocks <- list(run1 = fixture$effect_matrix, run2 = fixture$effect_matrix)
  relation_value <- relation(blocks, domain = bundle$domain)
  list(
    relation = relation_value,
    plan = plan_geometry(
      relation_value, family,
      cross_partitions(relation_value, independence = "independent")
    ),
    contrast = fixture$contrast
  )
}
