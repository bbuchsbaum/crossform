# Local reuse for one declared candidate family (layer 5) ------------------
# No process-global adaptive cache. Training and validation get separate
# objects; all form/path entries are discarded when their context changes.

.geometry_cache_pool_key <- function(basis_signature, weights) {
  .sha256_signature(list(model = basis_signature, weights = weights))
}

.geometry_cache_recipe_key <- function(basis_signature, weights, penalty) {
  .sha256_signature(list(model = basis_signature, weights = weights, penalty = penalty))
}

.new_geometry_training_cache <- function(candidates, role = c("training", "evaluation"), partner_bytes = 0) {
  role <- match.arg(role)
  if (!is.list(candidates) || !length(candidates)) .input_error("A training cache needs its declared candidate family.")
  .check_number(partner_bytes, "partner_bytes", nonnegative = TRUE)
  for (candidate in candidates) {
    .validate_model_basis(candidate$basis)
    .model_basis_weights(candidate$weights, candidate$basis$models$model)
    .check_number(candidate$penalty, "candidate penalty", nonnegative = TRUE)
  }
  cache <- new.env(parent = emptyenv())
  cache$candidates <- candidates
  cache$role <- role
  cache$partner_bytes <- partner_bytes
  cache$pools <- list(); cache$forms <- list(); cache$paths <- list()
  cache$context <- NULL; cache$form_id <- NULL; cache$bytes <- 0
  cache$counters <- stats::setNames(rep(0L, 7), c("context_resets", "form_hits", "form_misses",
    "pool_hits", "pool_misses", "path_hits", "path_misses"))
  cache$allowed_pools <- unique(vapply(candidates, function(x)
    .geometry_cache_pool_key(x$basis$signature, x$weights), ""))
  cache$allowed_recipes <- unique(vapply(candidates, function(x)
    .geometry_cache_recipe_key(x$basis$signature, x$weights, x$penalty), ""))
  class(cache) <- "effect_geometry_fit_cache"
  cache
}

.geometry_cache_check <- function(cache, role = NULL) {
  if (!is.environment(cache) || !inherits(cache, "effect_geometry_fit_cache") ||
      (!is.null(role) && !identical(cache$role, role))) {
    .input_error("Use a local cache with the matching training or evaluation role.")
  }
  invisible(cache)
}

.geometry_cache_bytes <- function(plan, candidates, role) {
  n <- as.double(plan$measurements)
  parts <- length(plan$task$left_relation$partitions)
  ids <- vapply(candidates, function(x) x$basis$signature, "")
  bytes <- 2 * as.double(utils::object.size(candidates))
  for (id in unique(ids)) {
    family <- candidates[ids == id]
    r <- as.double(family[[1]]$basis$dimension)
    paths <- length(unique(vapply(family, function(x)
      .geometry_cache_recipe_key(id, x$weights, x$penalty), "")))
    # Complete reduced forms, marginal/state overlap and their metadata.
    bytes <- bytes + n * (8 * (2 * r * (r + 1) + 4 * r * parts) + 1024)
    if (role == "training") {
      # Each path can retain its signed/shifted forms, eigensystem, support
      # and roots. This deliberately overcounts shared R references.
      bytes <- bytes + n * paths * (8 * (6 * r^2 + 24 * r) + 4096)
    }
  }
  if (!is.finite(bytes) || bytes > 2^53) .input_error("Training-cache dimensions exceed exact memory accounting.")
  1.25 * bytes
}

.geometry_cache_bind <- function(cache, prepared, row_block) {
  .geometry_cache_check(cache)
  key <- .sha256_signature(list(support = prepared$support$signature,
    target = prepared$target, compute = prepared$plan$compute,
    row_block = as.integer(row_block)), "geometry-cache-context-sha256:")
  if (!identical(cache$context, key)) {
    if (!is.null(cache$context)) cache$counters[["context_resets"]] <- cache$counters[["context_resets"]] + 1L
    cache$forms <- list(); cache$paths <- list(); cache$form_id <- NULL
    cache$context <- key
  }
  cache$bytes <- .geometry_cache_bytes(prepared$plan, cache$candidates, cache$role)
  cache$bytes + cache$partner_bytes
}

.geometry_cache_pool <- function(cache, basis, weights) {
  .geometry_cache_check(cache, "training")
  weights <- .model_basis_weights(weights, basis$models$model)
  key <- .geometry_cache_pool_key(basis$signature, weights)
  if (!key %in% cache$allowed_pools) .input_error("This model/weight recipe was not declared for the cache; create a new cache.")
  if (!is.null(cache$pools[[key]])) {
    cache$counters[["pool_hits"]] <- cache$counters[["pool_hits"]] + 1L
    return(cache$pools[[key]])
  }
  pool <- .model_basis_pool(basis, weights)
  cache$pools[[key]] <- pool
  cache$counters[["pool_misses"]] <- cache$counters[["pool_misses"]] + 1L
  pool
}

.geometry_cache_source <- function(cache, prepared) {
  .geometry_cache_check(cache)
  key <- .sha256_signature(list(context = cache$context, model = prepared$basis$signature,
    memory = prepared$workspace$geometry, row_block = prepared$workspace$row_block), "geometry-cached-form-sha256:")
  cache$form_id <- key
  if (!is.null(cache$forms[[key]])) {
    source <- cache$forms[[key]]
    .validate_effect_form(source, probe = FALSE)
    if (!identical(source$index, prepared$target$index) ||
        !identical(source$effect_space, prepared$basis$effect_space)) {
      .contract_error("Cached geometry does not match its ordered measurements and model coordinates.")
    }
    cache$counters[["form_hits"]] <- cache$counters[["form_hits"]] + 1L
    return(source)
  }
  source <- .geometry_fit_materialize(prepared, "memory", NULL)
  cache$forms[[key]] <- source
  cache$counters[["form_misses"]] <- cache$counters[["form_misses"]] + 1L
  source
}

.geometry_cache_path <- function(cache, row, S, pool, penalty) {
  .geometry_cache_check(cache, "training")
  .geometry_cache_require_recipe(cache, pool, penalty)
  key <- .sha256_signature(list(form = cache$form_id, row = row, signed_form = S,
    pool = pool$signature, penalty = penalty), "geometry-cached-spectrum-sha256:")
  if (!is.null(cache$paths[[key]])) {
    cache$counters[["path_hits"]] <- cache$counters[["path_hits"]] + 1L
    return(cache$paths[[key]])
  }
  path <- .geometry_fit_path(S, pool, penalty)
  cache$paths[[key]] <- path
  cache$counters[["path_misses"]] <- cache$counters[["path_misses"]] + 1L
  path
}

.geometry_cache_require_recipe <- function(cache, pool, penalty) {
  recipe <- .geometry_cache_recipe_key(pool$model_signature, pool$weights, penalty)
  if (!recipe %in% cache$allowed_recipes) .input_error("This penalty recipe was not declared for the cache; create a new cache.")
  invisible(TRUE)
}

.geometry_cache_summary <- function(cache) {
  .geometry_cache_check(cache)
  list(role = cache$role, context = cache$context, reserved_bytes = cache$bytes + cache$partner_bytes,
    counters = cache$counters, stored_forms = length(cache$forms),
    stored_pools = length(cache$pools), stored_paths = length(cache$paths))
}
