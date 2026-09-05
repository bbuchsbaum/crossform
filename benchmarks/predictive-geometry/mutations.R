# Mutations of the actual production closures, installed only in a disposable
# worker namespace. AST replacement must find its exact declared target.
pg_mutation_replace <- function(fun, from, to, expected = 1L) {
  hits <- 0L
  walk <- function(node) {
    if (identical(node, from)) { hits <<- hits + 1L; return(to) }
    if (is.call(node)) return(as.call(lapply(as.list(node), walk)))
    node
  }
  changed <- walk(body(fun))
  if (hits != expected) stop("Mutation definition drift: expected ", expected, " AST matches; found ", hits)
  body(fun) <- changed
  fun
}

pg_mutation_full_form <- function(clip) {
  prepare <- crossform:::.geometry_fit_prepare
  body(prepare) <- substitute({
    answer <- ORIGINAL
    answer$mutation_original_plan <- plan
    answer
  }, list(ORIGINAL = body(prepare)))
  materialize <- crossform:::.geometry_fit_materialize
  body(materialize) <- substitute({
    answer <- ORIGINAL
    full <- materialize_geometry(prepared$mutation_original_plan)
    packed <- geometry_component(full, prepared$component)
    small <- matrix(0, nrow(packed), prepared$plan$packed_width)
    for (i in seq_len(nrow(packed))) {
      G <- .unsvec_symmetric(packed[i, ], prepared$basis$n_conditions)
      if (CLIP) G <- .latent_rank_psd_form(G)$form
      small[i, ] <- .svec_symmetric(crossprod(prepared$basis$Q, G %*% prepared$basis$Q))
    }
    answer$total <- .memory_geometry_store(small, "symmetric_packed")
    answer
  }, list(ORIGINAL = body(materialize), CLIP = clip))
  list(.geometry_fit_prepare = prepare, .geometry_fit_materialize = materialize)
}

pg_mutation_catalog <- function() {
  entry <- function(id, symbol, transform, file, primary, description) {
    list(id = id, symbol = symbol, transform = transform, file = file,
      primary = primary, description = description)
  }
  one <- function(symbol, from, to, expected = 1L) {
    force(symbol); force(from); force(to); force(expected)
    function() stats::setNames(list(pg_mutation_replace(get(symbol, asNamespace("crossform")),
      from, to, expected)), symbol)
  }
  list(
    entry("inverse_as_kernel", ".geometry_fit_path",
      one(".geometry_fit_path", quote(penalty/d), quote(penalty * d)),
      "predictive-fit-numerics", "T12", "Use kernel eigenvalues instead of their inverse in the penalty"),
    entry("positive_penalty_shift", ".geometry_fit_path",
      one(".geometry_fit_path", quote(restricted/2 + t(restricted)/2 - diag(penalty/d, length(d))),
        quote(restricted/2 + t(restricted)/2 + diag(penalty/d, length(d)))),
      "predictive-fit-numerics", "T12", "Add rather than subtract the penalty"),
    entry("fit_pooled_null_space", ".geometry_fit_path", function() {
      fun <- pg_mutation_replace(crossform:::.geometry_fit_path,
        quote(restricted <- crossprod(U, S %*% U)), quote({
          inverse <- U %*% diag(1/d, length(d)) %*% t(U)
          restricted <- S; U <- diag(nrow(S)); d <- rep(1, nrow(S))
        }))
      fun <- pg_mutation_replace(fun, quote(diag(penalty/d, length(d))), quote(penalty * inverse))
      list(.geometry_fit_path = fun)
    }, "predictive-fit-numerics", "T09", "Fit the whole union after subtracting the pooled pseudoinverse"),
    entry("clip_compressed_first", ".geometry_fit_path",
      one(".geometry_fit_path", quote(S <- .geometry_fit_symmetric(S)),
        quote(S <- .latent_rank_psd_form(.geometry_fit_symmetric(S))$form)),
      "predictive-fit-numerics", "T15", "Clip the compressed signed form before shifting"),
    entry("clip_full_first", ".geometry_fit_prepare + .geometry_fit_materialize",
      function() pg_mutation_full_form(TRUE), "predictive-mutation-witnesses", "T14",
      "Materialize and clip original geometry before compression"),
    entry("normalize_before_truncation", ".model_basis_factor",
      one(".model_basis_factor", quote(retained_trace <- sum(e$values[keep])), quote(retained_trace <- positive_mass)),
      "predictive-kernel", "T05", "Normalize retained factors by the pre-truncation positive trace"),
    entry("drop_packing_weight", ".unsvec_symmetric",
      one(".unsvec_symmetric", quote(sqrt(2)), quote(1)),
      "predictive-score", "T34", "Decode packed off-diagonals without sqrt(2)"),
    entry("omit_score_factor_two", ".geometry_score_row",
      one(".geometry_score_row", quote(2 * weighted_evidence), quote(weighted_evidence)),
      "predictive-score", "T33", "Drop the leading factor two in predictive gain"),
    entry("trace_prediction_cost", ".geometry_score_row",
      one(".geometry_score_row", quote(amp^2), quote(amp), 2L),
      "predictive-score", "T33", "Charge trace instead of squared Frobenius prediction norm"),
    entry("omit_prediction_cost", ".geometry_score_row",
      one(".geometry_score_row", quote(2 * weighted_evidence - amp^2), quote(2 * weighted_evidence)),
      "predictive-score", "T33", "Omit the frozen prediction cost"),
    entry("double_prediction_cost", ".geometry_score_row",
      one(".geometry_score_row", quote(2 * weighted_evidence - amp^2), quote(2 * weighted_evidence - 2 * amp^2)),
      "predictive-score", "T33", "Subtract the prediction cost twice"),
    entry("prune_on_test_gain", ".geometry_score_row",
      one(".geometry_score_row", quote(gains <- 2 * weighted_evidence - amp^2),
        quote(gains <- pmax(2 * weighted_evidence - amp^2, 0))),
      "predictive-score", "T36", "Discard harmful modes using their final test gain"),
    entry("origins_as_labels", ".geometry_support_independence",
      one(".geometry_support_independence", quote(intersect(training$observations, evaluation$observations)),
        quote(intersect(names(training$endpoints), names(evaluation$endpoints)))),
      "predictive-training-support", "T27", "Treat local endpoint labels as complete observation identity"),
    entry("origins_as_content_hash", ".geometry_support_independence",
      one(".geometry_support_independence", quote(intersect(training$observations, evaluation$observations)),
        quote(intersect(vapply(training$endpoints, `[[`, "", "source_revision"),
          vapply(evaluation$endpoints, `[[`, "", "source_revision")))),
      "predictive-training-support", "T28", "Treat content hashes as complete observation identity"),
    entry("train_all_parent_partitions", ".geometry_training_support",
      one(".geometry_training_support",
        quote(products <- as.data.frame(over[over$weight > 0, c("left", "right", "weight"), drop = FALSE])),
        quote(products <- as.data.frame(cross_partitions(relation, independence = "independent",
          generalizes_over = attr(over, "generalizes_over")))[, c("left", "right", "weight"), drop = FALSE])),
      "predictive-training-support", "T26", "Replace selected training edges by all parent partition pairs"),
    entry("train_all_outer_folds", ".geometry_nested",
      one(".geometry_nested", quote(train <- budget_plan(fold$training)), quote(train <- train_plan)),
      "predictive-selection", "T50", "Use all outer-training data when fitting each inner candidate"),
    entry("positional_measurements", ".geometry_prediction_match",
      one(".geometry_prediction_match", quote(mapping <- match(training_target$measurement_ids, evaluation_target$measurement_ids)),
        quote(mapping <- seq_along(training_target$measurement_ids))),
      "predictive-score-preflight", "T40", "Match measurements by position instead of identity"),
    entry("stale_source_cache", ".geometry_cache_bind",
      one(".geometry_cache_bind", quote(prepared$support$signature),
        quote(.sha256_signature(list(products = prepared$support$products,
          observations = prepared$support$observations, upstream = prepared$support$upstream)))),
      "predictive-cache", "T53", "Cache ignores changed source revisions and extractors"),
    entry("stale_fold_cache", ".geometry_cache_bind",
      one(".geometry_cache_bind", quote(prepared$support$signature),
        quote(.sha256_signature(list(endpoints = prepared$support$endpoints,
          upstream = prepared$support$upstream, scope = prepared$support$scope_signature)))),
      "predictive-mutation-witnesses", "T53", "Cache ignores changed training edges and weights"),
    entry("materialize_full_geometry", ".geometry_fit_prepare + .geometry_fit_materialize",
      function() pg_mutation_full_form(FALSE), "predictive-lowering", "T25",
      "Compute original full geometry while returning only its reduced form")
  )
}

pg_mutation_instrument <- function(replacements) {
  counter <- new.env(parent = emptyenv()); counter$hits <- stats::setNames(integer(length(replacements)), names(replacements))
  for (name in names(replacements)) {
    fun <- replacements[[name]]
    scope <- new.env(parent = environment(fun)); scope$.mutation_counter <- counter
    body(fun) <- substitute({
      .mutation_counter$hits[[NAME]] <- .mutation_counter$hits[[NAME]] + 1L
      ORIGINAL
    }, list(NAME = name, ORIGINAL = body(fun)))
    environment(fun) <- scope
    replacements[[name]] <- fun
  }
  list(replacements = replacements, counter = counter)
}
