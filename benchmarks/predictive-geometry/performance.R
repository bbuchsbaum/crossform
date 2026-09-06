# Workload generation, independent dense checks and measured numerical work.
# Source this harness after loading crossform and performance-config.R.

pg_perf_raw <- function(Q, features, partition) {
  n <- nrow(Q); r <- ncol(Q)
  loadings <- outer(seq_len(r), features, function(a, b)
    sin(a * b * .071) + cos((a + 1) * b * .037))
  loadings <- sqrt(seq(r, 1) / r) * loadings
  noise <- .15 * outer(seq_len(n), features, function(a, b)
    sin(a * b * .113 + partition * 1.31) + cos(a * .39 + b * .07 * partition))
  Q %*% loadings + noise
}

pg_perf_fixture <- function(case, config = pg_performance_config()) {
  n <- case$conditions; r <- case$model_rank; p <- case$features
  # Analytic centered orthonormal DCT directions; no RNG state or dense source.
  Q <- sqrt(2/n) * cos(pi * outer(seq_len(n) - .5, seq_len(r)) / n)
  names <- paste0("condition", seq_len(n)); rownames(Q) <- names
  roots <- .9^(seq_len(r) - 1); roots <- roots / sum(roots)
  K <- Q %*% diag(roots) %*% t(Q)
  parts <- paste0("run", seq_len(case$partitions))
  reads <- new.env(parent = emptyenv())
  reads$calls <- stats::setNames(integer(length(parts)), parts)
  reads$columns <- reads$calls; reads$maximum_width <- 0L
  sources <- lapply(seq_along(parts), function(a) {
    force(a)
    function(features) {
      reads$calls[[a]] <- reads$calls[[a]] + 1L
      reads$columns[[a]] <- reads$columns[[a]] + length(features)
      reads$maximum_width <- max(reads$maximum_width, length(features))
      pg_perf_raw(Q, features, a)
    }
  }); names(sources) <- parts
  manifest <- list(id = paste0("performance-", case$id),
    partitions = stats::setNames(as.list(paste0("raw-", parts)), parts),
    independence = "independent", assumption = "Declared independent runs; deterministic computational witness.")
  domain <- abstract_domain(p, id = paste0("performance-", case$id))
  capabilities <- lapply(seq_along(parts), function(a) source_capabilities(TRUE,
    stable_revision = crossform:::.sha256_signature(list(generator = "pg-perf-v1", case = case, partition = a))))
  rel <- relation(sources, source_dims = rep(list(c(n, p)), length(parts)),
    effects = names, domain = domain, capabilities = capabilities,
    provenance = list(observation_origins = manifest))
  members <- lapply(seq_len(case$nodes), function(i)
    sort(as.integer(((i - 1L) * 17L + seq_len(case$support) - 1L) %% p + 1L)))
  at <- additive_frame(members = members, normalization = "local", domain = domain,
    measurements = paste0("node", seq_len(case$nodes)))
  half <- length(parts) / 2
  make_pairing <- function(ids) {
    edges <- utils::combn(ids, 2)
    pairing(edges[1, ], edges[2, ], weight = seq_len(ncol(edges)),
      independence = "independent", generalizes_over = "run")
  }
  train_pair <- make_pairing(parts[seq_len(half)])
  test_pair <- make_pairing(parts[half + seq_len(half)])
  policy <- compute_policy(block_features = config$block_features, workspace_bytes = config$workspace_bytes)
  train <- plan_geometry(rel, at, train_pair, compute = policy)
  test <- plan_geometry(rel, at, test_pair, compute = policy)
  basis <- model_basis(kernels = list(model = K), conditions = rel$effect_space, normalize = "trace")
  dense_train <- crossform:::.geometry_prediction_pair_plan(train, as.data.frame(train$pairing))
  dense_test <- crossform:::.geometry_prediction_pair_plan(test, as.data.frame(test$pairing))
  list(case = case, Q = Q, roots = roots, basis = basis, train = train, test = test,
    dense_train = dense_train, dense_test = dense_test,
    reads = reads, members = members, parts = parts, config = config)
}

pg_perf_admission <- function(f, route) {
  c <- f$config; rank <- f$case$response_rank; r <- f$basis$dimension
  train <- crossform:::.geometry_fit_prepare(f$train, f$basis)
  test <- crossform:::.geometry_fit_prepare(f$test, f$basis)
  if (route == "dense") {
    # Full original-space output must itself be admitted before execution.
    full <- lapply(list(f$dense_train, f$dense_test), crossform:::.compile_geometry_execution_plan, storage = "memory")
    # The second execution coexists with the retained first full form and
    # output summaries. Include that overlap in the declared dense baseline.
    retained <- 16 * f$case$nodes * f$case$conditions * (f$case$conditions + 1) / 2 +
      16 * f$case$nodes * (r^2 + 8 * r)
    planned <- max(vapply(full, function(x) x$memory$planned_workspace_bytes, 0)) + retained
    if (planned > c$workspace_bytes) stop("Dense comparison exceeds its prespecified workspace budget")
    return(list(route = route, planned_workspace_bytes = planned,
      retained_overlap_bytes = retained, geometry = lapply(full, `[[`, "memory")))
  }
  prototype <- list(values = list(representation = route,
    dim = c(f$case$nodes, r * rank + rank + 2L * r)))
  fitting <- crossform:::.geometry_prediction_workspace(train, rank, route, c$row_block, "fit")$workspace
  scoring <- crossform:::.geometry_prediction_workspace(test, rank, route, c$row_block,
    "score", fit = prototype, modes = TRUE)$workspace
  list(route = route, planned_workspace_bytes = max(fitting$planned_workspace_bytes,
    scoring$planned_workspace_bytes), fit = fitting, score = scoring)
}

pg_perf_instrument <- function(expected_dimension = NULL) {
  counters <- new.env(parent = emptyenv())
  counters$kernels <- 0L; counters$eigen <- 0L; counters$paths <- 0L
  counters$materializations <- 0L; counters$packed_widths <- integer()
  counters$expected_dimension <- expected_dimension
  # Isolated worker-owned instrumentation, never a production cache.
  assign(".pg_performance_counter", counters, envir = globalenv())
  namespace <- asNamespace("crossform")
  kernels <- c(".packed_effect_form_atoms_cpp", ".coherent_effect_form_atoms_cpp", ".support_metric_gram_pairs_cpp")
  increment <- function(field) substitute({
    counter <- get(".pg_performance_counter", envir = globalenv())
    counter[[FIELD]] <- counter[[FIELD]] + 1L
  }, list(FIELD = field))
  for (name in kernels) trace(name, increment("kernels"), where = namespace, print = FALSE)
  trace("eigen", increment("eigen"), where = baseenv(), print = FALSE)
  trace(".geometry_fit_path", increment("paths"), where = namespace, print = FALSE)
  trace(".geometry_fit_materialize", quote({
    counter <- get(".pg_performance_counter", envir = globalenv())
    counter$materializations <- counter$materializations + 1L
    counter$packed_widths <- c(counter$packed_widths, prepared$plan$packed_width)
    if (!is.null(counter$expected_dimension) &&
        prepared$plan$logical_shape[[1]] != counter$expected_dimension)
      stop("Forbidden original-space materialization in predictive route")
  }), where = namespace, print = FALSE)
  list(counters = counters, cleanup = function() {
    for (name in c(kernels, ".geometry_fit_path", ".geometry_fit_materialize")) untrace(name, where = namespace)
    untrace("eigen", where = baseenv())
    rm(".pg_performance_counter", envir = globalenv())
  })
}

pg_perf_prediction_rows <- function(fit, rows, row_block) {
  r <- fit$model$dimension
  answer <- matrix(0, length(rows), r * r)
  blocks <- split(seq_along(rows), ceiling(seq_along(rows) / row_block))
  for (indices in blocks) {
    predictions <- crossform:::.geometry_prediction_rows(fit, rows[indices])
    for (j in seq_along(indices)) answer[indices[j], ] <- as.numeric(tcrossprod(predictions[[j]]$factor))
  }
  answer
}

pg_perf_direct_form <- function(f, node, plan) {
  ids <- f$members[[node]]; edges <- as.data.frame(plan$pairing)
  blocks <- lapply(seq_along(f$parts), function(a) pg_perf_raw(f$Q, ids, a))
  names(blocks) <- f$parts
  G <- matrix(0, f$case$conditions, f$case$conditions)
  for (edge in seq_len(nrow(edges))) {
    product <- tcrossprod(blocks[[edges$left[edge]]], blocks[[edges$right[edge]]]) / length(ids)
    G <- G + edges$weight[edge] * (product + t(product)) / 2
  }
  G
}

# Independent base-R eigensolver, without production fit or packed codecs.
pg_perf_dense_fit <- function(S, J, rank, penalty) {
  model <- if (is.list(J)) J else eigen(J, symmetric = TRUE)
  keep <- which(model$values > 1e-10 * max(abs(model$values)))
  U <- model$vectors[, keep, drop = FALSE]; d <- model$values[keep]
  shifted <- crossprod(U, S %*% U) - diag(penalty/d, length(d))
  e <- eigen((shifted + t(shifted))/2, symmetric = TRUE)
  active <- head(which(e$values > 0), rank)
  Z <- U %*% e$vectors[, active, drop = FALSE]
  Z <- sweep(Z, 2, sqrt(e$values[active]), "*")
  list(form = tcrossprod(Z), amplitudes = e$values[active], factor = Z,
    signed = eigen(S, symmetric = TRUE, only.values = TRUE)$values, shifted = e$values)
}

pg_perf_reference <- function(f, rows) {
  J <- tcrossprod(f$basis$R)
  forms <- matrix(0, length(rows), f$basis$dimension^2); gain <- numeric(length(rows))
  trace <- numeric(length(rows))
  for (j in seq_along(rows)) {
    G <- pg_perf_direct_form(f, rows[j], f$train)
    S <- crossprod(f$basis$Q, G %*% f$basis$Q)
    fit <- pg_perf_dense_fit(S, J, f$case$response_rank, f$config$penalty)
    test <- pg_perf_direct_form(f, rows[j], f$test)
    F <- f$basis$Q %*% fit$form %*% t(f$basis$Q)
    forms[j, ] <- as.numeric(fit$form)
    gain[j] <- 2 * sum(F * test) - sum(F^2)
    trace[j] <- sum(diag(S))
  }
  list(forms = forms, gain = gain, signed_trace = trace)
}

pg_perf_error <- function(actual, expected) max(abs(actual - expected)) / max(1, max(abs(expected)))

pg_perf_rows <- function(f) {
  if (f$case$id != "large") seq_len(f$case$nodes) else
    unique(as.integer(round(seq(1, f$case$nodes, length.out = f$config$large_oracle_rows))))
}

pg_perf_public <- function(f, storage, directory) {
  c <- f$config; rows <- pg_perf_rows(f); n <- f$case$nodes; r <- f$basis$dimension
  fit <- fit_geometry(f$train, f$basis, rank = f$case$response_rank, penalty = c$penalty,
    storage = storage, storage_path = if (storage == "block") file.path(directory, "fit") else NULL,
    row_block = c$row_block)
  score <- score_geometry(fit, f$test, modes = TRUE, storage = storage,
    storage_path = if (storage == "block") file.path(directory, "score") else NULL, row_block = c$row_block)
  spectra <- matrix(0, n, 2 * r)
  mode_values <- matrix(0, n, 3 * fit$layout$budget)
  for (start in seq.int(1L, n, by = c$row_block)) {
    at <- seq.int(start, min(n, start + c$row_block - 1L))
    values <- crossform:::.geometry_prediction_rows(fit, at)
    for (j in seq_along(at)) spectra[at[j], ] <- c(values[[j]]$signed_spectrum, values[[j]]$shifted_spectrum)
    if (fit$layout$budget) mode_values[at, ] <- crossform:::.read_geometry_store(score$mode_values, at)[,
      seq_len(3 * fit$layout$budget), drop = FALSE]
  }
  files <- if (storage == "block") list.files(directory, recursive = TRUE, full.names = TRUE) else character()
  stored_bytes <- if (storage == "block") sum(file.info(files)$size) else
    as.double(object.size(fit) + object.size(score))
  list(table = score$table, modes = mode_values, spectra = spectra,
    diagnostics = fit$diagnostics, rows = rows,
    forms = pg_perf_prediction_rows(fit, rows, c$row_block),
    stored_bytes = stored_bytes, storage_kind = if (storage == "block") "durable_file_bytes" else "R_object_bytes",
    planned_workspace_bytes = max(fit$workspace$planned_workspace_bytes, score$workspace$planned_workspace_bytes),
    workspace = list(fit = fit$workspace, score = score$workspace),
    receipts = list(fit = fit$receipt, score = score$receipt))
}

pg_perf_full_space <- function(f) {
  c <- f$config; n <- f$case$nodes; r <- f$basis$dimension; s <- f$case$response_rank
  train <- materialize_geometry(f$dense_train)
  test <- materialize_geometry(f$dense_test)
  rows <- pg_perf_rows(f); forms <- matrix(0, n, r * r)
  spectra <- matrix(0, n, 2 * r); modes <- matrix(0, n, 3 * s)
  table <- data.frame(measurement = paste0("node", seq_len(n)), gain = numeric(n),
    rank = integer(n), inner_product = numeric(n), prediction_norm_sq = numeric(n))
  diagnostics <- data.frame(measurement = table$measurement, rank = integer(n),
    prediction_norm_sq = numeric(n), penalty_trace = numeric(n), compressed_residual_sq = numeric(n),
    objective = numeric(n), source_negative_mass = numeric(n), shifted_negative_mass = numeric(n),
    shifted_truncated_positive_mass = numeric(n))
  model <- eigen(tcrossprod(f$basis$R), symmetric = TRUE)
  for (start in seq.int(1L, n, by = c$row_block)) {
    at <- seq.int(start, min(n, start + c$row_block - 1L))
    train_values <- geometry_component(train, "total", rows = at)
    test_values <- geometry_component(test, "total", rows = at)
    for (j in seq_along(at)) {
      i <- at[j]
      G <- crossform:::.unsvec_symmetric(train_values[j, ], f$case$conditions)
      E <- crossform:::.unsvec_symmetric(test_values[j, ], f$case$conditions)
      S <- crossprod(f$basis$Q, G %*% f$basis$Q)
      prediction <- pg_perf_dense_fit(S, model, s, c$penalty)
      amp <- prediction$amplitudes; k <- length(amp); active <- seq_len(k)
      F <- f$basis$Q %*% prediction$form %*% t(f$basis$Q)
      weighted <- colSums(prediction$factor * (crossprod(f$basis$Q, E %*% f$basis$Q) %*% prediction$factor))
      evidence <- weighted / amp; gains <- 2 * weighted - amp^2
      forms[i, ] <- as.numeric(prediction$form)
      spectra[i, ] <- c(prediction$signed, prediction$shifted)
      modes[i, active] <- amp; modes[i, s + active] <- evidence; modes[i, 2 * s + active] <- gains
      table[i, -1L] <- list(2 * sum(F * E) - sum(F^2), k, sum(F * E), sum(F^2))
      penalty_trace <- sum(rowSums(crossprod(model$vectors, prediction$factor)^2) / model$values)
      residual <- sum((S - prediction$form)^2)
      truncated <- setdiff(which(prediction$shifted > 0), active)
      diagnostics[i, -1L] <- list(k, sum(F^2), penalty_trace, residual,
        .5 * residual + c$penalty * penalty_trace, -sum(pmin(prediction$signed, 0)),
        -sum(pmin(prediction$shifted, 0)), sum(prediction$shifted[truncated]))
    }
  }
  list(table = table, modes = modes, spectra = spectra, diagnostics = diagnostics,
    rows = rows, forms = forms[rows, , drop = FALSE],
    stored_bytes = as.double(object.size(list(table, modes, spectra, diagnostics, forms))),
    storage_kind = "R_numeric_reference_bytes", receipts = list(train = train$receipt, test = test$receipt))
}

pg_perf_allocation <- function(path) {
  con <- file(path, "r"); on.exit(close(con))
  count <- 0; total <- 0; largest <- 0
  repeat {
    lines <- readLines(con, n = 10000L, warn = FALSE)
    if (!length(lines)) break
    values <- suppressWarnings(as.numeric(sub(" .*", "", lines)))
    values <- values[is.finite(values)]
    count <- count + length(values); total <- total + sum(values)
    largest <- max(c(largest, values))
  }
  list(allocations = count, cumulative_R_bytes = total, largest_R_allocation = largest,
    scope = "Rprofmem allocations; native workspace additionally observed through RSS and compiler receipts")
}

pg_perf_legacy_overlap <- function(f) {
  old <- model_geometry(f$dense_train, f$basis, structure = "shared", rank = f$case$response_rank)
  unpenalized <- f; unpenalized$config$penalty <- 0
  reference <- pg_perf_reference(unpenalized, seq_len(f$case$nodes))
  forms <- t(vapply(seq_len(f$case$nodes), function(i)
    as.numeric(f$basis$R %*% old$coefficients[, , i] %*% t(f$basis$R)),
    numeric(f$basis$dimension^2)))
  errors <- c(signed_trace = pg_perf_error(old$addressable, reference$signed_trace),
    unpenalized_form = pg_perf_error(forms, reference$forms))
  stopifnot(all(errors <= f$config$relative_tolerance))
  list(errors = errors, nodes = f$case$nodes,
    scope = "Descriptive reader parity only for signed model trace and unpenalized shared form; its projector energy is not predictive gain.")
}

pg_perf_reuse_comparison <- function(f) {
  grid <- expand.grid(rank = c(0L, f$case$response_rank, f$basis$dimension),
    penalty = c(f$config$penalty, 4 * f$config$penalty))
  candidates <- crossform:::.geometry_candidates(lapply(seq_len(nrow(grid)), function(i)
    list(basis = f$basis, rank = grid$rank[i], penalty = grid$penalty[i])))
  run <- function(reuse) {
    caches <- crossform:::.geometry_selection_caches(f$train, f$test, candidates, reuse)
    f$reads$calls[] <- 0L; f$reads$columns[] <- 0L; f$reads$maximum_width <- 0L
    instrument <- pg_perf_instrument(f$basis$dimension)
    on.exit(instrument$cleanup(), add = TRUE)
    gains <- matrix(0, f$case$nodes, length(candidates))
    # Bounded numerical witnesses; the gain includes every measurement.
    rows <- unique(as.integer(round(seq(1, f$case$nodes, length.out = min(32L, f$case$nodes)))))
    forms <- array(0, c(length(rows), f$basis$dimension^2, length(candidates)))
    workspace <- numeric(length(candidates))
    started <- proc.time()[["elapsed"]]
    for (j in seq_along(candidates)) {
      fit <- crossform:::.geometry_fit_candidate(f$train, candidates[[j]], "total",
        f$config$row_block, cache = caches$training)
      score <- crossform:::.geometry_score_candidate(fit, f$test,
        f$config$row_block, cache = caches$evaluation)
      gains[, j] <- score$table$gain
      forms[, , j] <- pg_perf_prediction_rows(fit, rows, f$config$row_block)
      workspace[j] <- max(fit$workspace$planned_workspace_bytes, score$workspace$planned_workspace_bytes)
    }
    list(seconds = proc.time()[["elapsed"]] - started, gains = gains, forms = forms,
      calls = as.list(instrument$counters), reads = as.list(f$reads),
      maximum_planned_bytes = max(workspace),
      cache = crossform:::.geometry_selection_cache_summary(caches))
  }
  uncached <- run(FALSE); cached <- run(TRUE)
  errors <- c(gain = pg_perf_error(cached$gains, uncached$gains),
    form = pg_perf_error(cached$forms, uncached$forms))
  stopifnot(all(errors <= f$config$relative_tolerance),
    cached$calls$materializations == 2L, uncached$calls$materializations == 2L * length(candidates),
    cached$calls$paths == 2L * f$case$nodes,
    uncached$calls$paths == length(candidates) * f$case$nodes)
  uncached$gains <- uncached$forms <- NULL
  cached$gains <- cached$forms <- NULL
  list(grid = grid, uncached = uncached, cached = cached, errors = errors,
    scope = "Six fixed candidate recipes in one worker, uncached then cached; warm instrumented timings without an idle-machine or universal speedup claim.")
}
