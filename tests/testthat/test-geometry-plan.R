geometry_plan_fixture <- function(reads = NULL, effects = 4L, features = 12L) {
  domain <- abstract_domain(
    features,
    coordinates = cbind(x = seq_len(features), y = 0),
    id = "geometry-plan-domain"
  )
  effect_names <- paste0("effect", seq_len(effects))
  values <- list(
    run1 = matrix(seq_len(effects * features) / 10, effects, features,
      dimnames = list(effect_names, NULL)),
    run2 = matrix(rev(seq_len(effects * features)) / 11, effects, features,
      dimnames = list(effect_names, NULL))
  )
  if (is.null(reads)) {
    relation <- relation(values, domain = domain)
  } else {
    source <- function(value) function(features) {
      reads$count <- reads$count + 1L
      value[, features, drop = FALSE]
    }
    revision <- paste0("sha256:", paste(rep("a", 64), collapse = ""))
    relation <- relation(
      lapply(values, source),
      source_dims = rep(list(c(effects, features)), 2L),
      effects = effect_names,
      domain = domain,
      capabilities = source_capabilities(
        block_read = TRUE, stable_revision = revision
      )
    )
  }
  frame <- compile_frame(searchlights(radius = 1.01), domain)
  list(
    relation = relation,
    frame = frame,
    pairing = cross_partitions(relation),
    domain = domain
  )
}

test_that("geometry plans are lazy, canonical, and readable", {
  reads <- new.env(parent = emptyenv())
  reads$count <- 0L
  fixture <- geometry_plan_fixture(reads)
  plan <- plan_geometry(
    fixture$relation, fixture$frame, fixture$pairing,
    compute = compute_policy(block_features = 3)
  )

  expect_s3_class(plan, "effect_geometry_plan")
  expect_identical(reads$count, 0L)
  expect_identical(plan$lowering, "additive_contraction")
  expect_identical(plan$metric_schedule$feature_additive, TRUE)
  expect_identical(plan$logical_shape, c(4L, 4L))
  expect_identical(plan$packed_width, 10L)
  expect_gt(plan$dense_payload_bytes, 0)
  printed <- capture.output(print(plan))
  expect_true(any(grepl("query-first", printed, fixed = TRUE)))
  expect_false(any(grepl("<environment:", printed, fixed = TRUE)))
})

test_that("compatibility and plan calls execute one identical lowering", {
  fixture <- geometry_plan_fixture()
  plan <- plan_geometry(fixture$relation, fixture$frame, fixture$pairing,
    compute = compute_policy(block_features = 4))
  planned <- geometry(plan)
  compatible <- geometry(
    fixture$relation, fixture$frame, fixture$pairing,
    compute = compute_policy(block_features = 4)
  )

  for (component in c("total", "coherent", "configuration")) {
    expect_equal(
      geometry_component(planned, component),
      geometry_component(compatible, component),
      tolerance = 0
    )
  }
  expect_identical(
    planned$metadata$execution_plan$lowering,
    "additive_form_contraction"
  )
  expect_false(planned$metadata$execution_plan$query_fused)
  expect_identical(planned$receipt$kernel_version,
    "additive-effect-form-v2")
  expect_identical(planned$receipt$scientific_plan_id,
    planned$metadata$scientific_plan_id)
})

test_that("query-first execution fuses the query before spatial contraction", {
  fixture <- geometry_plan_fixture(effects = 8L, features = 16L)
  plan <- plan_geometry(fixture$relation, fixture$frame, fixture$pairing,
    compute = compute_policy(block_features = 4))
  query <- bilinear_query(tcrossprod(c(1, -1, rep(0, 6))))
  old_route_called <- FALSE
  testthat::local_mocked_bindings(
    .streamed_crossgram_contraction = function(...) {
      old_route_called <<- TRUE
      stop("legacy route must not execute")
    },
    .package = "effectagram"
  )

  direct <- evaluate_geometry(plan, query = query)
  complete <- geometry(plan)
  late <- query_geometry(complete, query)

  expect_false(old_route_called)
  expect_equal(direct$values, late$values, tolerance = 1e-12)
  expect_true(direct$metadata$execution_plan$query_fused)
  expect_identical(
    direct$metadata$execution_plan$lowering,
    "additive_query_fused_contraction"
  )
  expect_identical(direct$receipt$kernel_version,
    "additive-query-fused-v2")
  expect_identical(direct$metadata$requirements$materialization,
    "direct_total")
  expect_lt(
    direct$receipt$memory$categories[["atom_block"]],
    complete$receipt$memory$categories[["atom_block"]]
  )
})

test_that("RDM and RSA views execute only their compiled query coordinates", {
  fixture <- geometry_plan_fixture(effects = 5L, features = 18L)
  # The full feature block makes the per-block accounting comparison
  # asymptotically honest: the structured query's durable payload is a small
  # constant, while its per-block saving scales with the feature block.
  plan <- plan_geometry(fixture$relation, fixture$frame, fixture$pairing,
    compute = compute_policy(block_features = 18))
  complete <- geometry(plan)
  model <- as.matrix(stats::dist(seq_len(5)))
  effect_names <- fixture$relation$effect_space$coordinates
  dimnames(model) <- list(effect_names, effect_names)
  old_route_called <- FALSE
  testthat::local_mocked_bindings(
    .streamed_crossgram_contraction = function(...) {
      old_route_called <<- TRUE
      stop("legacy route must not execute")
    },
    .package = "effectagram"
  )

  direct_rdm <- rdm(plan, component = "configuration")
  late_rdm <- rdm(complete, component = "configuration")
  direct_rsa <- rsa(plan, models = list(distance = model),
    component = "configuration")
  late_rsa <- rsa(complete, models = list(distance = model),
    component = "configuration")

  expect_false(old_route_called)
  expect_equal(direct_rdm$values, late_rdm$values, tolerance = 1e-12)
  expect_equal(direct_rsa$coefficients, late_rsa$coefficients,
    tolerance = 1e-12)
  expect_identical(
    direct_rdm$receipt$kernel_version,
    "additive-query-fused-v2"
  )
  expect_identical(
    direct_rsa$receipt$kernel_version,
    "additive-query-fused-v2"
  )
  expect_lt(
    direct_rsa$receipt$memory$categories[["atom_block"]],
    complete$receipt$memory$categories[["atom_block"]]
  )
})

test_that("contrast plans fuse energy components and signed marginals in one pass", {
  reads <- new.env(parent = emptyenv())
  reads$count <- 0L
  fixture <- geometry_plan_fixture(reads, effects = 4L, features = 15L)
  plan <- plan_geometry(fixture$relation, fixture$frame, fixture$pairing,
    compute = compute_policy(block_features = 5))
  weights <- c(effect1 = 1, effect2 = -1, effect3 = 0, effect4 = 0)

  direct <- contrast(plan, weights)
  direct_reads <- reads$count
  complete <- geometry(plan)
  late <- contrast(complete, weights)

  expect_equal(direct$signed, late$signed, tolerance = 1e-12)
  expect_equal(direct$total, late$total, tolerance = 1e-12)
  expect_equal(direct$coherent, late$coherent, tolerance = 1e-12)
  expect_equal(direct$configuration, late$configuration, tolerance = 1e-12)
  expect_identical(direct$receipt$kernel_version,
    "additive-query-fused-v2")
  expect_equal(
    direct$receipt$observed$task_counts[["completed"]],
    ceiling(fixture$relation$n_features / 5),
    tolerance = 0
  )
  expect_equal(direct_reads, 2L * ceiling(15 / 5), tolerance = 0)
})

test_that("executed receipts are derived from the compiled execution plan", {
  fixture <- geometry_plan_fixture()
  plan <- plan_geometry(fixture$relation, fixture$frame, fixture$pairing,
    compute = compute_policy(block_features = 5))
  query <- matrix(1, plan$packed_width, 2,
    dimnames = list(NULL, c("first", "second")))
  result <- evaluate_geometry(plan, query = query,
    component = "configuration")
  execution <- result$metadata$execution_plan

  expect_identical(colnames(result$values), c("first", "second"))
  expect_identical(execution$signature,
    result$metadata$execution_plan$signature)
  expect_identical(execution$kernel_version, result$receipt$kernel_version)
  expect_identical(result$metadata$tiles,
    result$receipt$observed$tiles)
  expect_identical(result$receipt$reduction_plan_id,
    execution$stage_signature)
})

test_that("mutated plans fail before source reads", {
  reads <- new.env(parent = emptyenv())
  reads$count <- 0L
  fixture <- geometry_plan_fixture(reads)
  plan <- plan_geometry(fixture$relation, fixture$frame, fixture$pairing)
  forged <- plan
  forged$frame$weights[1, 1] <- forged$frame$weights[1, 1] + 1

  expect_error(
    evaluate_geometry(forged, query = bilinear_query(diag(4))),
    "Locally normalized|identity|inconsistent"
  )
  expect_identical(reads$count, 0L)
})

test_that("query-only plans cannot silently accept materialization controls", {
  fixture <- geometry_plan_fixture()
  plan <- plan_geometry(fixture$relation, fixture$frame, fixture$pairing)

  expect_error(
    evaluate_geometry(plan, query = bilinear_query(diag(4)),
      compute = compute_policy(block_features = 1)),
    "cannot be combined with raw plan inputs"
  )
  expect_error(
    geometry(plan, at = fixture$frame),
    "cannot be combined with raw plan inputs"
  )
})

metric_plan_oracle <- function(fixture, K, contrast = NULL) {
  weights <- as.matrix(fixture$frame$weights)
  effects <- fixture$relation$effects
  relation_values <- lapply(fixture$relation$partitions, function(partition) {
    relation_block(
      fixture$relation, partition, seq_len(fixture$relation$n_features)
    )
  })
  names(relation_values) <- fixture$relation$partitions
  pack <- function(value) {
    out <- numeric(length(effects) * (length(effects) + 1L) / 2L)
    coordinate <- 0L
    for (column in seq_along(effects)) {
      for (row in column:length(effects)) {
        coordinate <- coordinate + 1L
        out[[coordinate]] <- value[row, column] *
          if (row == column) 1 else sqrt(2)
      }
    }
    out
  }
  forms <- lapply(seq_len(nrow(weights)), function(node) {
    support <- which(weights[node, ] > 0)
    root <- sqrt(weights[node, support])
    local_K <- K[support, support, drop = FALSE] * tcrossprod(root)
    value <- matrix(0, length(effects), length(effects))
    for (edge in seq_len(nrow(fixture$pairing))) {
      left <- relation_values[[fixture$pairing$left[[edge]]]][
        , support, drop = FALSE
      ]
      right <- relation_values[[fixture$pairing$right[[edge]]]][
        , support, drop = FALSE
      ]
      cross <- left %*% local_K %*% t(right)
      value <- value + fixture$pairing$weight[[edge]] *
        (cross + t(cross)) / 2
    }
    value
  })
  if (is.null(contrast)) {
    do.call(rbind, lapply(forms, pack))
  } else {
    matrix(vapply(forms, function(value) {
      drop(contrast %*% value %*% contrast)
    }, numeric(1)), ncol = 1L)
  }
}

metric_component_oracle <- function(fixture, K, component, contrast = NULL) {
  component <- match.arg(component, c("total", "coherent", "configuration"))
  weights <- as.matrix(fixture$frame$weights)
  effects <- fixture$relation$effects
  relation_values <- lapply(fixture$relation$partitions, function(partition) {
    relation_block(
      fixture$relation, partition, seq_len(fixture$relation$n_features)
    )
  })
  names(relation_values) <- fixture$relation$partitions
  pack <- function(value) {
    out <- numeric(length(effects) * (length(effects) + 1L) / 2L)
    coordinate <- 0L
    for (column in seq_along(effects)) {
      for (row in column:length(effects)) {
        coordinate <- coordinate + 1L
        out[[coordinate]] <- value[row, column] *
          if (row == column) 1 else sqrt(2)
      }
    }
    out
  }
  forms <- lapply(seq_len(nrow(weights)), function(node) {
    support <- which(weights[node, ] > 0)
    weight <- weights[node, support]
    root <- sqrt(weight)
    local_K <- K[support, support, drop = FALSE] * tcrossprod(root)
    a <- weight / sum(weight)
    denominator <- drop(crossprod(a, solve(local_K, a)))
    total <- matrix(0, length(effects), length(effects))
    coherent <- matrix(0, length(effects), length(effects))
    for (edge in seq_len(nrow(fixture$pairing))) {
      left <- relation_values[[fixture$pairing$left[[edge]]]][
        , support, drop = FALSE
      ]
      right <- relation_values[[fixture$pairing$right[[edge]]]][
        , support, drop = FALSE
      ]
      cross_total <- left %*% local_K %*% t(right)
      left_mean <- drop(left %*% a)
      right_mean <- drop(right %*% a)
      cross_coherent <- tcrossprod(left_mean, right_mean) / denominator
      edge_weight <- fixture$pairing$weight[[edge]]
      total <- total + edge_weight * (cross_total + t(cross_total)) / 2
      coherent <- coherent +
        edge_weight * (cross_coherent + t(cross_coherent)) / 2
    }
    switch(component,
      total = total,
      coherent = coherent,
      configuration = total - coherent
    )
  })
  if (is.null(contrast)) {
    do.call(rbind, lapply(forms, pack))
  } else {
    matrix(vapply(forms, function(value) {
      drop(contrast %*% value %*% contrast)
    }, numeric(1)), ncol = 1L)
  }
}

test_that("fixed dense metrics stream supports without pair materialization", {
  fixture <- geometry_plan_fixture(effects = 3L, features = 10L)
  raw <- matrix(seq_len(100) / 137, 10, 10)
  K <- crossprod(raw) + diag(0.75, 10)
  metric <- neural_metric(K, fixture$domain)
  plan <- plan_geometry(
    fixture$relation, fixture$frame, fixture$pairing,
    metric = metric,
    compute = compute_policy(workspace_bytes = 64 * 1024^2)
  )
  contrast_weights <- c(1, -1, 0)
  direct <- evaluate_geometry(
    plan, query = bilinear_query(tcrossprod(contrast_weights)),
    component = "total"
  )
  complete <- geometry(plan)
  late <- query_geometry(
    complete, bilinear_query(tcrossprod(contrast_weights)), "total"
  )

  expect_identical(plan$lowering, "support_streamed_pair_contraction")
  expect_identical(
    direct$metadata$execution_plan$lowering,
    "support_streamed_metric_query_contraction"
  )
  expect_identical(direct$receipt$kernel_version,
    "support-streamed-metric-v1")
  expect_false(direct$metadata$diagnostics$total$pair_atoms_materialized)
  expect_false(direct$metadata$diagnostics$total$pair_frame_materialized)
  expect_equal(
    direct$values,
    metric_plan_oracle(fixture, K, contrast_weights),
    tolerance = 2e-11,
    ignore_attr = TRUE
  )
  expect_equal(direct$values, late$values, tolerance = 2e-11)
  expect_equal(
    geometry_component(complete, "total"),
    metric_plan_oracle(fixture, K),
    tolerance = 2e-11
  )
  for (component in c("coherent", "configuration")) {
    expected <- metric_component_oracle(fixture, K, component)
    expect_equal(
      geometry_component(complete, component), expected,
      tolerance = 2e-11
    )
    queried <- evaluate_geometry(
      plan,
      query = bilinear_query(tcrossprod(contrast_weights)),
      component = component
    )
    expect_equal(
      queried$values,
      metric_component_oracle(fixture, K, component, contrast_weights),
      tolerance = 2e-11,
      ignore_attr = TRUE
    )
  }
  expect_equal(
    geometry_component(complete, "coherent") +
      geometry_component(complete, "configuration"),
    geometry_component(complete, "total"),
    tolerance = 2e-12
  )
})

test_that("fixed native-diagonal metrics retain the additive query collapse", {
  fixture <- geometry_plan_fixture(effects = 3L, features = 10L)
  diagonal <- seq(0.5, 1.4, length.out = 10)
  metric <- neural_metric(diag(diagonal), fixture$domain)
  plan <- plan_geometry(
    fixture$relation, fixture$frame, fixture$pairing, metric = metric
  )
  contrast_weights <- c(1, -1, 0)
  direct <- evaluate_geometry(
    plan, query = bilinear_query(tcrossprod(contrast_weights)),
    component = "total"
  )

  expect_identical(plan$lowering, "additive_contraction")
  expect_identical(
    direct$metadata$execution_plan$lowering,
    "additive_metric_query_fused_contraction"
  )
  expect_identical(direct$receipt$kernel_version,
    "additive-fixed-metric-v1")
  expect_equal(
    direct$values,
    metric_plan_oracle(fixture, diag(diagonal), contrast_weights),
    tolerance = 2e-11,
    ignore_attr = TRUE
  )
})
