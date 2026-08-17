# Contract tests for R/coupling-plan.R: `coupling()`, the adjoint-side
# closure that reuses a geometry plan's relation, spatial frame, and pairing
# instead of asking for a second declaration of each.

coupling_plan_fixture <- local({
  built <- NULL
  function() {
    if (!is.null(built)) {
      return(built)
    }
    set.seed(20260819L)
    q <- 6L
    p <- 4L
    domain <- abstract_domain(p, id = "coupling-plan:neural:v1")
    effects <- effect_space(paste0("trial", seq_len(q)),
      basis_id = "coupling-plan:trials:v1")
    run1 <- matrix(rnorm(q * p), q, p, dimnames = list(effects$coordinates, NULL))
    run2 <- matrix(rnorm(q * p), q, p, dimnames = list(effects$coordinates, NULL))
    values <- list(run1 = run1, run2 = run2)
    relation <- relation(values, effects = effects, domain = domain)
    frame <- compile_frame(regions(c("r1", "r1", "r2", "r2")), domain)
    over <- pairing(
      relation$partitions, relation$partitions, directed = TRUE,
      self_pairs = "allow_biased", independence = "not_independent"
    )
    built <<- list(
      q = q, p = p, domain = domain, effects = effects, values = values,
      relation = relation, frame = frame, over = over,
      plan = plan_geometry(relation, frame, over),
      by = variation_query(
        (diag(q) - matrix(1 / q, q, q)) / (q - 1L), effects, "trial",
        "joint_covariance", provenance = list(estimator = "centered")
      )
    )
    built
  }
})

test_that("coupling reuses the plan's frame, relation, and pairing", {
  fixture <- coupling_plan_fixture()
  form <- coupling(fixture$plan, cbind("r1", "r2"), by = fixture$by)

  expect_s3_class(form, "effect_measurement_form")
  expect_identical(nrow(form$block_index), 1L)
  expect_identical(form$block_index$left, "r1")
  expect_identical(form$block_index$right, "r2")
  expect_identical(form$capabilities$sampling_axis, "trial")
  expect_true(form$capabilities$joint_covariance)

  # The block is the measured neural evidence form: L_{r1} B' H B L_{r2}',
  # aggregated over the plan's own partition pairs. The oracle rebuilds it
  # from the stored relation blocks and the frame weights.
  legs <- as.matrix(measurement_frame(fixture$frame)$legs$r1$operator)
  right_legs <- as.matrix(measurement_frame(fixture$frame)$legs$r2$operator)
  edges <- as.data.frame(fixture$over)
  oracle <- Reduce(`+`, lapply(seq_len(nrow(edges)), function(edge) {
    left <- fixture$values[[edges$left[[edge]]]]
    right <- fixture$values[[edges$right[[edge]]]]
    edges$weight[[edge]] *
      legs %*% crossprod(left, fixture$by$operator %*% right) %*% t(right_legs)
  }))

  expect_equal(effect_coupling(form)$values[[1L]], oracle, tolerance = 1e-10,
    ignore_attr = TRUE)
})

test_that("coupling closes the plan the other way from evaluate_geometry", {
  fixture <- coupling_plan_fixture()
  # A geometry plan closes the neural boundary and leaves the experimental
  # axes open; `coupling()` closes the experimental boundary instead. Read
  # the same cross-partition plan both ways and the scalar evidence agrees.
  plan <- plan_geometry(
    fixture$relation, fixture$frame,
    cross_partitions(fixture$relation, independence = "independent")
  )
  operator <- as.matrix(fixture$by$operator)
  by <- pair_query(operator, fixture$effects, fixture$effects)
  block <- effect_coupling(
    coupling(plan, cbind("r1", "r1"), by = by)
  )$values[[1L]]
  geometry <- evaluate_geometry(plan, query = bilinear_query(operator))
  node <- match("r1", measurement_frame(fixture$frame)$node_ids)

  # tr(L Q_H L') = <L'L, Q_H> = <H, B (L'L) B'>, which is the geometry the
  # forward closure reports at the same node.
  expect_equal(sum(diag(block)), drop(geometry$values)[[node]],
    tolerance = 1e-10)
  expect_gt(abs(sum(diag(block))), 1e-6)
})

test_that("coupling refuses a rectangular cross-axis plan", {
  fixture <- coupling_plan_fixture()
  right_effects <- effect_space(c("x", "y"),
    basis_id = "coupling-plan:right:v1")
  right <- relation(
    list(s1 = matrix(rnorm(2 * fixture$p), 2L, fixture$p,
      dimnames = list(c("x", "y"), NULL))),
    effects = right_effects, domain = fixture$domain
  )
  rectangular <- plan_geometry(
    fixture$relation, fixture$frame,
    pairing("run1", "s1", 1, directed = TRUE, independence = "independent"),
    right = right
  )
  refusal <- catch_refusal(
    coupling(rectangular, cbind("r1", "r2"), by = fixture$by)
  )

  expect_s3_class(refusal, "effect_capability_refusal")
  expect_identical(refusal$capability, "self_form_coupling")
  expect_identical(refusal$namespace, "coupling_plans")
  expect_identical(refusal$reasons, "rectangular_cross_axis_plan")
  expect_match(refusal$remedies, "omitting `right`")
})

test_that("both closure declarations are required rather than defaulted", {
  fixture <- coupling_plan_fixture()

  expect_error(coupling(fixture$plan, by = fixture$by),
    "`between` is required", class = "effect_input_error")
  expect_error(coupling(fixture$plan, cbind("r1", "r2")), "`by` is required",
    class = "effect_input_error")
  expect_error(coupling(fixture$relation, cbind("r1", "r2"), by = fixture$by),
    "geometry plan|effect_geometry_plan", class = "effect_input_error")
})

test_that("between resolves node names and indices against the plan frame", {
  fixture <- coupling_plan_fixture()
  named <- coupling(fixture$plan, cbind("r1", "r2"), by = fixture$by)
  indexed <- coupling(fixture$plan, cbind(1L, 2L), by = fixture$by)
  framed <- coupling(
    fixture$plan, data.frame(from = "r1", to = "r2"), by = fixture$by
  )

  expect_equal(
    effect_coupling(indexed)$values[[1L]],
    effect_coupling(named)$values[[1L]], tolerance = 0
  )
  expect_equal(
    effect_coupling(framed)$values[[1L]],
    effect_coupling(named)$values[[1L]], tolerance = 0
  )

  expect_error(
    coupling(fixture$plan, cbind("r1", "r9"), by = fixture$by),
    "`between` to names `r9`, which is not a measurement of the plan's frame"
  , class = "effect_input_error")
  expect_error(
    coupling(fixture$plan, cbind("r9", "r1"), by = fixture$by),
    "`between` from names `r9`"
  , class = "effect_input_error")
  expect_error(
    coupling(fixture$plan, cbind(1L, 7L), by = fixture$by),
    "`between` to indices must be whole numbers in 1:2"
  , class = "effect_input_error")
  expect_error(
    coupling(fixture$plan, cbind(1.5, 2), by = fixture$by),
    "whole numbers in 1:2"
  , class = "effect_input_error")
  expect_error(
    coupling(fixture$plan, c("r1", "r2"), by = fixture$by),
    "nonempty two-column matrix of frame node pairs"
  , class = "effect_input_error")
  expect_error(
    coupling(fixture$plan, matrix(character(), 0L, 2L), by = fixture$by),
    "nonempty two-column matrix of frame node pairs"
  , class = "effect_input_error")
  expect_error(
    coupling(fixture$plan, cbind("r1", "r2", "r1"), by = fixture$by),
    "nonempty two-column matrix of frame node pairs"
  , class = "effect_input_error")
})

test_that("mode, reducer, and over are recorded in the compiled identity", {
  fixture <- coupling_plan_fixture()
  between <- cbind("r1", "r2")
  total <- coupling(fixture$plan, between, by = fixture$by)
  coherent <- coupling(fixture$plan, between, by = fixture$by,
    mode = "coherent")
  aggregated <- coupling(fixture$plan, between, by = fixture$by,
    reducer = reduce_partitions())
  narrowed <- coupling(fixture$plan, between, by = fixture$by,
    over = pairing("run1", "run1", 1, directed = TRUE,
      self_pairs = "allow_biased", independence = "not_independent"))

  # A coherent frame keeps only the weighted-mean direction of each node, so
  # its block is 1-by-1 where the total-mode block is 2-by-2.
  expect_identical(dim(effect_coupling(total)$values[[1L]]), c(2L, 2L))
  expect_identical(dim(effect_coupling(coherent)$values[[1L]]), c(1L, 1L))
  expect_error(coupling(fixture$plan, between, by = fixture$by, mode = "raw"),
    "'arg' should be one of")

  # Each declaration changes the estimand and therefore the plan identity.
  signatures <- c(total$plan$signature, coherent$plan$signature,
    aggregated$plan$signature, narrowed$plan$signature)
  expect_identical(anyDuplicated(signatures), 0L)
  expect_identical(total$plan$stages$reduction$operation$order,
    "aggregate_first")
  expect_identical(aggregated$plan$stages$reduction$operation$order,
    "edge_first")

  # Omitting `over` reuses the plan's own pairing exactly.
  expect_identical(
    total$plan$signature,
    coupling(fixture$plan, between, by = fixture$by,
      over = fixture$over)$plan$signature
  )
  expect_false(identical(total$plan$signature, narrowed$plan$signature))
})
