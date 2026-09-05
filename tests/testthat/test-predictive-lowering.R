# Dense raw-product oracles and actual-read guards for the lowered plan.
pl_fixture <- function(frame_kind = "regions", metric_kind = "identity") {
  set.seed(784)
  n <- 4L; p <- 5L
  effects <- effect_space(paste0("c", 1:n), basis_id = "condition-estimates",
    units = "percent", scale = 2)
  E <- cbind(diag(n), c(.2, -.1, .3, 0), c(0, .1, -.2, .4))
  ext <- effect_extractor(E, effects = effects)
  raw <- stats::setNames(lapply(1:6, function(i) matrix(rnorm(6*p), 6)), letters[1:6])
  calls <- new.env(parent = emptyenv()); calls$parts <- character(); calls$width <- integer()
  sources <- lapply(names(raw), function(part) {
    force(part)
    function(features) {
      if (!part %in% c("a", "b", "c")) stop("unused neural values were read")
      calls$parts <- c(calls$parts, part); calls$width <- c(calls$width, length(features))
      raw[[part]][, features, drop = FALSE]
    }
  }); names(sources) <- names(raw)
  caps <- lapply(raw, function(B) source_capabilities(TRUE,
    stable_revision = crossform:::.sha256_signature(B)))
  origins <- list(id = "lowering-origins", partitions = as.list(stats::setNames(names(raw), names(raw))),
    independence = "independent", assumption = "Independent runs and externally fixed preprocessing.")
  domain <- abstract_domain(p, coordinates = cbind(seq_len(p), 0, 0), id = "lowered-plan")
  rel <- relation(sources, extract = ext, source_dims = rep(list(c(6L, p)), 6),
    domain = domain, capabilities = caps, provenance = list(observation_origins = origins))
  frame_spec <- if (frame_kind == "regions") regions(c("A", "A", "B", "B", "B")) else searchlights(1.1)
  at <- compile_frame(frame_spec, domain)
  M <- if (metric_kind == "identity") diag(p) else diag(seq(1, 2, length.out = p)) + matrix(.2, p, p)
  metric <- if (metric_kind == "identity") NULL else neural_metric(M, domain)
  over <- pairing(c("a", "a", "b", "a"), c("b", "c", "c", "f"),
    weight = c(2, 3, 5, 0), directed = FALSE, independence = "independent", generalizes_over = "run")
  plan <- plan_geometry(rel, at, over, metric = metric, compute = compute_policy(block_features = 2))
  F <- rbind(c(1, 0), c(-1, .2), c(.5, 1), c(0, -1))
  rownames(F) <- effects$coordinates
  basis <- model_basis(features = F, conditions = effects, normalize = "trace")
  list(raw = raw, E = E, rel = rel, at = at, M = M, plan = plan,
    basis = basis, calls = calls, metric = metric, domain = domain)
}
pl_dense <- function(f, node, component) {
  w <- as.numeric(f$at$weights[node, ])
  keep <- which(w > 0)
  L <- diag(sqrt(w[keep]), length(keep)) %*% f$M[keep, keep, drop = FALSE] %*%
    diag(sqrt(w[keep]), length(keep))
  c <- w[keep]/sum(w[keep])
  C <- tcrossprod(c) / sum(c * solve(L, c))
  M <- switch(component, total = L, coherent = C, configuration = L - C)
  G <- matrix(0, 4, 4)
  edges <- data.frame(left = c("a", "a", "b"), right = c("b", "c", "c"), weight = c(.2, .3, .5))
  for (i in 1:3) {
    B1 <- f$E %*% f$raw[[edges$left[i]]][, keep, drop = FALSE]
    B2 <- f$E %*% f$raw[[edges$right[i]]][, keep, drop = FALSE]
    product <- B1 %*% M %*% t(B2)
    G <- G + edges$weight[i] * (product + t(product))/2
  }
  crossprod(f$basis$Q, G %*% f$basis$Q)
}

test_that("[T24 T25 T26] non-square lowering matches dense weighted partition products", {
  for (frame_kind in c("regions", "searchlights")) for (metric_kind in c("identity", "dense")) {
    f <- pl_fixture(frame_kind, metric_kind)
    prepared <- crossform:::.geometry_fit_prepare(f$plan, f$basis)
    expect_length(f$calls$parts, 0L)
    expect_identical(prepared$plan$task$left_relation$partitions, c("a", "b", "c"))
    expect_identical(prepared$plan$logical_shape, c(2L, 2L))
    expect_identical(prepared$plan$packed_width, 3L)
    expect_identical(prepared$plan$metric_schedule$signature, f$plan$metric_schedule$signature)
    expect_identical(prepared$plan$compute, f$plan$compute)
    expect_false(isTRUE(attr(prepared$plan$pairing, "directed")))
    expect_identical(prepared$support$observations, c("a", "b", "c"))
    G <- crossform:::.geometry_fit_materialize(prepared)
    expect_true(all(f$calls$parts %in% c("a", "b", "c")))
    if (metric_kind == "identity") expect_lte(max(f$calls$width), 2L)
    for (component in c("total", "coherent", "configuration")) {
      packed <- geometry_component(G, component)
      for (i in seq_len(nrow(packed))) {
        actual <- crossform:::.unsvec_symmetric(packed[i, ], 2)
        expect_equal(unname(actual), unname(pl_dense(f, i, component)), tolerance = 1e-10)
      }
    }
  }
})

test_that("[T25] prediction preparation never requests original-space geometry or trace queries", {
  f <- pl_fixture()
  calls <- 0L
  original <- materialize_geometry
  testthat::local_mocked_bindings(
    materialize_geometry = function(x, ...) {
      calls <<- calls + 1L
      if (!identical(x$logical_shape, c(2L, 2L))) stop("forbidden full geometry")
      original(x, ...)
    },
    evaluate_geometry = function(...) stop("forbidden original trace bank"),
    .package = "crossform")
  prepared <- crossform:::.geometry_fit_prepare(f$plan, f$basis)
  expect_identical(calls, 0L)
  G <- crossform:::.geometry_fit_materialize(prepared)
  expect_identical(calls, 1L)
  expect_identical(G$total$dim[[2]], 3L)
})

test_that("[T31] explicit effect meaning and units are checked before reads", {
  f <- pl_fixture()
  F <- f$basis$Q %*% f$basis$R
  wrong_space <- effect_space(f$rel$effects, basis_id = "different-estimand", units = "percent", scale = 2)
  wrong <- model_basis(features = F, conditions = wrong_space)
  refusal <- catch_refusal(crossform:::.geometry_fit_prepare(f$plan, wrong))
  expect_identical(refusal$reasons[[1]], "model_effect_space_mismatch")
  wrong <- model_basis(features = F, conditions = effect_space(f$rel$effects, units = "meters"))
  expect_identical(catch_refusal(crossform:::.geometry_fit_prepare(f$plan, wrong))$reasons[[1]],
    "model_effect_space_mismatch")
  wrong <- model_basis(features = F, conditions = rev(f$rel$effects))
  expect_identical(catch_refusal(crossform:::.geometry_fit_prepare(f$plan, wrong))$reasons[[1]],
    "model_effect_space_mismatch")
  expect_length(f$calls$parts, 0L)
})

test_that("[T31 T55] unsupported neural metric routes refuse before execution", {
  f <- pl_fixture(metric_kind = "dense")
  frozen <- neural_metric(f$M, f$domain, estimation = "learned_frozen",
    provenance = list(frozen = TRUE, training_signature = paste0("sha256:", strrep("a", 64))))
  learned <- plan_geometry(f$rel, f$at, f$plan$pairing, metric = frozen)
  expect_identical(catch_refusal(crossform:::.geometry_fit_prepare(learned, f$basis))$reasons[[1]],
    "unadmitted_learned_neural_metric")
  # Whitened plan construction itself eagerly reads values. Build that
  # unsupported input from ordinary matrices, then check our preflight.
  matrix_rel <- relation(f$raw, extract = f$rel$extractors, domain = f$domain,
    provenance = f$rel$provenance)
  whitened <- plan_geometry(matrix_rel, f$at, f$plan$pairing, metric = f$metric, composition = "whitened")
  expect_identical(catch_refusal(crossform:::.geometry_fit_prepare(whitened, f$basis))$reasons[[1]],
    "metric_schedule_not_lowerable")
  semidefinite <- neural_metric(diag(c(1, 1, 1, 1, 0)), f$domain)
  semi <- plan_geometry(f$rel, f$at, f$plan$pairing, metric = semidefinite)
  expect_identical(catch_refusal(crossform:::.geometry_fit_prepare(semi, f$basis))$reasons[[1]],
    "complete_form_metric_requires_spd")
  expect_length(f$calls$parts, 0L)
})

test_that("[T31 T48] a directed full-form request refuses instead of reaching a raw kernel error", {
  f <- pl_fixture()
  directed <- plan_geometry(f$rel, f$at,
    pairing("a", "b", directed = TRUE, independence = "independent", generalizes_over = "run"))
  expect_identical(catch_refusal(crossform:::.geometry_fit_prepare(directed, f$basis))$reasons[[1]],
    "directed_full_form_not_admitted")
  expect_length(f$calls$parts, 0L)
})
