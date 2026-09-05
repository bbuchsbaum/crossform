# Contract court for `model-coordinate-v1` (ticket M7).
#
# `design/model-coordinate-geometry-contract.md` quotes numbers measured on
# two fixtures: the base-R oracle (`design/oracles/model-coordinate-geometry.R`,
# whose laws the per-ticket test files already re-derive in base R) and the
# package route on the README fixture (`design/oracles/model-coordinate-
# package-route.R`). This file pins the second: every package-route number
# the contract quotes is reproduced here through public entry points against
# base matrix algebra written in the test, and the contract's own text is
# checked for the vocabulary it makes normative. Nothing here shells out to
# either script.

contract_path <- testthat::test_path("..", "..", "design",
  "model-coordinate-geometry-contract.md")
terminology_path <- testthat::test_path("..", "..", "design", "terminology.md")

read_contract <- function() {
  testthat::skip_if_not(file.exists(contract_path),
    "source-tree design contracts are intentionally excluded from the tarball")
  paste(readLines(contract_path, warn = FALSE), collapse = "\n")
}

readme_route <- function() {
  example <- example_fmri_effects()
  rel <- example$fit$relation
  basis <- model_basis(list(category = example$model_rdm),
    distance = "squared_euclidean", rank = 2, conditions = rel$effect_space)
  over <- cross_partitions(rel, independence = "independent",
    generalizes_over = "run")
  plan <- plan_geometry(rel, at = example$frame, over = over)
  list(example = example, rel = rel, basis = basis, over = over, plan = plan,
    Q = unname(basis$Q), n = length(rel$effects))
}

unsvec <- function(v, q) crossform:::.unsvec_symmetric(v, q)

test_that("the README category basis clamps to rank one and stays centered (section 1.1)", {
  route <- readme_route()
  expect_identical(route$basis$models$rank_requested, 2L)
  expect_identical(route$basis$models$rank_effective, 1L)
  expect_true(route$basis$models$rank_clamped)
  expect_identical(route$basis$dimension, 1L)
  expect_lt(route$basis$centering_error, 1e-12)
  # The trap the contract quotes at 1.15: a thin QR of the two-column literal
  # factor of the rank-1 Gram returns an uncentered second column.
  D <- route$example$model_rdm
  n <- nrow(D)
  H <- diag(n) - 1 / n
  e <- eigen(-0.5 * H %*% D %*% H, symmetric = TRUE)
  literal <- e$vectors[, 1:2] %*% diag(sqrt(pmax(e$values[1:2], 0)), 2L)
  Q_trap <- qr.Q(qr(literal))
  expect_equal(max(abs(colSums(Q_trap))), 1.15, tolerance = 0.01)
})

test_that("lowering commutes with the estimator on the README fixture (section 2.2 table)", {
  route <- readme_route()
  G <- materialize_geometry(route$plan)
  rel_m <- relation(
    stats::setNames(lapply(route$rel$partitions, function(p) {
      relation_block(route$rel, p, seq_len(route$rel$n_features))
    }), route$rel$partitions),
    extract = route$basis, domain = route$example$domain)
  G_m <- materialize_geometry(plan_geometry(rel_m, at = route$example$frame,
    over = cross_partitions(rel_m, independence = "independent",
      generalizes_over = "run")))
  Q <- route$Q
  q <- route$basis$dimension
  for (component in c("total", "coherent", "configuration")) {
    full <- geometry_component(G, component)
    lowered <- geometry_component(G_m, component)
    projected <- apply(full, 1, function(v) drop(t(Q) %*% unsvec(v, route$n) %*% Q))
    expect_lt(max(abs(projected - as.numeric(lowered))), 1e-12)
  }
  # And one operator of the proposal's query bank reproduces one entry
  # (section 2.3): the bank is admissible, only secondary.
  operator <- Q[, 1] %*% t(Q[, 1])
  bank <- evaluate_geometry(route$plan, query = bilinear_query(operator))
  direct <- apply(geometry_component(G, "total"), 1, function(v) {
    drop(t(Q[, 1]) %*% unsvec(v, route$n) %*% Q[, 1])
  })
  expect_lt(max(abs(as.numeric(bank$values) - direct)), 1e-12)
})

test_that("the trace split at the README peak reproduces the contract's numbers (section 5)", {
  route <- readme_route()
  reading <- model_geometry(route$plan, route$basis)
  # The witness picks the peak by the packed lowered form's absolute row sum;
  # with q = 1 that is |tr(S_x)| = |addressable|, the same measurement.
  peak <- which.max(abs(reading$addressable))
  expect_identical(peak, 144L)
  expect_equal(reading$centered_total[[peak]], 4.10545, tolerance = 1e-5)
  expect_equal(reading$addressable[[peak]], 4.10297, tolerance = 1e-5)
  expect_equal(reading$orthogonal[[peak]], 0.00248, tolerance = 1e-3)
  expect_lt(reading$identity_error, 1e-12)
  # tr(S_x) read from the lowered form directly equals the addressable term.
  rel_m <- relation(
    stats::setNames(lapply(route$rel$partitions, function(p) {
      relation_block(route$rel, p, seq_len(route$rel$n_features))
    }), route$rel$partitions),
    extract = route$basis, domain = route$example$domain)
  S <- geometry_component(materialize_geometry(plan_geometry(rel_m,
    at = route$example$frame, over = cross_partitions(rel_m,
      independence = "independent", generalizes_over = "run"))), "total")
  expect_equal(as.numeric(S[peak, ]), reading$addressable[[peak]],
    tolerance = 1e-10)
  # The lowered form is read by the existing readers unchanged (section 11).
  G_m <- materialize_geometry(plan_geometry(rel_m, at = route$example$frame,
    over = cross_partitions(rel_m, independence = "independent",
      generalizes_over = "run")))
  latent <- latent_geometry(G_m)
  expect_s3_class(geometry_spectrum(G_m), "effect_spectrum_view")
  expect_identical(latent$projection$clipped, 79L)
  expect_identical(latent$projection$measurements, 280L)
})

test_that("the contract's vocabulary is normative in the terminology table (section 9)", {
  contract <- read_contract()
  testthat::skip_if_not(file.exists(terminology_path))
  terminology <- paste(readLines(terminology_path, warn = FALSE),
    collapse = "\n")
  for (token in c("model basis", "lowered form", "addressable", "orthogonal",
    "latent_rank", "cross_fit", "span_fraction")) {
    expect_match(terminology, paste0("| ", token, " |"), fixed = TRUE)
  }
  for (phrase in c("model-addressable component", "model-orthogonal component",
    "cross-fitted model energy", "a neural subspace", "a denoised geometry",
    "an unbiased estimate")) {
    expect_match(terminology, phrase, fixed = TRUE)
    expect_match(contract, phrase, fixed = TRUE)
  }
  # The contract typesets the rank; the table spells it.
  expect_match(terminology, "latent rank-s model energy", fixed = TRUE)
  expect_match(contract, "latent rank-\\(s\\) model energy", fixed = TRUE)
  # The refusals the contract names exist under the namespace it names.
  for (capability in c("euclidean_model_geometry", "identified_model_coordinates",
    "model_coordinate_lowering", "model_geometry_test",
    "cross_fitted_model_energy")) {
    expect_match(contract, capability, fixed = TRUE)
  }
  expect_match(contract, "model-coordinate-v1", fixed = TRUE)
})

test_that("the two layers are kept apart on the package route (sections 3 and 7)", {
  route <- readme_route()
  policy <- metric_training_policy("exclude_evaluation")
  reading <- model_geometry(route$plan, route$basis, structure = "shared",
    rank = 1, training = policy)
  expect_identical(reading$layer$estimation,
    c("addressable", "orthogonal", "centered_total", "cross_fit"))
  expect_identical(reading$layer$latent, c("fit", "coefficients", "weights"))
  expect_identical(reading$reading, "latent descriptive layer; not for inference")
  # The plug-in energy is a PSD projection of the addressable term and never
  # below zero; the cross-fitted energy is signed. That it goes negative on
  # 28% of these searchlights is a fact about this fixture, not a law.
  expect_true(all(reading$fit$fitted_energy >= 0))
  expect_true(any(reading$addressable < 0))
  expect_true(any(reading$cross_fit$energy < 0))
  expect_false(isTRUE(all.equal(reading$cross_fit$energy,
    reading$fit$fitted_energy)))
})
