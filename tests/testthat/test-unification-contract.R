# Normative vocabulary and implementation court for `unification-v1`.
#
# The oracle uses only base matrix algebra.  These tests separately pin the
# prose contract and compare the oracle with the public package route, so a
# source refactor cannot certify itself by calling the same private helper.

unification_contract_path <- testthat::test_path(
  "..", "..", "design", "unification-contract.md"
)
unification_oracle_path <- testthat::test_path(
  "..", "..", "design", "oracles", "unification-vocabulary.R"
)

read_unification_contract <- function() {
  testthat::skip_if_not(file.exists(unification_contract_path),
    "source-tree design contracts are intentionally excluded from the tarball")
  paste(readLines(unification_contract_path, warn = FALSE), collapse = "\n")
}

load_unification_oracle <- function() {
  testthat::skip_if_not(file.exists(unification_oracle_path),
    "source-tree design oracles are intentionally excluded from the tarball")
  court <- new.env(parent = baseenv())
  sys.source(unification_oracle_path, envir = court)
  court$unification_oracle
}

unification_svec <- function(value) {
  index <- which(lower.tri(value, diag = TRUE), arr.ind = TRUE)
  packed <- value[index]
  packed[index[, "row"] != index[, "col"]] <-
    sqrt(2) * packed[index[, "row"] != index[, "col"]]
  unname(packed)
}

test_that("the normative vocabulary and its scope boundary cannot disappear", {
  contract <- read_unification_contract()
  compact <- gsub("[[:space:]]+", " ", contract)
  terms <- c(
    "representation", "total_geometry", "coherent_component",
    "configurational_component", "frame_family", "scale_profile", "query"
  )
  for (term in terms) {
    marker <- paste0("<!-- contract-term: ", term, " -->")
    expect_identical(lengths(regmatches(contract,
      gregexpr(marker, contract, fixed = TRUE))), 1L, info = term)
  }

  mappings <- c(
    "relation()", "plan_geometry()", "geometry_component(..., \"total\")",
    "geometry_component(..., \"coherent\")",
    "geometry_component(..., \"configuration\")", "frame_family()",
    "coherence_spectrum()", "bilinear_query()", "pair_query()",
    "plan_population()"
  )
  for (mapping in mappings) {
    expect_match(contract, mapping, fixed = TRUE, info = mapping)
  }

  expect_match(compact,
    "not a claim that all representational statistics are linear queries",
    fixed = TRUE)
  expect_match(compact,
    "nonlinear rank correlations, locally data-adaptive queries, classification",
    fixed = TRUE)
  expect_match(compact,
    "a point-estimate theorem is not an uncertainty or coverage theorem",
    fixed = TRUE)
})

test_that("the independent vocabulary oracle passes all algebraic gates", {
  oracle <- load_unification_oracle()

  expect_named(oracle, c(
    "representation", "pairing", "frame_family", "frame_index", "geometry",
    "query", "first_query", "second_query", "scale_profile", "whole",
    "whole_query"
  ))
  expect_identical(dim(oracle$representation$run1), c(2L, 4L))
  expect_equal(colSums(oracle$frame_family), rep(1, 4L), tolerance = 1e-12)
  expect_equal(
    oracle$second_query$total,
    oracle$second_query$coherent + oracle$second_query$configuration,
    tolerance = 1e-12
  )
  expect_equal(
    oracle$scale_profile$total,
    oracle$scale_profile$alpha * oracle$whole_query,
    tolerance = 1e-12
  )
})

unification_public_fixture <- function(oracle) {
  domain <- abstract_domain(
    4L, coordinates = cbind(0:3, 0), feature_ids = paste0("v", 1:4),
    id = "unification-v1"
  )
  relation_value <- relation(oracle$representation, domain = domain)
  family <- compile_frame(
    searchlights(
      c(0.5, 1.01), "conservative", weights = c(0.5, 0.5)
    ),
    domain
  )
  plan <- plan_geometry(
    relation_value, family,
    cross_partitions(
      relation_value, independence = "independent", generalizes_over = "run"
    )
  )
  list(domain = domain, relation = relation_value, family = family, plan = plan)
}

test_that("public geometry implements the independent definitions", {
  oracle <- load_unification_oracle()
  fixture <- unification_public_fixture(oracle)
  form <- materialize_geometry(fixture$plan)

  expect_identical(as.character(fixture$family$index$family),
    oracle$frame_index$family)
  expect_equal(as.matrix(fixture$family$weights), oracle$frame_family,
    tolerance = 0, ignore_attr = TRUE)
  expect_equal(
    geometry_component(form, "total"),
    do.call(rbind, lapply(oracle$geometry,
      function(value) unification_svec(value$total))),
    tolerance = 1e-12, ignore_attr = TRUE
  )
  expect_equal(
    geometry_component(form, "coherent"),
    do.call(rbind, lapply(oracle$geometry,
      function(value) unification_svec(value$coherent))),
    tolerance = 1e-12, ignore_attr = TRUE
  )
  expect_equal(
    geometry_component(form, "configuration"),
    geometry_component(form, "total") - geometry_component(form, "coherent"),
    tolerance = 0
  )
})

test_that("public first, second, and scale queries match the oracle", {
  oracle <- load_unification_oracle()
  fixture <- unification_public_fixture(oracle)
  contrast <- oracle$query$contrast
  effect <- contrast_energy(fixture$plan, contrast)
  profile <- as.data.frame(coherence_spectrum(fixture$plan, contrast))

  expect_equal(effect$signed, oracle$first_query,
    tolerance = 1e-12, ignore_attr = TRUE)
  expect_equal(effect$total, oracle$second_query$total,
    tolerance = 1e-12, ignore_attr = TRUE)
  expect_equal(effect$coherent, oracle$second_query$coherent,
    tolerance = 1e-12, ignore_attr = TRUE)
  expect_equal(effect$configuration, oracle$second_query$configuration,
    tolerance = 1e-12, ignore_attr = TRUE)

  match_scale <- match(profile$scale, oracle$scale_profile$scale)
  expect_false(anyNA(match_scale))
  expected <- oracle$scale_profile[match_scale, ]
  expect_equal(profile$total, expected$total, tolerance = 1e-12)
  expect_equal(profile$coherent, expected$coherent, tolerance = 1e-12)
  expect_equal(profile$configuration, expected$configuration,
    tolerance = 1e-12)
  expect_equal(profile$coherence_fraction, expected$coherence_fraction,
    tolerance = 1e-12)

  query <- bilinear_query(oracle$query$operator)
  direct <- evaluate_geometry(fixture$plan, query, component = "total")
  expect_equal(drop(direct$values), oracle$second_query$total,
    tolerance = 1e-12, ignore_attr = TRUE)
})
