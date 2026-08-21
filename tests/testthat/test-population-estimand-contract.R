# Contract and current provenance court for `population-estimand-v1`.

pe_contract_path <- testthat::test_path(
  "..", "..", "design", "population-estimand-contract.md"
)
pe_oracle_path <- testthat::test_path(
  "..", "..", "design", "oracles", "population-estimand-targets.R"
)

pe_contract_text <- function() {
  testthat::skip_if_not(file.exists(pe_contract_path),
    "source-tree design contracts are intentionally excluded from the tarball")
  paste(readLines(pe_contract_path, warn = FALSE), collapse = "\n")
}

pe_load_oracle <- function() {
  testthat::skip_if_not(file.exists(pe_oracle_path),
    "source-tree design oracles are intentionally excluded from the tarball")
  court <- new.env(parent = baseenv())
  sys.source(pe_oracle_path, envir = court)
  court$population_estimand_oracle
}

test_that("both targets and every coverage definition are normative", {
  contract <- pe_contract_text()
  markers <- c(
    "<!-- population-target: all_planned -->",
    "<!-- population-target: available_at_node -->",
    "<!-- population-diagnostic: node_sample_size -->",
    "<!-- population-diagnostic: effective_sample_size -->",
    "<!-- population-diagnostic: sink_mass -->",
    "<!-- population-diagnostic: transport_quality -->"
  )
  for (marker in markers) {
    expect_identical(lengths(regmatches(contract,
      gregexpr(marker, contract, fixed = TRUE))), 1L, info = marker)
  }

  compact <- gsub("[[:space:]]+", " ", contract)
  required <- c(
    "conditional on the following realized objects",
    "do **not** propagate uncertainty from estimating registration",
    "must not infer availability from whether `y` happens to equal zero",
    "planned_subject_unavailable",
    "available_at_node",
    "Transport mass is not silently used as a precision weight",
    "a hash alone is insufficient",
    "computed subject-first from the underlying transported values and then refit",
    "Inverse-probability or transport-quality weighting requires a frozen"
  )
  for (phrase in required) expect_match(compact, phrase, fixed = TRUE)
})

test_that("the independent target oracle distinguishes selection policies", {
  oracle <- pe_load_oracle()

  expect_equal(oracle$coverage$n, c(5, 4, 4), tolerance = 0)
  expect_true(all(is.finite(oracle$all_planned$full$coefficients)))
  expect_true(all(is.na(
    oracle$all_planned$left_missing$coefficients
  )))
  expect_true(all(is.na(
    oracle$all_planned$right_missing$coefficients
  )))
  expect_identical(
    oracle$available_at_node$left_missing$subjects,
    paste0("s", 1:4)
  )
  expect_identical(
    oracle$available_at_node$right_missing$subjects,
    paste0("s", 2:5)
  )
  expect_true(all(oracle$coverage$mass_n_eff <= oracle$coverage$n))
  expect_lt(oracle$coverage_association[["left_missing"]], 0)
  expect_gt(oracle$coverage_association[["right_missing"]], 0)

  for (node in names(oracle$subject_set_id)) {
    key <- oracle$subject_set_id[[node]]
    expect_identical(
      oracle$subject_set_dictionary[[key]],
      oracle$available_at_node[[node]]$subjects,
      info = node
    )
  }
})

pe_contract_effects <- function() {
  effect_space(c("face", "house"), basis_id = "pe-contract-effects:v1")
}

pe_contract_subject <- function(id, features) {
  domain <- abstract_domain(
    features, coordinates = cbind(seq_len(features) - 1),
    feature_ids = paste0("f", seq_len(features)), id = id
  )
  values <- function(offset) {
    matrix(seq_len(2 * features) / features + offset, 2L, features,
      dimnames = list(c("face", "house"), NULL))
  }
  relation_value <- relation(
    list(run1 = values(0), run2 = values(0.25)),
    effects = pe_contract_effects(), domain = domain
  )
  plan_geometry(
    relation_value,
    compile_frame(voxelwise("conservative"), domain),
    cross_partitions(
      relation_value, independence = "independent", generalizes_over = "run"
    )
  )
}

pe_contract_fixture <- function() {
  sizes <- c(s01 = 4L, s02 = 5L, s03 = 6L)
  subjects <- stats::setNames(lapply(names(sizes), function(id) {
    pe_contract_subject(id, sizes[[id]])
  }), names(sizes))
  transport <- stats::setNames(lapply(names(sizes), function(id) {
    anatomical_transport(
      native_coords = cbind(seq_len(sizes[[id]]) - 1),
      group_coords = cbind(c(0, 3)), semantics = "budget"
    )
  }), names(sizes))
  list(subjects = subjects, transport = transport)
}

test_that("public plans bind planned subjects to realized transports", {
  fixture <- pe_contract_fixture()
  plan <- plan_population(fixture$subjects, fixture$transport)

  expect_identical(plan$subject_index$subject, names(fixture$subjects))
  expect_identical(
    plan$subject_index$plan_id,
    unname(vapply(fixture$subjects, `[[`, character(1), "scientific_plan_id"))
  )
  expect_identical(
    plan$subject_index$transport_signature,
    unname(vapply(fixture$transport, `[[`, character(1), "signature"))
  )

  changed <- fixture$transport
  changed$s03 <- anatomical_transport(
    native_coords = cbind(seq_len(6L) - 1),
    group_coords = cbind(c(0, 3)), semantics = "budget", radius = 1
  )
  moved <- plan_population(fixture$subjects, changed)
  expect_false(identical(plan$scientific_plan_id, moved$scientific_plan_id))
})

test_that("population receipts preserve the plan-level identity table", {
  fixture <- pe_contract_fixture()
  plan <- plan_population(fixture$subjects, fixture$transport)
  fit <- estimate_population(
    plan, rbind(`face-house` = c(face = 1, house = -1))
  )

  expect_identical(fit$receipt$population_plan_id, plan$scientific_plan_id)
  expect_identical(fit$receipt$subjects$subject, plan$subject_index$subject)
  expect_identical(
    fit$receipt$subjects$plan_id, plan$subject_index$plan_id
  )
  expect_identical(
    fit$receipt$subjects$transport_signature,
    plan$subject_index$transport_signature
  )
})
