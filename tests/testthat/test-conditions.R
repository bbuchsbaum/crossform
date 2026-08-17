# The condition taxonomy is a public contract: `?crossform_conditions` tells
# readers they may branch on `effect_input_error`, `effect_contract_error`,
# `effect_invariant_error`, and `effect_capability_refusal`, and that the
# first three share `effect_error`. These tests hold the constructors, the
# guards that raise them, and the entry points that reach both to that
# promise, so a later refactor cannot quietly reclassify a failure.

test_that("each constructor carries its class vector and its fields", {
  input <- tryCatch(
    .input_error("bad argument", arg = "at", received = "3.5",
      expected = "one whole number"),
    condition = function(condition) condition
  )
  expect_s3_class(input, "effect_input_error")
  expect_identical(class(input),
    c("effect_input_error", "effect_error", "error", "condition"))
  expect_identical(conditionMessage(input), "bad argument")
  expect_null(conditionCall(input))
  expect_identical(input$arg, "at")
  expect_identical(input$received, "3.5")
  expect_identical(input$expected, "one whole number")

  contract <- tryCatch(
    .contract_error("these do not belong together", arg = "receipt",
      received = "sha256:aaa", expected = "sha256:bbb"),
    condition = function(condition) condition
  )
  expect_identical(class(contract),
    c("effect_contract_error", "effect_error", "error", "condition"))
  expect_identical(contract$received, "sha256:aaa")
  expect_identical(contract$expected, "sha256:bbb")

  invariant <- tryCatch(.invariant_error("a norm came out negative."),
    condition = function(condition) condition)
  expect_identical(class(invariant),
    c("effect_invariant_error", "effect_error", "error", "condition"))
  # An invariant failure is never the caller's fault, so it must say so.
  expect_match(conditionMessage(invariant), "crossform bug", fixed = TRUE)
  expect_match(conditionMessage(invariant), "a norm came out negative.",
    fixed = TRUE)
})

test_that("all three are errors, and `effect_error` catches exactly them", {
  for (raise in list(
    function() .input_error("x"),
    function() .contract_error("x"),
    function() .invariant_error("x")
  )) {
    expect_error(raise())
    expect_identical(
      tryCatch(raise(), effect_error = function(condition) "effect"),
      "effect")
    expect_identical(
      tryCatch(raise(), error = function(condition) "error"), "error")
  }

  # A capability refusal is not an `effect_error`: it is not a mistake in the
  # input, and code that branches on the taxonomy must be able to tell them
  # apart.
  refusal <- tryCatch(
    .capability_refusal("no", capability = "c", namespace = "n"),
    condition = function(condition) condition)
  expect_s3_class(refusal, "effect_capability_refusal")
  expect_false(inherits(refusal, "effect_error"))
})

test_that("catch_refusal returns refusals and lets other errors through", {
  refusal <- catch_refusal(
    .capability_refusal("no", capability = "guaranteed_psd",
      namespace = "views", reasons = "r", remedies = "m"))
  expect_s3_class(refusal, "effect_capability_refusal")
  expect_identical(refusal$capability, "guaranteed_psd")
  expect_null(catch_refusal(1 + 1))
  expect_error(catch_refusal(.input_error("still an error")),
    class = "effect_input_error")
})

test_that("an effect_error prints as a record of what was asked", {
  condition <- tryCatch(
    .input_error("`at` must be one measurement index.", arg = "at",
      received = "3.5", expected = "one whole number"),
    condition = function(condition) condition)
  lines <- format(condition)
  expect_identical(lines[[1L]], "<effect_input_error>")
  expect_true(any(grepl("must be one measurement index", lines)))
  expect_true(any(grepl("^  arg:  +at$", lines)))
  expect_true(any(grepl("^  received:  +3\\.5$", lines)))
  expect_true(any(grepl("^  expected:  +one whole number$", lines)))
  expect_output(print(condition), "<effect_input_error>", fixed = TRUE)

  # Fields the caller did not supply are omitted rather than printed empty.
  bare <- tryCatch(.contract_error("two objects disagree."),
    condition = function(condition) condition)
  expect_false(any(grepl("arg:", format(bare), fixed = TRUE)))
})

test_that("the guard helpers state the expectation and the observed value", {
  expect_error(.check_string(c("a", "b"), "id"),
    class = "effect_input_error")
  expect_error(.check_string(c("a", "b"), "id"),
    "`id` must be one nonempty character string; received a character vector")
  expect_error(.check_string("", "id"), "received an empty string")
  expect_error(.check_string(NA_character_, "id"), "received NA")
  expect_identical(.check_string("", "id", allow_empty = TRUE), "")
  expect_identical(.check_string("x", "id"), "x")

  expect_error(.check_flag(NA, "sort"), "`sort` must be TRUE or FALSE")
  expect_error(.check_flag(c(TRUE, TRUE), "sort"), class = "effect_input_error")
  expect_true(.check_flag(TRUE, "sort"))

  expect_error(.check_count(2.5, "workers"),
    "`workers` must be one positive whole number; received 2.5")
  expect_error(.check_count(0, "workers"), class = "effect_input_error")
  expect_error(.check_count(Inf, "workers"), class = "effect_input_error")
  expect_identical(.check_count(3, "workers"), 3L)
  expect_identical(.check_count(0, "offset", min = 0L), 0L)

  expect_error(.check_number(-1, "tolerance", positive = TRUE),
    "`tolerance` must be one positive finite number; received -1")
  expect_error(.check_number(Inf, "tolerance"), class = "effect_input_error")
  expect_identical(.check_number(Inf, "tolerance", finite = FALSE), Inf)

  expect_error(.check_matrix(matrix(c(1, NA), 1L), "values"),
    "must be a numeric matrix with no missing or infinite entries")
  expect_error(.check_matrix(matrix(1, 1L), "values", nrow = 2L),
    "with 2 rows and any number of columns")
  expect_identical(.check_matrix(matrix(1, 1L), "values"), matrix(1, 1L))

  expect_error(.check_class(1L, "effect_frame", "at", from = "compile_frame()"),
    "`at` must be an object of class `effect_frame` from `compile_frame\\(\\)`")
  expect_error(.check_class(1L, "effect_frame", "at"),
    class = "effect_input_error")
})

test_that("the guard helpers carry arg, received, and expected as fields", {
  condition <- tryCatch(.check_count(-2, "workers"),
    condition = function(condition) condition)
  expect_identical(condition$arg, "workers")
  expect_identical(condition$received, "-2")
  expect_identical(condition$expected, "one positive whole number")
})

test_that("the shape predicates answer rather than raise", {
  expect_true(.is_string("a"))
  expect_false(.is_string(""))
  expect_true(.is_string("", allow_empty = TRUE))
  expect_false(.is_string(NA_character_))
  expect_false(.is_string(c("a", "b")))
  # Total on any input: a predicate is evaluated in any position of a chain.
  expect_false(.is_string(list(1)))
  expect_false(.is_string(NULL))

  expect_true(.is_strings(c("a", "b")))
  expect_false(.is_strings(c("a", "")))
  expect_false(.is_strings(c("a", NA)))
  expect_true(.is_strings(c("a", "a")))
  expect_false(.is_strings(c("a", "a"), unique = TRUE))

  expect_true(.is_flag(FALSE))
  expect_false(.is_flag(NA))
  expect_false(.is_flag(list(TRUE)))

  expect_true(.is_number(1.5))
  expect_false(.is_number(Inf))
  expect_true(.is_number(Inf, finite = FALSE))
  expect_false(.is_number("1"))

  expect_true(.is_count(3))
  expect_false(.is_count(3.5))
  expect_false(.is_count(0))
  expect_true(.is_count(0, min = 0L))
  expect_false(.is_count(Inf))

  expect_true(.is_finite_numeric(c(1, 2)))
  expect_false(.is_finite_numeric(c(1, NA)))
  expect_true(.is_finite_matrix(matrix(1, 1L)))
  expect_false(.is_finite_matrix(matrix("a", 1L)))
  expect_false(.is_finite_matrix(1))
})

test_that("a sealed record's prologue and signature checks disagree usefully", {
  record <- structure(list(kind = "x", signature = "sha256:aa"),
    class = "effect_demo_record")
  expect_true(.sealed_fields(record, "effect_demo_record",
    c("kind", "signature")))
  expect_false(.sealed_fields(record, "effect_other", c("kind", "signature")))
  expect_false(.sealed_fields(record, "effect_demo_record", c("signature")))
  expect_false(.sealed_fields(unclass(record), "effect_demo_record",
    c("kind", "signature")))

  expect_true(.check_signature("sha256:aa", "sha256:aa", "unused"))
  condition <- tryCatch(
    .check_signature(
      paste0("sha256:", strrep("a", 64L)),
      paste0("sha256:", strrep("b", 64L)),
      "Demo-record identity is inconsistent."),
    condition = function(condition) condition)
  expect_s3_class(condition, "effect_contract_error")
  expect_identical(conditionMessage(condition),
    "Demo-record identity is inconsistent.")
  # The digests are attached rather than pasted into the sentence, and are
  # truncated: the reader needs to see that two identities differ, not all 64
  # hexadecimal digits of each.
  expect_identical(condition$received, "sha256:aaaaaaaaaaaa...")
  expect_identical(condition$expected, "sha256:bbbbbbbbbbbb...")
})

test_that("public entry points raise the class their failure deserves", {
  expect_error(abstract_domain(-1L), class = "effect_input_error")
  expect_error(abstract_domain(4L, id = ""), class = "effect_input_error")
  expect_error(numerical_contract(-1, 1e-6), class = "effect_input_error")

  domain <- abstract_domain(4L, id = "conditions-entry")
  other <- abstract_domain(4L, id = "conditions-other")
  betas <- function() {
    value <- matrix(seq_len(8), 2L, 4L)
    rownames(value) <- c("a", "b")
    value
  }
  rel <- relation(list(run1 = betas(), run2 = betas()), domain = domain)
  expect_error(compile_frame(voxelwise(), 1L), class = "effect_input_error")

  # Two individually valid objects that were built against different domains
  # is a contract failure, not a shape error.
  expect_error(
    plan_geometry(rel, compile_frame(voxelwise(), other),
      cross_partitions(rel, independence = "independent")),
    class = "effect_contract_error")
})
