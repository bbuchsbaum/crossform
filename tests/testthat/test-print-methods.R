# Compact printing for user-facing objects.
#
# Every print method must produce a `<class>` header, stay inside one screen,
# and never leak a closure, an environment, a memory address, a whole digest,
# or a whole matrix. These tests assert the contract structurally rather than
# by snapshot, so they stay stable across BLAS and platform differences.

# Structural expectations ----------------------------------------------------

expect_compact_print <- function(object, class, max_lines = 20L) {
  output <- utils::capture.output(print(object))
  expect_gt(length(output), 0L)
  expect_identical(output[[1L]], paste0("<", class, ">"))
  expect_lte(length(output), max_lines)
  expect_false(any(grepl("<environment", output, fixed = TRUE)))
  expect_false(any(grepl("function (", output, fixed = TRUE)))
  expect_false(any(grepl("function(", output, fixed = TRUE)))
  expect_false(any(grepl("0x[0-9a-f]{6,}", output)))
  expect_false(any(grepl("[0-9a-f]{32,}", output)))
  expect_false(any(nchar(output) > 80L))
  invisible(output)
}

expect_compact_format <- function(object, class) {
  value <- format(object)
  expect_type(value, "character")
  expect_gte(length(value), 1L)
  expect_match(value[[1L]], paste0("^<", class), perl = TRUE)
  expect_false(any(grepl("[0-9a-f]{32,}", value)))
  invisible(value)
}

expect_prints_invisibly <- function(object) {
  result <- utils::capture.output(visible <- withVisible(print(object)))
  expect_false(visible$visible)
  expect_identical(visible$value, object)
  invisible(result)
}

# Fixtures -------------------------------------------------------------------
#
# Built once per file. Everything here is reachable through the public API.

print_fixture <- local({
  cached <- NULL
  function() {
    if (!is.null(cached)) {
      return(cached)
    }
    example <- example_fmri_effects()
    relation <- example$fit$relation
    pairing <- cross_partitions(relation, independence = "independent")
    plan <- plan_geometry(relation, example$frame, pairing)
    geometry <- materialize_geometry(plan)
    covariance <- rdm_sampling_covariance(plan, example$fit, target = "null")
    cached <<- list(
      example = example,
      relation = relation,
      pairing = pairing,
      plan = plan,
      geometry = geometry,
      covariance = covariance
    )
    cached
  }
})

coupling_fixture <- local({
  cached <- NULL
  function() {
    if (!is.null(cached)) {
      return(cached)
    }
    samples <- effect_space(paste0("time", seq_len(8)),
      basis_id = "print:time:v1")
    native <- abstract_domain(4, id = "print:native:v1")
    trend <- seq(-1.75, 1.75, length.out = 8)
    session <- function(shift) {
      cbind(
        trend + shift,
        -0.85 * trend + c(0.2, -0.1, 0.1, 0, -0.1, 0.1, -0.2, 0.1),
        sin(seq(shift, pi + shift, length.out = 8)),
        cos(seq(shift, pi + shift, length.out = 8))
      )
    }
    signals <- relation(
      list(session1 = session(0), session2 = session(0.1)),
      effects = samples, domain = native
    )
    nodes <- measurement_frame(
      list(anterior = matrix(c(1, 0, 0, 0), 1),
        posterior = matrix(c(0, 1, 0, 0), 1)),
      domain = native, id = "print:regional-means:v1"
    )
    node_pairs <- expand.grid(from = c("anterior", "posterior"),
      to = c("anterior", "posterior"), stringsAsFactors = FALSE)
    edges <- edge_frame(node_pairs$from, node_pairs$to, nodes)
    query <- variation_query(
      (diag(8) - matrix(1 / 8, 8, 8)) / 7, samples,
      sampling_axis = "time", construction = "joint_covariance",
      provenance = list(estimator = "centered within session")
    )
    products <- pairing(signals$partitions, signals$partitions,
      directed = TRUE, self_pairs = "allow_biased",
      independence = "not_independent")
    form <- measurement_form(left = signals, between = edges, by = query,
      over = products)
    cached <<- list(
      nodes = nodes, edges = edges, query = query, form = form,
      coupling = effect_coupling(form)
    )
    cached
  }
})

refusal_fixture <- function() {
  domain <- abstract_domain(2, id = "print:refusal:v1")
  relation <- relation(
    list(a = matrix(1:4, 2), b = matrix(2:5, 2)),
    effects = c("x", "y"), domain = domain
  )
  catch_refusal(rdm_sampling_covariance(
    plan_geometry(relation, compile_frame(whole_brain(), domain),
      cross_partitions(relation, independence = "independent")),
    relation, target = "null"
  ))
}

# The worked example ---------------------------------------------------------

test_that("printing the worked example stays inside one screen", {
  example <- example_fmri_effects()
  output <- utils::capture.output(print(example))
  expect_lt(length(output), 120L)
  expect_false(any(grepl("<environment", output, fixed = TRUE)))
  expect_false(any(grepl("function(", output, fixed = TRUE)))
  expect_compact_print(example, "effect_example_effects")
  expect_prints_invisibly(example)
  expect_match(paste(output, collapse = "\n"), "plan_geometry")
})

test_that("the worked example names its parts without expanding them", {
  output <- utils::capture.output(print(example_fmri_effects()))
  expect_match(paste(output, collapse = "\n"), "planted features")
  expect_match(paste(output, collapse = "\n"), "4 x 4 model dissimilarities")
  expect_false(any(grepl("^\\s*\\[1\\]", output)))
})

# Capability refusals --------------------------------------------------------

test_that("a capability refusal prints its capability, reasons, and remedies", {
  refusal <- refusal_fixture()
  expect_s3_class(refusal, "effect_capability_refusal")
  output <- expect_compact_print(refusal, "effect_capability_refusal")
  joined <- paste(output, collapse = "\n")
  expect_match(joined, "capability:\\s+sampling_covariance")
  expect_match(joined, "namespace:\\s+evidence_sampling")
  expect_match(joined, "reasons:")
  expect_match(joined, "remedies:")
  expect_match(joined, "- missing_error_channel")
  expect_match(joined, "lm_relation_fit")
  expect_prints_invisibly(refusal)
})

test_that("format() of a refusal returns one element per line", {
  refusal <- refusal_fixture()
  value <- format(refusal)
  expect_type(value, "character")
  expect_gt(length(value), 5L)
  expect_identical(value[[1L]], "<effect_capability_refusal>")
  expect_false(any(grepl("\n", value, fixed = TRUE)))
  expect_identical(value, utils::capture.output(print(refusal)))
})

test_that("a refusal with no reasons or remedies still prints", {
  bare <- structure(
    list(message = "refused", call = NULL, capability = "demo",
      namespace = "demo_namespace", reasons = character(),
      remedies = character()),
    class = c("effect_capability_refusal", "error", "condition")
  )
  output <- expect_compact_print(bare, "effect_capability_refusal")
  expect_match(paste(output, collapse = "\n"), "none recorded")
})

# Capability records ---------------------------------------------------------

test_that("capability records share one available/unavailable rendering", {
  fixture <- print_fixture()
  source_capability <- fixture$relation$capabilities[[1L]]
  error_capability <- fixture$example$fit$capabilities[[1L]]
  metric_capability <- metric_capabilities(
    identity_metric(fixture$example$domain))

  expect_compact_print(source_capability, "effect_source_capabilities")
  expect_compact_print(error_capability, "effect_error_capabilities")
  expect_compact_print(metric_capability, "effect_metric_capabilities")

  expect_match(format(source_capability), "granted")
  expect_match(format(error_capability), "granted")
  expect_match(format(metric_capability), "granted")

  expect_match(
    paste(utils::capture.output(print(source_capability)), collapse = "\n"),
    "block_read:\\s+yes"
  )
  expect_match(
    paste(utils::capture.output(print(metric_capability)), collapse = "\n"),
    "granted:"
  )
})

test_that("an adjudicated capability record reports its first reason", {
  unavailable <- structure(list(
    available = FALSE,
    capabilities = list(sampling_covariance = "unavailable",
      error_channel = "none"),
    reasons = data.frame(
      reason = "missing_error_channel",
      why = "the relation carries no residual channel",
      remedy = "refit with lm_relation_fit()",
      stringsAsFactors = FALSE
    )
  ), class = "effect_measurement_capabilities")
  output <- expect_compact_print(unavailable,
    "effect_measurement_capabilities")
  joined <- paste(output, collapse = "\n")
  expect_match(joined, "available:\\s+no")
  expect_match(joined, "first reason:")
  expect_match(joined, "missing_error_channel")
  expect_match(format(unavailable), "unavailable")
})

# Structural classes ---------------------------------------------------------

test_that("domains, frames, and relations print compactly", {
  fixture <- print_fixture()
  cases <- list(
    list(fixture$example$domain, "effect_domain"),
    list(fixture$example$domain$reference, "effect_domain_reference"),
    list(fixture$example$frame, "effect_frame"),
    list(fixture$example$frame$specification, "effect_frame_spec"),
    list(fixture$example$frame$support_index, "effect_support_index"),
    list(fixture$example$frame$support_index$cost, "effect_support_cost"),
    list(fixture$relation, "effect_relation"),
    list(fixture$relation$effect_space, "effect_space"),
    list(fixture$relation$extractors[[1L]], "effect_extractor"),
    list(fixture$relation$sources[[1L]], "effect_response_source"),
    list(fixture$relation$sources[[1L]]$descriptor,
      "effect_source_descriptor"),
    list(fixture$example$fit$error_models[[1L]], "effect_error_model"),
    list(fixture$example$fit$error_models[[1L]]$residual_source,
      "effect_residual_source")
  )
  for (case in cases) {
    expect_compact_print(case[[1L]], case[[2L]])
    expect_compact_format(case[[1L]], case[[2L]])
    expect_prints_invisibly(case[[1L]])
  }
})

test_that("domain and frame prints name the next call to make", {
  fixture <- print_fixture()
  domain_output <- paste(
    utils::capture.output(print(fixture$example$domain)), collapse = "\n")
  frame_output <- paste(
    utils::capture.output(print(fixture$example$frame)), collapse = "\n")
  relation_output <- paste(
    utils::capture.output(print(fixture$relation)), collapse = "\n")
  expect_match(domain_output, "compile_frame")
  expect_match(frame_output, "plan_geometry")
  expect_match(relation_output, "plan_geometry")
  expect_match(relation_output, "lazy")
})

test_that("a relation print never reads or shows its sources", {
  fixture <- print_fixture()
  output <- utils::capture.output(print(fixture$relation))
  expect_match(paste(output, collapse = "\n"), "unread")
  expect_false(any(grepl("[0-9]+\\.[0-9]{6,}", output)))
})

test_that("pairings print their contract rather than their rows", {
  fixture <- print_fixture()
  output <- expect_compact_print(fixture$pairing, "effect_pairing")
  joined <- paste(output, collapse = "\n")
  expect_match(joined, "independence:\\s+independent")
  expect_match(joined, "self pairs:\\s+forbid")
  expect_match(joined, "pairs:\\s+6")
  expect_prints_invisibly(fixture$pairing)
})

# Plans, policies, and schedules ---------------------------------------------

test_that("plans, schedules, and policies print compactly", {
  fixture <- print_fixture()
  domain <- fixture$example$domain
  cases <- list(
    list(fixture$plan$task, "effect_evidence_task"),
    list(fixture$plan$task$stages, "effect_evidence_stages"),
    list(fixture$plan$task$ordered_partition_products,
      "effect_ordered_edges"),
    list(fixture$plan$metric_schedule, "effect_metric_schedule"),
    list(fixture$plan$compute, "effect_compute_policy"),
    list(identity_metric(domain), "effect_metric_recipe"),
    list(neural_metric(diag(domain$n_features), domain),
      "effect_neural_metric"),
    list(metric_training_policy("exclude_evaluation"),
      "effect_metric_training_policy"),
    list(numerical_contract(), "effect_numerical_contract"),
    list(numerical_agreement(c(1, 2), c(1, 2)), "effect_numeric_agreement"),
    list(reduce_partitions(), "effect_partition_reducer"),
    list(gaussian_covariance_model(), "effect_gaussian_covariance_model"),
    list(frame_conservation(fixture$example$frame),
      "effect_frame_conservation")
  )
  for (case in cases) {
    expect_compact_print(case[[1L]], case[[2L]])
    expect_prints_invisibly(case[[1L]])
  }
  # `effect_ordered_edges` subclasses data.frame, so only print() is defined.
  for (case in cases[-3L]) {
    expect_compact_format(case[[1L]], case[[2L]])
  }
})

test_that("numbers print at four significant digits, never full precision", {
  agreement <- numerical_agreement(c(1, 2), c(1, 2) + 1e-12)
  output <- utils::capture.output(print(agreement))
  expect_false(any(grepl("[0-9]\\.[0-9]{7,}", output)))
  contract <- utils::capture.output(print(numerical_contract()))
  expect_false(any(grepl("[0-9]\\.[0-9]{7,}", contract)))
})

test_that("a metric never prints its matrix", {
  fixture <- print_fixture()
  domain <- fixture$example$domain
  output <- utils::capture.output(
    print(neural_metric(diag(domain$n_features), domain))
  )
  expect_lte(length(output), 20L)
  expect_match(paste(output, collapse = "\n"), "245 x 245 \\(not shown\\)")
})

# Results, receipts, and sampling --------------------------------------------

test_that("results, receipts, and sampling records print compactly", {
  fixture <- print_fixture()
  cases <- list(
    list(fixture$geometry$receipt, "effect_execution_receipt"),
    list(fixture$geometry$total, "effect_geometry_store"),
    list(fixture$covariance$plan, "effect_evidence_sampling_plan"),
    list(fixture$covariance$plan$evidence, "effect_sampling_record"),
    list(fixture$geometry$receipt$memory, "effect_memory_plan"),
    list(residual_pair_statistics(fixture$example$fit, fixture$example$frame),
      "effect_residual_pair_statistics")
  )
  for (case in cases) {
    expect_compact_print(case[[1L]], case[[2L]])
    expect_compact_format(case[[1L]], case[[2L]])
    expect_prints_invisibly(case[[1L]])
  }
})

test_that("an execution receipt prints no platform-dependent identity", {
  fixture <- print_fixture()
  output <- utils::capture.output(print(fixture$geometry$receipt))
  # No absolute filesystem path (the BLAS vendor string is deliberately
  # omitted; it stays reachable as `receipt$blas$vendor`).
  expect_false(any(grepl("(^|\\s)/[A-Za-z]", output)))
  expect_false(any(grepl("[Ff]ramework|dylib|[.]so\\b", output)))
  expect_match(paste(output, collapse = "\n"), "status:\\s+complete")
})

test_that("the base sampling covariance class has its own compact print", {
  fixture <- print_fixture()
  base <- fixture$covariance
  class(base) <- "effect_sampling_covariance"
  expect_compact_print(base, "effect_sampling_covariance")
  expect_compact_format(base, "effect_sampling_covariance")
})

# Study facts ----------------------------------------------------------------

test_that("study facts print compactly", {
  fixture <- bound_study_fixture()
  cases <- list(
    list(fixture$observations, "effect_observations"),
    list(fixture$fixture$indexes[[1L]], "effect_observation_index"),
    list(fixture$events, "effect_events"),
    list(fixture$confounds, "effect_observation_confounds"),
    list(fixture$hierarchy, "effect_partition_hierarchy")
  )
  for (case in cases) {
    expect_compact_print(case[[1L]], case[[2L]])
    expect_compact_format(case[[1L]], case[[2L]])
    expect_prints_invisibly(case[[1L]])
  }
})

test_that("study facts report their counts and censoring", {
  fixture <- bound_study_fixture()
  observations <- paste(
    utils::capture.output(print(fixture$observations)), collapse = "\n")
  confounds <- paste(
    utils::capture.output(print(fixture$confounds)), collapse = "\n")
  hierarchy <- paste(
    utils::capture.output(print(fixture$hierarchy)), collapse = "\n")
  expect_match(observations, "16 total")
  expect_match(confounds, "rows retained|dropped")
  expect_match(hierarchy, "subject")
})

# Design models --------------------------------------------------------------

test_that("design models and receipts print compactly", {
  fixture <- relation_plan_fixture()
  receipt <- relation_plan_receipts(fixture$plan)[[1L]]
  cases <- list(
    list(fixture$model, "effect_design_model"),
    list(receipt, "effect_design_receipt"),
    list(fixture$observation, "effect_observation_model")
  )
  for (case in cases) {
    expect_compact_print(case[[1L]], case[[2L]])
    expect_compact_format(case[[1L]], case[[2L]])
    expect_prints_invisibly(case[[1L]])
  }
  expect_match(
    paste(utils::capture.output(print(receipt)), collapse = "\n"),
    "rows retained"
  )
})

# Measurement forms and coupling ---------------------------------------------

test_that("measurement forms and coupling results print compactly", {
  fixture <- coupling_fixture()
  cases <- list(
    list(fixture$nodes, "effect_measurement_frame"),
    list(fixture$nodes$legs[[1L]], "effect_measurement_leg"),
    list(fixture$nodes$legs[[1L]]$output_space, "effect_measurement_axis"),
    list(fixture$edges, "effect_edge_frame"),
    list(fixture$edges$edges, "effect_measurement_edges"),
    list(fixture$query, "effect_pair_query"),
    list(fixture$form, "effect_measurement_form"),
    list(fixture$form$plan, "effect_measurement_plan"),
    list(fixture$form$plan$regularization,
      "effect_measurement_regularization"),
    list(fixture$form$store, "effect_measurement_store"),
    list(fixture$form$receipt, "effect_measurement_receipt"),
    list(fixture$form$diagnostics, "effect_measurement_diagnostics"),
    list(fixture$coupling, "effect_coupling_result"),
    list(measurement_space(3L, id = "print:measurements:v1"),
      "effect_measurement_space")
  )
  for (case in cases) {
    expect_compact_print(case[[1L]], case[[2L]])
    expect_compact_format(case[[1L]], case[[2L]])
    expect_prints_invisibly(case[[1L]])
  }
})

test_that("block indices print their shape, not their rows", {
  fixture <- coupling_fixture()
  output <- expect_compact_print(fixture$form$block_index,
    "effect_measurement_block_index")
  expect_match(paste(output, collapse = "\n"), "blocks:\\s+4")
  expect_prints_invisibly(fixture$form$block_index)
})

test_that("a measurement form states its claims and its next call", {
  fixture <- coupling_fixture()
  joined <- paste(utils::capture.output(print(fixture$form)), collapse = "\n")
  expect_match(joined, "symmetric")
  expect_match(joined, "effect_coupling\\(form\\)")
  coupling_output <- paste(
    utils::capture.output(print(fixture$coupling)), collapse = "\n")
  expect_match(coupling_output, "kind:\\s+effect_coupling")
  expect_match(coupling_output, "units:\\s+not claimed")
})

# Remaining constructible classes --------------------------------------------

test_that("crossnobis and raw-design objects print compactly", {
  fixture <- print_fixture()
  frozen <- plan_crossnobis(fixture$example$fit, at = fixture$example$frame,
    over = fixture$pairing)$metric_schedule
  design <- matrix(stats::rnorm(12L), 6L, 2L,
    dimnames = list(NULL, c("c1", "c2")))
  raw <- raw_design_model(list(a = design), list(a = seq_len(6L)), "qr")
  cases <- list(
    list(frozen, "effect_frozen_metric_schedule"),
    list(raw, "effect_raw_design_model")
  )
  for (case in cases) {
    expect_compact_print(case[[1L]], case[[2L]])
    expect_compact_format(case[[1L]], case[[2L]])
    expect_prints_invisibly(case[[1L]])
  }
  expect_match(
    paste(utils::capture.output(print(frozen)), collapse = "\n"),
    "frozen"
  )
  expect_match(
    paste(utils::capture.output(print(raw)), collapse = "\n"),
    "not claimed"
  )
})

test_that("queries and pair couplings print compactly", {
  fixture <- print_fixture()
  space <- fixture$relation$effect_space
  query <- bilinear_query(diag(4), effects = fixture$relation$effects)
  matches <- match_coupling(c("face", "body"), c("house", "tool"),
    space, space)
  controls <- control_coupling(matches)
  cases <- list(
    list(query, "effect_query"),
    list(matches, "effect_pair_coupling"),
    list(controls, "effect_pair_coupling")
  )
  for (case in cases) {
    expect_compact_print(case[[1L]], case[[2L]])
    expect_compact_format(case[[1L]], case[[2L]])
    expect_prints_invisibly(case[[1L]])
  }
  expect_match(
    paste(utils::capture.output(print(matches)), collapse = "\n"),
    "kind:\\s+match"
  )
})

test_that("hand-built value records print without a fixture", {
  # These classes are produced deep inside the kernel; the print methods are
  # pure formatters, so they are exercised against their documented shape.
  lowering <- structure(
    list(kind = "bilinear", collapsed = FALSE, reason = "rank above one"),
    class = "effect_lowering"
  )
  whitener <- structure(
    list(kind = "identity", dim = c(8L, 8L), signature = "sha256:abc123"),
    class = "effect_observation_whitener"
  )
  handle <- structure(
    list(metric = list(role = "same_space_metric"), diagnostics = list(),
      inverse_mode = "cholesky", factorizations = list(1L), shortcut = TRUE),
    class = "effect_metric_handle"
  )
  components <- structure(
    list(coherent_rank = 2L, configuration_psd = TRUE,
      denominator = "trace", inverse_quadratic_mode = "solve",
      factorization_count = 1L, signature = "sha256:aaa111"),
    class = "effect_metric_components"
  )
  inverse <- structure(list(kind = "none", signature = "sha256:bbb222"),
    class = "effect_metric_inverse_representation")
  for (object in list(lowering, whitener, handle, components, inverse)) {
    expect_compact_print(object, class(object)[[1L]])
    expect_compact_format(object, class(object)[[1L]])
    expect_prints_invisibly(object)
  }
})

test_that("a coupling partition policy prints its placement and weights", {
  policy <- crossform:::.coupling_partition_policy(
    c(0.5, 0.5), "within_partition_pair",
    structure(list(schema_version = 1L, kind = "identity", boundary = NULL,
      delta = NULL, ties = NULL), class = "effect_edge_transform")
  )
  expect_compact_print(policy, "effect_coupling_partition_policy")
  expect_compact_format(policy, "effect_coupling_partition_policy")
  expect_prints_invisibly(policy)
})

# Contrast views -------------------------------------------------------------

test_that("a contrast view prints the alignment its weights actually got", {
  fixture <- print_fixture()
  # Positional weights are accepted, so the print must show which effect each
  # weight landed on. A reader who mis-ordered them sees it immediately.
  view <- contrast_energy(fixture$plan, c(0.5, 0.5, -0.5, -0.5))
  output <- utils::capture.output(print(view))
  expect_identical(output[[1L]], "<effect_contrast_view>")
  expect_match(output[[3L]], "^  contrast: ")
  joined <- paste(output, collapse = "\n")
  for (effect in fixture$relation$effects) {
    expect_match(joined, effect, fixed = TRUE)
  }
  expect_match(joined, "face 0.5", fixed = TRUE)
  expect_match(joined, "tool -0.5", fixed = TRUE)
  # The header block stays aligned: `measurements` is the longest key.
  expect_identical(substr(output[[2L]], 1L, 16L), "  measurements: ")
  expect_identical(substr(output[[3L]], 1L, 16L), "  contrast:     ")
  expect_prints_invisibly(view)
})

test_that("named contrast weights print in relation order", {
  fixture <- print_fixture()
  view <- contrast_energy(fixture$plan,
    c(tool = -0.5, face = 0.5, house = -0.5, body = 0.5))
  contrast <- utils::capture.output(print(view))[[3L]]
  expect_match(contrast, "face 0.5, body 0.5, house -0.5, tool -0.5",
    fixed = TRUE)
})

test_that("the weight formatter caps long contrasts and handles odd input", {
  format_weights <- crossform:::.pf_weights
  expect_identical(format_weights(NULL), "none")
  expect_identical(format_weights(numeric()), "none")
  expect_identical(format_weights(c(a = 1, b = -1)), "a 1, b -1")
  # Unnamed weights fall back to positions rather than printing bare numbers.
  expect_identical(format_weights(c(1, -1)), "[1] 1, [2] -1")
  long <- stats::setNames(seq_len(9L), paste0("e", seq_len(9L)))
  expect_match(format_weights(long), "\\(\\+3 more\\)$")
  # Four significant digits, never LAPACK-length floats.
  expect_identical(format_weights(c(a = 1 / 3)), "a 0.3333")
})

# Small value records --------------------------------------------------------

test_that("small value records print compactly", {
  fixture <- print_fixture()
  domain <- fixture$example$domain
  recipe <- identity_metric(domain)
  cases <- list(
    list(neural_metric(diag(domain$n_features),
      domain)$inverse_representation,
      "effect_metric_inverse_representation"),
    list(recipe$capabilities, "effect_metric_capabilities")
  )
  for (case in cases) {
    expect_compact_print(case[[1L]], case[[2L]])
    expect_compact_format(case[[1L]], case[[2L]])
  }
})

test_that("every print method in the package returns its input invisibly", {
  namespace <- asNamespace("crossform")
  source_text <- function(name) {
    paste(deparse(body(get(name, envir = namespace))), collapse = " ")
  }
  # A method returns invisibly either directly or through the one shared
  # preview helper it delegates to.
  returns_invisibly <- function(name) {
    text <- source_text(name)
    if (grepl("invisible", text, fixed = TRUE)) {
      return(TRUE)
    }
    called <- unique(sub("[(]$", "",
      regmatches(text, gregexpr("[.][A-Za-z_.]+[(]", text))[[1L]]))
    called <- called[vapply(called, exists, logical(1),
      envir = namespace, inherits = FALSE)]
    any(vapply(called, function(helper) {
      grepl("invisible", source_text(helper), fixed = TRUE)
    }, logical(1)))
  }
  methods <- ls(namespace, pattern = "^print[.]effect_")
  expect_gt(length(methods), 80L)
  expect_setequal(methods[!vapply(methods, returns_invisibly, logical(1))],
    character())
})

test_that("the package registers a print method for every documented class", {
  namespace <- asNamespace("crossform")
  printed <- sub("^print[.]", "", ls(namespace, pattern = "^print[.]effect_"))
  # Regression guard: the classes a first-hour user is most likely to hold.
  expected <- c(
    "effect_relation", "effect_domain", "effect_frame", "effect_frame_spec",
    "effect_pairing", "effect_neural_metric", "effect_measurement_form",
    "effect_coupling_result", "effect_observations", "effect_events",
    "effect_observation_index", "effect_observation_confounds",
    "effect_partition_hierarchy", "effect_sampling_covariance",
    "effect_execution_receipt", "effect_design_receipt",
    "effect_design_model", "effect_evidence_sampling_plan",
    "effect_space", "effect_support_index", "effect_metric_schedule",
    "effect_metric_training_policy", "effect_compute_policy",
    "effect_numerical_contract", "effect_numeric_agreement",
    "effect_extractor", "effect_capability_refusal",
    "effect_source_capabilities", "effect_error_capabilities",
    "effect_metric_capabilities", "effect_measurement_capabilities"
  )
  expect_setequal(setdiff(expected, printed), character())
})
