# `materialize_population()` --- the streamed complete-form group executor of
# `population-form-v1`.
#
# E4 reads `K` declared contrasts; this route reads the whole `q`-by-`q` form at
# every group node without ever holding the `N`-by-`(m+1)`-by-`p` array that
# would be. Three things have to be true for that to be worth anything, and
# each has its own test below.
#
#   1. It is the same number the naive dense route gets. The oracle is written
#      out here in the *other* evaluation order --- materialize each
#      participant's complete packed geometry, transport all `p` coordinates,
#      fit every one of them --- out of the package's public verbs, touching
#      nothing this file is testing.
#   2. The assembled form answers a contrast the way `estimate_population()`
#      answers it. This is the commutation claim of `population-form-v1`
#      section 3 read across the two executors, and it is the test that makes
#      the form an object rather than a pile of coefficients: a form you cannot
#      query is not a form. It rests on the `sqrt(2)` off-diagonal packing of
#      `symmetric_packed` (section 5), so it is checked both through the packed
#      inner product and through the dense Frobenius one.
#   3. The tile bounds memory and changes nothing else.
#
# The fixture is E4's, deliberately: four participants on four different native
# frame sizes, a transport radius that leaves real sink mass, three effects so
# the packed width (six) is larger than any bank a reader would write.

pf_effects <- function() {
  effect_space(c("face", "house", "tool"), basis_id = "pop-form:v1")
}

pf_subject <- function(id, features, gain = 1, effects = pf_effects()) {
  q <- length(effects$coordinates)
  domain <- abstract_domain(features,
    coordinates = cbind(x = seq_len(features) - 1),
    feature_ids = paste0("f", seq_len(features)), id = id)
  # Signed by construction: a conservative ledger is signed, and a fixture of
  # nonnegative numbers would let a mass-law bug pass the budget certificate.
  values <- function(divisor) {
    set.seed(sum(as.integer(charToRaw(id))) + features +
      round(1000 * divisor))
    matrix(gain * stats::rnorm(q * features) / divisor, q, features,
      dimnames = list(effects$coordinates, NULL))
  }
  relation <- relation(
    list(run1 = values(1), run2 = values(1.3), run3 = values(0.8)),
    effects = effects, domain = domain
  )
  plan_geometry(relation, compile_frame(voxelwise(), domain),
    cross_partitions(relation))
}

# Group centres at 0, 4 and 9 with radius 1.5, so native nodes at x = 2, 6, 7
# or 11 land entirely in the sink and the sink column carries real mass. Every
# participant reaches every ordinary node so this file isolates streamed-form
# algebra; variable coverage has its own target tests.
pf_carrier <- function(features, semantics = "budget", radius = 1.5, ...) {
  anatomical_transport(
    native_coords = cbind(seq_len(features) - 1),
    group_coords = cbind(c(0, 4, 9)),
    semantics = semantics, radius = radius, ...
  )
}

pf_sizes <- c(s01 = 10L, s02 = 11L, s03 = 12L, s04 = 13L)
pf_gains <- c(s01 = 1, s02 = 1.6, s03 = 0.6, s04 = 1.2)

pf_plan <- function(sizes = pf_sizes, semantics = "budget", ...) {
  subjects <- stats::setNames(lapply(names(sizes), function(id)
    pf_subject(id, sizes[[id]], pf_gains[[id]])), names(sizes))
  plan_population(
    subjects,
    lapply(stats::setNames(names(sizes), names(sizes)), function(id)
      pf_carrier(sizes[[id]], semantics = semantics)),
    ...
  )
}

pf_packed <- function(bank) {
  width <- ncol(bank) * (ncol(bank) + 1L) / 2L
  matrix(vapply(seq_len(nrow(bank)), function(k) {
    crossform:::.svec_symmetric(tcrossprod(bank[k, ]))
  }, numeric(width)), nrow = width, dimnames = list(NULL, rownames(bank)))
}

# The dense oracle: transport every packed coordinate, then fit every one of
# them. This is the `N`-by-`(m+1)`-by-`p` array the streamed route refuses to
# build, built here on purpose so there is something to be equal to. Nothing
# below calls `materialize_population()`.
pf_dense_oracle <- function(plan, component = "total") {
  labels <- names(plan$subjects)
  nodes <- nrow(plan$group_index) + 1L
  width <- plan$subjects[[1L]]$packed_width
  stack <- matrix(NA_real_, length(labels), nodes * width,
    dimnames = list(labels, NULL))
  for (label in labels) {
    rows <- geometry_component(
      materialize_geometry(plan$subjects[[label]]), component
    )
    stack[label, ] <- as.numeric(transport_values(plan$transport[[label]], rows))
  }
  terms <- plan$model$columns
  finite <- apply(is.finite(stack), 2L, all)
  coefficients <- matrix(NA_real_, length(terms), ncol(stack))
  coefficients[, finite] <- qr.coef(plan$model$qr, stack[, finite, drop = FALSE])
  # `node x coordinate x term`, then permuted to the result's own axis order.
  dense <- array(t(coefficients), dim = c(nodes, width, length(terms)))
  aperm(dense, c(1L, 3L, 2L))
}


# The acceptance oracle --------------------------------------------------------

test_that("the streamed full-form executor equals the dense oracle", {
  plan <- pf_plan()
  form <- materialize_population(plan)
  oracle <- pf_dense_oracle(plan)

  expect_identical(dim(form$coefficient_forms), dim(oracle))
  expect_identical(names(dimnames(form$coefficient_forms)),
    c("node", "term", "coordinate"))
  expect_lt(max(abs(form$coefficient_forms - oracle)), 1e-12)

  # Not vacuously zero: the oracle carries real signal at every coordinate.
  expect_gt(min(apply(abs(oracle), 3L, max)), 1e-6)
})

test_that("the oracle agreement holds for the coherent and configuration ledgers", {
  plan <- pf_plan()
  for (component in c("coherent", "configuration")) {
    form <- materialize_population(plan, component = component)
    expect_lt(
      max(abs(form$coefficient_forms - pf_dense_oracle(plan, component))),
      1e-12
    )
  }
})

test_that("the additive ledger identity survives transport and the group fit", {
  plan <- pf_plan()
  total <- materialize_population(plan, component = "total")
  coherent <- materialize_population(plan, component = "coherent")
  configuration <- materialize_population(plan, component = "configuration")

  expect_lt(max(abs(
    coherent$coefficient_forms + configuration$coefficient_forms -
      total$coefficient_forms
  )), 1e-12)
  expect_identical(
    c(total$ledger, coherent$ledger, configuration$ledger),
    c("transported_total", "native_coherent_ledger",
      "native_configuration_ledger")
  )
})

test_that("the budget certificate closes on the transported total", {
  plan <- pf_plan()
  form <- materialize_population(plan, component = "total")

  budget <- form$receipt$budget
  expect_true(budget$asserted)
  expect_identical(budget$scale, "relative_to_ledger_l1_norm")
  expect_lt(budget$max_relative_deviation, 1e-12)

  # And it is a real closure, not an artefact of an empty sink: the sink
  # carries mass at every packed coordinate for at least one participant, and
  # the transported total including the sink is the native total.
  expect_gt(max(abs(form$receipt$sink_budget)), 1e-6)
  for (label in names(plan$subjects)) {
    rows <- geometry_component(
      materialize_geometry(plan$subjects[[label]]), "total"
    )
    carried <- transport_values(plan$transport[[label]], rows)
    expect_lt(
      max(abs(colSums(carried) - form$receipt$native_total[label, ])), 1e-12
    )
  }
})


# Commutation across the two executors -----------------------------------------

test_that("querying the assembled form equals the query-first executor", {
  plan <- pf_plan()
  bank <- rbind(`face-house` = c(1, -1, 0), `face-tool` = c(1, 0, -1),
    `all` = c(1, 1, 1))
  form <- materialize_population(plan)
  queried <- estimate_population(plan, bank)
  packed <- pf_packed(bank)
  terms <- dimnames(form$coefficient_forms)[[2L]]

  # (a) through the packed inner product: the assembled `(m+1)`-by-`p` block
  # per term times the packed query matrix.
  for (term in terms) {
    assembled <- form$coefficient_forms[, term, ] %*% packed
    expect_lt(
      max(abs(assembled - queried$coefficients[, , term, drop = TRUE])), 1e-12
    )
  }

  # (b) through the dense Frobenius inner product, which is the same number
  # only because `symmetric_packed` carries `sqrt(2)` off the diagonal
  # (`population-form-v1` section 5: a naive packing is off by an O(1) amount,
  # not by a tolerance).
  q <- length(form$effects)
  for (k in seq_len(nrow(bank))) {
    operator <- tcrossprod(bank[k, ])
    for (term in terms) {
      frobenius <- vapply(seq_len(nrow(form$index)), function(node) {
        sum(operator * crossform:::.unsvec_symmetric(
          form$coefficient_forms[node, term, ], q
        ))
      }, numeric(1))
      expect_lt(
        max(abs(frobenius - queried$coefficients[, k, term])), 1e-12
      )
    }
  }

  # Not vacuous: the queried coefficients are far from zero.
  expect_gt(max(abs(queried$coefficients)), 1e-3)
})

test_that("the commutation survives density semantics and its NA group node", {
  # Group node 9 is out of reach of the two smaller participants under any
  # radius, so density returns NA there rather than 0, and both executors must
  # solve around the same columns.
  partial_sizes <- c(s01 = 6L, s02 = 8L, s03 = 10L, s04 = 12L)
  plan <- pf_plan(sizes = partial_sizes, semantics = "density")
  bank <- rbind(`face-house` = c(1, -1, 0))
  form <- materialize_population(plan)
  queried <- estimate_population(plan, bank)
  packed <- pf_packed(bank)

  assembled <- form$coefficient_forms[, "(Intercept)", ] %*% packed
  expect_true(anyNA(assembled))
  expect_identical(is.na(assembled[, 1L]),
    is.na(queried$coefficients[, 1L, "(Intercept)"]))
  finite <- !is.na(assembled[, 1L])
  expect_lt(max(abs(
    assembled[finite, 1L] - queried$coefficients[finite, 1L, "(Intercept)"]
  )), 1e-12)

  expect_false(form$receipt$budget$asserted)
  expect_gt(form$receipt$unresolved_columns, 0L)
})


# The tile ---------------------------------------------------------------------

test_that("two coordinate tiles give identical coefficient forms", {
  plan <- pf_plan()
  width <- plan$subjects[[1L]]$packed_width
  wide <- materialize_population(plan, coordinate_tile = width)
  narrow <- materialize_population(plan, coordinate_tile = 1L)
  middling <- materialize_population(plan, coordinate_tile = 4L)

  expect_identical(wide$receipt$streaming$passes_per_subject, 1L)
  expect_identical(narrow$receipt$streaming$passes_per_subject, width)
  expect_identical(middling$receipt$streaming$passes_per_subject, 2L)

  # The tile is a blocking choice, and `numerical_contract()` declares
  # `bitwise_across_blocking = FALSE` with a `block_partition` guarantee of
  # "tolerance": the tile changes the column count of the fused contraction, so
  # BLAS may reassociate one sum. Asserting bit-identity here would be
  # asserting something the package's own contract says it does not promise.
  # `atol` is what it does promise, and what these agree to is ~1 ulp.
  atol <- numerical_contract()$atol
  expect_lt(max(abs(wide$coefficient_forms - narrow$coefficient_forms)), atol)
  expect_lt(max(abs(wide$coefficient_forms - middling$coefficient_forms)), atol)
  # Not vacuous: the forms are far from zero, so agreement at `atol` is a
  # statement about the route rather than about a table of zeros.
  expect_gt(max(abs(wide$coefficient_forms)), 1e-3)

  # The tile is physical, so it is not in the estimand.
  expect_identical(wide$scientific_plan_id, narrow$scientific_plan_id)
  expect_identical(wide$receipt$native_total, narrow$receipt$native_total)
  expect_identical(wide$receipt$sink_budget, narrow$receipt$sink_budget)
})

test_that("the tile comes from the compute policy when the caller names none", {
  plan <- pf_plan()
  width <- plan$subjects[[1L]]$packed_width
  # No declared workspace: the compiler's own default, `min(64, width)`.
  expect_identical(
    materialize_population(plan)$receipt$streaming$coordinate_tile,
    as.integer(min(64L, width))
  )

  # A declared workspace sizes the tile against the arrays that are linear in
  # it, and never returns zero.
  nodes <- nrow(plan$group_index) + 1L
  per_coordinate <- 8 * (length(plan$subjects) * nodes +
    max(vapply(plan$subjects, function(s) as.double(s$measurements),
      numeric(1))))
  tight <- plan_population(
    plan$subjects, plan$transport,
    compute = compute_policy(workspace_bytes = 3 * per_coordinate)
  )
  expect_identical(
    materialize_population(tight)$receipt$streaming$coordinate_tile, 3L
  )
  starved <- plan_population(
    plan$subjects, plan$transport,
    compute = compute_policy(workspace_bytes = per_coordinate / 10)
  )
  expect_identical(
    materialize_population(starved)$receipt$streaming$coordinate_tile, 1L
  )
})

test_that("peak vector memory tracks the coordinate tile", {
  skip_on_cran()
  # A case where the refused dense stack is worth refusing: five participants,
  # ten effects (fifty-five packed coordinates) and two thousand group nodes,
  # so `N * (m+1) * p` is 550,275 doubles --- 4.2 MiB --- while a tile of
  # eleven holds a fifth of that at a time.
  effects <- effect_space(paste0("e", 1:10), basis_id = "pop-form-mem:v1")
  sizes <- c(m01 = 90L, m02 = 100L, m03 = 110L, m04 = 120L, m05 = 130L)
  subjects <- stats::setNames(lapply(names(sizes), function(id)
    pf_subject(id, sizes[[id]], 1, effects)), names(sizes))
  group <- cbind(seq(0, 130, length.out = 2000))
  carriers <- lapply(sizes, function(n) anatomical_transport(
    native_coords = cbind(seq_len(n) - 1), group_coords = group,
    semantics = "budget"
  ))
  plan <- plan_population(subjects, carriers)

  peak_bytes <- function(expression) {
    gc(reset = TRUE, full = TRUE)
    before <- gc(full = TRUE)[["Vcells", "used"]]
    value <- force(expression)
    after <- gc(full = TRUE)[["Vcells", "max used"]]
    list(value = value, bytes = (after - before) * 8)
  }

  wide <- peak_bytes(materialize_population(plan, coordinate_tile = 55L))
  narrow <- peak_bytes(materialize_population(plan, coordinate_tile = 11L))

  # The numbers themselves, for the record.
  streaming <- wide$value$receipt$streaming
  expect_identical(streaming$packed_width, 55L)
  expect_identical(streaming$refused_dense_doubles, 5 * 2001 * 55)
  expect_identical(streaming$group_stack_doubles, 5 * 2001 * 55)
  expect_identical(
    narrow$value$receipt$streaming$group_stack_doubles, 5 * 2001 * 11
  )
  cat(sprintf(
    "\n  coordinate tile 55: %.2f MiB peak, 1 pass\n  coordinate tile 11: %.2f MiB peak, 5 passes\n",
    wide$bytes / 1024^2, narrow$bytes / 1024^2
  ))

  # The answer does not move, and the peak does: the fivefold reduction in the
  # group stack is 3.4 MiB, so a route whose peak did not fall by megabytes
  # would be holding the dense array somewhere.
  expect_identical(
    is.na(wide$value$coefficient_forms),
    is.na(narrow$value$coefficient_forms)
  )
  finite <- is.finite(wide$value$coefficient_forms) &
    is.finite(narrow$value$coefficient_forms)
  expect_true(any(finite))
  expect_lt(
    max(abs(wide$value$coefficient_forms[finite] -
      narrow$value$coefficient_forms[finite])),
    numerical_contract()$atol
  )
  expect_lt(narrow$bytes, wide$bytes - 2 * 1024^2)
})


# Refusals ---------------------------------------------------------------------

test_that("a complete form refuses a normalization that varies by coordinate", {
  plan <- pf_plan(normalization = "unit_budget")
  refusal <- catch_refusal(materialize_population(plan))
  expect_identical(refusal$capability, "complete_form_normalization")
  expect_identical(refusal$reasons, c(
    "normalization_varies_along_the_coordinate_axis",
    "declared_normalization:unit_budget"
  ))
  # The remedy is the other executor, where the divisor is defined.
  expect_true(any(grepl("estimate_population", refusal$remedies, fixed = TRUE)))

  # And `none` is admitted on the same plan shape.
  expect_s3_class(materialize_population(pf_plan()),
    "effect_population_result")
})

test_that("a group index carrying the result's own columns is refused", {
  sizes <- pf_sizes
  group <- data.frame(node = c("g1", "g2", "g3"), units = c("mm", "mm", "mm"),
    stringsAsFactors = FALSE)
  subjects <- stats::setNames(lapply(names(sizes), function(id)
    pf_subject(id, sizes[[id]], pf_gains[[id]])), names(sizes))
  plan <- plan_population(subjects,
    lapply(stats::setNames(names(sizes), names(sizes)), function(id)
      pf_carrier(sizes[[id]], group_index = group)))
  refusal <- catch_refusal(materialize_population(plan))
  expect_identical(refusal$capability, "reserved_group_index_columns")
})


# The record -------------------------------------------------------------------

test_that("the complete-form result is a sealed record with a stable contract", {
  plan <- pf_plan()
  form <- materialize_population(plan)

  expect_s3_class(form, "effect_population_result")
  expect_identical(form$basis, "complete_form")
  expect_silent(crossform:::.validate_population_result(form))
  expect_match(form$scientific_plan_id, "^population-sha256:")

  # The field contract the population views build on.
  expect_identical(names(form), c(
    "basis", "coefficient_forms", "residual_df", "index", "effects",
    "coordinates", "component", "ledger", "semantics", "normalization",
    "coverage", "uncertainty", "receipt", "scientific_plan_id"
  ))
  expect_null(form$values)
  expect_null(form$fitted)
  expect_null(form$residuals)

  # The packed coordinate table, and the `sqrt(2)` that makes the packed inner
  # product the Frobenius one.
  expect_identical(names(form$coordinates),
    c("coordinate", "row", "column", "scale"))
  expect_identical(nrow(form$coordinates), 6L)
  expect_identical(form$coordinates$coordinate[[1L]], "face:face")
  expect_identical(form$coordinates$scale,
    c(1, sqrt(2), sqrt(2), 1, sqrt(2), 1))
  expect_identical(form$effects, c("face", "house", "tool"))

  # Every `[node, term, ]` slice unpacks to a symmetric form.
  for (node in seq_len(nrow(form$index))) {
    unpacked <- crossform:::.unsvec_symmetric(
      form$coefficient_forms[node, "(Intercept)", ], 3L
    )
    expect_identical(unpacked, t(unpacked))
  }

  # The sink is a row of the result, materialized and marked.
  expect_identical(sum(form$index$sink), 1L)
  expect_identical(as.character(form$index$node[form$index$sink]), "<sink>")

  receipt <- form$receipt
  expect_identical(receipt$contract_version, "population-form-v1")
  expect_identical(receipt$basis, "complete_form")
  expect_identical(receipt$readout$width, 6L)
  expect_null(receipt$queries)
  expect_identical(receipt$population_plan_id, plan$scientific_plan_id)
  expect_identical(names(receipt$subject_receipts), names(plan$subjects))
  for (value in receipt$subject_receipts) {
    expect_s3_class(value, "effect_execution_receipt")
  }
  expect_identical(dim(receipt$native_total), c(4L, 6L))
  expect_identical(dim(receipt$sink_budget), c(4L, 6L))
  expect_identical(receipt$streaming$output_doubles, 4 * 1 * 6)
})

test_that("the identity separates the two bases and moves with the component", {
  plan <- pf_plan()
  form <- materialize_population(plan)

  # A complete form and a bank are two estimands, not two renderings of one.
  bank <- estimate_population(plan, rbind(`face-house` = c(1, -1, 0)))
  expect_false(identical(form$scientific_plan_id, bank$scientific_plan_id))
  expect_false(identical(
    form$scientific_plan_id,
    materialize_population(plan, component = "coherent")$scientific_plan_id
  ))
  expect_identical(form$scientific_plan_id,
    materialize_population(plan)$scientific_plan_id)
})

test_that("the query-bank result keeps its own field contract and basis", {
  fit <- estimate_population(pf_plan(), rbind(`face-house` = c(1, -1, 0)))
  expect_identical(fit$basis, "query_bank")
  expect_identical(names(fit), c(
    "basis", "coefficients", "values", "fitted", "residuals", "residual_df",
    "index", "queries", "component", "ledger", "semantics", "normalization",
    "coverage", "uncertainty", "receipt", "scientific_plan_id"
  ))
  expect_identical(fit$receipt$basis, "query_bank")
  expect_identical(fit$receipt$queries, "face-house")
  expect_null(fit$receipt$streaming)
})

test_that("the complete form prints its readout, its bound and what it omits", {
  form <- materialize_population(pf_plan())
  output <- utils::capture.output(print(form))

  expect_true(any(grepl("effect_population_result", output, fixed = TRUE)))
  expect_true(any(grepl("transported_total", output, fixed = TRUE)))
  expect_true(any(grepl("complete 3x3 form", output, fixed = TRUE)))
  expect_true(any(grepl("coordinate tile", output, fixed = TRUE)))
  expect_true(any(grepl("No participant-level arrays", output, fixed = TRUE)))
  expect_true(any(grepl("coefficient_forms", output, fixed = TRUE)))
  expect_match(format(form),
    "^<effect_population_result: transported_total, 4 subjects")
  expect_true(grepl("complete 3x3 form", format(form), fixed = TRUE))

  # A query-bank result still prints its bank under the word `queries`.
  bank <- utils::capture.output(
    print(estimate_population(pf_plan(), rbind(`face-house` = c(1, -1, 0))))
  )
  expect_true(any(grepl("queries:", bank, fixed = TRUE)))
  expect_false(any(grepl("streaming:", bank, fixed = TRUE)))
})

test_that("as.data.frame gives the coefficient forms in long form", {
  form <- materialize_population(pf_plan())
  table <- as.data.frame(form)

  expect_identical(nrow(table), 4L * 1L * 6L)
  expect_true(all(c("node", "sink", "units", "term", "coordinate", "row",
    "column", "scale", "estimate") %in% names(table)))
  expect_identical(table$estimate, as.numeric(form$coefficient_forms))
  # The coordinate metadata travels with the estimate, so a reader can rebuild
  # the symmetric form from the table alone.
  expect_identical(unique(table$scale), c(1, sqrt(2)))
  expect_true(all(table$row >= table$column))
})
