# A location transport is the population layer's one estimand-bearing input.
# These tests pin the laws `design/population-form-contract.md` §1, §2 and
# §11 make normative -- row-stochasticity including a required sink, exact
# budget preservation for SIGNED ledgers, density as a declared row-mass ratio
# with the sink still in budget units, and cross-fit provenance on any
# data-derived operator -- at the tolerances §11 states.

transport_block <- function() {
  # Node 3 splits 0.6/0.4 across the two group nodes, node 4 leaves 30% of its
  # territory unmapped, node 5 is unmapped entirely. The same shape as the
  # oracle fixture in `design/oracles/population-transport-contract.R` §P1.a.
  rbind(
    c(1, 0),
    c(1, 0),
    c(0.6, 0.4),
    c(0, 0.7),
    c(0, 0)
  )
}

transport_fixture <- function(semantics = "budget", ...) {
  location_transport(
    transport_block(),
    native_index = c("v1", "v2", "v3", "v4", "v5"),
    group_index = c("left", "right"),
    semantics = semantics,
    provenance = list(method = "anatomical", details = "fixed warp, atlas v1"),
    ...
  )
}

test_that("the constructor materializes the sink and closes every row", {
  transport <- transport_fixture()

  expect_s3_class(transport, "effect_location_transport")
  expect_identical(names(transport), c("matrix", "native_index",
    "group_index", "semantics", "row_mass", "provenance", "signature"))
  expect_true(methods::is(transport$matrix, "sparseMatrix"))

  # The sink is column m + 1 and is not indexed: `group_index` describes the
  # group nodes only.
  expect_identical(dim(transport$matrix), c(5L, 3L))
  expect_identical(nrow(transport$group_index), 2L)
  expect_identical(nrow(transport$native_index), 5L)

  operator <- as.matrix(transport$matrix)
  expect_equal(as.numeric(rowSums(operator)), rep(1, 5), tolerance = 1e-12)
  expect_equal(as.numeric(operator[, 3L]), c(0, 0, 0, 0.3, 1),
    tolerance = 1e-12)
  expect_true(all(operator >= 0))

  # The default declared row mass is the unit vector (§1.3, §14.6).
  expect_identical(transport$row_mass, rep(1, 5))
  expect_identical(transport$semantics, "budget")
  expect_true(.strong_sha256(transport$signature))
})

test_that("a sink column exists even when nothing is unmapped", {
  full <- location_transport(
    rbind(c(1, 0), c(0, 1)), c("a", "b"), c("g1", "g2"),
    semantics = "budget",
    provenance = list(method = "anatomical", details = "identity assignment")
  )
  expect_identical(dim(full$matrix), c(2L, 3L))
  expect_identical(sum(full$matrix[, 3L]), 0)
})

test_that("the constructor refuses an ill-formed operator", {
  # (a) a row claiming more than unit mass would need negative sink mass.
  expect_error(
    location_transport(rbind(c(0.8, 0.4)), "v1", c("left", "right"),
      semantics = "budget",
      provenance = list(method = "anatomical", details = "warp")),
    "sink mass would be negative",
    class = "effect_input_error"
  )

  # (b) a negative entry.
  expect_error(
    location_transport(rbind(c(1.2, -0.2)), "v1", c("left", "right"),
      semantics = "budget",
      provenance = list(method = "anatomical", details = "warp")),
    "must be nonnegative",
    class = "effect_input_error"
  )

  # (c) a non-finite weight declares nothing.
  expect_error(
    location_transport(rbind(c(NA_real_, 0.5)), "v1", c("left", "right"),
      semantics = "budget",
      provenance = list(method = "anatomical", details = "warp")),
    "no missing or infinite entries",
    class = "effect_input_error"
  )

  # (d) tolerance is a refusal threshold, not a renormalizing licence.
  expect_error(
    location_transport(rbind(c(1 + 1e-6, 0)), "v1", c("left", "right"),
      semantics = "budget", tolerance = 1e-9,
      provenance = list(method = "anatomical", details = "warp")),
    class = "effect_input_error"
  )
  expect_s3_class(
    location_transport(rbind(c(1 + 1e-12, 0)), "v1", c("left", "right"),
      semantics = "budget", tolerance = 1e-9,
      provenance = list(method = "anatomical", details = "warp")),
    "effect_location_transport"
  )
})

test_that("semantics is required and closed", {
  expect_error(
    location_transport(transport_block(), c("v1", "v2", "v3", "v4", "v5"),
      c("left", "right"),
      provenance = list(method = "anatomical", details = "warp")),
    "must be declared",
    class = "effect_input_error"
  )
  expect_error(
    transport_fixture(semantics = "mean"),
    class = "effect_input_error"
  )
  expect_error(
    transport_fixture(semantics = c("budget", "density")),
    class = "effect_input_error"
  )
})

test_that("a declared row mass must be positive, and is carried", {
  declared <- transport_fixture(semantics = "density",
    row_mass = c(1, 2, 3, 4, 5))
  expect_identical(declared$row_mass, c(1, 2, 3, 4, 5))

  expect_error(transport_fixture(row_mass = c(1, 0, 1, 1, 1)),
    class = "effect_input_error")
  expect_error(transport_fixture(row_mass = c(1, -1, 1, 1, 1)),
    class = "effect_input_error")
  expect_error(transport_fixture(row_mass = c(1, 1, 1)),
    class = "effect_input_error")
})

test_that("index tables must match the operator and name their rows", {
  block <- transport_block()
  provenance <- list(method = "anatomical", details = "warp")

  expect_error(
    location_transport(block, c("v1", "v2"), c("left", "right"),
      semantics = "budget", provenance = provenance),
    "one row per node",
    class = "effect_input_error"
  )
  expect_error(
    location_transport(block, c("v1", "v2", "v3", "v4", "v5"), "left",
      semantics = "budget", provenance = provenance),
    "one row per node",
    class = "effect_input_error"
  )
  expect_error(
    location_transport(block, c("v1", "v1", "v3", "v4", "v5"),
      c("left", "right"), semantics = "budget", provenance = provenance),
    "nonempty and unique",
    class = "effect_input_error"
  )
  expect_error(
    location_transport(block, data.frame(label = paste0("v", 1:5)),
      c("left", "right"), semantics = "budget", provenance = provenance),
    "must carry a `node` column",
    class = "effect_input_error"
  )
  expect_error(
    location_transport(block, c("v1", "v2", "v3", "v4", "v5"),
      data.frame(node = c("left", "right"), scale = c("a", "b")),
      semantics = "budget", provenance = provenance),
    "must be numeric",
    class = "effect_input_error"
  )
})

test_that("per-row frame metadata survives construction unchanged", {
  index <- data.frame(
    node = paste0("v", 1:5),
    family = c("fine", "fine", "fine", "coarse", "coarse"),
    scale = c(1, 1, 1, 3, 3),
    center = c("f1", "f2", "f3", "f9", "f12"),
    alpha = c(0.5, 0.5, 0.5, 0.5, 0.5),
    stringsAsFactors = FALSE
  )
  transport <- location_transport(transport_block(), index,
    c("left", "right"), semantics = "budget",
    provenance = list(method = "anatomical", details = "warp"))

  expect_identical(transport$native_index, index)
  expect_identical(transport$native_index$family, index$family)
  expect_identical(transport$native_index$scale, index$scale)
  expect_identical(transport$native_index$alpha, index$alpha)
})

test_that("a functional transport without cross-fit provenance is refused", {
  refusal <- catch_refusal(location_transport(
    transport_block(), c("v1", "v2", "v3", "v4", "v5"), c("left", "right"),
    semantics = "budget",
    provenance = list(method = "functional", details = "response clustering")
  ))
  expect_s3_class(refusal, "effect_capability_refusal")
  expect_identical(refusal$capability, "cross_fit_provenance")
  expect_identical(refusal$namespace, "location_transport")
  expect_identical(refusal$reasons, "cross_fit_partitions_not_declared")
  expect_true(nzchar(refusal$remedies))

  # Naming the partitions that built it is what the refusal asks for.
  honest <- location_transport(
    transport_block(), c("v1", "v2", "v3", "v4", "v5"), c("left", "right"),
    semantics = "budget",
    provenance = list(method = "functional", details = "response clustering",
      cross_fit = c("run-1", "run-2"))
  )
  expect_identical(honest$provenance$cross_fit, c("run-1", "run-2"))

  # An anatomical or external transport never saw the responses, so it owes
  # no cross-fit record.
  expect_s3_class(transport_fixture(), "effect_location_transport")
})

test_that("provenance must declare a closed method and how it was built", {
  block <- transport_block()
  nodes <- paste0("v", 1:5)
  expect_error(
    location_transport(block, nodes, c("left", "right"), semantics = "budget",
      provenance = list(method = "atlas", details = "warp")),
    class = "effect_input_error"
  )
  expect_error(
    location_transport(block, nodes, c("left", "right"), semantics = "budget",
      provenance = list(method = "anatomical")),
    "must be one nonempty string",
    class = "effect_input_error"
  )
  expect_error(
    location_transport(block, nodes, c("left", "right"), semantics = "budget",
      provenance = "anatomical"),
    class = "effect_input_error"
  )
  expect_error(
    location_transport(block, nodes, c("left", "right"), semantics = "budget"),
    "must declare how the transport was built",
    class = "effect_input_error"
  )

  # `method` and `details` lead the record, so a provenance typed in another
  # order is the same identity.
  a <- location_transport(block, nodes, c("left", "right"),
    semantics = "budget",
    provenance = list(method = "anatomical", details = "warp", atlas = "x"))
  b <- location_transport(block, nodes, c("left", "right"),
    semantics = "budget",
    provenance = list(details = "warp", atlas = "x", method = "anatomical"))
  expect_identical(a$signature, b$signature)
  expect_identical(names(a$provenance), c(
    "method", "details", "fitting_sample", "cross_fit_folds",
    "conditioning", "atlas"
  ))
})

test_that("budget preservation is exact, for signed ledgers, incl. the sink", {
  transport <- transport_fixture()
  ledger <- c(1.5, -0.5, 2, -1, 0.25)
  carried <- transport_values(transport, ledger)

  expect_identical(names(carried), c("left", "right", "<sink>"))
  # `population-form-v1` §11: tolerance scaled by the ledger's L1 norm.
  scale <- 1e-12 * sum(abs(ledger))
  expect_lt(abs(sum(carried) - sum(ledger)), scale)
  expect_lt(
    abs(sum(carried[1:2]) - (sum(ledger) - carried[["<sink>"]])), scale
  )

  # Hand-checkable: left = v1 + v2 + 0.6*v3, right = 0.4*v3 + 0.7*v4.
  expect_equal(carried[["left"]], 1.5 - 0.5 + 0.6 * 2, tolerance = 1e-12)
  expect_equal(carried[["right"]], 0.4 * 2 + 0.7 * -1, tolerance = 1e-12)
  expect_equal(carried[["<sink>"]], 0.3 * -1 + 1 * 0.25, tolerance = 1e-12)

  # Deleting the sink loses signed mass, which is the whole reason it is
  # required rather than optional.
  expect_gt(abs(sum(carried[1:2]) - sum(ledger)), 1e-3)
})

test_that("budget preservation holds over random transports", {
  set.seed(20260820)
  worst <- 0
  for (draw in seq_len(200L)) {
    n <- sample(4:16, 1L)
    m <- sample(2:6, 1L)
    block <- matrix(stats::rexp(n * m), n, m)
    block[stats::runif(length(block)) < 0.6] <- 0
    keep <- stats::runif(n, 0.2, 1)
    row_sums <- rowSums(block)
    block[row_sums > 0, ] <- block[row_sums > 0, , drop = FALSE] /
      row_sums[row_sums > 0] * keep[row_sums > 0]
    transport <- location_transport(
      block, as.character(seq_len(n)), paste0("g", seq_len(m)),
      semantics = "budget",
      provenance = list(method = "anatomical", details = "random draw")
    )
    ledger <- stats::rnorm(n, sd = 3)
    carried <- transport_values(transport, ledger)
    worst <- max(worst, abs(sum(carried) - sum(ledger)) / sum(abs(ledger)))
  }
  expect_lt(worst, 1e-12)
})

test_that("a matrix of values transports column by column", {
  transport <- transport_fixture()
  values <- cbind(a = c(1.5, -0.5, 2, -1, 0.25), b = c(1, 1, 1, 1, 1))
  carried <- transport_values(transport, values)

  expect_identical(dim(carried), c(3L, 2L))
  expect_identical(rownames(carried), c("left", "right", "<sink>"))
  expect_identical(colnames(carried), c("a", "b"))
  expect_equal(carried[, "a"], transport_values(transport, values[, "a"]),
    tolerance = 1e-12)
  expect_equal(colSums(carried), colSums(values), tolerance = 1e-12)
})

test_that("transport_values refuses values it cannot carry", {
  transport <- transport_fixture()
  expect_error(transport_values(transport, c(1, 2, 3)),
    class = "effect_input_error")
  expect_error(transport_values(transport, c(1, 2, 3, 4, NA)),
    "must be finite", class = "effect_input_error")
  expect_error(transport_values(transport, letters[1:5]),
    class = "effect_input_error")
  expect_error(transport_values(transport, matrix(1, 3, 2)),
    class = "effect_input_error")
})

test_that("density is the declared row-mass ratio, with the sink in budget", {
  row_mass <- c(1, 2, 3, 4, 5)
  transport <- transport_fixture(semantics = "density", row_mass = row_mass)
  ledger <- c(1.5, -0.5, 2, -1, 0.25)

  operator <- as.matrix(transport$matrix)
  numerator <- as.numeric(crossprod(operator, ledger))
  denominator <- as.numeric(crossprod(operator, row_mass))
  carried <- transport_values(transport, ledger)

  # `population-form-v1` §11: density equals the fixed linear map
  # D(1/(P^T mu)) P^T applied to the data, to 1e-12.
  linear <- diag(1 / denominator[1:2]) %*% t(operator)[1:2, , drop = FALSE]
  expect_equal(as.numeric(carried[1:2]), as.numeric(linear %*% ledger),
    tolerance = 1e-12)

  # The sink is an accounting column and stays in budget units under either
  # semantics.
  expect_equal(carried[["<sink>"]], numerator[[3L]], tolerance = 1e-12)
  expect_equal(carried[["<sink>"]],
    transport_values(transport_fixture(row_mass = row_mass),
      ledger)[["<sink>"]],
    tolerance = 1e-12)

  # Density satisfies no conservation law; that is what it trades away.
  expect_gt(abs(sum(carried) - sum(ledger)), 1e-6)
})

test_that("a group node reached by no native mass is NA, never zero", {
  transport <- location_transport(
    rbind(c(1, 0, 0), c(1, 0, 0), c(0, 0.5, 0)),
    c("v1", "v2", "v3"), c("reached", "partial", "empty"),
    semantics = "density",
    provenance = list(method = "anatomical", details = "warp")
  )
  carried <- transport_values(transport, c(2, 4, 6))
  expect_equal(carried[["reached"]], 3, tolerance = 1e-12)
  expect_equal(carried[["partial"]], 6, tolerance = 1e-12)
  expect_true(is.na(carried[["empty"]]))
  expect_false(isTRUE(carried[["empty"]] == 0))

  # The sink still reports the unmapped budget as a number.
  expect_equal(carried[["<sink>"]], 0.5 * 6, tolerance = 1e-12)
})

test_that("the signature is content-addressed over the whole estimand", {
  base <- transport_fixture()
  expect_identical(base$signature, transport_fixture()$signature)

  # Semantics changes the estimand, so it changes the identity.
  expect_false(identical(base$signature,
    transport_fixture(semantics = "density")$signature))

  # So does one entry of the operator itself.
  moved <- transport_block()
  moved[3L, ] <- c(0.5, 0.5)
  shifted <- location_transport(moved, paste0("v", 1:5), c("left", "right"),
    semantics = "budget",
    provenance = list(method = "anatomical", details = "fixed warp, atlas v1"))
  expect_false(identical(base$signature, shifted$signature))

  # And the declared row mass, the labels, and the provenance.
  expect_false(identical(base$signature,
    transport_fixture(row_mass = c(1, 1, 1, 1, 2))$signature))
  expect_false(identical(base$signature, location_transport(
    transport_block(), paste0("w", 1:5), c("left", "right"),
    semantics = "budget",
    provenance = list(method = "anatomical", details = "fixed warp, atlas v1")
  )$signature))
  expect_false(identical(base$signature, location_transport(
    transport_block(), paste0("v", 1:5), c("left", "right"),
    semantics = "budget",
    provenance = list(method = "anatomical", details = "a different warp")
  )$signature))
})

test_that("the validator catches a tampered transport", {
  transport <- transport_fixture()
  expect_identical(crossform:::.validate_location_transport(transport),
    transport)

  relabelled <- transport
  relabelled$semantics <- "density"
  expect_error(crossform:::.validate_location_transport(relabelled),
    class = "effect_contract_error")
  expect_error(transport_values(relabelled, rep(1, 5)),
    class = "effect_contract_error")

  truncated <- transport
  truncated$native_index <- truncated$native_index[1:3, , drop = FALSE]
  expect_error(crossform:::.validate_location_transport(truncated),
    "do not match the operator", class = "effect_contract_error")

  stripped <- unclass(transport)
  expect_error(
    crossform:::.validate_location_transport(structure(stripped[1:3],
      class = "effect_location_transport")),
    class = "effect_input_error"
  )
})

test_that("anatomical_transport assigns to the nearest centre", {
  # Two group centres 10 apart on a line, six native nodes placed so that
  # every branch of the rule is exercised by hand.
  native <- rbind(c(0, 0), c(1, 0), c(9, 0), c(10, 0), c(5, 0), c(0, 20))
  group <- rbind(c(0, 0), c(10, 0))

  bounded <- anatomical_transport(native, group, semantics = "budget",
    radius = 2)
  operator <- as.matrix(bounded$matrix)
  expect_identical(dim(operator), c(6L, 3L))
  # v1, v2 -> g1 (distance 0, 1); v3, v4 -> g2 (distance 1, 0);
  # v5 is equidistant at 5 and v6 is 20 away, so both exceed the radius.
  expect_identical(as.numeric(operator[, 1L]), c(1, 1, 0, 0, 0, 0))
  expect_identical(as.numeric(operator[, 2L]), c(0, 0, 1, 1, 0, 0))
  expect_identical(as.numeric(operator[, 3L]), c(0, 0, 0, 0, 1, 1))
  expect_equal(as.numeric(rowSums(operator)), rep(1, 6), tolerance = 1e-12)

  # Without a radius nothing is unmapped, the sink is materialized at zero,
  # and the equidistant node goes to the lowest group position.
  unbounded <- anatomical_transport(native, group, semantics = "budget")
  full <- as.matrix(unbounded$matrix)
  expect_identical(as.numeric(full[, 1L]), c(1, 1, 0, 0, 1, 1))
  expect_identical(sum(full[, 3L]), 0)

  # Every row is a hard assignment, so the operator is an indicator matrix.
  expect_true(all(full %in% c(0, 1)))

  # The centres are carried so displacement can be recomputed from the
  # transport alone.
  expect_identical(bounded$native_index$coord1, native[, 1L])
  expect_identical(bounded$native_index$coord2, native[, 2L])
  expect_identical(bounded$group_index$coord1, group[, 1L])
  expect_identical(bounded$provenance$method, "anatomical")
  expect_identical(bounded$provenance$rule, "nearest_centre")
  expect_identical(bounded$provenance$radius, 2)
})

test_that("anatomical_transport refuses inconsistent geometry", {
  native <- rbind(c(0, 0), c(1, 0))
  expect_error(
    anatomical_transport(native, rbind(c(0, 0, 0)), semantics = "budget"),
    "same coordinate system", class = "effect_input_error"
  )
  expect_error(anatomical_transport(native, rbind(c(0, 0))),
    "must be declared", class = "effect_input_error")
  expect_error(
    anatomical_transport(native, rbind(c(0, 0)), semantics = "budget",
      radius = -1),
    class = "effect_input_error"
  )
  expect_error(
    anatomical_transport(native, rbind(c(0, 0)), semantics = "budget",
      provenance = list(method = "functional", details = "clustering")),
    "cannot be declared", class = "effect_input_error"
  )
  expect_error(
    anatomical_transport(native, rbind(c(0, 0)), semantics = "budget",
      native_index = data.frame(node = c("a", "b"), coord1 = c(0, 1))),
    "already carries", class = "effect_input_error"
  )
})

test_that("external_transport admits a declared operator and its names", {
  P <- matrix(c(0.7, 0.2, 0, 0.3, 0.8, 0.5), nrow = 3,
    dimnames = list(c("n1", "n2", "n3"), c("A", "B")))
  transport <- external_transport(P, semantics = "budget",
    provenance = list(details = "partial-volume warp from atlas-tool 2.1"))

  expect_identical(transport$provenance$method, "external")
  expect_identical(transport$native_index$node, c("n1", "n2", "n3"))
  expect_identical(transport$group_index$node, c("A", "B"))
  expect_equal(as.numeric(transport$matrix[, 3L]), c(0, 0, 0.5),
    tolerance = 1e-12)

  # Without dimension names the nodes are named by position.
  bare <- external_transport(unname(P), semantics = "budget",
    provenance = list(details = "same operator, unnamed"))
  expect_identical(bare$native_index$node, c("native1", "native2", "native3"))
  expect_identical(bare$group_index$node, c("group1", "group2"))

  # "external" is not a route around the cross-fit obligation: an operator
  # fitted to responses is functional whoever fitted it.
  refusal <- catch_refusal(external_transport(P, semantics = "budget",
    provenance = list(method = "functional", details = "hyperalignment")))
  expect_s3_class(refusal, "effect_capability_refusal")
  expect_identical(refusal$capability, "cross_fit_provenance")

  expect_error(external_transport(P, semantics = "budget"),
    "must say how the external transport was built",
    class = "effect_input_error")
  expect_error(external_transport(P, provenance = list(details = "x")),
    "must be declared", class = "effect_input_error")
})

test_that("a transport prints its shape, semantics, sink and provenance", {
  transport <- transport_fixture()
  expect_snapshot(
    print(transport),
    transform = function(lines) {
      sub("(signature:\\s+sha256:)[0-9a-f]+", "\\1<digest>", lines)
    }
  )
  expect_snapshot(format(transport))

  # Density names its denominator, because the ratio is per unit of it.
  expect_snapshot(
    print(transport_fixture(semantics = "density",
      row_mass = c(1, 2, 3, 4, 5))),
    transform = function(lines) {
      sub("(signature:\\s+sha256:)[0-9a-f]+", "\\1<digest>", lines)
    }
  )
})
