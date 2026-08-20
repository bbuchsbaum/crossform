# The latent PSD descriptive layer (ticket D7).
#
# `design/conservative-geometry-contract.md` section 6 separates the signed
# estimation layer from the latent PSD descriptive layer and forbids the
# second's arithmetic on the first. Section 11.4 gap G7 adds that
# "nonnegativity projection" names no single operator, so the layer has to be
# built by a projection taken from a closed set with the mass it moved on the
# record.
#
# Three of the fixture's four measurements are DIAGONAL, and that is the point
# rather than laziness: the eigenvalues of a diagonal symmetric form are its
# diagonal entries, so the truth for the projected spectrum, the participation
# ratio and the whole cumulative curve is written down in the fixture and
# computed by hand in the assertion, not read back out of the thing under
# test. The fourth is block-diagonal with a 2x2 block whose eigenvalues are
# still exact (`[[2, 1], [1, 2]]` has roots 3 and 1), which is what makes the
# packing claim below real: a diagonal form has no off-diagonal entries, so a
# sqrt(2) scaling error in the codec would be invisible to a diagonal-only
# fixture whichever codec packed it.

latent_fixture <- function() {
  effects <- c("a", "b", "c", "d")
  # Row 1 is already PSD, so the projection must be the identity.
  # Row 2 is indefinite with two known negative roots.
  # Row 3 is negative definite, so nothing survives the projection at all.
  # Row 4 is not diagonal: roots 3, 1, 0.5, -1, one of them negative.
  block <- rbind(
    cbind(matrix(c(2, 1, 1, 2), 2, 2), matrix(0, 2, 2)),
    cbind(matrix(0, 2, 2), diag(c(0.5, -1)))
  )
  total_matrices <- list(
    diag(c(4, 2, 1, 0.5)),
    diag(c(3, 1, -0.5, -1.5)),
    diag(c(-0.25, -0.75, -1, -2)),
    block
  )
  coherent_matrices <- list(
    tcrossprod(c(1, 0.5, 0, 0)),
    tcrossprod(c(0, 1, 0.25, 0)),
    tcrossprod(c(0.5, 0, 0, 1)),
    tcrossprod(c(0.25, 0.25, 1, 0))
  )
  # Packed with the independent oracle codec, so a shared svec/unsvec scaling
  # error cannot cancel between the fixture and the reader.
  pack <- function(values) do.call(rbind, lapply(values, oracle_svec))
  endpoint <- matrix(c(
    1, 2, 3, 4,
    -1, 0, 1, 2,
    0.5, 0.5, 0.5, 0.5,
    2, -1, 0, 1
  ), 4, 4, byrow = TRUE, dimnames = list(NULL, effects))
  marginals <- structure(list(endpoint = endpoint),
    semantics = "undirected_endpoint", class = c("effect_marginals", "list"))
  geometry <- effect_geometry(pack(total_matrices), pack(coherent_matrices),
    marginals, effects, receipt_fixture(),
    index = c("psd", "indefinite", "negative", "block"))
  list(geometry = geometry, total = total_matrices,
    coherent = coherent_matrices)
}

# One measurement per fixture, so each exercises a whole-object branch of the
# print that the mixed fixture above never reaches.
latent_single_fixture <- function(form, id = "m1") {
  effects <- colnames(form)
  if (is.null(effects)) effects <- letters[seq_len(ncol(form))]
  endpoint <- matrix(seq_len(ncol(form)), 1L, ncol(form),
    dimnames = list(NULL, effects))
  effect_geometry(
    matrix(oracle_svec(form), nrow = 1L),
    matrix(oracle_svec(tcrossprod(rep(0, ncol(form)))), nrow = 1L),
    structure(list(endpoint = endpoint), semantics = "undirected_endpoint",
      class = c("effect_marginals", "list")),
    effects, receipt_fixture(), index = id
  )
}

latent_clean_fixture <- function() {
  effects <- c("a", "b")
  pack <- function(values) do.call(rbind, lapply(values, oracle_svec))
  endpoint <- matrix(c(1, 2, 3, 4), 2, 2, byrow = TRUE,
    dimnames = list(NULL, effects))
  effect_geometry(
    pack(list(diag(c(3, 1)), diag(c(2, 2)))),
    pack(list(tcrossprod(c(1, 0)), tcrossprod(c(0, 1)))),
    structure(list(endpoint = endpoint), semantics = "undirected_endpoint",
      class = c("effect_marginals", "list")),
    effects, receipt_fixture(), index = c("m1", "m2")
  )
}

latent_relation_fixture <- function() {
  domain <- abstract_domain(4, id = "latent-relation")
  relation <- relation(
    list(run1 = rbind(a = c(1, 0, 2, 1), b = c(0, 1, 1, 0)),
      run2 = rbind(a = c(1.1, 0.1, 1.9, 0.8), b = c(0.1, 0.9, 1.2, 0.2))),
    domain = domain
  )
  plan_geometry(
    relation, compile_frame(regions(c("v1", "v1", "it", "it")), domain),
    cross_partitions(relation, independence = "independent")
  )
}

test_that("a PSD source passes through the projection unchanged", {
  latent <- latent_geometry(latent_fixture()$geometry)

  # Identity on the already-nonnegative row: same roots, nothing moved.
  expect_equal(unname(latent$spectrum[1, ]), c(4, 2, 1, 0.5),
    tolerance = 1e-12)
  expect_identical(latent$moved_mass[[1L]], 0)
  expect_identical(latent$moved_share[[1L]], 0)

  # The participation ratio, by hand: (4 + 2 + 1 + 0.5)^2 = 56.25 over
  # (16 + 4 + 1 + 0.25) = 21.25.
  expect_equal(latent$n_eff[[1L]], 56.25 / 21.25, tolerance = 1e-12)
  expect_equal(latent$n_eff[[1L]], 2.647058823529412, tolerance = 1e-12)

  # C(k) = cumsum(lambda) / sum(lambda), by hand over a total of 7.5.
  expect_equal(unname(latent$cumulative[1, ]),
    c(4, 6, 7, 7.5) / 7.5, tolerance = 1e-12)
  expect_identical(unname(latent$cumulative[1, 4L]), 1)
  expect_identical(colnames(latent$cumulative), paste0("C", 1:4))
  expect_identical(colnames(latent$spectrum), paste0("root", 1:4))
})

test_that("an indefinite source reports exactly the mass it moved", {
  latent <- latent_geometry(latent_fixture()$geometry)

  # Row 2 is diag(3, 1, -0.5, -1.5): the truncation removes 0.5 + 1.5 = 2 of
  # a total absolute mass of 3 + 1 + 0.5 + 1.5 = 6.
  expect_equal(latent$moved_mass[[2L]], 2, tolerance = 1e-12)
  expect_equal(latent$moved_share[[2L]], 2 / 6, tolerance = 1e-12)
  expect_equal(unname(latent$spectrum[2, ]), c(3, 1, 0, 0), tolerance = 1e-12)
  expect_true(all(latent$spectrum >= 0))

  # n_eff on the PROJECTED spectrum: 4^2 / (9 + 1) = 1.6. On the signed
  # spectrum the same expression is 2^2 / 12.5, which counts nothing.
  expect_equal(latent$n_eff[[2L]], 1.6, tolerance = 1e-12)
  expect_equal(unname(latent$cumulative[2, ]), c(0.75, 1, 1, 1),
    tolerance = 1e-12)

  # The share is against total ABSOLUTE mass, not the signed trace. The
  # signed trace of this row is 2, so a share against it would read exactly
  # 1.0 -- "everything moved" -- for a form that kept two thirds of its mass.
  expect_lt(latent$moved_share[[2L]], 1)
  expect_equal(latent$moved_mass[[2L]] / sum(diag(latent_fixture()$total[[2L]])),
    1, tolerance = 1e-12)
})

test_that("a non-diagonal form is read through the packed codec correctly", {
  # Row 4 is block-diagonal with off-diagonal entries, so it exercises the
  # sqrt(2) scaling in the packed codec that rows 1-3 cannot. Its roots are
  # 3 and 1 from [[2, 1], [1, 2]], then 0.5 and -1.
  latent <- latent_geometry(latent_fixture()$geometry)

  expect_equal(unname(latent$spectrum[4, ]), c(3, 1, 0.5, 0),
    tolerance = 1e-12)
  expect_equal(latent$moved_mass[[4L]], 1, tolerance = 1e-12)
  expect_equal(latent$moved_share[[4L]], 1 / 5.5, tolerance = 1e-12)

  # (3 + 1 + 0.5)^2 = 20.25 over (9 + 1 + 0.25) = 10.25.
  expect_equal(latent$n_eff[[4L]], 20.25 / 10.25, tolerance = 1e-12)
  expect_equal(unname(latent$cumulative[4, ]),
    c(3, 4, 4.5, 4.5) / 4.5, tolerance = 1e-12)
})

test_that("a measurement with no nonnegative mass is masked, not zeroed", {
  latent <- latent_geometry(latent_fixture()$geometry)

  # Row 3 is negative definite: every root is truncated, so there is no
  # nonnegative partition to take a ratio or a share of.
  expect_equal(unname(latent$spectrum[3, ]), c(0, 0, 0, 0), tolerance = 1e-12)
  expect_equal(latent$moved_mass[[3L]], 0.25 + 0.75 + 1 + 2, tolerance = 1e-12)
  expect_equal(latent$moved_share[[3L]], 1, tolerance = 1e-12)
  expect_true(is.na(latent$n_eff[[3L]]))
  expect_true(all(is.na(latent$cumulative[3, ])))
  expect_false(any(is.nan(latent$n_eff)))
  expect_identical(latent$projection$masked, 1L)
})

test_that("an exactly zero form moves no mass and earns no share", {
  # The one place the layer departs from `.coherence_fraction()`'s masking on
  # purpose: a form with no absolute mass at all moved none of it, so the
  # share is 0 rather than NA. Everything else about the row is masked.
  latent <- latent_geometry(latent_single_fixture(matrix(0, 3, 3), "zero"))

  expect_identical(latent$moved_mass, 0)
  expect_identical(latent$moved_share, 0)
  expect_true(is.na(latent$n_eff))
  expect_true(all(is.na(latent$cumulative)))
  expect_identical(latent$projection$masked, 1L)
  expect_identical(latent$projection$clipped, 0L)

  # And the whole-object print branch for a layer that kept nothing anywhere.
  expect_snapshot(print(latent))
})

test_that("the masking guard is exact, and inherits gap G3's open question", {
  # `design/conservative-geometry-contract.md` section 11.4 gap G3 records
  # whether the package's nonnegativity guards should take a RELATIVE
  # tolerance as an open contract decision, shared with every per-node
  # fraction here. D7 inherits it rather than making a private one, so a form
  # that is zero only to within rounding is not masked and its functionals
  # are computed from the rounding. This test characterizes that rather than
  # asserting it is desirable: the roxygen `Masking` section says the same,
  # and the day G3 is decided this test is the one that has to change.
  noise <- latent_geometry(
    latent_single_fixture(diag(c(1e-17, 5e-18, -2e-17)), "noise")
  )
  signal <- latent_geometry(
    latent_single_fixture(diag(c(2, 1, -0.5)), "signal")
  )

  expect_false(is.na(noise$n_eff))
  expect_equal(noise$n_eff, signal$n_eff, tolerance = 1e-9)
  expect_equal(unname(noise$cumulative), unname(signal$cumulative),
    tolerance = 1e-9)
  expect_identical(noise$projection$masked, 0L)

  # What the object does say is the scale: the spectrum is reported, not only
  # the scale-free ratio computed from it.
  expect_lt(max(noise$spectrum), 1e-16)
  expect_gt(max(signal$spectrum), 1)
})

test_that("the projection receipt records moved mass per measurement", {
  fixture <- latent_fixture()
  latent <- latent_geometry(fixture$geometry)
  record <- latent$projection

  expect_s3_class(record, "effect_latent_projection_receipt")
  expect_identical(record$method, "psd_projection")
  expect_identical(record$operator, "eigenvalue_truncation_at_zero")
  expect_identical(record$component, "total")
  expect_identical(record$source, "effect_form")

  # Per measurement, not a summary: four rows in, four moved masses out,
  # aligned to `$index`. That the object and the receipt agree is enforced by
  # the validator and asserted in the forgery test below; what is asserted
  # here is that both equal the hand-computed truth.
  expect_identical(record$measurements, 4L)
  expect_equal(record$moved_mass, c(0, 2, 4, 1), tolerance = 1e-12)
  expect_equal(record$absolute_mass, c(7.5, 6, 4, 5.5), tolerance = 1e-12)
  expect_equal(record$moved_share, c(0, 2 / 6, 1, 1 / 5.5), tolerance = 1e-12)
  expect_identical(record$clipped, 3L)
  expect_identical(record$masked, 1L)
  expect_equal(record$total_moved_mass, 7, tolerance = 1e-12)
  expect_equal(record$max_moved_share, 1, tolerance = 1e-12)

  # The execution receipt names the projection twice over: a derived identity
  # so the latent layer never shares one with the signed source, and the
  # operator on the task partition so a reader does not have to decode a
  # digest to see that a projection happened.
  expect_false(identical(latent$receipt$scientific_plan_id,
    fixture$geometry$receipt$scientific_plan_id))
  expect_identical(latent$receipt$scientific_plan_id,
    record$scientific_plan_id)
  expect_identical(record$source_scientific_plan_id,
    fixture$geometry$receipt$scientific_plan_id)
  expect_match(latent$receipt$task_partition_id, "\\+psd_projection$")
  expect_identical(latent$metadata$scientific_plan_id,
    latent$receipt$scientific_plan_id)
})

test_that("constructing the latent layer leaves the signed source untouched", {
  fixture <- latent_fixture()
  geometry <- fixture$geometry

  before_packed <- serialize(geometry_component(geometry, "total"), NULL)
  before_spectrum <- serialize(geometry_spectrum(geometry)$values, NULL)
  before_signature <- geometry$contract_signature
  before_plan_id <- geometry$receipt$scientific_plan_id
  before_fields <- names(geometry)

  latent <- latent_geometry(geometry)
  expect_s3_class(latent, "effect_latent_geometry")

  # Byte-identical, not merely numerically close: the projection reads the
  # source and returns a new object, it never rewrites a stored geometry.
  expect_identical(serialize(geometry_component(geometry, "total"), NULL),
    before_packed)
  expect_identical(serialize(geometry_spectrum(geometry)$values, NULL),
    before_spectrum)
  expect_identical(geometry$contract_signature, before_signature)
  expect_identical(geometry$receipt$scientific_plan_id, before_plan_id)
  expect_identical(names(geometry), before_fields)

  # And the signed spectrum still is one: the negative roots the latent layer
  # truncated are still there to be read.
  signed <- geometry_spectrum(geometry)
  expect_true(isTRUE(signed$indefinite_estimates_preserved))
  expect_true(any(signed$values < 0))
})

test_that("fractions and cumulative curves live only on the latent layer", {
  # Section 6, normative clause 1: the functionals that need a nonnegative
  # partition are defined here and nowhere upstream. This asserts the
  # signed readouts did not quietly acquire one.
  fixture <- latent_fixture()
  latent <- latent_geometry(fixture$geometry)
  signed <- geometry_spectrum(fixture$geometry)
  energy <- contrast_energy(fixture$geometry, c(a = 1, b = -1, c = 0, d = 0))

  latent_only <- c("cumulative", "n_eff", "moved_mass", "moved_share")
  expect_true(all(latent_only %in% names(latent)))
  expect_identical(intersect(names(signed), latent_only), character(0))
  expect_identical(intersect(names(energy), latent_only), character(0))

  # The one fraction the signed layer does carry is the masked coherence
  # fraction, which is the template rather than the exception: it is still
  # masked, and it is still not a contribution share.
  expect_true(all(c("coherence_fraction", "coherence_fraction_valid") %in%
    names(energy)))
  expect_false("coherence_fraction" %in% names(latent))
})

test_that("the spectrum route carries the same identity as the form route", {
  # Not an independent-agreement check and not written as one: the form route
  # IS `geometry_spectrum()` followed by the spectrum route, so a PSD
  # projection of a symmetric form is the truncation of its eigenvalues by
  # construction and there is no second numerical path. What is asserted is
  # that the plumbing carries the component, the index and the derived
  # identity across the two entry points.
  fixture <- latent_fixture()
  from_form <- latent_geometry(fixture$geometry)
  from_spectrum <- latent_geometry(geometry_spectrum(fixture$geometry))

  expect_identical(from_spectrum$spectrum, from_form$spectrum)
  expect_identical(from_spectrum$n_eff, from_form$n_eff)
  expect_identical(from_spectrum$index, from_form$index)
  expect_identical(from_spectrum$receipt$scientific_plan_id,
    from_form$receipt$scientific_plan_id)
  expect_identical(from_spectrum$projection$source, "effect_spectrum_view")

  # Each named component reaches the projection receipt and gets its own
  # identity, because the three are three different estimands.
  for (component in c("coherent", "configuration")) {
    latent <- latent_geometry(fixture$geometry, component = component)
    expect_identical(latent$component, component)
    expect_identical(latent$projection$component, component)
    expect_false(identical(latent$receipt$scientific_plan_id,
      from_form$receipt$scientific_plan_id))
    expect_true(all(latent$spectrum >= 0))
  }
})

test_that("a spectrum that is not ordered largest root first is refused", {
  # `effect_spectrum_view` is unsealed, and C(k) is a statement about roots
  # taken largest first: read out of order it is still a monotone curve to 1
  # and is no longer a contribution curve.
  spectrum <- geometry_spectrum(latent_fixture()$geometry)
  reversed <- spectrum
  reversed$values <- reversed$values[, rev(seq_len(ncol(reversed$values))),
    drop = FALSE]

  expect_error(latent_geometry(reversed), "largest\\s+root first",
    class = "effect_input_error")
  expect_s3_class(latent_geometry(spectrum), "effect_latent_geometry")
})

test_that("the latent layer refuses everything it cannot project", {
  fixture <- latent_fixture()
  plan <- latent_relation_fixture()

  # A plan names an estimand and holds no geometry to decompose.
  plan_refusal <- catch_refusal(latent_geometry(plan))
  expect_s3_class(plan_refusal, "effect_capability_refusal")
  expect_identical(plan_refusal$capability, "complete_geometry")
  expect_match(conditionMessage(plan_refusal), "A latent geometry requires")
  expect_identical(plan_refusal$remedies,
    "Call `materialize_geometry(x)` and pass the result.")

  # A rectangular form has no symmetric self form to take a spectrum of.
  rectangular <- crossform:::effect_form(
    total = matrix(c(1, 4, -2, 3, 5, 0, 2, -1, 6, 4, 0, 3), nrow = 2,
      byrow = TRUE),
    left_space = effect_space(c("e_a", "e_b"), basis_id = "encoding:v1"),
    right_space = effect_space(c("r_a", "r_b", "r_c"),
      basis_id = "retrieval:v1"),
    receipt = receipt_fixture(), index = c("m1", "m2"),
    codec = "rectangular"
  )
  rectangular_refusal <- catch_refusal(latent_geometry(rectangular))
  expect_identical(rectangular_refusal$capability, "symmetric_self_form")

  # A contracted readout: clamping it would be a per-node total clamp, which
  # gap G7 names as a different operator moving different mass.
  energy <- contrast_energy(fixture$geometry, c(a = 1, b = -1, c = 0, d = 0))
  energy_refusal <- catch_refusal(latent_geometry(energy))
  expect_s3_class(energy_refusal, "effect_capability_refusal")
  expect_identical(energy_refusal$capability, "latent_projection_source")
  expect_identical(energy_refusal$namespace, "latent_geometry")
  expect_identical(energy_refusal$reasons,
    c("readout_is_not_a_symmetric_form",
      "per_node_clamp_is_a_different_projection"))
  expect_match(conditionMessage(energy_refusal), "per-node total clamp")

  # Every readout class the refusal names, not only the first: the roxygen
  # and gap G7 both list them, so each is exercised here.
  readouts <- list(
    rdm = rdm(fixture$geometry),
    rsa = rsa(fixture$geometry,
      models = list(model = abs(outer(1:4, 1:4, "-")))),
    view = query_geometry(fixture$geometry,
      bilinear_query(tcrossprod(c(1, -1, 0, 0))))
  )
  for (name in names(readouts)) {
    expect_identical(
      catch_refusal(latent_geometry(readouts[[name]]))$capability,
      "latent_projection_source", info = name
    )
  }

  # Anything else is a shape error, with the two accepted inputs named.
  expect_error(latent_geometry(diag(3)),
    "complete effect form|effect_spectrum_view", class = "effect_input_error")
})

test_that("the projection set is closed, and a reserved member refuses", {
  fixture <- latent_fixture()

  # A declared-but-unbuilt member refuses rather than falling back: it moves
  # different mass, so substituting `psd_projection` would silently answer a
  # different question (section 11.4, gap G7).
  reserved <- catch_refusal(
    latent_geometry(fixture$geometry, method = "nearest_psd")
  )
  expect_s3_class(reserved, "effect_capability_refusal")
  expect_identical(reserved$capability, "nearest_psd_projection")
  expect_identical(reserved$namespace, "latent_geometry")
  expect_identical(reserved$reasons,
    "projection_kind_declared_but_not_implemented")
  expect_match(conditionMessage(reserved), "not implemented")

  # A name outside the set is a typo, and the message lists the set.
  expect_error(latent_geometry(fixture$geometry, method = "clamp"),
    "psd_projection", class = "effect_input_error")
  expect_error(latent_geometry(fixture$geometry, method = c("a", "b")),
    "psd_projection", class = "effect_input_error")

  # A spectrum has already been decomposed; asking for another component of
  # it is a mistake rather than a second decomposition.
  spectrum <- geometry_spectrum(fixture$geometry, component = "coherent")
  expect_error(latent_geometry(spectrum, component = "total"),
    "fixed by the spectrum", class = "effect_input_error")
  expect_identical(
    latent_geometry(spectrum, component = "coherent")$component, "coherent"
  )
})

test_that("the latent layer is sealed against forgery", {
  latent <- latent_geometry(latent_fixture()$geometry)

  forged <- latent
  forged$spectrum[2, 3] <- -1
  expect_error(crossform:::.validate_effect_latent_geometry(forged),
    "nonnegative", class = "effect_contract_error")

  # `n_eff` and `C(k)` are functions of the spectrum, so a plausible value
  # beside a spectrum that does not produce it is caught. Range checks alone
  # would pass both of these.
  forged <- latent
  forged$n_eff[[1L]] <- 3
  expect_error(crossform:::.validate_effect_latent_geometry(forged),
    "not the one its own", class = "effect_contract_error")

  forged <- latent
  forged$cumulative[1, ] <- c(0.1, 0.2, 0.3, 1)
  expect_error(crossform:::.validate_effect_latent_geometry(forged),
    "not the one its own", class = "effect_contract_error")

  forged <- latent
  forged$cumulative[1, ] <- 2
  expect_error(crossform:::.validate_effect_latent_geometry(forged),
    "leaves \\[0, 1\\]", class = "effect_contract_error")

  # Zeroing the moved mass on BOTH the object and the receipt keeps them
  # consistent with each other, and is still caught: what was kept plus what
  # was moved has to add back to the source's absolute mass.
  forged <- latent
  forged$moved_mass <- rep(0, 4)
  forged$projection$moved_mass <- rep(0, 4)
  forged$projection$clipped <- 0L
  forged$projection$total_moved_mass <- 0
  expect_error(crossform:::.validate_effect_latent_geometry(forged),
    "do not add back", class = "effect_contract_error")

  # The receipt's summaries are the fields the print reads, so they are
  # re-derived rather than trusted: an honest per-measurement series under a
  # "nothing was clipped" headline is exactly the silent clipping section 6
  # forbids.
  forged <- latent
  forged$projection$clipped <- 0L
  expect_error(crossform:::.validate_effect_latent_geometry(forged),
    "summarizes moved mass it does not carry", class = "effect_contract_error")

  forged <- latent
  forged$projection$total_moved_mass <- 0
  expect_error(crossform:::.validate_effect_latent_geometry(forged),
    "summarizes moved mass it does not carry", class = "effect_contract_error")

  forged <- latent
  forged$projection$masked <- 0L
  expect_error(crossform:::.validate_effect_latent_geometry(forged),
    "miscounts the measurements", class = "effect_contract_error")

  forged <- latent
  forged$moved_mass <- rep(0, 4)
  expect_error(crossform:::.validate_effect_latent_geometry(forged),
    "do not add back", class = "effect_contract_error")

  forged <- latent
  forged$component <- "coherent"
  expect_error(crossform:::.validate_effect_latent_geometry(forged),
    "disagree about what was projected", class = "effect_contract_error")

  forged <- latent
  forged$extra <- TRUE
  expect_error(crossform:::.validate_effect_latent_geometry(forged),
    "canonical", class = "effect_input_error")
})

test_that("the latent layer prints what it is and what the clipping cost", {
  fixture <- latent_fixture()
  latent <- latent_geometry(fixture$geometry)

  # The acceptance sentence is asserted directly as well as by snapshot, so a
  # snapshot accepted in a hurry cannot drop it.
  printed <- capture.output(print(latent))
  expect_true(any(grepl("latent descriptive layer; not for inference", printed,
    fixed = TRUE)))
  expect_true(any(grepl("psd_projection", printed, fixed = TRUE)))
  expect_true(any(grepl("moved_mass", printed, fixed = TRUE)))
  expect_true(any(grepl("n_eff", printed, fixed = TRUE)))

  expect_snapshot(print(latent))
  expect_snapshot(format(latent))
  expect_snapshot(print(latent$projection))
  expect_snapshot(format(latent$projection))

  # A source that needed no clipping says so rather than staying quiet.
  clean <- latent_geometry(latent_clean_fixture())
  clean_printed <- capture.output(print(clean))
  expect_true(any(grepl("none", clean_printed, fixed = TRUE)))
  expect_true(any(grepl("no source root was negative", clean_printed,
    fixed = TRUE)))
  expect_snapshot(print(clean))
})

test_that("the latent layer coerces to one row per measurement", {
  latent <- latent_geometry(latent_fixture()$geometry)
  frame <- as.data.frame(latent)

  expect_identical(nrow(frame), 4L)
  expect_identical(frame$measurement,
    c("psd", "indefinite", "negative", "block"))
  expect_identical(names(frame), c("measurement", "n_eff", "moved_mass",
    "moved_share", paste0("root", 1:4)))
  expect_equal(frame$moved_mass, c(0, 2, 4, 1), tolerance = 1e-12)
  expect_true(is.na(frame$n_eff[[3L]]))
})
