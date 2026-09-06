# The model basis (ticket M2 of the model-coordinate epic).
#
# `design/model-coordinate-geometry-contract.md` sections 1, 1.5 and 6 are
# normative here, and `design/oracles/model-coordinate-geometry.R` section O2
# is the measurement behind the first test: a thin QR of a rank-deficient
# model factor is not centered, and the compressed form built on it moves
# under a condition-independent baseline shift. Every law below is checked
# against base matrix algebra written in the test, never against the value
# under test.

basis_conditions <- c("face", "body", "house", "tool", "chair", "scene")

category_rdm <- function(labels, conditions = basis_conditions) {
  D <- outer(labels, labels, function(a, b) as.numeric(a != b))
  dimnames(D) <- list(conditions, conditions)
  D
}

feature_matrix <- function(seed, columns, conditions = basis_conditions) {
  set.seed(seed)
  F <- matrix(rnorm(length(conditions) * columns), length(conditions), columns)
  dimnames(F) <- list(conditions, paste0("f", seq_len(columns)))
  F
}

squared_euclidean <- function(F) {
  D <- as.matrix(stats::dist(F))^2
  dimnames(D) <- list(rownames(F), rownames(F))
  D
}

oracle_gram <- function(D) {
  n <- nrow(D)
  H <- diag(n) - 1 / n
  -0.5 * H %*% D %*% H
}

oracle_factor <- function(D, rank, tolerance = 1e-10) {
  e <- eigen(oracle_gram(D), symmetric = TRUE)
  keep <- which(e$values > tolerance * max(abs(e$values)))
  keep <- keep[seq_len(min(rank, length(keep)))]
  e$vectors[, keep, drop = FALSE] %*%
    diag(sqrt(e$values[keep]), length(keep))
}

projector <- function(X) {
  X %*% solve(crossprod(X), t(X))
}

test_that("a rank request above a model's rank is clamped and the basis stays centered", {
  # The trap the contract's rule 1.1(a) exists for: the literal two-column
  # factor of a rank-1 category Gram has a zero column, and a thin QR of it
  # returns a second vector that is not centered.
  category <- category_rdm(c(1, 1, 1, 0, 0, 0))
  e <- eigen(oracle_gram(category), symmetric = TRUE)
  literal <- e$vectors[, 1:2] %*% diag(sqrt(pmax(e$values[1:2], 0)), 2L)
  Q_trap <- qr.Q(qr(literal))
  expect_gt(max(abs(colSums(Q_trap))), 0.5)

  basis <- model_basis(list(category = category),
    distance = "squared_euclidean", rank = 2)

  expect_s3_class(basis, "effect_model_basis")
  expect_identical(basis$models$rank_requested, 2L)
  expect_identical(basis$models$rank_effective, 1L)
  expect_true(basis$models$rank_clamped)
  expect_identical(basis$dimension, 1L)
  expect_true(basis$baseline_invariant)
  expect_lt(basis$centering_error, 1e-12)
  expect_lt(max(abs(colSums(basis$Q))), 1e-12)
  # The single direction is the centered category indicator.
  indicator <- c(1, 1, 1, 0, 0, 0) - 0.5
  expect_equal(unname(abs(drop(crossprod(basis$Q,
    indicator / sqrt(sum(indicator^2)))))), 1, tolerance = 1e-12)
  expect_equal(unname(basis$Q %*% basis$R),
    oracle_factor(category, 1L) * sign(sum(basis$Q * oracle_factor(category, 1L))),
    tolerance = 1e-12)
  expect_identical(basis$span_overlap_rank, 0L)
  expect_equal(basis$span_fraction, 1 / 5)
  expect_false(basis$saturated)
})

test_that("the basis factors the stacked model factor and spans the joint model span", {
  category <- category_rdm(c(1, 1, 1, 0, 0, 0))
  F <- feature_matrix(1, 3L)
  D <- squared_euclidean(F)
  basis <- model_basis(list(category = category, feature = D),
    distance = "squared_euclidean", rank = c(category = 2, feature = 2))

  T_oracle <- cbind(oracle_factor(category, 1L), oracle_factor(D, 2L))
  expect_identical(basis$models$rank_effective, c(1L, 2L))
  expect_identical(basis$models$rank_clamped, c(TRUE, FALSE))
  expect_identical(basis$dimension, 3L)
  expect_identical(dim(basis$R), c(3L, 3L))
  expect_identical(basis$columns, list(category = 1L, feature = 2:3))
  expect_identical(colnames(basis$R), c("category.1", "feature.1", "feature.2"))
  expect_equal(unname(crossprod(basis$Q)), diag(3), tolerance = 1e-12)
  # Each model's factor is recovered up to the sign convention on its
  # columns: the projectors onto the per-model column spans agree exactly.
  T_pkg <- basis$Q %*% basis$R
  expect_equal(unname(projector(T_pkg[, 1, drop = FALSE])),
    projector(T_oracle[, 1, drop = FALSE]), tolerance = 1e-12)
  expect_equal(unname(projector(T_pkg[, 2:3])), projector(T_oracle[, 2:3]),
    tolerance = 1e-12)
  expect_equal(unname(projector(T_pkg)), projector(T_oracle), tolerance = 1e-12)
  expect_equal(abs(unname(T_pkg)), abs(T_oracle), tolerance = 1e-12)
  s <- svd(T_oracle)$d
  expect_equal(basis$basis_condition, s[[1L]] / s[[3L]], tolerance = 1e-10)
  expect_equal(basis$span_fraction, 3 / 5)
  expect_identical(names(basis$spectra), c("category", "feature"))
  expect_equal(basis$spectra$feature,
    eigen(oracle_gram(D), symmetric = TRUE)$values, tolerance = 1e-12)
})

test_that("features and their squared Euclidean RDM declare the same span", {
  F <- feature_matrix(2, 2L)
  from_features <- model_basis(features = list(shape = F))
  from_rdm <- model_basis(list(shape = squared_euclidean(F)),
    distance = "squared_euclidean")

  expect_identical(from_features$models$kind, "features")
  expect_identical(from_rdm$models$kind, "rdm")
  expect_identical(from_features$distance, NA_character_)
  expect_identical(from_rdm$distance, "squared_euclidean")
  expect_identical(from_features$dimension, 2L)
  expect_equal(projector(from_features$Q), projector(from_rdm$Q),
    tolerance = 1e-12)
  # Same span, different declared kind: different objects.
  expect_false(identical(from_features$signature, from_rdm$signature))
})

test_that("the dissimilarity is declared explicitly from a closed set", {
  category <- category_rdm(c(1, 1, 1, 0, 0, 0))
  expect_error(model_basis(list(category = category)),
    "explicit", class = "effect_input_error")
  expect_error(model_basis(list(category = category), distance = "correlation"),
    "squared_euclidean", class = "effect_input_error")
  expect_error(model_basis(features = list(a = feature_matrix(3, 2L)),
    distance = "squared_euclidean"), "supply `models` too",
    class = "effect_input_error")
  expect_error(model_basis(), "at least one model", class = "effect_input_error")
  expect_error(model_basis(list(category), distance = "squared_euclidean"),
    "unique nonempty names", class = "effect_input_error")
  expect_error(model_basis(list(a = category), distance = "squared_euclidean",
    features = list(a = feature_matrix(3, 2L))), "unique across",
    class = "effect_input_error")
})

test_that("a non-Euclidean model is refused, and dropped negative mass is on the record", {
  F <- feature_matrix(4, 3L)
  D <- squared_euclidean(F)
  # Fourth powers of Euclidean distances are not squared Euclidean distances;
  # the centered Gram is indefinite with substantial negative mass.
  bent <- D^2
  values <- eigen(oracle_gram(bent), symmetric = TRUE)$values
  negative <- -sum(pmin(values, 0))
  positive <- sum(values[values > 1e-10 * max(abs(values))])
  expect_gt(negative / positive, 0.01)

  refusal <- catch_refusal(model_basis(list(bent = bent),
    distance = "squared_euclidean"))
  expect_s3_class(refusal, "effect_capability_refusal")
  expect_identical(refusal$capability, "euclidean_model_geometry")
  expect_identical(refusal$namespace, "model_coordinate")
  expect_identical(refusal$reasons[[1L]], "model_gram_indefinite")
  expect_true(any(grepl("negative_share", refusal$remedies)))
  expect_true(any(grepl("features", refusal$remedies)))

  admitted <- model_basis(list(bent = bent), distance = "squared_euclidean",
    negative_share = 1)
  expect_equal(admitted$models$dropped_negative_mass, negative,
    tolerance = 1e-10)
  expect_equal(admitted$models$negative_share, negative / positive,
    tolerance = 1e-10)
  expect_identical(admitted$models$rank_effective,
    sum(values > 1e-10 * max(abs(values))))
  expect_lt(admitted$centering_error, 1e-12)

  # A Euclidean model carries only round-off negative mass.
  clean <- model_basis(list(clean = D), distance = "squared_euclidean")
  expect_lt(clean$models$negative_share, 1e-12)
  expect_error(model_basis(list(zero = 0 * D), distance = "squared_euclidean"),
    "no geometry", class = "effect_input_error")
})

test_that("overlapping bases admit forms while the coefficient reader refuses", {
  coarse <- category_rdm(c(1, 1, 1, 0, 0, 0))
  fine <- category_rdm(c(1, 1, 1, 2, 2, 3))
  basis <- model_basis(list(coarse = coarse, fine = fine),
    distance = "squared_euclidean")
  expect_identical(basis$dimension, 2L)
  expect_identical(basis$span_overlap_rank, 1L)
  expect_identical(dim(basis$R), c(2L, 3L))
  expect_identical(crossform:::.validate_model_basis(basis), basis)
  refusal <- catch_refusal(crossform:::.model_geometry_model_side(basis))
  expect_identical(refusal$capability, "identified_model_coordinates")
  expect_identical(refusal$namespace, "model_coordinate")
  expect_identical(refusal$reasons[[1L]], "model_span_overlap")
  expect_match(refusal$reasons[[2L]], "span_overlap_rank = 1", fixed = TRUE)
  duplicate <- model_basis(list(a = coarse, b = coarse),
    distance = "squared_euclidean")
  expect_identical(duplicate$span_overlap_rank, 1L)
  expect_identical(duplicate$dimension, 1L)
  # The coefficient admission gate runs before planning or neural reads.
  expect_identical(catch_refusal(model_geometry(NULL, duplicate))$reasons[[1L]],
    "model_span_overlap")
  expect_identical(model_basis(list(fine = fine),
    distance = "squared_euclidean")$dimension, 2L)
})

test_that("the signature covers conditions, names, kinds, tolerance, ranks and factors", {
  category <- category_rdm(c(1, 1, 1, 0, 0, 0))
  D <- squared_euclidean(feature_matrix(5, 3L))
  build <- function(...) {
    model_basis(list(category = category, feature = D),
      distance = "squared_euclidean", ...)
  }
  base <- build(rank = 2)

  expect_identical(build(rank = 2), base)
  expect_true(.strong_sha256(sub("^model-basis-", "", base$signature)))
  variants <- list(
    build(rank = c(category = 2, feature = 1)),
    build(rank = NULL),
    build(tolerance = 1e-8, rank = 2),
    build(rank = 2, conditions = rev(basis_conditions)),
    model_basis(list(cat = category, feature = D),
      distance = "squared_euclidean", rank = 2),
    model_basis(list(feature = D, category = category),
      distance = "squared_euclidean", rank = 2)
  )
  for (variant in variants) {
    expect_false(identical(variant$signature, base$signature))
  }
  # Reversing the condition order permutes the rows of Q and nothing else.
  reversed <- variants[[4L]]
  expect_identical(rownames(reversed$Q), rev(basis_conditions))
  expect_equal(projector(reversed$Q[basis_conditions, ]), projector(base$Q),
    tolerance = 1e-12)

  # The model-coordinate effect space carries the signature: two relations
  # lowered through the same basis share one effect-space signature, and two
  # bases never do.
  expect_identical(base$effect_space$provenance$model_basis_signature,
    base$signature)
  expect_identical(build(rank = 2)$effect_space$signature,
    base$effect_space$signature)
  expect_false(identical(variants[[1L]]$effect_space$signature,
    base$effect_space$signature))
  expect_identical(base$extractor$effect_space, base$effect_space)
  expect_identical(base$extractor$estimator, "model_basis")
})

test_that("condition alignment follows rsa(): both axes by name, or unlabelled at the right size", {
  category <- category_rdm(c(1, 1, 1, 0, 0, 0))
  declared <- model_basis(list(category = category),
    distance = "squared_euclidean", conditions = basis_conditions)
  order <- c(3, 1, 6, 2, 5, 4)
  shuffled <- category[order, order]
  aligned <- model_basis(list(category = shuffled),
    distance = "squared_euclidean", conditions = basis_conditions)
  expect_identical(aligned, declared)

  bare <- unname(category)
  unlabelled <- model_basis(list(category = bare),
    distance = "squared_euclidean", conditions = basis_conditions)
  expect_identical(unlabelled, declared)
  expect_error(model_basis(list(category = bare), distance = "squared_euclidean"),
    "Supply `conditions`", class = "effect_input_error")
  expect_error(model_basis(list(category = bare), distance = "squared_euclidean",
    conditions = basis_conditions[1:5]), "unlabelled",
    class = "effect_input_error")
  half <- category
  colnames(half) <- NULL
  expect_error(model_basis(list(category = half), distance = "squared_euclidean"),
    "names only its rows", class = "effect_input_error")
  renamed <- category
  dimnames(renamed) <- list(basis_conditions, sub("face", "faces",
    basis_conditions))
  expect_error(model_basis(list(category = renamed),
    distance = "squared_euclidean"), "not a declared condition",
    class = "effect_input_error")
  expect_error(model_basis(list(category = category),
    distance = "squared_euclidean", conditions = basis_conditions[-1]),
    "not a declared condition", class = "effect_input_error")
  asymmetric <- category
  asymmetric[1, 2] <- 0.5
  expect_error(model_basis(list(category = asymmetric),
    distance = "squared_euclidean"), "not symmetric",
    class = "effect_input_error")
  expect_error(model_basis(list(category = category + diag(6)),
    distance = "squared_euclidean"), "nonzero diagonal",
    class = "effect_input_error")
  expect_error(model_basis(list(category = -category),
    distance = "squared_euclidean"), "negative entries",
    class = "effect_input_error")

  # An effect space lends its units to the model coordinates; a condition
  # space additionally yields the ingestion-path effect map.
  space <- effect_space(basis_conditions, basis_id = "condition-means",
    units = "percent-signal-change")
  typed <- model_basis(list(category = category),
    distance = "squared_euclidean", conditions = space)
  expect_identical(unique(unname(typed$effect_space$units)),
    "percent-signal-change")
  expect_null(typed$effect_map)
  conditions <- condition_space(basis_conditions, units = "arbitrary-BOLD")
  ingestion <- model_basis(list(category = category),
    distance = "squared_euclidean", conditions = conditions)
  expect_s3_class(ingestion$effect_map, "effect_condition_map")
  expect_identical(ingestion$effect_map$condition_space$signature,
    conditions$signature)
  expect_equal(unname(ingestion$effect_map$weights), unname(t(ingestion$Q)))
  expect_identical(ingestion$effect_map$effect_space, ingestion$effect_space)
  expect_identical(unique(unname(ingestion$effect_space$units)),
    "arbitrary-BOLD")
  mixed <- effect_space(basis_conditions, units = c(rep("a", 3), rep("b", 3)))
  expect_error(model_basis(list(category = category),
    distance = "squared_euclidean", conditions = mixed), "different units",
    class = "effect_input_error")
})

test_that("rank requests are one number, one per model, or named by model", {
  # Unequal category sizes give a resolved spectral gap. Positive tied cuts
  # now have their own refusal tests in test-predictive-kernel.R (T11).
  category <- category_rdm(c(1, 1, 1, 0, 0, 2))
  D <- squared_euclidean(feature_matrix(6, 3L))
  models <- list(category = category, feature = D)
  by_position <- model_basis(models, "squared_euclidean", rank = c(1, 2))
  by_name <- model_basis(models, "squared_euclidean",
    rank = c(feature = 2, category = 1))
  expect_identical(by_position, by_name)
  expect_identical(by_position$models$rank_effective, c(1L, 2L))
  expect_error(model_basis(models, "squared_euclidean", rank = c(1, 2, 3)),
    "one per model", class = "effect_input_error")
  expect_error(model_basis(models, "squared_euclidean", rank = 0),
    "positive whole number", class = "effect_input_error")
  expect_error(model_basis(models, "squared_euclidean",
    rank = c(category = 1, other = 2)), "must match the model names",
    class = "effect_input_error")
  expect_error(model_basis(models, "squared_euclidean", tolerance = 2),
    "below one", class = "effect_input_error")
})

test_that("saturation is recorded on the basis and refused for a shared fit", {
  conditions <- c("face", "body", "house", "tool")
  category <- category_rdm(c(1, 1, 0, 0), conditions)
  shape <- feature_matrix(7, 2L, conditions)
  saturated <- model_basis(list(category = category), "squared_euclidean",
    features = list(shape = shape))
  expect_identical(saturated$dimension, 3L)
  expect_equal(saturated$span_fraction, 1)
  expect_true(saturated$saturated)

  refusal <- catch_refusal(.model_basis_refuse_saturated(saturated))
  expect_s3_class(refusal, "effect_capability_refusal")
  expect_identical(refusal$capability, "model_geometry_test")
  expect_identical(refusal$namespace, "model_coordinate")
  expect_identical(refusal$reasons[[1L]], "model_span_saturated")
  expect_length(refusal$remedies, 3L)

  proper <- model_basis(list(category = category), "squared_euclidean")
  expect_false(proper$saturated)
  expect_identical(.model_basis_refuse_saturated(proper), proper)
})

test_that("the extractor lowers a relation into the projected geometry", {
  # Contract section 2.2: lowering commutes with the cross-partition
  # estimator for every component, so the lowered form equals Q' G Q read
  # from the full form. This is the package route of the oracle's law O1 on a
  # small fixture; the README fixture is the contract test's job (M7).
  set.seed(8)
  n <- 6L
  p <- 5L
  category <- category_rdm(c(1, 1, 1, 0, 0, 0))
  D <- squared_euclidean(feature_matrix(9, 2L))
  basis <- model_basis(list(category = category, feature = D),
    distance = "squared_euclidean", conditions = basis_conditions)
  Q <- basis$Q
  domain <- abstract_domain(p, id = "model-basis-lowering")
  betas <- lapply(1:3, function(run) {
    b <- matrix(rnorm(n * p), n, p)
    rownames(b) <- basis_conditions
    b
  })
  names(betas) <- paste0("run", 1:3)

  full <- relation(betas, effects = basis_conditions, domain = domain)
  lowered <- relation(betas, extract = basis$extractor, domain = domain)
  expect_identical(lowered$effect_space, basis$effect_space)
  expect_equal(relation_block(lowered, "run2", seq_len(p)),
    unname(t(Q) %*% betas$run2), tolerance = 1e-12, ignore_attr = TRUE)

  frame <- compile_frame(whole_brain(), domain)
  G <- materialize_geometry(plan_geometry(full, frame,
    cross_partitions(full, independence = "independent")))
  G_m <- materialize_geometry(plan_geometry(lowered, frame,
    cross_partitions(lowered, independence = "independent")))
  for (component in c("total", "coherent", "configuration")) {
    full_form <- .unsvec_symmetric(geometry_component(G, component)[1L, ], n)
    lowered_form <- .unsvec_symmetric(
      geometry_component(G_m, component)[1L, ], basis$dimension)
    expect_equal(lowered_form, t(Q) %*% full_form %*% Q, tolerance = 1e-12,
      ignore_attr = TRUE)
  }
  # And the lowered form is baseline invariant while the full form is not.
  shifted <- lapply(betas, function(b) b + rep(1, n) %*% t(rnorm(p)))
  full_shift <- relation(shifted, effects = basis_conditions, domain = domain)
  lowered_shift <- relation(shifted, extract = basis$extractor, domain = domain)
  G_shift <- materialize_geometry(plan_geometry(full_shift, frame,
    cross_partitions(full_shift, independence = "independent")))
  G_m_shift <- materialize_geometry(plan_geometry(lowered_shift, frame,
    cross_partitions(lowered_shift, independence = "independent")))
  expect_gt(max(abs(geometry_component(G_shift, "total") -
    geometry_component(G, "total"))), 1e-3)
  expect_lt(max(abs(geometry_component(G_m_shift, "total") -
    geometry_component(G_m, "total"))), 1e-10)
})

test_that("poorly conditioned models stay centered to round-off", {
  # The M2 review's blocker: a direct eigendecomposition of -HDH/2 leaves the
  # eigenvector of a small kept root uncentered by eps * cond^2, which at a
  # feature scale ratio of 1e-4 reached 1.1e-6 and tripped the invariant
  # assertion on ordinary inputs. Decomposing in an explicit basis of 1-perp
  # keeps |Q'1| at round-off whatever the conditioning.
  for (n in c(6L, 12L, 40L)) {
    conditions <- paste0("c", seq_len(n))
    for (ratio in c(1e-3, 1e-4, 1e-5)) {
      set.seed(n + round(-log10(ratio)))
      F <- cbind(rnorm(n), ratio * rnorm(n))
      rownames(F) <- conditions
      # The eigenvalue ratio is the square of the scale ratio, so the
      # tolerance is lowered to keep the small root rather than drop it.
      from_features <- model_basis(features = list(f = F), tolerance = 1e-12)
      from_rdm <- model_basis(list(f = squared_euclidean(F)),
        distance = "squared_euclidean", conditions = conditions,
        tolerance = 1e-12)
      for (basis in list(from_features, from_rdm)) {
        expect_identical(basis$dimension, 2L)
        expect_lt(basis$centering_error, 1e-12)
        expect_lt(max(abs(colSums(basis$Q))), 1e-12)
        expect_lt(max(abs(crossprod(basis$Q) - diag(2))), 1e-12)
        expect_gt(basis$basis_condition, 1 / ratio / 10)
      }
    }
  }
})

test_that("a model with no centered variation is an input error, not a bug", {
  constant <- cbind(rep(1, 6), rep(2, 6))
  rownames(constant) <- basis_conditions
  expect_error(model_basis(features = list(f = constant)),
    "constant across conditions", class = "effect_input_error")
  # A constant column beside a live one is dropped below tolerance.
  mixed <- cbind(rep(1, 6), feature_matrix(11, 1L)[, 1L])
  rownames(mixed) <- basis_conditions
  expect_identical(model_basis(features = list(f = mixed))$dimension, 1L)
})

test_that("feature rows and label orders align by name, and a clean model drops no mass", {
  F <- feature_matrix(12, 2L)
  declared <- model_basis(features = list(f = F), conditions = basis_conditions)
  permuted <- F[c(4, 2, 6, 1, 3, 5), , drop = FALSE]
  expect_identical(model_basis(features = list(f = permuted),
    conditions = basis_conditions), declared)
  # An RDM whose row and column labels are in different orders is symmetric
  # only once both axes follow the condition order.
  D <- squared_euclidean(F)
  crossed <- D[c(2, 1, 3, 4, 5, 6), c(6, 5, 4, 3, 2, 1)]
  expect_identical(model_basis(list(f = crossed), "squared_euclidean",
    conditions = basis_conditions),
    model_basis(list(f = D), "squared_euclidean", conditions = basis_conditions))
  expect_equal(projector(model_basis(list(f = crossed), "squared_euclidean",
    conditions = basis_conditions)$Q), projector(declared$Q), tolerance = 1e-12)
  # Roots within tolerance of zero are neither kept nor counted as dropped.
  clean <- model_basis(list(f = D), "squared_euclidean", negative_share = 0)
  expect_identical(clean$models$dropped_negative_mass, 0)
  expect_identical(clean$models$negative_share, 0)
  # A shared scale is carried; a per-coordinate scale is refused like units.
  scaled <- effect_space(basis_conditions, scale = 100)
  expect_identical(unique(unname(model_basis(list(f = D), "squared_euclidean",
    conditions = scaled)$effect_space$scale)), 100)
  uneven <- effect_space(basis_conditions, scale = 1:6)
  expect_error(model_basis(list(f = D), "squared_euclidean",
    conditions = uneven), "different scales", class = "effect_input_error")
})

test_that("derived fields are re-derived, so a forged basis cannot open a gate", {
  conditions <- c("face", "body", "house", "tool")
  category <- category_rdm(c(1, 1, 0, 0), conditions)
  full <- model_basis(list(category = category), "squared_euclidean",
    features = list(shape = feature_matrix(7, 2L, conditions)))
  proper <- model_basis(list(category = category), "squared_euclidean")

  forge <- function(basis, edit) {
    forged <- basis
    forged <- edit(forged)
    expect_error(format(forged), "derived fields|identity is inconsistent",
      class = "effect_contract_error")
  }
  forge(full, function(b) { b$saturated <- FALSE; b })
  forge(full, function(b) { b$span_fraction <- 9; b })
  forge(full, function(b) { b$basis_condition <- 1e300; b })
  forge(full, function(b) { b$baseline_invariant <- FALSE; b })
  forge(full, function(b) { b$centering_error <- 1; b })
  forge(full, function(b) { b$span_overlap_rank <- 1L; b })
  forge(full, function(b) { b$columns <- rev(b$columns); b })
  forge(full, function(b) { b$models$rank_effective <- c(2L, 1L); b })
  forge(full, function(b) { b$models$rank_clamped <- TRUE; b })
  forge(full, function(b) { colnames(b$R) <- rev(colnames(b$R)); b })
  forge(full, function(b) { colnames(b$Q) <- rev(colnames(b$Q)); b })
  # The saturation gate decides from the dimensions, and a saturated basis
  # with its flag flipped is rejected as forged rather than admitted.
  flipped <- full
  flipped$saturated <- FALSE
  expect_error(.model_basis_refuse_saturated(flipped), "derived fields",
    class = "effect_contract_error")
  expect_identical(.model_basis_refuse_saturated(proper), proper)
})

test_that("mutated model bases fail closed", {
  category <- category_rdm(c(1, 1, 1, 0, 0, 0))
  basis <- model_basis(list(category = category), distance = "squared_euclidean")

  tampered <- basis
  tampered$Q[1, 1] <- tampered$Q[1, 1] + 1e-3
  tampered$centering_error <- max(abs(colSums(tampered$Q)))
  expect_error(format(tampered), "not orthonormal and centered",
    class = "effect_contract_error")
  scaled <- basis
  scaled$R <- scaled$R * 2
  expect_error(format(scaled), "identity is inconsistent",
    class = "effect_contract_error")
  relabelled <- basis
  replacement <- if (grepl("0$", relabelled$signature)) "1" else "0"
  relabelled$signature <- sub("[0-9a-f]$", replacement, relabelled$signature)
  expect_error(format(relabelled), "identity is inconsistent",
    class = "effect_contract_error")
  truncated <- basis
  truncated$span_fraction <- NULL
  expect_error(format(truncated), "missing or noncanonical",
    class = "effect_input_error")
  swapped <- basis
  swapped$extractor$estimator <- "explicit"
  expect_error(format(swapped), "products are inconsistent",
    class = "effect_contract_error")
})

test_that("the basis prints its dimension, models, clamping and centering", {
  category <- category_rdm(c(1, 1, 1, 0, 0, 0))
  basis <- model_basis(list(category = category), distance = "squared_euclidean",
    rank = 2)
  expect_identical(format(basis), "model_basis<q = 1 of 5 centered; 1 model; span 0.20>")
  printed <- paste(capture.output(print(basis)), collapse = "\n")
  expect_match(printed, "rank_clamped", fixed = TRUE)
  expect_match(printed, "TRUE", fixed = TRUE)
  expect_match(printed, "baseline_invariant: TRUE", fixed = TRUE)
  expect_match(printed, "signature: model-basis-sha256:", fixed = TRUE)

  conditions <- c("face", "body", "house", "tool")
  full <- model_basis(list(category = category_rdm(c(1, 1, 0, 0), conditions)),
    "squared_euclidean", features = list(shape = feature_matrix(7, 2L, conditions)))
  expect_match(format(full), "(saturated)", fixed = TRUE)
})
