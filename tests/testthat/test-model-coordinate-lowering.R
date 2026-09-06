# Lowering a relation into model coordinates (ticket M3 of the
# model-coordinate epic).
#
# `design/model-coordinate-geometry-contract.md` section 2: the model basis is
# an extractor on the betas path and an effect map on the ingestion path; the
# ordinary plan and executor then produce the complete geometry in model
# coordinates, `S_x = Q' G_x Q`, and lowering commutes with every stage that
# is linear on the effect axis and with nothing else. The first test is the
# package-route law of section 2.2 on the README fixture, promoted from
# `design/oracles/model-coordinate-package-route.R` section R1; the others
# cover the two entries, the population signature rule of section 2.1(a), and
# the gate of section 2.2 against the internal evidence-task compiler, which
# is the only path that can reach it while the edge normalizers stay
# unexported.

lowering_conditions <- c("face", "body", "house", "tool")

lowering_category <- function(conditions = lowering_conditions) {
  animate <- c(1, 1, 0, 0)
  D <- outer(animate, animate, function(a, b) (a - b)^2)
  dimnames(D) <- list(conditions, conditions)
  D
}

unsvec <- function(v, q) crossform:::.unsvec_symmetric(v, q)
upper <- function(S) S[upper.tri(S, diag = TRUE)]

test_that("the lowered form is the projected full form on the README fixture", {
  example <- example_fmri_effects()
  rel <- example$fit$relation
  basis <- model_basis(list(category = example$model_rdm),
    distance = "squared_euclidean", rank = 2, conditions = rel$effect_space)
  # The README category RDM has a rank-1 Gram: the requested rank is clamped
  # and no uncentered direction can enter Q (contract section 1.1).
  expect_identical(basis$models$rank_effective, 1L)
  expect_true(basis$models$rank_clamped)
  expect_identical(basis$dimension, 1L)
  Q <- unname(basis$Q)
  q_full <- length(rel$effects)

  over <- cross_partitions(rel, independence = "independent",
    generalizes_over = "run")
  G <- materialize_geometry(plan_geometry(rel, at = example$frame, over = over))

  betas <- lapply(rel$partitions, function(p) {
    relation_block(rel, p, seq_len(rel$n_features))
  })
  names(betas) <- rel$partitions
  rel_m <- relation(betas, extract = basis, domain = example$domain)
  expect_identical(rel_m$effect_space, basis$effect_space)
  over_m <- cross_partitions(rel_m, independence = "independent",
    generalizes_over = "run")
  G_m <- materialize_geometry(plan_geometry(rel_m, at = example$frame,
    over = over_m))

  for (component in c("total", "coherent", "configuration")) {
    full <- geometry_component(G, component)
    lowered <- geometry_component(G_m, component)
    projected <- t(apply(full, 1, function(v) {
      upper(t(Q) %*% unsvec(v, q_full) %*% Q)
    }))
    read <- t(apply(lowered, 1, function(v) upper(unsvec(v, basis$dimension))))
    expect_lt(max(abs(projected - read)), 1e-12)
  }

  # The existing readers run unchanged on the lowered form.
  spectrum <- geometry_spectrum(G_m)
  latent <- latent_geometry(G_m)
  expect_s3_class(spectrum, "effect_spectrum_view")
  expect_s3_class(latent, "effect_latent_geometry")
  expect_identical(ncol(spectrum$values), basis$dimension)
})

test_that("relation() accepts the basis itself and lowers through its extractor", {
  set.seed(31)
  category <- lowering_category()
  basis <- model_basis(list(category = category), "squared_euclidean")
  domain <- abstract_domain(3L, id = "lowering-entry")
  betas <- lapply(1:2, function(run) {
    b <- matrix(rnorm(12), 4L, 3L)
    rownames(b) <- lowering_conditions
    b
  })
  names(betas) <- c("run1", "run2")

  by_basis <- relation(betas, extract = basis, domain = domain)
  by_extractor <- relation(betas, extract = basis$extractor, domain = domain)
  expect_identical(by_basis, by_extractor)
  expect_identical(by_basis$effects, "mc1")
  expect_identical(by_basis$effect_space$provenance$model_basis_signature,
    basis$signature)
  expect_true(crossform:::.is_model_coordinate_space(by_basis$effect_space))
  expect_false(crossform:::.is_model_coordinate_space(
    relation(betas, domain = domain)$effect_space))
  expect_equal(relation_block(by_basis, "run1", 1:3),
    unname(t(basis$Q) %*% betas$run1), tolerance = 1e-12, ignore_attr = TRUE)

  # Labelled betas are aligned to the basis's conditions by name; unlabelled
  # betas are taken by position; a label mismatch is refused.
  permuted <- lapply(betas, function(b) b[c(3, 1, 4, 2), , drop = FALSE])
  aligned <- relation(permuted, extract = basis, domain = domain)
  expect_equal(relation_block(aligned, "run1", 1:3),
    relation_block(by_basis, "run1", 1:3), tolerance = 1e-12)
  unlabelled <- lapply(permuted, unname)
  positional <- relation(unlabelled, extract = basis, domain = domain)
  expect_equal(relation_block(positional, "run1", 1:3),
    unname(t(basis$Q) %*% permuted$run1), tolerance = 1e-12,
    ignore_attr = TRUE)
  renamed <- lapply(betas, function(b) { rownames(b)[[1L]] <- "faces"; b })
  expect_error(relation(renamed, extract = basis, domain = domain),
    "does not carry the observations", class = "effect_input_error")

  # A declared effect space that is not the basis's is a contradiction.
  expect_error(relation(betas, extract = basis, domain = domain,
    effects = effect_space("other")), "incompatible",
    class = "effect_contract_error")
})

test_that("plan_relation() lowers through the basis's effect map on the ingestion path", {
  fixture <- relation_plan_fixture("qr", "cell")
  conditions <- fixture$conditions
  labels <- conditions$coordinates
  category <- outer(c(1, 0, 0), c(1, 0, 0), function(a, b) (a - b)^2)
  dimnames(category) <- list(labels, labels)
  basis <- model_basis(list(face = category), "squared_euclidean",
    conditions = conditions)
  expect_s3_class(basis$effect_map, "effect_condition_map")

  # The condition means themselves, to lower by hand afterwards.
  identity <- diag(3)
  dimnames(identity) <- list(labels, labels)
  means <- effect_map(identity, conditions)
  study <- fixture$plan$study
  plan_means <- plan_relation(study, fixture$model, means, fixture$observation)
  plan_lowered <- plan_relation(study, fixture$model, basis,
    fixture$observation)
  plan_by_map <- plan_relation(study, fixture$model, basis$effect_map,
    fixture$observation)
  expect_identical(plan_lowered$relation_plan_id, plan_by_map$relation_plan_id)
  expect_false(identical(plan_lowered$relation_plan_id,
    plan_means$relation_plan_id))

  fit_means <- estimate_relation(plan_means)
  fit_lowered <- estimate_relation(plan_lowered)
  expect_identical(fit_lowered$relation$effect_space, basis$effect_space)
  expect_true(crossform:::.is_model_coordinate_space(
    fit_lowered$relation$effect_space))
  for (partition in fit_means$relation$partitions) {
    full <- relation_block(fit_means, partition, 1:5)
    lowered <- relation_block(fit_lowered, partition, 1:5)
    expect_equal(lowered, unname(t(basis$Q) %*% full), tolerance = 1e-10,
      ignore_attr = TRUE)
  }

  # A basis built without a condition space has no ingestion-path product.
  bare <- model_basis(list(face = category), "squared_euclidean")
  expect_null(bare$effect_map)
  expect_error(plan_relation(study, fixture$model, bare, fixture$observation),
    "carries no effect map", class = "effect_input_error")
  # And a basis on a different condition space is the existing refusal.
  other <- model_basis(list(face = category), "squared_euclidean",
    conditions = condition_space(labels, basis_id = "other-basis"))
  refusal <- catch_refusal(plan_relation(study, fixture$model, other,
    fixture$observation))
  expect_identical(refusal$capability, "valid_effect_lowering")
})

test_that("two participants lowered through one basis share an effect space, two bases never do", {
  conditions <- c("face", "house")
  category <- outer(c(1, 0), c(1, 0), function(a, b) (a - b)^2)
  dimnames(category) <- list(conditions, conditions)
  basis <- model_basis(list(category = category), "squared_euclidean")
  other <- model_basis(list(category = category), "squared_euclidean",
    tolerance = 1e-8)
  expect_false(identical(basis$signature, other$signature))
  expect_false(identical(basis$effect_space$signature,
    other$effect_space$signature))

  subject <- function(id, features, basis) {
    domain <- abstract_domain(features,
      coordinates = cbind(x = seq_len(features) - 1),
      feature_ids = paste0("f", seq_len(features)), id = id)
    values <- function(divisor) {
      matrix(seq_len(2 * features) / (features * divisor), 2, features,
        dimnames = list(conditions, NULL))
    }
    relation <- relation(list(run1 = values(1), run2 = values(2)),
      extract = basis, domain = domain)
    plan_geometry(relation,
      compile_frame(voxelwise(normalization = "conservative"), domain),
      cross_partitions(relation))
  }
  carrier <- function(features) {
    anatomical_transport(native_coords = cbind(seq_len(features) - 1),
      group_coords = cbind(c(0, 3)), semantics = "budget")
  }
  subjects <- list(s01 = subject("s01", 5L, basis), s02 = subject("s02", 6L, basis))
  transport <- list(s01 = carrier(5L), s02 = carrier(6L))
  plan <- plan_population(subjects, transport)
  expect_s3_class(plan, "effect_population_plan")
  expect_identical(subjects$s01$task$left_relation$effect_space$signature,
    subjects$s02$task$left_relation$effect_space$signature)

  mixed <- list(s01 = subject("s01", 5L, basis), s02 = subject("s02", 6L, other))
  refusal <- catch_refusal(plan_population(mixed, transport))
  expect_s3_class(refusal, "effect_capability_refusal")
  expect_identical(refusal$capability, "common_experimental_space")
  expect_true(any(grepl("disagreeing_subject:s02", refusal$reasons,
    fixed = TRUE)))
})

test_that("nonlinear edge stages refuse model-coordinate lowering, linear ones admit it", {
  set.seed(32)
  category <- lowering_category()
  basis <- model_basis(list(category = category), "squared_euclidean")
  domain <- abstract_domain(3L, id = "lowering-gate")
  betas <- lapply(1:2, function(run) {
    b <- matrix(rnorm(12), 4L, 3L)
    rownames(b) <- lowering_conditions
    b
  })
  names(betas) <- c("run1", "run2")
  lowered <- relation(betas, extract = basis, domain = domain)
  plain <- relation(betas, domain = domain)
  over <- cross_partitions(lowered, independence = "independent")
  compile <- crossform:::.compile_effect_evidence_task

  refuses <- function(...) {
    refusal <- catch_refusal(compile(lowered, over, ...))
    expect_s3_class(refusal, "effect_capability_refusal")
    expect_identical(refusal$capability, "model_coordinate_lowering")
    expect_identical(refusal$namespace, "model_coordinate")
    expect_identical(refusal$reasons[[1L]], "nonlinear_edge_transform")
    invisible(refusal)
  }
  expect_match(refuses(normalizer = crossform:::correlation())$reasons[[2L]],
    "correlation", fixed = TRUE)
  expect_match(refuses(normalizer = crossform:::cosine())$reasons[[2L]],
    "cosine", fixed = TRUE)
  # A Fisher transform is admitted only over correlation-valued edges, so
  # the operation plan pairs it with `correlation()` and the gate names both.
  both <- refuses(normalizer = crossform:::correlation(),
    transform = crossform:::fisher_z())
  expect_length(both$reasons, 3L)
  expect_match(both$reasons[[3L]], "fisher_z", fixed = TRUE)
  expect_length(both$remedies, 2L)
  expect_match(refuses(transform = crossform:::rank_edges())$reasons[[2L]],
    "rank_edges", fixed = TRUE)

  # Linear stages are admitted, and the gate never touches an ordinary
  # relation whatever its stages.
  expect_s3_class(compile(lowered, over), "effect_evidence_task")
  expect_s3_class(compile(lowered, over, normalizer = crossform:::covariance()),
    "effect_evidence_task")
  expect_s3_class(compile(plain, cross_partitions(plain,
    independence = "independent"), normalizer = crossform:::correlation()),
    "effect_evidence_task")

  # The gate is keyed to the effect space's identity, not to a class: a space
  # that merely borrows the basis id without a basis signature is not lowered,
  # while one carrying a well-formed digest is treated as lowered. The check
  # is conservative by design (it only ever triggers refusals), and this pins
  # that it does not resolve the digest against a basis.
  impostor <- effect_space("mc1", basis_id = "model-coordinates")
  expect_false(crossform:::.is_model_coordinate_space(impostor))
  well_formed <- effect_space("mc1", basis_id = "model-coordinates",
    provenance = list(model_basis_signature = paste0("model-basis-sha256:",
      strrep("a", 64))))
  expect_true(crossform:::.is_model_coordinate_space(well_formed))
})

test_that("a relation lowered in place and one built from betas share one identity", {
  set.seed(33)
  category <- lowering_category()
  basis <- model_basis(list(category = category), "squared_euclidean")
  domain <- abstract_domain(3L, id = "lowering-identity")
  betas <- lapply(1:2, function(run) {
    b <- matrix(rnorm(12), 4L, 3L)
    rownames(b) <- lowering_conditions
    b
  })
  names(betas) <- c("run1", "run2")
  original <- relation(betas, effects = lowering_conditions, domain = domain)
  in_place <- crossform:::.relation_lowered(original, basis)
  from_betas <- relation(betas, extract = basis, domain = domain)
  expect_identical(in_place$extractors, from_betas$extractors)
  expect_identical(crossform:::.relation_family_identity(in_place),
    crossform:::.relation_family_identity(from_betas))
  frame <- compile_frame(whole_brain(), domain)
  public_plan <- plan_geometry(from_betas, frame,
    cross_partitions(from_betas, independence = "independent"))
  reading <- model_geometry(plan_geometry(original, frame,
    cross_partitions(original, independence = "independent")), basis)
  expect_identical(reading$lowered_scientific_plan_id,
    public_plan$scientific_plan_id)
  # A non-identity extractor composes, and the composed map keeps the
  # observation names of the map it composed with.
  fit_relation <- example_fmri_effects()$fit$relation
  composed_basis <- model_basis(list(category = example_fmri_effects()$model_rdm),
    "squared_euclidean", conditions = fit_relation$effect_space)
  composed <- crossform:::.relation_lowered(fit_relation, composed_basis)
  expect_identical(composed$extractors[[1L]]$estimator, "model_basis")
  expect_identical(rownames(composed$extractors[[1L]]$map), "mc1")
  expect_equal(unname(composed$extractors[[1L]]$map),
    unname(t(composed_basis$Q) %*% fit_relation$extractors[[1L]]$map),
    tolerance = 1e-12)
})
