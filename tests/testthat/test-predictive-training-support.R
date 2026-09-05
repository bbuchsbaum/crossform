# Ancestry tests use tiny plans, but never need to execute neural geometry.
ps_manifest <- function() {
  list(id = "prediction-origins-v1", partitions = list(
    a = "raw-a", b = "raw-b", c = "raw-c", d = "raw-d", unused = "raw-unused"),
    independence = "independent",
    assumption = "Origins are independent runs; design and preprocessing were fixed externally.")
}
ps_fixture <- function(manifest = ps_manifest(), same_values = FALSE) {
  nm <- letters[1:3]
  B <- matrix(seq_len(12)/10, 3, dimnames = list(nm, NULL))
  betas <- stats::setNames(lapply(seq_len(5), function(i) {
    if (same_values) B else B * i
  }), c("a", "b", "c", "d", "unused"))
  domain <- abstract_domain(4, id = "origin-fixture")
  provenance <- if (is.null(manifest)) list() else list(observation_origins = manifest)
  rel <- relation(betas, domain = domain, provenance = provenance)
  at <- compile_frame(whole_brain(), domain)
  plan <- function(rel, left, right, weight = NULL, independence = "independent") {
    plan_geometry(rel, at, pairing(left, right, weight,
      independence = independence, generalizes_over = "run"))
  }
  list(rel = rel, betas = betas, at = at, domain = domain, plan = plan,
    train = plan(rel, "a", "b"), test = plan(rel, "c", "d"))
}
ps_support <- function(plan, upstream = list()) crossform:::.geometry_training_support(plan, upstream)
ps_admit <- function(train, test) crossform:::.geometry_support_independence(train, test)

test_that("[T26 T32] actual support excludes unused and zero-weight endpoints", {
  f <- ps_fixture()
  support <- ps_support(f$train)
  expect_identical(names(support$endpoints), c("a", "b"))
  expect_identical(support$observations, c("raw-a", "raw-b"))
  expect_identical(support$products, data.frame(left = "a", right = "b", weight = 1))
  f$rel$sources$unused$read <- function(features) stop("unused source was read")
  with_zero <- f$plan(f$rel, c("a", "a"), c("b", "unused"), weight = c(1, 0))
  expect_identical(ps_support(with_zero), support)
  # Replace unused numeric content while keeping the parent origin declaration.
  altered <- f$betas; altered$unused <- altered$unused + 100
  rel2 <- relation(altered, domain = f$domain,
    provenance = list(observation_origins = ps_manifest()))
  expect_identical(ps_support(f$plan(rel2, "a", "b")), support)
  # A used source revision is a fitting dependency.
  altered$a <- altered$a + .5
  rel3 <- relation(altered, domain = f$domain,
    provenance = list(observation_origins = ps_manifest()))
  expect_false(identical(ps_support(f$plan(rel3, "a", "b"))$signature, support$signature))
})

test_that("[T27 T28] origin identity distinguishes copies from independent identical values", {
  f <- ps_fixture(same_values = TRUE)
  train <- ps_support(f$train); test <- ps_support(f$test)
  expect_identical(train$endpoints$a$source_revision, test$endpoints$c$source_revision)
  expect_true(ps_admit(train, test))
  copied <- crossform:::.relation_subset(f$rel, c("a", "b"), c("copy-c", "copy-d"))
  copy_plan <- f$plan(copied, "copy-c", "copy-d")
  copied_support <- ps_support(copy_plan)
  expect_identical(copied_support$observations, train$observations)
  expect_identical(copied$sources$`copy-c`$origin, f$rel$sources$a$origin)
  refused <- catch_refusal(ps_admit(train, copied_support))
  expect_identical(refused$reasons[[1]], "overlapping_training_evaluation_origins")
  # A reprocessed source has a new content revision and the same origin.
  m <- ps_manifest(); m$partitions$c <- "raw-a"
  reprocessed <- ps_fixture(m)
  expect_false(identical(ps_support(reprocessed$train)$endpoints$a$source_revision,
    ps_support(reprocessed$test)$endpoints$c$source_revision))
  expect_identical(catch_refusal(ps_admit(ps_support(reprocessed$train),
    ps_support(reprocessed$test)))$reasons[[1]], "overlapping_training_evaluation_origins")
})

test_that("[T27 T30] neural independence and model evaluation independence are separate gates", {
  m <- ps_manifest(); m$partitions$b <- "raw-a"
  f <- ps_fixture(m)
  expect_identical(catch_refusal(ps_admit(ps_support(f$train), ps_support(f$test)))$reasons[[1]],
    "overlapping_neural_origins")
  f <- ps_fixture()
  shared <- f$plan(f$rel, "a", "d")
  expect_identical(catch_refusal(ps_admit(ps_support(f$train), ps_support(shared)))$reasons[[1]],
    "overlapping_training_evaluation_origins")
  m <- ps_manifest(); m$independence <- "undeclared"
  f <- ps_fixture(m)
  expect_identical(catch_refusal(ps_admit(ps_support(f$train), ps_support(f$test)))$reasons[[1]],
    "independence_undeclared")
  f <- ps_fixture()
  undeclared <- f$plan(f$rel, "c", "d", independence = NULL)
  expect_identical(catch_refusal(ps_admit(ps_support(f$train), ps_support(undeclared)))$reasons[[1]],
    "independence_undeclared")
})

test_that("[T29 T30] unknown external lineage cannot become independent by a flag", {
  f <- ps_fixture(NULL)
  train <- ps_support(f$train); test <- ps_support(f$test)
  expect_false(train$known)
  expect_identical(train$observations, character())
  expect_identical(catch_refusal(ps_admit(train, test))$reasons[[1]], "unknown_observation_origins")
  forged <- train; forged$known <- TRUE
  expect_error(crossform:::.validate_geometry_training_support(forged),
    "origin status", class = "effect_contract_error")
  known <- ps_fixture()
  m <- ps_manifest(); m$id <- "unrelated-dataset"
  other <- ps_fixture(m)
  expect_identical(catch_refusal(ps_admit(ps_support(known$train), ps_support(other$test)))$reasons[[1]],
    "origin_scope_mismatch")
})

test_that("[T29 T32] upstream selection and preprocessing origins belong to fitting support", {
  f <- ps_fixture()
  selection <- ps_support(f$test)
  fitted <- ps_support(f$train, upstream = list(selection))
  expect_identical(fitted$observations, c("raw-a", "raw-b", "raw-c", "raw-d"))
  expect_identical(catch_refusal(ps_admit(fitted, ps_support(f$test)))$reasons[[1]],
    "overlapping_training_evaluation_origins")
  m <- ps_manifest(); m$dependencies <- list(a = "raw-c")
  f <- ps_fixture(m)
  expect_identical(ps_support(f$train)$observations, c("raw-a", "raw-b", "raw-c"))
  expect_identical(catch_refusal(ps_admit(ps_support(f$train), ps_support(f$test)))$reasons[[1]],
    "overlapping_training_evaluation_origins")
})

test_that("[T28 T29] lowering, serialization and row repacking preserve source ancestry", {
  f <- ps_fixture()
  F <- matrix(c(1, -1, 0), 3, dimnames = list(letters[1:3], "contrast"))
  basis <- model_basis(features = F, conditions = f$rel$effect_space)
  low <- crossform:::.relation_lowered(f$rel, basis)
  for (p in f$rel$partitions) expect_identical(low$sources[[p]]$origin, f$rel$sources[[p]]$origin)
  support <- ps_support(f$train)
  expect_identical(unserialize(serialize(support, NULL)), support)
  expect_identical(crossform:::.validate_geometry_training_support(support), support)
  # Repacking observation rows with labels does not rename their raw origins.
  repacked <- lapply(f$betas, function(B) B[c(3, 1, 2), , drop = FALSE])
  repacked <- relation(repacked, effects = f$rel$effect_space, domain = f$domain,
    provenance = f$rel$provenance)
  expect_identical(ps_support(f$plan(repacked, "a", "b")), support)
  modified <- support; modified$observations <- "raw-d"
  expect_error(crossform:::.validate_geometry_training_support(modified), class = "effect_contract_error")
  modified <- f$rel; modified$provenance$observation_origins <- NULL
  expect_error(crossform:::.validate_relation(modified), "cannot discard")
  modified <- f$rel; modified$sources$a$origin$observations <- "raw-c"
  expect_error(crossform:::.validate_relation(modified), "identity is inconsistent")
})

test_that("[T27 T30] malformed origin declarations fail before execution", {
  for (change in list(
    function(m) {m$partitions$a <- character(); m},
    function(m) {m$partitions$a <- c("raw-a", "raw-a"); m},
    function(m) {m$independence <- TRUE; m},
    function(m) {m$assumption <- ""; m},
    function(m) {m$dependencies <- list(absent = "raw-a"); m},
    function(m) {m$partitions$a <- NULL; m}
  )) expect_error(ps_fixture(change(ps_manifest())), class = "effect_input_error")
  f <- ps_fixture(); forged <- f$rel$provenance$observation_origins
  forged$partitions$a <- "raw-c"
  expect_error(crossform:::.validate_origin_manifest(forged), class = "effect_contract_error")
})
