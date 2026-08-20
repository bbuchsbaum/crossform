# Population views --- the reader verbs of `population-form-v1` over an
# estimated group form (ticket E6).
#
# The acceptance test is the additive identity of section 2 read through the
# group fit: transport is row-stochastic, so each participant's transported
# ledger sums --- over the group nodes AND the sink --- to its native total,
# and the fit is linear in the response, so the group coefficients sum the
# same way. `sum_u Theta_{j,u} = Theta_{j,Omega}` is therefore an *identity*,
# not an approximation, and the oracle below computes the right-hand side the
# long way: sum the transported responses over group nodes, refit the group
# model on those sums, and compare. Nothing in that route calls a view.
#
# The other half of the file is the refusal boundary. A population view is a
# combination in query space, and there are exactly two ways one can fail to
# be honest: the query is outside the span of what was estimated (nothing to
# combine), or the response is not linear along the basis axis so the
# combination does not estimate the combined query (`unit_budget`). Both are
# tested from the outside, through `catch_refusal()`, and both name a remedy.

pvw_effects <- function() {
  effect_space(c("face", "house", "tool"), basis_id = "pop-views:v1")
}

pvw_subject <- function(id, features) {
  domain <- abstract_domain(features,
    coordinates = cbind(x = seq_len(features) - 1),
    feature_ids = paste0("f", seq_len(features)), id = id)
  # Signed and participant-specific: a conservative ledger is signed by
  # construction, and a fixture of nonnegative numbers would let a sign bug in
  # the aggregate pass unnoticed.
  values <- function(divisor) {
    set.seed(sum(as.integer(charToRaw(id))) + features +
      round(1000 * divisor))
    matrix(stats::rnorm(3 * features) / divisor, 3, features,
      dimnames = list(c("face", "house", "tool"), NULL))
  }
  relation <- relation(
    list(run1 = values(1), run2 = values(1.3), run3 = values(0.8)),
    effects = pvw_effects(), domain = domain
  )
  plan_geometry(relation, compile_frame(voxelwise(), domain),
    cross_partitions(relation))
}

# Radius 1.5 around centres 0, 4 and 9 leaves native nodes at 2, 6, 7 and 11
# out of reach of every centre, so the sink carries real mass and the additive
# identity is a test of something rather than of an always-empty column.
pvw_carrier <- function(features, semantics = "budget") {
  anatomical_transport(
    native_coords = cbind(seq_len(features) - 1),
    group_coords = cbind(c(0, 4, 9)),
    semantics = semantics, radius = 1.5
  )
}

pvw_sizes <- c(s01 = 6L, s02 = 8L, s03 = 10L, s04 = 12L)

pvw_plan <- function(semantics = "budget", ...) {
  plan_population(
    stats::setNames(lapply(names(pvw_sizes), function(id)
      pvw_subject(id, pvw_sizes[[id]])), names(pvw_sizes)),
    lapply(pvw_sizes, pvw_carrier, semantics = semantics),
    ...
  )
}

# All three pairwise difference contrasts. This bank spans exactly the
# quadratic forms of zero-sum contrasts over three effects, which is what
# makes `rdm()` and `rsa()` derivable and every uncentred contrast refused.
pvw_bank <- function() {
  rbind(`face-house` = c(1, -1, 0), `face-tool` = c(1, 0, -1),
    `house-tool` = c(0, 1, -1))
}

pvw_fit <- function(plan = pvw_plan(), queries = pvw_bank(), ...) {
  estimate_population(plan, queries, ...)
}

pvw_animacy <- function() {
  matrix(c(0, 0, 1, 0, 0, 1, 1, 1, 0), 3, 3,
    dimnames = list(c("face", "house", "tool"),
      c("face", "house", "tool")))
}

# The oracle for the additive identity: the group coefficient of the summed
# response, computed without any view and without any population verb.
pvw_global <- function(plan, fit, term = "(Intercept)") {
  responses <- apply(fit$values, c(3L, 2L), sum)
  qr.coef(qr(plan$model$matrix), responses)[term, ]
}


# The additive identity ------------------------------------------------------

test_that("the group ledger adds up over the group nodes and the sink", {
  plan <- pvw_plan()
  fit <- pvw_fit(plan)

  # Section 2 through the fit. `apply(..., 2, sum)` runs over the node axis,
  # sink included -- deleting the sink row is exactly the failure the contract
  # measures at 48.1 % of the native total on its own fixture.
  summed <- apply(fit$coefficients[, , "(Intercept)"], 2L, sum)
  global <- pvw_global(plan, fit)
  expect_lt(max(abs(summed - global)), 1e-12)

  # The responses that were summed really are the participants' native
  # totals, so the identity is the contract's law and not an accident of the
  # fixture.
  expect_lt(max(abs(apply(fit$values, c(3L, 2L), sum) -
    fit$receipt$native_total)), 1e-12)

  # And through the verb, over a genuine partition of the group nodes.
  ledger <- contribution(fit, by = c("anterior", "anterior", "posterior"))
  expect_lt(max(abs(colSums(ledger$values) - global)), 1e-12)

  # Dropping the sink breaks it, which is why it is a row.
  without_sink <- ledger$values[!ledger$index$sink, , drop = FALSE]
  expect_gt(max(abs(colSums(without_sink) - global)), 1e-6)
})

test_that("the identity holds for a native ledger and on the form basis", {
  plan <- pvw_plan()

  for (component in c("coherent", "configuration")) {
    fit <- pvw_fit(plan, component = component)
    ledger <- contribution(fit, by = c("a", "a", "b"))
    expect_lt(max(abs(colSums(ledger$values) - pvw_global(plan, fit))), 1e-12)
  }

  # The complete-form route holds no `$values`, so the oracle here is the
  # query-first executor: the two bases estimate one population, and the
  # coordinate-form read of a contrast is the query-first coefficient of it.
  form <- materialize_population(plan)
  queried <- pvw_fit(plan)
  energy <- contrast_energy(form, c(1, -1, 0))
  expect_lt(max(abs(drop(energy$values) -
    queried$coefficients[, "face-house", "(Intercept)"])), 1e-12)

  aggregate <- contribution(energy, by = c("a", "a", "b"))
  expect_lt(abs(sum(aggregate$values) -
    sum(queried$coefficients[, "face-house", "(Intercept)"])), 1e-12)
})


# What the verbs read --------------------------------------------------------

test_that("a contrast in the bank is selected, and a rescaled one is exact", {
  fit <- pvw_fit()

  energy <- contrast_energy(fit, c(face = 1, house = -1, tool = 0))
  expect_identical(dim(energy$values), c(4L, 1L))
  expect_equal(drop(energy$values),
    fit$coefficients[, "face-house", "(Intercept)"],
    ignore_attr = TRUE)
  expect_identical(energy$receipt$basis$route, "selection")

  # A contrast energy is quadratic in the contrast, so `2c` is `4` times `c`.
  # The span solve reaches it, and reaching it is a *scaled* selection rather
  # than a selection because the weight is four rather than one.
  doubled <- contrast_energy(fit, c(2, -2, 0))
  expect_equal(drop(doubled$values), 4 * drop(energy$values))
  expect_identical(doubled$receipt$basis$route, "scaled_selection")
})

test_that("a genuine recombination is the estimate of the combined query", {
  # The load-bearing claim of the whole file, tested the only way a
  # commutation claim can be: estimate the combined query directly and
  # compare. `(1, 1, -2)` is zero-sum, so the pairwise bank spans it, but no
  # bank row is proportional to it -- this is a real three-column mixture.
  plan <- pvw_plan()
  fit <- pvw_fit(plan)

  combined <- contrast_energy(fit, c(1, 1, -2))
  expect_identical(combined$receipt$basis$route, "recombination")
  direct <- estimate_population(plan, rbind(`1+1-2` = c(1, 1, -2)))
  expect_lt(max(abs(drop(combined$values) -
    direct$coefficients[, "1+1-2", "(Intercept)"])), 1e-12)

  # And the same for RSA, whose every column is a mixture of the three pairs.
  fitted <- rsa(fit, models = list(animacy = pvw_animacy()))
  operators <- fitted$receipt$basis$coefficients
  expect_gt(max(colSums(operators != 0)), 1L)
  reference <- estimate_population(plan, pvw_bank())
  expect_lt(max(abs(fitted$values -
    reference$coefficients[, , "(Intercept)"] %*% operators)), 1e-12)
})

test_that("the pairwise bank reaches every distance and every RSA fit", {
  fit <- pvw_fit()

  distances <- rdm(fit)
  expect_identical(as.character(distances$columns$column),
    c("face - house", "face - tool", "house - tool"))
  expect_equal(unname(distances$values),
    unname(fit$coefficients[, , "(Intercept)"]))
  expect_identical(distances$receipt$basis$route, "selection")

  # RSA read independently: the coefficient at a node is the OLS fit of the
  # model design on that node's distances, computed here in plain R from the
  # group RDM rather than from anything the RSA verb touched.
  fitted <- rsa(fit, models = list(animacy = pvw_animacy()))
  design <- cbind(`(Intercept)` = 1,
    animacy = pvw_animacy()[lower.tri(pvw_animacy())])
  oracle <- t(qr.coef(qr(design), t(distances$values)))
  expect_equal(unname(fitted$values), unname(oracle))
  expect_identical(as.character(fitted$columns$role),
    c("intercept", "model"))

  # A pair the bank cannot reach is refused rather than approximated.
  narrow <- pvw_fit(queries = pvw_bank()[1L, , drop = FALSE])
  refusal <- catch_refusal(rdm(narrow))
  expect_identical(refusal$capability, "population_query_in_bank")
  expect_true(any(grepl("house - tool", refusal$reasons, fixed = TRUE)))
  expect_true(any(grepl("Re-estimate", refusal$remedies)))
})

test_that("an uncentred contrast is outside a zero-sum bank", {
  fit <- pvw_fit()

  # `(1, 1, -2)` is zero-sum, so the pairwise bank reaches it exactly even
  # though no bank row is proportional to it. This is the fact that makes the
  # span test the right test and set membership the wrong one.
  centred <- contrast_energy(fit, c(1, 1, -2))
  expect_identical(centred$receipt$basis$route, "recombination")
  expect_true(all(is.finite(centred$values)))

  refusal <- catch_refusal(contrast_energy(fit, c(1, 1, 0)))
  expect_identical(refusal$capability, "population_query_in_bank")
  expect_identical(refusal$namespace, "population_views")
  expect_true("query_outside_estimated_bank" %in% refusal$reasons)

  # The complete form spans every symmetric query, so the same contrast that
  # a bank refuses is read off a form without complaint.
  form <- materialize_population(pvw_plan())
  uncentred <- contrast_energy(form, c(1, 1, 0))
  expect_identical(uncentred$receipt$basis$route, "coordinate_form")
  expect_true(all(is.finite(uncentred$values)))
})

test_that("`unit_budget` reads one estimated query and refuses a mixture", {
  plan <- pvw_plan(normalization = "unit_budget")
  fit <- pvw_fit(plan)

  # Selection is exact under any normalization: it reads one estimated
  # column, and the column is what it is.
  energy <- contrast_energy(fit, c(1, -1, 0))
  expect_equal(drop(energy$values),
    fit$coefficients[, "face-house", "(Intercept)"], ignore_attr = TRUE)
  expect_identical(energy$receipt$basis$route, "selection")
  expect_s3_class(rdm(fit), "effect_population_view")

  # A multiple of an estimated query is the *same* share: the ledger and its
  # divisor scale together. So the weight is forced to one, and the values are
  # bit-identical to the selection rather than four times it.
  doubled <- contrast_energy(fit, c(2, -2, 0))
  expect_identical(doubled$values, energy$values)
  expect_identical(doubled$receipt$basis$route, "selection")

  # A mixture of estimated queries carries no denominator, and the numbers the
  # span solve would have returned are not close to anything -- so the refusal
  # is preventing a wrong answer, not a rounding difference.
  refusal <- catch_refusal(rsa(fit, models = list(animacy = pvw_animacy())))
  expect_identical(refusal$capability, "population_view_query_linearity")
  expect_true("unit_budget_divisor_is_query_specific" %in% refusal$reasons)
  expect_true(any(grepl("normalization = \\\"none\\\"", refusal$remedies)))

  mixed <- catch_refusal(contrast_energy(fit, c(1, 1, -2)))
  expect_identical(mixed$capability, "population_view_query_linearity")

  # How wrong the admitted alternative would be. The span coefficients are
  # exact for `normalization = "none"`, so applying them here measures the
  # gap between the two estimands rather than a numerical wobble.
  operators <- contrast_energy(pvw_fit(pvw_plan()),
    c(1, 1, -2))$receipt$basis$coefficients
  would_be <- fit$coefficients[, , "(Intercept)"] %*% operators
  honest <- estimate_population(plan,
    rbind(`1+1-2` = c(1, 1, -2)))$coefficients[, "1+1-2", "(Intercept)"]
  expect_gt(max(abs(drop(would_be) - honest)), 1)
})


# Refusals at the edges ------------------------------------------------------

test_that("estimation-time choices are not view-time arguments", {
  fit <- pvw_fit()

  expect_identical(catch_refusal(rdm(fit, component = "coherent"))$capability,
    "population_component_fixed_at_estimation")
  expect_identical(
    catch_refusal(rsa(fit, models = list(a = pvw_animacy()),
      component = "total"))$capability,
    "population_component_fixed_at_estimation")
  expect_identical(
    catch_refusal(rdm(fit, normalize = "correlation"))$capability,
    "guaranteed_psd")
  expect_identical(
    catch_refusal(contrast_energy(fit, c(1, -1, 0),
      remove_univariate = TRUE))$capability,
    "nondestructive_decomposition")
})

test_that("contribution refuses density and refuses to place the sink", {
  density <- pvw_fit(pvw_plan(semantics = "density"))
  refusal <- catch_refusal(contribution(density, by = c("a", "a", "b")))
  expect_identical(refusal$capability, "budget_semantics")
  expect_true("density_conserves_no_budget" %in% refusal$reasons)

  fit <- pvw_fit()
  named_sink <- catch_refusal(
    contribution(fit, by = c("a", "a", "<sink>"))
  )
  expect_identical(named_sink$capability, "sink_is_not_a_territory")

  # Every group node must belong to a territory: a missing label would drop
  # part of the budget the aggregate claims to partition.
  expect_error(contribution(fit, by = c("a", NA, "b")), "exactly one group")
  expect_error(contribution(fit), "`by` is required")

  # And the sink is not one of the group nodes the territories were built
  # from, however it is reported.
  ledger <- contribution(fit, by = c("a", "a", "b"))
  expect_identical(ledger$metadata$aggregation$measurements, 3L)
  expect_false(grepl("from 4 group node",
    paste(utils::capture.output(print(ledger)), collapse = " ")))
})

test_that("`using` is joined on `node`, never by row order", {
  fit <- pvw_fit()
  ordered <- contribution(fit, by = c("anterior", "anterior", "posterior"))

  # The same atlas with its rows permuted must give the same ledger: the join
  # is on the identifier, so row order carries no meaning.
  atlas <- data.frame(
    node = c("group3", "group1", "group2"),
    network = c("posterior", "anterior", "anterior"),
    stringsAsFactors = FALSE
  )
  joined <- contribution(fit, by = "network", using = atlas)
  expect_identical(joined$values, ordered$values)
  expect_identical(joined$index$node, ordered$index$node)

  # A table that cannot say which node each row is about is refused rather
  # than bound by position, which is how a permuted atlas becomes a wrong
  # ledger with nothing to see.
  expect_error(
    contribution(fit, by = "network",
      using = data.frame(network = c("anterior", "anterior", "posterior"))),
    "must name the group nodes"
  )
  expect_error(
    contribution(fit, by = c("a", "a", "b"), using = atlas),
    "only meaningful when `by` names one of its columns"
  )
})

test_that("`term` is named when the group model fits more than one", {
  data <- data.frame(subject = names(pvw_sizes),
    group = c("a", "a", "b", "b"), stringsAsFactors = FALSE)
  plan <- pvw_plan(model = ~ group, data = data)
  fit <- pvw_fit(plan)

  expect_error(contrast_energy(fit, c(1, -1, 0)), "`term` is required")
  by_name <- contrast_energy(fit, c(1, -1, 0), term = "groupb")
  by_position <- contrast_energy(fit, c(1, -1, 0), term = 2L)
  expect_identical(by_name$scientific_plan_id, by_position$scientific_plan_id)
  expect_equal(drop(by_name$values),
    fit$coefficients[, "face-house", "groupb"], ignore_attr = TRUE)
  expect_error(contrast_energy(fit, c(1, -1, 0), term = "age"),
    "names no column of the group model")

  # The identity is per term, so it has to hold for the slope too.
  ledger <- contribution(fit, by = c("a", "a", "b"), term = "groupb")
  expect_lt(max(abs(colSums(ledger$values) -
    pvw_global(plan, fit, "groupb"))), 1e-12)
})


# Identity -------------------------------------------------------------------

test_that("a view's identity is its parameters and nothing else", {
  plan <- pvw_plan()
  fit <- pvw_fit(plan)
  baseline <- contrast_energy(fit, c(1, -1, 0))

  # Same result, same arguments -- including the same contrast spelled with
  # names -- is the same view.
  expect_identical(baseline$scientific_plan_id,
    contrast_energy(fit, c(1, -1, 0))$scientific_plan_id)
  expect_identical(baseline$scientific_plan_id,
    contrast_energy(fit,
      c(face = 1, house = -1, tool = 0))$scientific_plan_id)
  expect_identical(baseline$scientific_plan_id,
    contrast_energy(pvw_fit(plan), c(1, -1, 0))$scientific_plan_id)

  # A contrast energy is quadratic in its contrast, so `-c` and `c` are one
  # view. The identity is over the packed operator for exactly this reason:
  # two spellings that produce the same numbers must not produce two ids.
  negated <- contrast_energy(fit, c(-1, 1, 0))
  expect_identical(negated$values, baseline$values)
  expect_identical(negated$scientific_plan_id, baseline$scientific_plan_id)

  # Different weights, a different verb, a different term and a different
  # ledger are four different views.
  expect_false(identical(baseline$scientific_plan_id,
    contrast_energy(fit, c(1, 0, -1))$scientific_plan_id))
  expect_false(identical(baseline$scientific_plan_id,
    rdm(fit, pairs = cbind("face", "house"))$scientific_plan_id))
  expect_false(identical(
    rdm(fit)$scientific_plan_id,
    rdm(fit, pairs = cbind("face", "house"))$scientific_plan_id))
  expect_false(identical(baseline$scientific_plan_id,
    contrast_energy(pvw_fit(plan, component = "coherent"),
      c(1, -1, 0))$scientific_plan_id))

  data <- data.frame(subject = names(pvw_sizes),
    group = c("a", "a", "b", "b"), stringsAsFactors = FALSE)
  covariate <- pvw_fit(pvw_plan(model = ~ group, data = data))
  expect_false(identical(
    contrast_energy(covariate, c(1, -1, 0), term = 1L)$scientific_plan_id,
    contrast_energy(covariate, c(1, -1, 0), term = 2L)$scientific_plan_id))

  # An aggregate is a view of a view, so the grouping moves the identity.
  expect_false(identical(
    contribution(baseline, by = c("a", "a", "b"))$scientific_plan_id,
    contribution(baseline, by = c("a", "b", "b"))$scientific_plan_id))
  expect_true(.strong_sha256(sub("^population-", "",
    contribution(baseline, by = c("a", "a", "b"))$scientific_plan_id)))
})


# Labelling ------------------------------------------------------------------

test_that("a transported component is named by its ledger, never bare", {
  fit <- pvw_fit(component = "coherent")
  view <- contrast_energy(fit, c(1, -1, 0))

  expect_identical(view$ledger, "native_coherent_ledger")
  expect_true(view$native_ledger)
  expect_identical(view$receipt$ledger, "native_coherent_ledger")

  bare <- c("coherent", "configuration")
  for (value in list(view, rdm(fit), contribution(fit, by = c("a", "a", "b")),
      contrast_energy(pvw_fit(component = "configuration"), c(1, -1, 0)))) {
    table <- as.data.frame(value)
    expect_false(any(names(table) %in% bare))
    expect_false(any(vapply(table, function(column) {
      any(as.character(column) %in% bare)
    }, logical(1))))
    expect_true("ledger" %in% names(table))
    expect_true(all(table$ledger %in% .population_ledger_names))
  }

  # The printed lines too. `native_coherent_ledger` contains the bare word as
  # a substring, so the assertion removes the ledger names first and then
  # looks for the word standing on its own.
  printed <- paste(utils::capture.output(print(view)), collapse = "\n")
  stripped <- gsub("native_(coherent|configuration)_ledger", "", printed)
  expect_false(grepl("\\bcoherent\\b", stripped))
  expect_false(grepl("\\bconfiguration\\b", stripped))
})

test_that("the print lines section 8.1 requires are all there", {
  fit <- pvw_fit(component = "coherent")
  printed <- utils::capture.output(print(contrast_energy(fit, c(1, -1, 0))))
  # The closing note is `strwrap()`ed, so a phrase may straddle two lines;
  # the assertions below are about the sentence, not about the line breaks.
  text <- gsub("\\s+", " ", paste(printed, collapse = " "))

  # The native frame family the ledger belongs to, and the transport identity
  # that carried it: semantics, provenance and cross-fit status.
  expect_true(any(grepl("^  frame:", printed)))
  expect_true(any(grepl("^  transport:", printed)))
  expect_true(grepl("budget", text, fixed = TRUE))
  expect_true(grepl("anatomical", text, fixed = TRUE))
  expect_true(grepl("cross-fit", text, fixed = TRUE))
  expect_true(grepl("native_coherent_ledger", text, fixed = TRUE))
  # The other half of the obligation: the additive identity does survive
  # transport, and only the common-mode reading fails.
  expect_true(grepl("transported_total holds exactly", text, fixed = TRUE))
  expect_true(grepl("not a group-node common mode", text, fixed = TRUE))
  # And the basis line, which says how the number was reached.
  expect_true(any(grepl("^  basis:", printed)))

  expect_true(grepl("effect_population_view",
    format(contrast_energy(fit, c(1, -1, 0))), fixed = TRUE))

  aggregate <- contribution(fit, by = c("a", "a", "b"))
  expect_true(any(grepl("^  aggregation:",
    utils::capture.output(print(aggregate)))))
})

test_that("the receipt carries transport semantics and the view identity", {
  fit <- pvw_fit(component = "coherent")
  view <- rdm(fit)
  receipt <- view$receipt

  expect_identical(receipt$transport$semantics, "budget")
  expect_identical(receipt$transport$provenance, "anatomical")
  expect_true(receipt$native_ledger)
  expect_identical(receipt$scientific_plan_id, view$scientific_plan_id)
  expect_identical(receipt$population_plan_id,
    fit$receipt$population_plan_id)
  expect_identical(receipt$population_result_id, fit$scientific_plan_id)
  expect_true(isTRUE(receipt$budget$asserted))
  expect_identical(receipt$normalization$mode, "none")
  expect_identical(receipt$basis$kind, "query_bank")
  expect_identical(receipt$basis$queries, rownames(fit$queries))

  density <- rdm(pvw_fit(pvw_plan(semantics = "density")))
  expect_identical(density$receipt$transport$semantics, "density")
  expect_false(isTRUE(density$receipt$budget$asserted))
})


# Dispatch registration ------------------------------------------------------
#
# Converting a plain function into a generic moves half of its behaviour into
# NAMESPACE, and an `S3method()` line that goes missing does not fail at load:
# the generic still resolves, and every call then dies with "no applicable
# method" at first use. That happened once on this branch, in the window
# between the generics landing in `R/views.R` and their registration landing
# in NAMESPACE, and what it broke was a concurrent workstream rather than a
# test. So the registration is pinned in both directions --- declared in the
# file, and live in the loaded namespace.

pvw_dispatch <- list(
  c("contrast_energy", "default"),
  c("contrast_energy", "effect_population_result"),
  c("contribution", "default"),
  c("contribution", "effect_population_result"),
  c("contribution", "effect_population_view"),
  c("rdm", "default"),
  c("rdm", "effect_population_result"),
  c("rsa", "default"),
  c("rsa", "effect_population_result"),
  c("as.data.frame", "effect_population_view"),
  c("format", "effect_population_view"),
  c("print", "effect_population_view")
)

test_that("every population and default method is registered for dispatch", {
  for (entry in pvw_dispatch) {
    method <- utils::getS3method(entry[[1L]], entry[[2L]], optional = TRUE,
      envir = asNamespace("crossform"))
    expect_true(is.function(method),
      info = paste0(entry[[1L]], ".", entry[[2L]]))
  }

  # Declared: the line is in NAMESPACE, so a roxygen run or a hand edit that
  # drops one fails here rather than in somebody else's script.
  namespace_file <- testthat::test_path("..", "..", "NAMESPACE")
  skip_if_not(file.exists(namespace_file), "NAMESPACE not available")
  declared <- sub("^S3method\\(\\s*([^,]+?)\\s*,\\s*(.+?)\\s*\\).*$", "\\1|\\2",
    grep("^S3method\\(", readLines(namespace_file, warn = FALSE), value = TRUE))
  expected <- vapply(pvw_dispatch, paste, character(1), collapse = "|")
  expect_identical(setdiff(expected, declared), character(0))
})


# The geometry-side generics ------------------------------------------------

test_that("dispatch leaves the per-participant verbs exactly as they were", {
  domain <- abstract_domain(8, id = "pop-views-geometry")
  labels <- rep(c("r1", "r2", "r3", "r4"), each = 2)
  run <- function(gain) rbind(
    face = gain * c(1, 0.2, 0, 0, 0.5, 0.1, 0, 0),
    house = gain * c(0, 1, 0.1, 0, 0, 0.6, 0.2, 0),
    tool = gain * c(0, 0.1, 1, 0.2, 0, 0.1, 0.7, 0)
  )
  relation <- relation(list(run1 = run(1), run2 = run(0.9)), domain = domain)
  plan <- plan_geometry(relation,
    compile_frame(regions(labels, "conservative"), domain),
    cross_partitions(relation, independence = "independent"))

  effect <- contrast_energy(plan, c(face = 1, house = -1, tool = 0))
  expect_s3_class(effect, "effect_contrast_view")
  expect_length(effect$total, 4L)

  expect_s3_class(rdm(plan), "effect_rdm_view")
  expect_s3_class(rdm(plan, component = "coherent"), "effect_rdm_view")
  expect_s3_class(rsa(plan, models = list(animacy = pvw_animacy())),
    "effect_rsa_view")

  ledger <- contribution(effect, by = c("early", "early", "late", "late"))
  expect_s3_class(ledger, "effect_contrast_view")
  expect_lt(abs(sum(ledger$total) - sum(effect$total)), 1e-12)

  # The refusals of the default methods still reach their default methods.
  expect_identical(
    catch_refusal(contrast_energy(plan, c(1, -1, 0),
      remove_univariate = TRUE))$capability,
    "nondestructive_decomposition")
  expect_identical(
    catch_refusal(rdm(plan, normalize = "correlation"))$namespace,
    "geometry_views")
  expect_identical(
    catch_refusal(contribution(rdm(plan), by = labels[1:4]))$capability,
    "additive_contribution")

  # A generic still refuses an input neither method understands, with the
  # default method's own message rather than a dispatch error.
  expect_error(contrast_energy(42, c(1, -1, 0)), "effect_geometry_plan")
  expect_error(rdm(42), "effect_geometry_plan")

  # An argument that exists to be refused must not be able to disappear into
  # the generic's `...`. A misspelt `normalize` that vanished would return a
  # number answering a different question than the one that was asked.
  expect_error(rdm(plan, normalise = "correlation"), "does not take")
  expect_error(
    contrast_energy(plan, c(1, -1, 0), remove_univarite = TRUE),
    "does not take")
  expect_error(
    rsa(plan, models = list(animacy = pvw_animacy()), intercpt = FALSE),
    "does not take")
  expect_error(contribution(effect, by = labels[1:4], usng = 1),
    "does not take")

  fit <- pvw_fit()
  expect_error(rdm(fit, pars = cbind("face", "house")), "does not take")
  expect_error(contrast_energy(fit, c(1, -1, 0), trm = 1L), "does not take")
  expect_error(rsa(fit, models = list(a = pvw_animacy()), nusance = NULL),
    "does not take")
  expect_error(
    contribution(contrast_energy(fit, c(1, -1, 0)), by = c("a", "a", "b"),
      term = "(Intercept)"),
    "does not take")
})
