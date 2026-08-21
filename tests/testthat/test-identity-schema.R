# The identity-schema consolidation is a visible migration, not a silent one.
#
# `fixtures/identity-schema.rds` records every evidence-task identity this
# package publishes, twice: `$before` was recorded on the tree that still
# carried the legacy `effect-form-v1` semantic alongside `evidence-pairing-v1`,
# and `$after` on the consolidated tree. The tests below prove three things
# together -- that the current tree reproduces `$after` byte for byte, that the
# set of ids that moved is exactly the bridged effect-form route, and that the
# routes which were already native to `evidence-pairing-v1` did not move at
# all. See `design/crossform-execution-design.md`, "Identity schema
# consolidation (B6)".

identity_schema_fixture <- function() {
  readRDS(test_path("fixtures", "identity-schema.rds"))
}

# The bridged effect-form route: an open experimental boundary over a
# bridge-closed neural boundary. This is the route whose identity used to be
# the legacy `.effect_task_semantic()` digest.
identity_schema_migrated_ids <- c(
  "effect_form_complete", "effect_form_complete_adapter",
  "effect_form_pair_query", "effect_form_pair_query_base",
  "effect_form_physical_query", "effect_form_reversed",
  "effect_task_plan_id", "geometry_plan_fixed_metric",
  "geometry_plan_learned_metric"
)

test_that("every recorded identity reproduces exactly on the current tree", {
  fixture <- identity_schema_fixture()
  expect_identical(fixture$schema_version, 1L)
  current <- identity_schema_ids()
  expect_identical(names(current), names(fixture$after))

  # The two geometry-plan ids digest a metric whose content is BLAS-computed,
  # and the package's numerical contract deliberately does not promise
  # bitwise equality across platforms
  # (`numerical_contract()$bitwise_across_platforms`). Those two are pinned
  # exactly on any machine that reproduces them and structurally elsewhere;
  # every other id hashes structure only and must reproduce everywhere.
  metric_bearing <- c("geometry_plan_fixed_metric",
    "geometry_plan_learned_metric")
  stable <- setdiff(names(current), metric_bearing)
  expect_identical(current[stable], fixture$after[stable])
  for (name in metric_bearing) {
    if (!identical(current[[name]], fixture$after[[name]])) {
      expect_match(current[[name]],
        "^geometry-sha256:[[:xdigit:]]{64}$")
      # The migration claim still holds off-platform: the current id is not
      # the retired pre-consolidation id either.
      expect_false(identical(current[[name]], fixture$before[[name]]))
    } else {
      expect_identical(current[[name]], fixture$after[[name]])
    }
  }
})

test_that("exactly the bridged effect-form route migrated identity", {
  fixture <- identity_schema_fixture()
  expect_identical(names(fixture$before), names(fixture$after))
  moved <- names(fixture$after)[fixture$after != fixture$before]
  expect_identical(sort(moved), sort(identity_schema_migrated_ids))

  # The two routes that were already native to `evidence-pairing-v1` keep the
  # ids they were recorded with: consolidating the schema retired a duplicate
  # naming rule, it did not rename the tasks that never used it.
  unmoved <- setdiff(names(fixture$after), identity_schema_migrated_ids)
  expect_identical(sort(unmoved),
    c("effect_form_neural_query", "measurement_form"))
  expect_identical(fixture$after[unmoved], fixture$before[unmoved])
})

test_that("one identity schema names every evidence task", {
  parts <- identity_schema_relations()
  native <- identity_schema_native_tasks(parts)
  bridged <- crossform:::.compile_effect_evidence_task(
    parts$self, cross_partitions(parts$self)
  )
  schemas <- vapply(
    list(bridged, native$measurement, native$query_closed),
    function(task) task$identity_schema, character(1)
  )

  expect_identical(unique(schemas), "evidence-pairing-v1")
  # The constructor no longer takes a schema or a legacy semantic, so there is
  # no second naming rule left to select.
  expect_false(any(
    c("identity_schema", "compatibility_semantic") %in%
      names(formals(crossform:::.new_evidence_task))
  ))
  expect_false(exists(".effect_task_semantic",
    envir = asNamespace("crossform"), inherits = FALSE))
})

test_that("the bridged route's identity is its evidence-pairing semantic", {
  parts <- identity_schema_relations()
  over <- cross_partitions(parts$self)
  task <- crossform:::.compile_effect_evidence_task(parts$self, over)
  expected <- crossform:::.evidence_task_general_semantic(
    task$left_relation_id, task$right_relation_id, task$spaces,
    task$ordered_partition_products, task$experimental_boundary,
    task$neural_boundary, task$stages, task$materialization
  )

  expect_identical(task$semantic, expected)
  expect_identical(task$task_id, crossform:::.effect_task_id(expected))
  expect_identical(task$semantic$transport, "evidence_pairing")
  expect_silent(crossform:::.validate_evidence_task(task))
})
