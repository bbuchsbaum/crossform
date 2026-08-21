# Representative tasks for the identity-schema fixture.
#
# `test-identity-schema.R` compares the identities these cases produce against
# `fixtures/identity-schema.rds`, which records the ids on both sides of the
# B6 consolidation (one identity schema, `evidence-pairing-v1`). The cases are
# deliberately deterministic and self-contained: every relation is built from
# literal matrices, so an id only moves when the identity schema moves.

identity_schema_relations <- function() {
  domain <- abstract_domain(3, id = "identity-schema:neural:v1")
  self_space <- effect_space(
    c("a", "b"),
    basis_id = "identity-schema:conditions:v1"
  )
  right_space <- effect_space(
    c("x", "y", "z"),
    basis_id = "identity-schema:retrieve:v1"
  )
  self <- relation(
    list(
      r1 = matrix(1, 2, 3), r2 = matrix(2, 2, 3), r3 = matrix(3, 2, 3)
    ),
    effects = self_space, domain = domain
  )
  left <- relation(
    list(e1 = matrix(1, 2, 3), e2 = matrix(2, 2, 3)),
    effects = self_space, domain = domain
  )
  right <- relation(
    list(s1 = matrix(1, 3, 3), s2 = matrix(2, 3, 3)),
    effects = right_space, domain = domain
  )
  list(
    domain = domain, self_space = self_space, right_space = right_space,
    self = self, left = left, right = right
  )
}

# The two routes that never carried the legacy `effect-form-v1` semantic: a
# measurement form (closed experimental boundary, open neural boundary) and an
# effect form whose neural boundary closes with a fixed query rather than a
# bridge. Their ids must not move.
identity_schema_native_tasks <- function(parts = identity_schema_relations()) {
  measurement_edges <- crossform:::.ordered_partition_edges(
    cross_partitions(parts$self), parts$self$partitions,
    parts$self$partitions, TRUE
  )
  frame <- crossform:::.measurement_frame_from_additive(
    additive_frame(diag(3), domain = parts$domain)
  )
  measurement <- crossform:::.new_evidence_task(
    parts$self, parts$self, TRUE, measurement_edges,
    crossform:::.closed_experimental_boundary(
      pair_query(diag(2), parts$self_space, parts$self_space),
      role = "variation", sampling_axis = "trial"
    ),
    crossform:::.open_neural_boundary(
      frame, frame,
      crossform:::.measurement_edges(
        c("measurement_1", "measurement_2"),
        c("measurement_2", "measurement_3"), frame
      )
    ),
    crossform:::.evidence_stage_plan(),
    crossform:::.evidence_materialization(
      "measurement_form", "complete_form"
    )
  )

  rectangular_edges <- crossform:::.ordered_partition_edges(
    pairing(c("e1", "e2"), c("s1", "s2"), c(0.5, 0.5), directed = TRUE),
    parts$left$partitions, parts$right$partitions, FALSE
  )
  neural_query <- crossform:::.neural_pair_query(
    matrix(c(1, 0, -1, 2, 3, 1, 0, 1, 2), 3, 3),
    parts$domain, parts$domain,
    provenance = list(bridge = "identity-schema-fixture")
  )
  query_closed <- crossform:::.new_evidence_task(
    parts$left, parts$right, FALSE, rectangular_edges,
    crossform:::.open_experimental_boundary(
      parts$self_space, parts$right_space
    ),
    crossform:::.closed_neural_query_boundary(neural_query),
    crossform:::.evidence_stage_plan(),
    crossform:::.evidence_materialization("effect_form", "complete_form")
  )

  list(measurement = measurement, query_closed = query_closed)
}

# Every identity this package records that depends on an evidence task, keyed
# by the route that produces it.
identity_schema_ids <- function() {
  parts <- identity_schema_relations()
  native <- identity_schema_native_tasks(parts)
  self_over <- cross_partitions(parts$self)
  rectangular_over <- pairing(
    c("e1", "e2"), c("s1", "s2"), c(0.5, 0.5), directed = TRUE
  )
  bridged_complete <- crossform:::.compile_effect_evidence_task(
    parts$self, self_over
  )
  bridged_pair_query <- crossform:::.compile_effect_evidence_task(
    parts$left, rectangular_over, parts$right,
    pair_query(
      matrix(as.numeric(1:6), 2, 3), parts$self_space, parts$right_space
    )
  )
  bridged_physical_query <- crossform:::.compile_effect_evidence_task(
    parts$self, self_over, query = matrix(1, 3, 2)
  )
  adapter <- crossform:::.compile_effect_task(parts$self, self_over)

  geometry <- identity_schema_geometry_ids()

  c(
    effect_form_complete = bridged_complete$task_id,
    effect_form_complete_adapter = adapter$task_id,
    effect_form_pair_query = bridged_pair_query$task_id,
    effect_form_pair_query_base =
      crossform:::.effect_task_base_id(bridged_pair_query),
    effect_form_physical_query = bridged_physical_query$task_id,
    effect_form_reversed =
      crossform:::.reverse_evidence_task(bridged_pair_query)$task_id,
    effect_task_plan_id = crossform:::.effect_task_plan_id(
      parts$self, additive_frame(diag(3), domain = parts$domain),
      self_over, "full_geometry", NULL, "full"
    ),
    measurement_form = native$measurement$task_id,
    effect_form_neural_query = native$query_closed$task_id,
    geometry_plan_fixed_metric = geometry[["fixed"]],
    geometry_plan_learned_metric = geometry[["learned"]]
  )
}

# The two plan-level identities, kept separate because they need the residual
# fixture (helper-residual-statistics.R) rather than literal matrices.
identity_schema_geometry_ids <- function() {
  setup <- metric_learning_setup()
  fit <- setup$fixture$fit
  fixed <- plan_geometry(fit$relation, setup$fixture$frame, setup$over)
  # A learned recipe needs the fit itself, for its residual channel.
  learned <- plan_geometry(
    fit, setup$fixture$frame, setup$over,
    metric = shrinkage_precision(0.2),
    residual_workspace_bytes = setup$budgets$wider
  )
  c(
    fixed = fixed$scientific_plan_id,
    learned = learned$scientific_plan_id
  )
}
