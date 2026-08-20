# The exported surface of crossform is a deliberate list, not a residue of
# where `@export` tags happened to land. `design/api-tiers.md` assigns every
# name below a tier (who it is for) and a disposition; this file is the gate
# that makes the ledger binding.
#
# Adding or removing an export is therefore a two-step edit, in this order:
#   1. record the decision in `design/api-tiers.md` (tier, user evidence,
#      disposition), then
#   2. edit `crossform_public_api` below to match.
# A diff that touches only NAMESPACE will fail here, which is the point: the
# failure is a request for the justification, not for a rubber stamp.
crossform_public_api <- c(
  "abstract_domain", "additive_frame", "aggregate_first",
  "anatomical_transport", "as_neurovol", "bids_study", "bilinear_query",
  "canonical_coupling", "catch_refusal",
  "coefficient_parameterization", "coherence_spectrum", "compile_frame",
  "compiler_conformance", "compute_policy", "condition_space",
  "connectivity", "contrast_energy", "contribution",
  "control_coupling", "coupling",
  "coupling_contrast", "covariance_coupling", "cross_partitions",
  "crossnobis", "design_model", "diagonal_precision", "edge_frame",
  "effect_coupling", "effect_extractor", "effect_map",
  "effect_space", "estimate_relation", "evaluate_geometry",
  "example_fmri_effects", "external_transport", "file_matrix_source",
  "fmridesign_design_model", "fmrireg_relation",
  "frame_conservation", "frame_family",
  "gaussian_covariance_model",
  "geometry_alignment", "geometry_component", "geometry_spectrum",
  "identity_metric", "lm_extractor", "lm_relation_fit",
  "location_transport", "lower_effect_map", "match_control", "match_coupling",
  "materialize_geometry", "measurement_components",
  "measurement_form", "measurement_frame", "metric_capabilities",
  "metric_training_policy", "neural_metric", "neuroim2_searchlights",
  "neuroim2_volume_domain", "noise_precision", "numerical_agreement",
  "numerical_contract", "observation_confounds",
  "observation_events", "observation_index", "observation_model",
  "observations", "pair_lm_query", "pair_query", "pairing",
  "partition_hierarchy", "plan_crossnobis", "plan_geometry",
  "plan_population", "plan_relation", "query_geometry", "raw_design_model",
  "raw_effect_map", "rdm", "rdm_sampling_covariance",
  "reconstruct_evidence", "reduce_partitions", "regions", "relation",
  "relation_block", "relation_fit", "relation_fit_capabilities",
  "relation_plan_receipts", "residual_block", "residual_df", "rsa",
  "sampling_capabilities", "sampling_covariance", "searchlights",
  "shrinkage_precision", "source_capabilities", "study",
  "study_axis", "study_capabilities", "transport_values", "variation_query",
  "volume_domain", "voxelwise", "whole_brain"
)

test_that("the exported surface is exactly the ledgered set", {
  exports <- sort(getNamespaceExports("crossform"))
  expected <- sort(crossform_public_api)

  # Report the drift in both directions before the identity check, so a
  # failure names the offending functions rather than printing two lists.
  expect_identical(setdiff(exports, expected), character(0))
  expect_identical(setdiff(expected, exports), character(0))
  expect_identical(exports, expected)

  # The count is stated separately because the ledger quotes it: the
  # subtraction release took 105 exports to 97, and WS-D has since added
  # `frame_family()`, `contribution()` and `coherence_spectrum()`, and WS-E
  # the four population transport names plus `plan_population()`
  # (`design/api-tiers.md`, "Additions after the subtraction release").
  expect_identical(length(exports), 105L)
})

test_that("NAMESPACE and the loaded namespace agree", {
  # Only meaningful from a source tree; an installed package is checked
  # through `getNamespaceExports()` above, which reads the same declarations.
  namespace_file <- testthat::test_path("..", "..", "NAMESPACE")
  skip_if_not(file.exists(namespace_file), "NAMESPACE not available")

  lines <- readLines(namespace_file, warn = FALSE)
  declared <- sort(sub(
    "^export\\(([^)]+)\\).*$", "\\1", grep("^export\\(", lines, value = TRUE)
  ))
  expect_identical(declared, sort(crossform_public_api))

  # NAMESPACE is hand-maintained, so a stray `exportPattern()` would widen the
  # surface without adding a single `export()` line.
  expect_identical(grep("^exportPattern", lines, value = TRUE), character(0))
})

test_that("the sanctioned developer entry points stay exported", {
  # `design/api-tiers.md`, developer tier: the extension-only surface an
  # adapter package needs and an end user never meets. It is five names, and
  # the ledger's standing rule is that it does not grow silently.
  developer_api <- c(
    "file_matrix_source", "source_capabilities", "effect_extractor",
    "relation_fit", "relation_block"
  )
  expect_true(all(developer_api %in% getNamespaceExports("crossform")))
  expect_true(all(developer_api %in% crossform_public_api))
})

test_that("the subtraction release's demotions stay demoted", {
  # Eight exports went internal in the subtraction release (tickets A3/A4).
  # Each remains a working function reachable with `crossform:::`; what was
  # removed is the promise that it is a stable public entry point.
  demoted <- c(
    "inner_product", "measurement_space", "measurement_bridge",
    "reverse_bridge", "effect_covariance", "residual_pair_statistics",
    "bids_events", "bids_confounds"
  )
  exports <- getNamespaceExports("crossform")

  expect_identical(intersect(demoted, exports), character(0))
  expect_identical(intersect(demoted, crossform_public_api), character(0))

  for (name in demoted) {
    expect_true(
      exists(name, envir = asNamespace("crossform"), inherits = FALSE),
      info = name
    )
    expect_true(
      is.function(get(name, envir = asNamespace("crossform"))),
      info = name
    )
  }
})
