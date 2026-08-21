# Map-scale admission coverage.
#
# A recorded gate is evidence for a map-scale compute verb: an export that,
# given a compiled plan and neural values, produces a scientific result at
# many spatial measurements. Constructors, printers, inspectors, and adapters
# do not count. A verb is certified only when its row's artifact is in the
# promote table and still binds.
#
# This file is the coverage list. `promote-artifacts.R` derives the shipped
# set from it. `tests/testthat/test-admission-coverage.R` refuses drift.
# `benchmarks/README.md` is the prose index of the same table.
#
# Roles:
#   certified  counted verb or in-scope contraction; artifact must ship
#   covered    counted verb already certified by another row's artifact
#   refused    shipped receipt of a path that is not admitted
#   local      recorded evidence that must stay out of the tarball
#   gap        counted verb with no promoted gate yet

.admission_row <- function(role, verbs, artifact = NA_character_,
                           runner = NA_character_, reader = NA_character_,
                           notes = "") {
  data.frame(
    role = role,
    verbs = verbs,
    artifact = artifact,
    runner = runner,
    reader = reader,
    notes = notes,
    stringsAsFactors = FALSE
  )
}

.crossform_admission_out_of_scope <- function() {
  c(
    "compile_frame",
    "plan_geometry",
    "relation",
    "pairing",
    "sampling_capabilities",
    "neuroim2_searchlights",
    "neuroim2_volume_domain"
  )
}

.crossform_admission_coverage <- function() {
  rbind(
    .admission_row(
      "certified", "rdm, contrast_energy",
      "public-map-scale-gate.rds", "run-public-map-scale-gate.R",
      "test-public-map-scale.R",
      "implicit vs explicit identity at 576 centers"
    ),
    .admission_row(
      "certified", "rdm, rsa, contrast_energy",
      "query-first-scale-gate.rds", "run-query-first-scale.R",
      "test-query-first-scale.R",
      "query-first rdm, rsa, and bilinear contrast at q = 100"
    ),
    .admission_row(
      "certified", "native_pair_query_allocation",
      "native-pair-allocation.rds", "run-native-pair-allocation.R",
      "test-certification-artifacts.R",
      paste0(
        "internal Rcpp admission court; cumulative allocation against the ",
        "retained two-pass R oracle"
      )
    ),
    .admission_row(
      "covered", "evaluate_geometry",
      "query-first-scale-gate.rds", "run-query-first-scale.R",
      "test-query-first-scale.R",
      "query-first public entry; same artifact as rdm/rsa"
    ),
    .admission_row(
      "covered", "materialize_geometry",
      "public-map-scale-gate.rds", "run-public-map-scale-gate.R",
      "test-public-map-scale.R",
      "late materialize-then-project path; do not add a third gate"
    ),
    .admission_row(
      "covered", "materialize_geometry",
      "query-first-scale-gate.rds", "run-query-first-scale.R",
      "test-query-first-scale.R",
      "materialized comparator in the query-first court"
    ),
    .admission_row(
      "certified", "crossnobis",
      "crossnobis-scale-gate.rds", "run-crossnobis-scale-gate.R",
      "test-certification-artifacts.R",
      "learned crossnobis on the 52,416-feature fixture"
    ),
    .admission_row(
      "certified", "rdm_sampling_covariance",
      "sampling-covariance-scale.rds", "run-sampling-covariance-scale.R",
      "test-evidence-sampling-scale.R",
      "factorized RDM-variance operations at volume scale"
    ),
    .admission_row(
      "certified", "plan_relation, estimate_relation, fmrireg_relation",
      "first-moment-vertical-slice.rds", "run-first-moment-vertical-slice.R",
      "test-first-moment-vertical-slice.R",
      "BIDS-shaped facts through a study-bound relation plan"
    ),
    .admission_row(
      "certified", "memory_contraction",
      "small-dense-memory-cold.rds", "run-memory-benchmarks.R",
      "test-benchmark.R",
      "sequential sparse additive-frame workspace; not a public export"
    ),
    .admission_row(
      "certified", "memory_contraction",
      "medium-sparse-memory-cold.rds", "run-memory-benchmarks.R",
      "test-benchmark.R",
      "sequential sparse additive-frame workspace; not a public export"
    ),
    .admission_row(
      "certified", "memory_contraction",
      "medium-sparse-block-cold.rds", "run-memory-benchmarks.R",
      "test-benchmark.R",
      "sequential sparse additive-frame workspace; not a public export"
    ),
    .admission_row(
      "certified", "memory_contraction",
      "medium-sparse-memory-warm.rds", "run-memory-benchmarks.R",
      "test-benchmark.R",
      "sequential sparse additive-frame workspace; not a public export"
    ),
    .admission_row(
      "certified", "measurement_form",
      "measurement-profile.rds", "run-measurement-profile.R",
      "test-certification-artifacts.R",
      paste0(
        "scalar and requested-multivariate route court; BLAS-backed path ",
        "retained after the native-admission rule declined Rcpp"
      )
    ),
    .admission_row(
      "refused", "shard_executor",
      "shard-admission.rds", "run-shard-admission.R",
      "test-certification-artifacts.R",
      "parity and cleanup hold; cold speedup did not; artifact is unbound"
    ),
    .admission_row(
      "local", "rdm_sampling_covariance",
      "sampling-covariance-validation.rds",
      "run-sampling-covariance-validation.R",
      "test-certification-artifacts.R",
      "10,000-rep Monte Carlo; exceeds the 64 KiB shipped cap"
    ),
    .admission_row(
      "local", "crossnobis",
      "learned-metric-policy-validation.rds",
      "run-learned-metric-policy-validation.R",
      "test-certification-artifacts.R",
      "500-rep statistical recovery; exceeds the 64 KiB shipped cap"
    ),
    .admission_row(
      "local", "population_uncertainty",
      "population-null-coverage.rds",
      "run-population-null-coverage.R",
      "test-population-uncertainty.R",
      paste0("2,000-rep null coverage of the between-subject SE; the record ",
        "is small but it certifies a reader verb rather than a map-scale ",
        "one, so the committed receipt is the summary CSV")
    )
  )
}

.crossform_promotable_artifacts <- function(
    coverage = .crossform_admission_coverage()) {
  keep <- coverage$role %in% c("certified", "refused") &
    !is.na(coverage$artifact) & nzchar(coverage$artifact)
  stats::setNames(coverage$runner[keep], coverage$artifact[keep])
}
