prospective_readiness <- function(
    checks, phase = c("discovery", "replication"),
    execution = NULL, discovery_artifact_root = NULL) {
  phase <- match.arg(phase)
  required <- c("eligible_data", "frozen_configuration",
    "environment_locked", "provenance_bound", "analyst_identified",
    "artifact_storage_reserved")
  if (!is.list(checks) || !identical(sort(names(checks)), sort(required)) ||
      !all(vapply(checks, is.logical, logical(1))) ||
      !all(lengths(checks) == 1L) || anyNA(unlist(checks))) {
    stop("checks must contain one logical value for every readiness item.",
         call. = FALSE)
  }
  blocked <- required[!unlist(checks[required])]
  if (length(blocked)) {
    return(list(state = "BLOCKED", phase = phase, reasons = blocked,
      ready = FALSE, executed = FALSE,
      evidence_state = "prospective_protocol"))
  }
  if (is.null(execution)) {
    return(list(state = "READY", phase = phase, reasons = character(),
      ready = TRUE, executed = FALSE,
      evidence_state = "prospective_protocol"))
  }
  execution_required <- c("manifest_verified", "protocol_hash_predates_outcomes",
                          "deviations_logged", "artifact_root")
  if (!is.list(execution) ||
      !all(execution_required %in% names(execution)) ||
      !all(vapply(execution[c("manifest_verified",
                             "protocol_hash_predates_outcomes",
                             "deviations_logged")], isTRUE, logical(1))) ||
      !is.character(execution$artifact_root) ||
      length(execution$artifact_root) != 1L ||
      !nzchar(execution$artifact_root)) {
    return(list(state = "BLOCKED", phase = phase,
      reasons = "execution_evidence_incomplete", ready = FALSE,
      executed = FALSE, evidence_state = "prospective_protocol"))
  }
  if (phase == "replication" &&
      (is.null(discovery_artifact_root) ||
       identical(execution$artifact_root, discovery_artifact_root))) {
    return(list(state = "BLOCKED", phase = phase,
      reasons = "replication_artifact_must_be_distinct", ready = FALSE,
      executed = FALSE, evidence_state = "prospective_protocol"))
  }
  list(state = "EXECUTED", phase = phase, reasons = character(),
    ready = TRUE, executed = TRUE,
    evidence_state = if (phase == "discovery")
      "completed_real_data_result" else "independent_replication",
    artifact_root = execution$artifact_root)
}
