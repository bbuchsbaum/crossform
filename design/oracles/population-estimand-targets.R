# Independent finite-sample court for `population-estimand-v1`.
#
# Base R only: this file distinguishes the all-planned target from the
# available-at-node target, computes coverage and effective sample sizes, and
# makes the exact subject sets recoverable.  It does not call crossform.

tol <- 1e-12
subjects <- paste0("s", seq_len(5L))
z <- c(-2, -1, 0, 1, 2)
X <- cbind(`(Intercept)` = 1, z = z)
rownames(X) <- subjects

# A nonlinear response makes coverage-associated selection visibly change the
# best linear coefficient rather than merely changing its standard error.
complete_response <- cbind(
  full = 1 + 2 * z + 0.5 * z^2,
  left_missing = -0.5 + z + 0.75 * z^2,
  right_missing = 2 - z + 0.25 * z^2
)
rownames(complete_response) <- subjects

# Transported native territory at each ordinary group node.  Positive mass is
# operator coverage; an exact zero is absence even when budget semantics would
# encode the transported response itself as numeric zero.
transported_mass <- cbind(
  full = c(1.0, 0.8, 1.1, 0.9, 1.0),
  left_missing = c(1.0, 0.8, 0.6, 0.5, 0.0),
  right_missing = c(0.0, 0.4, 0.7, 0.9, 1.0)
)
rownames(transported_mass) <- subjects
available <- transported_mass > 0
response <- complete_response
response[!available] <- NA_real_

fit_cell <- function(y, admitted) {
  ids <- subjects[admitted]
  design <- X[admitted, , drop = FALSE]
  rank <- qr(design)$rank
  estimable <- rank == ncol(design)
  coefficients <- rep(NA_real_, ncol(design))
  names(coefficients) <- colnames(design)
  if (estimable) coefficients <- qr.coef(qr(design), y[admitted])
  list(
    subjects = ids,
    n = length(ids),
    rank = rank,
    residual_df = length(ids) - rank,
    estimable = rep(estimable, ncol(design)),
    coefficients = coefficients
  )
}

all_planned <- lapply(seq_len(ncol(response)), function(j) {
  admitted <- available[, j]
  if (!all(admitted)) {
    return(list(
      subjects = subjects, n = length(subjects), rank = qr(X)$rank,
      residual_df = length(subjects) - qr(X)$rank,
      estimable = rep(FALSE, ncol(X)),
      coefficients = stats::setNames(rep(NA_real_, ncol(X)), colnames(X)),
      reason = "planned_subject_unavailable"
    ))
  }
  c(fit_cell(response[, j], rep(TRUE, length(subjects))),
    list(reason = NA_character_))
})
names(all_planned) <- colnames(response)

available_at_node <- lapply(seq_len(ncol(response)), function(j) {
  c(fit_cell(response[, j], available[, j]), list(reason = NA_character_))
})
names(available_at_node) <- colnames(response)

coverage <- data.frame(
  node = colnames(response),
  n = colSums(available),
  fraction = colMeans(available),
  regression_n_eff = colSums(available),
  mass_n_eff = vapply(seq_len(ncol(transported_mass)), function(j) {
    mass <- transported_mass[available[, j], j]
    sum(mass)^2 / sum(mass^2)
  }, numeric(1)),
  rank = vapply(available_at_node, `[[`, integer(1), "rank"),
  residual_df = vapply(available_at_node, `[[`, integer(1), "residual_df")
)

coverage_association <- vapply(seq_len(ncol(available)), function(j) {
  if (all(available[, j]) || !any(available[, j])) return(NA_real_)
  mean(z[available[, j]]) - mean(z[!available[, j]])
}, numeric(1))
names(coverage_association) <- colnames(available)

# A compact result may dictionary-encode sets, but the exact identifiers—not
# just counts or hashes—remain recoverable.
set_key <- vapply(available_at_node, function(value) {
  paste(value$subjects, collapse = "\r")
}, character(1))
unique_key <- unique(set_key)
subject_set_dictionary <- stats::setNames(lapply(unique_key, function(key) {
  if (!nzchar(key)) character() else strsplit(key, "\r", fixed = TRUE)[[1L]]
}), paste0("set", seq_along(unique_key)))
subject_set_id <- stats::setNames(
  names(subject_set_dictionary)[match(set_key, unique_key)], names(set_key)
)

transport_signature <- stats::setNames(
  paste0("transport-fixture-", subjects), subjects
)
sink_territory <- stats::setNames(c(0.00, 0.05, 0.10, 0.20, 0.35), subjects)

stopifnot(
  all(is.finite(all_planned$full$coefficients)),
  all(is.na(all_planned$left_missing$coefficients)),
  all(is.na(all_planned$right_missing$coefficients)),
  all(vapply(available_at_node, function(value) {
    all(is.finite(value$coefficients))
  }, logical(1))),
  identical(available_at_node$left_missing$subjects, subjects[1:4]),
  identical(available_at_node$right_missing$subjects, subjects[2:5]),
  coverage$n[[1L]] == 5L,
  coverage$n[[2L]] == 4L,
  coverage$n[[3L]] == 4L,
  all(coverage$mass_n_eff <= coverage$n + tol),
  coverage_association[["left_missing"]] < 0,
  coverage_association[["right_missing"]] > 0,
  identical(subject_set_dictionary[[subject_set_id[["left_missing"]]]],
    subjects[1:4])
)

population_estimand_oracle <- list(
  subjects = subjects,
  design = X,
  complete_response = complete_response,
  response = response,
  transported_mass = transported_mass,
  available = available,
  all_planned = all_planned,
  available_at_node = available_at_node,
  coverage = coverage,
  coverage_association = coverage_association,
  subject_set_id = subject_set_id,
  subject_set_dictionary = subject_set_dictionary,
  transport_signature = transport_signature,
  sink_territory = sink_territory
)

message(sprintf(
  "population-estimand-v1 PASS: planned n %d, node coverage %s",
  length(subjects), paste(coverage$n, collapse = "/")
))
