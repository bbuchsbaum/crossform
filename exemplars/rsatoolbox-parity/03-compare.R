#!/usr/bin/env Rscript
# 03-compare.R -- the agreement table.
#
# Joins crossform's and rsatoolbox's outputs on (measurement, pair) and
# (measurement, term) and records max absolute and max relative differences
# against a declared tolerance. Nothing here is rounded before comparison and
# nothing is dropped: every row written by either arm must find a partner.

exemplar_dir <- if (nzchar(Sys.getenv("EXEMPLAR_DIR"))) {
  Sys.getenv("EXEMPLAR_DIR")
} else {
  normalizePath(dirname(sub("^--file=", "",
    grep("^--file=", commandArgs(FALSE), value = TRUE)[1])))
}
source(file.path(exemplar_dir, "00-common.R"))
results <- file.path(exemplar_dir, "results")

need <- c("crossform-rdm.csv", "crossform-rsa.csv", "rsatoolbox-rdm.csv",
          "rsatoolbox-rsa.csv", "rsatoolbox-meta.csv", "fixture-meta.csv",
          "crossform-distances.csv", "rsatoolbox-distances.csv",
          "crossform-correlation.csv",
          "crossform-similarity.csv", "rsatoolbox-similarity.csv")
missing <- need[!file.exists(file.path(results, need))]
if (length(missing)) {
  stop("Missing ", paste(missing, collapse = ", "),
       ". Run 01-fixture.R then 02-rsatoolbox.py first.")
}

read_results <- function(name) {
  utils::read.csv(file.path(results, name), stringsAsFactors = FALSE)
}

cf_rdm <- read_results("crossform-rdm.csv")
py_rdm <- read_results("rsatoolbox-rdm.csv")
cf_rsa <- read_results("crossform-rsa.csv")
py_rsa <- read_results("rsatoolbox-rsa.csv")
py_meta <- read_results("rsatoolbox-meta.csv")
fixture_meta <- read_results("fixture-meta.csv")

summarise <- function(a, b) {
  stopifnot(length(a) == length(b), length(a) > 0L)
  scale <- pmax(abs(a), abs(b))
  relative <- ifelse(scale > 0, abs(a - b) / scale, 0)
  list(n = length(a), max_abs = max(abs(a - b)), max_rel = max(relative))
}

rows <- list()
record <- function(quantity, comparator, n_values, max_abs, max_rel,
                   tolerance, note) {
  rows[[length(rows) + 1L]] <<- data.frame(
    quantity = quantity, comparator = comparator, n_values = n_values,
    max_abs_diff = max_abs, max_rel_diff = max_rel, tolerance = tolerance,
    passes = max_abs <= tolerance, note = note, stringsAsFactors = FALSE
  )
}

## ---- Crossnobis RDM -----------------------------------------------------
joined <- merge(cf_rdm, py_rdm,
                by = c("measurement", "pair", "left", "right"))
stopifnot(nrow(joined) == nrow(cf_rdm), nrow(joined) == nrow(py_rdm))
agreement <- summarise(joined$crossform, joined$rsatoolbox)
record("crossnobis_rdm", "rsatoolbox::calc_rdm_crossnobis",
       agreement$n, agreement$max_abs, agreement$max_rel, TOLERANCE,
       "fixed noise precision, cross-run, 4 measurements x 15 pairs")

oracle <- summarise(joined$crossform, joined$allpairs_oracle)
record("crossnobis_rdm", "explicit all-pairs numpy oracle",
       oracle$n, oracle$max_abs, oracle$max_rel, TOLERANCE,
       "third-party check that LOO folding equals uniform C(P,2) pairing")

# Per-measurement breakdown: an aggregate can hide one bad region.
for (m in sort(unique(joined$measurement))) {
  part <- joined[joined$measurement == m, ]
  s <- summarise(part$crossform, part$rsatoolbox)
  record(paste0("crossnobis_rdm[", m, "]"),
         "rsatoolbox::calc_rdm_crossnobis", s$n, s$max_abs, s$max_rel,
         TOLERANCE, "per-measurement breakdown")
}

## ---- Linear RSA ---------------------------------------------------------
lstsq <- py_rsa[py_rsa$route == "numpy_lstsq", ]
rsa_joined <- merge(cf_rsa, lstsq, by = c("measurement", "term"))
stopifnot(nrow(rsa_joined) == nrow(cf_rsa), nrow(rsa_joined) == nrow(lstsq))
s <- summarise(rsa_joined$crossform, rsa_joined$rsatoolbox)
record("linear_rsa_coefficients", "numpy least squares on vectorised RDMs",
       s$n, s$max_abs, s$max_rel, TOLERANCE,
       "crossform::rsa() is OLS in RDM space; same design, same response")

with_intercept <- rsa_joined[!grepl("^nointercept:", rsa_joined$term), ]
s <- summarise(with_intercept$crossform, with_intercept$rsatoolbox)
record("linear_rsa_coefficients[intercept]",
       "numpy least squares on vectorised RDMs", s$n, s$max_abs, s$max_rel,
       TOLERANCE, "(Intercept) + category + animacy")

regress <- py_rsa[py_rsa$route == "fit_regress_cosine_rescaled", ]
regress_joined <- merge(cf_rsa, regress, by = c("measurement", "term"))
stopifnot(nrow(regress_joined) == nrow(regress))
s <- summarise(regress_joined$crossform, regress_joined$rsatoolbox)
record("linear_rsa_coefficients[no intercept]",
       "rsatoolbox ModelWeighted + fit_regress(cosine), rescaled",
       s$n, s$max_abs, s$max_rel, 1e-8,
       paste0("fit_regress optimises a cosine-normalised objective: theta = ",
              "beta_OLS / sqrt(mean(d^2)); the scale factor is undone here"))

## ---- Non-crossvalidated distances ---------------------------------------
# crossform reaches these through a declared biased self-pairing; rsatoolbox
# through calc_rdm(method = ...). Correlation distance is read downstream from
# a guaranteed-PSD within-sample self form, which is the case the
# correlation-distance policy licenses; `rdm(normalize=)` stays refused.
cf_dist <- read_results("crossform-distances.csv")
py_dist <- read_results("rsatoolbox-distances.csv")

distance_notes <- c(
  euclidean = paste0("biased self-pairing + identity metric; both sides ",
                     "divide by the channel count"),
  mahalanobis = paste0("biased self-pairing + the same fixed noise ",
                       "precision restricted to each support"),
  correlation = paste0("downstream from a PSD within-sample self form: ",
                       "G_ij = (G_ii + G_jj - d_ij)/2, then ",
                       "1 - G_ij/sqrt(G_ii G_jj)")
)
cf_corr <- read_results("crossform-correlation.csv")
for (method in names(distance_notes)) {
  source_table <- if (method == "correlation") cf_corr else cf_dist
  left <- source_table[, c("measurement", "pair", "left", "right", method)]
  names(left)[names(left) == method] <- "crossform"
  right <- py_dist[py_dist$method == method,
                   c("measurement", "pair", "left", "right", "rsatoolbox")]
  joined_d <- merge(left, right, by = c("measurement", "pair", "left", "right"))
  stopifnot(nrow(joined_d) == nrow(left), nrow(joined_d) == nrow(right))
  s <- summarise(joined_d$crossform, joined_d$rsatoolbox)
  record(paste0(method, "_rdm"), paste0("rsatoolbox::calc_rdm(method=\"",
                                        method, "\")"),
         s$n, s$max_abs, s$max_rel, TOLERANCE, distance_notes[[method]])
}

## ---- Downstream comparison statistics -----------------------------------
# Eight rsatoolbox compare() methods, computed in base R from
# `rdm()$values`. No crossform export is involved; see 00-common.R.
cf_sim <- read_results("crossform-similarity.csv")
py_sim <- read_results("rsatoolbox-similarity.csv")
sim_joined <- merge(cf_sim, py_sim, by = c("measurement", "model", "method"))
stopifnot(nrow(sim_joined) == nrow(cf_sim), nrow(sim_joined) == nrow(py_sim))
s <- summarise(sim_joined$crossform, sim_joined$rsatoolbox)
record("rdm_comparison_statistics", "rsatoolbox::compare (all methods)",
       s$n, s$max_abs, s$max_rel, TOLERANCE,
       paste0("cosine, corr, spearman, kendall, tau-a, rho-a, cosine_cov, ",
              "corr_cov over 3 models x 4 measurements, downstream in base R"))

for (method in sort(unique(sim_joined$method))) {
  part <- sim_joined[sim_joined$method == method, ]
  s <- summarise(part$crossform, part$rsatoolbox)
  record(paste0("rdm_comparison[", method, "]"),
         paste0("rsatoolbox::compare(method=\"",
                sub("_", "-", sub("^(tau|rho)_a$", "\\1-a", method)), "\")"),
         s$n, s$max_abs, s$max_rel, TOLERANCE, "per-method breakdown")
}

# The tie-sensitive pairs must actually be distinguished by the fixture, or
# four of the eight rows above would be testing the same number twice.
tie_check <- function(a, b) {
  x <- sim_joined[sim_joined$method == a, ]
  y <- sim_joined[sim_joined$method == b, ]
  m <- merge(x, y, by = c("measurement", "model"))
  max(abs(m$rsatoolbox.x - m$rsatoolbox.y))
}
tau_separation <- tie_check("tau_a", "kendall")
rho_separation <- tie_check("rho_a", "spearman")
stopifnot(tau_separation > 1e-6, rho_separation > 1e-6)
comparison_meta <- c(
  tau_a_vs_kendall_separation = tau_separation,
  rho_a_vs_spearman_separation = rho_separation,
  n_similarity_values = nrow(sim_joined),
  n_models = length(unique(sim_joined$model)),
  n_methods = length(unique(sim_joined$method))
)
utils::write.csv(
  data.frame(key = names(comparison_meta), value = unname(comparison_meta),
             stringsAsFactors = FALSE),
  file.path(results, "comparison-meta.csv"), row.names = FALSE)
message("tie separation: tau-a vs kendall ",
        format(tau_separation, digits = 3), ", rho-a vs spearman ",
        format(rho_separation, digits = 3),
        " (both must be nonzero or the tied fixture is not doing its job)")

## ---- Noise estimators ---------------------------------------------------
env_all <- setNames(py_meta$value, py_meta$key)
record("noise_precision[diagonal]",
       "rsatoolbox::prec_from_residuals(method=\"diag\")", 1L,
       as.numeric(env_all[["prec_from_residuals_diag_max_abs_diff"]]), NA_real_,
       TOLERANCE, "diagonal precision estimator, same pooled residuals")
record("noise_precision[shrinkage_eye]",
       "rsatoolbox::cov_from_residuals(method=\"shrinkage_eye\")", 1L,
       as.numeric(env_all[["shrinkage_eye_reconstruction_max_abs_diff"]]),
       NA_real_, TOLERANCE,
       paste0("conditional parity: lambda = ",
              format(as.numeric(env_all[["shrinkage_eye_implied_lambda"]]),
                     digits = 6),
              " derived analytically by rsatoolbox; the identity-target ",
              "shrinkage formula reproduces it"))
record("noise_precision[shrinkage]",
       "rsatoolbox::cov_from_residuals(method=\"shrinkage_diag\")", 1L,
       as.numeric(env_all[["shrinkage_diag_reconstruction_max_abs_diff"]]),
       NA_real_, TOLERANCE,
       paste0("conditional parity: rsatoolbox derives lambda = ",
              format(as.numeric(env_all[["shrinkage_diag_implied_lambda"]]),
                     digits = 6),
              " analytically; crossform's (1-lambda)S + lambda diag(S) ",
              "reproduces it. crossform declares lambda rather than tuning ",
              "it on evaluation data -- a policy difference, not a gap"))

agreement_table <- do.call(rbind, rows)
utils::write.csv(agreement_table, file.path(results, "agreement.csv"),
                 row.names = FALSE)

## ---- Coverage ledger ----------------------------------------------------
# One disposition per rsatoolbox capability. This is the artifact that answers
# "where there is overlap and there should be parity, can we show it?" -- and,
# just as importantly, where there is deliberately no overlap and why.
#
#   parity                  a direct numerical agreement row
#   parity_conditional      agreement under a stated condition
#   reproducible_downstream reproduced from crossform outputs in base R,
#                           with no crossform export involved
#   refused_by_design       a written contract refuses it; capability named
#   absent_no_claim         not implemented, and not claimed
#   out_of_scope            outside what this exemplar undertakes to match
coverage <- rbind(
  data.frame(area = "calc_rdm", capability = "euclidean",
             disposition = "parity", evidence = "euclidean_rdm"),
  data.frame(area = "calc_rdm", capability = "mahalanobis",
             disposition = "parity", evidence = "mahalanobis_rdm"),
  data.frame(area = "calc_rdm", capability = "crossnobis",
             disposition = "parity", evidence = "crossnobis_rdm"),
  data.frame(area = "calc_rdm", capability = "correlation",
             disposition = "reproducible_downstream",
             evidence = "correlation_rdm"),
  data.frame(area = "calc_rdm", capability = "poisson",
             disposition = "out_of_scope",
             evidence = "not a bilinear or PSD-normalised distance"),
  data.frame(area = "calc_rdm", capability = "poisson_cv",
             disposition = "out_of_scope",
             evidence = "not a bilinear or PSD-normalised distance"),
  data.frame(area = "calc_rdm", capability = "calc_rdm_unbalanced",
             disposition = "out_of_scope",
             evidence = "balanced fixture; LOO-equals-uniform-pairing needs it"),
  data.frame(area = "calc_rdm", capability = "calc_rdm_movie",
             disposition = "out_of_scope", evidence = "no temporal datasets"),
  data.frame(area = "noise", capability = "prec_from_residuals(full)",
             disposition = "parity", evidence = "crossnobis_rdm"),
  data.frame(area = "noise", capability = "prec_from_residuals(diag)",
             disposition = "parity", evidence = "noise_precision[diagonal]"),
  data.frame(area = "noise", capability = "cov_from_residuals(shrinkage_diag)",
             disposition = "parity_conditional",
             evidence = "noise_precision[shrinkage]"),
  data.frame(area = "noise", capability = "cov_from_residuals(shrinkage_eye)",
             disposition = "parity_conditional",
             evidence = "noise_precision[shrinkage_eye]"),
  data.frame(area = "compare", capability = "cosine",
             disposition = "reproducible_downstream",
             evidence = "rdm_comparison[cosine]"),
  data.frame(area = "compare", capability = "corr",
             disposition = "reproducible_downstream",
             evidence = "rdm_comparison[corr]"),
  data.frame(area = "compare", capability = "spearman",
             disposition = "reproducible_downstream",
             evidence = "rdm_comparison[spearman]"),
  data.frame(area = "compare", capability = "kendall/tau-b",
             disposition = "reproducible_downstream",
             evidence = "rdm_comparison[kendall]"),
  data.frame(area = "compare", capability = "tau-a",
             disposition = "reproducible_downstream",
             evidence = "rdm_comparison[tau_a]"),
  data.frame(area = "compare", capability = "rho-a",
             disposition = "reproducible_downstream",
             evidence = "rdm_comparison[rho_a]"),
  data.frame(area = "compare", capability = "cosine_cov",
             disposition = "reproducible_downstream",
             evidence = "rdm_comparison[cosine_cov]"),
  data.frame(area = "compare", capability = "corr_cov",
             disposition = "reproducible_downstream",
             evidence = "rdm_comparison[corr_cov]"),
  data.frame(area = "compare", capability = "neg_riem_dist",
             disposition = "out_of_scope",
             evidence = "Riemannian metric on PSD forms; separate contract"),
  data.frame(area = "compare", capability = "bures / bures_metric",
             disposition = "out_of_scope",
             evidence = "Riemannian metric on PSD forms; separate contract"),
  data.frame(area = "model", capability = "ModelFixed",
             disposition = "parity", evidence = "linear_rsa_coefficients"),
  data.frame(area = "model", capability = "ModelWeighted + fit_regress",
             disposition = "parity",
             evidence = "linear_rsa_coefficients[no intercept]"),
  data.frame(area = "model", capability = "ModelSelect / ModelInterpolate",
             disposition = "absent_no_claim",
             evidence = "data-adaptive model selection; excluded from the ",
             stringsAsFactors = FALSE),
  data.frame(area = "rdm_view", capability = "rdm(normalize=)",
             disposition = "refused_by_design",
             evidence = "capability guaranteed_psd; correlation-distance policy"),
  data.frame(area = "inference", capability = "eval_fixed / eval_bootstrap*",
             disposition = "absent_no_claim",
             evidence = "crossform supplies no inference layer, and says so"),
  data.frame(area = "inference", capability = "noise ceilings",
             disposition = "absent_no_claim",
             evidence = "crossform supplies no inference layer, and says so"),
  stringsAsFactors = FALSE
)
coverage$evidence[coverage$capability == "ModelSelect / ModelInterpolate"] <-
  "data-adaptive; excluded by common-geometry-equivalence.md:156"

# A disposition that claims evidence must name an agreement row that exists
# and passes. Without this the ledger could drift into decoration.
claimed <- coverage$disposition %in%
  c("parity", "parity_conditional", "reproducible_downstream")
unknown <- setdiff(coverage$evidence[claimed], agreement_table$quantity)
if (length(unknown)) {
  stop("Coverage ledger cites agreement rows that do not exist: ",
       paste(unknown, collapse = ", "))
}
failing <- agreement_table$quantity[!agreement_table$passes]
if (any(coverage$evidence[claimed] %in% failing)) {
  stop("Coverage ledger cites a failing agreement row.")
}
utils::write.csv(coverage, file.path(results, "coverage.csv"),
                 row.names = FALSE)
message("Coverage ledger: ", nrow(coverage), " capabilities -- ",
        paste(sprintf("%s %d", names(table(coverage$disposition)),
                      as.integer(table(coverage$disposition))),
              collapse = ", "))

## ---- Report -------------------------------------------------------------
env <- setNames(py_meta$value, py_meta$key)
message("Environment: python ", env[["python"]], ", rsatoolbox ",
        env[["rsatoolbox"]], ", numpy ", env[["numpy"]])
message("Fixture: seed ", fixture_meta$value[fixture_meta$key == "seed"],
        ", ", fixture_meta$value[fixture_meta$key == "n_conditions"],
        " conditions x ",
        fixture_meta$value[fixture_meta$key == "n_runs"], " runs x ",
        fixture_meta$value[fixture_meta$key == "n_voxels"], " voxels; ",
        "residual covariance condition number ",
        format(as.numeric(fixture_meta$value[
          fixture_meta$key == "covariance_condition_number"]), digits = 4))
message("Noise precision: rsatoolbox prec_from_residuals(method='full', ",
        "dof=", fixture_meta$value[fixture_meta$key == "residual_df"],
        ") vs the R pooled estimate: max abs diff = ",
        format(as.numeric(env[["prec_from_residuals_max_abs_diff"]]),
               digits = 3))
message("")
print(agreement_table, row.names = FALSE, digits = 4)
message("")

failures <- agreement_table[!agreement_table$passes, ]
if (nrow(failures)) {
  print(failures, row.names = FALSE)
  stop("Parity failed for ", nrow(failures), " comparison(s).")
}
message("All ", nrow(agreement_table), " comparisons within tolerance. ",
        "Worst max abs diff = ",
        format(max(agreement_table$max_abs_diff), digits = 3))
