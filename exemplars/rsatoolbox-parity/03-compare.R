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
          "rsatoolbox-rsa.csv", "rsatoolbox-meta.csv", "fixture-meta.csv")
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

agreement_table <- do.call(rbind, rows)
utils::write.csv(agreement_table, file.path(results, "agreement.csv"),
                 row.names = FALSE)

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
