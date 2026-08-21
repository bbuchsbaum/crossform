#!/usr/bin/env Rscript
# 04-extension.R -- the strict extension.
#
# Everything below comes out of the SAME relation fit, the SAME fixed noise
# precision, and the SAME cross-run pairing that 01-fixture.R handed to
# rsatoolbox in 02. Nothing is refitted. The point of the arm is that the RDM
# rsatoolbox returns is one linear functional of an object that also answers
# other questions, rather than the analysis product itself.
#
# Recorded here:
#   A. signed contrast energy (rsatoolbox's RDM entries are squared and
#      unsigned by construction)
#   B. the exact coherent / configuration / total partition of that energy
#   C. the crossnobis RDM re-derived pair by pair as `crossnobis()` on the
#      difference contrasts, showing the RDM is one query among many
#   D. the linear RSA coefficient compiled as one query over the same plan
#   E. analytic sampling covariance of the RDM under the admitted
#      fixed-metric law, and its exact transport to the RSA coefficient
#   F. the refusals that bound E

exemplar_dir <- if (nzchar(Sys.getenv("EXEMPLAR_DIR"))) {
  Sys.getenv("EXEMPLAR_DIR")
} else {
  normalizePath(dirname(sub("^--file=", "",
    grep("^--file=", commandArgs(FALSE), value = TRUE)[1])))
}
source(file.path(exemplar_dir, "00-common.R"))
invisible(load_crossform(exemplar_dir))
results <- file.path(exemplar_dir, "results")

state_file <- file.path(results, "fixture.rds")
if (!file.exists(state_file)) stop("Run 01-fixture.R first.")
state <- readRDS(state_file)
fit <- state$fit
plans <- state$plans
models <- state$models
pairs <- state$pairs
conditions <- CONDITIONS
n_pairs <- nrow(pairs)

records <- list()
add <- function(...) records[[length(records) + 1L]] <<-
  data.frame(..., stringsAsFactors = FALSE)

## ---- A/B. Signed contrast energy and its exact partition ----------------
message("A/B. contrast_energy(face - house) on the same plan")
energy <- lapply(plans, function(plan) contrast_energy(plan, CONTRAST))

for (nm in names(energy)) {
  e <- energy[[nm]]
  measurements <- as.character(e$index)
  # crossnobis() is the named Mahalanobis reading of exactly $total.
  named <- crossnobis(plans[[nm]], CONTRAST)
  stopifnot(identical(named$estimand,
                      "crossvalidated_squared_mahalanobis_contrast"))
  same_estimand <- max(abs(named$values - e$total))
  recomposition <- max(abs(e$coherent + e$configuration - e$total))
  message("  ", nm, ": crossnobis() vs contrast_energy()$total = ",
          format(same_estimand, digits = 3),
          "; coherent + configuration - total = ",
          format(recomposition, digits = 3))
  stopifnot(same_estimand < 1e-12, recomposition < 1e-12)
  for (i in seq_along(measurements)) {
    add(frame = nm, measurement = measurements[i], quantity = "signed",
        value = e$signed[[i]])
    add(frame = nm, measurement = measurements[i], quantity = "coherent",
        value = e$coherent[[i]])
    add(frame = nm, measurement = measurements[i], quantity = "configuration",
        value = e$configuration[[i]])
    add(frame = nm, measurement = measurements[i], quantity = "total",
        value = e$total[[i]])
    add(frame = nm, measurement = measurements[i],
        quantity = "coherence_fraction",
        value = if (isTRUE(e$coherence_fraction_valid[[i]])) {
          e$coherence_fraction[[i]]
        } else {
          NA_real_
        })
    add(frame = nm, measurement = measurements[i], quantity = "crossnobis",
        value = named$values[[i]])
  }
}

## ---- C. The RDM is one linear functional of the same evidence ----------
message("C. rdm() entries re-derived as crossnobis() of difference contrasts")
rdm_total <- lapply(plans, function(plan) as.matrix(rdm(plan)$values))
rdm_coherent <- lapply(plans, function(plan) {
  as.matrix(rdm(plan, component = "coherent")$values)
})
rdm_configuration <- lapply(plans, function(plan) {
  as.matrix(rdm(plan, component = "configuration")$values)
})

difference_gap <- 0
for (nm in names(plans)) {
  redone <- vapply(seq_len(n_pairs), function(k) {
    weights <- setNames(rep(0, length(conditions)), conditions)
    weights[[pairs$left[[k]]]] <- 1
    weights[[pairs$right[[k]]]] <- -1
    crossnobis(plans[[nm]], weights)$values
  }, numeric(nrow(rdm_total[[nm]])))
  if (is.null(dim(redone))) redone <- matrix(redone, 1L, n_pairs)
  difference_gap <- max(difference_gap, max(abs(redone - rdm_total[[nm]])))
  split_gap <- max(abs(
    rdm_coherent[[nm]] + rdm_configuration[[nm]] - rdm_total[[nm]]
  ))
  message("  ", nm, ": rdm() vs crossnobis(e_i - e_j) = ",
          format(max(abs(redone - rdm_total[[nm]])), digits = 3),
          "; coherent + configuration - total = ",
          format(split_gap, digits = 3))
  stopifnot(split_gap < 1e-12)
}
stopifnot(difference_gap < 1e-12)

for (nm in names(plans)) {
  measurements <- as.character(rdm(plans[[nm]])$index)
  for (i in seq_along(measurements)) {
    for (k in seq_len(n_pairs)) {
      label <- paste0(pairs$left[[k]], "|", pairs$right[[k]])
      add(frame = nm, measurement = measurements[i],
          quantity = paste0("rdm_total[", label, "]"),
          value = rdm_total[[nm]][i, k])
      add(frame = nm, measurement = measurements[i],
          quantity = paste0("rdm_coherent[", label, "]"),
          value = rdm_coherent[[nm]][i, k])
      add(frame = nm, measurement = measurements[i],
          quantity = paste0("rdm_configuration[", label, "]"),
          value = rdm_configuration[[nm]][i, k])
    }
  }
}

## ---- D. RSA coefficient, compiled over the same plan -------------------
rsa_fits <- lapply(plans, function(plan) rsa(plan, models = models))
for (nm in names(rsa_fits)) {
  coefficients <- as.matrix(rsa_fits[[nm]]$coefficients)
  measurements <- as.character(rsa_fits[[nm]]$index)
  for (i in seq_along(measurements)) {
    for (term in colnames(coefficients)) {
      add(frame = nm, measurement = measurements[i],
          quantity = paste0("rsa[", term, "]"), value = coefficients[i, term])
    }
  }
}

## ---- E. Analytic sampling covariance under the admitted law ------------
message("E. rdm_sampling_covariance() under the fixed-metric law")
design <- cbind(`(Intercept)` = 1,
                category = rdm_pair_vector(models$category),
                animacy = rdm_pair_vector(models$animacy))
transport <- solve(crossprod(design), t(design))["category", , drop = FALSE]

uncertainty <- list()
for (nm in names(plans)) {
  plan <- plans[[nm]]
  measurements <- as.character(rdm(plan)$index)
  for (target in c("null", "plugin")) {
    for (i in seq_along(measurements)) {
      cv <- rdm_sampling_covariance(plan, fit, target = target, at = i)
      se <- sqrt(sampling_covariance(cv))
      rsa_var <- as.numeric(sampling_covariance(
        cv, operation = "transport", query = transport
      ))
      uncertainty[[length(uncertainty) + 1L]] <- data.frame(
        frame = nm, measurement = measurements[i], target = target,
        pair = seq_len(n_pairs), left = pairs$left, right = pairs$right,
        rdm = rdm_total[[nm]][i, ], se = se,
        z = rdm_total[[nm]][i, ] / se,
        rsa_category = as.matrix(rsa_fits[[nm]]$coefficients)[i, "category"],
        rsa_category_se = sqrt(rsa_var),
        residual_df = cv$source$residual_df,
        residual_effective_dimension =
          cv$source$residual_effective_dimension,
        stringsAsFactors = FALSE
      )
      if (identical(target, "plugin") && identical(nm, "whole")) {
        add(frame = nm, measurement = measurements[i],
            quantity = "rsa[category] SE (plugin)", value = sqrt(rsa_var))
        for (k in seq_len(n_pairs)) {
          add(frame = nm, measurement = measurements[i],
              quantity = paste0("rdm SE (plugin)[", pairs$left[[k]], "|",
                                pairs$right[[k]], "]"),
              value = se[[k]])
        }
      }
    }
  }
}
uncertainty <- do.call(rbind, uncertainty)
message("  ", nrow(uncertainty), " (measurement, target, pair) rows; ",
        "median RDM |z| (plugin) = ",
        format(median(abs(uncertainty$z[uncertainty$target == "plugin"])),
               digits = 3))
message("  RSA[category] whole-brain estimate ",
        format(as.matrix(rsa_fits$whole$coefficients)[1, "category"],
               digits = 4), " +/- ",
        format(unique(uncertainty$rsa_category_se[
          uncertainty$frame == "whole" & uncertainty$target == "plugin"]),
          digits = 3), " (plugin)")

## ---- F. The refusals that bound E --------------------------------------
capture_refusal <- function(label, expr) {
  message <- tryCatch({ force(expr); NA_character_ },
                      error = function(e) conditionMessage(e))
  data.frame(label = label, refused = !is.na(message),
             message = if (is.na(message)) "" else message,
             stringsAsFactors = FALSE)
}

refusals <- rbind(
  capture_refusal(
    "sampling covariance target is never inferred",
    rdm_sampling_covariance(plans$whole, fit, at = 1L)),
  capture_refusal(
    "correlation-style normalization of the RDM",
    rdm(plans$whole, normalize = "correlation")),
  capture_refusal(
    "crossnobis without a declared noise-precision metric",
    crossnobis(
      plan_geometry(fit$relation, at = state$frames$whole, over = state$over),
      CONTRAST))
)
message("F. refusals captured: ", sum(refusals$refused), " of ",
        nrow(refusals))
stopifnot(all(refusals$refused))

## ---- Write --------------------------------------------------------------
extension <- do.call(rbind, records)
utils::write.csv(extension, file.path(results, "extension.csv"),
                 row.names = FALSE)
utils::write.csv(uncertainty, file.path(results, "extension-uncertainty.csv"),
                 row.names = FALSE)
utils::write.csv(refusals, file.path(results, "extension-refusals.csv"),
                 row.names = FALSE)

## ---- The comparison table ----------------------------------------------
w <- as.character(rdm(plans$whole)$index)[1]
pick <- function(quantity) {
  extension$value[extension$frame == "whole" &
                    extension$measurement == w &
                    extension$quantity == quantity][1]
}
plugin_se <- unique(uncertainty$rsa_category_se[
  uncertainty$frame == "whole" & uncertainty$target == "plugin"])
null_se <- unique(uncertainty$rsa_category_se[
  uncertainty$frame == "whole" & uncertainty$target == "null"])
first_pair <- paste0(pairs$left[[1]], "|", pairs$right[[1]])

# The parity number is read from the recorded agreement table, never retyped.
agreement_file <- file.path(results, "agreement.csv")
rdm_parity <- if (file.exists(agreement_file)) {
  agreement <- utils::read.csv(agreement_file, stringsAsFactors = FALSE)
  row <- agreement[agreement$quantity == "crossnobis_rdm" &
                     agreement$comparator ==
                       "rsatoolbox::calc_rdm_crossnobis", ]
  sprintf("yes, agrees to %.2e", row$max_abs_diff[[1]])
} else {
  "yes (run 03-compare.R for the recorded tolerance)"
}
rsa_parity <- if (file.exists(agreement_file)) {
  row <- agreement[agreement$quantity == "linear_rsa_coefficients", ]
  sprintf("yes, one compiled query, agrees to %.2e", row$max_abs_diff[[1]])
} else {
  "yes, one compiled query over the same plan"
}

table <- data.frame(
  quantity = c(
    "crossnobis RDM (15 pairs)",
    "linear RSA coefficient (category)",
    "signed contrast energy, face - house",
    "coherent energy (common spatial mode)",
    "configuration energy (pattern remainder)",
    "total energy = coherent + configuration",
    "coherence fraction",
    "coherent/configuration split of every RDM entry",
    "analytic SE of one RDM entry",
    "analytic SE of the RSA coefficient",
    "refusal when the estimand is not admitted"
  ),
  rsatoolbox = c(
    "yes (calc_rdm_crossnobis)",
    "yes (fit_regress / external OLS)",
    "no -- RDM entries are squared, unsigned",
    "no", "no", "no", "no", "no",
    "bootstrap / subject-level inference",
    "bootstrap / subject-level inference",
    "no typed refusal channel"
  ),
  crossform_same_fit = c(
    rdm_parity,
    rsa_parity,
    sprintf("%+.6f", pick("signed")),
    sprintf("%.6f", pick("coherent")),
    sprintf("%.6f", pick("configuration")),
    sprintf("%.6f", pick("total")),
    sprintf("%.4f", pick("coherence_fraction")),
    sprintf("yes, e.g. %s: %.6f + %.6f", first_pair,
            pick(paste0("rdm_coherent[", first_pair, "]")),
            pick(paste0("rdm_configuration[", first_pair, "]"))),
    sprintf("%.6f (plugin, %s)",
            pick(paste0("rdm SE (plugin)[", first_pair, "]")), first_pair),
    sprintf("%.6f (plugin) / %.6f (null)", plugin_se, null_se),
    sprintf("%d of %d provoked refusals are classed conditions",
            sum(refusals$refused), nrow(refusals))
  ),
  stringsAsFactors = FALSE
)
utils::write.csv(table, file.path(results, "extension-table.csv"),
                 row.names = FALSE)

message("")
message("what rsatoolbox gives / what crossform gives from the SAME fit")
message("(whole-brain measurement, contrast = mean(face) - mean(house))")
print(table, row.names = FALSE, right = FALSE)
message("")
message("Wrote extension.csv (", nrow(extension), " rows), ",
        "extension-uncertainty.csv (", nrow(uncertainty), " rows), ",
        "extension-refusals.csv, extension-table.csv")
