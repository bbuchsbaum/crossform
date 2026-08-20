#!/usr/bin/env Rscript
# 03-recover.R -- planted versus estimated, and the match-versus-control
# tests. Reads only the committed CSV from 02-analyze.R plus the planted
# truth, and does its own statistics; it calls no crossform function, so the
# recovery check is external to the thing being checked.
#
# Writes:
#   results/planted-vs-estimated.csv   every query: planted, estimate, bias
#   results/coupling-levels.csv        the headline match/control levels
#   results/recovery-verdicts.csv      the qualitative claims, each with a test

exemplar_dir <- if (nzchar(Sys.getenv("EXEMPLAR_DIR"))) {
  Sys.getenv("EXEMPLAR_DIR")
} else {
  normalizePath(dirname(sub("^--file=", "",
    grep("^--file=", commandArgs(FALSE), value = TRUE)[1])))
}
source(file.path(exemplar_dir, "00-common.R"))
results <- file.path(exemplar_dir, "results")

estimates <- read.csv(file.path(results, "estimates-by-subject.csv"),
  stringsAsFactors = FALSE)
n_subjects <- length(unique(estimates$subject))
message("subjects: ", n_subjects, "  |  rows: ", nrow(estimates))

## ---- Planted versus estimated -------------------------------------------
summarise <- function(block) {
  value <- block$estimate
  planted <- unique(block$planted)
  stopifnot(length(planted) == 1L)
  mean_value <- mean(value)
  se <- stats::sd(value) / sqrt(length(value))
  data.frame(
    plan = block$plan[[1L]], query = block$query[[1L]],
    region = block$region[[1L]], readout = block$readout[[1L]],
    planted = planted, estimate = mean_value, se = se,
    bias = mean_value - planted,
    bias_t = (mean_value - planted) / se,
    bias_p = 2 * stats::pt(abs((mean_value - planted) / se),
      df = length(value) - 1L, lower.tail = FALSE),
    t_vs_zero = mean_value / se,
    p_vs_zero = 2 * stats::pt(abs(mean_value / se),
      df = length(value) - 1L, lower.tail = FALSE),
    ci_low = mean_value - stats::qt(0.975, length(value) - 1L) * se,
    ci_high = mean_value + stats::qt(0.975, length(value) - 1L) * se,
    stringsAsFactors = FALSE
  )
}
blocks <- split(estimates, list(estimates$plan, estimates$query,
  estimates$region), drop = TRUE)
recovery <- do.call(rbind, lapply(blocks, summarise))
recovery <- recovery[order(recovery$plan, recovery$query, recovery$region), ]
rownames(recovery) <- NULL
write.csv(
  cbind(recovery[, c("plan", "query", "region", "readout")],
    round(recovery[, c("planted", "estimate", "se", "bias", "bias_t",
      "bias_p", "t_vs_zero", "p_vs_zero", "ci_low", "ci_high")], 6)),
  file.path(results, "planted-vs-estimated.csv"), row.names = FALSE
)

pick <- function(plan, query, region) {
  row <- recovery[recovery$plan == plan & recovery$query == query &
    recovery$region == region, ]
  stopifnot(nrow(row) == 1L)
  row
}

cross <- recovery[recovery$plan == "cross", ]
covered <- sum(cross$planted >= cross$ci_low & cross$planted <= cross$ci_high)
message("\nplanted value inside the 95% subject CI for ", covered, " of ",
        nrow(cross), " cross-run readouts")
message("largest |bias| over cross-run readouts: ",
        format(max(abs(cross$bias)), digits = 3),
        "  (largest |bias| / se = ", format(max(abs(cross$bias_t)),
          digits = 3), ")")

## ---- The headline coupling levels ---------------------------------------
levels_table <- do.call(rbind, lapply(names(REGION_SIZES), function(region) {
  match_level <- pick("cross", "level_match", region)
  control_all <- pick("cross", "level_control_all", region)
  control_category <- pick("cross", "level_control_category", region)
  contrast_category <- pick("cross", "contrast_category", region)
  contrast_all <- pick("cross", "contrast_all", region)
  data.frame(
    region = region,
    match_coupling = match_level$estimate,
    control_coupling_all = control_all$estimate,
    control_coupling_category = control_category$estimate,
    item_specific_contrast = contrast_category$estimate,
    item_specific_planted = contrast_category$planted,
    item_specific_t = contrast_category$t_vs_zero,
    naive_contrast = contrast_all$estimate,
    naive_planted = contrast_all$planted,
    stringsAsFactors = FALSE
  )
}))
write.csv(cbind(levels_table[, "region", drop = FALSE],
  round(levels_table[, -1L], 4)),
  file.path(results, "coupling-levels.csv"), row.names = FALSE)
message("\nmatch and control couplings, cross-run plan:")
print(cbind(levels_table[, "region", drop = FALSE],
  round(levels_table[, -1L], 4)), row.names = FALSE)

## ---- Qualitative verdicts, each with an explicit test -------------------
verdict <- function(claim, statistic, criterion, passed, detail) {
  data.frame(claim = claim, statistic = statistic, criterion = criterion,
    passed = passed, detail = detail, stringsAsFactors = FALSE)
}
significant <- function(row) row$p_vs_zero < 0.001 && row$estimate > 0
null_result <- function(row) row$ci_low <= 0 && row$ci_high >= 0

a_item <- pick("cross", "contrast_category", "regionA")
b_item <- pick("cross", "contrast_category", "regionB")
c_item <- pick("cross", "contrast_category", "regionC")
b_naive <- pick("cross", "contrast_all", "regionB")
c_same <- pick("cross", "contrast_category", "regionC")
c_same_run <- pick("same", "contrast_category", "regionC")
a_depth <- pick("cross", "lm_match_depth", "regionA")
a_lm <- pick("cross", "lm_match", "regionA")
a_mc <- pick("cross", "match_control_category", "regionA")

report <- function(row) {
  sprintf("%.4f (planted %.4f, se %.4f, t = %.2f)", row$estimate, row$planted,
    row$se, row$t_vs_zero)
}

verdicts <- rbind(
  verdict("regionA carries item-specific reinstatement",
    "coupling_contrast, category-matched controls", "t > 0, p < 0.001",
    significant(a_item), report(a_item)),
  verdict("regionB carries category structure but no item reinstatement",
    "coupling_contrast, category-matched controls", "95% CI contains 0",
    null_result(b_item), report(b_item)),
  verdict("regionC carries nothing",
    "coupling_contrast, category-matched controls", "95% CI contains 0",
    null_result(c_item), report(c_item)),
  verdict("an unrestricted control set reports reinstatement in regionB",
    "coupling_contrast, all eligible controls", "t > 0, p < 0.001",
    significant(b_naive), report(b_naive)),
  verdict("same-run pairing invents reinstatement in the null region",
    "coupling_contrast under the same-run plan", "t > 0, p < 0.001",
    significant(c_same_run),
    paste0(report(c_same_run), " vs cross-run ", report(c_same))),
  verdict("the encoding study-duration covariate is recovered in regionA",
    "pair_lm_query match_depth coefficient", "t > 0, p < 0.001",
    significant(a_depth), report(a_depth)),
  verdict("the regression and contrast routes agree in regionA",
    "pair_lm_query match vs match_control vs coupling_contrast",
    "spread < 0.05",
    max(c(a_lm$estimate, a_mc$estimate, a_item$estimate)) -
      min(c(a_lm$estimate, a_mc$estimate, a_item$estimate)) < 0.05,
    sprintf("pair_lm_query %.4f, match_control %.4f, coupling_contrast %.4f",
      a_lm$estimate, a_mc$estimate, a_item$estimate)),
  verdict("every cross-run readout recovers its planted value",
    "planted value vs 95% subject CI",
    sprintf("all %d readouts covered", nrow(cross)),
    covered == nrow(cross),
    sprintf("%d of %d covered; largest |bias| = %.4f", covered, nrow(cross),
      max(abs(cross$bias))))
)
write.csv(verdicts, file.path(results, "recovery-verdicts.csv"),
  row.names = FALSE)

message("\nverdicts:")
for (i in seq_len(nrow(verdicts))) {
  message(sprintf("  [%s] %s\n        %s",
    if (verdicts$passed[[i]]) "PASS" else "FAIL",
    verdicts$claim[[i]], verdicts$detail[[i]]))
}
if (!all(verdicts$passed)) {
  message("\nNOTE: at least one verdict failed; the tables above are the ",
          "record, not the claim.")
}
message("\nwrote ", file.path(results, "planted-vs-estimated.csv"))
