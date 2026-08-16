#!/usr/bin/env Rscript
# 05-crossnobis-uncertainty.R -- the error-channel arm (crossform only).
#
# This is NOT a comparison row. It is a different, better-defined estimand
# that rMVPA does not expose: the same cross-run geometry carried together
# with an identified error channel, so a distance can be reported with an
# analytic standard error rather than as a bare point estimate.
#
# The precomputed condition means used by 02 and 03 cannot support this --
# averaging discards the residuals. So this arm refits the SAME effects from
# the raw labelled volumes with `lm_relation_fit()`. The indicator design
# means the refitted coefficients are numerically identical to the condition
# means (asserted below), so the point estimand is unchanged; what is added
# is the residual channel that was previously thrown away.
#
# Four things are recorded:
#   A. fixed-metric RDM point estimate + analytic sampling standard errors
#   B. transport of that covariance to the RSA coefficient (a linear
#      functional of the RDM, so the transport is exact)
#   C. learned-metric crossnobis contrast (face vs house)
#   D. the two capability REFUSALS, captured verbatim. They are the point of
#      the arm as much as the numbers are: the package declines to calibrate
#      what it cannot calibrate, instead of returning a plausible number.

suppressMessages(library(neuroim2))

exemplar_dir <- if (nzchar(Sys.getenv("EXEMPLAR_DIR"))) {
  Sys.getenv("EXEMPLAR_DIR")
} else {
  normalizePath(dirname(sub("^--file=", "",
    grep("^--file=", commandArgs(FALSE), value = TRUE)[1])))
}
source(file.path(exemplar_dir, "models.R"))
invisible(load_crossform(exemplar_dir))
paths <- exemplar_paths(exemplar_dir)

prep <- readRDS(paths$prepared)
if (is.null(prep$response_runs)) {
  stop("prepared-smoke.rds carries no response blocks; delete it and rerun ",
       "01-prepare-data.R so the error channel has raw observations.")
}

# How many frame nodes get an analytic covariance, spread evenly over the VT
# centres. Bounded by default because the cost remains material. A 2026-08-14
# 25-node run on the validation-memo path measured 0.65 s/node; the script
# records every node rather than treating that machine-specific rate as a
# permanent guarantee. It is not a correctness limit -- `UNCERTAINTY_NODES=all`
# runs every centre.
node_arg <- Sys.getenv("UNCERTAINTY_NODES", "120")

## ---- Domain, frame ------------------------------------------------------
vt_vol <- neuroim2::read_vol(prep$mask_file)
vt_mask <- neuroim2::LogicalNeuroVol(as.array(vt_vol) > 0,
                                     neuroim2::space(vt_vol))
domain <- neuroim2_volume_domain(vt_mask, id = "haxby-subj1-vt")
frame <- neuroim2_searchlights(vt_mask, radius = RADIUS_MM, domain = domain,
                               normalization = "local")
n_features <- length(domain$feature_ids)
n_centers <- nrow(frame$weights)

## ---- Refit the same effects, keeping residuals --------------------------
effect_targets <- diag(length(CONDITIONS))
dimnames(effect_targets) <- list(CONDITIONS, CONDITIONS)

fit <- lm_relation_fit(prep$response_runs, prep$response_design,
                       effect_targets, effect_names = CONDITIONS,
                       sampling_unit = "scan", domain = domain)
message("lm_relation_fit: ", length(fit$relation$partitions), " partitions, ",
        nrow(prep$response_runs[[1]]), " observations x ", n_features,
        " features per partition")

over <- cross_partitions(fit$relation)

## ---- Consistency with the point-estimate arm ----------------------------
# The refit must reproduce 02's RDM exactly, or this arm is describing
# different effects and its error bars would not belong to that estimate.
nominal_plan <- plan_geometry(fit$relation, at = frame, over = over)
D_fit <- as.matrix(rdm(nominal_plan)$values)
arm2 <- file.path(paths$results, "smoke-crossform.rds")
consistency <- NA_real_
if (file.exists(arm2)) {
  consistency <- max(abs(D_fit - readRDS(arm2)$rdm))
  message("refit RDM vs 02's precomputed-mean RDM: max abs difference = ",
          format(consistency, digits = 3))
  if (consistency > TOLERANCE) {
    stop("The error-channel refit does not reproduce the point-estimate arm.")
  }
}

## ---- A. Fixed-metric RDM with analytic sampling standard errors ---------
metric <- noise_precision(diag(n_features), domain,
                          covariance = diag(n_features),
                          provenance = list(source = "fixed-identity-smoke"))
fixed_plan <- plan_geometry(fit$relation, at = frame, over = over,
                            metric = metric)
# One evaluation only: with a domain-wide metric this view is the expensive
# part of the arm, so hold on to both the values and the centre index.
fixed_view <- rdm(fixed_plan)
D_fixed <- as.matrix(fixed_view$values)
fixed_centers <- as.integer(fixed_view$index)

nodes <- if (identical(node_arg, "all")) {
  seq_len(n_centers)
} else {
  k <- as.integer(node_arg)
  unique(round(seq(1, n_centers, length.out = min(k, n_centers))))
}
message("Analytic sampling covariance for ", length(nodes), " of ",
        n_centers, " frame nodes (timed individually) ...")

n_pairs <- ncol(D_fixed)
se <- matrix(NA_real_, length(nodes), n_pairs)

# B. transport operator: the OLS RSA coefficient is a fixed linear functional
# of the 28 distances, so its variance is an exact transport of the RDM
# covariance. (A Spearman score is not linear and gets no such transport --
# that boundary is recorded in the report.)
model_vec <- model_vector(MODEL_RDMS$animacy)
X <- cbind(`(Intercept)` = 1, animacy = model_vec)
transport <- solve(crossprod(X), t(X))[2L, , drop = FALSE]
rsa_var <- rep(NA_real_, length(nodes))

# Per-node cost and process RSS are recorded so the arm's cost is measured
# rather than guessed. An earlier run appeared to slow down badly, which
# turned out to be contention from other jobs on the machine; keeping the
# instrumentation attached to the run is what settled that.
rss_kb <- function() {
  out <- suppressWarnings(try(
    system2("ps", c("-o", "rss=", "-p", Sys.getpid()), stdout = TRUE),
    silent = TRUE))
  if (inherits(out, "try-error") || !length(out)) return(NA_real_)
  as.numeric(trimws(out[1]))
}

node_seconds <- rep(NA_real_, length(nodes))
node_rss_kb <- rep(NA_real_, length(nodes))

t0 <- Sys.time()
for (k in seq_along(nodes)) {
  tk <- Sys.time()
  cv <- rdm_sampling_covariance(fixed_plan, fit, target = "plugin",
                                at = nodes[k])
  se[k, ] <- sqrt(sampling_covariance(cv))
  rsa_var[k] <- as.numeric(sampling_covariance(cv, operation = "transport",
                                               query = transport))
  node_seconds[k] <- as.numeric(difftime(Sys.time(), tk, units = "secs"))
  node_rss_kb[k] <- rss_kb()
  if (k %% 20 == 0) {
    el <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
    message("  node ", k, "/", length(nodes), "  (", round(el, 0), " s cum, ",
            round(mean(tail(node_seconds[seq_len(k)], 20)), 2),
            " s/node recent, RSS ", round(node_rss_kb[k] / 1024), " MB)")
  }
}
elapsed_cov <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
message("sampling covariance: ", round(elapsed_cov, 1), " s total, ",
        round(elapsed_cov / length(nodes), 2), " s/node mean")
message("  first 10 nodes: ", round(mean(head(node_seconds, 10)), 2),
        " s/node, RSS ", round(head(node_rss_kb, 1) / 1024), " MB")
message("  last 10 nodes : ", round(mean(tail(node_seconds, 10)), 2),
        " s/node, RSS ", round(tail(node_rss_kb, 1) / 1024), " MB")

rsa_se <- sqrt(rsa_var)
rsa_point <- as.numeric(as.matrix(
  rsa(fixed_plan, models = list(animacy = MODEL_RDMS$animacy))$coefficients
)[nodes, "animacy"])
rsa_z <- rsa_point / rsa_se

message("RSA coefficient with analytic SE: median estimate ",
        round(median(rsa_point), 4), ", median SE ", round(median(rsa_se), 4),
        ", median |z| ", round(median(abs(rsa_z)), 3))

## ---- C. Learned-metric crossnobis contrast ------------------------------
weights <- setNames(rep(0, length(CONDITIONS)), CONDITIONS)
weights[["face"]] <- 1
weights[["house"]] <- -1

cn_plan <- plan_crossnobis(fit, frame, over, metric = shrinkage_precision(0.2))
t0 <- Sys.time()
cn <- crossnobis(cn_plan, weights)
elapsed_cn <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
cn_df <- as.data.frame(cn)
message("learned-metric crossnobis (face - house) over ", nrow(cn_df),
        " centres in ", round(elapsed_cn, 1), " s; median ",
        round(median(cn_df$crossnobis), 4))

## ---- D. Capability refusals, captured verbatim --------------------------
capture_refusal <- function(label, expr) {
  msg <- tryCatch({ force(expr); NA_character_ },
                  error = function(e) conditionMessage(e))
  if (is.na(msg)) {
    message("NOT REFUSED: ", label, " -- expected a refusal and got a value")
  } else {
    message("refusal captured: ", label)
  }
  list(label = label, message = msg)
}

refusals <- list(
  capture_refusal(
    "analytic sampling covariance on a learned-metric crossnobis plan",
    rdm_sampling_covariance(cn_plan, fit, target = "plugin", at = 1L)),
  capture_refusal(
    "analytic sampling covariance from precomputed effects (no error channel)",
    {
      rel0 <- relation(prep$runs,
                       effects = effect_space(CONDITIONS,
                                              basis_id = "haxby-condition-means:v1"),
                       domain = domain)
      p0 <- plan_geometry(rel0, at = frame, over = cross_partitions(rel0),
                          metric = metric)
      rdm_sampling_covariance(p0, rel0, target = "plugin", at = 1L)
    }),
  capture_refusal(
    "sampling covariance target not inferred when `target` is omitted",
    rdm_sampling_covariance(fixed_plan, fit, at = 1L))
)

## ---- Save ---------------------------------------------------------------
out <- list(
  arm = "crossform error channel",
  note = paste("A different, better-defined estimand than the parity",
               "comparison -- not a replacement-map row."),
  centers = fixed_centers,
  nodes = nodes,
  node_centers = fixed_centers[nodes],
  rdm_fixed_metric = D_fixed,
  rdm_se = se,
  pairs = rdm_pair_index()[, c("left", "right")],
  rsa = data.frame(node = nodes, estimate = rsa_point, se = rsa_se, z = rsa_z),
  crossnobis = cn_df,
  crossnobis_weights = weights,
  crossnobis_metric = "shrinkage_precision(0.2), learned per node",
  refusals = refusals,
  consistency_vs_arm2 = consistency,
  timing = list(sampling_covariance_seconds = elapsed_cov,
                seconds_per_node = elapsed_cov / length(nodes),
                node_seconds = node_seconds,
                node_rss_kb = node_rss_kb,
                crossnobis_seconds = elapsed_cn),
  session = list(crossform = as.character(utils::packageVersion("crossform")),
                 when = Sys.time())
)
saveRDS(out, file.path(paths$results, "smoke-uncertainty.rds"))
message("Wrote ", file.path(paths$results, "smoke-uncertainty.rds"))
