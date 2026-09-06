# Run from any directory with an installed crossform:
# Rscript exemplars/predictive-geometry/fixed-split.R
library(crossform)

# The neural/model directions are centered and orthonormal. All numerical
# values below are deterministic: independent origins describe the intended
# sampling design, which here has zero noise to expose the arithmetic.
condition_names <- c("face", "body", "house", "tool")
V <- cbind(c(1, -1, 0, 0), c(0, 0, 1, -1)) / sqrt(2)
rownames(V) <- condition_names
model_values <- c(.6, .4)
model_features <- V %*% diag(sqrt(model_values))
K <- tcrossprod(model_features)
dimnames(K) <- list(condition_names, condition_names)
D <- outer(diag(K), diag(K), "+") - 2 * K # explicitly squared distances

penalty <- .1
B_train <- V %*% diag(sqrt(c(2, 1) + penalty / model_values))
B_test <- V %*% diag(sqrt(c(2, 0)))
domain <- abstract_domain(2, id = "predictive-example")
origin_manifest <- list(id = "four-independent-runs",
  partitions = list(train1 = "raw-run1", train2 = "raw-run2",
    test1 = "raw-run3", test2 = "raw-run4"),
  independence = "independent",
  assumption = "Independent runs; effect extraction and preprocessing were fixed externally.")
rel <- relation(list(train1 = B_train, train2 = B_train,
  test1 = B_test, test2 = B_test), domain = domain,
  provenance = list(observation_origins = origin_manifest))
# Sum over the two neural features; local normalization would divide by two
# and therefore change both the geometry scale and the appropriate penalty.
at <- compile_frame(whole_brain(normalization = "none"), domain)
train_plan <- plan_geometry(rel, at, pairing("train1", "train2",
  independence = "independent", generalizes_over = "run"))
test_plan <- plan_geometry(rel, at, pairing("test1", "test2",
  independence = "independent", generalizes_over = "run"))

bases <- list(
  kernel = model_basis(kernels = list(model = K), conditions = rel$effect_space, normalize = "trace"),
  features = model_basis(features = list(model = model_features), conditions = rel$effect_space, normalize = "trace"),
  squared_rdm = model_basis(models = list(model = D), distance = "squared_euclidean",
    conditions = rel$effect_space, normalize = "trace"))
fits <- lapply(bases, function(basis) fit_geometry(train_plan, basis, rank = 2, penalty = penalty))
scores <- lapply(fits, function(fit) score_geometry(fit, test_plan, modes = TRUE, components = TRUE))
fit <- fits$kernel
evidence <- scores$kernel
mode_evidence <- as.data.frame(evidence, view = "modes")
zero <- score_geometry(fit_geometry(train_plan, bases$kernel, rank = 0, penalty = penalty), test_plan)

# A reduced rank may be DECLARED before evaluation. This example reports both
# declared fits; it does not choose rank by inspecting final test gains.
rank_one <- score_geometry(fit_geometry(train_plan, bases$kernel, rank = 1, penalty = penalty), test_plan)
stopifnot(abs(evidence$table$gain - 3) < 1e-10,
  max(abs(mode_evidence$amplitude - c(2, 1))) < 1e-10,
  max(abs(mode_evidence$evidence - c(2, 0))) < 1e-10,
  max(abs(mode_evidence$gain - c(4, -1))) < 1e-10,
  zero$table$gain == 0, abs(rank_one$table$gain - 4) < 1e-10,
  all(vapply(scores, function(x) abs(x$table$gain - 3) < 1e-10, FALSE)))

print(evidence)
print(mode_evidence, row.names = FALSE)
cat("The second frozen mode costs 1 and has no replicated evidence: its gain is -1.\n")
cat("Total gain is 3 squared geometry units; the zero predictor has exactly zero gain.\n")
cat("This evaluates independent runs on the same conditions. It does not test new-condition prediction.\n")

# These are refusals, not ways to bypass independent evaluation.
overlap <- catch_refusal(score_geometry(fit, train_plan))
unsupported_loss <- catch_refusal(fit_geometry(train_plan, bases$kernel,
  rank = 2, penalty = penalty, loss = "rdm_gls"))
stopifnot("overlapping_training_evaluation_origins" %in% overlap$reasons,
  "predictive_loss_not_implemented" %in% unsupported_loss$reasons)
