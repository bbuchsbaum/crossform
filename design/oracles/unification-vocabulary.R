# Independent executable court for `design/unification-contract.md`.
#
# This script deliberately does not load crossform.  It constructs the common
# geometry, its first moments, the coherent/configuration split, a fixed query,
# a conservative frame family, and a scale profile from base matrix algebra.
# `tests/testthat/test-unification-contract.R` compares these values with the
# public package route.

tol <- 1e-12

# A representation is a partition-indexed map from experimental coordinates
# (rows) to neural features (columns).  The two partitions are deliberately
# not identical, so the court exercises cross-generalization rather than a
# self-product shortcut.
representation <- list(
  run1 = rbind(a = c(2.0, -2.0, 1.0, 1.0), b = c(0, 0, 0, 0)),
  run2 = rbind(a = c(1.0, -3.0, 1.0, 1.0), b = c(0, 0, 0, 0))
)
effects <- rownames(representation[[1L]])
features <- paste0("v", seq_len(ncol(representation[[1L]])))
for (value in representation) colnames(value) <- features

# One undirected cross-partition edge is the two half-weight ordered edges.
Gamma <- matrix(c(0, 0.5, 0.5, 0), 2L, 2L,
  dimnames = list(names(representation), names(representation)))

node_geometry <- function(weight) {
  mass <- sum(weight)
  stopifnot(length(weight) == length(features), mass > 0)
  total <- coherent <- matrix(0, length(effects), length(effects),
    dimnames = list(effects, effects))
  first <- matrix(0, length(effects), length(representation),
    dimnames = list(effects, names(representation)))
  for (r in seq_along(representation)) {
    first[, r] <- drop(representation[[r]] %*% weight) / mass
    for (s in seq_along(representation)) {
      total <- total + Gamma[r, s] *
        representation[[r]] %*% diag(weight) %*% t(representation[[s]])
      coherent <- coherent + Gamma[r, s] *
        tcrossprod(representation[[r]] %*% weight,
          representation[[s]] %*% weight) / mass
    }
  }
  endpoint_first <- drop(first %*% rowSums(Gamma))
  names(endpoint_first) <- effects
  list(
    mass = mass,
    first = endpoint_first,
    total = total,
    coherent = coherent,
    configuration = total - coherent
  )
}

# Each member is conservative on its own.  Alpha scales complete members; the
# stacked family therefore has unit column mass and retains alpha as part of
# the estimand.
fine <- diag(4L)
coarse <- rbind(
  v1 = c(1 / 2, 1 / 3, 0, 0),
  v2 = c(1 / 2, 1 / 3, 1 / 3, 0),
  v3 = c(0, 1 / 3, 1 / 3, 1 / 2),
  v4 = c(0, 0, 1 / 3, 1 / 2)
)
alpha <- c(`radius-0.5` = 0.5, `radius-1.01` = 0.5)
frame_family <- rbind(
  `radius-0.5` = alpha[["radius-0.5"]] * fine,
  `radius-1.01` = alpha[["radius-1.01"]] * coarse
)
family <- c(rep("radius-0.5", nrow(fine)),
  rep("radius-1.01", nrow(coarse)))
scale <- c(rep(0.5, nrow(fine)), rep(1.01, nrow(coarse)))
node <- rep(features, 2L)
rownames(frame_family) <- paste(family, node, sep = "::")

geometry <- lapply(seq_len(nrow(frame_family)), function(i) {
  node_geometry(frame_family[i, ])
})
names(geometry) <- rownames(frame_family)

# A fixed query is chosen before any geometry is inspected.  H = c c' reads
# the signed a-b contrast energy from every second-moment component; c reads
# the corresponding first moment.
contrast <- c(a = 1, b = -1)
H <- tcrossprod(contrast)
query <- function(G) drop(sum(H * G))
first_query <- vapply(geometry,
  function(value) drop(crossprod(contrast, value$first)), numeric(1))
second_query <- data.frame(
  family = family,
  scale = scale,
  node = node,
  coherent = vapply(geometry, function(value) query(value$coherent), numeric(1)),
  configuration = vapply(geometry,
    function(value) query(value$configuration), numeric(1)),
  total = vapply(geometry, function(value) query(value$total), numeric(1)),
  row.names = rownames(frame_family), check.names = FALSE
)

scale_profile <- do.call(rbind, lapply(names(alpha), function(label) {
  rows <- second_query$family == label
  coherent <- sum(second_query$coherent[rows])
  configuration <- sum(second_query$configuration[rows])
  total <- sum(second_query$total[rows])
  data.frame(
    family = label,
    scale = unique(second_query$scale[rows]),
    alpha = unname(alpha[[label]]),
    coherent = coherent,
    configuration = configuration,
    total = total,
    coherence_fraction = if (total > 0 && coherent >= 0 && configuration >= 0) {
      coherent / total
    } else {
      NA_real_
    }
  )
}))
rownames(scale_profile) <- NULL

whole <- node_geometry(rep(1, length(features)))
whole_query <- query(whole$total)

# Algebraic gates.  These are definitions and exact finite-dimensional laws,
# not simulation claims.
stopifnot(
  identical(dim(representation[[1L]]), c(2L, 4L)),
  max(abs(rowSums(Gamma) - c(0.5, 0.5))) < tol,
  max(abs(colSums(frame_family) - 1)) < tol,
  all(vapply(geometry, function(value) {
    max(abs(value$total - value$coherent - value$configuration)) < tol
  }, logical(1))),
  max(abs(second_query$total - second_query$coherent -
    second_query$configuration)) < tol,
  max(abs(scale_profile$total - alpha * whole_query)) < tol,
  max(abs(scale_profile$total - scale_profile$coherent -
    scale_profile$configuration)) < tol,
  abs(first_query[[1L]] -
    sum(contrast * geometry[[1L]]$first)) < tol,
  abs(second_query$total[[1L]] -
    drop(crossprod(contrast, geometry[[1L]]$total %*% contrast))) < tol,
  scale_profile$coherence_fraction[[1L]] >
    scale_profile$coherence_fraction[[2L]]
)

unification_oracle <- list(
  representation = representation,
  pairing = Gamma,
  frame_family = frame_family,
  frame_index = data.frame(
    measurement = rownames(frame_family), family = family, scale = scale,
    node = node, alpha = unname(alpha[family]), stringsAsFactors = FALSE
  ),
  geometry = geometry,
  query = list(contrast = contrast, operator = H),
  first_query = first_query,
  second_query = second_query,
  scale_profile = scale_profile,
  whole = whole,
  whole_query = whole_query
)

message(sprintf(
  paste0("unification-v1 PASS: %d representations, %d family nodes, ",
    "max recomposition error %.3g"),
  length(representation), nrow(frame_family),
  max(abs(second_query$total - second_query$coherent -
    second_query$configuration))
))
