pg_nested_fixture <- function(test_scale = 1) {
  f <- pg_fixture()
  first <- c(1, 1.1, .9, .8, 1.2, 1, .9, 1.05)
  second <- c(1, .8, 1.2, -.2, .4, .9, .6, .8)
  parts <- stats::setNames(lapply(1:8, function(i) f$B %*% diag(c(first[i], second[i]))), letters[1:8])
  parts$g <- test_scale * parts$g
  manifest <- list(id = "nested-eight-runs", partitions = as.list(stats::setNames(paste0("raw-", letters[1:8]), letters[1:8])),
    independence = "independent", assumption = "Independent runs; fixed externally defined effects and preprocessing.")
  f$rel <- relation(parts, domain = f$domain, provenance = list(observation_origins = manifest))
  edges <- t(combn(letters[1:6], 2))
  f$train <- plan_geometry(f$rel, f$at, pairing(edges[, 1], edges[, 2], weight = seq_len(nrow(edges)),
    independence = "independent", generalizes_over = "run"))
  f$test <- pg_plan(f, "g", "h")
  f$inner <- list(list(train = letters[1:4], validate = letters[5:6], weight = .2),
    list(train = letters[3:6], validate = letters[1:2], weight = .3),
    list(train = letters[c(1, 2, 5, 6)], validate = letters[3:4], weight = .5))
  f$candidates <- list(list(basis = f$basis, rank = 0L, penalty = .2),
    list(basis = f$basis, rank = 1L, penalty = 0),
    list(basis = f$basis, rank = 1L, penalty = .2),
    list(basis = f$basis, rank = 2L, penalty = .2),
    list(basis = f$basis, rank = 2L, penalty = .5))
  f
}
