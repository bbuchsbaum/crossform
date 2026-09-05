pg_fixture <- function(n = 4L, origins = TRUE) {
  V <- qr.Q(qr(cbind(1, diag(n))))[, 2:3, drop = FALSE]
  nm <- paste0("c", seq_len(n)); rownames(V) <- nm
  B <- V %*% diag(sqrt(c(6, 3)))
  parts <- list(a = B, b = B, c = B * .8, d = B * 1.1)
  domain <- abstract_domain(2L, id = "frozen-fit")
  origins <- if (origins) list(observation_origins = list(id = "frozen-fit-origins",
    partitions = list(a = "raw-a", b = "raw-b", c = "raw-c", d = "raw-d"),
    independence = "independent", assumption = "Independent runs with fixed external preprocessing.")) else list()
  rel <- relation(parts, domain = domain, provenance = origins)
  at <- additive_frame(members = list(1:2, 1L),
    domain = domain, measurements = c("first", "second"))
  plan <- plan_geometry(rel, at, pairing("a", "b", independence = "independent", generalizes_over = "run"))
  K <- V %*% diag(c(.6, .4)) %*% t(V); dimnames(K) <- list(nm, nm)
  basis <- model_basis(kernels = list(m = K), conditions = rel$effect_space, normalize = "trace")
  list(plan = plan, basis = basis, V = V, K = K, B = B, rel = rel, at = at, domain = domain)
}
pg_forms <- function(fit) {
  lapply(crossform:::.geometry_prediction_rows(fit, seq_len(length(fit$target$index))),
    function(row) fit$model$Q %*% tcrossprod(row$factor) %*% t(fit$model$Q))
}

pg_plan <- function(f, left = "c", right = "d", rel = f$rel, at = f$at,
                    metric = NULL, independence = "independent", generalizes_over = "run",
                    compute = compute_policy()) {
  plan_geometry(rel, at, pairing(left, right,
    independence = independence, generalizes_over = generalizes_over), metric = metric,
    compute = compute)
}
