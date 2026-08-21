pdec_fixture <- function(coverage_policy = "available_at_node",
                         normalization = "none", outlier_subject = NULL) {
  effects <- effect_space(c("face", "house", "tool"), basis_id = "pdec:v1")
  sizes <- c(s01 = 6L, s02 = 7L, s03 = 8L, s04 = 9L,
             s05 = 10L, s06 = 11L, s07 = 12L, s08 = 13L)
  gains <- stats::setNames(seq(0.7, 1.4, length.out = length(sizes)),
                           names(sizes))
  if (!is.null(outlier_subject)) gains[[outlier_subject]] <- 8
  subjects <- stats::setNames(lapply(names(sizes), function(id) {
    n <- sizes[[id]]
    domain <- abstract_domain(n, coordinates = cbind(x = seq_len(n) - 1),
      feature_ids = paste0("f", seq_len(n)), id = id)
    base <- matrix(gains[[id]] * seq_len(3L * n), 3L, n,
      dimnames = list(c("face", "house", "tool"), NULL)) / n
    relation <- relation(list(run1 = base, run2 = base * 1.13,
      run3 = base * 0.91, run4 = base * 1.04),
      effects = effects, domain = domain)
    plan_geometry(relation, compile_frame(voxelwise(), domain),
                  cross_partitions(relation))
  }), names(sizes))
  transports <- lapply(sizes, function(n) anatomical_transport(
    cbind(seq_len(n) - 1), cbind(c(0, 5, 12)), semantics = "budget",
    radius = 1.5
  ))
  covariates <- data.frame(age = seq(-1, 1, length.out = length(sizes)),
                           row.names = names(sizes))
  plan_population(subjects, transports, model = ~ age, data = covariates,
    coverage_policy = coverage_policy, normalization = normalization)
}

pdec_results <- function(plan = pdec_fixture()) {
  bank <- rbind(`face-house` = c(1, -1, 0),
                `face-tool` = c(1, 0, -1))
  list(
    total = estimate_population(plan, bank, component = "total"),
    coherent = estimate_population(plan, bank, component = "coherent"),
    configuration = estimate_population(plan, bank,
      component = "configuration")
  )
}
