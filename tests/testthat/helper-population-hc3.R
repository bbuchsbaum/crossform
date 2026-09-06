# Shared public-route fixture for the production and standalone HC3 courts.

ph_subject <- function(id, features, gain) {
  effects <- effect_space(c("face", "house"), basis_id = "population-hc3:v1")
  domain <- abstract_domain(features,
    coordinates = cbind(x = seq_len(features) - 1),
    feature_ids = paste0("f", seq_len(features)), id = id)
  base <- matrix(seq_len(2L * features), 2L, features,
    dimnames = list(c("face", "house"), NULL)) / features
  rel <- relation(list(run1 = gain * base, run2 = gain * base * 1.2,
    run3 = gain * base * 0.8), effects = effects, domain = domain)
  plan_geometry(rel, compile_frame(voxelwise(), domain),
    cross_partitions(rel))
}
ph_fit <- function(covariates, gains, coverage_policy = "all_planned",
                   model = ~ x1 + x2) {
  labels <- rownames(covariates)
  sizes <- stats::setNames(10L + seq_along(labels) - 1L, labels)
  subjects <- stats::setNames(lapply(labels, function(id) {
    ph_subject(id, sizes[[id]], gains[[id]])
  }), labels)
  transport <- lapply(sizes, function(features) {
    anatomical_transport(
      native_coords = cbind(seq_len(features) - 1),
      group_coords = cbind(c(0, 4, 9)), semantics = "budget"
    )
  })
  plan <- plan_population(subjects, transport, model = model,
    data = covariates, coverage_policy = coverage_policy)
  estimate_population(plan, rbind(`face-house` = c(1, -1)))
}
