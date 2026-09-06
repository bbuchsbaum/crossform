# Workloads fixed before measuring implementation performance.
pg_performance_config <- function() list(
  schema = "predictive-performance-v1", seed = 2026090501L,
  cases = data.frame(id = c("tiny", "medium", "large"),
    conditions = c(12L, 64L, 200L), model_rank = c(4L, 12L, 20L),
    partitions = c(4L, 8L, 8L), features = c(64L, 2048L, 20000L),
    nodes = c(32L, 2000L, 10000L), support = c(8L, 32L, 32L),
    response_rank = c(2L, 6L, 8L)),
  penalty = .005, block_features = 128L, row_block = 64L,
  workspace_bytes = 512 * 1024^2, large_oracle_rows = 64L,
  relative_tolerance = 1e-8,
  scope = "same-condition independent partitions; total form; identity metric; local frame normalization",
  outputs = "prediction form, signed/shifted spectra, per-mode amplitudes/evidence/gain and total gain",
  timing = "instrumented serial execution; no universal speedup threshold")
