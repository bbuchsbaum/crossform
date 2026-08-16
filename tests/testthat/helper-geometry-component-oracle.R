# An independent first-principles oracle for the total / coherent /
# configuration decomposition of a spatial-frame geometry. It calls no
# production crossform function, so an assertion made against it cannot be
# satisfied by the implementation's own definition of `configuration`.
#
# For one frame node with support weights w over its features and a fixed
# neural metric K, write
#
#   local_K = diag(sqrt(w)) K diag(sqrt(w))          weighted local metric
#   a       = w / sum(w)                             normalized mass profile
#
# The total cross product of one ordered partition edge (L, R) is
# `T = L local_K R'`. The coherent component is the part of `local_K` carried
# by the frame-weighted mean. That functional is `x -> a'x`; its Riesz
# representative under the local metric is `v` with `local_K v = a`, and the
# coherent operator is `local_K` restricted to that one direction:
#
#   P       = v v' / (v' local_K v)                  local_K-orthogonal onto v
#   K_coh   = local_K P local_K                      rank one, PSD
#   C       = L K_coh R'
#
# `configuration = T - C` then annihilates `v` exactly, which is the property
# that makes the split canonical. This oracle builds `P` explicitly rather
# than reusing the closed form the compiler evaluates.
#
# design/effect-form-contract.md section 7 states the K = I specialization:
# there `v = w`, `v' local_K v = sum(w)`, and `C = (L w)(R w)' / sum(w)`.
# test-crossnobis-known.R checks both that specialization and the annihilation
# property, so the metric-general form is anchored to the written contract
# rather than assumed.
#
# `relation_values` is a named list of effect-by-feature matrices, one per
# partition. `partition_edges` is a data frame of `left`, `right`, `weight`.
# Set `symmetrize = FALSE` for a rectangular cross-axis form, where the two
# experimental axes differ and no half-edge expansion applies.

geometry_component_oracle <- function(relation_values, frame_weights,
                                      partition_edges, metric = NULL,
                                      symmetrize = TRUE,
                                      right_values = relation_values) {
  frame_weights <- as.matrix(frame_weights)
  features <- ncol(frame_weights)
  metric <- if (is.null(metric)) diag(features) else as.matrix(metric)
  partition_edges <- as.data.frame(partition_edges)
  q_left <- nrow(relation_values[[1L]])
  q_right <- nrow(right_values[[1L]])
  nodes <- lapply(seq_len(nrow(frame_weights)), function(node) {
    support <- which(frame_weights[node, ] > 0)
    weight <- frame_weights[node, support]
    root <- sqrt(weight)
    local_metric <- metric[support, support, drop = FALSE] * tcrossprod(root)
    profile <- weight / sum(weight)
    representative <- drop(solve(local_metric, profile))
    projector <- tcrossprod(representative) /
      drop(crossprod(representative, local_metric %*% representative))
    coherent_metric <- local_metric %*% projector %*% local_metric
    total <- matrix(0, q_left, q_right)
    coherent <- matrix(0, q_left, q_right)
    for (edge in seq_len(nrow(partition_edges))) {
      left <- relation_values[[partition_edges$left[[edge]]]][
        , support, drop = FALSE
      ]
      right <- right_values[[partition_edges$right[[edge]]]][
        , support, drop = FALSE
      ]
      edge_total <- left %*% local_metric %*% t(right)
      edge_coherent <- left %*% coherent_metric %*% t(right)
      if (isTRUE(symmetrize)) {
        edge_total <- (edge_total + t(edge_total)) / 2
        edge_coherent <- (edge_coherent + t(edge_coherent)) / 2
      }
      mass <- partition_edges$weight[[edge]]
      total <- total + mass * edge_total
      coherent <- coherent + mass * edge_coherent
    }
    list(total = total, coherent = coherent,
      configuration = total - coherent)
  })
  list(
    total = lapply(nodes, `[[`, "total"),
    coherent = lapply(nodes, `[[`, "coherent"),
    configuration = lapply(nodes, `[[`, "configuration")
  )
}

# The same decomposition read through one fixed contrast, giving one scalar
# per frame node: the quantity `contrast_energy()` and `crossnobis()` report.
geometry_contrast_oracle <- function(weights, ...) {
  components <- geometry_component_oracle(...)
  lapply(components, function(component) {
    vapply(component, function(value) {
      drop(weights %*% value %*% weights)
    }, numeric(1))
  })
}

# Read a relation's partition blocks without going through any geometry
# machinery. This is a convenience for fit-backed relations, whose stored
# values are only reachable through the extractor; where the raw source
# matrices are available, pass those instead, since they share nothing at all
# with the code under test.
geometry_oracle_relation_values <- function(relation) {
  values <- lapply(relation$partitions, function(partition) {
    relation_block(relation, partition, seq_len(relation$n_features))
  })
  stats::setNames(values, relation$partitions)
}
