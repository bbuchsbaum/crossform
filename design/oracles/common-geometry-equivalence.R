# Independent slow oracle for design/common-geometry-equivalence.md.
# Base matrix arithmetic only; no crossform implementation helper is called.

oracle_bilinear_inputs <- function(operator, left, right, gamma, metric) {
  if (!is.list(left) || !length(left) || !is.list(right) || !length(right) ||
      !all(vapply(left, is.matrix, logical(1))) ||
      !all(vapply(right, is.matrix, logical(1)))) {
    stop("oracle blocks must be nonempty lists of matrices")
  }
  q_left <- nrow(left[[1L]])
  d_left <- ncol(left[[1L]])
  q_right <- nrow(right[[1L]])
  d_right <- ncol(right[[1L]])
  if (any(vapply(left, nrow, integer(1)) != q_left) ||
      any(vapply(left, ncol, integer(1)) != d_left) ||
      any(vapply(right, nrow, integer(1)) != q_right) ||
      any(vapply(right, ncol, integer(1)) != d_right) ||
      any(!is.finite(unlist(left, use.names = FALSE))) ||
      any(!is.finite(unlist(right, use.names = FALSE)))) {
    stop("oracle blocks must have stable finite dimensions")
  }
  if (!is.matrix(operator) || !identical(dim(operator), c(q_left, q_right)) ||
      any(!is.finite(operator))) {
    stop("oracle operator is not aligned to the effect axes")
  }
  if (!is.matrix(metric) || !identical(dim(metric), c(d_left, d_right)) ||
      any(!is.finite(metric))) {
    stop("oracle metric is not aligned to the feature axes")
  }
  if (!is.matrix(gamma) ||
      !identical(dim(gamma), c(length(left), length(right))) ||
      any(!is.finite(gamma))) {
    stop("oracle Gamma is not aligned to the partition axes")
  }

  left_names <- names(left)
  right_names <- names(right)
  gamma_names <- dimnames(gamma)
  named <- !is.null(left_names) || !is.null(right_names) ||
    !is.null(gamma_names[[1L]]) || !is.null(gamma_names[[2L]])
  if (named) {
    valid_names <- function(x, n) is.character(x) && length(x) == n &&
      !anyNA(x) && !any(!nzchar(x)) && !anyDuplicated(x)
    if (!valid_names(left_names, length(left)) ||
        !valid_names(right_names, length(right)) ||
        !valid_names(gamma_names[[1L]], nrow(gamma)) ||
        !valid_names(gamma_names[[2L]], ncol(gamma)) ||
        !setequal(left_names, gamma_names[[1L]]) ||
        !setequal(right_names, gamma_names[[2L]])) {
      stop("oracle partition labels do not bind Gamma to the blocks")
    }
    gamma <- gamma[left_names, right_names, drop = FALSE]
  }
  list(operator = operator, left = left, right = right, gamma = gamma,
    metric = metric)
}

# Explicit matrix formula. It supports square self forms and rectangular
# left/right forms; symmetry is a property of the former, not assumed here.
oracle_bilinear <- function(operator, left, right, gamma, metric) {
  input <- oracle_bilinear_inputs(operator, left, right, gamma, metric)
  total <- 0
  for (r in seq_along(input$left)) {
    for (s in seq_along(input$right)) {
      total <- total + input$gamma[r, s] * sum(input$operator *
        (input$left[[r]] %*% input$metric %*% t(input$right[[s]])))
    }
  }
  as.numeric(total)
}

# The deliberately slow court: every trace term is expanded over effect and
# feature indices. This shares neither a matrix contraction nor a lowering
# with the explicit formula above.
oracle_bilinear_direct <- function(operator, left, right, gamma, metric) {
  input <- oracle_bilinear_inputs(operator, left, right, gamma, metric)
  total <- 0
  for (r in seq_along(input$left)) {
    for (s in seq_along(input$right)) {
      for (i in seq_len(nrow(input$operator))) {
        for (j in seq_len(ncol(input$operator))) {
          for (a in seq_len(nrow(input$metric))) {
            for (b in seq_len(ncol(input$metric))) {
              total <- total + input$gamma[r, s] * input$operator[i, j] *
                input$left[[r]][i, a] * input$metric[a, b] *
                input$right[[s]][j, b]
            }
          }
        }
      }
    }
  }
  as.numeric(total)
}

oracle_effective_condition <- function(value,
                                       tolerance = sqrt(.Machine$double.eps)) {
  singular <- svd(value, nu = 0, nv = 0)$d
  if (!length(singular) || max(singular) == 0) return(Inf)
  retained <- singular[singular > tolerance * max(singular)]
  if (!length(retained)) return(Inf)
  max(retained) / min(retained)
}

oracle_numeric_comparison <- function(reference, candidate, matrices = list(),
                                      operations = 1L, atol = 1e-12,
                                      rtol = 64 * .Machine$double.eps) {
  if (!is.numeric(reference) || !is.numeric(candidate) ||
      !identical(length(reference), length(candidate)) ||
      any(!is.finite(reference)) || any(!is.finite(candidate))) {
    stop("oracle comparison requires aligned finite numeric values")
  }
  conditions <- if (length(matrices)) {
    vapply(matrices, oracle_effective_condition, numeric(1))
  } else {
    1
  }
  condition <- max(conditions)
  amplification <- min(condition, 1 / sqrt(.Machine$double.eps))
  scale <- pmax(1, abs(reference))
  bound <- atol + rtol * max(1L, operations) * amplification * scale
  error <- abs(candidate - reference)
  list(
    pass = all(error <= bound),
    error = error,
    bound = bound,
    atol = atol,
    rtol = rtol,
    condition = condition,
    operations = as.integer(operations)
  )
}

oracle_common_geometry <- function(blocks, gamma, metric) {
  q <- nrow(blocks[[1L]])
  geometry <- matrix(0, q, q)
  for (r in seq_along(blocks)) {
    for (s in seq_along(blocks)) {
      geometry <- geometry + gamma[r, s] *
        blocks[[r]] %*% metric %*% t(blocks[[s]])
    }
  }
  geometry
}

oracle_energy <- function(operator, blocks, gamma, metric) {
  oracle_bilinear(operator, blocks, blocks, gamma, metric)
}

oracle_pairs <- function(q) t(utils::combn(seq_len(q), 2L))

oracle_rdm <- function(geometry) {
  pairs <- oracle_pairs(nrow(geometry))
  values <- apply(pairs, 1L, function(pair) {
    contrast <- numeric(nrow(geometry))
    contrast[pair] <- c(1, -1)
    as.numeric(crossprod(contrast, geometry %*% contrast))
  })
  names(values) <- apply(pairs, 1L, paste, collapse = "-")
  values
}

oracle_rdm_adjoint <- function(weights, q) {
  pairs <- oracle_pairs(q)
  operator <- matrix(0, q, q)
  for (edge in seq_len(nrow(pairs))) {
    contrast <- numeric(q)
    contrast[pairs[edge, ]] <- c(1, -1)
    operator <- operator + weights[[edge]] * tcrossprod(contrast)
  }
  operator
}

common_geometry_equivalence_oracle <- function() {
  blocks <- list(
    matrix(c(1, 0, 0, 1, 1, 1), 3L, 2L, byrow = TRUE),
    matrix(c(2, 0, 0, 2, 1, 3), 3L, 2L, byrow = TRUE)
  )
  gamma <- matrix(c(0, 0.5, 0.5, 0), 2L, 2L)
  metric <- diag(c(2, 1))
  geometry <- oracle_common_geometry(blocks, gamma, metric)
  rdm <- oracle_rdm(geometry)

  direct_pairs <- vapply(seq_len(3L), function(edge) {
    pair <- oracle_pairs(3L)[edge, ]
    contrast <- numeric(3L)
    contrast[pair] <- c(1, -1)
    oracle_energy(tcrossprod(contrast), blocks, gamma, metric)
  }, numeric(1))
  names(direct_pairs) <- names(rdm)

  design <- cbind(`(Intercept)` = 1, model = c(1, 0, 1))
  coefficient_map <- solve(crossprod(design), t(design))
  rsa <- as.numeric(coefficient_map %*% rdm)
  names(rsa) <- colnames(design)
  adjoint <- lapply(seq_len(nrow(coefficient_map)), function(term) {
    oracle_rdm_adjoint(coefficient_map[term, ], 3L)
  })
  adjoint_rsa <- vapply(adjoint, oracle_energy, numeric(1),
    blocks = blocks, gamma = gamma, metric = metric)
  names(adjoint_rsa) <- colnames(design)

  centred <- c(1, -2, 1)
  centred_weights <- numeric(3L)
  pairs <- oracle_pairs(3L)
  for (edge in seq_len(nrow(pairs))) {
    centred_weights[[edge]] <- -prod(centred[pairs[edge, ]])
  }
  list(
    blocks = blocks,
    gamma = gamma,
    metric = metric,
    geometry = geometry,
    rdm = rdm,
    direct_pairs = direct_pairs,
    contrast_ab = oracle_energy(tcrossprod(c(1, -1, 0)),
      blocks, gamma, metric),
    rsa_design = design,
    coefficient_map = coefficient_map,
    rsa = rsa,
    rsa_adjoint = adjoint_rsa,
    centred_operator = tcrossprod(centred),
    centred_rdm_adjoint = oracle_rdm_adjoint(centred_weights, 3L)
  )
}

if (sys.nframe() == 0L) {
  court <- common_geometry_equivalence_oracle()
  stopifnot(
    identical(court$geometry, matrix(c(
      4, 0, 3, 0, 2, 2.5, 3, 2.5, 5
    ), 3L, 3L, byrow = TRUE)),
    identical(unname(court$rdm), c(6, 3, 2)),
    identical(court$rdm, court$direct_pairs),
    identical(court$contrast_ab, 6),
    max(abs(court$rsa - c(3, 1))) < 1e-14,
    max(abs(court$rsa - court$rsa_adjoint)) < 1e-14,
    max(abs(court$centred_operator - court$centred_rdm_adjoint)) < 1e-14
  )
  cat("common-geometry-equivalence-oracle PASS\n")
}
