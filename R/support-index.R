# Canonical support topology ------------------------------------------------

.support_cell_keys <- function(cells) {
  if (!is.matrix(cells) || !is.numeric(cells) || ncol(cells) < 1L) {
    .input_error("Cell coordinates must be a numeric matrix.")
  }
  do.call(paste, c(unname(as.data.frame(cells)), sep = "\034"))
}

.support_index_pattern_identity <- function(x) {
  list(
    dim = unname(x@Dim),
    p = unname(x@p),
    i = unname(x@i),
    uplo = x@uplo
  )
}

.support_index_signature <- function(domain, node_ids, ptr, members,
                                     construction) {
  semantic <- list(
    schema_version = 1L,
    domain = domain$signature,
    node_ids = node_ids,
    ptr = ptr,
    members = members,
    construction = construction
  )
  .sha256_signature(semantic)
}

.support_index_summary <- function(x) {
  values <- as.double(x)
  stats::setNames(
    c(min(values), stats::median(values), mean(values), max(values)),
    c("min", "median", "mean", "max")
  )
}

.support_index_cost <- function(domain, node_ids, ptr, members,
                                pair_pattern, construction) {
  sizes <- diff(ptr)
  payload <- list(
    domain = domain,
    node_ids = node_ids,
    ptr = ptr,
    members = members,
    construction = construction
  )
  if (!is.null(pair_pattern)) {
    payload$pair_pattern <- pair_pattern
    diagonal <- as.numeric(Matrix::diag(pair_pattern) != 0)
    degree <- as.numeric(Matrix::colSums(pair_pattern)) - diagonal
    pair_pattern_nnz <- as.double(Matrix::nnzero(pair_pattern))
    pair_pattern_stored_nnz <- as.double(length(pair_pattern@i))
    union_degree <- .support_index_summary(degree)
  } else {
    pair_pattern_nnz <- NA_real_
    pair_pattern_stored_nnz <- NA_real_
    union_degree <- stats::setNames(
      rep(NA_real_, 4L), c("min", "median", "mean", "max")
    )
  }
  structural_bytes <- sum(vapply(payload, function(value) {
    as.double(utils::object.size(value))
  }, numeric(1)))
  structure(list(
    nodes = as.double(length(node_ids)),
    features = as.double(domain$n_features),
    support_memberships = as.double(length(members)),
    support_size = .support_index_summary(sizes),
    pair_pattern_nnz = pair_pattern_nnz,
    pair_pattern_stored_nnz = pair_pattern_stored_nnz,
    union_degree = union_degree,
    diagonal_metric_entries = sum(sizes),
    dense_metric_entries = sum(sizes^2),
    one_dense_factorization_pass_units = sum(sizes^3),
    largest_local_metric_bytes = 8 * max(sizes)^2,
    materialized_dense_metric_bytes = 8 * sum(sizes^2),
    estimated_structural_bytes = structural_bytes
  ), class = "effect_support_cost")
}

# The pair pattern marks feature pairs that share at least one support. The
# sparse product costs about sum(support sizes^2) operations, the dense product
# n_features^3 flops in BLAS; the dense route wins once supports cover more than
# a few percent of the domain and the dense square fits comfortably in memory.
# Both routes yield the identical pattern object, so signatures do not depend on
# the route taken.
.support_pair_pattern <- function(membership) {
  n <- ncol(membership)
  density <- Matrix::nnzero(membership) / (as.double(nrow(membership)) * n)
  if (density > 0.05 && n <= 12000L) {
    dense <- crossprod(as.matrix(membership) * 1)
    pattern <- methods::as(dense, "sparseMatrix")
  } else {
    pattern <- Matrix::crossprod(membership)
  }
  methods::as(methods::as(pattern, "symmetricMatrix"), "nMatrix")
}

.support_index_from_membership <- function(membership, domain, node_ids,
                                           construction,
                                           pair_pattern_route = c("lazy",
                                                                  "eager")) {
  .validate_domain(domain)
  pair_pattern_route <- match.arg(pair_pattern_route)
  if (!inherits(membership, "Matrix") || length(dim(membership)) != 2L ||
      any(dim(membership) < 1L) || ncol(membership) != domain$n_features) {
    .input_error(paste0(
      "Support membership must be a nonempty sparse Matrix with one column ",
      "per feature."
    ))
  }
  if (length(node_ids) != nrow(membership) || anyNA(node_ids) ||
      anyDuplicated(node_ids)) {
    .input_error("Support nodes must uniquely identify every membership row.")
  }
  if (!is.list(construction)) {
    .input_error("Support construction metadata must be a list.")
  }
  membership <- methods::as(membership != 0, "nMatrix")
  if (any(Matrix::rowSums(membership) < 1L)) {
    .input_error("Every support node must contain at least one feature.")
  }
  row_membership <- methods::as(membership, "RsparseMatrix")
  ptr <- as.double(row_membership@p)
  members <- as.integer(row_membership@j + 1L)
  pair_pattern <- if (identical(pair_pattern_route, "eager")) {
    .support_pair_pattern(membership)
  } else {
    NULL
  }
  domain_reference <- domain$reference
  cost <- .support_index_cost(
    domain_reference, node_ids, ptr, members, pair_pattern, construction
  )
  signature <- .support_index_signature(
    domain_reference, node_ids, ptr, members, construction
  )
  structure(list(
    domain = domain_reference,
    node_ids = node_ids,
    ptr = ptr,
    members = members,
    pair_pattern = pair_pattern,
    construction = construction,
    cost = cost,
    signature = signature
  ), class = "effect_support_index", pair_pattern_route = pair_pattern_route)
}

.support_index_pair_pattern_route <- function(index) {
  route <- attr(index, "pair_pattern_route", exact = TRUE)
  if (is.null(route)) {
    if (is.null(index$pair_pattern)) "lazy" else "eager"
  } else {
    route
  }
}

.support_index_materialize_pair_pattern <- function(index) {
  if (!is.null(index$pair_pattern)) {
    return(index)
  }
  index <- .validate_support_index(index)
  membership <- .support_index_membership(index)
  pair_pattern <- .support_pair_pattern(membership)
  index$pair_pattern <- pair_pattern
  index$cost <- .support_index_cost(
    index$domain, index$node_ids, index$ptr, index$members,
    pair_pattern, index$construction
  )
  index
}

.support_index_pair_pattern <- function(index) {
  .support_index_materialize_pair_pattern(index)$pair_pattern
}

.support_index_from_members <- function(members, domain, node_ids,
                                        construction,
                                        pair_pattern_route = c("lazy",
                                                               "eager")) {
  .validate_domain(domain)
  if (!is.list(members) || length(members) < 1L) {
    .input_error("Support members must be a nonempty list.")
  }
  members <- lapply(members, function(value) {
    if (!.is_finite_numeric(value) || length(value) < 1L || anyNA(value) ||
        any(value %% 1 != 0) || any(value < 1L) || any(value > domain$n_features)) {
      .input_error("Every support must contain valid finite feature positions.")
    }
    sort(unique(as.integer(value)))
  })
  sizes <- lengths(members)
  membership <- Matrix::sparseMatrix(
    i = rep.int(seq_along(members), sizes),
    j = unlist(members, use.names = FALSE),
    dims = c(length(members), domain$n_features)
  )
  .support_index_from_membership(
    membership, domain, node_ids, construction,
    pair_pattern_route = pair_pattern_route
  )
}

.volume_ball_membership <- function(domain, radius) {
  .validate_domain(domain)
  if (!identical(domain$kind, "volume")) {
    .input_error("A volume stencil requires a volume domain.")
  }
  dims <- domain$metadata$dim
  spacing <- domain$metadata$spacing
  voxel <- domain$metadata$voxel
  if (!is.integer(dims) || length(dims) != 3L ||
      !is.numeric(spacing) || length(spacing) != 3L ||
      !is.matrix(voxel) || !identical(dim(voxel), c(domain$n_features, 3L))) {
    .input_error("The volume domain lacks canonical grid metadata.")
  }
  # An offset can only pair two grid voxels when it is shorter than the grid
  # itself, so the stencil never needs to reach past `dims - 1` per axis. Without
  # this clamp the enumeration grows cubically in the radius regardless of the
  # volume, and a large radius stalls on work that maps nowhere.
  extent <- pmin(floor(radius / spacing + 1e-12), as.double(dims) - 1)
  stencil_size <- prod(2 * extent + 1)
  if (stencil_size > domain$n_features / 4) {
    # The ball is large relative to the volume: one vectorized pass per stencil
    # offset costs about four direct pairwise tests, so past this point the
    # pairwise route is cheaper. Membership, and therefore the support-index
    # signature, is identical.
    return(.volume_ball_membership_pairwise(voxel, spacing, radius))
  }
  offsets <- as.matrix(expand.grid(lapply(extent, function(value) {
    seq.int(-value, value)
  }), KEEP.OUT.ATTRS = FALSE))
  offsets <- offsets[
    rowSums(sweep(offsets, 2L, spacing, `*`)^2) <=
      radius^2 * (1 + 1e-12),
    , drop = FALSE
  ]
  full_size <- prod(as.double(dims))
  if (!is.finite(full_size) || full_size > .Machine$integer.max) {
    .input_error("The volume grid is too large for the compact stencil index.")
  }
  lookup <- integer(as.integer(full_size))
  lookup[domain$feature_ids] <- seq_len(domain$n_features)
  center_chunks <- vector("list", nrow(offsets))
  member_chunks <- vector("list", nrow(offsets))
  strides <- c(1, dims[[1L]], dims[[1L]] * dims[[2L]])
  for (offset_index in seq_len(nrow(offsets))) {
    shifted <- sweep(voxel, 2L, offsets[offset_index, ], `+`)
    valid <- rowSums(shifted < 1L) == 0L &
      rowSums(sweep(shifted, 2L, dims, `>`)) == 0L
    centers <- which(valid)
    if (!length(centers)) next
    full_indices <- 1 + drop((shifted[centers, , drop = FALSE] - 1) %*% strides)
    mapped <- lookup[full_indices]
    present <- mapped > 0L
    center_chunks[[offset_index]] <- as.integer(centers[present])
    member_chunks[[offset_index]] <- as.integer(mapped[present])
  }
  Matrix::sparseMatrix(
    i = unlist(center_chunks, use.names = FALSE),
    j = unlist(member_chunks, use.names = FALSE),
    dims = c(domain$n_features, domain$n_features)
  )
}

# Direct pairwise membership on the spacing-scaled voxel grid. Used when the
# offset stencil would be larger than the volume it indexes; produces exactly
# the membership the stencil would.
.volume_ball_membership_pairwise <- function(voxel, spacing, radius) {
  n <- nrow(voxel)
  limit <- radius^2 * (1 + 1e-12)
  block <- max(1L, min(n, as.integer(floor(2^22 / max(n, 1L)))))
  starts <- seq.int(1L, n, by = block)
  center_chunks <- vector("list", length(starts))
  member_chunks <- vector("list", length(starts))
  for (chunk in seq_along(starts)) {
    rows <- seq.int(starts[[chunk]], min(starts[[chunk]] + block - 1L, n))
    # Integer grid differences times spacing, squared and summed per axis:
    # the same arithmetic as the stencil test, so boundary membership agrees.
    squared <- 0
    for (axis in seq_len(ncol(voxel))) {
      squared <- squared +
        (outer(voxel[rows, axis], voxel[, axis], `-`) * spacing[[axis]])^2
    }
    hits <- which(squared <= limit, arr.ind = TRUE, useNames = FALSE)
    center_chunks[[chunk]] <- as.integer(rows[hits[, 1L]])
    member_chunks[[chunk]] <- as.integer(hits[, 2L])
  }
  Matrix::sparseMatrix(
    i = unlist(center_chunks, use.names = FALSE),
    j = unlist(member_chunks, use.names = FALSE),
    dims = c(n, n)
  )
}

.coordinate_ball_members <- function(domain, radius) {
  .validate_domain(domain)
  coordinates <- domain$coordinates
  if (is.null(coordinates)) {
    .input_error(sprintf(paste0(
      "Searchlights need feature coordinates, and domain `%s` has none. ",
      "Build it with `abstract_domain(n, coordinates = )`, `volume_domain()`, ",
      "or `neuroim2_volume_domain()`; `regions()` and `whole_brain()` need no ",
      "geometry."
    ), domain$id))
  }
  dimensions <- ncol(coordinates)
  neighbor_cells <- 3^dimensions
  if (!is.finite(neighbor_cells) || neighbor_cells > 100000L) {
    .input_error(paste0(
      "High-dimensional coordinates require a domain-native neighborhood ",
      "provider."
    ))
  }
  cells <- floor(coordinates / radius)
  offsets <- as.matrix(expand.grid(
    rep(list(-1:1), dimensions), KEEP.OUT.ATTRS = FALSE
  ))
  cell_min <- apply(cells, 2L, min)
  cell_max <- apply(cells, 2L, max)
  radix <- cell_max - cell_min + 3
  code_space <- prod(radix)
  exact_mixed_radix <- is.finite(code_space) && code_space <= 2^52

  if (exact_mixed_radix) {
    strides <- c(1, cumprod(radix[-length(radix)]))
    shifted <- sweep(cells, 2L, cell_min - 1, `-`)
    cell_codes <- 1 + drop(shifted %*% strides)
    occupied_codes <- sort(unique(cell_codes))
    groups <- split(seq_len(domain$n_features),
      match(cell_codes, occupied_codes), drop = TRUE)
    offset_delta <- drop(offsets %*% strides)
    candidate_groups <- matrix(
      match(c(outer(cell_codes, offset_delta, `+`)), occupied_codes,
        nomatch = 0L),
      nrow = domain$n_features,
      ncol = nrow(offsets)
    )
    lookup <- function(center) {
      ids <- candidate_groups[center, ]
      unlist(groups[ids[ids > 0L]], use.names = FALSE)
    }
    lookup_kind <- "mixed_radix_vectorized"
  } else {
    keys <- .support_cell_keys(cells)
    groups <- split(seq_len(domain$n_features), keys, drop = TRUE)
    table <- new.env(
      hash = TRUE, parent = emptyenv(), size = max(29L, length(groups))
    )
    list2env(groups, envir = table)
    lookup <- function(center) {
      candidate_keys <- .support_cell_keys(
        sweep(offsets, 2L, cells[center, ], `+`)
      )
      unlist(mget(candidate_keys, envir = table, inherits = FALSE,
        ifnotfound = list(integer())), use.names = FALSE)
    }
    lookup_kind <- "hashed_tuple"
  }
  members <- vector("list", domain$n_features)
  for (center in seq_len(domain$n_features)) {
    candidates <- sort(unique(as.integer(lookup(center))))
    squared <- rowSums((coordinates[candidates, , drop = FALSE] -
      matrix(coordinates[center, ], nrow = length(candidates),
        ncol = dimensions, byrow = TRUE))^2)
    members[[center]] <- candidates[
      squared <= radius^2 * (1 + 1e-12)
    ]
  }
  structure(members, lookup = lookup_kind)
}

.euclidean_support_index <- function(domain, radius,
                                     provider = c("auto", "volume_stencil",
                                                  "coordinate_cells")) {
  .validate_domain(domain)
  if (!.is_number(radius) || radius <= 0) {
    .input_error("Support radius must be one positive finite number.")
  }
  provider <- match.arg(provider)
  if (provider == "auto") {
    provider <- if (identical(domain$kind, "volume")) {
      "volume_stencil"
    } else {
      "coordinate_cells"
    }
  }
  construction <- list(
    kind = "euclidean_ball",
    provider = provider,
    radius = as.numeric(radius),
    coordinate_units = domain$coordinate_units
  )
  if (provider == "volume_stencil") {
    construction$lookup <- "offset_stencil"
    membership <- .volume_ball_membership(domain, radius)
    return(.support_index_from_membership(
      membership, domain, domain$feature_ids, construction
    ))
  }
  members <- .coordinate_ball_members(domain, radius)
  construction$lookup <- attr(members, "lookup", exact = TRUE)
  .support_index_from_members(
    members, domain, domain$feature_ids, construction
  )
}

.support_index_membership <- function(index) {
  .validate_support_index(index)
  sizes <- diff(index$ptr)
  Matrix::sparseMatrix(
    i = rep.int(seq_along(index$node_ids), sizes),
    j = index$members,
    dims = c(length(index$node_ids), index$domain$n_features)
  )
}

.support_index_support <- function(index, nodes) {
  .validate_support_index(index)
  if (!.is_finite_numeric(nodes) || length(nodes) < 1L || anyNA(nodes) ||
      any(nodes %% 1 != 0) || any(nodes < 1L) ||
      any(nodes > length(index$node_ids))) {
    .input_error("Support nodes must be valid integer positions.")
  }
  lapply(as.integer(nodes), function(node) {
    start <- index$ptr[[node]] + 1
    end <- index$ptr[[node + 1L]]
    index$members[seq.int(start, end)]
  })
}

# Trusted hot-loop accessor. The containing plan has already validated the
# complete support index, so repeating global domain, CSR, and signature scans
# for every node would turn one linear execution into quadratic validation.
.support_index_support_trusted <- function(index, node) {
  start <- index$ptr[[node]] + 1
  end <- index$ptr[[node + 1L]]
  index$members[seq.int(start, end)]
}

.support_index_blocks <- function(index, block_size) {
  .validate_support_index(index)
  block_size <- .validate_tile_size(block_size, "block_size")
  starts <- seq.int(1L, length(index$node_ids), by = block_size)
  data.frame(
    block = seq_along(starts),
    start = starts,
    end = pmin(starts + block_size - 1L, length(index$node_ids)),
    stringsAsFactors = FALSE
  )
}

.support_index_gather <- function(index, x, node, margin = 2L) {
  support <- .support_index_support(index, node)[[1L]]
  if (is.atomic(x) && is.null(dim(x))) {
    if (length(x) != index$domain$n_features) {
      .contract_error(
        "A gathered vector must match the support feature domain."
      )
    }
    return(x[support])
  }
  if (!is.matrix(x) || !is.numeric(x) ||
      !is.numeric(margin) || length(margin) != 1L ||
      !margin %in% c(1L, 2L)) {
    .input_error(
      "Support gathering requires a numeric vector or matrix and margin 1 or 2."
    )
  }
  if (margin == 1L) {
    if (nrow(x) != index$domain$n_features) {
      .contract_error("Matrix rows do not match the support feature domain.")
    }
    return(x[support, , drop = FALSE])
  }
  if (ncol(x) != index$domain$n_features) {
    .contract_error("Matrix columns do not match the support feature domain.")
  }
  x[, support, drop = FALSE]
}

.support_index_preflight <- function(index, workspace_bytes = 4 * 1024^3,
                                     evaluation_edges = 1L) {
  .validate_support_index(index)
  if (!.is_number(workspace_bytes) || workspace_bytes <= 0) {
    .input_error("Support preflight requires one positive finite byte budget.")
  }
  if (!.is_number(evaluation_edges) || evaluation_edges < 1 ||
      evaluation_edges %% 1 != 0) {
    .input_error(
      "Support preflight requires a positive integer evaluation-edge count."
    )
  }
  cost <- index$cost
  evaluation_edges <- as.double(evaluation_edges)
  dense_bytes_per_edge <- cost$materialized_dense_metric_bytes
  dense_bytes_total <- dense_bytes_per_edge * evaluation_edges
  factorization_units_per_edge <-
    cost$one_dense_factorization_pass_units
  factorization_units_total <-
    factorization_units_per_edge * evaluation_edges
  structure(list(
    workspace_bytes = as.double(workspace_bytes),
    evaluation_edges = evaluation_edges,
    structural_fits = cost$estimated_structural_bytes <= workspace_bytes,
    materialized_dense_metric_bytes_per_edge = dense_bytes_per_edge,
    materialized_dense_metric_bytes_total = dense_bytes_total,
    materialized_dense_metric_per_edge_fits =
      dense_bytes_per_edge <= workspace_bytes,
    materialized_dense_metric_fits = dense_bytes_total <= workspace_bytes,
    largest_local_metric_fits =
      cost$largest_local_metric_bytes <= workspace_bytes,
    one_dense_factorization_pass_units_per_edge =
      factorization_units_per_edge,
    one_dense_factorization_pass_units_total =
      factorization_units_total,
    factorization_passes = "metric_capability_dependent",
    dense_metric_lowering = if (
      dense_bytes_total <= workspace_bytes
    ) "materializable" else "support_streamed_only",
    metric_schedule_storage = "derive_on_demand_from_frozen_recipe",
    cost = cost
  ), class = "effect_support_preflight")
}

.validate_support_index <- function(index, deep = FALSE) {
  expected <- c("domain", "node_ids", "ptr", "members", "pair_pattern",
    "construction", "cost", "signature")
  if (!.sealed_fields(index, "effect_support_index", expected)) {
    .input_error("Support-index fields are missing or noncanonical.")
  }
  domain <- .validate_domain_reference(index$domain)
  nodes <- length(index$node_ids)
  if (nodes < 1L || anyNA(index$node_ids) || anyDuplicated(index$node_ids) ||
      !is.double(index$ptr) || length(index$ptr) != nodes + 1L ||
      index$ptr[[1L]] != 0 || any(diff(index$ptr) < 1) ||
      utils::tail(index$ptr, 1L) != length(index$members) ||
      !is.integer(index$members) || any(index$members < 1L) ||
      any(index$members > domain$n_features)) {
    .input_error("Support-index nodes or CSR membership are invalid.")
  }
  if (!is.null(index$pair_pattern) &&
      (!inherits(index$pair_pattern, "symmetricMatrix") ||
       !inherits(index$pair_pattern, "nMatrix") ||
       !identical(as.integer(dim(index$pair_pattern)),
         c(domain$n_features, domain$n_features)))) {
    .input_error("Support-index pair pattern is invalid.")
  }
  if (!is.list(index$construction) ||
      !inherits(index$cost, "effect_support_cost") ||
      !.strong_sha256(index$signature)) {
    .input_error("Support-index identity or cost metadata are invalid.")
  }
  if (isTRUE(deep)) {
    sizes <- diff(index$ptr)
    starts <- index$ptr[-length(index$ptr)] + 1
    canonical <- vapply(seq_len(nodes), function(node) {
      values <- index$members[seq.int(starts[[node]],
        starts[[node]] + sizes[[node]] - 1)]
      length(values) == 1L || all(diff(values) > 0L)
    }, logical(1))
    if (!all(canonical)) {
      .input_error(
        "Support-index CSR rows must be strictly increasing and unique."
      )
    }
    if (!is.null(index$pair_pattern)) {
      membership <- .support_index_membership(index)
      expected_pattern <- methods::as(
        Matrix::crossprod(membership), "nMatrix"
      )
      if (!identical(index$pair_pattern, expected_pattern)) {
        .contract_error(
          "Support-index pair pattern does not match its supports."
        )
      }
    }
    expected_cost <- .support_index_cost(
      domain, index$node_ids, index$ptr, index$members,
      index$pair_pattern, index$construction
    )
    expected_signature <- .support_index_signature(
      domain, index$node_ids, index$ptr, index$members,
      index$construction
    )
    if (!identical(index$cost, expected_cost) ||
        !identical(index$signature, expected_signature)) {
      .contract_error("Support-index cached cost or identity is inconsistent.")
    }
  }
  invisible(index)
}
