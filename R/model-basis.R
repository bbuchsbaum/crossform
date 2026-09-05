# Model bases: the model family as a declared basis of the effect axis --------
#
# Layer 2 (values), beside `effect-map.R`, `extractor.R` and `metric.R`; it
# calls sideways into those three and downward into layer 1, and nothing
# below it calls back. `design/predictive-geometry-contract.md` and its
# migration table govern the current admission; the older model-coordinate
# contract records the descriptive reader's original derivation.
#
# A model basis is one data-independent value: an orthonormal, centered basis
# `Q` of the joint span of a family of model geometries, the factor `R` that
# recovers each model's own coordinates (`T = Q R`), the per-model spectra,
# the requested and effective ranks, the negative model mass dropped, the
# span fraction, and a signature. It is built before any neural data is read.
# Lowering a relation through it (`B~ = Q' B`) yields the complete geometry
# in model coordinates, `S_x = Q' G_x Q`, on which the existing readers
# apply.
#
# Four construction rules are the reason this file exists rather than a
# five-line helper in a vignette, and each was measured before it was written:
#
#  1. EVERY TRUNCATION IS RANK-REVEALING. The literal factor
#     `U_r Lambda_r^{1/2}` with `r` above a model's true rank has a zero
#     column, and a thin unpivoted QR of it returns a second "basis vector"
#     that is not centered (`|Q'1| = 1.10` on the oracle fixture, `1.15` on
#     the README fixture); under a condition-independent baseline shift the
#     compressed form then moves by `8.12`. So each model's Gram is truncated
#     at a relative eigenvalue tolerance, the requested rank is clamped to the
#     effective rank and the clamp is recorded, and the stacked factor is
#     factored by a thin SVD at the same tolerance.
#
#  2. EVERY DECOMPOSITION RUNS IN THE CENTERED SUBSPACE. The Gram `-HDH/2` is
#     only centered to round-off, and the eigenvector of a small kept root
#     `lambda_k` picks up a component along `1` of order
#     `eps * lambda_max / lambda_k`: on a two-feature model whose second
#     feature is scaled by `1e-4`, a direct eigendecomposition left
#     `|Q'1| = 1.1e-6` (fresh-context review of M2, 2026-09-04). So every
#     Gram is formed and decomposed in an explicit orthonormal basis `V` of
#     `1-perp`, and `Q = V U~` is centered to `eps` whatever the conditioning.
#
#  3. OVERLAP IS ADMITTED FOR FORM PREDICTION. Rectangular R retains each
#     model's factor while Q spans their union once. Unique coefficients in
#     redundant model coordinates are a stronger requirement, checked by
#     readers that report them. Pooling kernels never requires solve(R).
#
#  4. CENTERING IS A PROPERTY OF THE BASIS. Every operator the basis induces
#     is `sym(q_j q_k')`, whose marginals vanish iff `Q'1 = 0`; so baseline
#     invariance of the whole lowered form is one assertion on `Q`, made here
#     and recorded as `baseline_invariant`, rather than a diagnosis of
#     `q(q+1)/2` operators downstream. By rule 2 a correctly built basis
#     cannot fail it, which is why its failure is the package's invariant
#     error and not a refusal.
#
# What the basis does NOT do: it never subtracts a noise floor, never clips
# the neural side, and is not a neural estimate. Its model spectra and factors
# retain the geometry used by `fit_geometry()`; the descriptive PSD projection
# remains a separate reading in `latent_geometry()`.

.model_basis_distances <- "squared_euclidean"
.model_basis_coordinate_prefix <- "mc"

# Flip each column so that its largest-magnitude entry (the first, on ties)
# is positive. Eigenvectors and singular vectors are defined up to sign, and a
# signature over `Q` and `R` should not depend on which sign LAPACK chose.
.model_basis_fix_signs <- function(x) {
  for (j in seq_len(ncol(x))) {
    pivot <- which.max(abs(x[, j]))
    if (x[pivot, j] < 0) x[, j] <- -x[, j]
  }
  x
}

# An orthonormal basis of 1-perp: the Householder QR of `[1, I]` puts `1`
# first and the remaining columns orthogonal to it, deterministically.
.model_basis_centered_frame <- function(n) {
  V <- qr.Q(qr(cbind(1, diag(n))))[, -1L, drop = FALSE]
  .model_basis_fix_signs(V)
}

# One named list of matrices, or a single matrix wrapped under `default_name`.
.model_basis_entries <- function(value, label, default_name) {
  if (is.null(value)) return(list())
  if (is.matrix(value)) value <- stats::setNames(list(value), default_name)
  if (!is.list(value) || !length(value)) {
    .input_error(sprintf(paste0(
      "`%s` must be one matrix or a nonempty named list of them; received %s."
    ), label, .msg_value(value)),
      arg = label, received = .msg_value(value),
      expected = "one matrix, or a named list of matrices")
  }
  if (!.is_strings(names(value), unique = TRUE) ||
      length(names(value)) != length(value)) {
    .input_error(sprintf(paste0(
      "`%s` must be a list with unique nonempty names; the names become the ",
      "model labels of the basis."
    ), label),
      arg = label, received = .msg_names(names(value)),
      expected = "a list with unique nonempty names")
  }
  for (entry in names(value)) {
    x <- value[[entry]]
    where <- sprintf("`%s` entry `%s`", label, entry)
    if (!is.matrix(x) || !is.numeric(x)) {
      .input_error(sprintf("%s must be a numeric matrix; received %s.",
        where, .msg_value(x)),
        arg = label, received = .msg_value(x), expected = "a numeric matrix")
    }
    if (any(!is.finite(x))) {
      .input_error(sprintf("%s contains %s non-finite %s.", where,
        sum(!is.finite(x)),
        if (sum(!is.finite(x)) == 1L) "entry" else "entries"),
        arg = label,
        received = sprintf("%d non-finite of %d", sum(!is.finite(x)),
          length(x)),
        expected = "all entries finite")
    }
    if (any(dim(x) < 1L)) {
      .input_error(sprintf("%s is empty.", where), arg = label,
        received = sprintf("%d x %d", nrow(x), ncol(x)),
        expected = "a nonempty matrix")
    }
  }
  value
}

# The condition axis the basis is declared over: an `effect_space()`, a
# `condition_space()`, condition names, or NULL (then read from the labels the
# first labelled model carries). Returns the ordered names, the shared unit
# and scale, where the names came from, and the condition space when one was
# supplied.
.model_basis_conditions <- function(conditions, entries) {
  labels_of <- function(x, kind) {
    if (kind != "features" && is.null(rownames(x)) && is.null(colnames(x))) {
      return(NULL)
    }
    rownames(x)
  }
  shared_of <- function(values, what) {
    shared <- unique(unname(values))
    if (length(shared) != 1L) {
      .input_error(sprintf(paste0(
        "The condition coordinates carry different %ss (%s), so no unit-norm ",
        "combination of them is well typed. Declare one shared %s on the ",
        "space before building a model basis."
      ), what, .msg_names(as.character(shared)), what),
        arg = "conditions", received = .msg_names(as.character(shared)),
        expected = sprintf("one shared %s", what))
    }
    shared
  }
  units <- "arbitrary"
  scale <- 1
  space <- NULL
  declaration <- list(kind = "names", signature = NULL)
  origin <- "`conditions`"
  if (inherits(conditions, "effect_condition_space")) {
    space <- .validate_condition_space(conditions)
    declaration <- list(kind = "condition_space", signature = space$signature)
    names <- space$coordinates
    units <- shared_of(space$units, "unit")
    scale <- shared_of(space$scale, "scale")
  } else if (inherits(conditions, "effect_space")) {
    reference <- .validate_effect_space(conditions)
    declaration <- list(kind = "effect_space", signature = reference$signature)
    names <- reference$coordinates
    units <- shared_of(reference$units, "unit")
    scale <- shared_of(reference$scale, "scale")
  } else if (is.null(conditions)) {
    labelled <- Filter(Negate(is.null), Map(labels_of,
      lapply(entries, `[[`, "value"), vapply(entries, `[[`, "", "kind")))
    if (!length(labelled)) {
      .input_error(paste0(
        "Supply `conditions` (condition names, an `effect_space()`, or a ",
        "`condition_space()`), or label every model's condition axis: an ",
        "unlabelled model cannot be aligned to the relation's effects."
      ), arg = "conditions", received = "NULL and unlabelled models",
        expected = "condition names or labelled models")
    }
    names <- labelled[[1L]]
    origin <- sprintf("the labels of model `%s`, the first labelled model",
      names(labelled)[[1L]])
  } else {
    if (!is.character(conditions) || length(conditions) < 2L) {
      .input_error(sprintf(paste0(
        "`conditions` must name at least two conditions, or be an ",
        "`effect_space()` or `condition_space()`; received %s."
      ), .msg_value(conditions)),
        arg = "conditions", received = .msg_value(conditions),
        expected = "at least two condition names")
    }
    names <- conditions
  }
  names <- .validate_effect_names(names, length(names))
  if (length(names) < 2L) {
    .input_error(sprintf(paste0(
      "A model basis needs at least two conditions: with %s there is no ",
      "centered direction to span."
    ), .msg_count(length(names), "condition")),
      arg = "conditions", received = .msg_count(length(names), "condition"),
      expected = "at least two")
  }
  list(names = names, units = units, scale = scale, space = space,
    origin = origin, declaration = declaration)
}

# Align one model to the declared condition order, as `rsa()` aligns a model
# RDM: both axes named with the declared conditions (any order), or neither
# axis named and the dimension equal to the condition count. The numeric
# checks on an RDM run AFTER alignment, because a matrix whose rows and
# columns carry the same labels in different orders is symmetric only once
# both axes follow one order.
.model_basis_align <- function(x, kind, conditions, where, origin) {
  n <- length(conditions)
  square <- kind != "features"
  axes <- if (square) list(rownames(x), colnames(x)) else {
    list(rownames(x))
  }
  named <- !vapply(axes, is.null, logical(1))
  if (square) {
    if (nrow(x) != ncol(x)) {
      .input_error(sprintf(
        "%s is %d x %d; a model kernel or RDM must be square.", where, nrow(x), ncol(x)),
        received = sprintf("%d x %d", nrow(x), ncol(x)), expected = "square")
    }
    if (xor(named[[1L]], named[[2L]])) {
      .input_error(sprintf(paste0(
        "%s names only its %s; name both axes with the conditions (%s), or ",
        "neither."
      ), where, if (named[[1L]]) "rows" else "columns", .msg_names(conditions)),
        received = if (named[[1L]]) "row names only" else "column names only",
        expected = "both axes named, or neither")
    }
  }
  if (!any(named)) {
    if (nrow(x) != n) {
      .input_error(sprintf(paste0(
        "%s is unlabelled and has %s, but the basis declares %s (%s, from ",
        "%s). An unlabelled model is admitted only when its dimension equals ",
        "the condition count; otherwise name its axes."
      ), where, .msg_count(nrow(x), "row"), .msg_count(n, "condition"),
        .msg_names(conditions), origin),
        received = .msg_count(nrow(x), "row"),
        expected = .msg_count(n, "condition"))
    }
    dimnames(x) <- if (square) list(conditions, conditions) else {
      list(conditions, colnames(x))
    }
  } else {
    for (labels in axes) {
      if (anyNA(labels) || anyDuplicated(labels) ||
          !setequal(labels, conditions)) {
        unknown <- setdiff(unique(labels), conditions)
        absent <- setdiff(conditions, labels)
        detail <- c(
          if (length(unknown)) sprintf("%s is not a declared condition",
            .msg_names(unknown)),
          if (length(absent)) sprintf("%s is missing from an axis",
            .msg_names(absent))
        )
        if (!length(detail)) detail <- "an axis repeats a condition"
        .input_error(sprintf(
          "%s axis names do not match the conditions (%s, from %s): %s.",
          where, .msg_names(conditions), origin,
          paste(detail, collapse = "; ")),
          received = paste(detail, collapse = "; "),
          expected = .msg_names(conditions))
      }
    }
    x <- if (square) {
      x[conditions, conditions, drop = FALSE]
    } else {
      x[conditions, , drop = FALSE]
    }
  }
  if (square) {
    reference <- max(abs(x))
    asymmetry <- max(abs(x - t(x)))
    if (asymmetry > 1e-12 * reference) {
      .input_error(sprintf(paste0(
        "%s is not symmetric once both axes follow the condition order; the ",
        "largest difference between `m[i, j]` and `m[j, i]` is %g."
      ), where, asymmetry),
        received = sprintf("asymmetry %g", asymmetry),
        expected = "a symmetric matrix")
    }
    if (kind == "rdm" && max(abs(diag(x))) > 1e-12 * reference) {
      .input_error(sprintf(paste0(
        "%s has a nonzero diagonal (largest |m[i, i]| is %g); a condition is ",
        "at squared distance zero from itself."
      ), where, max(abs(diag(x)))),
        received = sprintf("largest |m[i, i]| is %g", max(abs(diag(x)))),
        expected = "a zero diagonal")
    }
    if (kind == "rdm" && any(x < 0)) {
      .input_error(sprintf(
        "%s has negative entries; a squared Euclidean distance is nonnegative.",
        where), received = sprintf("min %g", min(x)),
        expected = "nonnegative entries")
    }
  }
  x
}

# Admission before centering: a negative constant direction cannot be hidden
# by restricting a malformed kernel to 1-perp. An explicitly allowed repair
# projects the input kernel first, with its cost retained separately from the
# later centered-model spectrum.
.model_basis_admit_kernel <- function(x, name, tolerance, negative_share) {
  e <- eigen(x/2 + t(x)/2, symmetric = TRUE)
  scale <- max(abs(e$values))
  if (!is.finite(scale)) {
    .input_error(sprintf("Model kernel `%s` overflows its spectrum; rescale its input.", name))
  }
  positive <- e$values > tolerance * scale
  negative <- e$values < -tolerance * scale
  mass <- -sum(e$values[negative])
  numerical <- -sum(e$values[e$values < 0 & !negative])
  share <- if (sum(e$values[positive]) > 0) {
    mass / sum(e$values[positive])
  } else if (mass > 0) Inf else 0
  if (share > negative_share) {
    .capability_refusal(sprintf(
      "Model kernel `%s` is not positive semidefinite before centering (negative mass %.6g).",
      name, mass), capability = "euclidean_model_geometry",
      namespace = "model_coordinate", reasons = "model_gram_indefinite",
      remedies = c("Supply a positive semidefinite kernel or its generating features.",
        "Declare a larger `negative_share` only to admit an explicit positive-part model replacement."))
  }
  if (any(negative)) {
    repaired <- e$vectors %*% (pmax(e$values, 0) * t(e$vectors))
    dimnames(repaired) <- dimnames(x)
    x <- repaired
  }
  list(value = x, negative_mass = mass, numerical_negative_mass = numerical)
}

# The centered Gram of one model, formed directly in the centered frame `V`:
# `-V' D V / 2` for a squared Euclidean RDM, `(V'F)(V'F)'` for features. The
# two agree exactly when `D` is the squared Euclidean distance of the rows
# of `F`, and both equal `V' K V` for the full-space Gram `K`.
.model_basis_gram <- function(x, kind, V) {
  if (kind == "rdm") {
    -0.5 * t(V) %*% x %*% V
  } else if (kind == "kernel") {
    # Subtract the common offsets before multiplying by V. This avoids
    # roundoff in V'1 magnifying a large condition-independent baseline.
    centered <- sweep(x, 1L, rowMeans(x), "-")
    centered <- sweep(centered, 2L, colMeans(centered), "-")
    t(V) %*% centered %*% V
  } else {
    tcrossprod(t(V) %*% sweep(x, 2L, colMeans(x), "-"))
  }
}

# The rank-revealing factor `T_i = U Lambda^{1/2}` of one model's Gram, with
# `U = V U~` lifted from the centered frame. Roots within `tolerance` of zero
# relative to the largest absolute root are not model geometry on either
# side: positive ones are not kept, negative ones are not counted as dropped
# mass, so a Euclidean model records a dropped mass of exactly zero. The
# requested rank is clamped to what survives and the clamp is reported.
.model_basis_factor <- function(K, V, x, name, kind, rank, tolerance,
                                negative_share, normalize) {
  n <- nrow(V)
  if (!.is_finite_matrix(K)) {
    .input_error(sprintf("Model `%s` overflows while forming its Gram; rescale its input.", name))
  }
  e <- eigen(K, symmetric = TRUE)
  scale <- max(abs(e$values))
  if (!is.finite(scale) || scale <= 0) {
    .input_error(sprintf(paste0(
      "Model `%s` has no geometry: %s, so it distinguishes no conditions."
    ), name, if (kind == "rdm") {
      "every squared distance is zero"
    } else if (kind == "features") {
      "every feature is constant across conditions"
    } else {
      "the centered kernel is zero"
    }), received = "a model with no centered variation",
      expected = "at least one positive Gram root")
  }
  positive <- e$values > tolerance * scale
  negative <- e$values < -tolerance * scale
  positive_mass <- sum(e$values[positive])
  dropped <- -sum(e$values[negative])
  share <- if (positive_mass > 0) dropped / positive_mass else Inf
  if (share > negative_share) {
    .capability_refusal(
      sprintf(paste0(
        "Model `%s` is not a Euclidean geometry: its centered Gram carries ",
        "negative mass %.3g against positive mass %.3g (share %.3g, above ",
        "the declared limit %.3g), so no model factor `T` reproduces it."
      ), name, dropped, positive_mass, share, negative_share),
      capability = "euclidean_model_geometry",
      namespace = "model_coordinate",
      reasons = c(
        "model_gram_indefinite",
        sprintf("Model `%s`: negative Gram mass %.6g of positive mass %.6g.",
          name, dropped, positive_mass)
      ),
      remedies = c(
        paste0("Supply a Euclidean embedding of the model (for example a ",
          "classical multidimensional scaling of it) as `features`."),
        "Supply the feature matrix that generated the dissimilarities.",
        sprintf(paste0("Raise `negative_share` above %.3g to admit the ",
          "positive part alone, knowing the dropped mass is recorded."), share)
      )
    )
  }
  n_positive <- sum(positive)
  if (n_positive == 0L || !is.finite(positive_mass)) {
    .input_error(sprintf("Model `%s` has no finite positive centered geometry; rescale or replace it.", name))
  }
  effective <- if (is.na(rank)) n_positive else min(rank, n_positive)
  .model_basis_rank_cut(e$values, effective, tolerance, sprintf("model `%s`", name))
  keep <- seq_len(effective)
  retained_trace <- sum(e$values[keep])
  if (normalize == "trace" && !is.finite(1 / retained_trace)) {
    .input_error(sprintf("Model `%s` has an unrepresentable normalization scale; rescale its input.", name))
  }
  U <- .model_basis_fix_signs(V %*% e$vectors[, keep, drop = FALSE])
  factor <- U %*% diag(sqrt(e$values[keep]), effective)
  if (normalize == "trace") factor <- factor / sqrt(retained_trace)
  list(
    factor = factor,
    record = data.frame(
      model = name,
      kind = kind,
      rank_requested = as.integer(rank),
      rank_effective = as.integer(effective),
      rank_clamped = !is.na(rank) && rank > effective,
      positive_mass = positive_mass,
      raw_trace = sum(e$values),
      retained_trace = retained_trace,
      normalization_scale = if (normalize == "trace") 1 / retained_trace else 1,
      truncated_positive_mass = positive_mass - retained_trace,
      numerical_positive_mass = sum(e$values[e$values > 0 & !positive]),
      numerical_negative_mass = -sum(e$values[e$values < 0 & !negative]),
      dropped_negative_mass = dropped,
      negative_share = share,
      stringsAsFactors = FALSE
    ),
    # The full-space spectrum: the centered frame's n - 1 roots plus the one
    # exact zero that centering removes.
    spectrum = c(e$values, 0)
  )
}

# A cut through a positive tie chooses a direction the model did not choose.
# Keeping the whole tie is well defined even though its individual factors
# are free to rotate. The same policy is used by the predictive rank path.
.model_basis_rank_cut <- function(values, rank, tolerance, context) {
  if (rank > 0L && rank < length(values) && values[[rank + 1L]] > 0 &&
      abs(values[[rank]] - values[[rank + 1L]]) <= tolerance * max(abs(values))) {
    .capability_refusal(sprintf(
      "The requested rank of %s splits a positive tied eigenspace; retain or omit the whole group.",
      context), capability = "identified_rank_projection",
      namespace = "model_coordinate", reasons = "ambiguous_rank_cut",
      remedies = "Choose a rank at a gap in the positive spectrum.")
  }
  invisible(rank)
}

.model_basis_ranks <- function(rank, models) {
  k <- length(models)
  if (is.null(rank)) return(rep(NA_integer_, k))
  if (!is.numeric(rank) || !length(rank) %in% c(1L, k) || anyNA(rank) ||
      any(!is.finite(rank)) || any(rank < 1) ||
      any(rank > .Machine$integer.max) || any(rank %% 1 != 0)) {
    .input_error(sprintf(paste0(
      "`rank` must be NULL, one positive whole number, or one per model ",
      "(%s); received %s."
    ), .msg_names(models), .msg_value(rank)),
      arg = "rank", received = .msg_value(rank),
      expected = "NULL, one positive whole number, or one per model")
  }
  if (!is.null(names(rank))) {
    if (!setequal(names(rank), models) || anyDuplicated(names(rank))) {
      .input_error(sprintf(
        "`rank` names (%s) must match the model names (%s).",
        .msg_names(names(rank)), .msg_names(models)),
        arg = "rank", received = .msg_names(names(rank)),
        expected = .msg_names(models))
    }
    rank <- rank[models]
  }
  if (length(rank) == 1L) rank <- rep(rank, k)
  as.integer(unname(rank))
}

# The names of the basis coordinates and of the model coordinates, derived
# from the model names and effective ranks so the validator can rebuild them.
.model_basis_coordinate_names <- function(q) {
  paste0(.model_basis_coordinate_prefix, seq_len(q))
}

.model_basis_model_coordinate_names <- function(models, ranks) {
  unlist(Map(function(nm, r) paste0(nm, ".", seq_len(r)), models, ranks),
    use.names = FALSE)
}

.model_basis_columns <- function(models, ranks) {
  ends <- cumsum(ranks)
  starts <- ends - ranks + 1L
  stats::setNames(Map(function(a, b) seq.int(a, b), starts, ends), models)
}

# The two products of a basis: the extractor that lowers a beta block on the
# betas path, and the effect map that lowers through a study's condition
# space on the ingestion path. Both carry the model-coordinate effect space,
# whose identity carries the basis signature (contract section 2.1).
.model_basis_products <- function(Q, conditions, units, scale, space,
                                  signature, models) {
  q <- ncol(Q)
  coordinates <- .model_basis_coordinate_names(q)
  effects <- effect_space(
    coordinates,
    basis_id = "model-coordinates",
    units = units,
    scale = scale,
    provenance = list(
      model_basis_signature = signature,
      models = models,
      conditions = conditions
    )
  )
  map <- t(Q)
  dimnames(map) <- list(coordinates, conditions)
  extractor <- effect_extractor(
    map,
    effects = effects,
    estimator = "model_basis",
    diagnostics = list(model_basis_signature = signature)
  )
  effect_map <- if (is.null(space)) NULL else {
    effect_map(map, conditions = space, effects = effects,
      provenance = list(model_basis_signature = signature))
  }
  list(effect_space = effects, extractor = extractor, effect_map = effect_map)
}

# Is this effect space one a model basis produced? The lowering entries and
# the evidence-task gate ask this of a relation's effect space, so the answer
# is read from the space's identity and not from any object that happens to
# be nearby: the basis id names the vocabulary, and the provenance carries the
# basis signature that section 2.1(a) puts into every plan downstream. The
# signature is checked for shape, not resolved against a basis (none is at
# hand where this is asked); the answer only ever triggers refusals, so a
# hand-built space that borrows the id and a well-formed digest is treated as
# lowered and refused with it, never admitted by it.
.is_model_coordinate_space <- function(space) {
  if (!inherits(space, "effect_space") ||
      !identical(space$basis_id, "model-coordinates")) {
    return(FALSE)
  }
  signature <- space$provenance$model_basis_signature
  is.character(signature) && length(signature) == 1L &&
    .strong_sha256(sub("^model-basis-", "", signature))
}

.model_basis_fields <- c(
  "Q", "R", "conditions", "n_conditions", "dimension", "models", "columns",
  "spectra", "distance", "tolerance", "factor_tolerance", "negative_share", "span_fraction",
  "saturated", "span_overlap_rank", "basis_condition", "baseline_invariant",
  "centering_error", "normalize", "centering", "effect_space", "extractor",
  "effect_map", "signature"
)

.model_basis_semantic <- function(x) {
  list(
    schema_version = 2L,
    contract = "predictive-geometry-v1",
    conditions = x$conditions,
    models = x$models,
    distance = x$distance,
    tolerance = x$tolerance,
    factor_tolerance = x$factor_tolerance,
    negative_share = x$negative_share,
    normalize = x$normalize,
    centering = x$centering,
    spectra = x$spectra,
    rank_requested = x$models$rank_requested,
    rank_effective = x$models$rank_effective,
    Q = unname(x$Q),
    R = unname(x$R)
  )
}

# Centering is asserted at a fixed tolerance: `Q = V U~` with `V'1` at
# round-off, so `|Q'1|` is a few `eps` whatever the model's conditioning
# (rule 2 of the header). A failure is therefore a construction bug.
.model_basis_centering_tolerance <- function(n) 1e-10 * sqrt(n)

#' Declare a model family as a basis of the effect axis
#'
#' A model basis turns one or more model geometries over the relation's
#' conditions into a single declared value: an orthonormal, centered basis
#' `Q` of their joint span, and the factor `R` that recovers each model's own
#' coordinates (`T = Q R`, with `T = [T_1, ..., T_k]` the stacked model
#' factors). Lowering a relation through it, `B~_r = Q' B_r`, is an ordinary
#' relation on the model-coordinate effect space, and the ordinary geometry
#' plan then yields the complete geometry in model coordinates,
#' `S_x = Q' G_x Q`. [geometry_spectrum()] and [latent_geometry()] read that
#' form as they read any other; two participants lowered through one basis
#' share one effect-space signature, which is what the population layer
#' pools on. The basis is built before any neural data is read and its
#' signature enters every plan downstream.
#'
#' @section Construction:
#' For each model the centered Gram is formed, `K_i = -H D_i H / 2` for a
#' squared Euclidean RDM, `K_i = H F_i F_i' H` for a feature matrix, or
#' `K_i = H K_input H` for an explicitly supplied kernel, and
#' truncated at a relative eigenvalue tolerance. The requested rank is clamped
#' to the effective rank and the clamp is recorded, never silent: the literal
#' factor with a rank above the model's true rank has a zero column, and a
#' thin QR of it returns a direction that is not centered and leaks a
#' condition-independent baseline into the lowered form. Every
#' decomposition runs in an explicit orthonormal basis of the centered
#' subspace, so a poorly conditioned model (a feature scaled far below the
#' others) cannot leak an uncentered direction either. The stacked factor is
#' then factored by a thin SVD at the same tolerance. Negative Gram
#' eigenvalues are dropped and their mass recorded; when it exceeds
#' `negative_share` of the positive mass the RDM is not a Euclidean geometry
#' and the constructor refuses. Overlapping model spans are admitted for
#' invariant form prediction; `$R` can be rectangular. Readers that require
#' unique redundant coefficients refuse there. `$basis_condition` reports
#' the condition number on the retained union span. Centering, `Q'1 = 0`, is asserted and recorded, because
#' it is exactly the condition under which every reading of the lowered form
#' cancels additive baselines.
#'
#' @section Alignment:
#' Models are aligned to the declared conditions by name, both axes of an
#' kernel or RDM and the rows of a feature matrix; an unlabelled model is admitted only
#' when its dimension equals the condition count. The `$extractor` product
#' lowers a beta block by **position**: [relation()] aligns a labelled beta
#' matrix to the basis's `$conditions` by row name and refuses a mismatch,
#' but an unlabelled block is taken in the basis's condition order, so build
#' the basis with `conditions = relation$effect_space` when in doubt.
#'
#' @section What it is not:
#' The basis is not a geometry, and building it subtracts nothing from the
#' neural side. The rank-limited fit of a geometry to the model span is a
#' latent projection and lives on the latent layer; the signed compressed
#' form is the estimation layer.
#'
#' @param models One squared Euclidean model RDM over the conditions, or a
#'   named list of them. Each is condition-by-condition, symmetric, zero on
#'   the diagonal, and either labelled on both axes with the condition names
#'   (any order) or unlabelled with the condition count as its dimension.
#' @param distance The dissimilarity `models` declare. Required whenever
#'   `models` is supplied, and the only admitted value is
#'   `"squared_euclidean"`: only a squared Euclidean RDM has a centered Gram
#'   whose factor reproduces it. This is a narrower admission than [rsa()]'s,
#'   which reads any symmetric zero-diagonal matrix for a different purpose.
#' @param features One condition-by-feature model matrix, or a named list of
#'   them; rows are conditions, labelled or not under the same rule. A model
#'   given as features enters through `H F F' H`, which equals the Gram of its
#'   squared Euclidean RDM exactly.
#' @param kernels One positive semidefinite condition-by-condition kernel,
#'   or a named list of them. Both axes follow the same alignment rules as
#'   RDMs. The input kernel is checked for positive semidefiniteness before
#'   centering, so centering cannot hide an indefinite constant direction.
#' @param rank `NULL` to keep every model's full effective rank, one positive
#'   whole number shared by all models, or one per model (named by model or in
#'   model order). Ranks are clamped to each model's effective rank and the
#'   clamp is recorded.
#' @param conditions The condition axis the basis is declared over: condition
#'   names, an [effect_space()] (typically the relation's), or a
#'   [condition_space()]. `NULL` reads the order from the first labelled
#'   model. A space lends its shared unit and scale to the model coordinates
#'   and must carry one of each. Supplying a condition space additionally
#'   yields `$effect_map`, the ingestion-path product.
#' @param tolerance Positive relative tolerance under which an eigenvalue or
#'   singular value counts as zero, against the largest absolute one.
#' @param negative_share The largest admitted ratio of dropped negative Gram
#'   mass to positive Gram mass. Above it the model is refused as
#'   non-Euclidean. Roots within `tolerance` of zero count on neither side,
#'   so a Euclidean model records a dropped mass of exactly zero and
#'   `negative_share = 0` (the default) admits it. A larger value explicitly
#'   admits replacement by the positive part; input-kernel repair mass is
#'   recorded separately from centered-Gram mass. Numerical negative roots
#'   are also recorded separately and do not consume this allowance.
#' @param normalize `"none"` preserves each retained model's scale. `"trace"`
#'   divides each retained kernel by its trace, after rank truncation. Use
#'   this to declare comparable scales before pooling model kernels. The
#'   original trace, retained trace and normalization scale are recorded.
#' @return An `effect_model_basis`: a sealed list with the condition-by-basis
#'   matrix `$Q` (orthonormal columns, each summing to zero), the
#'   basis-by-model-coordinate factor `$R` such that `Q %*% R` is the stacked
#'   model factor, the ordered `$conditions`, `$n_conditions`, the basis
#'   `$dimension` `q`, a per-model data frame `$models` (kind, requested and
#'   effective rank, whether the rank was clamped, positive and dropped
#'   negative Gram mass and their ratio), `$columns` naming each model's
#'   columns of `R`, the full Gram `$spectra`, the declared `$distance` and
#'   `$tolerance` and `$negative_share`, the `$span_fraction`
#'   `q / (n - 1)` with `$saturated` true at one, `$span_overlap_rank` (the
#'   model-coordinate count minus union dimension) and `$basis_condition`
#'   (the condition number on the retained union span), `$baseline_invariant` with the measured `$centering_error`, the
#'   model-coordinate `$effect_space`, the `$extractor` that lowers a beta
#'   block (map `t(Q)`), the `$effect_map` over the supplied condition space
#'   or `NULL`, and the `$signature`.
#'   `$normalize` and `$centering` record the scale convention and the
#'   condition-mean centering origin, units and scale. A supplied effect space
#'   retains that space's signature in the centering
#'   declaration, so identical labels cannot silently change their meaning.
#'   A rank cut through a
#'   positive tied eigenspace refuses with reason `ambiguous_rank_cut`;
#'   retaining the whole group is allowed.
#' @section Refusals:
#' Each is an `effect_capability_refusal` (see [catch_refusal()]) in namespace
#' `"model_coordinate"`. A model whose centered Gram carries more negative
#' mass than `negative_share` admits refuses capability
#' `"euclidean_model_geometry"` with reason `model_gram_indefinite`. A
#' positive rank-cut tie refuses `"identified_rank_projection"` with reason
#' `ambiguous_rank_cut`. Overlap itself is admitted: the fitted form does not
#' require a unique coefficient for every redundant model coordinate.
#' @section Identity:
#' The signature covers the ordered conditions, the model names, kinds and
#' declared distance, the tolerance, the requested and effective ranks, and
#' the unrounded `Q` and `R`. Signs of the computed factors are fixed
#' deterministically, but a model with tied Gram eigenvalues has a factor
#' determined only up to a rotation inside the tie, so two platforms whose
#' LAPACK resolves the tie differently can produce equal spans under
#' different signatures. Compare declarations by signature on one platform.
#' Across platforms compare each reconstructed model kernel using its
#' columns of `Q %*% R`; comparing only `Q %*% t(Q)` establishes equal spans
#' but discards the model's directional strengths.
#' @seealso [effect_extractor()] and [effect_map()], the two forms the basis
#'   is lowered through; [relation()] on the betas path and [plan_relation()]
#'   on the ingestion path; [latent_geometry()] for the rank-limited reading
#'   of the lowered form; [rsa()] for the fixed linear reading of model RDMs
#'   that does not lower anything.
#'   Use [fit_geometry()] and [score_geometry()] for regularized prediction
#'   and independent signed gain; `vignette("predictive-geometry")` gives the
#'   complete workflow.
#' @examples
#' conditions <- c("face", "body", "house", "tool")
#'
#' # A category model: animate against inanimate. Its Gram has rank one, so a
#' # requested rank of two is clamped to one and the clamp is on the record.
#' animate <- c(1, 1, 0, 0)
#' category <- outer(animate, animate, function(a, b) (a - b)^2)
#' dimnames(category) <- list(conditions, conditions)
#' basis <- model_basis(list(category = category),
#'   distance = "squared_euclidean", rank = 2)
#' basis
#' basis$models[, c("rank_requested", "rank_effective", "rank_clamped")]
#' max(abs(colSums(basis$Q)))
#'
#' # A second model given as features spans two more centered directions;
#' # with four conditions the family then saturates the effect axis.
#' shape <- cbind(elongation = c(0.2, 0.9, 0.1, 0.8), size = c(1, 2, 3, 1))
#' rownames(shape) <- conditions
#' both <- model_basis(list(category = category), "squared_euclidean",
#'   features = list(shape = shape))
#' both$dimension
#' both$span_fraction
#'
#' # The lowering products: the extractor for betas, and the effect map when
#' # a condition space is supplied.
#' both$extractor$map
#' with_space <- model_basis(list(category = category), "squared_euclidean",
#'   conditions = condition_space(conditions))
#' with_space$effect_map
#'
#' # Duplicate models retain separate factors in one common span.
#' duplicate <- model_basis(
#'   list(a = category, b = category), distance = "squared_euclidean"
#' )
#' duplicate$span_overlap_rank
#' dim(duplicate$R)
#' @export
model_basis <- function(models = NULL, distance = NULL, features = NULL,
                        rank = NULL, conditions = NULL, tolerance = 1e-10,
                        negative_share = 0, kernels = NULL,
                        normalize = c("none", "trace")) {
  normalize <- match.arg(normalize)
  .check_number(tolerance, "tolerance", positive = TRUE)
  if (tolerance >= 1) {
    .input_error("`tolerance` is relative and must be below one.",
      arg = "tolerance", received = .msg_value(tolerance),
      expected = "a positive number below one")
  }
  .check_number(negative_share, "negative_share", nonnegative = TRUE)
  rdms <- .model_basis_entries(models, "models", "model")
  feature_sets <- .model_basis_entries(features, "features", "features")
  kernel_sets <- .model_basis_entries(kernels, "kernels", "kernel")
  if (!length(rdms) && !length(feature_sets) && !length(kernel_sets)) {
    .input_error(
      "Supply at least one model, as `models` (an RDM), `features`, or `kernels`.",
      arg = "models", received = "nothing", expected = "at least one model")
  }
  if (length(rdms)) {
    if (is.null(distance)) {
      .input_error(paste0(
        "`distance` is required with `models`: a model RDM enters the Gram ",
        "construction only under an explicit `distance = ",
        "\"squared_euclidean\"` declaration."
      ), arg = "distance", received = "NULL",
        expected = "\"squared_euclidean\"")
    }
    .check_string(distance, "distance", what = "one dissimilarity name")
    if (!distance %in% .model_basis_distances) {
      .input_error(sprintf(paste0(
        "`distance` must be %s; received %s. Only a squared Euclidean RDM ",
        "has a centered Gram whose factor reproduces it; for another ",
        "dissimilarity supply the features or a Euclidean embedding."
      ), .msg_names(.model_basis_distances), .msg_value(distance)),
        arg = "distance", received = .msg_value(distance),
        expected = .msg_names(.model_basis_distances))
    }
    distance <- unname(distance)
  } else if (!is.null(distance)) {
    .input_error(
      "`distance` declares the dissimilarity of `models`; supply `models` too.",
      arg = "distance", received = .msg_value(distance),
      expected = "NULL without `models`")
  } else {
    distance <- NA_character_
  }
  entries <- c(
    lapply(names(rdms), function(nm) list(name = nm, kind = "rdm",
      value = rdms[[nm]])),
    lapply(names(feature_sets), function(nm) list(name = nm, kind = "features",
      value = feature_sets[[nm]])),
    lapply(names(kernel_sets), function(nm) list(name = nm, kind = "kernel",
      value = kernel_sets[[nm]]))
  )
  model_names <- vapply(entries, `[[`, "", "name")
  names(entries) <- model_names
  if (anyDuplicated(model_names)) {
    .input_error(sprintf(
      "Model names must be unique across `models`, `features` and `kernels`; %s repeats.",
      .msg_names(unique(model_names[duplicated(model_names)]))),
      arg = "features", received = .msg_names(model_names),
      expected = "unique names across the input lists")
  }
  ranks <- .model_basis_ranks(rank, model_names)
  axis <- .model_basis_conditions(conditions, entries)
  n <- length(axis$names)
  V <- .model_basis_centered_frame(n)

  factors <- vector("list", length(entries))
  records <- vector("list", length(entries))
  spectra <- vector("list", length(entries))
  for (i in seq_along(entries)) {
    entry <- entries[[i]]
    where <- sprintf("`%s` entry `%s`",
      switch(entry$kind, rdm = "models", features = "features", kernel = "kernels"),
      entry$name)
    aligned <- .model_basis_align(entry$value, entry$kind, axis$names, where,
      axis$origin)
    admission <- if (entry$kind == "kernel") {
      .model_basis_admit_kernel(aligned, entry$name, tolerance, negative_share)
    } else list(value = aligned, negative_mass = 0, numerical_negative_mass = 0)
    K <- .model_basis_gram(admission$value, entry$kind, V)
    part <- .model_basis_factor(K, V, aligned, entry$name, entry$kind,
      ranks[[i]], tolerance, negative_share, normalize)
    part$record$input_negative_mass <- admission$negative_mass
    part$record$input_numerical_negative_mass <- admission$numerical_negative_mass
    factors[[i]] <- part$factor
    records[[i]] <- part$record
    spectra[[i]] <- part$spectrum
  }
  names(spectra) <- model_names
  records <- do.call(rbind, records)
  rownames(records) <- NULL
  columns <- .model_basis_columns(model_names, records$rank_effective)

  # The stacked factor lives in col(V); factor its centered coordinates and
  # lift the left singular vectors back, so Q'1 = 0 by construction.
  stacked <- do.call(cbind, factors)
  total <- ncol(stacked)
  s <- svd(t(V) %*% stacked)
  keep <- s$d > tolerance * s$d[[1L]]
  q <- sum(keep)
  Q <- V %*% s$u[, keep, drop = FALSE]
  R <- diag(s$d[keep], q) %*% t(s$v[, keep, drop = FALSE])
  flip <- vapply(seq_len(q), function(j) {
    pivot <- which.max(abs(Q[, j]))
    if (Q[pivot, j] < 0) -1 else 1
  }, numeric(1))
  Q <- Q * rep(flip, each = n)
  R <- R * flip
  basis_condition <- s$d[[1L]] / s$d[[q]]
  centering_error <- max(abs(colSums(Q)))
  if (!is.finite(centering_error) ||
      centering_error > .model_basis_centering_tolerance(n)) {
    .invariant_error(sprintf(paste0(
      "The model basis is not centered: max |Q'1| = %.3g exceeds %.3g, so a ",
      "lowered form would not cancel additive baselines."
    ), centering_error, .model_basis_centering_tolerance(n)))
  }
  dimnames(Q) <- list(axis$names, .model_basis_coordinate_names(q))
  dimnames(R) <- list(colnames(Q),
    .model_basis_model_coordinate_names(model_names, records$rank_effective))

  value <- list(
    Q = Q,
    R = R,
    conditions = axis$names,
    n_conditions = as.integer(n),
    dimension = as.integer(q),
    models = records,
    columns = columns,
    spectra = spectra,
    distance = distance,
    tolerance = as.numeric(tolerance),
    factor_tolerance = as.numeric(tolerance),
    negative_share = as.numeric(negative_share),
    span_fraction = q / (n - 1),
    saturated = q == n - 1L,
    span_overlap_rank = as.integer(ncol(R) - q),
    basis_condition = basis_condition,
    baseline_invariant = TRUE,
    centering_error = centering_error,
    normalize = normalize,
    centering = list(method = "condition_mean", conditions = axis$names,
      units = axis$units, scale = axis$scale, declaration = axis$declaration)
  )
  signature <- .sha256_signature(.model_basis_semantic(value),
    "model-basis-sha256:")
  products <- .model_basis_products(Q, axis$names, axis$units, axis$scale,
    axis$space, signature, model_names)
  structure(c(value, products, list(signature = signature)),
    class = "effect_model_basis")
}

# Fail closed on any record whose fields, factors, products or identity
# disagree with one another. The Gram inputs are not retained, so this
# validates the stored value against itself: every field derivable from `Q`,
# `R` and the per-model ranks is re-derived and compared, the products are
# rebuilt from `Q`, and the signature is recomputed over the semantic fields.
.validate_model_basis <- function(x) {
  if (!.sealed_fields(x, "effect_model_basis", .model_basis_fields)) {
    .input_error("Model-basis fields are missing or noncanonical.")
  }
  Q <- x$Q
  R <- x$R
  models <- x$models
  if (!.is_finite_matrix(Q) || !.is_finite_matrix(R) || !is.data.frame(models) ||
      !.is_strings(x$conditions, unique = TRUE) ||
      !identical(nrow(Q), x$n_conditions) || !identical(ncol(Q), x$dimension) ||
      !identical(nrow(R), x$dimension) ||
      !identical(rownames(Q), x$conditions) ||
      !is.integer(models$rank_effective) || any(models$rank_effective < 1L) ||
      !identical(ncol(R), sum(models$rank_effective)) ||
      !identical(names(x$spectra), models$model) ||
      !all(lengths(x$spectra) == x$n_conditions) ||
      !identical(x$factor_tolerance, x$tolerance) ||
      !identical(x$normalize, "none") && !identical(x$normalize, "trace")) {
    .input_error("Model-basis factors or bookkeeping are inconsistent.")
  }
  n <- x$n_conditions
  q <- x$dimension
  d <- svd(R, nu = 0L, nv = 0L)$d
  derived <- list(
    columns = .model_basis_columns(models$model, models$rank_effective),
    Q_names = list(x$conditions, .model_basis_coordinate_names(q)),
    R_names = list(.model_basis_coordinate_names(q),
      .model_basis_model_coordinate_names(models$model, models$rank_effective)),
    span_fraction = q / (n - 1),
    saturated = q == n - 1L,
    span_overlap_rank = as.integer(ncol(R) - q),
    basis_condition = d[[1L]] / d[[length(d)]],
    centering_error = max(abs(colSums(Q)))
  )
  if (!identical(x$columns, derived$columns) ||
      !identical(dimnames(Q), derived$Q_names) ||
      !identical(dimnames(R), derived$R_names) ||
      !identical(x$span_fraction, derived$span_fraction) ||
      !identical(x$saturated, derived$saturated) ||
      !identical(x$span_overlap_rank, derived$span_overlap_rank) ||
      !isTRUE(x$baseline_invariant) ||
      !isTRUE(all.equal(x$basis_condition, derived$basis_condition,
        tolerance = 1e-8)) ||
      !identical(x$centering_error, derived$centering_error) ||
      any(models$rank_clamped != (!is.na(models$rank_requested) &
        models$rank_requested > models$rank_effective))) {
    .contract_error(
      "Model-basis derived fields are not the ones its factors produce."
    )
  }
  if (max(abs(crossprod(Q) - diag(q))) > 1e-8 ||
      derived$centering_error > .model_basis_centering_tolerance(n)) {
    .contract_error("Model-basis columns are not orthonormal and centered.")
  }
  expected <- .sha256_signature(.model_basis_semantic(x), "model-basis-sha256:")
  .check_signature(x$signature, expected,
    "Model-basis identity is inconsistent with its fields.")
  effects <- .validate_effect_space(x$effect_space)
  declaration <- x$centering$declaration
  if (!is.list(declaration) || !identical(names(declaration), c("kind", "signature")) ||
      !.is_string(declaration$kind) ||
      !declaration$kind %in% c("names", "effect_space", "condition_space") ||
      (identical(declaration$kind, "names") && !is.null(declaration$signature)) ||
      (!identical(declaration$kind, "names") &&
       !.strong_sha256(sub("^condition-space-", "", declaration$signature)))) {
    .contract_error("Model-basis source-coordinate declaration is inconsistent.")
  }
  if (!identical(x$centering, list(method = "condition_mean",
      conditions = x$conditions, units = unique(unname(effects$units)),
      scale = unique(unname(effects$scale)), declaration = declaration))) {
    .contract_error("Model-basis centering origin or units are inconsistent.")
  }
  for (i in seq_len(nrow(models))) {
    values <- x$spectra[[i]]
    threshold <- x$tolerance * max(abs(values))
    retained <- sum(values[seq_len(models$rank_effective[[i]])])
    positive <- sum(values[values > threshold])
    expected_record <- c(raw_trace = sum(values), retained_trace = retained,
      normalization_scale = if (x$normalize == "trace") 1 / retained else 1,
      positive_mass = positive, truncated_positive_mass = positive - retained,
      numerical_positive_mass = sum(values[values > 0 & values <= threshold]),
      numerical_negative_mass = -sum(values[values < 0 & values >= -threshold]),
      dropped_negative_mass = -sum(values[values < -threshold]))
    recorded <- unlist(models[i, names(expected_record), drop = FALSE],
      use.names = FALSE)
    block <- R[, x$columns[[i]], drop = FALSE]
    expected_trace <- if (x$normalize == "trace") 1 else retained
    if (!isTRUE(all.equal(recorded, unname(expected_record), tolerance = 1e-10)) ||
        !isTRUE(all.equal(sum(block * block), expected_trace, tolerance = 1e-8))) {
      .contract_error("Model-basis normalization or spectral mass record is inconsistent.")
    }
  }
  if (!identical(effects$provenance$model_basis_signature, x$signature) ||
      !identical(x$extractor$diagnostics$model_basis_signature, x$signature)) {
    .contract_error("Model-basis products do not carry the basis signature.")
  }
  space <- if (is.null(x$effect_map)) NULL else x$effect_map$condition_space
  units <- unique(unname(effects$units))
  scale <- unique(unname(effects$scale))
  rebuilt <- .model_basis_products(Q, x$conditions, units, scale, space,
    x$signature, models$model)
  if (!identical(rebuilt$effect_space, x$effect_space) ||
      !identical(rebuilt$extractor, .validate_effect_extractor(x$extractor)) ||
      !identical(rebuilt$effect_map, x$effect_map)) {
    .contract_error("Model-basis products are inconsistent with its basis.")
  }
  x
}

# The saturation refusal of contract section 6(b). A shared (learned-rotation)
# fit on a basis that spans all of 1-perp reproduces the best rank-s PSD
# approximation of H G H whatever models produced the basis, so it tests
# nothing about them. Isotropic and diagonal structures learn no rotation and
# are never refused on this ground; the layer-5 reader decides which
# structures call this. Saturation is decided from the dimensions, not read
# from the flag, so a forged flag cannot open the gate.
.model_basis_refuse_saturated <- function(basis, structure = "shared") {
  basis <- .validate_model_basis(basis)
  if (basis$dimension < basis$n_conditions - 1L) return(invisible(basis))
  .capability_refusal(
    sprintf(paste0(
      "The model basis spans every centered direction (span_fraction = 1, ",
      "q = %d of n - 1 = %d), so a %s fit reproduces the best rank-limited ",
      "PSD approximation of the centered geometry whatever models were ",
      "supplied, and tests none of them."
    ), basis$dimension, basis$n_conditions - 1L, structure),
    capability = "model_geometry_test",
    namespace = "model_coordinate",
    reasons = c(
      "model_span_saturated",
      sprintf("span_fraction = %g over models %s.", basis$span_fraction,
        .msg_names(basis$models$model))
    ),
    remedies = c(
      "Reduce the model ranks so the family spans a proper subspace.",
      "Use the diagonal or block structure, which learns no rotation.",
      paste0("Evaluate the learned readout under the edge-disjoint cross-fit ",
        "instead of the plug-in fit.")
    )
  )
}

# Pool model kernels in fixed union coordinates. The declared sum and its
# numerically admitted positive support are distinct records. Zero weights
# cannot leave an unpenalized direction available to the estimator.
.model_basis_weights <- function(weights, models) {
  if (is.null(weights) && length(models) == 1L) {
    return(stats::setNames(1, models))
  }
  if (!is.numeric(weights) || !is.null(dim(weights)) ||
      length(weights) != length(models) || any(!is.finite(weights)) ||
      any(weights < 0) || !.is_strings(names(weights), unique = TRUE) ||
      !setequal(names(weights), models)) {
    .input_error("`weights` must supply one finite nonnegative weight per named model.",
      arg = "weights", expected = .msg_names(models))
  }
  if (!is.finite(sum(weights)) || abs(sum(weights) - 1) > 1e-12) {
    .input_error("Model `weights` must sum to one; they are not silently normalized.",
      arg = "weights", expected = "a named simplex vector")
  }
  stats::setNames(as.numeric(weights[models]), models)
}

.model_basis_pool <- function(basis, weights = NULL) {
  basis <- .validate_model_basis(basis)
  weights <- .model_basis_weights(weights, basis$models$model)
  J <- matrix(0, basis$dimension, basis$dimension)
  for (i in seq_along(weights)) {
    # Skip zero-weight factors before multiplying: their support is absent.
    if (weights[[i]] > 0) {
      factor <- basis$R[, basis$columns[[i]], drop = FALSE] * sqrt(weights[[i]])
      J <- J + tcrossprod(factor)
    }
  }
  if (!.is_finite_matrix(J)) {
    .input_error("The pooled model kernel overflows; rescale or trace-normalize its inputs.")
  }
  e <- eigen(J, symmetric = TRUE)
  scale <- max(abs(e$values))
  keep <- e$values > basis$tolerance * scale
  if (!any(keep) || !is.finite(scale)) {
    .capability_refusal("The pooled kernel has no finite positive support.",
      capability = "predictive_model_support", namespace = "predictive_geometry",
      reasons = "empty_model_support", remedies = "Supply a nonzero model with positive weight.")
  }
  U <- .model_basis_fix_signs(e$vectors[, keep, drop = FALSE])
  values <- e$values[keep]
  effective <- tcrossprod(sweep(U, 2L, sqrt(values), "*"))
  result <- list(kernel = J, effective_kernel = effective,
    vectors = U, values = values, spectrum = e$values,
    dimension = as.integer(sum(keep)), union_dimension = basis$dimension,
    condition = max(values) / min(values),
    discarded_positive_mass = sum(e$values[e$values > 0 & !keep]),
    numerical_negative_mass = -sum(e$values[e$values < 0]),
    tolerance = basis$tolerance, weights = weights,
    model_signature = basis$signature)
  result$signature <- .sha256_signature(result, "model-pool-sha256:")
  result
}

#' @export
format.effect_model_basis <- function(x, ...) {
  x <- .validate_model_basis(x)
  sprintf(
    "model_basis<q = %d of %d centered; %s; span %.2f%s>",
    x$dimension, x$n_conditions - 1L,
    .msg_count(nrow(x$models), "model"),
    x$span_fraction, if (isTRUE(x$saturated)) " (saturated)" else ""
  )
}

#' @export
print.effect_model_basis <- function(x, ...) {
  cat(format(x), "\n", sep = "")
  shown <- x$models[, c("model", "kind", "rank_requested", "rank_effective",
    "rank_clamped", "dropped_negative_mass")]
  print(shown, row.names = FALSE)
  cat(sprintf(
    "baseline_invariant: %s (max |Q'1| = %.1e); basis_condition = %.3g\n",
    x$baseline_invariant, x$centering_error, x$basis_condition))
  cat(sprintf("signature: %s\n", .msg_signature(x$signature)))
  invisible(x)
}
