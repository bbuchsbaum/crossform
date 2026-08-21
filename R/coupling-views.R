# Capability-gated views over completed measurement forms ------------------

#' Coupling results and the readings they may carry
#'
#' Every view over a completed `effect_measurement_form` returns one
#' `effect_coupling_result`. There are no subclasses: the seven readings
#' below are distinguished by the `$kind` string, which is a **closed**
#' enumeration. A result of any kind carries the same fourteen fields, is
#' validated the same way, and prints the same way, so code that reads one
#' kind reads them all; only the shape of `$values` differs, and it differs
#' two ways rather than seven.
#'
#' Every kind carries `$values`, `$edge_index` naming the edges those values
#' are keyed to, `$source_plan` and `$source_receipt` identifying the forms it
#' was read from, `$normalization_axis` and `$summary_axis` stating what was
#' divided out and what was summarized over, `$stage_order` listing the
#' pipeline that produced it, `$units` and `$terminology` stating what the
#' numbers may be called, `$edge_completeness`, and a `$signature` over all of
#' it. `$regularization` and `$partition_policy` are present exactly when the
#' kind applied one; they are `NULL` otherwise, never silently defaulted.
#'
#' @section Kinds:
#' `$kind` is one of the following seven values and no others. Constructing a
#' result with an unlisted kind, or with a `$values` shape that disagrees with
#' its kind, is a contract error rather than an accepted record.
#'
#' Two kinds report **blocks**: `$values` is a named list holding one numeric
#' matrix per edge, in `$edge_index$edge_id` order.
#' \itemize{
#'   \item `"effect_coupling"` -- the raw measurement block, with no
#'     covariance claim attached. `$normalization_axis` is `"none"` and
#'     `$units` is `NULL`, because nothing has been divided out and nothing is
#'     claimed about what the numbers mean.
#'   \item `"covariance_coupling"` -- the same blocks, now certified as
#'     repeated-sample covariance: the form established repeated variation of
#'     effective rank above one and symmetric positive self-blocks.
#' }
#'
#' Five kinds report a **table**: `$values` is a data frame keyed by `edge_id`
#' whose remaining columns are fixed by the kind.
#' \itemize{
#'   \item `"pearson_correlation"` -- columns `edge_id`, `correlation`. One
#'     signed scalar per edge, requiring rank-one measurement axes.
#'   \item `"partitioned_pearson_coupling"` -- columns `edge_id`, `value`,
#'     `transform`. A weighted reduction across several source forms; carries
#'     a `$partition_policy` recording the weights, the placement, and the
#'     edge transform, and `$units` is `"correlation"` or `"fisher_z"`
#'     according to that transform.
#'   \item `"canonical_coupling"` -- columns `edge_id`, `mode`,
#'     `canonical_correlation`. Several rows per edge, one per canonical mode
#'     in descending order; carries the `$regularization` whose ridge changed
#'     the values.
#'   \item `"geometry_alignment"` -- columns `edge_id`, `geometry_alignment`.
#'     Static linear CKA/RV-like alignment, normalized over form entries
#'     rather than over experimental samples.
#'   \item `"gaussian_mutual_information"` -- columns `edge_id`,
#'     `information`, `units`. Carries the `$regularization` used for the
#'     underlying canonical spectrum, and a `$terminology` naming the
#'     signature of the [gaussian_covariance_model()] declaration it rests on.
#' }
#'
#' @seealso [effect_coupling()], [covariance_coupling()],
#'   [canonical_coupling()], [geometry_alignment()], and [connectivity()],
#'   which produce these results, and [gaussian_covariance_model()] for the
#'   declaration Gaussian information requires.
#' @family coupling and connectivity views
#' @name effect_coupling_result
NULL

# The closed set of readings an `effect_coupling_result` may carry.
#
# One class, seven kinds -- and this table is what makes that honest rather
# than lazy. Every kind builds the same fourteen sealed fields under the same
# signature scheme, is checked by the same validator, and is read by the same
# two readers (`format()` and `print()`), neither of which branches on `kind`:
# they report it as a recorded string, the way they report `$terminology` and
# `$units`. What genuinely varies is the shape of `$values`, and that varies
# two ways, not seven -- a named list holding one matrix per edge
# ("edge_blocks"), or one data frame keyed by `edge_id` ("edge_table").
# Splitting the class per kind would mint seven types to express a two-way
# distinction.
#
# `columns` is the exact `names(values)` an edge_table kind must produce, so
# the column shape a caller reads is a checked contract rather than whatever
# the view happened to `data.frame()` together. `regularization` and
# `partition_policy` record whether that optional field must be present, so a
# view cannot drop the ridge or the weights that changed its numbers.
#
# See `?effect_coupling_result` for the user-facing statement of this table
# and design/decisions/2026-08-17-coupling-result-kinds.md for why it is one
# class.
.coupling_kinds <- list(
  effect_coupling = list(
    shape = "edge_blocks",
    columns = NULL,
    regularization = FALSE,
    partition_policy = FALSE
  ),
  covariance_coupling = list(
    shape = "edge_blocks",
    columns = NULL,
    regularization = FALSE,
    partition_policy = FALSE
  ),
  pearson_correlation = list(
    shape = "edge_table",
    columns = c("edge_id", "correlation"),
    regularization = FALSE,
    partition_policy = FALSE
  ),
  partitioned_pearson_coupling = list(
    shape = "edge_table",
    columns = c("edge_id", "value", "transform"),
    regularization = FALSE,
    partition_policy = TRUE
  ),
  canonical_coupling = list(
    shape = "edge_table",
    columns = c("edge_id", "mode", "canonical_correlation"),
    regularization = TRUE,
    partition_policy = FALSE
  ),
  geometry_alignment = list(
    shape = "edge_table",
    columns = c("edge_id", "geometry_alignment"),
    regularization = FALSE,
    partition_policy = FALSE
  ),
  gaussian_mutual_information = list(
    shape = "edge_table",
    columns = c("edge_id", "information", "units"),
    regularization = TRUE,
    partition_policy = FALSE
  )
)

.coupling_kind_contract <- function(kind) {
  if (!.is_string(kind) || !kind %in% names(.coupling_kinds)) {
    .input_error(sprintf(
      "Coupling-view kind must be one of %s; received %s.",
      .msg_names(names(.coupling_kinds)),
      if (.is_string(kind)) .msg_names(kind) else
        paste0("a ", class(kind)[[1L]], " value")
    ))
  }
  .coupling_kinds[[kind]]
}

# The one branch any reader of `$values` needs. Two shapes, not seven kinds.
.coupling_value_shape <- function(x) {
  .coupling_kind_contract(x$kind)$shape
}

# Checks the parts of the kind contract that are cheap enough to re-run on
# every validation: the shape of `$values`, its exact column names when it is
# a table, and whether the fields that record what changed the numbers are
# present. The finite-value sweep stays in the constructor; the signature
# already protects the values themselves from drifting after construction.
.check_coupling_kind_contract <- function(kind, values, regularization,
                                          partition_policy) {
  contract <- .coupling_kind_contract(kind)
  if (contract$shape == "edge_table") {
    if (!is.data.frame(values) || !identical(names(values),
        contract$columns)) {
      .contract_error(sprintf(
        "Coupling kind %s reports the columns %s; received %s.",
        .msg_names(kind), .msg_names(contract$columns),
        if (is.data.frame(values)) .msg_names(names(values)) else
          "a non-table value"
      ))
    }
  } else if (is.data.frame(values) || !is.list(values)) {
    .contract_error(sprintf(
      "Coupling kind %s reports one matrix block per edge, not a table.",
      .msg_names(kind)
    ))
  }
  if (contract$regularization != !is.null(regularization)) {
    .contract_error(sprintf(
      "Coupling kind %s %s record the regularization it applied.",
      .msg_names(kind),
      if (contract$regularization) "must" else "must not"
    ))
  }
  if (contract$partition_policy != !is.null(partition_policy)) {
    .contract_error(sprintf(
      "Coupling kind %s %s record a partition policy.",
      .msg_names(kind),
      if (contract$partition_policy) "must" else "must not"
    ))
  }
  invisible(kind)
}

.coupling_result_signature <- function(fields) {
  .sha256_signature(c(
    list(schema_version = 1L), fields
  ))
}

.new_coupling_result <- function(kind, values, x, normalization_axis,
                                 summary_axis, stage_order,
                                 regularization = NULL, units = NULL,
                                 terminology = NULL,
                                 partition_policy = NULL,
                                 source_forms = list(x)) {
  .validate_measurement_form(x, probe = FALSE)
  if (!is.list(source_forms) || length(source_forms) < 1L) {
    .input_error("Coupling results require at least one completed source form.")
  }
  lapply(source_forms, .validate_measurement_form, probe = FALSE)
  if (!is.null(partition_policy)) {
    .validate_coupling_partition_policy(partition_policy)
  }
  .coupling_kind_contract(kind)
  if (!.is_string(normalization_axis, allow_empty = TRUE) ||
      !normalization_axis %in% c("none", "experimental_samples", "neural_features", "partition_pairs", "form_entries") ||
      !.is_string(summary_axis, allow_empty = TRUE) ||
      !summary_axis %in% c("measurement_coordinates", "form_entries", "canonical_modes") ||
      !is.character(stage_order) || length(stage_order) < 1L ||
      anyNA(stage_order)) {
    .input_error("Coupling-view identity or axes are invalid.")
  }
  valid_values <- if (is.data.frame(values)) {
    nrow(values) >= 1L && "edge_id" %in% names(values) &&
      all(values$edge_id %in% x$block_index$edge_id) &&
      all(vapply(values, function(value) {
        if (is.numeric(value)) all(is.finite(value)) else
          all(!is.na(value))
      }, logical(1)))
  } else if (is.list(values)) {
    length(values) == nrow(x$block_index) &&
      identical(names(values), x$block_index$edge_id) &&
      all(vapply(values, function(value) {
        is.matrix(value) && is.numeric(value) && all(is.finite(value))
      }, logical(1)))
  } else {
    FALSE
  }
  if (!valid_values) {
    .contract_error(
      "Coupling-view values do not match their completed edge set."
    )
  }
  .check_coupling_kind_contract(
    kind, values, regularization, partition_policy
  )
  fields <- list(
    kind = kind,
    values = values,
    edge_index = x$block_index,
    source_plan = vapply(source_forms, function(form) {
      form$plan$signature
    }, character(1)),
    source_receipt = vapply(source_forms, function(form) {
      form$receipt$signature
    }, character(1)),
    normalization_axis = normalization_axis,
    summary_axis = summary_axis,
    stage_order = stage_order,
    regularization = regularization,
    units = units,
    terminology = terminology,
    partition_policy = partition_policy,
    edge_completeness = x$edge_completeness
  )
  structure(c(fields, list(
    signature = .coupling_result_signature(fields)
  )), class = "effect_coupling_result")
}

.validate_coupling_result <- function(x) {
  expected <- c("kind", "values", "edge_index", "source_plan",
    "source_receipt", "normalization_axis", "summary_axis", "stage_order",
    "regularization", "units", "terminology", "partition_policy",
    "edge_completeness", "signature")
  if (!.sealed_fields(x, "effect_coupling_result", expected)) {
    .input_error("Coupling result is missing or noncanonical.")
  }
  .check_coupling_kind_contract(
    x$kind, x$values, x$regularization, x$partition_policy
  )
  fields <- x[setdiff(names(x), "signature")]
  .check_signature(
    x$signature, .coupling_result_signature(fields),
    "Coupling-result identity is inconsistent."
  )
  if (!is.character(x$source_plan) || length(x$source_plan) < 1L ||
      anyNA(x$source_plan) || !is.character(x$source_receipt) ||
      length(x$source_receipt) != length(x$source_plan) ||
      anyNA(x$source_receipt)) {
    .input_error("Coupling-result source identities are inconsistent.")
  }
  if (!is.null(x$partition_policy)) {
    .validate_coupling_partition_policy(x$partition_policy)
  }
  invisible(x)
}

.coupling_partition_policy <- function(
    weights, placement = c("within_partition_pair",
                           "after_partition_aggregation"),
    transform = NULL) {
  placement <- match.arg(placement)
  if (!.is_finite_numeric(weights) || length(weights) < 1L || anyNA(weights) ||
      any(weights < 0) || sum(weights) <= 0) {
    .input_error(
      "Partitioned coupling requires finite nonnegative positive-mass weights."
    )
  }
  weights <- as.numeric(weights / sum(weights))
  if (is.null(transform)) transform <- .identity_edge_transform()
  transform <- .validate_edge_transform(transform)
  if (!transform$kind %in% c("identity", "fisher_z")) {
    .input_error(
      "Partitioned Pearson coupling supports identity or Fisher transforms."
    )
  }
  semantic <- list(
    schema_version = 1L,
    placement = placement,
    weights = weights,
    transform = unclass(transform)
  )
  structure(c(semantic[-1L], list(
    signature = .sha256_signature(semantic)
  )), class = "effect_coupling_partition_policy")
}

.validate_coupling_partition_policy <- function(x) {
  expected <- c("placement", "weights", "transform", "signature")
  if (!.sealed_fields(x, "effect_coupling_partition_policy", expected)) {
    .input_error("Coupling partition policy is missing or noncanonical.")
  }
  transform <- structure(x$transform, class = "effect_edge_transform")
  rebuilt <- .coupling_partition_policy(x$weights, x$placement, transform)
  if (!identical(x, rebuilt)) {
    .contract_error("Coupling partition-policy identity is inconsistent.")
  }
  rebuilt
}

.coupling_blocks <- function(x) {
  stats::setNames(lapply(seq_len(nrow(x$block_index)), function(edge) {
    .measurement_block(x, edge)
  }), x$block_index$edge_id)
}

.coupling_self_blocks <- function(x, edge, tolerance = 1e-10,
                                  ridge_left = 0, ridge_right = 0,
                                  require_positive = TRUE) {
  .validate_measurement_form(x, probe = FALSE)
  position <- .measurement_edge_position(x$block_index, edge)
  row <- x$block_index[position, , drop = FALSE]
  if (!identical(x$left_frame$signature, x$right_frame$signature)) {
    .input_error(
      "Normalized coupling requires one compatible self-measurement frame."
    )
  }
  locate <- function(left, right) {
    which(x$block_index$left == left & x$block_index$right == right)
  }
  left_position <- locate(row$left[[1L]], row$left[[1L]])
  right_position <- locate(row$right[[1L]], row$right[[1L]])
  if (length(left_position) != 1L || length(right_position) != 1L) {
    .input_error(
      "Normalized coupling requires explicit left and right self-blocks."
    )
  }
  cross <- .measurement_block(x, position)
  left_self <- .measurement_block(x, left_position)
  right_self <- .measurement_block(x, right_position)
  if (max(abs(left_self - t(left_self))) > tolerance ||
      max(abs(right_self - t(right_self))) > tolerance) {
    .input_error("Coupling self-blocks must be symmetric.")
  }
  left_eigen <- eigen((left_self + t(left_self)) / 2,
    symmetric = TRUE, only.values = TRUE)$values
  right_eigen <- eigen((right_self + t(right_self)) / 2,
    symmetric = TRUE, only.values = TRUE)$values
  if (require_positive &&
      (min(left_eigen) < -tolerance || min(right_eigen) < -tolerance)) {
    .input_error("Indefinite self-blocks cannot normalize coupling.")
  }
  list(
    edge_id = row$edge_id[[1L]],
    cross = cross,
    left_self = left_self,
    right_self = right_self,
    left_regularized = left_self + diag(ridge_left, nrow(left_self)),
    right_regularized = right_self + diag(ridge_right, nrow(right_self)),
    left_eigenvalues = left_eigen,
    right_eigenvalues = right_eigen
  )
}

# The one edge-reduction loop behind the four block-reducing coupling views.
#
# Each view walks the same completed edge set, pulls the same cross- and
# self-blocks for each edge, and rbinds one row per edge. They differ in
# exactly two places, which are the two parameters here:
#
#   * `statistic` -- what one edge reduces to. A view that returns `NULL`
#     reduces to a validation sweep (`.covariance_coupling()` proves every
#     edge has symmetric positive self-blocks and then reports the raw
#     blocks); a view that returns a one-row data frame reports one number per
#     edge; a view that returns several rows reports a spectrum.
#   * the ridge, which is where the normalization is *placed*. Pearson and
#     geometry alignment divide by a scalar built from the untouched
#     self-blocks, so their ridge is zero; canonical coupling whitens by the
#     inverse square root of the self-blocks, so its ridge has to be applied
#     inside `.coupling_self_blocks()`, before this loop sees them.
.coupling_edge_reduce <- function(x, statistic, tolerance = 1e-10,
                                  ridge_left = 0, ridge_right = 0) {
  rows <- lapply(seq_len(nrow(x$block_index)), function(edge) {
    statistic(.coupling_self_blocks(
      x, edge, tolerance, ridge_left, ridge_right
    ))
  })
  if (all(vapply(rows, is.null, logical(1)))) return(invisible(NULL))
  do.call(rbind, rows)
}

.require_nondegenerate_variation <- function(x) {
  if (!isTRUE(x$capabilities$repeated_variation) ||
      is.null(x$capabilities$sampling_axis)) {
    .capability_refusal(paste0(
      "Normalized connectivity requires certified repeated variation along a ",
      "named sampling axis: a correlation divides by a variance, and this ",
      "form has not established that its variation is repeated sampling ",
      "rather than a single fixed pattern."
    ),
      capability = "certified_repeated_variation",
      namespace = "coupling_views",
      reasons = "repeated_variation_not_certified",
      remedies = paste0(
        "Build the form with a `variation_query()` that declares its ",
        "sampling axis, over a pairing whose partitions repeat that axis."
      )
    )
  }
  if (x$diagnostics$experimental_effective_rank <= 1L) {
    .capability_refusal(sprintf(paste0(
      "Normalized connectivity requires an effective sampling rank above ",
      "one; this variation query has effective rank %s. A rank-one variation ",
      "direction carries no independent repeats to normalize by, so every ",
      "edge would report a correlation of plus or minus one by construction."
    ), format(x$diagnostics$experimental_effective_rank)),
      capability = "nondegenerate_variation",
      namespace = "coupling_views",
      reasons = "rank_one_variation_axis",
      remedies = paste0(
        "Supply a variation query of rank two or more, for example a ",
        "centered covariance over several repeated observations, or read the ",
        "unnormalized block with `effect_coupling()`."
      )
    )
  }
  invisible(x)
}

.effect_coupling <- function(x) {
  .require_measurement_view_capabilities(x, "effect_coupling")
  .new_coupling_result(
    "effect_coupling",
    .coupling_blocks(x),
    x,
    normalization_axis = "none",
    summary_axis = "measurement_coordinates",
    stage_order = c(x$plan$stages$order,
      "effect_coupling_materialization"),
    terminology = "algebraic_effect_coupling_no_covariance_claim"
  )
}

.covariance_coupling <- function(x, tolerance = 1e-10) {
  .require_nondegenerate_variation(x)
  .require_measurement_view_capabilities(
    x, "canonical_coupling", positive_self_blocks = TRUE
  )
  .coupling_edge_reduce(x, function(blocks) NULL, tolerance)
  .new_coupling_result(
    "covariance_coupling",
    .coupling_blocks(x),
    x,
    normalization_axis = "none",
    summary_axis = "measurement_coordinates",
    stage_order = c(x$plan$stages$order, "joint_covariance_validation",
      "covariance_materialization"),
    terminology = "repeated_sample_covariance_coupling"
  )
}

.pearson_coupling <- function(x, tolerance = 1e-10) {
  .require_nondegenerate_variation(x)
  .require_measurement_view_capabilities(
    x, "canonical_coupling", positive_self_blocks = TRUE
  )
  values <- .coupling_edge_reduce(x, function(blocks) {
    if (!identical(dim(blocks$cross), c(1L, 1L))) {
      .input_error("Pearson coupling requires rank-one measurement axes.")
    }
    denominator <- sqrt(blocks$left_self[[1L]] * blocks$right_self[[1L]])
    if (!is.finite(denominator) || denominator <= tolerance) {
      # A correlation divides by this; a node that does not vary has no
      # correlation to report, which is a capability refusal and not an
      # arithmetic accident.
      .capability_refusal(sprintf(paste0(
        "Pearson coupling requires strictly positive scalar self-variance; ",
        "edge %s has self-variances %s and %s, whose geometric mean is %s at ",
        "tolerance %s. A node with no measured variation has no correlation."
      ), .msg_names(blocks$edge_id), format(blocks$left_self[[1L]]),
        format(blocks$right_self[[1L]]), format(denominator),
        format(tolerance)),
        capability = "nondegenerate_self_variance",
        namespace = "coupling_views",
        reasons = "self_variance_not_strictly_positive",
        remedies = paste0(
          "Drop the constant node from the edge frame, or read the ",
          "unnormalized block with `effect_coupling()`."
        )
      )
    }
    value <- blocks$cross[[1L]] / denominator
    if (abs(value) > 1 + 10 * tolerance) {
      .input_error(
        "Scalar covariance blocks do not define a valid correlation."
      )
    }
    data.frame(edge_id = blocks$edge_id,
      correlation = max(-1, min(1, value)), stringsAsFactors = FALSE)
  }, tolerance)
  .new_coupling_result(
    "pearson_correlation",
    values,
    x,
    normalization_axis = "experimental_samples",
    summary_axis = "measurement_coordinates",
    stage_order = c(x$plan$stages$order, "joint_covariance_validation",
      "experimental_sample_normalization", "signed_scalar_summary"),
    units = "correlation",
    terminology = "signed_functional_connectivity"
  )
}

.validate_partitioned_coupling_forms <- function(forms) {
  if (!is.list(forms) || length(forms) < 1L) {
    .input_error(
      "Partitioned coupling requires a nonempty list of measurement forms."
    )
  }
  lapply(forms, .validate_measurement_form, probe = FALSE)
  reference <- forms[[1L]]
  compatible <- vapply(forms, function(form) {
    identical(form$block_index, reference$block_index) &&
      identical(form$left_frame$signature, reference$left_frame$signature) &&
      identical(form$right_frame$signature, reference$right_frame$signature) &&
      identical(form$plan$experimental_query,
        reference$plan$experimental_query) &&
      identical(form$capabilities$sampling_axis,
        reference$capabilities$sampling_axis)
  }, logical(1))
  if (!all(compatible)) {
    .contract_error(paste0(
      "Partitioned coupling forms must share exact measurement edges, frames, ",
      "experimental query, and sampling axis."
    ))
  }
  lapply(forms, function(form) {
    .require_nondegenerate_variation(form)
    .require_measurement_view_capabilities(
      form, "canonical_coupling", positive_self_blocks = TRUE
    )
  })
  forms
}

.scalar_correlations_from_blocks <- function(form, blocks,
                                              tolerance = 1e-10) {
  if (!is.list(blocks) || length(blocks) != nrow(form$block_index)) {
    .input_error("Scalar coupling blocks are absent or incomplete.")
  }
  vapply(seq_len(nrow(form$block_index)), function(edge) {
    row <- form$block_index[edge, , drop = FALSE]
    left_self <- which(form$block_index$left == row$left[[1L]] &
      form$block_index$right == row$left[[1L]])
    right_self <- which(form$block_index$left == row$right[[1L]] &
      form$block_index$right == row$right[[1L]])
    if (length(left_self) != 1L || length(right_self) != 1L ||
        !identical(dim(blocks[[edge]]), c(1L, 1L)) ||
        !identical(dim(blocks[[left_self]]), c(1L, 1L)) ||
        !identical(dim(blocks[[right_self]]), c(1L, 1L))) {
      .input_error(paste0(
        "Partitioned Pearson coupling requires rank-one axes and explicit ",
        "self-blocks."
      ))
    }
    denominator <- sqrt(
      blocks[[left_self]][[1L]] * blocks[[right_self]][[1L]]
    )
    if (!is.finite(denominator) || denominator <= tolerance) {
      .input_error(
        "Partitioned Pearson coupling requires positive self-variance."
      )
    }
    value <- blocks[[edge]][[1L]] / denominator
    if (!is.finite(value) || abs(value) > 1 + 10 * tolerance) {
      .input_error(
        "Partitioned covariance blocks do not define valid correlations."
      )
    }
    max(-1, min(1, value))
  }, numeric(1))
}

.partitioned_pearson_coupling <- function(
    forms, weights = rep(1, length(forms)),
    placement = c("within_partition_pair", "after_partition_aggregation"),
    transform = NULL, tolerance = 1e-10) {
  forms <- .validate_partitioned_coupling_forms(forms)
  policy <- .coupling_partition_policy(
    weights, match.arg(placement), transform
  )
  if (length(policy$weights) != length(forms)) {
    .contract_error("Partition coupling weights must match the source forms.")
  }
  transform <- structure(policy$transform, class = "effect_edge_transform")
  reference <- forms[[1L]]
  if (policy$placement == "within_partition_pair") {
    correlations <- vapply(forms, function(form) {
      result <- .pearson_coupling(form, tolerance)
      result$values$correlation[
        match(reference$block_index$edge_id, result$values$edge_id)
      ]
    }, numeric(nrow(reference$block_index)))
    transformed <- apply(correlations, 2L, function(value) {
      drop(.apply_edge_transform(matrix(value, ncol = 1L), transform))
    })
    if (is.null(dim(transformed))) {
      transformed <- matrix(transformed,
        nrow = nrow(reference$block_index), ncol = length(forms))
    }
    values <- drop(transformed %*% policy$weights)
    stage_order <- c(
      "partition_pair_covariance",
      "joint_covariance_validation",
      "experimental_sample_normalization",
      "edge_transform",
      "partition_reduction"
    )
  } else {
    blocks_by_form <- lapply(forms, .coupling_blocks)
    aggregate <- lapply(seq_len(nrow(reference$block_index)), function(edge) {
      Reduce(`+`, Map(function(blocks, weight) {
        blocks[[edge]] * weight
      }, blocks_by_form, policy$weights))
    })
    correlations <- .scalar_correlations_from_blocks(
      reference, aggregate, tolerance
    )
    values <- drop(.apply_edge_transform(
      matrix(correlations, ncol = 1L), transform
    ))
    stage_order <- c(
      "partition_pair_covariance",
      "partition_reduction",
      "joint_covariance_validation",
      "experimental_sample_normalization",
      "edge_transform"
    )
  }
  value_table <- data.frame(
    edge_id = reference$block_index$edge_id,
    value = values,
    transform = transform$kind,
    stringsAsFactors = FALSE
  )
  .new_coupling_result(
    "partitioned_pearson_coupling",
    value_table,
    reference,
    normalization_axis = "experimental_samples",
    summary_axis = "measurement_coordinates",
    stage_order = stage_order,
    units = if (transform$kind == "identity") "correlation" else "fisher_z",
    terminology = "ordered_signed_functional_connectivity",
    partition_policy = policy,
    source_forms = forms
  )
}

.canonical_coupling <- function(x, regularization, tolerance = 1e-10) {
  .require_nondegenerate_variation(x)
  .require_measurement_view_capabilities(
    x, "canonical_coupling", positive_self_blocks = TRUE
  )
  regularization <- .validate_measurement_regularization(regularization)
  if (regularization$kind != "ridge" || regularization$applied) {
    .input_error(paste0(
      "Canonical coupling requires an explicit unapplied ridge policy; the ",
      "view applies it to its self-blocks."
    ))
  }
  values <- .coupling_edge_reduce(x, function(blocks) {
    left_inverse <- .decomposition_inverse_sqrt(
      blocks$left_regularized, tolerance
    )
    right_inverse <- .decomposition_inverse_sqrt(
      blocks$right_regularized, tolerance
    )
    values <- svd(left_inverse %*% blocks$cross %*% right_inverse,
      nu = 0L, nv = 0L)$d
    if (any(values > 1 + 10 * tolerance)) {
      .input_error(
        "Regularized covariance blocks yield canonical values above one."
      )
    }
    data.frame(
      edge_id = blocks$edge_id,
      mode = seq_along(values),
      canonical_correlation = pmin(1, pmax(0, values)),
      stringsAsFactors = FALSE
    )
  }, tolerance, regularization$lambda_left, regularization$lambda_right)
  .new_coupling_result(
    "canonical_coupling",
    values,
    x,
    normalization_axis = "experimental_samples",
    summary_axis = "canonical_modes",
    stage_order = c(x$plan$stages$order, "joint_covariance_validation",
      "ridge_regularization", "experimental_sample_whitening", "svd"),
    regularization = regularization,
    units = "canonical_correlation",
    terminology = "multivariate_covariance_coupling_spectrum"
  )
}

.geometry_alignment <- function(x, tolerance = 1e-10) {
  .require_nondegenerate_variation(x)
  if (!isTRUE(x$capabilities$joint_covariance)) {
    .input_error(
      "Geometry alignment requires a coherent positive joint covariance."
    )
  }
  values <- .coupling_edge_reduce(x, function(blocks) {
    denominator <- sqrt(sum(blocks$left_self^2) * sum(blocks$right_self^2))
    if (!is.finite(denominator) || denominator <= tolerance) {
      .input_error("Geometry alignment requires nonzero self-geometries.")
    }
    value <- sum(blocks$cross^2) / denominator
    if (value > 1 + 10 * tolerance) {
      .input_error(
        "Covariance blocks do not define bounded geometry alignment."
      )
    }
    data.frame(edge_id = blocks$edge_id,
      geometry_alignment = min(1, max(0, value)), stringsAsFactors = FALSE)
  }, tolerance)
  .new_coupling_result(
    "geometry_alignment",
    values,
    x,
    normalization_axis = "form_entries",
    summary_axis = "form_entries",
    stage_order = c(x$plan$stages$order, "joint_covariance_validation",
      "form_entry_frobenius_normalization", "geometry_alignment"),
    units = "linear_cka",
    terminology = paste0(
      "static_geometry_alignment_not_dynamic_informational_connectivity"
    )
  )
}

.gaussian_covariance_model <- function(provenance = list()) {
  .validate_effect_provenance(provenance, "Gaussian-model provenance")
  semantic <- list(
    schema_version = 1L,
    family = "joint_gaussian_covariance",
    fixed = TRUE,
    provenance = provenance
  )
  structure(c(semantic[-1L], list(
    signature = .sha256_signature(semantic)
  )), class = "effect_gaussian_covariance_model")
}

.validate_gaussian_covariance_model <- function(x) {
  expected <- c("family", "fixed", "provenance", "signature")
  if (!.sealed_fields(x, "effect_gaussian_covariance_model", expected) ||
      !identical(x$family, "joint_gaussian_covariance") ||
      !identical(x$fixed, TRUE)) {
    .input_error("Gaussian covariance-model declaration is missing or invalid.")
  }
  rebuilt <- .gaussian_covariance_model(x$provenance)
  if (!identical(x, rebuilt)) {
    .contract_error("Gaussian covariance-model identity is inconsistent.")
  }
  rebuilt
}

.gaussian_information <- function(x, regularization, model,
                                  units = c("nats", "bits"),
                                  tolerance = 1e-10) {
  units <- match.arg(units)
  model <- .validate_gaussian_covariance_model(model)
  canonical <- .canonical_coupling(x, regularization, tolerance)
  split_values <- split(
    canonical$values$canonical_correlation,
    canonical$values$edge_id
  )
  rows <- lapply(names(split_values), function(edge_id) {
    rho <- split_values[[edge_id]]
    if (any(rho >= 1 - tolerance)) {
      .input_error(
        "Gaussian information is infinite at canonical correlation one."
      )
    }
    value <- -0.5 * sum(log1p(-(rho^2)))
    if (units == "bits") value <- value / log(2)
    data.frame(edge_id = edge_id, information = value,
      units = units, stringsAsFactors = FALSE)
  })
  .new_coupling_result(
    "gaussian_mutual_information",
    do.call(rbind, rows),
    x,
    normalization_axis = "experimental_samples",
    summary_axis = "canonical_modes",
    stage_order = c(canonical$stage_order,
      "joint_gaussian_information_transform"),
    regularization = regularization,
    units = units,
    terminology = paste0("Gaussian_mutual_information_model:", model$signature)
  )
}

.connectivity_view <- function(
    x, view = c("correlation", "canonical", "geometry_alignment",
                "gaussian_information"),
    regularization = NULL, model = NULL, units = c("nats", "bits"),
    tolerance = 1e-10) {
  view <- match.arg(view)
  switch(view,
    correlation = .pearson_coupling(x, tolerance),
    canonical = {
      if (is.null(regularization)) {
        .capability_refusal(paste0(
          "Canonical connectivity requires an explicit regularization ",
          "policy: the canonical correlations invert both self-blocks, and ",
          "the ridge that makes that inversion well posed changes the ",
          "numbers, so it is never chosen for you."
        ),
          capability = "declared_regularization",
          namespace = "coupling_views",
          reasons = "regularization_not_declared",
          remedies = "Pass `ridge = ` (and `ridge_right = ` when it differs)."
        )
      }
      .canonical_coupling(x, regularization, tolerance)
    },
    geometry_alignment = .geometry_alignment(x, tolerance),
    gaussian_information = {
      if (is.null(regularization) || is.null(model)) {
        missing_parts <- c(
          if (is.null(regularization)) "`ridge`",
          if (is.null(model)) "`model`"
        )
        .capability_refusal(sprintf(paste0(
          "Gaussian-information connectivity requires explicit ",
          "regularization and an explicit Gaussian model declaration; %s ",
          "%s missing. Mutual information is a modeling claim about the ",
          "joint distribution, not a rescaling of the coupling block, so ",
          "crossform records the declaration instead of assuming it."
        ), paste(missing_parts, collapse = " and "),
          if (length(missing_parts) == 1L) "is" else "are"),
          capability = "declared_gaussian_model",
          namespace = "coupling_views",
          reasons = c(
            if (is.null(model)) "gaussian_model_not_declared",
            if (is.null(regularization)) "regularization_not_declared"
          ),
          remedies = paste0(
            "Pass `model = gaussian_covariance_model(...)` recording the ",
            "assumption, together with an explicit `ridge`."
          )
        )
      }
      .gaussian_information(
        x, regularization, model, match.arg(units), tolerance
      )
    }
  )
}
