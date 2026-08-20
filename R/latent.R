# The latent PSD descriptive layer -------------------------------------------
#
# `design/conservative-geometry-contract.md` section 6 names two layers and
# says they must not be confused. The signed estimation layer holds the
# crossvalidated estimates: unbiased precisely because no within-partition
# term is ever formed (`design/effect-form-contract.md` section 8), signed as
# the visible cost of that, and conserving (section 2) as a statement about a
# signed sum. Fractions, cumulative-contribution curves, effective counts and
# every other functional that needs a nonnegative partition are invalid on it
# -- on the contract's own pure-noise fixture the implied shares span
# `[-0.304, 0.599]`, and clipping to nonnegativity inflates the conserved
# budget by `+53.45%`.
#
# The latent PSD layer is where those functionals are legal. It is not the
# estimation layer with a cosmetic repair: it is a declared, biased,
# non-conserving projection of it, and this file is the only place in the
# package that produces one.
#
# Three things follow from section 6's normative clauses, and each is a design
# decision here rather than a caveat in the documentation:
#
#  1. THE PROJECTION IS NAMED, FROM A CLOSED SET. Section 11.4 gap G7 says
#     "nonnegativity projection" is underdetermined -- a per-node total clamp
#     and an eigenvalue truncation of the form are different operators moving
#     different mass. `.latent_projection_methods` is that closed set.
#     `"psd_projection"` -- eigenvalue truncation at zero of the symmetric
#     form -- is the one implemented; `"nearest_psd"` is reserved and refuses
#     by capability rather than by a base-R `match.arg()` error, so a caller
#     asking for it learns that the name is real and unbuilt rather than that
#     it is a typo.
#
#  2. CLIPPING IS NEVER SILENT. Every measurement carries the absolute mass
#     the projection moved and that mass as a share of the source spectrum's
#     total absolute mass, the object carries them, the projection receipt
#     carries them beside the operator that moved them, and the print shows
#     them. A measurement whose source spectrum was already nonnegative says
#     so with a zero rather than by omission.
#
#  3. THE SIGNED SOURCE IS NOT TOUCHED. Constructing a latent layer reads the
#     source and returns a new object; it never rewrites a stored geometry,
#     and no fraction, cumulative curve or effective count is added to any
#     signed view. `geometry_spectrum()` keeps its
#     `$indefinite_estimates_preserved = TRUE` and keeps returning negative
#     roots.
#
# WHY THE EXECUTION RECEIPT IS NOT WIDENED. `effect_execution_receipt` is a
# sealed layer-2 record whose canonical field list every receipt in the
# package shares; adding a latent-only field to it would put an empty slot on
# every receipt the package has ever written. The projection instead enters
# the receipt two ways it already has room for -- the derived
# `$scientific_plan_id`, so a latent layer never shares an identity with the
# signed source it came from, and `$task_partition_id`, which gains
# `+psd_projection` so the operator is readable off the record without
# decoding a digest -- and the moved mass lives in the projection receipt the
# object carries at `$projection`, which is a record of exactly this
# operation and of nothing else.

# The sentence every latent-layer record in the package prints, written once.
# Two layers carry it -- this file's projection of a signed spectrum and
# `population-prevalence.R`'s fractions over participants -- and they are
# different operations that fail the same way, so a reader has to meet the
# same words on both. A second copy would be a second chance to soften one.
.latent_reading_line <- "latent descriptive layer; not for inference"

# The closed set of section 11.4 G7, and the subset built today. A method is
# in `.latent_projection_methods` because it is a projection the latent layer
# admits; it is in `.latent_projection_implemented` because it exists.
.latent_projection_methods <- c("psd_projection", "nearest_psd")
.latent_projection_implemented <- "psd_projection"

.latent_projection_operators <- c(
  psd_projection = "eigenvalue_truncation_at_zero",
  nearest_psd = "nearest_psd_in_frobenius_norm"
)

.latent_operator_phrase <- function(method) {
  switch(method,
    psd_projection = "eigenvalue truncation at zero",
    nearest_psd = "nearest PSD form in Frobenius norm",
    method
  )
}

.latent_method_list <- function() {
  paste(sprintf("\"%s\"", .latent_projection_methods), collapse = ", ")
}

.latent_projection_method <- function(method) {
  # The formal default is the vector of implemented names, so the first entry
  # is the default the same way `match.arg()` would take it. The membership
  # test is written out rather than delegated to `match.arg()` because the
  # closed set is wider than the implemented set, and the two failures are
  # different: an unknown name is a typo, a reserved name is a refusal.
  if (is.character(method) && length(method) > 1L) method <- method[[1L]]
  .check_string(method, "method")
  if (!method %in% .latent_projection_methods) {
    .input_error(sprintf(paste0(
      "`method` must name one of the latent layer's projections (%s); ",
      "received %s."
    ), .latent_method_list(), .msg_value(method)),
      arg = "method", received = .msg_value(method),
      expected = sprintf("one of %s", .latent_method_list()))
  }
  if (!method %in% .latent_projection_implemented) {
    .capability_refusal(sprintf(paste0(
      "`method = \"%s\"` is a declared member of the latent layer's ",
      "projection set and is not implemented. It moves different mass than ",
      "`\"psd_projection\"` does, so it cannot be silently substituted for ",
      "it."
    ), method),
      capability = paste0(method, "_projection"),
      namespace = "latent_geometry",
      reasons = "projection_kind_declared_but_not_implemented",
      remedies = paste0(
        "Use `method = \"psd_projection\"`, which truncates the symmetric ",
        "form's eigenvalues at zero and records the mass it moved."
      )
    )
  }
  method
}

# What the latent layer accepts, and the reasoned refusal for everything else.
#
# Two routes, both reading the same object: the per-measurement symmetric form
# `geometry_spectrum()` decomposes, and the signed spectrum it returns. They
# agree by construction, because the form route is the spectrum route with
# `geometry_spectrum()` called first -- a PSD projection of a symmetric form
# *is* the truncation of its eigenvalues, so there is no second numerical path
# to keep in agreement.
.latent_source <- function(x, component, component_given, row_block) {
  if (inherits(x, "effect_spectrum_view")) {
    if (component_given && !identical(component, x$component)) {
      .input_error(sprintf(paste0(
        "`component` is fixed by the spectrum that was passed: this view ",
        "decomposed `%s`, and `component = \"%s\"` asks for a different one. ",
        "A spectrum cannot be re-decomposed."
      ), x$component, component),
        arg = "component", received = component,
        expected = sprintf("\"%s\", or no `component` argument", x$component))
    }
    if (!.is_finite_matrix(x$values) || NROW(x$index) != nrow(x$values)) {
      .input_error("Spectrum values or measurement index are invalid.")
    }
    # `effect_spectrum_view` is an unsealed record, and `C(k)` is a statement
    # about roots taken largest first: read out of order it is still a
    # monotone curve to 1 and is no longer a contribution curve. Checked
    # rather than sorted, because silently reordering a view that says it is
    # ordered would repair the label and not the object.
    if (ncol(x$values) > 1L &&
        any(x$values[, -1L, drop = FALSE] >
            x$values[, -ncol(x$values), drop = FALSE] + 1e-12)) {
      .input_error(paste0(
        "A spectrum passed to `latent_geometry()` must be ordered largest ",
        "root first, as `geometry_spectrum()` returns it: the cumulative ",
        "contribution curve is defined on that order."
      ),
        arg = "x", received = "a spectrum whose roots are not descending",
        expected = "an `effect_spectrum_view` from `geometry_spectrum()`")
    }
    return(list(kind = "effect_spectrum_view", spectrum = x$values,
      component = x$component, index = x$index, receipt = x$receipt))
  }
  if (inherits(x, "effect_geometry_plan") || inherits(x, "effect_form")) {
    # Raised here rather than left to `geometry_spectrum()` so that the
    # refusal names the operation the caller asked for.
    .require_self_form_view(x, "A latent geometry")
    spectrum <- geometry_spectrum(x, component = component,
      row_block = row_block)
    return(list(kind = "effect_form", spectrum = spectrum$values,
      component = spectrum$component, index = spectrum$index,
      receipt = spectrum$receipt))
  }
  .latent_unsupported(x)
}

.latent_readouts <- c("effect_contrast_view", "effect_view", "effect_rdm_view",
  "effect_rsa_view", "effect_crossnobis_view")

.latent_unsupported <- function(x) {
  if (inherits(x, .latent_readouts)) {
    .capability_refusal(sprintf(paste0(
      "A latent geometry projects a symmetric form, and %s is a query ",
      "readout: its values are already contracted against fixed weights, so ",
      "there is no spectrum left to truncate. Clamping such a value at zero ",
      "is a per-node total clamp, which is a different projection moving ",
      "different mass, and it is not in the latent layer's set."
    ), .msg_value(x)),
      capability = "latent_projection_source",
      namespace = "latent_geometry",
      reasons = c("readout_is_not_a_symmetric_form",
        "per_node_clamp_is_a_different_projection"),
      remedies = paste0(
        "Pass the geometry the readout came from -- ",
        "`latent_geometry(materialize_geometry(plan))` -- or the signed ",
        "spectrum `geometry_spectrum()` reads from it."
      )
    )
  }
  .input_error(sprintf(paste0(
    "`x` must be a complete effect form from `materialize_geometry()` or an ",
    "`effect_spectrum_view` from `geometry_spectrum()`; received %s."
  ), .msg_value(x)),
    arg = "x", received = .msg_value(x),
    expected = "a complete `effect_geometry` or an `effect_spectrum_view`")
}

# The projection itself, and its accounting.
#
# `pmax(lambda, 0)` preserves the descending order `geometry_spectrum()`
# returns, so the projected spectrum needs no re-sort. Three quantities come
# out of it, and the denominators are the whole argument:
#
#   * MOVED MASS is `sum(|lambda|)` over the negative roots -- the absolute
#     mass the truncation removed. It is not the change in the trace with a
#     sign attached, because the reader has to be able to compare it against
#     something, which is:
#   * TOTAL ABSOLUTE MASS, `sum(|lambda|)` over the *signed* source spectrum.
#     The signed trace is the wrong denominator: it can be small, zero or
#     negative while the form is far from PSD, and a share against it would be
#     unbounded or undefined exactly where the warning matters most.
#   * n_eff is the participation ratio `(sum lambda)^2 / sum lambda^2` on the
#     PROJECTED spectrum, per section 6. On the signed spectrum it is not a
#     count of anything: the numerator is the square of a signed sum.
#
# A measurement whose projected spectrum sums to zero -- every root at or
# below zero, or an exactly zero form -- has no nonnegative partition at all,
# so `n_eff` and `C(k)` are masked rather than reported as `0/0`. That is the
# same discipline `.coherence_fraction()` applies to a fraction it cannot
# earn, and the count of masked measurements is on the projection receipt.
# The guard is `> 0` exactly, not `> tol`, for the same reason
# `.coherence_fraction()`'s is: section 11.4 gap G3 records the question of
# whether the package's nonnegativity guards take a relative tolerance as an
# open contract decision, shared with every per-node fraction here, and D7
# inherits that decision rather than making a private one. The consequence is
# stated in the `Masking` section of `latent_geometry()` so a reader of a form
# that is zero only to within rounding is not surprised by it.

# `$n_eff` and `$cumulative` are functions of the projected spectrum alone, so
# they are computed in one place and RE-DERIVED by the validator rather than
# trusted. Two of the record's twelve fields would otherwise be forgeable
# without contradiction -- a plausible `n_eff` beside a spectrum that does not
# produce it is exactly the kind of number this layer exists to make
# impossible.
.latent_functionals <- function(projected) {
  projected <- unname(projected)
  kept <- rowSums(projected)
  square <- rowSums(projected^2)
  supported <- kept > 0 & square > 0

  n_eff <- rep(NA_real_, length(kept))
  n_eff[supported] <- kept[supported]^2 / square[supported]

  cumulative <- projected
  if (ncol(projected) > 1L) {
    cumulative <- matrix(
      apply(projected, 1L, cumsum), nrow(projected), ncol(projected),
      byrow = TRUE
    )
  }
  cumulative[supported, ] <- cumulative[supported, , drop = FALSE] /
    kept[supported]
  cumulative[!supported, ] <- NA_real_
  dimnames(cumulative) <- list(NULL, paste0("C", seq_len(ncol(projected))))

  list(n_eff = n_eff, cumulative = cumulative, kept = kept,
    supported = supported)
}

.latent_psd_projection <- function(signed) {
  projected <- pmax(signed, 0)
  moved_mass <- rowSums(pmax(-signed, 0))
  absolute_mass <- rowSums(abs(signed))

  # A form with no absolute mass at all moved none of it. The share is `0`
  # rather than masked, which is the one place this file departs from
  # `.coherence_fraction()`'s discipline and does so deliberately: `NA` here
  # would read as "we could not tell you how much was clipped" about a
  # measurement where the answer is known and is none.
  moved_share <- rep(0, length(moved_mass))
  positive <- absolute_mass > 0
  moved_share[positive] <- moved_mass[positive] / absolute_mass[positive]

  functionals <- .latent_functionals(projected)
  dimnames(projected) <- list(NULL, paste0("root", seq_len(ncol(projected))))

  list(spectrum = projected, cumulative = functionals$cumulative,
    n_eff = functionals$n_eff, moved_mass = moved_mass,
    moved_share = moved_share, absolute_mass = absolute_mass,
    supported = functionals$supported)
}

.latent_scientific_id <- function(parent, method, component) {
  .sha256_signature(list(
    schema_version = 1L,
    role = "latent_geometry",
    parent = parent,
    method = method,
    operator = .latent_projection_operators[[method]],
    component = component
  ), "latent-sha256:")
}

# `.projection_receipt()` marks a view as projected from a parent execution.
# The latent layer names *which* projection on top of that, because section 6
# requires the operation to be visible on the record and not only inside a
# digest a reader cannot decode.
.latent_execution_receipt <- function(parent, scientific_plan_id, method) {
  receipt <- .projection_receipt(parent, scientific_plan_id)
  receipt$task_partition_id <- paste0(receipt$task_partition_id, "+", method)
  .validate_execution_receipt(receipt)
  receipt
}

.new_latent_projection_receipt <- function(method, component, projection,
                                           source_kind, source_plan_id,
                                           scientific_plan_id) {
  structure(
    list(
      method = method,
      operator = .latent_projection_operators[[method]],
      component = component,
      measurements = length(projection$moved_mass),
      clipped = sum(projection$moved_mass > 0),
      masked = sum(!projection$supported),
      moved_mass = projection$moved_mass,
      moved_share = projection$moved_share,
      absolute_mass = projection$absolute_mass,
      total_moved_mass = sum(projection$moved_mass),
      max_moved_share = if (length(projection$moved_share)) {
        max(projection$moved_share)
      } else {
        0
      },
      source = source_kind,
      source_scientific_plan_id = source_plan_id,
      scientific_plan_id = scientific_plan_id
    ),
    class = "effect_latent_projection_receipt"
  )
}

.validate_latent_projection_receipt <- function(x) {
  expected <- c("method", "operator", "component", "measurements", "clipped",
    "masked", "moved_mass", "moved_share", "absolute_mass",
    "total_moved_mass", "max_moved_share", "source",
    "source_scientific_plan_id", "scientific_plan_id")
  if (!.sealed_fields(x, "effect_latent_projection_receipt", expected)) {
    .input_error("`x` must be a canonical latent projection receipt.")
  }
  if (!x$method %in% .latent_projection_implemented ||
      !identical(x$operator, .latent_projection_operators[[x$method]])) {
    .contract_error("A latent projection receipt names an unbuilt operator.")
  }
  if (!.is_finite_numeric(x$moved_mass) || any(x$moved_mass < 0) ||
      !.is_finite_numeric(x$moved_share) || any(x$moved_share < 0) ||
      any(x$moved_share > 1) ||
      !.is_finite_numeric(x$absolute_mass) || any(x$absolute_mass < 0) ||
      length(x$moved_share) != length(x$moved_mass) ||
      length(x$absolute_mass) != length(x$moved_mass) ||
      !identical(length(x$moved_mass), as.integer(x$measurements))) {
    .input_error("Latent projection moved mass is missing or inconsistent.")
  }
  # The three summaries are functions of the per-measurement series, and they
  # are the fields the print reads. Left unchecked, a receipt could carry the
  # honest per-measurement mass and a `clipped: 0 of 3` headline over it,
  # which is precisely the silent clipping section 6 forbids.
  if (!identical(x$clipped, sum(x$moved_mass > 0)) ||
      !isTRUE(all.equal(x$total_moved_mass, sum(x$moved_mass),
        tolerance = 1e-12)) ||
      !isTRUE(all.equal(x$max_moved_share,
        if (length(x$moved_share)) max(x$moved_share) else 0,
        tolerance = 1e-12))) {
    .contract_error(paste0(
      "A latent projection receipt summarizes moved mass it does not carry."
    ))
  }
  invisible(x)
}

.new_effect_latent_geometry <- function(projection, component, method, index,
                                        projection_receipt, receipt,
                                        metadata) {
  if (!is.list(metadata)) .input_error("`metadata` must be a list.")
  metadata$scientific_plan_id <- receipt$scientific_plan_id
  value <- structure(
    list(
      spectrum = projection$spectrum,
      cumulative = projection$cumulative,
      n_eff = projection$n_eff,
      moved_mass = projection$moved_mass,
      moved_share = projection$moved_share,
      component = component,
      method = method,
      index = index,
      projection = projection_receipt,
      receipt = receipt,
      metadata = metadata,
      contract_signature = .latent_contract_signature(
        method, component, dim(projection$spectrum),
        receipt$scientific_plan_id
      )
    ),
    class = "effect_latent_geometry"
  )
  .validate_effect_latent_geometry(value)
  value
}

.latent_contract_signature <- function(method, component, shape,
                                       scientific_plan_id) {
  .sha256_signature(list(
    schema_version = 1L,
    result_capability = "latent_descriptive_layer",
    method = method,
    operator = .latent_projection_operators[[method]],
    component = component,
    shape = as.integer(shape),
    scientific_plan_id = scientific_plan_id
  ))
}

.validate_effect_latent_geometry <- function(x) {
  expected <- c("spectrum", "cumulative", "n_eff", "moved_mass", "moved_share",
    "component", "method", "index", "projection", "receipt", "metadata",
    "contract_signature")
  if (!.sealed_fields(x, "effect_latent_geometry", expected)) {
    .input_error("`x` must be a canonical `effect_latent_geometry`.")
  }
  if (!.is_finite_matrix(x$spectrum) || any(x$spectrum < 0)) {
    .contract_error(paste0(
      "A latent spectrum must be nonnegative; this one is not, so the ",
      "projection that produced it did not run."
    ))
  }
  if (!is.matrix(x$cumulative) || !identical(dim(x$cumulative),
      dim(x$spectrum))) {
    .input_error("A latent cumulative curve must match its spectrum's shape.")
  }
  finite_curve <- x$cumulative[is.finite(x$cumulative)]
  if (any(finite_curve < -1e-12) || any(finite_curve > 1 + 1e-12)) {
    .contract_error("A latent cumulative curve leaves [0, 1].")
  }
  if (!is.numeric(x$n_eff) || length(x$n_eff) != nrow(x$spectrum) ||
      any(x$n_eff[is.finite(x$n_eff)] < 1 - 1e-12) ||
      any(is.nan(x$n_eff))) {
    .contract_error(paste0(
      "A latent effective count must be at least one wherever it is defined, ",
      "and masked where it is not."
    ))
  }
  # Re-derived, not range-checked. `n_eff` and `C(k)` are functions of the
  # projected spectrum, so a plausible value beside a spectrum that does not
  # produce it is detectable and must be detected: these are the two numbers
  # the layer exists to supply, and a reader has no other way to check them.
  # Both sides run the same computation on the same input, so an honest
  # object matches bit for bit.
  derived <- .latent_functionals(x$spectrum)
  if (!identical(derived$n_eff, unname(x$n_eff)) ||
      !identical(derived$cumulative, x$cumulative)) {
    .contract_error(paste0(
      "A latent effective count or cumulative curve is not the one its own ",
      "spectrum produces."
    ))
  }
  if (!x$method %in% .latent_projection_implemented ||
      !x$component %in% c("total", "coherent", "configuration")) {
    .input_error("A latent layer names an unknown projection or component.")
  }
  if (NROW(x$index) != nrow(x$spectrum) || anyNA(x$index)) {
    .input_error("`index` must have one entry per measurement.")
  }
  .validate_latent_projection_receipt(x$projection)
  # The one identity that ties the moved mass to the spectrum that survived
  # it: what was kept plus what was moved is the source's total absolute
  # mass. Without it `moved_mass` is unanchored, and a record claiming to
  # have clipped nothing sits consistently beside a spectrum that lost four
  # units. It is a floating-point identity between a one-pass sum and a
  # two-pass one, so it is asserted at the numerical contract's tolerance
  # rather than exactly.
  scale <- max(1, max(x$projection$absolute_mass))
  if (any(abs(rowSums(x$spectrum) + x$moved_mass -
      x$projection$absolute_mass) > 1e-10 * scale)) {
    .contract_error(paste0(
      "A latent layer's kept and moved mass do not add back to the absolute ",
      "mass of the spectrum it was projected from."
    ))
  }
  if (!identical(x$projection$masked, sum(!derived$supported)) ||
      !identical(x$projection$measurements, nrow(x$spectrum))) {
    .contract_error(
      "A latent projection receipt miscounts the measurements it describes."
    )
  }
  if (!identical(x$projection$moved_mass, x$moved_mass) ||
      !identical(x$projection$moved_share, x$moved_share) ||
      !identical(x$projection$method, x$method) ||
      !identical(x$projection$component, x$component)) {
    .contract_error(paste0(
      "A latent layer and its projection receipt disagree about what was ",
      "projected or how much mass it moved."
    ))
  }
  .validate_execution_receipt(x$receipt)
  if (!identical(x$projection$scientific_plan_id,
      x$receipt$scientific_plan_id)) {
    .contract_error(
      "A latent projection receipt and its execution receipt disagree."
    )
  }
  if (!is.list(x$metadata) || !identical(x$metadata$scientific_plan_id,
      x$receipt$scientific_plan_id)) {
    .contract_error(
      "Latent metadata and receipt identify different scientific plans."
    )
  }
  .check_signature(
    x$contract_signature,
    .latent_contract_signature(x$method, x$component, dim(x$spectrum),
      x$receipt$scientific_plan_id),
    "Latent contract signature is inconsistent with its claims."
  )
  invisible(x)
}

#' Construct the latent PSD descriptive layer of a signed geometry
#'
#' Crossvalidated contributions are **signed**, and the arithmetic that treats
#' a node's value as part of a nonnegative whole is invalid on them
#' (`design/conservative-geometry-contract.md` section 6): the implied shares
#' can leave `[0, 1]`, and clipping to nonnegativity both reintroduces the
#' bias the cross-partition estimator removes and destroys conservation.
#' Effective counts, cumulative-contribution curves and fractions are legal
#' only on a declared nonnegative projection of the estimates, and
#' `latent_geometry()` is the named operation that builds one.
#'
#' The projection is biased and does not conserve. That is not a defect to be
#' worked around; it is what makes the layer's arithmetic well defined, and it
#' is why the moved mass is reported per measurement rather than absorbed.
#' Read the layer descriptively -- which modes carry a form's nonnegative mass
#' and how concentrated they are -- and read the signed geometry for anything
#' that has to be unbiased or conserved.
#'
#' @param x A complete symmetric effect form from [materialize_geometry()], or
#'   the signed `effect_spectrum_view` [geometry_spectrum()] reads from one.
#'   A query readout has already been contracted against fixed weights and
#'   carries no spectrum, so it is refused.
#' @param method The projection to apply, from the latent layer's closed set.
#'   `"psd_projection"` truncates the symmetric form's eigenvalues at zero.
#'   `"nearest_psd"` is a declared member of the set and is not implemented;
#'   it refuses rather than falling back, because it moves different mass.
#' @param component Geometry component to project. Fixed by the view when `x`
#'   is an `effect_spectrum_view`.
#' @param row_block Positive number of measurement rows read per block when
#'   `x` is a form.
#' @return An `effect_latent_geometry`.
#' @section Structure:
#' One nonnegative spectrum per spatial measurement, with the functionals that
#' only a nonnegative partition admits and the record of what the projection
#' cost.
#'
#' - `$spectrum`: the projected spectrum, one row per measurement, one column
#'   per root, ordered `root1` (largest) through `rootq`. Every entry is
#'   nonnegative.
#' - `$cumulative`: the cumulative contribution curve
#'   \eqn{C(k) = \sum_{i \le k} \lambda_i / \sum_i \lambda_i}{C(k) = sum_(i <= k) lambda_i / sum_i lambda_i}, columns `C1`
#'   through `Cq`, so `Cq` is `1` wherever the curve is defined.
#' - `$n_eff`: the participation ratio
#'   \eqn{(\sum_i \lambda_i)^2 / \sum_i \lambda_i^2}{(sum_i lambda_i)^2 / sum_i lambda_i^2} of the projected
#'   spectrum, one per measurement.
#' - `$moved_mass`: the absolute mass the projection removed from each
#'   measurement -- the sum of the magnitudes of its negative roots. Zero for
#'   a measurement whose source spectrum was already nonnegative.
#' - `$moved_share`: that mass as a fraction of the source spectrum's total
#'   absolute mass, so a measurement that moved a lot in a large form is not
#'   confused with one that moved a little in a small one.
#' - `$component`, `$method`: what was projected, and by which named operator.
#' - `$index`: the measurement identifiers, carried from the source.
#' - `$projection`: the projection receipt -- the operator, the per-measurement
#'   moved mass and share, the counts of clipped and masked measurements, and
#'   the source identity it was derived from.
#' - `$receipt`: the execution receipt. Its `$scientific_plan_id` is derived
#'   from the source's and the projection's name, so a latent layer never
#'   shares an identity with the signed estimates behind it, and its
#'   `$task_partition_id` ends in `+psd_projection`.
#'
#' Any element not listed here is internal and may change.
#' @section Masking:
#' A measurement whose projected spectrum sums to zero -- every root at or
#' below zero -- has no nonnegative partition, so its `$n_eff` and its whole
#' `$cumulative` row are `NA` rather than `0/0`. `$projection$masked` counts
#' them. This is the discipline the coherence fraction already uses: a number
#' that was not earned is withheld, not clamped.
#'
#' The guard is exact (`> 0`), not a tolerance, and it is the same guard every
#' per-node fraction in the package uses. A form that is zero only to within
#' rounding therefore returns `$n_eff` and `$cumulative` computed from that
#' rounding: a row of roots at `1e-17` gives the same `$n_eff` as a row a
#' hundred million times larger, because the ratio is scale free.
#' `$moved_share` will read near `0.5` for such a row and `$spectrum` will
#' show the scale, so the object says what happened -- but the effective count
#' of a form that is numerically zero is noise, and no guard here will tell
#' you so. Whether these guards should take a relative tolerance is an open
#' contract decision (`design/conservative-geometry-contract.md` §11.4, gap
#' G3, where a singleton scale's aggregated configuration was measured at
#' `-2.8e-17`); this layer inherits it rather than making a private one.
#'
#' `$moved_share` is the one quantity that is *not* masked when its
#' denominator vanishes: a form with no absolute mass at all moved none of it,
#' so the share is `0`. `NA` there would read as "we could not tell you how
#' much was clipped" about a measurement where the answer is known.
#' @section Refusal:
#' A query-first `effect_geometry_plan` signals an
#' `effect_capability_refusal` with capability `"complete_geometry"`; a
#' rectangular or non-symmetric form signals `"symmetric_self_form"`; a query
#' readout (`contrast_energy()`, `rdm()`, `rsa()`, `crossnobis()`, or a
#' `query_geometry()` view) signals `"latent_projection_source"`, because
#' clamping a contracted value is a per-node total clamp and a different
#' projection. `method = "nearest_psd"` signals
#' `"nearest_psd_projection"`. See [catch_refusal()].
#' @seealso [geometry_spectrum()] for the signed spectrum this layer projects,
#'   which never truncates a negative root, and [contribution()] for the
#'   additive ledger reading that stays on the signed layer.
#' @family geometry plans and views
#' @examples
#' domain <- abstract_domain(4, id = "latent-example")
#' relation <- relation(
#'   list(
#'     run1 = rbind(face = c(1, 0, 2, 1), house = c(0, 1, 1, 0),
#'       tool = c(0.5, -1, 0, 1)),
#'     run2 = rbind(face = c(0.9, 0.2, 1.7, 0.6), house = c(0.2, 0.8, 1.3, 0.1),
#'       tool = c(-0.4, 1, 0.2, -0.8))
#'   ),
#'   domain = domain
#' )
#' geometry <- materialize_geometry(plan_geometry(
#'   relation, compile_frame(regions(c("v1", "v1", "it", "it")), domain),
#'   cross_partitions(relation, independence = "independent")
#' ))
#'
#' # The signed spectrum keeps its negative root; both regions have one.
#' geometry_spectrum(geometry)$values
#'
#' # The latent layer truncates it and says what that cost, per region.
#' latent <- latent_geometry(geometry)
#' latent
#' latent$cumulative
#' latent$projection
#'
#' # A contracted readout has no spectrum to project, and clamping it would be
#' # a per-node total clamp: a different operator moving different mass.
#' effect <- contrast_energy(geometry, c(face = 1, house = -1, tool = 0))
#' catch_refusal(latent_geometry(effect))$capability
#' @export
latent_geometry <- function(x,
                            method = c("psd_projection"),
                            component = c("total", "coherent",
                              "configuration"),
                            row_block = 1024L) {
  method <- .latent_projection_method(method)
  component_given <- !missing(component)
  component <- match.arg(component)
  row_block <- .validate_tile_size(row_block, "row_block")

  source <- .latent_source(x, component, component_given, row_block)
  projection <- .latent_psd_projection(source$spectrum)

  scientific_plan_id <- .latent_scientific_id(
    source$receipt$scientific_plan_id, method, source$component
  )
  receipt <- .latent_execution_receipt(
    source$receipt, scientific_plan_id, method
  )
  projection_receipt <- .new_latent_projection_receipt(
    method, source$component, projection, source$kind,
    source$receipt$scientific_plan_id, scientific_plan_id
  )
  .new_effect_latent_geometry(
    projection, source$component, method, source$index, projection_receipt,
    receipt, list(projected_from = source$receipt$scientific_plan_id)
  )
}

# Printing -------------------------------------------------------------------

.latent_moved_items <- function(x) {
  moved <- x$moved_mass
  clipped <- sum(moved > 0)
  if (!clipped) {
    return(c("none", "no source root was negative"))
  }
  c(
    paste0(.pf_num(sum(moved)), " moved"),
    sprintf("%d of %s clipped", clipped,
      .msg_count(length(moved), "measurement")),
    paste0("max share ", .pf_num(max(x$moved_share), 3L)),
    if (clipped < length(moved)) sprintf("%d moved none", length(moved) - clipped)
  )
}

.latent_n_eff_line <- function(n_eff) {
  finite <- n_eff[is.finite(n_eff)]
  if (!length(finite)) {
    return("undefined; no measurement kept positive mass")
  }
  masked <- length(n_eff) - length(finite)
  paste0(
    sprintf("median %s (range %s to %s)",
      .pf_num(stats::median(finite), 3L), .pf_num(min(finite), 3L),
      .pf_num(max(finite), 3L)),
    if (masked) sprintf(", %s masked", .msg_count(masked, "measurement")) else ""
  )
}

#' @export
as.data.frame.effect_latent_geometry <- function(x, row.names = NULL,
                                                 optional = FALSE, ...) {
  .validate_effect_latent_geometry(x)
  .bind_result_values(x$index, cbind(
    n_eff = x$n_eff, moved_mass = x$moved_mass, moved_share = x$moved_share,
    x$spectrum
  ))
}

#' @export
format.effect_latent_geometry <- function(x, ...) {
  .format_counted_result("effect_latent_geometry", nrow(x$spectrum),
    paste0(ncol(x$spectrum), " nonnegative roots"))
}

#' @export
print.effect_latent_geometry <- function(x, ...) {
  .validate_effect_latent_geometry(x)
  .format_result_preview(x, "effect_latent_geometry", fields = list(
    component = x$component,
    projection = sprintf("%s (%s)", x$method,
      .latent_operator_phrase(x$method)),
    moved_mass = .pf_wrap(.latent_moved_items(x)),
    n_eff = .latent_n_eff_line(x$n_eff),
    reading = .latent_reading_line
  ), ...)
  cat(sprintf("  %-14s%s\n", "next:", "x$cumulative, x$projection"))
  invisible(x)
}

#' @export
format.effect_latent_projection_receipt <- function(x, ...) {
  .pf_inline("effect_latent_projection_receipt", x$method,
    sprintf("%d of %s clipped", x$clipped,
      .msg_count(x$measurements, "measurement")))
}

#' @export
print.effect_latent_projection_receipt <- function(x, ...) {
  .validate_latent_projection_receipt(x)
  .pf_emit("effect_latent_projection_receipt", list(
    method = x$method,
    operator = .latent_operator_phrase(x$method),
    component = x$component,
    clipped = sprintf("%d of %s", x$clipped,
      .msg_count(x$measurements, "measurement")),
    moved_mass = sprintf("%s total, max share %s",
      .pf_num(x$total_moved_mass), .pf_num(x$max_moved_share, 3L)),
    masked = if (x$masked) .msg_count(x$masked, "measurement") else "none",
    source = .pf_sig(x$source_scientific_plan_id)
  ))
  invisible(x)
}
