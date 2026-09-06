# Location transports --------------------------------------------------------
#
# Layer 2 (values). A transport is a declaration: it says how one participant's
# *native* conservative frame nodes are carried onto a shared set of *group*
# nodes, and it is an input to crossform, never something crossform estimates
# (`design/population-form-contract.md` §9.2). Nothing here reads data, opens a
# source, or fits anything; the file calls only layer-1 guards and the
# provenance walker in `effect-space.R`.
#
# The normative shape is `population-form-v1` §1: a sparse nonnegative
# `n_native x (m + 1)` operator whose last column is the sink, row-stochastic
# *including* the sink. The sink is what turns partial coverage into a visible
# number instead of a silent loss, so it is materialized even when it is
# identically zero, and the constructor computes it rather than trusting a
# caller to have appended one --- a supplied "sink" that did not close the row
# sums would defeat the whole point of the column.

.transport_sink_label <- "<sink>"

#' Declare a location transport onto a shared group node set
#'
#' A location transport carries each participant's native conservative frame
#' nodes onto a common set of group nodes. It is a **declaration**, not an
#' estimate: `crossform` accepts a typed, sparse, provenance-bearing transport
#' and evaluates it, and refuses to learn one (registration and
#' functional-transport learning stay outside the package).
#'
#' @section The operator:
#' `matrix` supplies only the `n_native` by `m` **group** block. The
#' constructor appends the required sink column itself, with row `x` receiving
#' sink mass `1 - sum_j matrix[x, j]`, so that the assembled
#' `n_native x (m + 1)` operator is row-stochastic by construction. A row whose
#' group mass exceeds one would need negative sink mass and is refused. The
#' sink column always exists, even when it is identically zero: unmapped
#' territory has to be a number a reader can see, not an absence.
#'
#' @section Budget and density:
#' `semantics` has no default, because the two choices estimate different
#' things and the difference is not a display option.
#'
#' * `"budget"` transports the signed sum: group value `j` is
#'   `sum_x P[x, j] * values[x]`. Each participant's transported total over the
#'   group nodes equals its native total minus the sink mass, exactly.
#' * `"density"` divides that by the transported row mass,
#'   `sum_x P[x, j] * row_mass[x]`, giving transported budget per unit
#'   transported territory. A group node reached by no native mass is `NA`,
#'   never `0`. The sink is reported in budget units under either semantics,
#'   because a "density" of unmapped territory has no referent.
#'
#' `row_mass` is the declared positive territory measure of each native row.
#' Its default is the unit vector --- density is then budget per native node
#' count --- and it is recorded either way, because it is the denominator of
#' the estimand.
#'
#' @section Provenance:
#' `provenance` must be a list carrying `method`, one of `"anatomical"`
#' (built from a registration or parcellation that never saw the response
#' data), `"functional"` (built from data), or `"external"` (built outside
#' `crossform` by a procedure the caller names), and `details`, one string
#' saying how the operator was built. A `"functional"` transport additionally
#' requires `cross_fit`: the runs, sessions or tasks whose data built it. That
#' field is not advisory. A circular transport --- one fitted on a partition
#' that is also used to evaluate it --- reports a transport gain roughly three
#' times the honest one and is otherwise indistinguishable from it, so the
#' field that would let a plan exclude those partitions is required rather
#' than encouraged. `fitting_sample` and `cross_fit_folds` may distinguish the
#' operator's full fitting sample from the named fold assignments; for
#' compatibility they default to the required `cross_fit` declaration.
#'
#' Cross-fitting is a circularity control, not uncertainty propagation. Every
#' transport receives a derived `$provenance$conditioning` record naming its
#' source, fixed-versus-estimated status, fitting sample, folds, conditional
#' inference scope, excluded uncertainty, and the future propagation
#' capability. Current population inference is always conditional on the
#' realized operator: `uncertainty_propagated` and
#' `marginal_over_transport` are immutable `FALSE`. Any further non-reserved
#' keys are kept as declared.
#'
#' @param matrix The `n_native` by `m` nonnegative group block, as a base
#'   matrix or a `Matrix`. Row sums may be at most one; the shortfall becomes
#'   sink mass.
#' @param native_index One row per native node: a character vector of node
#'   identifiers, or a data frame with a `node` column plus any per-row frame
#'   metadata (`family`, `scale`, `center`, `alpha` are type-checked when
#'   present and carried through unchanged).
#' @param group_index One row per group node: a character vector of node
#'   identifiers, or a data frame with a `node` column plus any group-node
#'   metadata, such as the `coord1`, `coord2`, ... centres the displacement
#'   diagnostics need. It describes the `m` group nodes only; the sink is
#'   column `m + 1` of the operator and has no index row.
#' @param semantics `"budget"` or `"density"`. Required; there is no default.
#' @param provenance A list carrying at least `method` and `details`, plus
#'   `cross_fit` when `method` is `"functional"`. Optional `fitting_sample`
#'   and `cross_fit_folds` refine that declaration. Conditioning fields are
#'   derived and cannot be supplied by the caller.
#' @param row_mass Optional declared positive territory measure, one entry per
#'   native node. Defaults to the unit vector.
#' @param tolerance Positive row-sum tolerance. Group mass above `1 +
#'   tolerance` on any row is refused rather than renormalized.
#' @return An `effect_location_transport` carrying `$matrix` (the assembled
#'   sparse operator, sink included), `$native_index`, `$group_index`,
#'   `$semantics`, `$row_mass`, `$provenance`, and a content-addressed
#'   `$signature` that binds all of them --- including the operator's own
#'   contents --- so the transport can enter a scientific plan identity.
#' @seealso [anatomical_transport()] and [external_transport()] for the two
#'   convenience constructors, and [transport_values()] for applying one.
#' @family population transports
#'
#' @section Refusal:
#' A `"functional"` transport declared without `cross_fit` signals an
#' `effect_capability_refusal` carrying capability `"cross_fit_provenance"` in
#' namespace `"location_transport"`. Branch on it with [catch_refusal()]
#' rather than on the message text.
#'
#' @examples
#' # Four native nodes onto two group nodes. Node 3 splits its mass; node 4
#' # keeps 30% of its territory unmapped, so 0.3 lands in the sink.
#' block <- rbind(c(1, 0), c(1, 0), c(0.6, 0.4), c(0, 0.7))
#' transport <- location_transport(
#'   block,
#'   native_index = c("v1", "v2", "v3", "v4"),
#'   group_index = c("left", "right"),
#'   semantics = "budget",
#'   provenance = list(method = "anatomical", details = "hand-declared")
#' )
#' transport
#'
#' # The sink column is materialized and closes every row sum to one.
#' as.matrix(transport$matrix)
#'
#' # Budget preservation: the transported total, sink included, is the native
#' # total, for signed ledgers as well as nonnegative ones.
#' ledger <- c(1.5, -0.5, 2, -1)
#' carried <- transport_values(transport, ledger)
#' carried
#' abs(sum(carried) - sum(ledger)) < 1e-12
#'
#' # A row that claims more than unit mass would need a negative sink, and is
#' # refused rather than renormalized.
#' bad <- try(location_transport(
#'   rbind(c(0.8, 0.4)), "v1", c("left", "right"),
#'   semantics = "budget",
#'   provenance = list(method = "anatomical", details = "hand-declared")
#' ), silent = TRUE)
#' conditionMessage(attr(bad, "condition"))
#'
#' # A data-derived transport must name the data that built it.
#' refusal <- catch_refusal(location_transport(
#'   block, c("v1", "v2", "v3", "v4"), c("left", "right"),
#'   semantics = "budget",
#'   provenance = list(method = "functional", details = "clustered responses")
#' ))
#' refusal$capability
#' refusal$reasons
#' @export
location_transport <- function(matrix, native_index, group_index, semantics,
                               provenance, row_mass = NULL,
                               tolerance = 1e-9) {
  if (missing(semantics)) {
    .input_error(paste0(
      "`semantics` must be declared as \"budget\" or \"density\". The two are ",
      "different estimands built from the same operator --- budget transports ",
      "the signed sum, density divides it by the transported row mass --- so ",
      "there is no default to fall back on."
    ), arg = "semantics", expected = "\"budget\" or \"density\"")
  }
  if (missing(provenance)) {
    .input_error(paste0(
      "`provenance` must declare how the transport was built. A transport ",
      "without provenance cannot be told apart from a circular one."
    ), arg = "provenance", expected = "a list with `method` and `details`")
  }
  semantics <- .transport_semantics(semantics)
  tolerance <- .check_number(tolerance, "tolerance", positive = TRUE)
  provenance <- .transport_provenance(provenance)

  group <- .transport_sparse(matrix, "matrix")
  n_native <- nrow(group)
  n_group <- ncol(group)
  if (n_native < 1L || n_group < 1L) {
    .input_error(paste0(
      "`matrix` must have at least one native row and one group column; ",
      "received a ", n_native, " x ", n_group, " operator."
    ), arg = "matrix")
  }
  native_index <- .transport_index(native_index, "native_index", n_native)
  group_index <- .transport_index(group_index, "group_index", n_group)
  row_mass <- .transport_row_mass(row_mass, n_native)

  if (length(group@x) && any(group@x < 0)) {
    worst <- which.min(group@x)
    .input_error(sprintf(paste0(
      "A transport must be nonnegative, but `matrix` holds %s. A negative ",
      "entry is not a small transport --- it moves evidence the wrong way and ",
      "breaks the row-stochastic law the sink column closes."
    ), format(group@x[[worst]])), arg = "matrix",
      received = sprintf("minimum entry %s", format(min(group@x))),
      expected = "every entry at least zero")
  }

  deficit <- 1 - as.numeric(Matrix::rowSums(group))
  if (any(deficit < -tolerance)) {
    worst <- which.min(deficit)
    .input_error(sprintf(paste0(
      "Native row %d carries group mass %s, which is more than one, so its ",
      "sink mass would be negative. Rows are stochastic including the sink, ",
      "and crossform refuses to renormalize a row rather than silently ",
      "changing what the transport declares."
    ), worst, format(1 - deficit[[worst]])), arg = "matrix",
      received = sprintf("worst row mass %s", format(max(1 - deficit))),
      expected = "row mass at most one")
  }
  deficit[deficit < 0] <- 0

  nonzero <- which(deficit > 0)
  sink_column <- Matrix::sparseMatrix(
    i = nonzero, j = rep(1L, length(nonzero)), x = deficit[nonzero],
    dims = c(n_native, 1L)
  )
  operator <- Matrix::drop0(cbind(group, sink_column))
  row_error <- max(abs(as.numeric(Matrix::rowSums(operator)) - 1))
  if (!is.finite(row_error) || row_error > tolerance) {
    .input_error(sprintf(paste0(
      "Transport rows must sum to one within `tolerance`; the worst row is ",
      "off by %s."
    ), format(row_error)), arg = "matrix",
      received = sprintf("worst row deviation %s", format(row_error)),
      expected = sprintf("at most %s", format(tolerance)))
  }

  .transport_seal(operator, native_index, group_index, semantics, row_mass,
    provenance)
}

# The sealed record and its identity. The operator's own contents enter the
# signature through a digest of its triplet form: two transports that differ in
# one entry are two different estimands, and a signature that hashed only the
# labels would call them the same.
.transport_seal <- function(operator, native_index, group_index, semantics,
                            row_mass, provenance) {
  semantic <- list(
    schema_version = 1L,
    operator = .transport_operator_digest(operator),
    native_index = native_index,
    group_index = group_index,
    semantics = semantics,
    row_mass = row_mass,
    provenance = provenance
  )
  structure(list(
    matrix = operator,
    native_index = native_index,
    group_index = group_index,
    semantics = semantics,
    row_mass = row_mass,
    provenance = provenance,
    signature = .sha256_signature(semantic)
  ), class = "effect_location_transport")
}

.transport_operator_digest <- function(operator) {
  triplet <- Matrix::summary(Matrix::drop0(operator))
  .sha256_signature(list(
    dim = as.integer(dim(operator)),
    i = as.integer(triplet$i),
    j = as.integer(triplet$j),
    x = as.numeric(triplet$x)
  ), prefix = "sha256:")
}

.validate_location_transport <- function(x) {
  expected <- c("matrix", "native_index", "group_index", "semantics",
    "row_mass", "provenance", "signature")
  if (!.sealed_fields(x, "effect_location_transport", expected) ||
      !methods::is(x$matrix, "sparseMatrix") ||
      !is.data.frame(x$native_index) || !is.data.frame(x$group_index) ||
      !.is_string(x$semantics) || !is.numeric(x$row_mass) ||
      !is.list(x$provenance)) {
    .input_error("Location-transport fields are missing or noncanonical.")
  }
  provenance <- x$provenance
  if (!.is_string(provenance$method) ||
      !provenance$method %in% c("anatomical", "functional", "external") ||
      !.is_string(provenance$details) ||
      (identical(provenance$method, "functional") &&
        !.is_strings(provenance[["cross_fit", exact = TRUE]]))) {
    .contract_error("Location-transport provenance is noncanonical.")
  }
  expected_conditioning <- .transport_conditioning(
    provenance$method, provenance$details, provenance$fitting_sample,
    provenance$cross_fit_folds
  )
  if (!is.character(provenance$fitting_sample) ||
      !is.character(provenance$cross_fit_folds) ||
      !identical(provenance$conditioning, expected_conditioning)) {
    .contract_error(paste0(
      "Location-transport conditioning is missing or inconsistent. ",
      "Transport inference must remain explicitly conditional on the ",
      "realized operator."
    ))
  }
  n_native <- nrow(x$matrix)
  if (nrow(x$native_index) != n_native ||
      nrow(x$group_index) != ncol(x$matrix) - 1L ||
      length(x$row_mass) != n_native) {
    .contract_error(paste0(
      "Location-transport index tables do not match the operator: the ",
      "operator is ", n_native, " x ", ncol(x$matrix),
      " (group nodes plus the sink) against ", nrow(x$native_index),
      " native and ", nrow(x$group_index), " group rows."
    ))
  }
  .check_signature(
    x$signature,
    .sha256_signature(list(
      schema_version = 1L,
      operator = .transport_operator_digest(x$matrix),
      native_index = x$native_index,
      group_index = x$group_index,
      semantics = x$semantics,
      row_mass = x$row_mass,
      provenance = x$provenance
    )),
    "Location-transport identity is inconsistent."
  )
  invisible(x)
}

#' Carry native node values onto the group nodes
#'
#' `transport_values()` applies a [location_transport()] to one value per
#' native node --- a queried ledger, or a matrix of packed forms with one row
#' per native node --- and returns the group-node values with the sink last.
#'
#' Under `"budget"` semantics the result is `t(P) %*% values`, and the
#' transported total including the sink equals the native total exactly. The
#' argument is `P 1 = 1` and never `values >= 0`, so this is a signed-sum law:
#' it holds for the signed, crossvalidated ledgers a conservative geometry
#' actually produces. Under `"density"` semantics each group entry is divided
#' by the transported row mass, group nodes reached by no native mass come back
#' `NA` rather than `0`, and the sink row stays in budget units.
#'
#' @param transport A [location_transport()].
#' @param values One finite value per native node: a numeric vector, or a
#'   matrix with one row per native node.
#' @return A named numeric vector of length `m + 1` when `values` is a vector,
#'   otherwise an `(m + 1)` by `ncol(values)` matrix. The last entry or row is
#'   the sink, labelled `<sink>`.
#' @seealso [location_transport()].
#' @family population transports
#' @examples
#' transport <- anatomical_transport(
#'   native_coords = cbind(c(0, 1, 4, 5, 9)),
#'   group_coords = cbind(c(0, 5)),
#'   semantics = "budget", radius = 2
#' )
#'
#' # Native node 5 sits further than the radius from either group centre, so
#' # its whole budget lands in the sink and the total still closes.
#' ledger <- c(2, -1, 0.5, 1.5, 3)
#' transport_values(transport, ledger)
#' sum(transport_values(transport, ledger)) - sum(ledger)
#'
#' # A matrix of packed forms transports column by column.
#' forms <- cbind(a = ledger, b = rev(ledger))
#' transport_values(transport, forms)
#' @export
transport_values <- function(transport, values) {
  transport <- .validate_location_transport(transport)
  operator <- transport$matrix
  n_native <- nrow(operator)
  n_group <- ncol(operator) - 1L

  vector_input <- is.null(dim(values))
  block <- if (vector_input) {
    if (!is.numeric(values) || length(values) != n_native) {
      .check_failed("values",
        sprintf("a numeric vector with %s", .msg_count(n_native, "entry",
          "entries")), values)
    }
    cbind(as.numeric(values))
  } else {
    given <- as.matrix(values)
    if (!is.numeric(given) || nrow(given) != n_native) {
      .check_failed("values",
        sprintf("a numeric matrix with %s", .msg_count(n_native, "row")),
        values)
    }
    given
  }
  if (!all(is.finite(block))) {
    .input_error(paste0(
      "`values` must be finite. A missing native value would spread over ",
      "every group node the row touches, so the transport refuses it rather ",
      "than producing a ledger whose absences are invisible."
    ), arg = "values")
  }

  carried <- as.matrix(Matrix::crossprod(operator, block))
  if (identical(transport$semantics, "density")) {
    mass <- as.numeric(Matrix::crossprod(operator, transport$row_mass))
    rows <- seq_len(n_group)
    empty <- mass[rows] <= 0
    scale <- mass[rows]
    scale[empty] <- 1
    carried[rows, ] <- carried[rows, , drop = FALSE] / scale
    carried[rows[empty], ] <- NA_real_
  }
  labels <- c(as.character(transport$group_index$node), .transport_sink_label)
  rownames(carried) <- labels
  if (vector_input) {
    return(stats::setNames(as.numeric(carried), labels))
  }
  colnames(carried) <- colnames(block)
  carried
}

#' Build a hard nearest-centre transport from node coordinates
#'
#' `anatomical_transport()` is the transport a registration or a parcellation
#' induces: each native node is assigned whole to one group node, and nothing
#' is estimated from the response data.
#'
#' @section The assignment rule:
#' Every native node is assigned, with mass exactly one, to the group node
#' whose centre is closest in Euclidean distance in the coordinate system both
#' matrices are already expressed in. No registration is performed and no
#' image is resampled: the transport maps nodes, not voxels. Ties go to the
#' lowest group-node position, so the operator is a deterministic function of
#' its inputs. When `radius` is supplied and the closest group centre is
#' further away than `radius`, the native node is left unassigned and its
#' entire mass goes to the sink; without `radius` every native node is
#' assigned and the sink is identically zero but still materialized.
#'
#' The result is a hard assignment, so every row has entropy zero. A partial
#' volume warp does not, which is why the population diagnostics report row
#' entropy rather than assuming it.
#'
#' @param native_coords Native node centres, one row per native node.
#' @param group_coords Group node centres, one row per group node, with the
#'   same number of coordinate columns as `native_coords`.
#' @param semantics `"budget"` or `"density"`. Required; see
#'   [location_transport()].
#' @param radius Optional positive assignment radius. A native node whose
#'   closest group centre is further than this goes entirely to the sink.
#' @param native_index,group_index Optional index tables or node-identifier
#'   vectors, as in [location_transport()]. The coordinates are appended to
#'   them as `coord1`, `coord2`, ... columns, so the displacement diagnostics
#'   can be recomputed from the transport alone.
#' @param provenance Optional extra provenance entries, such as the atlas or
#'   warp identity. `method` is fixed at `"anatomical"` and `details` records
#'   the assignment rule unless the caller supplies its own.
#' @param row_mass Optional declared positive territory measure, one entry per
#'   native node.
#' @return An `effect_location_transport`.
#' @seealso [location_transport()] for the general constructor and the field
#'   semantics.
#' @family population transports
#' @examples
#' # Five native nodes on a line, two group centres at 0 and 5, radius 2.
#' transport <- anatomical_transport(
#'   native_coords = cbind(c(0, 1, 4, 5, 9)),
#'   group_coords = cbind(c(0, 5)),
#'   semantics = "budget", radius = 2,
#'   provenance = list(atlas = "example:line-atlas:v1")
#' )
#' as.matrix(transport$matrix)
#'
#' # The unassigned node is visible as sink territory, not as absence.
#' transport
#' @export
anatomical_transport <- function(native_coords, group_coords, semantics,
                                 radius = NULL, native_index = NULL,
                                 group_index = NULL, provenance = list(),
                                 row_mass = NULL) {
  if (missing(semantics)) {
    .input_error(paste0(
      "`semantics` must be declared as \"budget\" or \"density\"; an ",
      "anatomical transport is as much an estimand choice as any other."
    ), arg = "semantics", expected = "\"budget\" or \"density\"")
  }
  native_coords <- .transport_coords(native_coords, "native_coords")
  group_coords <- .transport_coords(group_coords, "group_coords")
  if (ncol(native_coords) != ncol(group_coords)) {
    .input_error(sprintf(paste0(
      "`native_coords` and `group_coords` must live in the same coordinate ",
      "system, but they have %d and %d columns. Transport maps nodes that ",
      "are already registered; it does not register them."
    ), ncol(native_coords), ncol(group_coords)), arg = "group_coords")
  }
  if (!is.null(radius)) {
    radius <- .check_number(radius, "radius", positive = TRUE)
  }

  n_native <- nrow(native_coords)
  n_group <- nrow(group_coords)
  nearest <- .transport_nearest(native_coords, group_coords)
  assigned <- if (is.null(radius)) {
    rep(TRUE, n_native)
  } else {
    nearest$distance <= radius
  }
  rows <- which(assigned)
  block <- Matrix::sparseMatrix(
    i = rows, j = nearest$group[rows], x = rep(1, length(rows)),
    dims = c(n_native, n_group)
  )

  native_index <- .transport_with_coords(
    .transport_index(native_index, "native_index", n_native,
      default_prefix = "native"),
    native_coords, "native_index"
  )
  group_index <- .transport_with_coords(
    .transport_index(group_index, "group_index", n_group,
      default_prefix = "group"),
    group_coords, "group_index"
  )

  details <- if (is.null(radius)) {
    "nearest group centre, ties to the lowest group position"
  } else {
    sprintf(paste0("nearest group centre within radius %s, ties to the lowest ",
      "group position; otherwise sink"), format(radius))
  }
  provenance <- .transport_fixed_method(provenance, "anatomical", details,
    list(rule = "nearest_centre",
      radius = if (is.null(radius)) NA_real_ else as.numeric(radius),
      coordinates = ncol(native_coords)))

  location_transport(block, native_index, group_index, semantics = semantics,
    provenance = provenance, row_mass = row_mass)
}

#' Admit a transport built outside crossform
#'
#' `crossform` does not learn transports. `external_transport()` is the door a
#' transport built elsewhere --- by a registration pipeline, a parcellation
#' tool, or a functional alignment method --- comes through: a thin wrapper
#' over [location_transport()] that takes the group block as it stands, reads
#' node identifiers from its dimension names when they are there, and requires
#' the provenance that makes the operator auditable.
#'
#' Declaring `method = "external"` says only that the operator was built
#' outside this package. A transport that was fitted to response data is
#' `"functional"` whichever program fitted it, and must name its `cross_fit`
#' partitions; `"external"` is not a way around that.
#'
#' @param P The `n_native` by `m` nonnegative group block, as a base matrix or
#'   a `Matrix`. Dimension names, when present, become the default node
#'   identifiers.
#' @param semantics `"budget"` or `"density"`. Required; see
#'   [location_transport()].
#' @param provenance A list carrying at least `details`; `method` defaults to
#'   `"external"`.
#' @param native_index,group_index Optional index tables or node-identifier
#'   vectors overriding the dimension names.
#' @param row_mass Optional declared positive territory measure, one entry per
#'   native node.
#' @param tolerance Positive row-sum tolerance, as in [location_transport()].
#' @return An `effect_location_transport`.
#' @seealso [location_transport()], [anatomical_transport()].
#' @family population transports
#' @examples
#' # A partial-volume operator produced by some other program, with its node
#' # names already on the dimensions.
#' P <- matrix(c(0.7, 0.2, 0, 0.3, 0.8, 0.5), nrow = 3,
#'   dimnames = list(c("n1", "n2", "n3"), c("A", "B")))
#' transport <- external_transport(
#'   P, semantics = "density",
#'   provenance = list(details = "partial-volume warp from atlas-tool 2.1")
#' )
#' transport
#'
#' # Density divides the transported budget by the transported row mass; with
#' # the default unit row mass that is the transported node count.
#' transport_values(transport, c(1, 1, 1))
#' @export
external_transport <- function(P, semantics, provenance, native_index = NULL,
                               group_index = NULL, row_mass = NULL,
                               tolerance = 1e-9) {
  if (missing(semantics)) {
    .input_error(paste0(
      "`semantics` must be declared as \"budget\" or \"density\"; crossform ",
      "does not guess the estimand of a transport it did not build."
    ), arg = "semantics", expected = "\"budget\" or \"density\"")
  }
  if (missing(provenance)) {
    .input_error(paste0(
      "`provenance` must say how the external transport was built. An ",
      "operator admitted without provenance cannot be told apart from a ",
      "circular one."
    ), arg = "provenance", expected = "a list with `details`")
  }
  if (!is.list(provenance)) {
    .input_error("`provenance` must be a list.", arg = "provenance")
  }
  if (is.null(provenance$method)) provenance$method <- "external"

  names_of <- dimnames(P)
  if (is.null(native_index) && !is.null(names_of) &&
      .is_strings(names_of[[1L]])) {
    native_index <- names_of[[1L]]
  }
  if (is.null(group_index) && !is.null(names_of) && length(names_of) > 1L &&
      .is_strings(names_of[[2L]])) {
    group_index <- names_of[[2L]]
  }
  block <- .transport_sparse(P, "P")
  if (is.null(native_index)) {
    native_index <- .transport_default_nodes(nrow(block), "native")
  }
  if (is.null(group_index)) {
    group_index <- .transport_default_nodes(ncol(block), "group")
  }
  location_transport(block, native_index, group_index, semantics = semantics,
    provenance = provenance, row_mass = row_mass, tolerance = tolerance)
}

# Helpers ---------------------------------------------------------------------

.transport_semantics <- function(semantics) {
  if (!.is_string(semantics) || !semantics %in% c("budget", "density")) {
    .check_failed("semantics", "either \"budget\" or \"density\"", semantics)
  }
  semantics
}

.transport_sparse <- function(x, arg) {
  if (is.data.frame(x)) x <- as.matrix(x)
  dense <- is.matrix(x) && is.numeric(x)
  if (!dense && !methods::is(x, "Matrix")) {
    .check_failed(arg, "a numeric matrix or a Matrix object", x)
  }
  if (dense) x <- Matrix::Matrix(x, sparse = TRUE)
  x <- methods::as(
    methods::as(methods::as(x, "dMatrix"), "generalMatrix"), "CsparseMatrix"
  )
  if (length(x@x) && !all(is.finite(x@x))) {
    .input_error(sprintf(paste0(
      "`%s` must have no missing or infinite entries. A transport is a ",
      "declaration of where evidence goes, and a missing weight declares ",
      "nothing."
    ), arg), arg = arg)
  }
  Matrix::drop0(x)
}

.transport_default_nodes <- function(n, prefix) {
  paste0(prefix, seq_len(n))
}

# Index tables are typed but open: `node` is required and the four per-row
# frame metadata columns are checked when present, because a transport maps
# *named* rows and a result that cannot say which instrument produced a number
# is not auditable. Anything else a caller carries is kept as declared.
.transport_index <- function(index, arg, n, default_prefix = NULL) {
  if (is.null(index) && !is.null(default_prefix)) {
    index <- .transport_default_nodes(n, default_prefix)
  }
  if (is.character(index) || is.factor(index)) {
    index <- data.frame(node = as.character(index), stringsAsFactors = FALSE)
  }
  if (!is.data.frame(index)) {
    .check_failed(arg, paste0("a data frame with a `node` column, or a ",
      "character vector of node identifiers"), index)
  }
  index <- as.data.frame(index, stringsAsFactors = FALSE)
  if (!"node" %in% names(index)) {
    .input_error(sprintf(paste0(
      "`%s` must carry a `node` column naming each row. Transport maps named ",
      "nodes; positions alone cannot survive a group result."
    ), arg), arg = arg, expected = "a `node` column")
  }
  index$node <- as.character(index$node)
  if (!.is_strings(index$node, unique = TRUE)) {
    .input_error(sprintf(
      "`%s$node` must be nonempty and unique identifiers.", arg
    ), arg = arg)
  }
  if (nrow(index) != n) {
    .input_error(sprintf(
      "`%s` must have one row per node; received %d for %s.",
      arg, nrow(index), .msg_count(n, "node")
    ), arg = arg, received = .msg_count(nrow(index), "row"),
      expected = .msg_count(n, "row"))
  }
  for (column in c("family", "center")) {
    if (column %in% names(index)) {
      index[[column]] <- as.character(index[[column]])
    }
  }
  for (column in c("scale", "alpha")) {
    if (column %in% names(index) && !is.numeric(index[[column]])) {
      .input_error(sprintf(
        "`%s$%s` must be numeric; it is the row's own %s.", arg, column,
        if (identical(column, "scale")) "scale parameter" else "family weight"
      ), arg = arg)
    }
  }
  if (any(vapply(index, is.list, logical(1)))) {
    .input_error(sprintf(
      "`%s` must hold atomic columns; a list column is not portable.", arg
    ), arg = arg)
  }
  rownames(index) <- NULL
  index
}

.transport_row_mass <- function(row_mass, n) {
  if (is.null(row_mass)) {
    return(rep(1, n))
  }
  if (!is.numeric(row_mass) || length(row_mass) != n ||
      !all(is.finite(row_mass)) || any(row_mass <= 0)) {
    .check_failed("row_mass",
      sprintf(paste0("%s, every one positive and finite --- it is the ",
        "denominator of density semantics"), .msg_count(n, "declared mass",
        "declared masses")), row_mass)
  }
  as.numeric(row_mass)
}

# Transport uncertainty is deliberately a future capability, not a label a
# caller can assert. Cross-fitting limits reuse of the response data, but it
# does not integrate over the fitted operator or its fold assignment. This
# canonical record follows the operator into every population receipt.
.transport_conditioning <- function(method, details, fitting_sample,
                                    cross_fit_folds) {
  estimated <- identical(method, "functional")
  list(
    source = details,
    operator_status = if (estimated) "estimated" else "fixed",
    fitting_sample = as.character(fitting_sample),
    cross_fit_folds = as.character(cross_fit_folds),
    circularity_control = if (estimated) {
      "cross_fit_partitions_declared"
    } else {
      "fixed_operator"
    },
    inference_scope = "conditional_on_realized_transport",
    uncertainty_propagated = FALSE,
    marginal_over_transport = FALSE,
    excluded_uncertainty = if (estimated) {
      c("transport_operator_estimation", "cross_fit_fold_assignment")
    } else {
      "transport_source_choice"
    },
    future = list(
      capability = "transport_uncertainty_propagation",
      status = "not_implemented",
      requires = c(
        "transport_sampling_law",
        "joint_transport_response_resampling",
        "validated_propagation_operator"
      )
    )
  )
}

# `method` and `details` lead the record so that two logically identical
# provenances hash the same regardless of the order they were typed in; the
# caller's remaining keys keep the order they were given.
.transport_provenance <- function(provenance) {
  if (!is.list(provenance)) {
    .input_error("`provenance` must be a list.", arg = "provenance")
  }
  .validate_effect_provenance(provenance, "provenance")
  method <- provenance$method
  if (!.is_string(method) ||
      !method %in% c("anatomical", "functional", "external")) {
    .check_failed("provenance$method",
      "one of \"anatomical\", \"functional\", or \"external\"", method)
  }
  details <- provenance$details
  if (!.is_string(details)) {
    .input_error(paste0(
      "`provenance$details` must be one nonempty string saying how the ",
      "operator was built --- the warp, the atlas, or the alignment method. ",
      "The transport is part of the estimand, so how it was made belongs in ",
      "its identity."
    ), arg = "provenance$details", received = .msg_value(details),
      expected = "one nonempty character string")
  }
  reserved <- intersect(names(provenance), c(
    "conditioning", "uncertainty_propagated", "marginal_over_transport"
  ))
  if (length(reserved)) {
    .input_error(paste0(
      "Transport-conditioning fields are derived by crossform and cannot ",
      "be supplied in provenance. Cross-fitting does not make inference ",
      "marginal over transport estimation. Remove ", .msg_names(reserved),
      "; the sealed record will declare the supported conditional scope."
    ), arg = paste0("provenance$", reserved[[1L]]))
  }
  cross_fit <- provenance[["cross_fit", exact = TRUE]]
  if (identical(method, "functional")) {
    if (!.is_strings(cross_fit)) {
      .capability_refusal(paste0(
        "A functional transport must name the runs, sessions, or tasks whose ",
        "data built it, in `provenance$cross_fit`. A transport fitted on a ",
        "partition that is also used to evaluate it reports roughly three ",
        "times the honest transport gain and lands at the oracle ceiling, and ",
        "in a result object the two numbers are indistinguishable. The field ",
        "is what lets a population plan exclude those partitions."
      ),
        capability = "cross_fit_provenance",
        namespace = "location_transport",
        reasons = "cross_fit_partitions_not_declared",
        remedies = paste0(
          "Pass `provenance = list(method = \"functional\", details = ..., ",
          "cross_fit = c(\"run-1\", \"run-2\"))` naming every partition that ",
          "went into fitting the operator, or declare `method = ",
          "\"anatomical\"` if the response data was never seen."
        ))
    }
    cross_fit <- as.character(cross_fit)
  } else if (!is.null(cross_fit)) {
    if (!.is_strings(cross_fit)) {
      .input_error(paste0(
        "`provenance$cross_fit`, when given, must name partitions as nonempty ",
        "strings."
      ), arg = "provenance$cross_fit")
    }
    cross_fit <- as.character(cross_fit)
  }
  fitting_sample <- provenance[["fitting_sample", exact = TRUE]]
  if (is.null(fitting_sample) && identical(method, "functional")) {
    fitting_sample <- cross_fit
  }
  if (is.null(fitting_sample)) fitting_sample <- character()
  if (!is.character(fitting_sample) || anyNA(fitting_sample) ||
      any(!nzchar(fitting_sample))) {
    .input_error(paste0(
      "`provenance$fitting_sample`, when given, must identify the sample ",
      "that built the operator using nonempty strings."
    ), arg = "provenance$fitting_sample")
  }
  cross_fit_folds <- provenance[["cross_fit_folds", exact = TRUE]]
  if (is.null(cross_fit_folds) && identical(method, "functional")) {
    cross_fit_folds <- cross_fit
  }
  if (is.null(cross_fit_folds)) cross_fit_folds <- character()
  if (!is.character(cross_fit_folds) || anyNA(cross_fit_folds) ||
      any(!nzchar(cross_fit_folds))) {
    .input_error(paste0(
      "`provenance$cross_fit_folds`, when given, must name each declared ",
      "cross-fitting fold using nonempty strings."
    ), arg = "provenance$cross_fit_folds")
  }
  conditioning <- .transport_conditioning(
    method, details, fitting_sample, cross_fit_folds
  )
  leading <- c(
    "method", "details", "fitting_sample", "cross_fit", "cross_fit_folds",
    "conditioning"
  )
  rest <- provenance[setdiff(names(provenance), leading)]
  c(list(
      method = method,
      details = details,
      fitting_sample = as.character(fitting_sample)
    ),
    if (is.null(cross_fit)) NULL else list(cross_fit = cross_fit),
    list(
      cross_fit_folds = as.character(cross_fit_folds),
      conditioning = conditioning
    ),
    rest)
}

.transport_fixed_method <- function(provenance, method, details, extra) {
  if (!is.list(provenance)) {
    .input_error("`provenance` must be a list.", arg = "provenance")
  }
  declared <- provenance$method
  if (!is.null(declared) && !identical(declared, method)) {
    .input_error(sprintf(paste0(
      "This constructor builds a `%s` transport, so `provenance$method` ",
      "cannot be declared as %s. Use `location_transport()` or ",
      "`external_transport()` to declare a different method."
    ), method, .msg_value(declared)), arg = "provenance$method",
      received = .msg_value(declared), expected = sprintf("\"%s\"", method))
  }
  provenance$method <- method
  if (is.null(provenance$details)) provenance$details <- details
  for (name in names(extra)) {
    if (is.null(provenance[[name]])) provenance[[name]] <- extra[[name]]
  }
  provenance
}

.transport_coords <- function(coords, arg) {
  if (is.data.frame(coords)) coords <- as.matrix(coords)
  if (is.numeric(coords) && is.null(dim(coords))) {
    coords <- cbind(coords)
  }
  if (!is.matrix(coords) || !is.numeric(coords) || !all(is.finite(coords)) ||
      nrow(coords) < 1L || ncol(coords) < 1L) {
    .check_failed(arg,
      "a numeric matrix of node centres with no missing or infinite entries",
      coords)
  }
  unname(coords)
}

.transport_with_coords <- function(index, coords, arg) {
  labels <- paste0("coord", seq_len(ncol(coords)))
  clash <- intersect(labels, names(index))
  if (length(clash)) {
    .input_error(sprintf(paste0(
      "`%s` already carries %s, which is where the node centres are recorded. ",
      "Rename the column, or drop it and let the constructor write the ",
      "coordinates it was given."
    ), arg, paste0("`", clash, "`", collapse = ", ")), arg = arg)
  }
  for (position in seq_len(ncol(coords))) {
    index[[labels[[position]]]] <- as.numeric(coords[, position])
  }
  index
}

# Nearest group centre, blocked so that the n_native x m distance table is
# never materialized: a whole-brain native frame against a fine group node set
# would be gigabytes, and only the argmin and its distance survive the block.
.transport_nearest <- function(native_coords, group_coords, block = 4096L) {
  n <- nrow(native_coords)
  group_square <- rowSums(group_coords^2)
  group <- integer(n)
  distance <- numeric(n)
  start <- 1L
  while (start <= n) {
    last <- min(n, start + block - 1L)
    rows <- start:last
    chunk <- native_coords[rows, , drop = FALSE]
    squared <- outer(rowSums(chunk^2), group_square, "+") -
      2 * tcrossprod(chunk, group_coords)
    squared[squared < 0] <- 0
    nearest <- max.col(-squared, ties.method = "first")
    group[rows] <- nearest
    distance[rows] <- sqrt(squared[cbind(seq_along(rows), nearest)])
    start <- last + 1L
  }
  list(group = group, distance = distance)
}

# The two sink numbers population diagnostics report: unmapped native
# *territory*, which is a property of the operator alone and computable before
# any data is read, and the raw sink row mass it comes from.
.transport_sink_territory <- function(transport) {
  operator <- transport$matrix
  sink <- as.numeric(operator[, ncol(operator), drop = TRUE])
  mass <- transport$row_mass
  list(
    mass = sum(sink),
    share = sum(mass * sink) / sum(mass),
    rows = sum(sink > 0),
    all_sink = sum(sink >= 1)
  )
}
