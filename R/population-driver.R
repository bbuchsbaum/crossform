# Population execution -------------------------------------------------------
#
# Layer 4 (execution, and the records execution produces). `plan_population()`
# names the group estimand; this file is the only thing that computes a number
# from one. It is the population counterpart of `execution-driver.R`, and it
# calls sideways into that file's runner rather than up into
# `evaluate_geometry()`: the public entry points sit in layer 5, above the
# records an executor emits, and `crossnobis-driver.R` takes the same route for
# the same reason.
#
# `design/population-form-contract.md` (`population-form-v1`) is normative:
#
#   * the evaluation order is query -> transport -> OLS, and it is *forced*
#     rather than chosen: with heterogeneous native frames there is no common
#     node axis, so fit-then-transport is not even definable (section 3,
#     "scope of the tensor argument"). Query and transport act on different
#     tensor factors, so the answer does not depend on which of the two runs
#     first --- but which one *can* run first does (section 3);
#   * the transported total, sink included, is the native total exactly, for
#     signed ledgers. That is a certificate asserted per subject at fit time
#     with a tolerance scaled by the ledger's L1 norm, not by its signed total
#     (section 2);
#   * the budget normalization is a closed set of three that differ in argmax,
#     so it is applied here exactly as section 4.1 writes it and never as a
#     display option;
#   * transported components carry the names of section 8.1
#     (`transported_total`, `native_coherent_ledger`,
#     `native_configuration_ledger`, `sink_budget`), and the bare names
#     `coherent` and `configuration` are forbidden on a transported result.
#
# Why query-first is not an optimization. Complete geometry at a native node is
# `p = q(q+1)/2` packed coordinates; a bank of `K` queries is `K` numbers. The
# naive route materializes `n_i x p` per subject and transports `(m+1) x p`;
# this one holds `n_i x K` and transports `(m+1) x K`, and by section 3 the two
# agree exactly. For a four-condition study read through four contrasts that is
# ten coordinates against four, and the ratio grows like `q^2 / K`.
#
# Measured against the naive route written out in
# `tests/testthat/test-population-executor.R` (`pex_dense_oracle()`), on
# `N = 8` participants, `m = 500` group nodes, `K = 4` queries, native frames
# of 560-900 nodes, as the growth in R's vector high-water mark over the run:
#
#                      query-first   naive
#   q  =  4 (p = 10):    24.0 MiB   24.4 MiB
#   q  =  8 (p = 36):    27.4 MiB   39.1 MiB
#   q  = 12 (p = 78):    30.9 MiB   49.9 MiB
#
# The gap is the packed width, and it is bounded below by the fact that both
# routes run the same kernel: at four conditions the saving is 1.5 % and at
# twelve it is 38 %. The retained arrays are the cleaner statement of the same
# thing --- the transported group stack is `N(m+1)K` doubles here against
# `N(m+1)p` there (125 KiB vs 2.4 MiB at `q = 12`), and the largest native
# block held at once is `n_i K` against `2 n_i p`, the factor two being the
# separate `total` and `coherent` stores a materialized geometry keeps.

# The three admitted components and the transported ledger name each one takes.
# `population-form-v1` section 8.1 is the whole reason this table exists: the
# transported total *is* the group node's own total under a feature-additive
# metric, while the transported coherent part is a ledger of native-node
# coherence carried to a group location and is a different object from any
# group-node common mode. Naming them apart is what stops a reader inferring
# the second from the first.
.population_ledger_names <- c(
  total = "transported_total",
  coherent = "native_coherent_ledger",
  configuration = "native_configuration_ledger"
)

# Section 11: the budget-preservation tolerance is `1e-12 * sum_x |c_x|`. The
# L1 norm and not the total, because the total is a signed sum that may sit
# arbitrarily close to zero while its summands are large, and a
# relative-to-total tolerance would be unbounded exactly where the estimate is
# weakest.
.population_budget_tolerance <- 1e-12

# The query bank ---------------------------------------------------------------
#
# One bank parser for the whole package. `.sampling_query_bank()` is D8's, and
# reusing it is what makes `estimate_population(plan, queries)` and
# `sampling_covariance(x, queries = )` agree on what a bank *is* --- a numeric
# contrast vector, a K-by-q matrix of contrast rows, or a list of them, aligned
# against the relation's declared effect order and named once. A second parser
# here would be a second dialect.
.population_query_bank <- function(plan, queries) {
  if (missing(queries)) {
    .input_error(paste0(
      "`queries` is required: pass the bank of contrasts the group estimate ",
      "is about, for example `queries = rbind(faces = c(1, -1))`. A ",
      "population fit is a fit *of* a query --- there is no default readout ",
      "of a transported geometry, and guessing one would name an estimand ",
      "the plan never sealed."
    ),
      arg = "queries", received = "no argument",
      expected = "a contrast vector, a K-by-q contrast matrix, or a list")
  }
  effects <- plan$subjects[[1L]]$task$left_relation$effect_space$coordinates
  bank <- .sampling_query_bank(queries, effects)
  width <- length(effects) * (length(effects) + 1L) / 2L
  packed <- matrix(vapply(seq_len(nrow(bank$matrix)), function(k) {
    .svec_symmetric(tcrossprod(bank$matrix[k, ]))
  }, numeric(width)), nrow = width,
    dimnames = list(NULL, bank$labels))
  # The packed operator is exactly what `contrast_energy()` lowers a contrast
  # to (`bilinear_query(tcrossprod(weights))`), so a one-query bank and
  # `contrast_energy()` read the same physical query and not merely the same
  # arithmetic.
  list(matrix = bank$matrix, labels = bank$labels, packed = packed,
    effects = effects)
}

# Normalization ----------------------------------------------------------------
#
# Section 4.1's group functionals, written out rather than templated: `none` is
# the mean subject ledger and `unit_budget` is the *unweighted* mean of
# *unit-normalized* ledgers, and no single weighted-mean template covers both
# (the contract's own section 4.1 records an earlier draft getting this wrong).
# Both are per-subject scalar rescalings of the response, so both commute with
# transport (section 4.4, measured at 1.1e-16) --- which is why they are
# applied here, to the transported ledger, and give the same answer as
# rescaling the native one.
#
# The scalar is per subject *and per query*: `T_i` is the native total of one
# queried ledger, and a bank of K queries has K of them.
.population_scale <- function(plan, native, budget_floor) {
  totals <- colSums(native)
  if (identical(plan$normalization, "none")) {
    return(list(totals = totals, scale = rep(1, length(totals)),
      admitted = rep(TRUE, length(totals))))
  }
  # Section 4.3: `unit_budget`'s divisor is a signed estimate, so it can sit
  # near zero or be negative, and a subject whose total is not bounded away
  # from zero is marked NA rather than emitting a divergent share. The
  # criterion is scale-free --- `|T_i| >= floor * sum_x |c_ix|` --- because a
  # signed total may be a near-cancellation of large summands, and that is
  # exactly the case an absolute floor cannot see.
  norm <- colSums(abs(native))
  admitted <- abs(totals) > budget_floor * norm & totals != 0
  scale <- ifelse(admitted, 1 / totals, NA_real_)
  list(totals = totals, scale = scale, admitted = admitted)
}

# Section 4.1 admits `unit_budget` because it conserves: `sum_j g_j = 1`
# exactly. That is a statement about a *budget*, and density conserves nothing
# at all (section 1.3, measured: the density columns of the contract's own
# fixture satisfy no conservation law). Dividing a density profile by a native
# budget total therefore produces a number whose documented property is false,
# which is the one thing a receipt may not do.
.population_admit_normalization <- function(plan) {
  if (!identical(plan$normalization, "unit_budget") ||
      identical(plan$semantics, "budget")) {
    return(invisible(NULL))
  }
  .capability_refusal(paste0(
    "`normalization = \"unit_budget\"` divides each participant's ledger by ",
    "its own native budget total, and it is admitted because the result ",
    "conserves: every participant contributes budget one, so the group ",
    "ledger sums to one exactly. Density semantics transports budget per unit ",
    "transported territory, which satisfies no conservation law, so the ",
    "quotient conserves nothing and the share it looks like is a share of ",
    "nothing. crossform will not print a number under a name whose stated ",
    "property it does not have."
  ),
    capability = "budget_normalization_semantics",
    namespace = "population_execution",
    reasons = c(
      "unit_budget_requires_budget_semantics",
      paste0("declared_semantics:", plan$semantics)
    ),
    remedies = c(
      paste0(
        "Plan the population with `normalization = \"none\"`, which is the ",
        "mean subject density profile and makes no conservation claim."
      ),
      paste0(
        "Declare `semantics = \"budget\"` on every transport if the ",
        "attribution *share* is the estimand you want."
      )
    ))
}

# Budget preservation ----------------------------------------------------------
#
# Section 2's certificate, asserted per subject at fit time on the raw
# transported ledger --- before any normalization, because a per-subject scalar
# rescales both sides and would make the check vacuous.
#
# The argument is `P 1 = 1` and never `c >= 0`, so this holds for the signed
# crossvalidated ledgers a conservative geometry actually produces, and it
# holds for every component: the coherent and configuration ledgers are
# ledgers too. It is checked rather than assumed because deleting the sink
# column from a well-formed `P` loses 48.1 % of the native total on the
# contract's own fixture and leaves no trace anywhere else.
.population_budget_certificate <- function(label, native, carried) {
  native_total <- colSums(native)
  carried_total <- colSums(carried)
  norm <- pmax(colSums(abs(native)), .Machine$double.eps)
  deviation <- abs(carried_total - native_total) / norm
  worst <- which.max(deviation)
  if (deviation[[worst]] > .population_budget_tolerance) {
    .contract_error(sprintf(paste0(
      "Budget preservation failed for participant `%s` on query `%s`: the ",
      "transported total including the sink is %s against a native total of ",
      "%s, a relative deviation of %s. A row-stochastic transport preserves ",
      "the signed total exactly, so a failure here means the operator is ",
      "ill-formed rather than imprecise, and crossform refuses rather than ",
      "renormalizing it."
    ), label, colnames(native)[[worst]], format(carried_total[[worst]]),
      format(native_total[[worst]]), format(deviation[[worst]])))
  }
  max(deviation)
}

# The group index --------------------------------------------------------------
#
# The plan's `m` group nodes plus the sink, which is a row of the result and
# not a footnote: section 3.3 fits it and section 1.1 requires it materialized
# even when it is identically zero, because unmapped territory has to be a
# number a reader can see. `units` is what stops a density result reading as if
# its sink were a density: the sink is reported in budget units under either
# semantics, since a density of unmapped territory has no referent
# (section 1.3).
.population_result_index <- function(plan) {
  group <- plan$group_index
  # A group index is an open table --- a caller may carry any column on it ---
  # so the two the result adds have to be checked rather than assumed free.
  # Overwriting a caller's own `sink` or `units` column would replace declared
  # metadata with the result's bookkeeping and leave nothing saying so;
  # `location_transport()` refuses a `coord1` clash for the same reason.
  clash <- intersect(c("sink", "units"), names(group))
  if (length(clash)) {
    .capability_refusal(sprintf(paste0(
      "The group index already carries %s, which is where a population ",
      "result records whether a row is the sink and which units its values ",
      "are in. crossform will not overwrite a declared column with its own ",
      "bookkeeping."
    ), paste0("`", clash, "`", collapse = " and ")),
      capability = "reserved_group_index_columns",
      namespace = "population_execution",
      reasons = paste0("group_index_column_reserved:", clash),
      remedies = paste0(
        "Rename ", paste0("`", clash, "`", collapse = " and "),
        " in the `group_index` every transport was built against, then ",
        "re-plan the population."
      ))
  }
  sink <- group[1L, , drop = FALSE]
  for (column in names(sink)) {
    sink[[column]] <- group[[column]][NA_integer_]
  }
  sink$node <- .transport_sink_label
  index <- rbind(group, sink)
  index$sink <- c(rep(FALSE, nrow(group)), TRUE)
  index$units <- c(rep(plan$semantics, nrow(group)), "budget")
  rownames(index) <- NULL
  index
}

# Uncertainty ------------------------------------------------------------------
#
# What the transported covariance would need, and why it is refused.
#
# The transported value at group node `j` is `b_j = sum_x P_xj z_x`, so
#
#   Cov(b_j) = sum_x sum_x' P_xj P_x'j Cov(z_x, z_x')
#
# and every term with `x != x'` is a *cross-node* sampling covariance. D8
# refuses exactly that object under capability `cross_node_sampling_covariance`
# (reasons `spatial_covariance_model_unavailable` and
# `conservation_gives_no_uncertainty`), and `population-form-v1` section 4.5
# records that the refusal stands. A transport with more than one native row
# per group column --- which is every transport that does anything --- needs
# those terms, so the transported covariance cannot be assembled from what
# exists. What *can* be carried honestly is the per-native-node marginal block
# D8 does supply, untransported and labelled as such.
#
# The one case where the cross terms vanish is a group column fed by exactly
# one native row, where `Cov(b_j) = P_xj^2 Cov(z_x)`. It is not implemented
# here: a real transport is a mixture of that case and the general one, so
# shipping it means deciding what a result looks like when some nodes carry a
# covariance and others carry an absence, and that adjudication belongs with
# the uncertainty layer rather than to a passthrough. The condition is recorded
# in the refusal so the successor ticket does not have to rediscover it.
.population_transport_covariance_refusal <- function() {
  list(
    capability = "transported_sampling_covariance",
    namespace = "population_execution",
    reasons = c(
      "cross_node_sampling_covariance_unavailable",
      "transported_value_mixes_native_nodes"
    ),
    remedies = c(
      paste0(
        "Read `$uncertainty$native` as what it is: one K-by-K marginal ",
        "covariance per native node per participant, before transport."
      ),
      paste0(
        "The refusal lifts for a group column fed by exactly one native row, ",
        "where the cross-node terms carry weight zero; a general transport ",
        "has no such column."
      )
    ),
    message = paste0(
      "The sampling covariance of a transported value is refused. A group ",
      "node's value is a weighted sum over many native nodes, so its ",
      "covariance needs the covariance *between* native nodes, and crossform ",
      "has no route that produces one: the analytic sampling law is local to ",
      "one measurement, and a conservative frame supplies no substitute, ",
      "because conservation is a law about estimates and not about their ",
      "uncertainty. The per-native-node blocks are carried untransported."
    )
  )
}

.population_uncertainty <- function(plan, bank, uncertainty) {
  if (is.null(uncertainty)) return(NULL)
  labels <- names(plan$subjects)
  if (!is.list(uncertainty) || inherits(uncertainty, "effect_sampling_covariance") ||
      inherits(uncertainty, "effect_sampling_covariance_batch")) {
    .input_error(paste0(
      "`uncertainty` must be a named list of `rdm_sampling_covariance()` ",
      "results, one per participant. A single covariance belongs to one ",
      "participant's own frame, so pairing it with the population by ",
      "position would attach one person's error bars to another's estimate."
    ), arg = "uncertainty", received = .msg_value(uncertainty),
      expected = "a named list of sampling covariances, one per participant")
  }
  given <- names(uncertainty)
  if (!.is_strings(given, unique = TRUE) || !setequal(given, labels)) {
    .capability_refusal(paste0(
      "Every participant needs exactly one sampling covariance and every ",
      "covariance needs a participant. The names are the binding, exactly as ",
      "they are for the transports."
    ),
      capability = "aligned_subject_uncertainty",
      namespace = "population_execution",
      reasons = c(
        paste0("subject_without_uncertainty:",
          setdiff(labels, as.character(given))),
        paste0("uncertainty_without_subject:",
          setdiff(as.character(given), labels))
      ),
      remedies = paste0(
        "Name `uncertainty` with the participant identifiers of the plan: ",
        .msg_names(labels), "."
      ))
  }
  # D8's own route, per subject. It refuses an uncentred contrast --- the
  # crossvalidated distance basis carries no information about one --- and that
  # refusal propagates rather than being caught, because a bank the covariance
  # cannot express is a bank the caller has to know about.
  native <- lapply(uncertainty[labels], function(value) {
    sampling_covariance(value, queries = bank$matrix)
  })
  names(native) <- labels
  list(
    native = native,
    basis = "query_bank",
    scope = "native_node_marginal",
    transported = .population_transport_covariance_refusal()
  )
}

# Identity ---------------------------------------------------------------------
#
# The result's own scientific identity: its parent estimand, the bank it was
# read through, the component it names, and the normalization criterion that
# decided which participants contributed a share. `budget_floor` is in it
# because it decides which subjects are marked NA, and a different set of
# contributing participants is a different population.
.population_result_id <- function(plan, bank, component, budget_floor) {
  .sha256_signature(list(
    schema_version = 1L,
    contract_version = "population-form-v1",
    role = "population_result",
    parent = plan$scientific_plan_id,
    component = component,
    ledger = unname(.population_ledger_names[[component]]),
    queries = unname(bank$matrix),
    query_labels = bank$labels,
    normalization = plan$normalization,
    budget_floor = if (identical(plan$normalization, "unit_budget")) {
      as.numeric(budget_floor)
    } else {
      NULL
    }
  ), "population-sha256:")
}

# The verb ---------------------------------------------------------------------

#' Estimate a planned population form through a bank of queries
#'
#' `estimate_population()` is the authorized group-level execution verb: it
#' takes the estimand [plan_population()] sealed, reads each participant's
#' geometry through one bank of `K` contrasts, carries the queried ledgers onto
#' the shared group nodes with that participant's declared transport, and fits
#' the group model once at every group node and query.
#'
#' @section The order, and why it is not a choice:
#' The evaluation order is **query, then transport, then fit**. Query and
#' transport act on different axes of the same array --- the query mixes
#' experimental coordinates, the transport mixes spatial nodes --- so they
#' commute exactly, and under a subject-constant weight operator the fit
#' commutes with both. What is *not* free is which may run first: participants
#' have different native frames, so there is no common node axis to fit at
#' before transporting, and fit-then-transport is not definable at all. The
#' order is recorded on the plan rather than selected here.
#'
#' Reading the query first is what makes the route cheap. Complete geometry at
#' a native node is \eqn{q(q+1)/2}{q(q+1)/2} packed coordinates; a bank of
#' `K` queries is `K` numbers, and complete geometry is never allocated.
#'
#' @section What the result is a ledger of:
#' `component` selects which part of each participant's conservative geometry
#' is carried, and the transported object takes its own name
#' (`population-form-v1` section 8.1) because it is not the group-node quantity
#' the bare name would suggest.
#'
#' | `component` | `$ledger` | what it is |
#' |---|---|---|
#' | `"total"` | `transported_total` | a genuine group-node total for that participant |
#' | `"coherent"` | `native_coherent_ledger` | native-node coherent evidence carried to a group node |
#' | `"configuration"` | `native_configuration_ledger` | native-node configuration evidence carried to a group node |
#'
#' `native_coherent_ledger + native_configuration_ledger = transported_total`
#' holds exactly; what fails is reading the coherent ledger as a group-node
#' common mode. Transport carries *nodes*, so there is no group frame and no
#' group-node weight vector for such a mode to be defined against.
#'
#' The sink is a row of the result, labelled `<sink>`, not a discarded
#' remainder: it is fitted like any other column because covariate-linked sink
#' mass is how a differentially failing transport becomes visible. It is never
#' a value at a location, and it is always reported in budget units.
#'
#' @section Budget preservation:
#' Under `"budget"` semantics each participant's transported total, sink
#' included, equals its native total exactly. The certificate is asserted per
#' participant and per query at fit time, with a tolerance of `1e-12` times the
#' ledger's \eqn{L^1}{L1} norm --- not its total, which is a signed sum that
#' may sit near zero while its summands are large. A participant that fails it
#' has an ill-formed transport, and the run refuses rather than renormalizing.
#'
#' @section Normalization:
#' `"none"` fits the ledgers as they are, so participants with more evidence
#' weigh more. `"unit_budget"` divides each participant's ledger by its own
#' native total first, so each contributes budget one; the divisor is a
#' *signed* estimate, so a participant whose total is not bounded away from
#' zero by `budget_floor` is marked `NA` for that query rather than emitting a
#' divergent share. The default `budget_floor = 0` admits every nonzero total,
#' which is the weakest honest criterion and is recorded as such in the
#' receipt: the constant is an open maintainer decision, and it needs a
#' per-participant standard error that does not exist yet. Both normalizations
#' are computed against a total read from the same data as the ledger, which
#' the receipt declares (`budget_estimate = "same_data_ratio"`).
#'
#' @section Refusals:
#' Each is an `effect_capability_refusal` (see [catch_refusal()]).
#'
#' * `budget_normalization_semantics` --- `"unit_budget"` on a `"density"`
#'   plan. Density conserves nothing, so the share would be a share of
#'   nothing.
#' * `aligned_subject_uncertainty` --- `uncertainty` names that do not match
#'   the participants.
#' * `distance_basis_query_bank` --- an uncentred contrast in the bank when
#'   `uncertainty` is supplied. Point estimates admit uncentred contrasts;
#'   their sampling law does not.
#'
#' @param plan An `effect_population_plan` from [plan_population()].
#' @param queries The bank of contrasts to read: a numeric contrast vector, a
#'   `K`-by-`q` matrix with one contrast per row, or a list of contrast
#'   vectors, over the participants' shared experimental effects. Named rows
#'   name the queries; unnamed ones are numbered. This is the same bank
#'   [sampling_covariance()] takes.
#' @param component Which component of each participant's conservative
#'   geometry to carry: `"total"` (the default), `"coherent"`, or
#'   `"configuration"`. See the table above for the name the transported
#'   object takes.
#' @param uncertainty Optional named list of [rdm_sampling_covariance()]
#'   results, one per participant, lowered onto the query bank and carried
#'   **untransported**. The covariance of a transported value is refused, and
#'   the refusal record travels with the result.
#' @param budget_floor Nonnegative relative floor for the `"unit_budget"`
#'   divisor: a participant contributes a share for query `k` only when
#'   \eqn{|T_i| > \mathrm{floor} \cdot \sum_x |c_{ix}|}{|T_i| > floor * sum_x
#'   |c_ix|}. Ignored under `"none"`.
#' @return An `effect_population_result`: a sealed record carrying
#'   `$coefficients` (a `node`-by-`query`-by-`term` array), `$values`,
#'   `$fitted` and `$residuals` (`node`-by-`query`-by-`subject`), the
#'   `$index` of group nodes plus the sink, `$queries`, `$ledger`,
#'   `$uncertainty`, and a `$receipt` recording every participant's read, its
#'   transport signature, the budget certificate, and the normalization.
#'   `as.data.frame()` returns the coefficient table in long form.
#' @references `design/population-form-contract.md` (`population-form-v1`),
#'   sections 2, 3, 4 and 8.
#' @family population transports
#' @seealso [plan_population()] for the estimand this executes,
#'   [transport_values()] for the carrying step, and [contrast_energy()] for
#'   the single-participant reading of the same query.
#' @examples
#' # Four participants on different native frames, two shared group nodes.
#' effects <- effect_space(c("face", "house"), basis_id = "population:v1")
#' subject <- function(id, n, gain) {
#'   domain <- abstract_domain(n, coordinates = cbind(x = seq_len(n) - 1),
#'     feature_ids = paste0("f", seq_len(n)), id = id)
#'   values <- function(divisor) matrix(
#'     gain * seq_len(2 * n) / (n * divisor), 2, n,
#'     dimnames = list(c("face", "house"), NULL)
#'   )
#'   rel <- relation(list(run1 = values(1), run2 = values(2)),
#'     effects = effects, domain = domain)
#'   plan_geometry(rel, compile_frame(voxelwise(), domain),
#'     cross_partitions(rel))
#' }
#' carrier <- function(n) anatomical_transport(
#'   native_coords = cbind(seq_len(n) - 1), group_coords = cbind(c(0, 3)),
#'   semantics = "budget"
#' )
#' sizes <- c(s01 = 5L, s02 = 6L, s03 = 7L, s04 = 8L)
#' gains <- c(s01 = 1, s02 = 1.4, s03 = 0.7, s04 = 1.1)
#' subjects <- stats::setNames(lapply(names(sizes), function(id)
#'   subject(id, sizes[[id]], gains[[id]])), names(sizes))
#' plan <- plan_population(subjects, lapply(sizes, carrier))
#'
#' # One bank of two contrasts, read at every group node in one pass.
#' fit <- estimate_population(plan, rbind(
#'   `face-house` = c(1, -1), `face+house` = c(1, 1)
#' ))
#' fit
#'
#' # The group mean at each node and query, with the sink as its own row.
#' fit$coefficients[, , "(Intercept)"]
#'
#' # Budget preservation: each participant's transported total, sink
#' # included, is its native total.
#' fit$receipt$budget$max_relative_deviation < 1e-12
#'
#' # The transported object is named for what it is.
#' fit$ledger
#' @export
estimate_population <- function(plan, queries,
                                component = c("total", "coherent",
                                  "configuration"),
                                uncertainty = NULL, budget_floor = 0) {
  .validate_population_plan(plan, deep = FALSE)
  component <- match.arg(component)
  budget_floor <- .check_number(budget_floor, "budget_floor",
    nonnegative = TRUE,
    what = "one nonnegative relative floor on the unit-budget divisor")
  .population_admit_normalization(plan)
  bank <- .population_query_bank(plan, queries)

  labels <- names(plan$subjects)
  n_subject <- length(labels)
  n_query <- length(bank$labels)
  index <- .population_result_index(plan)
  nodes <- nrow(index)
  node_labels <- as.character(index$node)

  # The stack is `N` by `(m + 1) * K`, and it is the only object here whose
  # size grows with the group node count. Nothing native survives the loop:
  # each participant's `n_i` by `K` view is contracted to `(m + 1)` by `K` and
  # dropped before the next participant is opened, so peak memory is set by
  # the group side rather than by the largest native frame.
  stack <- matrix(NA_real_, n_subject, nodes * n_query,
    dimnames = list(labels, NULL))
  receipts <- vector("list", n_subject)
  names(receipts) <- labels
  native_totals <- matrix(NA_real_, n_subject, n_query,
    dimnames = list(labels, bank$labels))
  sink_budget <- native_totals
  admitted <- matrix(TRUE, n_subject, n_query,
    dimnames = list(labels, bank$labels))
  deviations <- numeric(n_subject)

  for (position in seq_len(n_subject)) {
    label <- labels[[position]]
    carrier <- plan$transport[[label]]
    # Sideways into the executor's own runner, exactly as `crossnobis()` does:
    # `evaluate_geometry()` is the layer-5 entry point for the same call, and
    # a driver must not reach up into the views it feeds.
    view <- .run_geometry_compiler(
      plan$subjects[[label]], query = bank$packed, component = component
    )
    native <- view$values
    colnames(native) <- bank$labels
    receipts[[label]] <- view$receipt

    carried <- transport_values(carrier, native)
    if (identical(plan$semantics, "budget")) {
      deviations[[position]] <-
        .population_budget_certificate(label, native, carried)
    }
    native_totals[position, ] <- colSums(native)
    sink_budget[position, ] <- as.numeric(
      Matrix::crossprod(carrier$matrix[, ncol(carrier$matrix), drop = FALSE],
        native)
    )

    scaling <- .population_scale(plan, native, budget_floor)
    admitted[position, ] <- scaling$admitted
    carried <- carried * rep(scaling$scale, each = nodes)
    stack[position, ] <- as.numeric(carried)
    rm(view, native, carried)
  }

  fit <- .population_ols(plan$model, stack)
  shape <- function(values, third, third_name) {
    array(t(values), dim = c(nodes, n_query, length(third)),
      dimnames = stats::setNames(
        list(node_labels, bank$labels, third),
        c("node", "query", third_name)
      ))
  }

  value <- structure(list(
    coefficients = shape(fit$coefficients, plan$model$columns, "term"),
    values = shape(stack, labels, "subject"),
    fitted = shape(fit$fitted, labels, "subject"),
    residuals = shape(fit$residuals, labels, "subject"),
    residual_df = as.integer(n_subject - plan$model$rank),
    index = index,
    queries = bank$matrix,
    component = component,
    ledger = unname(.population_ledger_names[[component]]),
    semantics = plan$semantics,
    normalization = plan$normalization,
    uncertainty = .population_uncertainty(plan, bank, uncertainty),
    receipt = .population_receipt(
      plan, bank, component, budget_floor, receipts, native_totals,
      sink_budget, admitted, deviations, fit$unresolved
    ),
    scientific_plan_id = .population_result_id(
      plan, bank, component, budget_floor
    )
  ), class = "effect_population_result")
  .validate_population_result(value)
  value
}

# The group fit ----------------------------------------------------------------
#
# One backsolve for the whole result. The plan stored the design's pivoted QR
# because the design varies along neither the node axis nor the coordinate
# axis, so the same factorization serves every group node and every query; the
# per-(node, query) loop section 3 describes is the column loop inside LAPACK
# rather than an R one, and `Theta` comes back for all `(m + 1) * K` responses
# at once.
#
# A column carrying a non-finite response is solved around rather than through.
# Two things produce one, and both are measurements: density semantics returns
# `NA` at a group node reached by no native mass (never `0`, because zero is a
# measurement and absence is not), and `unit_budget` marks a participant whose
# signed native total is not bounded away from zero. Dropping the participant
# instead would silently change which population the estimate is about.
.population_ols <- function(model, stack) {
  finite <- apply(is.finite(stack), 2L, all)
  coefficients <- matrix(NA_real_, ncol(model$matrix), ncol(stack),
    dimnames = list(model$columns, NULL))
  fitted <- matrix(NA_real_, nrow(stack), ncol(stack))
  residuals <- fitted
  if (any(finite)) {
    responses <- stack[, finite, drop = FALSE]
    coefficients[, finite] <- qr.coef(model$qr, responses)
    fitted[, finite] <- qr.fitted(model$qr, responses)
    residuals[, finite] <- qr.resid(model$qr, responses)
  }
  list(coefficients = coefficients, fitted = fitted, residuals = residuals,
    unresolved = as.integer(sum(!finite)))
}

# The receipt ------------------------------------------------------------------
#
# What a reader would have to know to say what these numbers mean, and what
# they would have to check to disbelieve them: who was read and under which
# plan identity, how each participant was carried and under what provenance,
# whether the budget closed, which normalization made incommensurable native
# totals commensurable and on what criterion, and --- section 8.1's required
# print line --- which native frame family the ledger belongs to.
.population_receipt <- function(plan, bank, component, budget_floor,
                                receipts, native_totals, sink_budget,
                                admitted, deviations, unresolved) {
  # The per-query numbers stay as `subject`-by-`query` matrices rather than
  # becoming list columns of the audit table: a list column is not portable,
  # `location_transport()` refuses one in an index for exactly that reason, and
  # a matrix reads the same way whether the bank holds one query or twenty.
  subjects <- plan$subject_index
  subjects$budget_deviation <- if (identical(plan$semantics, "budget")) {
    deviations
  } else {
    rep(NA_real_, nrow(subjects))
  }
  subjects$contributes <- rowSums(admitted) == ncol(admitted)
  list(
    contract_version = "population-form-v1",
    population_plan_id = plan$scientific_plan_id,
    population_signature = plan$signature,
    subjects = subjects,
    native_total = native_totals,
    sink_budget = sink_budget,
    subject_receipts = receipts,
    component = component,
    ledger = unname(.population_ledger_names[[component]]),
    frame = .population_frame_identity(plan),
    semantics = plan$semantics,
    evaluation_order = plan$fit$evaluation_order,
    fit = plan$fit,
    queries = bank$labels,
    normalization = list(
      mode = plan$normalization,
      # Section 4.3(a): the divisor is read from the same data as the ledger,
      # so it is declared descriptively rather than presented as an
      # independent estimate.
      budget_estimate = if (identical(plan$normalization, "unit_budget")) {
        "same_data_ratio"
      } else {
        NA_character_
      },
      budget_floor = if (identical(plan$normalization, "unit_budget")) {
        as.numeric(budget_floor)
      } else {
        NA_real_
      },
      floor_criterion = if (identical(plan$normalization, "unit_budget")) {
        "abs(native_total) > budget_floor * sum(abs(native_ledger))"
      } else {
        NA_character_
      },
      # One logical per participant and query: `TRUE` where that
      # participant's ledger was admitted a share. Reported as the admissions
      # rather than the exclusions so the empty case reads as "everyone
      # contributed" instead of as an absent field.
      admitted = admitted
    ),
    budget = list(
      asserted = identical(plan$semantics, "budget"),
      tolerance = .population_budget_tolerance,
      scale = "relative_to_ledger_l1_norm",
      max_relative_deviation = if (identical(plan$semantics, "budget")) {
        max(deviations)
      } else {
        NA_real_
      }
    ),
    unresolved_columns = unresolved
  )
}

# Section 8.1's required print line needs the native frame family the ledger
# belongs to. The transports carry per-row frame metadata when the frame was
# built by `frame_family()`; when they do not, the plan's own audit table still
# names the declared normalization, and saying "undeclared" is the honest
# answer rather than inventing one.
.population_frame_identity <- function(plan) {
  families <- unlist(lapply(plan$transport, function(carrier) {
    if ("family" %in% names(carrier$native_index)) {
      unique(as.character(carrier$native_index$family))
    } else {
      NULL
    }
  }), use.names = FALSE)
  scales <- unlist(lapply(plan$transport, function(carrier) {
    if ("scale" %in% names(carrier$native_index)) {
      unique(as.numeric(carrier$native_index$scale))
    } else {
      NULL
    }
  }), use.names = FALSE)
  list(
    family = if (length(families)) sort(unique(families)) else NA_character_,
    scale = if (length(scales)) sort(unique(scales)) else NA_real_,
    normalization = sort(unique(plan$subject_index$declared_normalization)),
    row_metadata = all(plan$subject_index$row_metadata)
  )
}

# The record -------------------------------------------------------------------

.population_result_fields <- c(
  "coefficients", "values", "fitted", "residuals", "residual_df", "index",
  "queries", "component", "ledger", "semantics", "normalization",
  "uncertainty", "receipt", "scientific_plan_id"
)

.validate_population_result <- function(x) {
  if (!inherits(x, "effect_population_result")) {
    .input_error(sprintf(paste0(
      "Expected an `effect_population_result` from `estimate_population()`; ",
      "received %s."
    ), .msg_value(x)), arg = "x", received = .msg_value(x),
      expected = "an `effect_population_result` from `estimate_population()`")
  }
  if (!.sealed_fields(x, "effect_population_result",
      .population_result_fields) ||
      !is.array(x$coefficients) || length(dim(x$coefficients)) != 3L ||
      !is.array(x$values) || length(dim(x$values)) != 3L ||
      !is.data.frame(x$index) || !is.matrix(x$queries) ||
      !.is_string(x$component) || !.is_string(x$ledger) ||
      !.is_string(x$semantics) || !.is_string(x$normalization) ||
      !is.integer(x$residual_df) || !is.list(x$receipt) ||
      !.strong_sha256(sub("^population-", "", x$scientific_plan_id))) {
    .input_error("Population-result fields are missing or noncanonical.")
  }
  shape <- dim(x$coefficients)
  if (shape[[1L]] != nrow(x$index) || shape[[2L]] != nrow(x$queries) ||
      !identical(dim(x$values), dim(x$fitted)) ||
      !identical(dim(x$values), dim(x$residuals)) ||
      dim(x$values)[[1L]] != shape[[1L]] ||
      dim(x$values)[[2L]] != shape[[2L]] ||
      dim(x$values)[[3L]] != nrow(x$receipt$subjects)) {
    .contract_error(paste0(
      "Population-result arrays do not agree with their group node index, ",
      "query bank and participant list."
    ))
  }
  if (!identical(unname(.population_ledger_names[[x$component]]), x$ledger)) {
    .contract_error(paste0(
      "A transported component must carry the ledger name ",
      "`population-form-v1` section 8.1 assigns it; this result names a ",
      "different one."
    ))
  }
  if (sum(x$index$sink) != 1L ||
      !identical(as.character(x$index$node[x$index$sink]),
        .transport_sink_label)) {
    .contract_error(paste0(
      "A population result carries exactly one sink row, materialized even ",
      "when it is empty."
    ))
  }
  invisible(x)
}
