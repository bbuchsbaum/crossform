# `crossform` API tier ledger

Status: ticket A1 of WS-A (subtraction release)
Date: 2026-08-17
Companion to: [architecture.md](architecture.md) (source layering),
[crossform-package-design.md](crossform-package-design.md) §4 (the 0.1 public API list),
[ingestion-contract.md](ingestion-contract.md) §13 (compatibility promises)
Input inventory: [`.planning/2026-08-17-feedback-assessment.md`](../.planning/2026-08-17-feedback-assessment.md)

## Why this document exists

`NAMESPACE` exports 105 functions. Nothing in the tree says which of them a
first-time user is supposed to meet, which are a specialist's tools, and which
exist only so an extension package can hand data in. Without that statement the
package cannot be subtracted from: every export looks equally load-bearing, and
the only available evidence — "some test calls it" — is satisfied by all 105.

This ledger assigns every export a **tier** (who it is for) and a
**disposition** (what happens to it in the subtraction release). It is the
authority the rest of WS-A works against: `_pkgdown.yml` regrouping, the
`\keyword{internal}` demotions, and the README narrative split all read their
tier from here.

## Method

The user columns are mechanical, not remembered. For each export name `N`, a
script matched the call pattern `(^|[^A-Za-z0-9._])N\s*\(` against four corpora,
excluding `N`'s own defining file:

| Corpus | Files |
|---|---|
| vignette | `vignettes/*.Rmd` (9 files) |
| README | `README.md` |
| exemplar | `exemplars/**/*.R` (8 files) |
| tests | `tests/testthat/{test,helper}-*.R` |

Two derived corpora were also measured and are cited in justifications where
they change a call: **foreign Rd examples** (an export called in the
`\examples{}` block of a *different* topic — these run under `R CMD check`, so
demoting the callee breaks the caller), and **exported-consumer** (does any
exported function accept this value as an argument?). Scripts:
`inventory.R`, `titles.R`, `pkgd.R`, `testfiles.R`, `rdex.R` (session
scratchpad; regenerable from the method above).

Notation in the *Users* column: `V:` vignette stems, `README`, `E:` exemplar
scripts, `T:n` = number of test/helper files that call it. A `†` marks a
vignette/README credit that is a **prose mention only** — see the caveat below.

**Method caveat (added by the 2026-08-17 fresh-context review).** The call
pattern above is matched against the *whole* `.Rmd`/`.md` file, so a name
written in backticks in running prose — `` `coupling(plan, between, by)` `` in a
narrative table cell — counts as a "user" indistinguishably from a call inside
an executed chunk. Two vignettes, `novelty.Rmd` and
`correlation-distance-policy.Rmd`, contain **no code fences at all**; every
`V:novelty` and `V:corr-dist` credit in this ledger is therefore a prose
mention. Restricting the corpora to fenced R chunks (comments stripped) moves
eleven exports into the tests-only set, taking it from 35 to **46**. The
per-row Users columns are left in the original notation for continuity with the
tier tables, with `†` marking prose-only credits; the review section at the
foot of this document lists the eleven and what follows from them.

Likewise, the **foreign Rd examples** counts as originally measured included
`\seealso` `\link{}` cross-references, `\usage{}` defaults, `\description{}`
prose, and `#` comments inside example blocks. Re-measured with `tools::Rd2ex`
(R's own example extractor) and comments stripped, 29 of the 35 counts were
overstated; the corrected values are in the tests-only table. This matters
because the count is used as a *cost of demotion* — and a `\seealso` link to a
`\keyword{internal}` topic costs nothing, since the topic still exists.

## Tier definitions

| Tier | Who it is for | Doc obligation |
|---|---|---|
| **core** | Anyone doing a representational analysis. Meeting all of these is the price of entry. | Every one must appear in a vignette or the README, with a worked example. |
| **core-ingestion** | The typed-facts spine: study → observations → design → effects → relation plan. A user who brings raw data rather than a fitted relation meets all of it; a user who starts from betas meets none of it. | `from-observations.Rmd` is its single narrative. Separately titled so the core count is not inflated for the beta-first user. |
| **advanced** | A specialist reaching past the four view functions: explicit queries, frame algebra, metrics, crossnobis, error-channel readers, numerical evidence. | Reference page + runnable Rd examples. Vignette coverage desirable, not required. |
| **advanced (experimental: connect)** | Measurement forms, coupling closures, tomography. The connect narrative. Scale-gated, and the package refuses brain-scale claims here. | Must carry the experimental label at the top of its reference section and in `_pkgdown.yml`. |
| **adapters** | Labeled sub-tier of advanced: morphisms into the typed core from BIDS, fmridesign, fmrireg, neuroim2. | Each adapter's dependency is `Suggests:`; the core is not shaped by any of them. |
| **developer** | Extension and adapter *packages*, not end users: the **extension-*only*** protocol for handing data, a design, and an error channel into the core. It is not the whole surface an extension package meets — an adapter *also* meets the core-ingestion tier (the four in-tree adapters use thirteen exported crossform functions between them: fourteen before ticket A3 applied decision 1, twelve after it, and thirteen once ticket D3 gave `neuroim2_searchlights()` a multiscale route that delegates its refusals to `searchlights()` and stacks the result with `frame_family()` — of which only two are developer-tier). | Target: **≤ 5 sanctioned entry points remain exported.** Everything else that was developer-tier becomes internal. |

**On the developer tier's boundary** (maintainer decision 7, 2026-08-17). Read
literally — "extension and adapter *packages*, not end users" — the definition
also captures `compiler_conformance()` and `relation_plan_receipts()`, since the
admission test and receipt reader for an external design compiler are by
construction called by another package. Both are nonetheless **kept exported**
and tiered core-ingestion, on the ground that `from-observations.Rmd` genuinely
*executes* them (lines 254, 265, 350 — real chunks, not prose). The tier is
therefore defined by *what the surface is for* (an extension-only seam) rather
than by *who happens to call it*; an export a published vignette runs is not an
extension-only seam whatever its caller's package.

Precedent for demotion: 13 `man/*.Rd` topics already carry `\keyword{internal}`
(`compile_metric_schedule`, `effect_view`, `nonlinear_query`,
`coherent_functional`, `compile_lowering`, `factor_frame`, `metric_components`,
`execution_receipt`, `memory_plan`, `materialize_metric`, `pairing_marginals`,
`effect_geometry`, `crossform-package`). `pairing_marginals` is named in
`crossform-package-design.md` §4's public API list and is nonetheless internal —
that list is a 0.1 sketch, not a binding contract.

## Summary

### Counts per tier

| Tier | Exports | Keep | Removed from `NAMESPACE` |
|---|---:|---:|---:|
| core | 21 | 21 | 0 |
| core-ingestion | 20 | 20 | 0 |
| advanced | 33 | 32 | 1 |
| advanced (experimental: connect) | 16 | 13 | 3 |
| adapters | 8 | 6 | 2 |
| developer | 7 | 5 | 2 |
| **total** | **105** | **97** | **8** |

(Advanced is 33 and developer 7 because maintainer decision 2 re-tiers
`lm_extractor` from developer to advanced; adapters lose 2 to decision 1.)

### Counts per disposition

| Disposition | Count | Names |
|---|---:|---|
| keep | 97 | — |
| demote-internal | 7 | `inner_product`, `measurement_space`, `measurement_bridge`, `effect_covariance`, `residual_pair_statistics`, `bids_events`, `bids_confounds` |
| merge-into | 1 | `reverse_bridge` → `measurement_bridge` |
| deprecate | 0 | — |

**Post-subtraction: 97 exports.** The arithmetic: 105 exports − 7
`demote-internal` = 98, then − 1 for the merge = **97**. The merge subtracts a
full name rather than folding two names into one surviving export, because
`reverse_bridge` merges into `measurement_bridge` — and `measurement_bridge` is
itself one of the seven demotions. Both names leave `NAMESPACE`; the merged
behaviour survives as the internal `measurement_bridge(reverse = TRUE)`.
Equivalently, the "Removed from `NAMESPACE`" column above sums to 8 and
105 − 8 = 97.

Of those 97, the developer tier is exactly **5** (`file_matrix_source`,
`source_capabilities`, `effect_extractor`, `relation_fit`, `relation_block`).
WS-A's target was "≤ 6 sanctioned entry points"; maintainer decision 2 tightened
it to **≤ 5** by re-tiering `lm_extractor`, and the tier definition now states
the tighter figure. The set meets it exactly.

### Tier × disposition

| Tier | keep | demote-internal | merge-into |
|---|---:|---:|---:|
| core | 21 | 0 | 0 |
| core-ingestion | 20 | 0 | 0 |
| advanced | 32 | 1 | 0 |
| advanced (experimental: connect) | 13 | 2 | 1 |
| adapters | 6 | 2 | 0 |
| developer | 5 | 2 | 0 |
| **total** | **97** | **7** | **1** |

### Where this differs from the input inventory

The assessment's independent count was 36 core / 51 advanced / 18 developer.
This ledger reaches 41 core-plus-core-ingestion / 49 advanced-plus-connect /
8 adapters / 7 developer (41 + 49 + 8 + 7 = 105; the figure previously printed
here, "56 advanced-plus-connect / 8 developer", did not sum to 105 and predates
the `lm_extractor` re-tier). The material difference is the developer tier: ten
names the assessment counted as developer are re-tiered upward here because
they are **accepted arguments of exported functions** or are **named in
error-message remedies**, which makes them reachable by an ordinary user by
design rather than by accident:

`bilinear_query`, `aggregate_first`, `reduce_partitions` (arguments of
`query_geometry`/`evaluate_geometry`, `measurement_form`, `coupling`);
`query_geometry`, `evaluate_geometry`, `materialize_geometry` (named as
remedies in `R/compiler.R:33`, `R/views.R:166,180,215`); `numerical_contract`,
`numerical_agreement` (mutually parameterizing — `numerical_contract()` is the
`contract =` default of `numerical_agreement()`); `compute_policy` (the
`compute =` budget knob, used in `from-rmvpa.Rmd`).

*Review correction (2026-08-17):* `geometry_component` was listed above as
named in an error-message remedy. It is not — none of `R/compiler.R:33` or
`R/views.R:166,180,215` mentions it. Its user-facing pointer is
`print.effect_geometry_store`'s `next` hint (`R/print-methods.R:1134`), which
is a weaker but still real ground; and the `numerical_*` pair's "§4 public
list" ground is dropped, since this ledger elsewhere establishes that §4 is a
0.1 sketch rather than a binding contract.

The assessment also reports **34** exports with no vignette/README/exemplar
user; the mechanical count here is **35**. The difference is `residual_df`,
which the assessment classified as exemplar-only. The only exemplar occurrence
is `exemplars/haxby2001/results/smoke-report.md:338` — generated *output*
recording that a call **failed** ("requires a `partition` argument that its
usage was not"), not exemplar source. `vignettes/interpreting-results.Rmd:678`
uses `residual_df` as a list element *name*, not a call. `noise_precision` is
correctly exemplar-only (a real call at
`exemplars/haxby2001/05-crossnobis-uncertainty.R:91`).

---

## Ledger

### Tier: core (21)

| Name | Defining file | Users | Justification | Disposition |
|---|---|---|---|---|
| `abstract_domain()` | `R/domain.R` | V:pairing/from-obs/from-rmvpa/intro, README, T:48 | the non-volumetric neural domain constructor; every non-imaging entry starts here | keep |
| `catch_refusal()` | `R/conditions.R` | V:failures/from-rmvpa/interp/intro/novelty, README, T:31 | the reader for the package's central idiom — a refusal is a value, not a crash | keep |
| `compile_frame()` | `R/frame.R` | V:from-obs/from-rmvpa/interp/intro/neuroim2, README, E:6-coherent-configuration, T:34 | binds a spatial specification to a domain; the only route from `searchlights()`/`regions()` to a usable frame | keep |
| `compute_policy()` | `R/compute-policy.R` | V:from-rmvpa, T:15 | the declared memory/threads budget; `compute =` on `plan_geometry`, `measurement_form`, `coupling`, `plan_crossnobis` | keep |
| `contrast_energy()` | `R/views.R` | V:failures/from-obs/from-rmvpa/interp/intro/neuroim2/novelty, README, E:6-coherent-configuration, T:17 | the signed/coherent/configuration/total readout — the package's headline view | keep |
| `cross_partitions()` | `R/pairing.R` | V:failures/from-obs/from-rmvpa/interp/intro/neuroim2/novelty, README, E:2/5/6, T:39 | the default crossvalidated pairing; unbiasedness of every view depends on it | keep |
| `effect_space()` | `R/effect-space.R` | V:pairing/from-obs/from-rmvpa/interp/intro/neuroim2, README, E:2/5/6, T:26 | the experimental coordinate system every query and view is typed against | keep |
| `example_fmri_effects()` | `R/example-data.R` | V:failures/from-rmvpa/interp/intro, README, T:5 | the reproducible dataset with known spatial truth that all teaching material uses | keep |
| `lm_relation_fit()` | `R/relation-fit.R` | V:failures/from-obs/from-rmvpa/interp/intro/novelty, README, E:5-crossnobis-uncertainty, T:19 | the one-call route from betas + design to a relation with an error channel | keep |
| `plan_geometry()` | `R/geometry-plan.R` | V:failures/from-obs/from-rmvpa/interp/intro/neuroim2/novelty, README, E:2/5/6, T:33 | the query-first compiler entry; the plan is the scientific-identity object | keep |
| `rdm()` | `R/views.R` | V:corr-dist/failures/from-obs/from-rmvpa/interp/intro/neuroim2, README, E:2/5/models, T:24 | squared crossvalidated distances; the interchange format with the RSA literature | keep |
| `rdm_sampling_covariance()` | `R/evidence-sampling-product.R` | V:failures/from-obs/from-rmvpa/interp/intro, README, E:5-crossnobis-uncertainty, T:12 | exact analytic sampling covariance for distances — a first-class claim of the package | keep |
| `regions()` | `R/frame.R` | V:from-obs/from-rmvpa/neuroim2, README, E:6-coherent-configuration, T:7 | one of four frame specifications; the label-driven one | keep |
| `relation()` | `R/relation.R` | V:pairing/from-obs/from-rmvpa/interp/intro/neuroim2, README, E:2/5/6, T:33 | the lazy experimental–neural object everything downstream consumes | keep |
| `rsa()` | `R/views.R` | V:corr-dist/from-obs/from-rmvpa/interp/intro/neuroim2, README, E:2/4/5, T:12 | multiple-regression RSA compiled as one geometry query; the comparison point with rsatoolbox | keep |
| `sampling_capabilities()` | `R/evidence-sampling-product.R` | V:failures/from-obs/from-rmvpa/interp/intro/novelty, README, T:2 | ask-before-you-provoke guard for the analytic sampling law; the idiom the failure gallery teaches | keep |
| `sampling_covariance()` | `R/evidence-sampling-product.R` | V:from-obs/from-rmvpa/interp/intro, README, E:5-crossnobis-uncertainty, T:7 | the factorized sampling-covariance query, incl. the `"transport"` route WS-D/WS-E build on | keep |
| `searchlights()` | `R/frame.R` | V:from-rmvpa/interp/intro/neuroim2, README, T:13 | the frame specification the whole searchlight literature arrives expecting | keep |
| `volume_domain()` | `R/domain.R` | V:interp/intro/neuroim2, README, T:6 | the native volumetric domain constructor; no neuroim2 dependency required | keep |
| `voxelwise()` | `R/frame.R` | README†, T:13 | fourth frame specification; dropping it would leave an asymmetric set of three. The README credit is prose only (`README.md:343` lists it as an alternative to swap in); no corpus executes it | keep |
| `whole_brain()` | `R/frame.R` | V:from-obs/from-rmvpa/intro, README, T:17 | the single-node frame; the degenerate case every conservation argument is checked against | keep |

### Tier: core-ingestion (20)

The `from-observations.Rmd` spine. Not core for a user who arrives with fitted
betas; unavoidable for a user who arrives with raw data.

| Name | Defining file | Users | Justification | Disposition |
|---|---|---|---|---|
| `observations()` | `R/study-facts.R` | V:from-obs, T:4 | binds raw sources to indexes and a domain — the spine's first node | keep |
| `observation_index()` | `R/study-facts.R` | V:from-obs, T:4 | declares one partition's observation axis; row lineage starts here | keep |
| `observation_events()` | `R/study-facts.R` | V:from-obs, T:5 | the typed event record; onset/duration semantics are checked once, here | keep |
| `observation_confounds()` | `R/study-facts.R` | V:from-obs, T:4 | confound and censor facts, separated from events so censoring is auditable | keep |
| `partition_hierarchy()` | `R/study-facts.R` | V:from-obs, T:4 | declares nested partition axes; the generalization-axis claims depend on it | keep |
| `study()` | `R/study.R` | V:from-obs, README, T:4 | binds facts + clocks + axes into the one object a design compiler consumes | keep |
| `study_axis()` | `R/study.R` | T:1 | typed accessor for one declared axis. Stronger support than first recorded: the package itself directs users to it by name in two `print` methods' `next` hints (`R/print-methods.R:1237,1286`), which is the same class of package-authored instruction as an error-message remedy. 0 foreign Rd examples (corrected from 1 — `partition_hierarchy.Rd:48` is a `\seealso`) | keep |
| `study_capabilities()` | `R/study.R` | V:from-obs, T:1 | the ask-before-you-provoke guard on a study, matching the same idiom elsewhere | keep |
| `condition_space()` | `R/effect-map.R` | V:from-obs, T:3 | the semantic vocabulary that makes the coding-invariance claim expressible | keep |
| `effect_map()` | `R/effect-map.R` | V:from-obs, T:5 | declares effects in condition-space terms; the input to lowering | keep |
| `coefficient_parameterization()` | `R/effect-map.R` | V:from-obs, T:3 | names one concrete coefficient basis; the other half of the lowering pair | keep |
| `lower_effect_map()` | `R/effect-map.R` | T:1 | the inspectable step of the semantic route; the observable behind the coding-invariance claim; `plan_relation()` calls it at `R/relation-plan.R:121` (verified). 0 foreign Rd examples (corrected from 3: `raw_effect_map.Rd:49` is a `#` comment inside an example, the other two are `\seealso`) | keep — see reviewer note |
| `design_model()` | `R/design-model.R` | V:from-obs/intro, T:4 | declares a semantic design with compiled routes; the coding-invariant route | keep |
| `raw_design_model()` | `R/design-model.R` | V:from-obs†, T:3 | the escape hatch for pre-built designs; reports `coding_invariant = FALSE` — merging into `design_model()` would hide exactly that capability difference. Vignette credit is prose only (`from-observations.Rmd:601`) | keep |
| `raw_effect_map()` | `R/effect-map.R` | V:from-obs†, T:3 | the paired escape hatch for a raw target matrix; same capability-honesty argument. Vignette credit is prose only | keep |
| `observation_model()` | `R/observation-model.R` | V:from-obs, T:5 | declares the first-moment observation model; separates estimand from estimator | keep |
| `plan_relation()` | `R/relation-plan.R` | V:from-obs/intro, README, T:5 | plans the experimental–neural relation; produces the portable design receipt | keep |
| `relation_plan_receipts()` | `R/relation-plan.R` | V:from-obs, T:2 | reads the portable receipts; the evidence that a third-party compiler conformed | keep |
| `compiler_conformance()` | `R/compiler-conformance.R` | V:from-obs, T:3 | the conformance court for external design compilers; the ingestion contract's admission test | keep |
| `estimate_relation()` | `R/relation-plan.R` | V:from-obs/interp, README, T:4 | executes a planned relation family; the spine's terminal node | keep |

> **Reviewer note (2026-08-17) — `lower_effect_map()`: proposed
> `demote-internal`, or document it.** It satisfies none of the four standing
> keep-rules: it is not an accepted argument value of any exported function, it
> is named in no error-message remedy, it is not a compatibility promise in
> `ingestion-contract.md`, and it is not a developer entry point. `T:1`, and
> its foreign-Rd support was 0 once `\seealso` and comments are excluded. The
> real argument for it — that it is the *observable* of the coding-invariance
> claim — is a good one, but nothing in the tree demonstrates it. Either
> `from-observations.Rmd` gains a chunk that lowers one effect map two ways and
> shows the invariance (the better fix, since coding invariance is a headline
> claim), or the function goes internal and `plan_relation()` keeps calling it.
> Keeping it exported *and* undemonstrated is the one option a subtraction
> release should not take.
>
> **Maintainer decision 4 (2026-08-17): keep — the documentation fix, not the
> demotion.** The reviewer's framing is accepted in full, including that the
> third option is not available. The obligation is therefore attached to a
> ticket: **A4 adds a runnable vignette chunk demonstrating the
> coding-invariance observable** — lowering one effect map under two coefficient
> parameterizations and showing the invariant result — in
> `from-observations.Rmd` or `interpreting-results.Rmd`. Coding invariance is a
> headline claim of the package and currently has no executed demonstration
> anywhere; this row is where that gap is paid off. If A4 ships without the
> chunk, `lower_effect_map()` reverts to `demote-internal` and `plan_relation()`
> keeps calling it at `R/relation-plan.R:121`.

> **Reviewer note (2026-08-17) — tier boundary.** `compiler_conformance()` and
> `relation_plan_receipts()` are tiered core-ingestion but read as **developer**
> under this document's own tier definition ("extension and adapter *packages*,
> not end users"): they are the admission test and receipt reader for an
> *external design compiler*, which is by construction another package. They
> are kept here because `from-observations.Rmd` genuinely executes them
> (lines 254, 265, 350 — real chunks, not prose). Flagged because if the
> developer tier is defined by *who calls it*, the tier is larger than its
> stated count and the WS-A target is met by where the line was drawn rather
> than by subtraction.
>
> **Maintainer decision 7 (2026-08-17): both stay exported, and the tier
> definition is reworded rather than the rows re-tiered.** The developer tier is
> now defined as the extension-***only*** surface — what an extension package
> meets *and an end user never does* — and the tier definitions section states
> explicitly that an extension package also meets the core-ingestion tier
> (about fourteen functions in total for the four in-tree adapters, of which
> five are developer-tier). Under that definition `compiler_conformance()` and
> `relation_plan_receipts()` are correctly outside it: `from-observations.Rmd`
> executes both, so they are not extension-only, whatever the caller's package.
> The published number is now true of the tier *and* answers the reader's
> question, which was the reviewer's actual complaint.

### Tier: advanced (33)

| Name | Defining file | Users | Justification | Disposition |
|---|---|---|---|---|
| `additive_frame()` | `R/scope.R` | V:pairing/interp, T:17 | the general diagonal frame constructor; WS-D's frame families are built on it | keep |
| `aggregate_first()` | `R/operations.R` | T:4 | the documented default `reducer =` of `measurement_form()` and `coupling()`; an exported argument value | keep |
| `reduce_partitions()` | `R/operations.R` | T:5 | the other admissible `reducer =`; the edge-first alternative that changes the estimand | keep |
| `inner_product()` | `R/operations.R` | T:3 | **no exported function accepts a normalizer**; its siblings `covariance()`/`cosine()`/`correlation()` (`R/operations.R:48,50,54`) are already private for exactly this reason — only `inner_product` appears in `NAMESPACE:253`; 0 foreign Rd examples | **demote-internal** |
| `bilinear_query()` | `R/scope.R` | T:11 | the only way a user supplies a custom operator `H` to `query_geometry()`/`evaluate_geometry()`; called in 4 foreign Rd examples | keep |
| `pair_query()` | `R/scope.R` | V:pairing/intro/novelty, T:14 | the axis-bound pair query — the rectangular (unequal-axes) entry point | keep |
| `pairing()` | `R/pairing.R` | V:pairing/novelty, T:29 | explicit partition pairing when `cross_partitions()` is not the intended design | keep |
| `materialize_geometry()` | `R/geometry-entry.R` | V:novelty†, T:13 | builds the complete geometry when a user wants many queries against one run; named in four error messages (`R/compiler.R:33`, `R/views.R:126,166,196`), one of them a `remedies` field. Vignette credit is prose only — `novelty.Rmd` has no code | keep |
| `evaluate_geometry()` | `R/geometry-entry.R` | T:7 | the query-first evaluation that never materializes; named as a remedy in `R/views.R:180,215`; 2 foreign Rd examples (corrected from 5) | keep |
| `query_geometry()` | `R/result.R` | T:11 | applies a linear query to an already-materialized geometry; `R/views.R:134-135` and `R/compiler.R:33` both route users here by name | keep |
| `geometry_component()` | `R/result.R` | T:11 | reads one component (total/coherent/configuration) of a complete geometry. **Not** named in any error-message remedy (corrected); its user-facing pointer is the `next` hint of `print.effect_geometry_store` (`R/print-methods.R:1134`), plus 1 foreign Rd example (`materialize_geometry.Rd`) | keep — see reviewer note |
| `geometry_spectrum()` | `R/views.R` | V:failures†, T:5 | the signed eigenvalue spectrum; the diagnostic the failure gallery *describes* to show indefiniteness. The vignette credit is prose only (`failure-gallery.Rmd:84`); no chunk runs it | keep |
| `frame_conservation()` | `R/frame.R` | V:novelty†, T:2 | the conservation diagnostic on a compiled frame; WS-D's certificate extends it. Vignette credit is prose only — `novelty.Rmd:399` is a status-table cell | keep — see reviewer note |
| `neural_metric()` | `R/metric.R` | V:failures, T:6 | the fixed support-local metric constructor | keep |
| `noise_precision()` | `R/crossnobis.R` | E:5-crossnobis-uncertainty, T:7 | the fixed noise-precision metric; fixed-metric crossnobis is exactly `plan_geometry()` + this | keep |
| `identity_metric()` | `R/metric-learning.R` | T:2 | the null-metric recipe, and an admissible `metric =` value of `plan_crossnobis()`/`plan_geometry()` (`effect_metric_recipe`, `R/metric.R:363,436`). The original justification had the dependency backwards: it is `identity_metric`'s own example that calls `metric_capabilities()`, not the reverse; 0 foreign Rd examples (corrected from 1) | keep |
| `diagonal_precision()` | `R/metric-learning.R` | T:1 | univariate noise normalization recipe with an explicit variance floor; an admissible `metric =` value alongside `identity_metric()`/`shrinkage_precision()`. 0 foreign Rd examples (corrected from 3 — all three were `\seealso` links) | keep |
| `shrinkage_precision()` | `R/metric-learning.R` | V:from-rmvpa, E:5-crossnobis-uncertainty, T:3 | the default learned recipe; `plan_crossnobis()`'s `metric =` default | keep |
| `metric_training_policy()` | `R/metric-learning.R` | V:from-rmvpa†, T:5 | declares which residual partitions may train a metric; the leakage guard. Vignette credit is prose only (`from-rmvpa.Rmd:378`) | keep |
| `metric_capabilities()` | `R/metric.R` | T:3 | the metric-side capability inspector: the fourth member of the `*_capabilities()` inspector family (`study_`, `relation_fit_`, `sampling_`), which is one idiom and should not be demoted piecemeal. It is *also* called in the `\examples{}` of `identity_metric`, `neural_metric`, `shrinkage_precision` (3 foreign examples, verified) — but per the standing rule that is a cost of demotion, not a reason to keep | keep |
| `plan_crossnobis()` | `R/crossnobis.R` | V:from-rmvpa, E:5-crossnobis-uncertainty, T:4 | the learned-metric plan | keep — WS-B may delete it in favour of `plan_geometry(metric = shrinkage_precision())`; that is a separate ticket and a breaking change, not a tiering decision |
| `crossnobis()` | `R/crossnobis-driver.R` | V:from-rmvpa, E:5-crossnobis-uncertainty, T:6 | the signed local crossnobis contrast; WS-B keeps it as a validating view | keep |
| `lm_extractor()` | `R/extractor.R` | T:5 | compiles a supplied design into the extraction map `E` with pivoted-QR / SVD rank handling. **Re-tiered from developer to advanced (maintainer decision 2, 2026-08-17):** it is an accepted `extract =` argument value of `relation()` and the ordinary user path in `crossform-package-design.md`'s "Raw response example" — a specialist's tool, not an extension-only seam. Its `ingestion-contract.md:481` (§13) compatibility promise binds regardless of tier; 0 foreign Rd examples | keep |
| `relation_fit_capabilities()` | `R/relation-fit.R` | T:3 | the "does this object have an error channel" guard; 5 foreign Rd examples, and it is the public alternative to poking at the fit's internals | keep |
| `residual_block()` | `R/relation-fit.R` | T:7 | the documented user path to residuals, and the carrier of the "divide by `residual_df()`, not `nrow()`" instruction (`R/relation-fit.R:721`); user-facing science, not extension plumbing | keep |
| `residual_df()` | `R/relation-fit.R` | T:6 | the correct divisor for any noise-variance estimate; also the canonical refusal demo in `relation()`'s examples (`R/relation.R:152`) | keep |
| `numerical_contract()` | `R/numerics.R` | T:4 | the `contract =` argument of `numerical_agreement()` and its default (verified, `R/numerics.R:93`); the package's reproducibility promise as a value | keep |
| `numerical_agreement()` | `R/numerics.R` | T:3 | the user-facing verifier that two runs agree under a named guarantee, and the only public reader of the contract; carries its own worked `\examples{}` block (`man/numerical_agreement.Rd:39-46`). The §4 citation is dropped: this ledger already established that §4 is a 0.1 sketch, not a binding contract | keep |
| `match_coupling()` | `R/pair-query.R` | T:2 | marks matched encoding–retrieval pairs — the first step of the ER-RSA family | keep, conditional on F2 — see below |
| `control_coupling()` | `R/pair-query.R` | T:2 | marks the eligible control cells against those matches | keep, conditional on F2 — see below |
| `coupling_contrast()` | `R/pair-query.R` | T:1 | builds the balanced matched-minus-control operator with its additive-baseline diagnostics | keep, conditional on F2 — see below |
| `match_control()` | `R/pair-query.R` | T:1 | the nuisance-adjusted LM form of the same contrast (item/probe nuisance families on by default) | keep, conditional on F2 — see below |
| `pair_lm_query()` | `R/pair-query.R` | T:1 | the general weighted pair-space LM coefficient, compiled once and reusable | keep, conditional on F2 — see below |

The five `pair-query.R` exports are a single coherent ER-RSA API with no
vignette, README or exemplar user. Their designated user is the **rectangular
ER-RSA exemplar named as pending** in `vignettes/novelty.Rmd:34,395,421` and
scheduled in this program as **ticket F2**. Five exports resting on one
unwritten exemplar is the weakest justification in this ledger, which is why
the keep is conditional rather than settled.

> **Decision (maintainer decision 3, 2026-08-17) — the ER-RSA five: keep,
> conditional on F2.** All five stay exported for this release, justified by
> F2's rectangular ER-RSA exemplar rather than by any current user. The
> condition is binding and mechanical: **any of the five that the delivered F2
> exemplar does not exercise is demoted to internal at F2's close.** That is an
> acceptance criterion of both **A4** and **F2** — F2 does not close until the
> exemplar is runnable and the un-exercised subset has been demoted, and A4's
> documentation pass must not present any of the five as settled public API
> before then. The reviewer's alternative below (demote `match_control()` and
> `pair_lm_query()` now) is declined only because F2 is scheduled inside this
> program; if F2 leaves the program, that alternative is the fallback and this
> block is the first to reconsider.

> **Reviewer note (2026-08-17) — the ER-RSA five: proposed `demote-internal`
> for `match_control()` and `pair_lm_query()` now, and for the remaining three
> if the exemplar does not land in this release. DECLINED in that form; the
> substance is carried by the conditional keep above (maintainer decision 3),
> which applies the reviewer's rule to all five at F2's close instead of to two
> of them now.** The ledger's own assessment ("the
> weakest justification in this ledger") is correct and, on the corrected
> numbers, understated: `coupling_contrast` has **0** foreign Rd examples, not
> 4; `match_control` **0**, not 3; `pair_lm_query` **0**, not 2. All five fail
> every standing keep-rule — none is an accepted argument value of an exported
> function (they are arguments to *each other*, which is self-referential), none
> is named in a remedy, none is in the ingestion contract, none is a developer
> entry. `novelty.Rmd` has no code fences, so the "pending exemplar" is named
> in prose in a vignette that demonstrates nothing. A pending exemplar is a
> plan, not a user; under a subtraction posture a plan does not hold an export
> open. If the maintainer wants to preserve the family's shape pending F2,
> the minimum defensible subtraction is `match_control()` and `pair_lm_query()`
> — two alternative terminal forms of a contrast that already has a terminal in
> `coupling_contrast()`, each at `T:1` with no other support.

> **Reviewer note (2026-08-17) — `geometry_component()`: keep, justification
> rebuilt.** The original justification was wrong in both halves. It is named in
> **no** error-message remedy (`R/compiler.R:33` and `R/views.R:166,180,215`
> name `materialize_geometry`, `query_geometry`, `evaluate_geometry` and
> `pair_query`, never `geometry_component`), and the "3 foreign Rd examples
> incl. `as_neurovol`" is 1 (`materialize_geometry.Rd`); the `as_neurovol`
> mention is an `@seealso` at `R/neuroim2-adapter.R:343`. I still say keep,
> because `print.effect_geometry_store` emits
> `next = geometry_component(geometry, "total")` (`R/print-methods.R:1134`) —
> a package-authored instruction naming the function, which is the same kind of
> obligation as a remedy. If the maintainer declines to extend the standing
> rule to `print` hints, this row becomes a `demote-internal` candidate and
> `study_axis()` goes with it, since both rest on exactly that ground.
>
> **Maintainer decision 5 (2026-08-17): the rule is extended. `print` `next:`
> hints count as remedies, and both rows are kept.** The `next:` hint is the
> package's on-ramp mechanism — the designed route from an object a user is
> holding to the function that reads it — so a name the package prints as the
> user's next call is package-authored instruction in exactly the sense an
> error-message remedy is. Demoting a function the package tells the user to
> call would leave `print` emitting a pointer to a `\keyword{internal}` topic.
> This settles `geometry_component()` and `study_axis()` together.

> **Reviewer note (2026-08-17) — `frame_conservation()`: subtraction
> candidate.** With the vignette credit corrected to prose-only, its entire
> support is `T:2` and a forward reference to WS-D. Not proposing a disposition
> flip, but it is now in the same evidentiary position as the ER-RSA five and
> should be judged by the same rule.
>
> **Maintainer decision 6 (2026-08-17): keep.** The distinction from the ER-RSA
> five is that WS-D's conservation certificate gives `frame_conservation()` real
> in-tree users inside this program — it is the diagnostic the certificate is
> built on, not a function waiting for an exemplar to justify it. A4 must give
> it an executed example.

### Tier: advanced (experimental: connect) (16)

Measurement forms, coupling closures, and tomography. `R/evidence-api.R:10-31`
and `R/tomography.R:115-122` gate this tier: a 256 MiB dense small-node ceiling
and an explicit refusal to claim brain-scale coupling. Every reference page in
this tier must carry the experimental label.

| Name | Defining file | Users | Justification | Disposition |
|---|---|---|---|---|
| `measurement_frame()` | `R/evidence-api.R` | V:pairing, T:5 | declares the identified small-node measurements a form is built over | keep |
| `edge_frame()` | `R/measurement.R` | V:pairing, T:4 | declares which measurement edges are requested; the cost bound is computed from it | keep |
| `measurement_form()` | `R/evidence-api.R` | V:pairing, T:4 | materializes the requested measurement-space forms; the tier's compiler entry | keep |
| `measurement_components()` | `R/evidence-api.R` | V:pairing, T:1 | summarizes the crossed node decomposition on one edge | keep |
| `variation_query()` | `R/evidence-api.R` | V:pairing, T:5 | declares a repeated-variation experimental query; the `by =` of `coupling()` | keep |
| `coupling()` | `R/evidence-api.R` | V:novelty†, T:3 | the adjoint-side coupling closure of a geometry plan; the tier's headline construction. Its only vignette credit is prose (`novelty.Rmd:396`, a status-table cell) — the tier's headline construction is nowhere executed outside tests | keep — see reviewer note |
| `effect_coupling()` | `R/evidence-api.R` | V:pairing, T:5 | reads a completed form as an effect coupling | keep |
| `covariance_coupling()` | `R/evidence-api.R` | V:pairing†, T:2 | reads it as a covariance; a scientifically distinct reading, not a `kind` flag. Vignette credit is prose only (`evidence-pairing.Rmd:509`, a capability-table cell) | keep — see reviewer note |
| `canonical_coupling()` | `R/evidence-api.R` | V:pairing, T:2 | reads it as a ridge-regularized canonical coupling; takes its own `ridge` arguments | keep |
| `geometry_alignment()` | `R/evidence-api.R` | V:pairing, T:1 | reads it as an alignment with a tolerance; fourth distinct reading | keep |
| `gaussian_covariance_model()` | `R/evidence-api.R` | V:pairing, T:3 | declares the joint Gaussian interpretation that licenses the covariance reading. (`T:2` corrected to `T:3`: `test-coupling-views.R`, `test-evidence-api.R`, `test-print-methods.R`) | keep |
| `connectivity()` | `R/evidence-api.R` | V:pairing, T:2 | the normalized connectivity view, gated behind explicit capabilities | keep |
| `reconstruct_evidence()` | `R/evidence-api.R` | V:pairing, T:1 | the only public entry to tomography (`R/evidence-api.R:819-832`); returns `effect_tomography_result` or a budget refusal | keep |
| `measurement_space()` | `R/bridge.R` | T:5 | **no exported function accepts a measurement space or a bridge** — bridges enter the IR only through internal `.effect_evidence_task()` constructors (`R/evidence-task.R:148,179,864,946`); tests already reach `crossform:::.apply_measurement_bridge()` for the consuming half | **demote-internal** |
| `measurement_bridge()` | `R/bridge.R` | T:4 | same: an exported value type with no exported consumer. Keep the type and its validator; re-export when WS-E's population form gives it a consumer | **demote-internal** |
| `reverse_bridge()` | `R/bridge.R` | T:2 | three lines that rebuild the bridge with its legs swapped (`R/bridge.R:236-244`); as an internal it is an argument, not a name. 0 foreign Rd examples (corrected from 1) | **merge-into:`measurement_bridge`** (as `reverse = TRUE`) |

> **Reviewer note (2026-08-17) — the connect tier's headline is
> undemonstrated.** `coupling()` and `covariance_coupling()` were credited with
> a vignette user; both credits are prose table cells. Combined with
> `measurement_space`/`measurement_bridge`/`reverse_bridge` already going
> internal, this tier's evidentiary base is thinner than the counts suggest:
> of its 16 exports, the two the ledger calls "the tier's headline
> construction" and a "scientifically distinct reading" are executed nowhere
> outside `tests/`. I am not proposing disposition flips — the tier is labelled
> experimental and that label does real work — but the experimental label
> should be read as *the reason* these keeps are admissible, not as an
> afterthought, and `_pkgdown.yml` should not present them at the same
> confidence as `contrast_energy()` or `rdm()`.
>
> **Maintainer decision 6 (2026-08-17): `coupling()` and `covariance_coupling()`
> are kept, and the reviewer's reading is adopted as the ground.** They are the
> connect tier's entry points — demoting the tier's headline construction and
> its second distinct reading would leave a tier a user cannot enter — and
> **ticket A9 marks the whole tier experimental**, which is the label that makes
> a thin evidentiary base admissible here and nowhere else in the ledger.
> Because the keep rests on the label rather than on demonstrated use, it comes
> with an obligation: **A4 must give `coupling()` and `covariance_coupling()`
> each an executed example** (Rd `\examples{}` at minimum), so that neither
> remains a function no corpus outside `tests/` ever runs.

> **Reviewer note (2026-08-17) — the three demotions in this tier are correct
> and their reasons verify.** Bridges enter the IR only through internal
> `.effect_evidence_task()` constructors (`R/evidence-task.R:58,483,879,990`);
> the two exported functions that mention the class, `metric_capabilities()`
> (`R/metric.R:427`) and `metric_components()` (`R/metric.R:916`), *reject* a
> bridge with an error rather than consuming one. "No exported consumer" holds.

### Tier: adapters (8)

Labeled sub-tier of advanced. Each is a morphism into the typed core; none of
their dependencies is mandatory. Note that the three `neuroim2-adapter.R`
entries are de-facto core for volumetric work (3 vignettes + README + 3
exemplar scripts each) and should be presented alongside core in `_pkgdown.yml`
even though they are tiered here by defining file.

| Name | Defining file | Users | Justification | Disposition |
|---|---|---|---|---|
| `neuroim2_volume_domain()` | `R/neuroim2-adapter.R` | V:from-rmvpa/intro/neuroim2, README, E:2/5/6, T:2 | the mask → domain morphism; every real-data path begins here | keep |
| `neuroim2_searchlights()` | `R/neuroim2-adapter.R` | V:from-rmvpa/intro/neuroim2, README, E:2/5/6, T:1 | compiles neuroim2 searchlight indices into a crossform frame | keep |
| `as_neurovol()` | `R/neuroim2-adapter.R` | V:from-rmvpa/intro/neuroim2, README, T:2 | maps a compact result vector back to a volume — the last step of every imaging workflow | keep — becomes an S3 generic so other packages can hook it (separate ticket; currently a plain function at `R/neuroim2-adapter.R:274`) |
| `bids_study()` | `R/adapter-bids.R` | V:from-obs, T:1 | binds BIDS files into a generic study; composes the two file adapters below | keep |
| `bids_events()` | `R/adapter-bids.R` | T:1 | the events-only adapter. Fails all four standing keep-rules: not an accepted argument value (`bids_study()` takes *file paths*, not the record it returns), named in no error-message remedy, absent from `ingestion-contract.md`'s compatibility section, not a developer entry. 0 foreign Rd examples (corrected from 3 — all three are `\seealso`), leaving `T:1` as its entire support. `bids_study()` composes it at `R/adapter-bids.R:267` (verified), so every documented workflow reaches it only through `bids_study()` | **demote-internal** (maintainer decision 1) — cost is one test file; re-export the day a partial-BIDS workflow is written |
| `bids_confounds()` | `R/adapter-bids.R` | T:1 | the fMRIPrep confound adapter; identical position on all four rules. 0 foreign Rd examples (corrected from 3), `T:1`, composed by `bids_study()` at `R/adapter-bids.R:272` | **demote-internal** (maintainer decision 1) — same cost, same re-export condition |
| `fmridesign_design_model()` | `R/adapter-fmridesign.R` | V:from-obs/intro†, T:1 | compiles a semantic design model from fmridesign; an admitted external compiler. Both vignette credits are prose only | keep |
| `fmrireg_relation()` | `R/adapter-fmrireg.R` | V:from-obs†, T:2 | executes a relation plan with fmrireg. Vignette credit is prose only | keep |

> **Reviewer note (2026-08-17) — `bids_events()` / `bids_confounds()`:
> proposed `demote-internal` for both. ACCEPTED (maintainer decision 1).** Their foreign-Rd support was 3 each;
> corrected, it is **0 each** — every one of those six references is a
> `\seealso` link. That removes the only evidence beyond `T:1`. They satisfy no
> standing keep-rule: not argument values (`bids_study()` takes *file paths*,
> not the records they return), not named in any remedy, not in
> `ingestion-contract.md`'s compatibility section, not developer entries. The
> ledger's own open question concedes "the partial-BIDS case, which no test or
> vignette currently exercises" — that is a hypothesis about a user, and this
> release is supposed to subtract exactly those. `bids_study()` calls both
> internally (`R/adapter-bids.R:267,272`), so demoting costs one test file's
> switch to `crossform:::` and nothing else; re-export the day a partial-BIDS
> workflow is written. This is the clearest available −2 in the whole ledger.

### Tier: developer (7 → 5 exported)

The extension-*only* protocol for handing data, a design, and an error channel
into the core. The sanctioned set is the "hand something in" surface plus the
bounded block read; the two accumulators that only the compiler consumes go
internal. `lm_extractor()` was tiered here until maintainer decision 2
(2026-08-17) moved it to **advanced**: it is an `extract =` value of
`relation()`, which makes it reachable by an ordinary user by design, so it is
not extension-only. That re-tier takes the tier from 8 rows to 7 and the
sanctioned set from 6 to **5**.

| Name | Defining file | Users | Justification | Disposition |
|---|---|---|---|---|
| `file_matrix_source()` | `R/source.R` | T:4 | the out-of-memory source constructor — the reason `relation()` is lazy at all; §4 public list | keep (sanctioned 1/5) |
| `source_capabilities()` | `R/receipt.R` | T:12 | the capability declaration every custom source must supply; the compiler refuses on it, so an extension author cannot avoid it | keep (sanctioned 2/5) |
| `effect_extractor()` | `R/extractor.R` | T:2 | the declared map `E` in `B = E Y`; named as a compatibility promise in `ingestion-contract.md:481` | keep (sanctioned 3/5) |
| `relation_fit()` | `R/relation-fit.R` | T:3 | attaches an externally computed error channel to a relation; the seam for a package that fits its own GLM | keep (sanctioned 4/5) |
| `relation_block()` | `R/relation.R` | T:15 | the bounded block-read protocol — "where neural values are actually read" (`R/relation.R:398`); §4 public list; 3 foreign Rd examples | keep (sanctioned 5/5) |
| `effect_covariance()` | `R/relation-fit.R` | T:2 | the design-side factor of the separable error model. Stronger than originally stated: `effect_covariance()` has **no call site anywhere in `R/` outside its own defining file** — `rdm_sampling_covariance()` does not call it (`R/evidence-sampling-kernel.R:219` takes `effect_covariance` as a *parameter name*, not a call). Its only callers are three test lines. `relation_fit_capabilities()` already answers the only question a user asks of it. 0 foreign Rd examples (corrected from 2) | **demote-internal** |
| `residual_pair_statistics()` | `R/residual-statistics.R` | T:5 | the residual sufficient-statistic accumulator that `plan_crossnobis()` runs internally (verified: `R/evidence-sampling-product.R:156`, `R/crossnobis.R:340`); a compiler stage exposed for testing, not an API. 0 foreign Rd examples (corrected from 1) | **demote-internal** |

> **Reviewer note (2026-08-17) — developer-set verdict: the sanctioned entries
> are the right *seams*, but the count does not mean what the tier definition
> says it means.** (Written against the pre-decision set of six; both of its
> recommendations were adopted — see maintainer decisions 2 and 7. Counts below
> are restated for the post-decision set of five.)
> I checked the claim against what an adapter actually calls. The four in-tree
> adapters call, in total:
>
> | Adapter | Core functions it calls |
> |---|---|
> | `R/adapter-bids.R` | `observation_events()`, `observation_confounds()`, `study()` |
> | `R/adapter-fmridesign.R` | `condition_space()`, `coefficient_parameterization()`, `design_model()` |
> | `R/adapter-fmrireg.R` | `effect_extractor()`, `relation()`, `relation_fit()` |
> | `R/neuroim2-adapter.R` | `abstract_domain()` |
>
> Only **2 of the 5** sanctioned developer entry points (`effect_extractor`,
> `relation_fit`) are exercised by any adapter in this tree. The rest of what an
> adapter needs is tiered core-ingestion. Two consequences:
>
> 1. **The set is defensible.** The five cover three genuinely
>    extension-only seams — data-in (`file_matrix_source`, `source_capabilities`,
>    `relation_block`), design-in (`effect_extractor`),
>    error-channel-in (`relation_fit`) — and no in-tree adapter is stranded by
>    the choice. The three unexercised ones are unexercised because no in-tree
>    adapter supplies its own out-of-memory source, not because the seam is
>    imaginary. I would cut none of the five.
> 2. **The tier definition oversells it.** "≤ 5 sanctioned entry points" reads
>    as *an extension package meets 5 crossform functions*. It meets about 14.
>    The sentence should say the developer tier is the extension-*only* surface,
>    and that adapters additionally use the core-ingestion tier — otherwise the
>    subtraction release publishes a number that is true of the tier and false
>    of the reader's question. **Adopted:** the tier definition now says exactly
>    this (maintainer decision 7).
>
> One re-tier worth making — **adopted as maintainer decision 2**:
> **`lm_extractor()` is not a developer function.** It
> is an accepted `extract =` argument value of `relation()`, `§4` lists it under
> *Relations*, and `crossform-package-design.md`'s own "Raw response example"
> uses it as the ordinary user path
> (`relation(bold_runs, extract = lm_extractor(...))`). Moving it to **advanced**
> costs nothing (it stays exported, and the `ingestion-contract.md` §13
> compatibility promise binds regardless of tier) and leaves a developer tier of
> **5** that is honestly extension-only. Its stated justification here — "named
> in `ingestion-contract.md:481`" — verifies (§13, confirmed), but that promise
> is about the API *remaining valid*, not about who calls it.

---

## Tests-only exports (35): dispositions

Thirty-five exports have no call in any vignette, the README, or any exemplar
script. Being tests-only is a *flag*, not a verdict: five of these are argument
values of exported functions, several are called in another topic's checked
`\examples{}`, and one family is waiting on a scheduled exemplar. Each is
dispositioned explicitly below with the reason it survives or does not.

**Two corrections from the 2026-08-17 review apply to this whole table.**
First, the *Foreign Rd examples* column has been re-measured with
`tools::Rd2ex` and comments stripped; 29 of the 35 original counts were
overstated because they counted `\seealso` `\link{}` cross-references,
`\usage{}` defaults, `\description{}` prose, and `#` comments inside example
blocks. The corrected column below is the number of *other* topics whose
`\examples{}` actually call the name. Second, on the corrected numbers
**22 of the 35 have zero foreign-example use**, so the `inner_product` row's
claim to be "the only export with zero foreign-example use" is withdrawn — it
was already contradicted by this table's own `relation_fit` row.

Third, and separately: this table's membership is itself understated. Under a
code-chunk-only reading of the corpora (see the Method caveat), eleven further
exports have no executable user outside `tests/` — `coupling`,
`covariance_coupling`, `fmridesign_design_model`, `fmrireg_relation`,
`frame_conservation`, `geometry_spectrum`, `materialize_geometry`,
`metric_training_policy`, `raw_design_model`, `raw_effect_map`, `voxelwise` —
making the true count **46**. They are not added as rows here because their
dispositions are settled in their tier tables, but a subtraction release should
know the number is 46.

| Name | Tests | Foreign Rd examples | Disposition | Reason |
|---|---:|---:|---|---|
| `aggregate_first` | 4 | 1 | keep | default `reducer =` of `measurement_form()`/`coupling()` — an exported argument value (verified `R/evidence-api.R:375,993`) |
| `reduce_partitions` | 5 | 0 | keep | the other admissible `reducer =` (`.validate_partition_reducer`, `R/pairing.R:287`); changes the estimand, so it must be nameable |
| `inner_product` | 3 | 0 | **demote-internal** | no exported function accepts a normalizer; its three siblings `covariance()`/`cosine()`/`correlation()` are already private (`R/operations.R:48,50,54`) |
| `bilinear_query` | 11 | 3 | keep | the `query =` argument of `query_geometry()`/`evaluate_geometry()` |
| `evaluate_geometry` | 7 | 2 | keep | named as a remedy in two error messages (`R/views.R:180,215`) |
| `query_geometry` | 11 | 0 | keep | `R/compiler.R:33` routes users to it by name — an error-message remedy |
| `geometry_component` | 11 | 1 | keep | *not* a remedy (correction); kept on `print.effect_geometry_store`'s `next` hint, `R/print-methods.R:1134` — **decision 5** rules that `print` `next:` hints count as remedies |
| `numerical_contract` | 4 | 0 | keep | the `contract =` argument and default of `numerical_agreement()` (`R/numerics.R:93`) |
| `numerical_agreement` | 3 | 0 | keep | the only user-facing reader of the reproducibility contract; carries its own worked `\examples{}` |
| `metric_capabilities` | 3 | 3 | keep | fourth member of the `*_capabilities()` inspector family; demoting one of four would break the idiom |
| `identity_metric` | 2 | 0 | keep | an admissible `metric =` recipe value (the original "`metric_capabilities()`'s example contrasts against it" was reversed) |
| `diagonal_precision` | 1 | 0 | keep | third member of the recipe family; an admissible `metric =` value of `plan_crossnobis()`/`plan_geometry()` |
| `relation_fit_capabilities` | 3 | 4 | keep | the public "is there an error channel" question; member of the inspector family |
| `residual_block` | 7 | 1 | keep | user-facing residual read; carries the "divide by `residual_df()`" instruction (`R/relation-fit.R:721`) |
| `residual_df` | 6 | 4 | keep | the correct noise-variance divisor; the refusal demo in `relation()`'s examples (`R/relation.R:152`) |
| `effect_covariance` | 2 | 0 | **demote-internal** | no call site anywhere in `R/` outside its own file; `rdm_sampling_covariance()` does not call it |
| `residual_pair_statistics` | 5 | 0 | **demote-internal** | compiler-stage accumulator run inside `plan_crossnobis()` (`R/crossnobis.R:340`) |
| `effect_extractor` | 2 | 0 | keep | sanctioned developer entry; compatibility promise in `ingestion-contract.md` §13 (verified) |
| `lm_extractor` | 5 | 0 | keep | compatibility promise in `ingestion-contract.md` §13; **re-tiered developer → advanced (decision 2)** — it is an `extract =` value of `relation()`, not an extension-only seam |
| `relation_fit` | 3 | 0 | keep | sanctioned developer entry — the only way an external fitter attaches an error channel; the sole developer entry an in-tree adapter actually calls (`R/adapter-fmrireg.R:115`) |
| `relation_block` | 15 | 2 | keep | sanctioned developer entry; the bounded read protocol |
| `file_matrix_source` | 4 | 0 | keep | sanctioned developer entry; the out-of-memory source |
| `source_capabilities` | 12 | 0 | keep | sanctioned developer entry; the compiler refuses without it |
| `measurement_space` | 5 | 3 | **demote-internal** | no exported consumer of a measurement space or bridge exists (verified) |
| `measurement_bridge` | 4 | 2 | **demote-internal** | same; re-export when WS-E gives it a consumer |
| `reverse_bridge` | 2 | 0 | **merge-into:`measurement_bridge`** | a leg swap, not a concept |
| `match_coupling` | 2 | 3 | keep, conditional | ER-RSA family; user is the rectangular exemplar scheduled as **ticket F2** — **decision 3**: demoted at F2's close if the delivered exemplar does not exercise it |
| `control_coupling` | 2 | 1 | keep, conditional | ER-RSA family; same F2 condition (decision 3) |
| `coupling_contrast` | 1 | 0 | keep, conditional | ER-RSA family; same F2 condition (foreign count was 4, is 0) |
| `match_control` | 1 | 0 | keep, conditional | ER-RSA family; same F2 condition — reviewer's "demote now" declined in favour of the F2 test |
| `pair_lm_query` | 1 | 0 | keep, conditional | ER-RSA family; same F2 condition — same |
| `bids_events` | 1 | 0 | **demote-internal** | **decision 1**: fails all four standing rules, 0 foreign Rd examples, composed by `bids_study()` (`R/adapter-bids.R:267`); cost is one test file |
| `bids_confounds` | 1 | 0 | **demote-internal** | **decision 1**: identical position; `bids_study()` composes it at `R/adapter-bids.R:272` |
| `lower_effect_map` | 1 | 0 | keep, with obligation | the observable of the coding-invariance claim — **decision 4**: kept on condition that **A4** adds a runnable vignette chunk demonstrating it |
| `study_axis` | 1 | 0 | keep | typed axis accessor named by the package itself in two `print` `next` hints (`R/print-methods.R:1237,1286`) — **decision 5** makes those hints a keep-rule |

## What demotion costs

Every demotion below is a `\keyword{internal}` + un-export, not a deletion. No
vignette calls any of the eight affected names — `from-observations.Rmd` reaches
BIDS only through `bids_study()` — so **no vignette needs to change** and no
vignette needs `crossform:::` (which would not be permitted in a vignette in any
case). The costs are confined to tests and to roxygen `@examples`.

| Demoted | Tests that switch to `crossform:::` | Roxygen `@examples` that must change |
|---|---|---|
| `inner_product` | `test-effect-form-laws.R`, `test-evidence-task-ir.R`, `test-operations.R` | `R/operations.R` — the `reduce_partitions()` and `aggregate_first()` example blocks call it; replace with the reducer objects themselves |
| `measurement_space` | `test-bridge.R`, `test-effect-form-laws.R`, `test-evidence-pairing-laws.R`, `test-metric.R`, `test-print-methods.R` | `R/bridge.R` (own block, becomes internal), `R/metric.R:419-421` — the `neural_metric()` refusal example builds a cross-space bridge to show it is refused; rewrite that example to provoke the refusal without constructing a bridge |
| `measurement_bridge` | `test-bridge.R`, `test-effect-form-laws.R`, `test-evidence-pairing-laws.R`, `test-metric.R` | `R/bridge.R`, `R/metric.R:419` — same example |
| `reverse_bridge` | `test-bridge.R`, `test-evidence-pairing-laws.R` | `R/bridge.R:206-234` — the topic is absorbed into `measurement_bridge`'s internal doc |
| `effect_covariance` | `test-evidence-sampling-product.R`, `test-relation-fit.R` | `R/relation-fit.R:780` (`effect_covariance`'s own block) and `R/relation-fit.R:384` `@seealso`; `residual_df`'s block also names it |
| `residual_pair_statistics` | `helper-residual-statistics.R`, `test-metric-learning.R`, `test-print-methods.R`, `test-relation-fit.R`, `test-residual-statistics.R` | `R/residual-statistics.R:361,370` (own block), `R/crossnobis.R:258` `@seealso` |
| `bids_events` | `test-adapter-bids.R` only (lines 39, 75) | `R/adapter-bids.R:52-68` — its own `@examples` block, which becomes an internal topic's; the `@seealso` mentions at `R/adapter-bids.R:126,229` and `R/study-facts.R:426` need no edit (a `\seealso` link to a `\keyword{internal}` topic still resolves) |
| `bids_confounds` | `test-adapter-bids.R` only (lines 44, 83) | `R/adapter-bids.R:127-148` — its own block, same treatment; `@seealso` at `R/adapter-bids.R:49,230` and `R/study-facts.R:552` unchanged. `bids_study()`'s own example does **not** call either adapter (verified), so it needs no edit |

`test-first-moment-vertical-slice.R:254` names `"bids_events"` and
`"bids_confounds"` as timing-path *strings* in an `expect_setequal()`, not as
calls; it needs no change. So the whole cost of decision 1 is one test file.

Tests may use `crossform:::`; the helper `helper-residual-statistics.R` is
shared by four test files, so the switch is one edit there plus the direct
callers listed. `test-architecture.R` and `test-api.R` should gain a check that
the exported set matches this ledger, so the next export is a deliberate
decision rather than an accident of an `@export` tag.

## Open questions for review

Questions 1–4 were **closed by the maintainer decisions recorded in the next
section**; they are kept here with their resolutions so the reasoning that fed
each decision stays attached to it.

1. **`study_axis()`** (T:1, zero foreign examples on the corrected count). It is
   a typed accessor whose only alternative is reaching into `study$` directly.
   Kept for symmetry with `study_capabilities()`, but it is the thinnest
   core-ingestion export. — *Closed: keep (decision 5); the two `print` `next:`
   hints that name it are a keep-rule, not a courtesy.*
2. **`lower_effect_map()`** (T:1). `plan_relation()` calls it internally, so no
   user is obliged to. It is kept because it is the *observable* of the
   coding-invariance claim and because `raw_effect_map()`'s documentation
   contrasts against it by name. A reviewer may reasonably call it internal.
   — *Closed: keep with an obligation (decision 4); A4 owes it a runnable
   vignette chunk.*
3. **The five `pair-query.R` ER-RSA exports.** Justified entirely by a pending
   exemplar. — *Closed: keep, conditional (decision 3); the exemplar is ticket
   F2, and anything F2 does not exercise is demoted at its close.*
4. **`bids_events()` / `bids_confounds()`.** Fully composed by `bids_study()` in
   every documented workflow. Kept for the partial-BIDS case, which no test or
   vignette currently exercises. — *Closed: `demote-internal` (decision 1); a
   hypothesised user is exactly what this release subtracts.*
5. **Five capability functions, one shape.** `study_capabilities()`,
   `relation_fit_capabilities()`, `metric_capabilities()`,
   `sampling_capabilities()` are four inspectors with the same
   `capabilities(x)` signature, and `source_capabilities()` is a *constructor*
   that shares the name pattern. A single `capabilities()` S3 generic would
   remove three exports and one naming collision. Not proposed here because it
   is an API redesign rather than a tiering decision, and because collapsing the
   constructor into the generic would be wrong.
6. **`plan_crossnobis()`.** Tiered advanced and kept, but WS-B's executor
   unification proposes deleting it in favour of
   `plan_geometry(metric = shrinkage_precision())`. If WS-B lands before the
   subtraction release, this row becomes `merge-into:plan_geometry`.
7. **`relation_fit()` has zero foreign Rd examples.** It is kept as a sanctioned
   developer entry point on the strength of its role, not its usage. If the
   developer protocol document (WS-A) does not give it a worked example, the
   justification is weak.

## Maintainer decisions (2026-08-17)

The fresh-context review at the foot of this document challenged eight
dispositions and flagged three more without proposing a flip. The maintainer
ruled as follows; every table above is consistent with these rulings, and each
decision is also recorded inline at the row or note it governs.

| # | Subject | Decision | Reason (one line) |
|---:|---|---|---|
| 1 | `bids_events`, `bids_confounds` | **demote-internal** | They fail all four standing keep-rules and `bids_study()` composes both, so the demotion costs one test file and nothing else. |
| 2 | `lm_extractor` | **re-tier developer → advanced**, keep | It is an accepted `extract =` value of `relation()`, i.e. an ordinary specialist's tool rather than an extension-only seam. |
| 3 | The ER-RSA five (`match_coupling`, `control_coupling`, `coupling_contrast`, `match_control`, `pair_lm_query`) | **keep, conditional** | Their justification is the rectangular ER-RSA exemplar scheduled as ticket **F2**; any of the five that the delivered exemplar does not exercise is demoted at F2's close (an A4/F2 acceptance criterion). |
| 4 | `lower_effect_map` | **keep, with obligation** | It is the observable of the coding-invariance claim, so **A4** must add a runnable vignette chunk demonstrating it (`from-observations.Rmd` or `interpreting-results.Rmd`); undemonstrated, it reverts to demote-internal. |
| 5 | `geometry_component`, `study_axis` | **keep** | A `print` method's `next:` hint is the package's on-ramp mechanism, so a name the package tells the user to call next counts as a remedy under the standing rule. |
| 6 | `coupling`, `covariance_coupling`, `frame_conservation` | **keep** | `frame_conservation` gains real users in WS-D's conservation certificate; `coupling`/`covariance_coupling` are the connect tier's entry points and **A9** marks that tier experimental — with **A4** owing each of the three an executed example. |
| 7 | Tier definition wording | **reword, do not re-tier** | The developer tier is the extension-*only* surface (an extension package additionally meets the core-ingestion tier), and `compiler_conformance()`/`relation_plan_receipts()` stay exported because `from-observations.Rmd` executes them. |

**What the decisions change in the counts.** Decision 1 moves two adapters from
keep to demote-internal (adapters: 8 keep → 6 keep, 2 removed). Decision 2 moves
one row from developer to advanced (developer 8 → 7 rows, sanctioned set 6 → 5;
advanced 32 → 33 rows). Nothing else changes a row's tier or disposition:
decisions 3–7 confirm keeps already recorded, three of them with obligations on
tickets A4, A9 and F2 rather than on this ledger. Net: **105 exports − 7
demote-internal − 1 merge = 97 exports** after subtraction, of which the
developer tier is 5.


**F2 discharge (2026-08-20).** The rectangular ER-RSA exemplar
(`exemplars/er-rsa`) exercised all five pair-query exports as load-bearing
steps — `match_coupling()`, `control_coupling()`, `coupling_contrast()`,
`match_control()`, `pair_lm_query()` — so maintainer decision 3's condition is
met and the conditional keep on the ER-RSA five resolves to **keep**.

**Obligations this ledger now hands to other tickets.** A4: a coding-invariance
chunk for `lower_effect_map()`, executed examples for `coupling()`,
`covariance_coupling()` and `frame_conservation()`, and no presentation of the
ER-RSA five as settled API. A9: the experimental label on the connect tier in
`_pkgdown.yml` and at the head of its reference section. F2: the rectangular
ER-RSA exemplar, plus the demotion of whichever of the five it does not
exercise. WS-D: the conservation certificate that gives `frame_conservation()`
its users.

**Discharge record (2026-08-17, A4/A11).** Decision 4 is paid off: the
`coding-invariance` chunk in `vignettes/from-observations.Rmd` lowers one
`effect_map()` under a cell-means and a treatment parameterization, shows the
effect-map identity surviving while the lowering identity does not, and checks
that both extract the same amplitudes from their own design's coefficients
(0.07 s). Decision 6 is paid off by executed `\examples{}` blocks:
`man/coupling.Rd` (1.57 s), `man/coupling_views.Rd`, which is the topic
`covariance_coupling()` is documented under (1.99 s), and
`man/frame_conservation.Rd` (0.04 s) — none under `\dontrun{}`. Both
obligations are pinned by tests rather than by this paragraph:
`tests/testthat/test-from-observations-vignette.R` fails if the chunk is
removed, and `tests/testthat/test-api-surface.R` (ticket A11) holds the exported
set to the 97 names ledgered here, the 5 sanctioned developer entry points, and
the 8 demotions. Editing the export list without editing this document now
breaks the suite, which is the ratchet the subtraction release needed.

**Discharge record (2026-08-17, A9).** The two obligations this ledger handed
to A9 are paid off in `_pkgdown.yml`, with `README.md` repeating both.

- *The experimental label on the connect tier, in `_pkgdown.yml` and at the
  head of its reference section* (decision 6). The section is titled
  "Coupling and measurement (experimental, small-node only)" and its `desc`
  opens with **Experimental**, says the API and the returned kinds may change,
  and quotes the scale gate in the package's own terms: a request whose dense
  payload exceeds the version 0.1 small-node limit of 256 MiB
  (`.public_measurement_dense_limit_bytes`, `R/evidence-api.R:10-31`) is
  refused outright rather than approximated, with the remedy "use the
  support-local geometry plan; brain-scale measurement tomography is not yet
  exported". The README's tier map repeats the label and the ceiling, and the
  `evidence-pairing` entries in both the articles index and the README's
  "Where next" list carry the label too.
- *No presentation of the ER-RSA five as settled public API* (decision 3).
  They are indexed in their own section, "Advanced — matched-pair contrasts
  (provisional API)", whose `desc` states that the family is held open pending
  ticket F2 and that whichever of the five the delivered exemplar does not
  exercise becomes internal.

The reference index is now grouped by this ledger's tiers rather than by
defining file: package overview; core (the beta-first spine, then the views);
core-ingestion; advanced (queries, frame algebra, metrics and crossnobis,
extraction and error channels, numerical evidence, provisional matched-pair);
coupling (experimental); adapters; and a "For extension packages" section
holding exactly the five sanctioned developer entry points and pointing at
`vignette("crossform-extending")`. All 98 exports are indexed exactly once
across 95 topics; no demoted or `\keyword{internal}` topic appears apart from
the package overview itself. `README.md`'s first-hour path is held to the core
and core-ingestion tiers by `tests/testthat/test-readme-surface.R`, which
parses the README's fenced `r` chunks and fails on any crossform call outside
those two vectors — so an example that reaches into advanced or connect
material is a deliberate edit of that test, and of a tier row here.

## Cross-check

The ledger is verified against `NAMESPACE` by a script that parses the tier
tables out of this file and compares the name set:

```r
ns <- readLines("NAMESPACE")
exports <- sort(sub("^export\\((.*)\\)$", "\\1", grep("^export\\(", ns, value = TRUE)))
doc <- readLines("design/api-tiers.md")

# 1. Ledger rows are table rows whose first cell is `name()`.
rows   <- grep("^\\| `[a-z_0-9]+\\(\\)` \\| `R/", doc, value = TRUE)
listed <- sub("^\\| `([a-z_0-9]+)\\(\\)`.*$", "\\1", rows)
stopifnot(
  !anyDuplicated(listed),        # each export appears exactly once
  setequal(listed, exports),     # and the sets agree
  length(listed) == length(exports)
)

# 2. Every tests-only export carries an explicit disposition.
a   <- grep("^## Tests-only exports", doc)
b   <- grep("^## What demotion costs", doc)
sec <- doc[a:(b - 1L)]
tests_only_rows <- sub("^\\| `([a-z_0-9]+)`.*$", "\\1",
                       grep("^\\| `[a-z_0-9]+` \\|", sec, value = TRUE))
stopifnot(!anyDuplicated(tests_only_rows), length(tests_only_rows) == 35L)

# 3. Tier blocks partition the ledger.
tiers <- grep("^### Tier: ", doc)
per   <- vapply(seq_along(tiers), function(i) {
  hi <- if (i < length(tiers)) tiers[i + 1L] - 1L else length(doc)
  length(grep("^\\| `[a-z_0-9]+\\(\\)` \\| `R/", doc[tiers[i]:hi]))
}, integer(1))
stopifnot(sum(per) == length(exports))
```

**Check result (2026-08-17, re-run after the maintainer decisions):** all three
assertions pass.

- 105 exports in `NAMESPACE`, 105 ledger rows, 0 duplicates, 0 in `NAMESPACE`
  but missing from the ledger, 0 in the ledger but not exported.
- 35 measured tests-only exports, 35 documented, 0 duplicates, 0 missing,
  0 extra. (Membership is unchanged by the decisions: decision 1 flips two of
  these rows' *dispositions*, not their tests-only status, and decision 2 moves
  `lm_extractor` between tiers without touching this table's membership.)
- Tier blocks 21 + 20 + 33 + 16 + 8 + 7 = 105, matching each tier heading;
  dispositions 97 keep + 7 demote-internal + 1 merge-into = 105, and
  105 − 7 − 1 = **97 exports** post-subtraction.

---

## Fresh-context review (2026-08-17)

An independent reviewer re-derived the inventory from the corpora rather than
from this document, re-measured the foreign-Rd column with `tools::Rd2ex`, read
`R/adapter-*.R` and `R/neuroim2-adapter.R` to test the developer-tier claim, and
re-ran the cross-check. Scope: read-only except this file. **No disposition was
flipped** — challenges are recorded for the maintainer to decide.

### What held up

The ledger's spine is sound. An independent mechanical scan of all 105 exports
against `vignettes/*.Rmd`, `README.md`, `exemplars/**/*.R` and
`tests/testthat/{test,helper}-*.R` reproduced the Users column exactly for
**104 of 105 rows** (vignette stems, README flags, exemplar names and `T:n`
all matching), reproduced the tests-only set as exactly the same 35 names, and
confirmed every `Defining file` cell. Spot-verified and correct: the
`reducer =`/`contract =`/`metric =`/`query =` argument-value claims; the
remedy citations for `materialize_geometry`, `query_geometry`, `evaluate_geometry`
and `pair_query`; the `ingestion-contract.md` §13 promise covering
`effect_extractor`/`lm_extractor`; `plan_relation()` → `lower_effect_map()`
(`R/relation-plan.R:121`); `bids_study()` → `bids_events()`/`bids_confounds()`
(`R/adapter-bids.R:267,272`); `residual_block`'s `residual_df()` instruction
(`R/relation-fit.R:721`); `relation()`'s refusal demo (`R/relation.R:152`); the
privacy of `inner_product`'s three siblings; and "no exported consumer" for
`measurement_space`/`measurement_bridge`. All five demote-internal decisions and
the one merge are correct.

### Rows corrected

**Counting method — vignette prose read as use (11 rows).** The Users pattern
was matched against whole `.Rmd` files, so a name in backticks in narrative
prose scored identically to a call in an executed chunk. `novelty.Rmd` and
`correlation-distance-policy.Rmd` contain **no code fences at all**, so every
`V:novelty` and `V:corr-dist` credit is prose. Eleven exports have no
executable user outside `tests/` despite carrying a vignette or README credit:
`coupling`, `covariance_coupling`, `fmridesign_design_model`, `fmrireg_relation`,
`frame_conservation`, `geometry_spectrum`, `materialize_geometry`,
`metric_training_policy`, `raw_design_model`, `raw_effect_map`, `voxelwise`.
Their rows are now marked `†`. **The un-demonstrated surface is 46 exports, not
35.**

**Foreign Rd examples — 29 of 35 counts overstated.** The original measurement
counted `\seealso` `\link{}` cross-references, `\usage{}` defaults,
`\description{}` prose and `#` comments inside example blocks. Re-measured with
`tools::Rd2ex`. Largest corrections: `coupling_contrast` 4→0, `bids_events` 3→0,
`bids_confounds` 3→0, `diagonal_precision` 3→0, `match_control` 3→0,
`lower_effect_map` 3→0, `evaluate_geometry` 5→2, `residual_block` 4→1,
`geometry_component` 3→1, `effect_extractor` 2→0, `file_matrix_source` 2→0,
`pair_lm_query` 2→0, `reduce_partitions` 2→0, `effect_covariance` 2→0. This
matters because the column is used as a *cost of demotion*, and a `\seealso`
link to a `\keyword{internal}` topic costs nothing — the topic still exists.
Consequence: **22 of the 35 have zero foreign-example use**, so
`inner_product`'s claim to be the only such export is withdrawn (it was already
contradicted by the table's own `relation_fit` row).

**`T:n` (1 row).** `gaussian_covariance_model` T:2 → T:3.

**Justifications rewritten (5 rows).** `geometry_component` — named in no
remedy, and its "3 foreign examples incl. `as_neurovol`" is 1; the `as_neurovol`
mention is an `@seealso`. `identity_metric` — the dependency was stated
backwards. `study_axis` — support is *stronger* than recorded (two `print`
`next` hints). `effect_covariance` — stronger than recorded: it has no call site
in `R/` at all; `rdm_sampling_covariance()` does not call it. `numerical_agreement`
— the "§4 public list" ground is dropped as inconsistent with this ledger's own
treatment of §4.

### Dispositions challenged

The *Outcome* column records the maintainer's ruling; the full statements are in
[Maintainer decisions (2026-08-17)](#maintainer-decisions-2026-08-17) above, and
every table in this ledger reflects them.

| Name | Ledger | Proposed | Reason | Outcome |
|---|---|---|---|---|
| `bids_events` | keep | **demote-internal** | foreign-Rd support was the evidence and it is 0; fails all four standing rules; `bids_study()` calls it internally; cost is one test file | **accepted** (decision 1) |
| `bids_confounds` | keep | **demote-internal** | identical; together the clearest −2 available | **accepted** (decision 1) |
| `match_control` | keep | **demote-internal** | ER-RSA; `T:1`, foreign-Rd 0 (was 3); an alternative terminal for a contrast that already has one | deferred to F2 (decision 3) |
| `pair_lm_query` | keep | **demote-internal** | ER-RSA; `T:1`, foreign-Rd 0 (was 2); same | deferred to F2 (decision 3) |
| `match_coupling`, `control_coupling`, `coupling_contrast` | keep | **demote-internal if the exemplar does not land this release** | justified solely by an exemplar that does not exist, named in a vignette that executes nothing | rule accepted, applied at F2's close to all five (decision 3) |
| `lower_effect_map` | keep | **demote-internal, or add a vignette chunk** | satisfies no standing rule; the "observable of coding invariance" argument is good but undemonstrated. Documenting it is the better fix | chunk chosen; A4 owes it (decision 4) |
| `geometry_component` | keep | keep *if* `print` hints count as remedies | its stated remedy ground is false; the replacement ground is weaker. If rejected, `study_axis` falls with it | `print` hints **do** count; both kept (decision 5) |
| `lm_extractor` | developer | **re-tier advanced** | an `extract =` value of `relation()` and §4's *Relations* entry, not an extension-only seam | **accepted** (decision 2) |

Not challenged but newly exposed by the prose-mention correction:
`frame_conservation` (`T:2`), `covariance_coupling` (`T:2`) and `coupling` now
rest on the same evidentiary footing as the ER-RSA five and should be judged by
the same rule. — *Ruled on as decision 6: all three kept, `frame_conservation`
because WS-D's conservation certificate gives it real users and the other two
because they are the connect tier's entry points under A9's experimental label,
with A4 owing each an executed example.*

### Verdict on the developer set

**The sanctioned entries are the right seams; the number is oversold.** Checked
against what an adapter actually calls: the four in-tree adapters use
`observation_events()`, `observation_confounds()`, `study()`,
`condition_space()`, `coefficient_parameterization()`, `design_model()`,
`effect_extractor()`, `relation()`, `relation_fit()` and `abstract_domain()` —
only **2 of the 5** sanctioned entry points among them. The other three are
unexercised because no in-tree adapter brings its own out-of-memory source, not
because the seam is imaginary; the set coherently covers data-in, design-in and
error-channel-in, and no adapter is stranded. I would cut none of them.

Two caveats the release notes must carry — **both adopted**. (1) "≤ 5 sanctioned
entry points" reads as *an extension package meets 5 crossform functions*; it
meets about 14, because what adapters actually need is tiered core-ingestion.
The tier definition should say the developer tier is the extension-*only*
surface — *done, decision 7*. (2) Applied literally, the tier definition
("extension and adapter packages, not end users") also captures
`compiler_conformance()` and `relation_plan_receipts()` — the admission test and
receipt reader for an external design compiler. They are correctly kept exported
(`from-observations.Rmd` genuinely executes them), but on the old wording the
tier's own definition swept them in — *resolved by rewording rather than
re-tiering, decision 7*. Re-tiering `lm_extractor` to advanced leaves an honest
**5** — *done, decision 2*.

### Cross-check

Re-run after the maintainer decisions were applied (scratch R session reading
`NAMESPACE` and this file): **all three assertions pass** — 105 exports / 105
ledger rows / 0 duplicates / 0 in `NAMESPACE` missing from the ledger / 0 in the
ledger not exported; 35 tests-only rows with explicit dispositions and no
duplicates; tier blocks 21 + 20 + 33 + 16 + 8 + 7 = 105, matching every tier
heading. Post-subtraction the ledger stands at 97 exports
(105 − 7 demote-internal − 1 merge).

---

## Additions after the subtraction release

Everything above audits the 105 exports the subtraction release inherited and
the 97 it left. This section is the append-only register of exports added
*after* that release, so the historical arithmetic above stays readable as the
record of what it actually decided. The binding total is the sum of the two:
**97 + 10 = 107**, which is what `tests/testthat/test-api-surface.R` pins. A
row short of it is a missing justification, which is exactly the failure this
ledger exists to force.

| Name | Defining file | Tier | Program | Justification | Disposition |
|---|---|---|---|---|---|
| `frame_family()` | `R/frame.R` | advanced | ws-d (ticket D2) | the named constructor for α-weighted conservative frame families, closing gaps G1, G2 and G11 of `design/conservative-geometry-contract.md` §11.4. Before it, the only working route was `rbind()` of member weight matrices into `additive_frame(..., "conservative")`: numerically sound, undocumented, and provenance-blind — it dropped `$index` and `$specification`, so a result's `measurement` column degraded to row positions and no row could say which scale produced it. `frame_family()` is the same operator with the per-row `family` / `scale` / `center` / `alpha` metadata WS-E's transport layer requires (§7.1), and it enforces the two preconditions the per-block law needs rather than assuming them (α sums to one; every member column-normalized on its own). Users: T:1 (`test-frame-family.R`), plus `frame_conservation()`'s per-block certificate | keep |

| `contribution()` | `R/views.R` | core | ws-d (ticket D4) | the aggregation half of the attribution instrument, closing gap G4 of `design/conservative-geometry-contract.md` §11.4. A conservative frame exists so that node values can be *added up over a territory* (§1.2); before this there was no verb that did it, so every user reinvented `tapply(view$total, region, sum)` — which silently accepts a locally normalized frame whose sum estimates nothing (§1.1), averages coherence fractions instead of recomputing them, and reports a coherent budget without saying it is frame-relative (§4). `contribution()` is that arithmetic with the three refusals attached: it groups by row (budget-exact, never splitting an overlapping node), labels coherent and configuration budgets `frame_relative`, and recomputes the group fraction under the same nonnegative-partition mask a node's fraction gets. Core rather than advanced because reading a ledger is the point of asking for a conservative frame at all, not a specialist follow-up. Users: T:1 (`test-contribution.R`) | keep |

| `coherence_spectrum()` | `R/views.R` | core | ws-d (ticket D5) | the named object `design/conservative-geometry-contract.md` §3.2 declares to be the scientific product of a conservative frame family, and the resolution of gap G3 of §11.4. The contract is normative in both directions here and neither direction is discoverable from the data: per-scale *energy* is `α_s·G_Ω` by construction (§3.1), so a panel of it plots the analyst's own weight vector and may not be presented as evidence about spatial scale, while the coherent *share* is exactly invariant to α (§3.2, as strengthened by the 2026-08-20 review) and is the only scale-resolved quantity that is a finding. A user reaching for `tapply(view$coherent, scale, sum) / tapply(view$total, scale, sum)` gets the right number and none of that: no conservative-frame refusal, no mask on a negative coherent budget, no `frame_relative` label, and no statement of which of the two columns they just computed is the one that means something. It also fixes the object's granularity — G3 says the share is a function of (location, scale), not a number, so `by_location = TRUE` returns that table and no collapse is offered. Core rather than advanced for the same reason `contribution()` is: it is the readout a multiscale conservative family exists to produce, and `searchlights(c(...), "conservative")` is already core. Users: T:1 (`test-coherence-spectrum.R`) | keep |

| `location_transport()` | `R/transport.R` | advanced (experimental: population) | ws-e (ticket E2) | the typed transport object `design/population-form-contract.md` §1 makes normative, and the first WS-E symbol to exist in `R/` at all (§10 grep-verified that none did). A population layer cannot be assembled out of a bare matrix: the four things that decide what a group number *means* — the sink column, the budget-versus-density semantics, the declared row mass, and the provenance — are exactly the four a bare matrix drops. Dropping the sink silently loses `48.1 %` of a subject's native total on the contract's own fixture (§2), and a functional transport arriving without cross-fit provenance reports `3.15×` the honest transport gain while being indistinguishable from an honest one in a result object (§7.2). The constructor is where those become refusals rather than caveats, and its content-addressed signature is what lets a transport enter a scientific plan identity — two analyses differing in `P` estimate different things (§1.5). Users: T:1 (`test-transport.R`) | keep |

| `anatomical_transport()` | `R/transport.R` | advanced (experimental: population) | ws-e (ticket E2) | the registration-or-parcellation transport, as a named rule instead of a hand-built indicator matrix. Every native node goes whole to the nearest group centre, ties to the lowest group position, and — when a `radius` is given — an out-of-range node goes entirely to the sink rather than to its nearest neighbour anyway. The rule is written down because a hard assignment has row entropy `0` by construction while a partial-volume warp does not (`design/population-form-contract.md` §7.5), so which one a result was built on is a fact a reader needs and cannot recover from the operator alone. It attaches the node centres to both index tables, which is what makes the displacement diagnostics recomputable from the transport by itself. Users: T:1 (`test-transport.R`) | keep |

| `external_transport()` | `R/transport.R` | advanced (experimental: population) | ws-e (ticket E2) | the door §9.2 requires: `crossform` refuses image registration, functional-transport learning, and image resampling, and instead "accepts a typed, sparse, provenance-bearing transport as an input". Without this wrapper that refusal has no counterpart and the boundary reads as a gap. It takes the operator as it stands, lifts node identifiers off the dimension names, and defaults `method` to `"external"` — which says only that some other program built it, and is deliberately *not* a route around the cross-fit obligation: an operator fitted to responses is `"functional"` whoever fitted it, and the refusal fires the same way. Users: T:1 (`test-transport.R`) | keep |

| `transport_values()` | `R/transport.R` | advanced (experimental: population) | ws-e (ticket E2) | the contraction the transport exists to perform, `Pᵀc`, with the two semantics and the sink rule attached. Core to the object rather than a helper: the budget-preservation law (§2) is a property of *this* verb — transported total plus sink equals native total, exactly, for the signed cross-validated ledgers a conservative geometry actually produces — and the density branch is where a group node reached by no native mass comes back `NA` instead of `0`, and where the sink stays in budget units because a density of unmapped territory has no referent (§1.3). Written by hand at the call site each of those is a coin flip. Users: T:1 (`test-transport.R`) | keep |

| `plan_population()` | `R/population-plan.R` | advanced (experimental: population) | ws-e (ticket E3) | the group-level estimand `design/population-form-contract.md` §§1–4 specify, and the object that turns four separate declarations into one sealed identity. Without it a group analysis is a script, and the three choices that decide what its number *means* — which transport carried each participant (§1.5), which of the closed set of three budget normalizations made incommensurable native totals commensurable (§4, measured to disagree by up to `94 %` and to move the argmax), and which fit and evaluation order produced it (§3.3) — live nowhere a reader can check them. Each is a plan-identity field here, so two analyses differing in any one of them have different identities rather than the same number twice. The constructor is also where the preconditions become refusals: a subject frame that is not conservative has no budget behind its node values, a learned-metric plan carries a data-dependent metric that differs per participant, and a rank-deficient group design has unidentified coefficients — the refusal names the aliased columns rather than dropping one and changing what the rest mean. The model matrix is built once with a pivoted QR stored on the plan, which is the whole numerical content of an OLS group fit whose design varies along neither the node nor the coordinate axis. `precision_weighted` stays in the closed set and is refused (§4.5, §14.1): its precision is the variance of a conserved budget, which needs the cross-node sampling covariance gap G8 records as missing. D8's query bank, which landed alongside this ticket, does not close it — it reads one measurement at a time and refuses `scope = "cross_measurement"` for want of a spatial model, which is precisely the route a conserved-budget variance would take. Users: T:1 (`test-population-plan.R`) | keep |

| `estimate_population()` | `R/population-driver.R` | advanced (experimental: population) | ws-e (ticket E4) | the authorized group-level execution verb, standing to `plan_population()` exactly as `estimate_relation()` stands to `plan_relation()`: the plan names the estimand and reads nothing, and this is the one place a group number is computed. It is the query-first path `design/population-form-contract.md` §3 forces — query, then transport, then OLS — and the order is not a performance choice: query and transport act on different tensor factors so they commute exactly, but with heterogeneous native frames there is no common node axis, so fit-then-transport is not even definable and the transport *must* precede the fit. Reading the query first is what makes the route affordable: complete geometry at a native node is `q(q+1)/2` packed coordinates and a bank of `K` contrasts is `K` numbers, so a four-condition study read through four contrasts never allocates the ten-wide packed rows the naive route transports. Three things become facts of the result rather than caveats around it: the §2 budget certificate is asserted per participant and per query at fit time against an L¹-scaled tolerance and refuses rather than renormalizing; the §4 normalization is applied as the group functional it is, with `unit_budget` refused on density semantics (which conserves nothing, so the share would be a share of nothing) and a signed divisor that marks a participant `NA` rather than emitting a divergent share; and the transported components carry the §8.1 names, so a reader is never handed a `coherent` column that invites reading as a group-node common mode. The sink is a row of the output, fitted like any other node, because covariate-linked sink mass is how a differentially failing transport becomes visible. The sampling-covariance passthrough carries D8's per-native-node `K`-by-`K` blocks untransported and attaches a refusal for the transported covariance, which needs the cross-node terms `cross_node_sampling_covariance` refuses. Users: T:1 (`test-population-executor.R`) | keep |

| `latent_geometry()` | `R/latent.R` | advanced | ws-d (ticket D7) | the named object `design/conservative-geometry-contract.md` §6 requires and the resolution of gap G7 of §11.4. §6 is normative that effective counts, cumulative-contribution curves and contribution fractions live **only** on a declared nonnegative projection of the signed estimates, and until this the package had the prohibition and not the layer: `coherence_fraction` masking was the sole instance of the discipline, and a user wanting `n_eff` or "the top two modes explain 90 %" had `pmax(geometry_spectrum(g)$values, 0)` and no record that anything had been clipped. That expression is wrong in three ways this function is right in, and none of them is discoverable from the number it returns. First, the projection is a *choice*: G7 says "nonnegativity projection" names no single operator — a per-node total clamp and an eigenvalue truncation of the form move different mass — so `method` is a closed set, `"psd_projection"` is the built member, and `"nearest_psd"` is declared and refuses rather than silently substituting. Second, the clipping is *never silent*: `$moved_mass` is the absolute mass removed per measurement and `$moved_share` is it against total **absolute** mass (the signed trace is the wrong denominator — it can be small, zero or negative while the form is far from PSD), both are on the projection receipt at `$projection` beside the operator that moved them, and the print states them. Third, `n_eff` is the participation ratio of the **projected** spectrum; on the signed one the numerator is the square of a signed sum and counts nothing. The layer also refuses what cannot be projected — a query readout has been contracted against fixed weights, so clamping it is the per-node clamp G7 distinguishes — and it leaves the signed source byte-identical, which is asserted rather than assumed. Advanced rather than core because the signed layer is what a conservative analysis reports and this is the descriptive follow-up, and because the object exists partly to make the *prohibition* legible: the fractions are here because they are not allowed anywhere else. Users: T:1 (`test-latent-geometry.R`) | keep |

Tier arithmetic after this section: core 21 → **23**; advanced 33 → **35**;
advanced (experimental) 16 → **22**, of which the six new ones are the
population sub-label; total 97 → **107**.
