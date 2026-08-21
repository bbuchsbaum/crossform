# Package index

## Package overview

What crossform is for, and the boundaries it will not cross.

- [`crossform-package`](https://bbuchsbaum.github.io/crossform/reference/crossform-package.md)
  [`crossform`](https://bbuchsbaum.github.io/crossform/reference/crossform-package.md)
  : crossform: cross-generalized effect geometry for task fMRI

## Core — the beta-first spine

Start here. Condition effects, the measurements you want answers at, and
the partitions that must generalize compile once into a geometry plan.
Nothing below this heading is optional for an ordinary analysis, and
nothing above the advanced sections is needed to complete one.

- [`example_fmri_effects()`](https://bbuchsbaum.github.io/crossform/reference/example_fmri_effects.md)
  : Generate a small task-fMRI relation with known spatial truth
- [`effect_space()`](https://bbuchsbaum.github.io/crossform/reference/effect_space.md)
  : Define an experimental coordinate space
- [`relation()`](https://bbuchsbaum.github.io/crossform/reference/relation.md)
  : Construct a lazy experimental-neural relation
- [`lm_relation_fit()`](https://bbuchsbaum.github.io/crossform/reference/lm_relation_fit.md)
  : Fit a partitioned linear-model relation with a residual error
  channel
- [`abstract_domain()`](https://bbuchsbaum.github.io/crossform/reference/abstract_domain.md)
  : Construct an abstract neural feature domain
- [`volume_domain()`](https://bbuchsbaum.github.io/crossform/reference/volume_domain.md)
  : Construct a native volumetric neural feature domain
- [`searchlights()`](https://bbuchsbaum.github.io/crossform/reference/searchlights.md)
  : Specify Euclidean searchlights
- [`regions()`](https://bbuchsbaum.github.io/crossform/reference/regions.md)
  : Specify region measurements
- [`voxelwise()`](https://bbuchsbaum.github.io/crossform/reference/voxelwise.md)
  : Specify point measurements
- [`whole_brain()`](https://bbuchsbaum.github.io/crossform/reference/whole_brain.md)
  : Specify a whole-brain additive measurement
- [`compile_frame()`](https://bbuchsbaum.github.io/crossform/reference/compile_frame.md)
  : Compile a spatial specification against a neural domain
- [`cross_partitions()`](https://bbuchsbaum.github.io/crossform/reference/cross_partitions.md)
  : Pair every distinct partition once
- [`plan_geometry()`](https://bbuchsbaum.github.io/crossform/reference/plan_geometry.md)
  : Compile a reusable query-first geometry plan

## Core — views, uncertainty, and refusals

The plan is an identity, not a matrix. These read it: signed and
decomposed contrast energy, squared crossvalidated distances, fixed
linear RSA, the ledger reading that adds a conservative map up over a
territory, the coherent share of a multiscale family against scale, the
analytic sampling law, and the refusal a view returns when the evidence
for it was never earned.

- [`contrast_energy()`](https://bbuchsbaum.github.io/crossform/reference/contrast_energy.md)
  : Read a contrast as signed, coherent, configuration, and total
  evidence
- [`contribution()`](https://bbuchsbaum.github.io/crossform/reference/contribution.md)
  : Add a conservative attribution map up over a territory
- [`coherence_spectrum()`](https://bbuchsbaum.github.io/crossform/reference/coherence_spectrum.md)
  : Read the coherent share of a conservative frame family against scale
- [`rdm()`](https://bbuchsbaum.github.io/crossform/reference/rdm.md) :
  Read squared experimental distances from geometry
- [`rsa()`](https://bbuchsbaum.github.io/crossform/reference/rsa.md) :
  Fit multiple-regression RSA as one compiled geometry query
- [`plot(`*`<effect_contrast_view>`*`)`](https://bbuchsbaum.github.io/crossform/reference/plot_views.md)
  [`plot(`*`<effect_rdm_view>`*`)`](https://bbuchsbaum.github.io/crossform/reference/plot_views.md)
  [`plot(`*`<effect_rsa_view>`*`)`](https://bbuchsbaum.github.io/crossform/reference/plot_views.md)
  [`plot(`*`<effect_crossnobis_view>`*`)`](https://bbuchsbaum.github.io/crossform/reference/plot_views.md)
  [`plot(`*`<effect_sampling_covariance>`*`)`](https://bbuchsbaum.github.io/crossform/reference/plot_views.md)
  : Plot geometry views
- [`sampling_capabilities()`](https://bbuchsbaum.github.io/crossform/reference/sampling_capabilities.md)
  : Ask whether the analytic sampling law is available before provoking
  it
- [`rdm_sampling_covariance()`](https://bbuchsbaum.github.io/crossform/reference/rdm_sampling_covariance.md)
  : Construct exact analytic sampling covariance for crossvalidated
  distances
- [`sampling_covariance()`](https://bbuchsbaum.github.io/crossform/reference/sampling_covariance.md)
  : Query an exact factorized sampling-covariance form
- [`catch_refusal()`](https://bbuchsbaum.github.io/crossform/reference/catch_refusal.md)
  : Inspect a capability refusal
- [`crossform_conditions`](https://bbuchsbaum.github.io/crossform/reference/crossform_conditions.md)
  : Conditions raised by crossform
- [`compute_policy()`](https://bbuchsbaum.github.io/crossform/reference/compute_policy.md)
  : Construct the version 0.1 compute policy

## Core (ingestion) — from observations to a fitted relation

The typed-facts spine, for arriving with scan responses and an event
table rather than fitted betas: study, observations, design, relation
plan, receipts. A user who starts from betas meets none of it.
[`vignette("from-observations")`](https://bbuchsbaum.github.io/crossform/articles/from-observations.md)
is its single narrative.

- [`observations()`](https://bbuchsbaum.github.io/crossform/reference/observations.md)
  : Bind raw observation sources to indexes and a neural domain
- [`observation_index()`](https://bbuchsbaum.github.io/crossform/reference/observation_index.md)
  : Define one partition's observation axis
- [`observation_events()`](https://bbuchsbaum.github.io/crossform/reference/observation_events.md)
  : Declare a typed event record
- [`observation_confounds()`](https://bbuchsbaum.github.io/crossform/reference/observation_confounds.md)
  : Declare observation-level confounds and censor facts
- [`partition_hierarchy()`](https://bbuchsbaum.github.io/crossform/reference/partition_hierarchy.md)
  : Declare nested partition axes
- [`study()`](https://bbuchsbaum.github.io/crossform/reference/study.md)
  : Bind observations, events, confounds, clocks, and partition axes
- [`study_axis()`](https://bbuchsbaum.github.io/crossform/reference/study_axis.md)
  : Select a typed study axis
- [`study_capabilities()`](https://bbuchsbaum.github.io/crossform/reference/study_capabilities.md)
  : Inspect factual study capabilities
- [`condition_space()`](https://bbuchsbaum.github.io/crossform/reference/condition_space.md)
  : Define a semantic condition space
- [`effect_map()`](https://bbuchsbaum.github.io/crossform/reference/effect_map.md)
  : Declare effects in a semantic condition vocabulary
- [`coefficient_parameterization()`](https://bbuchsbaum.github.io/crossform/reference/coefficient_parameterization.md)
  : Identify one concrete coefficient parameterization
- [`lower_effect_map()`](https://bbuchsbaum.github.io/crossform/reference/lower_effect_map.md)
  : Lower semantic effects into a concrete coefficient basis
- [`raw_effect_map()`](https://bbuchsbaum.github.io/crossform/reference/raw_effect_map.md)
  : Bind a raw coefficient-space target
- [`design_model()`](https://bbuchsbaum.github.io/crossform/reference/design_model.md)
  : Declare a semantic design model with compiled routes
- [`raw_design_model()`](https://bbuchsbaum.github.io/crossform/reference/raw_design_model.md)
  : Bind already compiled raw design matrices
- [`observation_model()`](https://bbuchsbaum.github.io/crossform/reference/observation_model.md)
  : Declare the first-moment observation model
- [`plan_relation()`](https://bbuchsbaum.github.io/crossform/reference/plan_relation.md)
  : Plan a first-moment experimental–neural relation
- [`compiler_conformance()`](https://bbuchsbaum.github.io/crossform/reference/compiler_conformance.md)
  : Inspect first-moment compiler conformance
- [`relation_plan_receipts()`](https://bbuchsbaum.github.io/crossform/reference/relation_plan_receipts.md)
  : Inspect portable design receipts
- [`estimate_relation()`](https://bbuchsbaum.github.io/crossform/reference/estimate_relation.md)
  : Estimate a planned relation family

## Advanced — explicit queries and geometries

Reaching past the view functions: materialize a geometry once and query
it many times, or evaluate a query without materializing anything at
all; custom bilinear operators, unequal axes, and explicit partition
pairings and reductions.

- [`materialize_geometry()`](https://bbuchsbaum.github.io/crossform/reference/materialize_geometry.md)
  : Materialize a complete cross-generalized geometry
- [`evaluate_geometry()`](https://bbuchsbaum.github.io/crossform/reference/evaluate_geometry.md)
  : Evaluate a fixed query without materializing complete geometry
- [`query_geometry()`](https://bbuchsbaum.github.io/crossform/reference/query_geometry.md)
  : Apply a linear query to a complete geometry
- [`geometry_component()`](https://bbuchsbaum.github.io/crossform/reference/geometry_component.md)
  : Read one component of a complete geometry
- [`geometry_spectrum()`](https://bbuchsbaum.github.io/crossform/reference/geometry_spectrum.md)
  : Read the signed eigenvalue spectrum of cross-generalized geometry
- [`bilinear_query()`](https://bbuchsbaum.github.io/crossform/reference/bilinear_query.md)
  : Describe a bilinear geometry query
- [`pair_query()`](https://bbuchsbaum.github.io/crossform/reference/pair_query.md)
  : Describe an axis-bound pair query
- [`pairing()`](https://bbuchsbaum.github.io/crossform/reference/pairing.md)
  : Construct a partition pairing
- [`aggregate_first()`](https://bbuchsbaum.github.io/crossform/reference/aggregate_first.md)
  : Aggregate edge sufficient statistics before normalization
- [`reduce_partitions()`](https://bbuchsbaum.github.io/crossform/reference/reduce_partitions.md)
  : Reduce normalized and transformed partition edges

## Advanced — frame algebra and conservation

Diagonal and α-weighted frame families, and the certificate that a frame
conserves the mass it claims to.

- [`additive_frame()`](https://bbuchsbaum.github.io/crossform/reference/additive_frame.md)
  : Describe an additive diagonal spatial frame
- [`frame_family()`](https://bbuchsbaum.github.io/crossform/reference/frame_family.md)
  : Combine conservative frames into one alpha-weighted family
- [`frame_conservation()`](https://bbuchsbaum.github.io/crossform/reference/frame_conservation.md)
  : Diagnose local-to-global conservation of a compiled frame

## Advanced — the latent PSD layer

Crossvalidated contributions are signed, and effective counts,
cumulative curves and contribution fractions are invalid on them: the
implied shares can leave `[0, 1]`, and clipping reintroduces the bias
the cross-partition estimator removes. Those functionals are legal only
on a declared nonnegative projection, which this builds and labels — a
named projection from a closed set, with the mass it moved reported per
measurement rather than absorbed.

- [`latent_geometry()`](https://bbuchsbaum.github.io/crossform/reference/latent_geometry.md)
  : Construct the latent PSD descriptive layer of a signed geometry

## Advanced — metrics and crossnobis

Fixed and learned neural metrics, the leakage guard on metric training,
and the separately qualified crossnobis path.

- [`neural_metric()`](https://bbuchsbaum.github.io/crossform/reference/neural_metric.md)
  : Construct a support-local same-space neural metric
- [`noise_precision()`](https://bbuchsbaum.github.io/crossform/reference/noise_precision.md)
  : Construct a fixed neural noise-precision metric
- [`identity_metric()`](https://bbuchsbaum.github.io/crossform/reference/identity_metric.md)
  : Specify an on-demand identity metric
- [`diagonal_precision()`](https://bbuchsbaum.github.io/crossform/reference/diagonal_precision.md)
  : Specify on-demand diagonal residual-variance precision
- [`shrinkage_precision()`](https://bbuchsbaum.github.io/crossform/reference/shrinkage_precision.md)
  : Specify on-demand shrinkage-to-diagonal precision
- [`metric_training_policy()`](https://bbuchsbaum.github.io/crossform/reference/metric_training_policy.md)
  : Declare which residual partitions may train a metric
- [`metric_capabilities()`](https://bbuchsbaum.github.io/crossform/reference/metric_capabilities.md)
  : Inspect exact neural-metric capabilities
- [`plan_crossnobis()`](https://bbuchsbaum.github.io/crossform/reference/plan_crossnobis.md)
  : Compile an on-demand learned-metric crossnobis plan
- [`crossnobis()`](https://bbuchsbaum.github.io/crossform/reference/crossnobis.md)
  : Evaluate a signed local crossnobis contrast

## Advanced — extraction maps, error channels, and residuals

Where the estimates and their error channel come from: the extraction
map in `B = E Y`, whether a fit carries residuals at all, and the
bounded residual read with its correct divisor.

- [`lm_extractor()`](https://bbuchsbaum.github.io/crossform/reference/lm_extractor.md)
  : Compile a supplied linear model into an effect extractor
- [`relation_fit_capabilities()`](https://bbuchsbaum.github.io/crossform/reference/relation_fit_capabilities.md)
  : Inspect statistical capabilities of a relation or relation fit
- [`residual_block()`](https://bbuchsbaum.github.io/crossform/reference/residual_block.md)
  : Read fitted residuals for a neural feature block
- [`residual_df()`](https://bbuchsbaum.github.io/crossform/reference/residual_df.md)
  : Read residual degrees of freedom from a fitted partition

## Advanced — numerical evidence

The reproducibility contract as a value, and the verifier that two runs
agree under a named guarantee.

- [`numerical_contract()`](https://bbuchsbaum.github.io/crossform/reference/numerical_contract.md)
  : Define crossform's numerical reproducibility contract
- [`numerical_agreement()`](https://bbuchsbaum.github.io/crossform/reference/numerical_agreement.md)
  : Assess results under a named numerical guarantee

## Advanced — matched-pair contrasts (provisional API)

Matched encoding–retrieval pairs, eligible controls, and the balanced
matched-minus-control operator. **Provisional, not settled public API:**
this family is held open pending the rectangular ER-RSA exemplar, and
whichever of these the delivered exemplar does not exercise becomes
internal. Do not build against it yet.

- [`match_coupling()`](https://bbuchsbaum.github.io/crossform/reference/match_coupling.md)
  : Mark matched encoding-retrieval effect pairs
- [`control_coupling()`](https://bbuchsbaum.github.io/crossform/reference/control_coupling.md)
  : Mark eligible control pairs
- [`coupling_contrast()`](https://bbuchsbaum.github.io/crossform/reference/coupling_contrast.md)
  : Contrast matched and control pair couplings
- [`match_control()`](https://bbuchsbaum.github.io/crossform/reference/match_control.md)
  : Compile a matched-versus-control pair-space coefficient
- [`pair_lm_query()`](https://bbuchsbaum.github.io/crossform/reference/pair_lm_query.md)
  : Compile a weighted pair-space linear-model coefficient to an effect
  query

## Coupling and measurement (experimental, small-node only)

**Experimental.** Measurement forms, coupling closures, and tomography —
results that relate two sets of measurements to each other rather than
scoring one. The API and the returned kinds may change. This tier is
**scale-gated and the package refuses brain-scale claims here**: a
request whose dense payload exceeds the version 0.1 small-node limit of
256 MiB (268,435,456 bytes) is refused outright rather than
approximated, with the remedy “use the support-local geometry plan;
brain-scale measurement tomography is not yet exported.”
[`vignette("evidence-pairing")`](https://bbuchsbaum.github.io/crossform/articles/evidence-pairing.md)
is the narrative.

- [`measurement_frame()`](https://bbuchsbaum.github.io/crossform/reference/measurement_frame.md)
  : Declare identified neural measurements
- [`edge_frame()`](https://bbuchsbaum.github.io/crossform/reference/edge_frame.md)
  : Declare requested neural measurement edges
- [`measurement_form()`](https://bbuchsbaum.github.io/crossform/reference/measurement_form.md)
  : Materialize requested measurement-space forms
- [`measurement_components()`](https://bbuchsbaum.github.io/crossform/reference/measurement_components.md)
  : Summarize a crossed node decomposition on one measurement edge
- [`variation_query()`](https://bbuchsbaum.github.io/crossform/reference/variation_query.md)
  : Declare a repeated-variation experimental query
- [`coupling()`](https://bbuchsbaum.github.io/crossform/reference/coupling.md)
  : Take the adjoint-side coupling closure of a geometry plan
- [`effect_coupling()`](https://bbuchsbaum.github.io/crossform/reference/coupling_views.md)
  [`covariance_coupling()`](https://bbuchsbaum.github.io/crossform/reference/coupling_views.md)
  [`canonical_coupling()`](https://bbuchsbaum.github.io/crossform/reference/coupling_views.md)
  [`geometry_alignment()`](https://bbuchsbaum.github.io/crossform/reference/coupling_views.md)
  : Interpret completed measurement forms
- [`effect_coupling_result`](https://bbuchsbaum.github.io/crossform/reference/effect_coupling_result.md)
  : Coupling results and the readings they may carry
- [`gaussian_covariance_model()`](https://bbuchsbaum.github.io/crossform/reference/gaussian_covariance_model.md)
  : Declare a joint Gaussian covariance interpretation
- [`connectivity()`](https://bbuchsbaum.github.io/crossform/reference/connectivity.md)
  : Request a validated connectivity view
- [`reconstruct_evidence()`](https://bbuchsbaum.github.io/crossform/reference/reconstruct_evidence.md)
  : Reconstruct or project the global neural evidence operator

## Population (experimental)

**Experimental.** Carrying one participant’s conservative geometry onto
a shared group node set. A transport is an *input*: crossform accepts a
typed, sparse, provenance-bearing operator and refuses to learn one —
image registration and functional-transport learning stay outside the
package. The sink column is required and always materialized, so partial
coverage is a number you can read rather than budget that quietly went
missing, and `semantics` has no default because budget and density are
two different estimands. A transport built from the response data must
name the runs that built it.
[`plan_population()`](https://bbuchsbaum.github.io/crossform/reference/plan_population.md)
assembles the group estimand from them: the transport, the budget
normalization, and the fit’s evaluation order all enter its scientific
identity, because each of them changes what the group number is rather
than how it is computed.
[`estimate_population()`](https://bbuchsbaum.github.io/crossform/reference/estimate_population.md)
executes one through a bank of contrasts, reading the query before the
transport so complete packed geometry is never allocated, and asserting
each participant’s budget certificate at fit time.
[`materialize_population()`](https://bbuchsbaum.github.io/crossform/reference/materialize_population.md)
is its complete-form counterpart, for when the whole group form is
wanted rather than a bank: it streams the packed coordinate axis in
tiles, so the group stack stays linear in the tile and the dense
participant-by-node-by-coordinate array is never built.
`population_views` are the reader verbs over either result —
[`contrast_energy()`](https://bbuchsbaum.github.io/crossform/reference/contrast_energy.md),
[`rdm()`](https://bbuchsbaum.github.io/crossform/reference/rdm.md),
[`rsa()`](https://bbuchsbaum.github.io/crossform/reference/rsa.md) and
[`contribution()`](https://bbuchsbaum.github.io/crossform/reference/contribution.md)
again, reading the group coefficients that run already produced rather
than executing a second one; a query the estimated basis cannot reach is
refused with the re-estimation remedy instead of being projected onto
it.
[`population_uncertainty()`](https://bbuchsbaum.github.io/crossform/reference/population_uncertainty.md)
reports the two error bars section 7 requires to be kept apart: the
between-subject scatter of the participants about the group fit, and —
only where a group column is fed by exactly one native row, so the
cross-node terms carry weight zero — the transported within-subject
sampling variance. They are separate blocks with no field holding their
sum, and the `t` is labelled uncalibrated whatever the recorded null
simulation measured, because a real fit adds misspecification and
transport heterogeneity that simulation did not contain.
[`heterogeneity()`](https://bbuchsbaum.github.io/crossform/reference/heterogeneity.md)
reads the other half of section 6’s decomposition: the N-by-N subject
Gram of participants’ deviations from the group fit, its spectrum, the
subject loadings on its modes and — at nodes you name — the geometry of
those modes. It defaults to the cross-fitted estimator, because the
plug-in one books every participant’s within-subject sampling noise as
between-participant heterogeneity and inflates the trace by a measured
62.7%; the plug-in route stays available for loading directions, named
as such and carrying no effective mode count.
[`population_prevalence()`](https://bbuchsbaum.github.io/crossform/reference/population_prevalence.md)
is the descriptive counterpart: how many participants stand behind a
group value, at each group node and query and across the bank as a
whole. It lives on the latent layer section 6.5 confines fractions to,
prints the same “not for inference” sentence
[`latent_geometry()`](https://bbuchsbaum.github.io/crossform/reference/latent_geometry.md)
does, and refuses to carry a field that looks like an error bar — a cell
at which nothing reproduces reports a fraction near 0.5, not near 0. See
[`design/population-form-contract.md`](https://github.com/bbuchsbaum/crossform/blob/main/design/population-form-contract.md).

- [`location_transport()`](https://bbuchsbaum.github.io/crossform/reference/location_transport.md)
  : Declare a location transport onto a shared group node set
- [`anatomical_transport()`](https://bbuchsbaum.github.io/crossform/reference/anatomical_transport.md)
  : Build a hard nearest-centre transport from node coordinates
- [`external_transport()`](https://bbuchsbaum.github.io/crossform/reference/external_transport.md)
  : Admit a transport built outside crossform
- [`transport_values()`](https://bbuchsbaum.github.io/crossform/reference/transport_values.md)
  : Carry native node values onto the group nodes
- [`plan_population()`](https://bbuchsbaum.github.io/crossform/reference/plan_population.md)
  : Plan a population form over transported conservative geometry
- [`estimate_population()`](https://bbuchsbaum.github.io/crossform/reference/estimate_population.md)
  : Estimate a planned population form through a bank of queries
- [`materialize_population()`](https://bbuchsbaum.github.io/crossform/reference/materialize_population.md)
  : Materialize a planned population form at every group node
- [`contrast_energy(`*`<effect_population_result>`*`)`](https://bbuchsbaum.github.io/crossform/reference/population_views.md)
  [`rdm(`*`<effect_population_result>`*`)`](https://bbuchsbaum.github.io/crossform/reference/population_views.md)
  [`rsa(`*`<effect_population_result>`*`)`](https://bbuchsbaum.github.io/crossform/reference/population_views.md)
  [`contribution(`*`<effect_population_result>`*`)`](https://bbuchsbaum.github.io/crossform/reference/population_views.md)
  [`contribution(`*`<effect_population_view>`*`)`](https://bbuchsbaum.github.io/crossform/reference/population_views.md)
  : Reader verbs on an estimated population form
- [`heterogeneity()`](https://bbuchsbaum.github.io/crossform/reference/heterogeneity.md)
  : Decompose a population form into consensus and heterogeneity
- [`population_uncertainty()`](https://bbuchsbaum.github.io/crossform/reference/population_uncertainty.md)
  : Read the group uncertainty layers of an estimated population form
- [`population_prevalence()`](https://bbuchsbaum.github.io/crossform/reference/population_prevalence.md)
  : Count how many participants a transported group value stands on

## Adapters

Morphisms into the typed core from another package’s vocabulary. Every
dependency here is `Suggests:`, and the core is not shaped by any of
them. The three `neuroim2` entries are de-facto core for volumetric work
— every real imaging path begins and ends there.

- [`neuroim2_volume_domain()`](https://bbuchsbaum.github.io/crossform/reference/neuroim2_volume_domain.md)
  : Construct a crossform domain from a neuroim2 volume mask
- [`neuroim2_searchlights()`](https://bbuchsbaum.github.io/crossform/reference/neuroim2_searchlights.md)
  : Compile neuroim2 searchlight indices into a crossform frame
- [`as_neurovol()`](https://bbuchsbaum.github.io/crossform/reference/as_neurovol.md)
  : Map a compact result vector back to a neuroim2 volume
- [`bids_study()`](https://bbuchsbaum.github.io/crossform/reference/bids_study.md)
  : Bind BIDS files into a generic study
- [`fmridesign_design_model()`](https://bbuchsbaum.github.io/crossform/reference/fmridesign_design_model.md)
  : Compile a semantic design model from fmridesign
- [`fmrireg_relation()`](https://bbuchsbaum.github.io/crossform/reference/fmrireg_relation.md)
  : Execute a relation plan with fmrireg

## For extension packages

The five sanctioned entry points an extension or adapter *package* uses
to hand data, a design, and an error channel into the core, plus the
version certificate every adapter is obliged to obtain before it touches
the package it adapts. Not an end user’s surface, and deliberately
small: everything else that once lived here is internal.
[`vignette("crossform-extending")`](https://bbuchsbaum.github.io/crossform/articles/crossform-extending.md)
states what is open, what is closed, and what an extension may rely on
remaining true.

- [`file_matrix_source()`](https://bbuchsbaum.github.io/crossform/reference/file_matrix_source.md)
  : Describe a read-only column-major matrix file
- [`source_capabilities()`](https://bbuchsbaum.github.io/crossform/reference/source_capabilities.md)
  : Declare source execution capabilities
- [`effect_extractor()`](https://bbuchsbaum.github.io/crossform/reference/effect_extractor.md)
  : Construct an explicit linear effect extractor
- [`relation_fit()`](https://bbuchsbaum.github.io/crossform/reference/relation_fit.md)
  : Attach statistical error information to a relation
- [`relation_block()`](https://bbuchsbaum.github.io/crossform/reference/relation_block.md)
  : Read one experimental-neural relation block
- [`adapter_version_certificate()`](https://bbuchsbaum.github.io/crossform/reference/adapter_version_certificate.md)
  : Certify the installed version of an adapter's upstream package
