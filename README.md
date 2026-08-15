# effectagram

`effectagram` asks a concrete task-fMRI question: **where does an experimental
effect or representational geometry reproduce across independent data
partitions?**

You declare the effects, the spatial measurement, and which partitions must
generalize. The package compiles that scientific request once. Contrasts,
squared-distance RDMs, and linear RSA then become views of the same fitted
geometry rather than separate analysis pipelines.

This is an experimental package. It preserves negative crossvalidated
estimates, runs sequentially, and does not preprocess fMRI, register images, or
perform group inference.

## A complete first analysis

The generated fixture contains four conditions in four independent runs over a
small volume. It retains raw trial responses and a known central multivariate
signal, so both the point estimate and its admitted uncertainty path can run.

```r
library(effectagram)

example <- example_fmri_effects()

plan <- plan_geometry(
  example$fit$relation,
  at = example$frame,
  over = cross_partitions(
    example$fit$relation,
    independence = "independent",
    generalizes_over = "run"
  )
)
```

`plan` is the estimand-bearing object. It records the named conditions, exact
neural domain, searchlight frame, partition pairs, metric, and units before
reading a neural block.

Ask where the animate-versus-inanimate contrast reproduces:

```r
effect <- contrast(plan, example$contrast)
peak <- which.max(effect$total)

c(
  peak_measurement = peak,
  total_energy = effect$total[peak],
  planted_signal = peak %in% example$truth$signal_measurements
)
```

`effect$signed` is the familiar signed contrast at each measurement — a
retained first-moment marginal, not an energy. The three energies decompose
its reproducibility: `total` is reproducible contrast energy; `coherent` is
the part carried by the locally averaged contrast (a quadratic quantity, not
the signed effect itself); `configuration` is reproducible spatial departure
from that weighted mean, orthogonal in the frame-weighted inner product; and
they sum to `total` up to floating-point tolerance. These crossvalidated
quantities can be negative and are not truncated.

Change the question without refitting the geometry:

```r
distances <- rdm(plan)
category_fit <- rsa(plan, models = list(category = example$model_rdm))

dim(distances$values)       # 245 searchlights by 6 condition pairs
dim(category_fit$coefficients)

# Or ask for exactly the pairs you mean; the rest is never materialized.
rdm(plan, pairs = rbind(c("face", "house")))
```

For conditions `i` and `j`, `rdm()` reports

\[
d_{ij} = G_{ii} + G_{jj} - 2G_{ij}.
\]

Under cross-partition pairing this is a signed crossvalidated squared
Euclidean distance; a fixed neural metric makes it squared Mahalanobis. It is
not `1 - Pearson correlation`. The deliberate boundary is documented in the
[correlation-distance policy](correlation-distance-policy.md).

## Start before the beta matrices

When raw observation-by-feature responses are available, the fitted relation
can be planned without hiding the first-level model in an analysis script:

```r
facts <- study(
  observations(response_runs, scan_indexes, native_domain),
  events(event_table),
  observation_confounds(confound_table, censor = "retained"),
  partition_hierarchy(run_table)
)

relation_request <- plan_relation(
  facts,
  model = declared_design,
  effects = named_condition_effects,
  observation_model = observation_model(
    "ols",
    sampling_unit = "scan",
    independence = "runs independently acquired conditional on the model"
  )
)
fit <- estimate(relation_request)
```

`plan_relation()` is the estimand-bearing request. A concrete design matrix,
coding, censor realization, rank, aliases, and solver live in portable design
receipts; the resulting fit binds those receipts and the exact observation
source revisions. Effects are declared against named condition coordinates,
so supported cell-means and treatment codings retain one plan identity while
their execution receipts remain distinct.

The complete executable journey—including timed events, scan-level confounds,
censoring, relation estimation, geometry, RDM, RSA, and admitted uncertainty—is
in `vignette("from-observations", package = "effectagram")`. BIDS,
`fmridesign`, and `fmrireg` are optional adapters into this typed core; they do
not define the core object model.

## Add an error bar only when it is earned

The fixture was built with `lm_relation_fit()`, so its residual channel is
still available. For the selected searchlight:

```r
distance_covariance <- rdm_sampling_covariance(
  plan,
  example$fit,
  target = "plugin",
  at = peak
)
distance_se <- sqrt(sampling_covariance(distance_covariance))
distance_se
```

This is a within-measurement covariance law under the declared equal-partition,
fixed-metric, separable plug-in model. It is not a spatial random-field model,
an automatic confidence interval, or a population analysis. `target = "null"`
is also available, but no target is chosen silently.

A relation made only from precomputed beta matrices remains valid for point
geometry. It cannot reconstruct discarded residual uncertainty, and
`effectagram` refuses to pretend that edge spread supplies it. Ask before
provoking: `sampling_capabilities(plan, example$fit)` reports whether the
analytic law is admitted and, if not, every unmet requirement with its
remedy. Refusals themselves are classed conditions — `catch_refusal(expr)`
returns the missing capability and remedies as data. The
[failure gallery](failure-gallery.md) shows six realistic errors the package
guards against. Unsupported callable interpretations return classed refusals;
other errors are prevented by the absence of a misleading API or by distinct
estimand identities.

## The small API most analyses need

Start with these functions:

| Task | Functions |
|---|---|
| Bind raw facts | `observations()`, `events()`, `observation_confounds()`, `study()` |
| Plan and estimate a relation | `condition_space()`, `effect_map()`, `design_model()`, `observation_model()`, `plan_relation()`, `estimate()` |
| Use precomputed effects or a raw `(X, T)` | `relation()`, `lm_relation_fit()`, `raw_design_model()`, `raw_effect_map()` |
| Define the spatial measurement | `volume_domain()`, `compile_frame()`, `searchlights()` or `regions()` |
| Declare generalization and compile | `cross_partitions()`, `plan_geometry()` |
| Read scientific results | `contrast()`, `rdm()`, `rsa()` |
| Calibrate an admitted fixed-metric RDM | `rdm_sampling_covariance()`, `sampling_covariance()` |

Beyond the self-form core, the same plan vocabulary reaches the other
closures of the calculus: `plan_geometry(x, at, over, right = )` compiles a
rectangular cross-axis plan read with `pair_query()`s (encoding-retrieval
similarity is the canonical use), `coupling(plan, between, by)` takes the
adjoint neural-side closure between named frame measurements, and
`cross_partitions(relation, independence = "independent",
generalizes_over = "run")` binds the declared
generalization axis into the estimand's identity. Metric learning, bridges,
source adapters, and execution controls are advanced surfaces. The generated
reference site groups all exported functions by these roles in
[`_pkgdown.yml`](_pkgdown.yml).

## What is actually novel?

Not RSA, crossnobis, or the analytic covariance formula. `effectagram` advances
a broader **evidence-pairing calculus**: self- and cross-experimental forms and
self- and cross-neural measurements are different boundary closures of one
typed second-order relation. In that architecture:

- crossvalidated contrast energy, squared-distance RDMs, fixed linear RSA, and
  ordered pair hypotheses are queries against one cross-partition form;
- coherent spatial level and configuration are additive components inherited
  by every fixed linear query;
- square RDMs are optional views, while unequal-axis forms and query-first
  execution retain structure that an RDM-first workflow may collapse; and
- generalization is an identity-bound part of the estimand, not a fold count.

The executable estimand contract is the proof mechanism: it separates plan
identity from execution receipts, keeps the point relation connected to its
error channel, and uses capabilities and refusals to prevent silent changes of
interpretation.

The strongest scientific claims are deliberately gated, and three gates have
landed: the one-plan coherent/configuration family, the large-condition
query-first benchmark (at 100 conditions over 1,080 searchlights, one hundred
selected pairs run in 0.23 s and the fused full RDM beats
materialize-then-project fourfold, with a `4.4e-16` oracle), and the executable
[failure gallery](failure-gallery.md). `rsatoolbox`
parity, a real rectangular exemplar, and an operational conservation example
remain work to be earned. The full claim ledger—established precedents,
current evidence with its machine-checkable artifacts, and promotion gates—is
in [What is novel in effectagram?](novelty.md).

## Real-data evidence

The [Haxby 2001 exemplar](exemplars/haxby2001/) runs the public effectagram and
rMVPA paths on 12 conditions and 577 VT searchlights. On the matched
crossvalidated squared-Euclidean/crossnobis estimand, effectagram agrees with an
independent reference loop to `1.33e-15` and with rMVPA to `8.88e-16`. Refitting
the raw responses to retain the error channel reproduces the point RDM to
`4.44e-16` before analytic covariance is transported to a fixed linear RSA
coefficient.

Those results demonstrate numerical parity and an integrated uncertainty path.
They do **not** demonstrate a matched-estimator speed advantage, and they say
nothing about correlation-distance RSA. The scripts, assumptions, refusals,
and timing records live with the exemplar. Map-scale runtime and storage claims
are qualified separately under [`benchmarks`](benchmarks/).

## Installation and deeper guides

The package is not yet published. Install a local checkout with:

```r
install.packages(".", repos = NULL, type = "source")
```

Then read:

- `vignette("from-observations", package = "effectagram")` for the complete
  facts-to-fit-to-geometry journey;
- `vignette("introduction", package = "effectagram")` for the continuous
  question-first workflow;
- `vignette("evidence-pairing", package = "effectagram")` for measurement
  forms and coupling views, including one bounded cross-domain contraction
  and the Parseval reconstruction law;
- [the effect-form contract](effect-form-contract.md) for normative algebra and
  execution identity;
- [the evidence-sampling contract](evidence-sampling-contract.md) for the exact
  admitted uncertainty specialization and its refusals.

For neuroim2 volumes, `neuroim2_volume_domain()` and
`neuroim2_searchlights()` preserve stable full-volume indices, and
`as_neurovol()` maps compact results back without interpolation or smoothing.
