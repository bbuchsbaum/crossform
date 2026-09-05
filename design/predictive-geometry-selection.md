# Nested geometry selection

The internal `.geometry_nested()` orchestration takes separate outer training
and test plans, explicit candidate tuples and explicit inner splits. A tuple
declares a model basis, named model weights, response rank and penalty.
Different model-rank/normalization declarations are distinct candidates.
Candidate IDs are derived from these recipes; presentation order cannot break
a tie. Exact inner-gain ties prefer lower rank, then larger penalty, then the
canonical candidate ID. A small nonzero gain difference is not an exact tie.

Each inner split declares `train`, `validate` and a positive `weight`; weights
sum to one. Both subsets must contain independent products in the declared
outer-training pairing. No synthetic undeclared edge is introduced, and no
outer-test observation may appear in any inner split. All split and target
checks happen before fitting begins. Eight runs give the reference example:
six outer-training runs, two outer-test runs, and four/two training/validation
runs within each inner fold.

`scope="global"` selects one candidate using declared measurement weights
(uniform by default). `scope="per_measurement"` selects independently for
each measurement and records its candidate mapping. The full candidate set,
all inner gains, fold weights, tie rule and actual origin dependencies are
frozen before outer outcomes are read. The chosen candidates are refitted on
outer training; their frozen fit records include all inner training and
validation dependencies. Only chosen candidates receive outer scores.

The result retains each selected prediction, its score and the fixed mapping
from measurements to candidates. A data-dependent collection of selected
recipes is not disguised as a single globally parameterized `fit_geometry`
record. Validators rederive choices solely from the retained inner gains and
verify outer row alignment without fitting or accessing sources.

This initial internal orchestrator uses memory outputs and reserves candidate,
fold and retained prediction space before the existing per-fit/per-score
memory gates. It does not add a public model-selection framework. The primary
fixed-split workflow remains a model declaration followed by `fit_geometry()`
and `score_geometry()`; the orchestrator makes a fully specified nested
workflow testable while its API remains private.

The statistical court repeats the entire inner-selection/refit/outer-score
procedure under null and planted signal. Its error is measured against each
outer fit's conditional risk improvement. A procedure that selects candidates
using outer outcomes is a separately labelled leakage control, not an
alternative tuning option.

The adaptive calibration protocol is frozen in `benchmarks/predictive-geometry/selection-config.R`: null and aligned generators, nine rank/penalty tuples (ranks 0/1/2; penalties 0/0.15/0.5), three inner folds weighted 0.2/0.3/0.5, six outer-training and two test partitions, one measurement. The independent 1,000-dataset pilot uses seed 2026090414; production uses seed 2026090404 and at least 2,000 datasets per arm. The same two-percent generator-based equivalence margin and five-MCSE bias/precision rules apply, with 50% pilot variance headroom. The outer-test-selected negative control must have a null lower five-MCSE bias bound above 0.001 of the reference scale. Production counts are fixed from the separate pilot before the run. Global and per-measurement orchestration are both covered by the eight-run deterministic court.

The separate pilot fixed production counts at 2,000 null and 5,000 aligned-signal datasets, recorded in `benchmark-results/predictive-geometry-selection-pilot.rds`. These counts were set before production outcomes were read.
