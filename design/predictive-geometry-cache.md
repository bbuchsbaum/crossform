# Local reuse for predictive geometry

The internal nested-selection procedure creates separate training and
evaluation caches for each inner fold. The objects are owned by that
invocation; there is no process-global cache. Public `fit_geometry()` and
`score_geometry()` retain their ordinary execution contracts.

The candidate family is declared before execution. Within a fold, candidates
with the same model declaration reuse its signed reduced form. A pooled
kernel is reused for each named weight vector, and its shifted eigensystem
is reused for each penalty. Changing response rank reads the existing
eigensystem. Evaluation caches hold signed forms only: they never fit a
model or choose modes from evaluation outcomes.

## Identity and invalidation

The neural context binds the actual positive-edge observation support,
source revisions, extractor signatures, upstream training dependencies,
ordered effect and measurement target, component, neural metric, frame,
normalization, generalization axis, compute policy and requested row block.
A changed context clears every neural form and fitted spectral path. A
change confined to unused held-out source values leaves the training context
unchanged.

Within that context, form keys also bind the complete model declaration and
admitted executor budget. Model-rank and normalization changes use separate
preparations; there is no maximal-basis shortcut. Spectral keys additionally
bind the pooled support, weights, penalty, row identity and signed contents.
Model-only pools may survive context changes because they contain no neural
data. Undeclared model, weight or penalty recipes are rejected before neural
execution, rather than extending an unbudgeted cache.

Physical row order is part of the cache context. Scientific prediction IDs
still use the normal canonical target identity and are unchanged by choosing
cached execution. Cache counters and reservations are execution metadata,
separate from the selection identity.

## Memory admission

The declared candidate family bounds retained forms and per-penalty spectral
paths. The byte model includes matrix payloads, marginal/state overlap,
metadata and a conservative allowance for R objects. Both simultaneously
live training and evaluation caches are reserved before the existing
executor receives its remaining budget. Full model-dimension row buffers
keep the physical execution budget stable across response-rank candidates.
A budget that cannot hold this reservation refuses before source reads.

Caches currently require memory storage. Ordinary fitting, scoring and
edge cross-fitting retain their admitted block-storage routes. Large tuning
grids that exceed the local-cache budget can use uncached orchestration;
disk-backed cache eviction is a separate feature.

## Evidence

`test-predictive-cache.R` instruments actual materialization and spectral
calls, compares cached predictions, IDs, choices and scores with uncached
execution, exercises source/fold/target/model/extractor/upstream changes,
checks training/evaluation role separation, and proves refusal before reads
for undeclared recipes and insufficient memory. The matched performance
runner records the cost of reuse under its explicit workload; these tests
alone make no speedup claim.
