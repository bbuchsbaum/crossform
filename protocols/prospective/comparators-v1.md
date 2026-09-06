# Prospective ambiguity and comparator protocol

Status: frozen decision rule, no empirical result

Version: prospective-comparators-v1

Frozen: 2026-08-21

All comparators use the identical admitted participants, runs, condition coding,
regions, features, cross-validation partitions, transport, coverage policy, and
normalization as the Crossform component analysis. A mismatch is a protocol
failure, not a sensitivity result.

## Comparator estimands

1. Activation: the population coefficient of the signed regional mean contrast,
   in the original response units.
2. Aggregate MVPA: the all-unordered-partition-pairs crossnobis coefficient,
   summed with the declared conservative frame weights and reported in squared
   response units.
3. Fixed linear RSA: the prespecified bilinear model contrast on the same
   crossvalidated condition geometry, with no correlation normalization,
   adaptive model selection, or outcome-fitted weights.
4. Crossform: directly fitted total, coherent, and configuration population
   coefficients at each frozen scale. A share is derived from population
   coefficients only where the total denominator and subject-set identity are
   valid; participant ratios are not averaged.

## Quantitative decision rule

Conventional ambiguity is established only if both activation and aggregate
MVPA differences meet their two-sided equivalence tests: the entire 95 percent pointwise HC3
interval lies within plus or minus 15 percent of the prespecified
reference magnitude. Fixed linear RSA must either meet its own frozen 15
percent equivalence margin or be reported separately as resolving ambiguity.

Interpretive gain requires all of the following: conventional ambiguity;
coherent-share profile separation of at least 0.20 at two adjacent nonpoint scales;
matching total effects within the 15 percent equivalence margin; the
same component ordering under HC3 and the frozen wild-bootstrap sensitivity;
and no coverage, subject-set, sink, transport-quality, or influential-subject
warning that crosses its frozen threshold. Statistical significance alone,
a smaller p-value, or a narrower interval is not new organization information.

Crossform fails the primary criterion if any required conventional comparator
is not ambiguous, the total is not matched, profile separation is below 0.20,
ordering changes across supported uncertainty methods, a required cell is
unresolved, or support diagnostics fail. Failure is retained and reported; the
criterion is never redefined after inspection.

## Negative controls and alternatives

The activation comparator must detect the frozen synthetic mean-shift control.
Condition-label and organization-label permutations must not create profile
separation. Alternative explanations reported beside the primary decision are
SNR, spatial smoothness, unequal coverage, sink loss, transport quality,
motion/censoring, leverage, one-subject influence, and model/ROI misspecification.
These are diagnostic explanations, not automatic causal corrections.
