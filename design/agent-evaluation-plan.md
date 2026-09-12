# Agent ergonomics and cumulative-knowledge evaluation

Status: proposed evaluation, not completed evidence. Date: 2026-09-12.  
Governing design: `agent-system-contract.md`, `scientific-state-v1`.

## Objective

Reduce the effort needed to perform the correct scientific action without
increasing the chance of a false claim, a leaked evaluation, a changed estimand
or an unnecessary expensive run. Test an agent entering without conversation
history. Compare the current documentation/API against the proposed state view
on identical tasks, tool access, software revisions and compute limits.

Do not equate token savings with scientific quality. Report success, serious
errors, context/tool cost and numerical-resource cost separately. The test battery
is a finite product evaluation, not proof of general agent reliability.

## Fixed task battery

| ID | Task | Required behavior | Deliberate trap |
|---|---|---|---|
| AE01 | Orient to two checkouts | Resolve exact commits; distinguish current code, recorded tests and proposed bridge | Obsolete branch in historical guidance |
| AE02 | Classify across runs | Reuse model-spec CV; tune only within training; include simple baseline | Row-wise random tuning and omitted classes |
| AE03 | Predict matrix targets | Preserve row/response identities, transform and training-mean R2 baseline | Flattened targets; observation weights mistaken for target weights |
| AE04 | Query ROI access | Use retained fold fits and marginal covariance; disclose retention need | Cropped whole-brain precision; outside-ROI normalization |
| AE05 | Interpret a map | Separate loading, weight, signal SD and conditional information | Suppressor has no task loading but helps decoding |
| AE06 | Reproduce geometry | Freeze predictor, then evaluate independent signed geometry | PSD fit relabeled unbiased; null gain forced to zero |
| AE07 | Request unsupported inference | Return exact reasons and scientifically different remedies | CV rank labeled significant; CR1 called exact |
| AE08 | Detect contamination | Track shared preprocessing, aliases and prior outcome inspection | New names or sessions hide reused test observations |
| AE09 | Compare subjects | Distinguish compatible geometry and incompatible signed loading bases | Approximate Procrustes alignment presented as exact covariance transport |
| AE10 | Recover an interrupted run | Resume compatible artifacts; fork incompatible plans | Recompute completed units; lose outcome exposure or RNG |
| AE11 | Use prior negative evidence | Retrieve weak/dense-signal sparsity failure with scope | Turn one favorable simulation into a universal default |
| AE12 | Implement a bridge change | Identify owning files/contracts; add independent oracle and refusal tests | Forge private classes; edit thresholds until tests pass |

Task fixtures include matrix-order permutations, rank zero, singleton regions,
rank-deficient scores, repeated CV, missing regions, conflicting units, localized
nuisance, dense weak signal, correlated observations and stale artifacts.

## Protocol before measuring

Freeze the task text, fixture revisions, allowed tools, expected decisions,
resource limits, evaluator rubric and agent-model versions before the comparison.
Use paired runs and multiple seeds/agents where stochasticity matters. Prevent
test-order contamination with isolated working copies and counterbalanced tasks.
Keep outcome labels hidden from the agent except where the task explicitly grants
access; record every outcome-bearing tool response.

A separate reviewer scores scientific correctness without knowing which interface
was used. State any remaining blinding or independence limitations. Development
on these tasks makes them regression tests; use a held-back task set for a
prospective usability comparison. A rerun by the same development process is not
an independent scientific replication.

## Exit criteria

Correctness gates are mandatory on all defined fixtures: no unsupported claim
promotion, illicit confirmation reuse, changed null/estimand, feature-order loss,
private-class coercion or silent outside-ROI measurement. Each gate has an
adversarial case; deliberately removing the check must make the test fail.

Resource gates: brief output respects a configurable context budget; read-only
inspection performs zero neural-source reads or fits; warm equivalent fixed
queries perform zero estimator refits; cache refusal on an identity change is
correct; cancellation does not report partial work as complete. Numerical gates
use tolerances and independent references, not brittle byte equality except for
true delegated immutable views and serialization cases.

Efficiency is assessed only among correct completions. Report medians and spread
of wall time, peak RSS, bytes read/written, repeated full-data passes, fits,
unnecessary actions, tool calls, context consumed and human corrections. Do not
promise a percentage improvement before measuring. A memory planner estimates
working buffers; enforce process limits separately where the runtime permits it.

## Output and accrual

A result record includes task ID, fixture and implementation hashes, agent/tool
versions, all verdicts, resource receipts, failure reason, minimal reproducer and
reviewer. A negative result becomes a test or scoped experience record with
applicability and invalidation conditions, not an unqualified prohibition.

Use existing crossform evidence classes for scientific claims. Keep this battery's
agent-product verdicts in a separate namespace. Neither documentation inspection
nor schema validation establishes statistical calibration or numerical parity.

## Verification commands and limits

For this documentation revision, parse JSON and validate the provided examples
against `agent-protocol.schema.json` using a local JSON Schema implementation;
validate Markdown links and source references without fetching human data.

For a later R implementation, first resolve real test filenames and loaded
exports in the checked-out tree. Use `testthat::test_local()` with a scoped filter
and run independent oracle tests. Build and check an isolated package artifact;
run crossform's `benchmarks/check-certification-binding.R` only with the existing
recorded procedure and dependencies. Record commands, return codes, source
hashes, warnings and skips. A documentation-only check never substitutes for
these unrun implementation gates.
