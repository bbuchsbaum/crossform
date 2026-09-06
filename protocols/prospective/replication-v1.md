# Independent replication and external-analyst protocol

Status: READY AS PROTOCOL; no replication has occurred

Version: prospective-replication-v1

Frozen: 2026-08-21

The primary estimand and success/failure rule are exactly discovery-v1 and
comparators-v1. Replication succeeds only when an eligible route independently
meets the frozen conventional-ambiguity and component-profile interpretive-gain
rule with the same direction and scale ordering. A failure, unresolved result,
or discrepancy is retained; discovery significance cannot substitute for it.

## Route A: independent eligible dataset

The dataset has no participant, acquisition session, or derived-image overlap
with discovery; satisfies eligibility-v1 independently; receives its own
immutable version/checksum manifest; and maps conditions to frozen roles before
outcomes are read. Dataset-specific preprocessing may differ only where
eligibility requires it and is disclosed as a protocol deviation. The primary
analysis code and thresholds remain locked.

## Route B: external analyst locked bundle

An analyst outside the discovery analysis team receives a read-only bundle
containing configuration, environment lock, synthetic rehearsal, eligible
inputs, one command, expected file inventory, and hash verifier. Before running,
the analyst records identity/affiliation, prior exposure to discovery outcomes,
conflicts, platform, and bundle digest. They may ask operational questions but
cannot change estimands, exclusions, thresholds, or adjudication rules.

Optional blinding replaces condition labels and hides discovery estimates until
the external execution manifest is signed. The mapping key is held by a named custodian
not conducting the run.

## Independence and shared dependencies

Route A supplies independent data; Route B supplies independent execution and
decision handling but may use the discovery data. They are distinct evidence
routes and cannot be described interchangeably. Both necessarily share
Crossform source, the frozen JSON, comparator definitions, and ordinary R/system
dependencies. Shared code is disclosed and source-hashed; it limits
implementation independence even when analyst independence is achieved.

## Adjudication and discrepancies

Two named reviewers compare manifests before unblinding. Hash mismatch, missing
output, environment drift, or typed failure blocks adjudication. Numerical
discrepancies within declared tolerances are logged; larger discrepancies are
reproduced from both bundles without overwriting either result. Protocol
questions, external feedback, exact analyst wording where permission permits,
responses, deviations, dates, and decisions are appended to
external-feedback-log.md. The original message or immutable attachment hash is
preserved; summaries never replace it.

Discovery execution and replication write different artifact roots and
different evidence-ledger rows. This protocol becoming ready is not a completed
replication and does not license the word replicated.
