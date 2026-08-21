# Prospective dataset eligibility and exclusion protocol

Status: frozen eligibility screen, not executed evidence

Version: prospective-eligibility-v1

Frozen: 2026-08-21

The screen is applied from metadata, data dictionaries, licenses, and file
inventories before any desired contrast or component result is inspected.

## Admission requirements

An eligible discovery dataset must satisfy all of these requirements:

1. At least 24 analyzable participants before outcome-dependent exclusions.
   This is the smallest subject count in the completed population calibration
   courts; it is a workflow floor, not a universal power guarantee. A frozen
   design-specific power analysis may require more.
2. At least four independent evaluation partitions per participant, with at
   least two observations per condition in each partition. Transport fitting
   uses different data or leave-one-participant-out information and names every
   fitting/evaluation partition.
3. At least three prespecified conditions supporting one primary contrast and
   negative-control contrasts under one common full-rank coding.
4. A spatial domain with coordinates, stable feature identifiers, declared
   support/mask construction, row mass, and a group-space transport with an
   explicit sink. At every primary group node, at least 80 percent of planned
   subjects must contribute and mean retained territory must be at least 70 percent.
5. Participant metadata needed for the frozen population model, plus motion,
   censoring, run, coverage, transport-quality, and exclusion diagnostics.
6. Documented preprocessing and first-level modeling, including software
   versions, confounds, censoring, HRF/design choices, and deviations.
7. A license permitting analysis and derived provenance, reproducible access,
   immutable dataset/version identity, and file checksums.

## Exclusions

The dataset or affected target is excluded for circular or outcome-selected
ROIs, unavailable independent partitions, incompatible condition coding,
undocumented preprocessing, unidentified participant/run mappings, missing
license or immutable version, response-derived non-cross-fitted transport,
unrecoverable subject sets, or failure of the frozen coverage/transport floor.
Zero filling and silent condition, run, node, or participant substitution are
forbidden.

Participant exclusions must be decidable from frozen metadata and QC rules.
Node exclusions are reported as failures of the primary all-planned target;
they do not silently create an available-at-node headline estimand.

## Existing illustrations

Haxby 2001 has independent runs and accessible provenance, but Crossform's
population illustration has six participants, includes already inspected
outcomes, and was not selected under this frozen protocol. It is ineligible for
prospective promotion and remains existing illustration/external parity.

OpenNeuro ds003745 v2.1.1 is CC0, has five trust-task evaluation runs and a
separate two-run transport task, explicit checksums, and full 12-of-12 reported
group-node coverage. It nevertheless has only 12 participants, uses an already
inspected retrospective analysis and outcome-informed design decisions, and
was not admitted before outcome access. It remains an existing illustration.

## Amendments

This version is immutable. Any deviation creates a new dated file and machine
configuration, states the reason and affected degrees of freedom, preserves
this version, and is frozen before newly eligible outcomes are inspected. An
amendment cannot retroactively make Haxby or ds003745 prospective.
