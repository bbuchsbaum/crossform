# Prospective real-data readiness gate

Status: normative state and evidence-promotion gate

Version: prospective-readiness-v1

The executable state machine in readiness.R returns exactly BLOCKED, READY, or
EXECUTED, with reasons. Its six pre-execution checks are eligible data, frozen
configuration, locked environment, bound provenance, identified analyst or
reviewer, and reserved artifact storage.

READY means only that execution may begin. It is still prospective_protocol and
does not imply a result, validation, confirmation, or replication. EXECUTED
additionally requires a verified execution manifest, proof that the protocol
hash predates outcome access, a deviation log, and a nonempty artifact root.

Discovery EXECUTED may promote its own ledger row to
completed_real_data_result. Replication EXECUTED may promote a separate row to
independent_replication only through replication-v1 and a different artifact
root. Discovery and replication cannot share a completion artifact.

Haxby and ds003745 remain existing illustrations: neither passed eligibility-v1
prospectively, neither has a pre-outcome frozen discovery hash, and neither can
be promoted by declaring this gate READY after the fact.

The current report is readiness-current.json. Future execution belongs to two
separate successor tracker tickets, one for prospective discovery execution and
one for independent replication execution. Protocol-ticket completion never
marks either successor executed.

- discovery execution: `bd-01M0HGAFB4MEQTFNWG4TMM55S9`;
- independent replication execution: `bd-01M0HGAPYHSFQ1SN4N1GEHHM59`.
