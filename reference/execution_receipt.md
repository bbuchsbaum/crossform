# Record the execution that produced a geometry or view

A receipt separates scientific-plan identity, numerical execution
policy, source revisions, and reduction identity. Progress reporters and
other nonsemantic observers are intentionally absent.

## Usage

``` r
execution_receipt(
  scientific_plan_id,
  compute,
  sources,
  memory,
  kernel_version,
  task_partition_id,
  reduction_plan_id,
  precision = "double",
  numeric_contract = numerical_contract(),
  completion_status = "complete",
  task_count = 1L,
  completed_task_count = NULL,
  elapsed_seconds = 0,
  blas = list(vendor = "unknown", requested_threads = 1L, observed_threads = NA_integer_),
  domain_signature = NULL,
  observed = .empty_execution_observations()
)
```

## Arguments

- scientific_plan_id:

  Stable identity of relation, frame, pairing, and query/materialization
  semantics.

- compute:

  An `effect_compute_policy`.

- sources:

  A nonempty list of `effect_source_capabilities` values.

- memory:

  An `effect_memory_plan`.

- kernel_version, task_partition_id, reduction_plan_id, precision:

  Nonempty execution identity strings.

- numeric_contract:

  An `effect_numerical_contract`.

- completion_status:

  One of `complete`, `planned`, `failed`, or `interrupted`.

- task_count, completed_task_count:

  Nonnegative exact whole task counts.

- elapsed_seconds:

  Nonnegative finite elapsed wall time.

- blas:

  A list naming the BLAS vendor and positive thread count.

- domain_signature:

  Optional exact neural-domain signature. Compiler receipts always
  provide it; standalone receipts may omit it.

- observed:

  Canonical observed execution facts. Planned receipts use an empty
  observation record that is populated by the coordinator.

## Value

An immutable-by-convention execution receipt.
