# Construct a crossform-owned workspace plan

The hard budget covers package-owned live objects and conservative
temporary overlap. Process baseline RSS and absolute RSS are
observations, not hidden additions to the workspace model.

## Usage

``` r
memory_plan(
  frame_bytes = 0,
  resident_source_bytes = 0,
  source_handle_bytes = 0,
  source_block_bytes = 0,
  relation_block_bytes = 0,
  atom_block_bytes = 0,
  local_state_bytes = 0,
  output_bytes = 0,
  contraction_bytes = 0,
  replacement_copy_bytes = 0,
  serialization_overlap_bytes = 0,
  reorder_buffer_bytes = 0,
  checkpoint_buffer_bytes = 0,
  workers = 1L,
  n_active = workers,
  safety_factor = 1.25,
  budget_bytes = NULL,
  measured_workspace_bytes = NULL,
  baseline_rss_bytes = NULL,
  peak_rss_bytes = NULL
)
```

## Arguments

- frame_bytes:

  Resident dense or sparse frame storage.

- resident_source_bytes:

  Resident source objects owned by crossform.

- source_handle_bytes:

  Execution-owned source handles and bookkeeping.

- source_block_bytes, relation_block_bytes, atom_block_bytes:

  Per-active-task data blocks.

- local_state_bytes:

  Durable component-dependent local relation state.

- output_bytes:

  Durable in-memory result components.

- contraction_bytes, replacement_copy_bytes:

  Per-task contraction and copy-on-modify temporaries.

- serialization_overlap_bytes, reorder_buffer_bytes,
  checkpoint_buffer_bytes:

  Optional per-task execution buffers.

- workers, n_active:

  Positive worker and active-task counts.

- safety_factor:

  Multiplicative headroom, at least one.

- budget_bytes:

  Optional hard owned-workspace budget.

- measured_workspace_bytes:

  Optional observed package-owned live bytes.

- baseline_rss_bytes, peak_rss_bytes:

  Optional process RSS observations.

## Value

An immutable-by-convention `effect_memory_plan`.
