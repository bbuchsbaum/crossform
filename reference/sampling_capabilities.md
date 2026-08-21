# Ask whether the analytic sampling law is available before provoking it

A scientist should be able to ask what the uncertainty channel will
grant instead of discovering refusals one at a time. This inspection
compiles the evidence-sampling admission for a plan and error source and
reports every unmet requirement with its remedy, without computing
anything.

## Usage

``` r
sampling_capabilities(x, fit = NULL)
```

## Arguments

- x:

  A compiled `effect_geometry_plan` (or crossnobis plan).

- fit:

  Optional `effect_relation_fit` supplying the error channel; omitting
  it probes the plan's bare relation, which has no channel.

## Value

An `effect_sampling_capabilities` list: `available`, the full
`capabilities` record, and a `reasons` data frame with one row per unmet
requirement (`reason`, `why`, `remedy`).

## What this cannot answer in advance

Every requirement reported here is a property of the plan and its error
channel, so it can be checked without touching neural values. One
requirement of
[`rdm_sampling_covariance()`](https://bbuchsbaum.github.io/crossform/reference/rdm_sampling_covariance.md)
is not: whether a *particular* measurement has enough residual degrees
of freedom for the number of residual directions its own support spends
variance on (capability `"sufficient_residual_df"`). That depends on the
local residual spectrum and can only be known once it is computed, so
`available = TRUE` here does not promise that every measurement will be
answerable.

## See also

[`rdm_sampling_covariance()`](https://bbuchsbaum.github.io/crossform/reference/rdm_sampling_covariance.md),
the call this inspection describes, and
[`catch_refusal()`](https://bbuchsbaum.github.io/crossform/reference/catch_refusal.md)
for branching on a refusal that has already happened.

Other sampling uncertainty:
[`rdm_sampling_covariance()`](https://bbuchsbaum.github.io/crossform/reference/rdm_sampling_covariance.md),
[`reconstruct_evidence()`](https://bbuchsbaum.github.io/crossform/reference/reconstruct_evidence.md),
[`sampling_covariance()`](https://bbuchsbaum.github.io/crossform/reference/sampling_covariance.md)

## Examples

``` r
# Ask before provoking: the fixture's fit retains a residual channel, so
# the analytic law is admitted.
example <- example_fmri_effects()
plan <- plan_geometry(
  example$fit$relation, example$frame,
  cross_partitions(example$fit$relation, independence = "independent")
)
capabilities <- sampling_capabilities(plan, example$fit)
capabilities$available
#> [1] TRUE
capabilities
#> <effect_sampling_capabilities>
#>   analytic sampling law: available 
#>   metric: fixed | partitions: equal | error channel: relation_fit 

# Probing the bare relation instead reports every unmet requirement with
# its remedy, rather than failing one refusal at a time.
without_fit <- sampling_capabilities(plan)
without_fit$reasons
#>                                  reason
#> 1                 missing_error_channel
#> 2 sampling_axis_missing_or_inconsistent
#>                                                                                                                                                                                                                                              why
#> 1 this evidence plan has only a precomputed relation and no error channel. Refit raw observations with `lm_relation_fit()` or supply a validated, identity-bound external error channel; beta matrices alone cannot recover residual uncertainty
#> 2                                                                                                                            no single sampling axis is declared, or the declared axis conflicts with the error channel's recorded sampling unit
#>                                                                        remedy
#> 1                            Refit raw observations with `lm_relation_fit()`.
#> 2 Declare one `sampling_axis` that matches the error channel's sampling unit.
```
