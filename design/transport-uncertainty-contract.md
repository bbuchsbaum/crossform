# Transport uncertainty boundary

Population inference in crossform is conditional on the realized transport.
This is true for fixed anatomical operators, externally supplied operators,
and functionally estimated operators. A cross-fitted functional transport
limits circular reuse of response data; it does not integrate over the
estimated operator, the fitting sample, or the fold assignment.

Every sealed transport therefore carries a canonical `conditioning` record:

- `source`, `operator_status`, `fitting_sample`, and `cross_fit_folds` state
  what was supplied and how it was fit;
- `inference_scope = "conditional_on_realized_transport"` and
  `marginal_over_transport = FALSE` prevent a cross-fitted result from being
  presented as marginal over alignment estimation;
- `uncertainty_propagated = FALSE` and `excluded_uncertainty` state what the
  current intervals and bootstrap do not cover.

For a concrete example, an HC3 population interval includes the observed
between-subject residual heteroskedasticity conditional on the transported
values. If those values came from a learned functional alignment, the interval
excludes variation from learning that alignment and from choosing its
cross-fitting folds. The subject-level wild bootstrap has the same boundary:
it resamples participant residuals while holding the realized transport fixed.

## Future extension point

The record reserves capability `transport_uncertainty_propagation`, currently
`not_implemented`. Implementing it requires all three declared ingredients:

1. a transport sampling law;
2. joint transport-response resampling that respects held-out fitting; and
3. a validated propagation operator.

Only an implementation satisfying that interface may set
`uncertainty_propagated` or `marginal_over_transport` to true. Until then those
fields are derived, immutable false values; caller-supplied claims are refused.
