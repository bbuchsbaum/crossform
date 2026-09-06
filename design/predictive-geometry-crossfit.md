# Cross-fitting frozen geometry predictions

The internal `.geometry_crossfit()` orchestrator is deliberately built from
`fit_geometry()` and `score_geometry()`. For each positive-weight edge of the
declared undirected pairing, it trains on declared edges whose endpoints are
disjoint from that evaluation edge. It preflights every fold's actual origin
ancestry before any source values are read. A star pairing refuses even with
many partition labels; a disconnected pairing can be valid. Reversing the
spelling of an undirected edge changes no product.

The aggregate is the original edge-weighted mean of fold gains:

```
sum_e w_e [2 <F_e, G_test,e> - ||F_e||^2].
```

Every squared norm stays with its own prediction. The squared norm of an
average prediction is a different quantity. Fold ranks and full score
receipts remain available. Overlapping evaluation edges are not independent
replicates for standard errors; the result supplies none.

This estimates performance of the specified fitting procedure at those
fold-specific training sizes. It is not a score for a later all-data refit.
Model weights, rank and penalty are fixed here; nested training-only selection
is a separate orchestration step. This helper is internal while the package
keeps the primary user workflow to a model declaration and two verbs.

Sequential folds reserve retained prediction/score/schedule memory before per-fold
fit/score admission. Block execution uses a newly owned directory, removes
temporary signed geometry after fitting/scoring, retains each frozen fit,
and publishes the final cross-fit record
only after all folds succeed. Unsupported schedules, targets and storage
routes preserve the same refusals as the fixed-split path.
