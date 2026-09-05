# Production mutation evidence

The standalone algebraic oracle distinguishes sixteen wrong formulas. That
is useful fixture evidence, but it is not a test of the production program.
This gate additionally mutates the actual production function bodies and
runs their named primary tests in disposable R processes:

```sh
Rscript benchmarks/predictive-geometry/run-mutations.R . benchmark-results/predictive-geometry-mutations all
```

The catalog covers twenty variants: inverse versus direct kernel penalty;
penalty sign; fitting pooled null directions; compressed and full-form
preclipping; normalization before truncation; missing packed-coordinate
weights; missing score factor two; trace, missing or doubled prediction
cost; pruning on final test gain; endpoint-label and content-hash origin
shortcuts; training on all parent partitions or outer folds; positional
measurement matching; stale source or fold caches; and full geometry
materialization concealed behind a reduced result.

For each variant, a fresh worker loads the current package and runs the
unmodified primary test file. It then runs a no-op instrumented copy of the
same production closures, to check that instrumentation itself changes no
test result. Finally, exact AST replacements install the wrong production
closure in the worker namespace and rerun the file. Every replaced function
must actually be called, and at least one test naming the declared primary
T-ID must fail. A parser error, missing AST target, failed baseline, failed
no-op control, or unused mutation is not a kill.

The workspace's production files are never rewritten. Namespace replacements
are scoped and restored, and each variant's process exits before the next
starts. The artifact records the original body digests, complete changed
bodies, hit counts, primary semantic failures and test-file digest. The
combined matrix binds the R source and all three mutation harness files.

Two small production witnesses accompany the catalog. The first fits a
signed full form with compressed amplitude 1; clipping the full form first
would instead give 1.5. The second changes training edge weights while
holding endpoint identities, source revisions and all dimensions fixed.
It catches a stale fold cache even when its storage-size key happens to be
unchanged. Both pass on the intended implementation.

The machine-readable matrix and per-variant logs live under
`benchmark-results/predictive-geometry-mutations/`. The gate requires every
critical row to be `KILLED`; there is no coverage-percentage substitute for a
surviving scientific error. New failures exposed by a survivor require a
deterministic regression witness before the gate can pass.
