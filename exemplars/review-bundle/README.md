# Crossform scientific review bundle

This is the reviewer navigation page for Crossform's population-estimation and
interpretive-validation program. It separates proof, software agreement,
synthetic operating evidence, retrospective illustration, prospective plans,
and independent replication because those are different kinds of evidence.

The bounded public claim is:

> Crossform gives univariate contrasts, multivariate distances and fixed linear
> RSA, and population summaries one declared crossvalidated bilinear geometry,
> then decomposes each scale's reproducible effect into coherent and
> configurational estimand components.

The bundle does not claim to unify preprocessing, nonlinear or adaptive
statistics, every RSA method, transport learning, or all downstream inference.

## Start here

1. Verify every payload byte with `shasum -a 256 -c SHA256SUMS`, or run
   `Rscript build-bundle.R --verify=<tarball>` from any directory.
2. Read `PUBLIC-SCOPE.md` for the public claim and non-goals.
3. Open `review/articles/derivation.html` for the common-geometry theorem.
4. Open `review/articles/matched_simulation.html` for the equal-total,
   different-organization simulation.
5. Open `review/articles/population.html` for the executable population API,
   uncertainty, decomposition, diagnostics, scale profile, and influence
   workflow.
6. Inspect `SOURCE-TO-ARTIFACT.csv` for source hashes, artifact hashes, commands,
   seeds or replication rules, versions, and the boundary of each result.
7. Read `certification/artifacts/evidence-claim-registry.csv` for one owner and
   current status per claim, then `design/evidence-status-ledger.md` for every
   earned evidence row.

`MANIFEST.txt` records the source commit, branch, whether the working tree was
clean, the aggregate `R/` digest, build environment, every included path, byte
count, and SHA-256. A dirty-tree build is not hidden: its modified paths are
listed and the commit is explicitly insufficient to reconstruct it.

## Evidence map

| label | question answered | principal bundle paths | current boundary |
|---|---|---|---|
| **THEOREM** | Are contrast energy, crossnobis distance, and fixed linear RSA queries of one declared geometry? | `review/articles/derivation.html`; `design/common-geometry-equivalence.md`; `design/oracles/common-geometry-equivalence.R` | Fixed linear operators and declared metrics only. |
| **INTERNAL ORACLE** | Does production code agree with independent matrix constructions? | `design/oracles/`; package test paths named by the ledger | Bounded generated fixtures and tolerances, not external validation. |
| **EXTERNAL PARITY** | Does the mapped estimator agree with an independent implementation? | `exemplars/rsatoolbox-parity/`; `certification/artifacts/common-geometry-external-parity.csv` | Version-pinned fixed-linear synthetic case and bounded Haxby route; no broad RSA or scientific parity. |
| **MATCHED SIMULATION** | Can equal-total effects with different organization be distinguished, and how do intervals behave in declared regimes? | `review/articles/matched_simulation.html`; `certification/artifacts/matched-interpretability-*`; `population-calibration-*`; `population-interpretability-*` | Synthetic evidence only. Informative coverage supports no marginal population claim. |
| **RETROSPECTIVE ILLUSTRATION** | Do the typed routes run on versioned public data and expose diagnostics? | `exemplars/haxby2001/`; `exemplars/population-slice2/` | Descriptive/regression evidence; no frozen hypothesis or transport-superiority result. |
| **PROSPECTIVE PROTOCOL** | What eligible analysis would be run next? | `protocols/prospective/` | Frozen and synthetically rehearsed, but readiness is `BLOCKED`; no real-data result exists. |
| **COMPLETED REAL DATA** | Did an eligible frozen discovery protocol produce a result? | none | Not executed. |
| **INDEPENDENT REPLICATION** | Was the frozen target repeated through an eligible independent route? | none | Not executed. |

The labels are not a single ladder. A theorem does not become a simulation;
external parity does not become scientific validation; a rehearsed protocol is
not a result.

## What the synthetic courts establish

### Matched interpretability

The certification contains 48 paired seeds across three planted organizations,
three noise regimes, two sample sizes, three SNRs, and four scales. All
prespecified recovery, false-separation, conservation, ambiguity, and negative
control gates pass. The signed crossvalidated coherent share uses no positivity
screen. This licenses interpretation of those planted regimes only.

Key files:

- `certification/artifacts/matched-interpretability-parameters.csv`
- `certification/artifacts/matched-interpretability-metrics.csv`
- `certification/artifacts/matched-interpretability-checksums.csv`
- `benchmarks/matched-interpretability/05-certify.R`

### Population interval calibration

The interval court compares classical, HC3, and null-imposed wild bootstrap on
500 paired datasets in each of eight regimes. Gaussian calibration, HC3's
heteroskedastic improvement, heavy-tail behavior, and computational failure
rates are explicit. The informative-coverage regime is a failure boundary:
conditional intervals do not become marginal population intervals.

Key files:

- `certification/artifacts/population-calibration-parameters.csv`
- `certification/artifacts/population-calibration-results.csv`
- `certification/artifacts/population-calibration-verdicts.csv`
- `benchmarks/population-calibration/01-certify.R`

### Hierarchical interpretability

The population interpretability court uses 200 paired 24-subject replications.
It recovers the planted coherent/configuration ordering, exact total
conservation, component-aware interval behavior, scale profiles, and the
planted influential subject under realized transport. It remains synthetic and
does not license an empirical neuroscience interpretation.

Key files:

- `certification/artifacts/population-interpretability-results.csv`
- `certification/artifacts/population-interpretability-verdicts.csv`
- `benchmarks/population-interpretability/00-certify.R`

## Prospective boundary

`protocols/prospective/discovery-v1.json` and the matching Markdown freeze data
eligibility, exclusions, estimands, comparators, thresholds, deviations, and
failure handling. `replication-v1.md` requires either independent eligible data
or an external locked-bundle analyst route. The rehearsal executes public APIs
and records both successful outputs and expected refusal states, but uses
synthetic inputs.

Read `protocols/prospective/readiness-current.json` literally:

- `BLOCKED` means required real-data or independent-route prerequisites are
  absent.
- `READY` would authorize execution; it would not be a scientific result.
- `EXECUTED` requires immutable inputs, receipts, results, and deviations under
  the appropriate discovery or replication artifact root.

Discovery and replication cannot share a completion artifact. No bundle file
uses the terms demonstrated, validated, confirmed, or replicated for the
prospective rows.

## Rebuilding

From a source checkout with package dependencies, `pkgload`, `rmarkdown`,
`knitr`, and Pandoc available:

```sh
Rscript exemplars/review-bundle/build-bundle.R --force
Rscript exemplars/review-bundle/build-bundle.R \
  --verify=dist/crossform-review-bundle-<date>.tar.gz
```

The builder fails when a required source or artifact is absent, when an article
does not render, or when verification finds a missing, extra, or hash-mismatched
payload. It writes only under `dist/`. `ENVIRONMENT.txt` records the actual
build environment; `RUNNING-TESTS.md` gives source-checkout test commands.

Rebuild individual certification artifacts with the exact commands in
`SOURCE-TO-ARTIFACT.csv`. Those commands overwrite governed artifacts and
should therefore be run only on a frozen candidate tree. Their parameter and
checksum CSVs are included beside the results.

## Tamper check

Verification is intentionally content-addressed. To confirm that it detects a
substitution, extract a copy, alter any listed payload byte, and run
`shasum -a 256 -c SHA256SUMS`; the changed path must fail. The test suite also
constructs an altered bundle copy and requires the R verifier to exit nonzero.

## What is not claimed

- No empirical coherent/configurational brain-region conclusion is established.
- The ds003745 transport-efficiency ratio failed its retrospective
  interpretation audit; no functional-over-anatomical superiority is claimed.
- No interval marginal over informative coverage or transport learning is
  provided.
- No simultaneous/maxT scale-profile band, cross-node covariance estimator,
  random-field inference, or general multiple-comparison procedure is provided.
- No completed prospective real-data result or independent replication exists.
- No cross-package speed advantage is claimed.
- Version-pinned external parity is not universal parity with every RSA method.

The authoritative current-state surfaces are
`certification/artifacts/evidence-claim-registry.csv` and
`design/evidence-promotion-history.md`. If prose and an artifact appear to
disagree, the narrower artifact boundary governs and the inconsistency should
be reported.
