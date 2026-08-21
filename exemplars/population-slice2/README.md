# population-slice2 — searchlight-level transport, and what `η_transport` is worth

**WS-E slice 2.** Twelve subjects of OpenNeuro
[`ds003745`](https://openneuro.org/datasets/ds003745) v2.1.1 (CC0), carried
from ~65,000 native conservative searchlights each onto **1,241 shared
consensus grey-matter group nodes**, under **two declared transports**, with
one group model fitted under each — and then the question slice 1 could not
ask: *does the functionally-informed transport carry the population better than
the anatomical one, and how would we know?*

Slice 1 does population form at region level on six Haxby subjects, where the
transport is a 3×3 identity: no displacement, no spreading, no sink, and
nothing for `population-form-v1` §7 to diagnose. Slice 2 is the hard half. Here
the transport moves mass a median of **5.3–6.7 mm**, spreads it over **2.4–2.6
effective group nodes**, and sends **18–25 % of every subject's territory to
the sink** — every §7.5 diagnostic is a number that had to be computed rather
than a zero that could be asserted, and the exercise turned up a way of gaming
`η` that those six diagnostics do not catch.

- **Dataset:** ds003745, *An fMRI Dataset on Social Reward Processing and
  Decision Making in Younger and Older Adults* (Smith, Ludwig, Dennison, Reeck
  & Fareri). **CC0.**
- **Derivatives:** fMRIPrep **21.0.2**, space `MNI152NLin2009cAsym`.
- **Cohort:** 12 subjects — 6 younger (ages 20–34), 6 older (63–80). See
  [`DECISION.md`](DECISION.md) §2 for why the split is deliberate and §7 risk 9
  for why **no younger-versus-older contrast is reported anywhere**.
- **Evidence:** task `trust`, 9 conditions × 5 runs.
- **Transport source:** task `sharedreward`, 2 runs — a **different task**, so
  `P^F` never sees the data `η` is evaluated on.

Read [`DECISION.md`](DECISION.md) for how the dataset was chosen and for the
eleven-item risk list this code is the answer to. See
[Fetching](#fetching) for the downloader and [Running it](#running-it) for the
pipeline.

---

## The headline: `η` said one thing and its components said another

```
η_transport = R(P^F) − R(P^A) = +0.926794        signed, unclamped

  P^A   V^C = +475.8643   V^W = +1408.7529 ± 3167.5   R = +0.337791
  P^F   V^C = +442.6487   V^W =  +350.0347 ± 1964.0   R = +1.264585

  null band, 200 permuted-fingerprint draws:
      mean +0.014251   sd 0.104664   q95 +0.194488   max +0.404403
      49.5 % of draws negative
  η exceeds 200 of 200 draws        one-sided p = 0.0050
```

Read as a headline, that is a functional transport beating an anatomical one
by nearly a whole unit of consensus share, past every draw of its own null.
**It is nothing of the kind, and the numbers that say so are `V^C` and `V^W`
reported separately** — which is exactly why `population-form-v1` §7.3 requires
them and forbids reporting the ratio alone.

| | `P^A` | `P^F` | change | where `P^F` sits in the null |
|---|---|---|---|---|
| **`V^C`** cross-participant | +475.86 | +442.65 | **−7.0 %** | **47.5th percentile** — 95 of 200 draws below it |
| **`V^W`** cross-partition | +1408.75 | +350.03 | **−75.2 %** | **0th percentile** — below *all* 200 draws |

The functional transport built **no more cross-participant consensus than a
scrambled fingerprint would have**: its `V^C` lands at the null median. What it
did do, uniquely and by a wide margin, is destroy within-subject cross-partition
energy — its `V^W` is lower than every one of the 200 randomized transports.
`η` rose because the denominator collapsed.

The mechanism is not mysterious. The true atlas sends functionally *similar*
voxels to the *same* group node; averaging them together is exactly what
removes the local geometric diversity that carried the cross-run reproducible
signal. A scrambled atlas mixes dissimilar voxels and preserves more of it. So
the transport that looks best on `η` is the one that has thrown away the most
of what `η`'s denominator is supposed to measure.

§7.4 measures an adversarial transport that wins `η` by **discarding**
territory into the sink. Slice 2 blocks that attack by construction — `P^F`'s
sink is `P^A`'s, row for row, to 5.6e-16. What it found instead is a **second
attack the six diagnostics do not cover: winning `η` by suppressing `V^W`.**
Every one of the six is clean here. The sink is pinned. Coverage is 12 of 12 at
every node. Displacement is bounded by the support radius. Entropy is healthy.
And `η` is still meaningless as stated.

Two further reasons not to read `R(P^F) = +1.26` as "126 % consensus":

- **`R > 1` is out of range for what the ratio claims to be.** `V^C` is
  unbiased for `‖μ‖²` and `V^W` for `‖μ‖² + tr Σ_B`, so `R ≤ 1` whenever the
  quantities mean what their names say. §7.1 permits a share outside `[0,1]`
  because cross-fitting makes both terms noisy — that permission is about
  noise, not a licence to read 1.26 as a percentage.
- **`V^W` is 0.4 standard errors from zero** (+1408.75 ± 3167.5 over the twelve
  per-participant cross-partition products; `P^F`'s is 0.2 SE). It clears this
  slice's declared floor of `V^W > 0`, and a bare sign test is a weak reading
  of §7.1's "bounded away from zero by a declared criterion". The share has a
  denominator that is barely determined, and the receipts say so: `V_W_PA_se`
  and `V_W_PF_se` are recorded beside every `V^W`.

**What survives.** The permutation test itself is valid — it compares the same
statistic under label permutation on the same data, so `p = 0.005` genuinely
says the true fingerprint is doing *something* no scrambled one does. The
component table says what: it suppresses `V^W`. That is a real, reproducible
property of the operator, and it is not the property `η` is named after.

**A calibration check that did pass.** The null band's median `R` is
**+0.3422** against `R(P^A) = +0.3378`. A permuted-fingerprint transport
recovers the anatomical consensus share to within half a percent — which is
what a well-calibrated null should do, and is good evidence that the machinery
is measuring what it claims and the finding above is not an artefact of it.

For scale, `population-form-v1` §7.2's own fixture reports `R(P^A) = +0.3573`;
the anatomical transport here lands at **+0.3378**. The contract's synthetic
baseline and a real twelve-subject MNI cohort agree to about two points of
consensus share.

The second independence axis is discussed under
[Two cross-fit axes](#two-cross-fit-axes); it returns `NA`, and that is also a
result.

## What the transport actually is

Every subject's functional data already lives on **one shared 66 × 78 × 61 MNI
lattice** (2.973 × 2.973 × 3.22 mm), so unlike slice 1 there *is* voxel
correspondence between subjects and a searchlight-level transport is
meaningful. What differs between subjects is **coverage**, and it differs a
lot: over the twelve seven-run intersection masks the union is 77,064 voxels
and the consensus is 51,136 — **33.6 % of the union is covered for some
subjects and not others.**

| | |
|---|---|
| **native nodes** | one conservative searchlight per covered voxel, **radius 8 mm**, mean 75.5 voxels. 60,443–69,867 per subject. The frame is conservative in `conservative-geometry-v1` §2's sense: each voxel's mass is split among the searchlights containing it, so the column masses sum to 1 (max deviation ≤ 3e-15) |
| **group nodes** | **1,241**. Every 3rd lattice voxel in each direction (8.92 × 8.92 × 9.66 mm), restricted to the 12-subject consensus coverage **and** to grey matter — the across-subject mean fMRIPrep GM probseg, box-averaged from the 1 mm anatomical grid onto the functional grid, > 0.25 |
| **P^A** | `anatomical_transport()`: hard nearest-centre assignment, mass 1 on the nearest group node **within 9 mm**, otherwise all sink. Row entropy is identically 0 |
| **P^F** | softmax over a wider **12 mm** neighbourhood weighted by fingerprint similarity — *restricted to the rows `P^A` placed*, so its sink is `P^A`'s sink, row by row and bit for bit |

Radii are in **millimetres throughout**, measured in the domain's own
coordinates. The voxels are anisotropic (2.973 × 2.973 × 3.22 mm), so a radius
in voxels would build spheres squashed along z and would bias `P^A` in the same
direction ([`DECISION.md`](DECISION.md) risk 3).

### The sink is pinned, deliberately

`population-form-v1` §7.4 measures an adversarial transport that reports
η = +0.167 — 3.3× the honest functional gain in the contract's own fixture —
purely by sending **83 % of native territory to the sink**, discarding whatever
disagrees. Slice 2 removes that degree of freedom by construction: a native row
is carried if and only if `P^A` carries it, and `max |P^A_sink − P^F_sink|` is
asserted at ≤ 1e-12 for every subject. The two transports carry **exactly the
same territory**, so η can only be responding to where the mass goes among the
nodes they both reach.

The 12 mm support is nonetheless a real asymmetry — `P^F` may redistribute
where `P^A` may not, and spreading smooths. At 9 mm the mean support is only
1.8 group nodes (the grid is grey-matter gated and consensus restricted, so it
is sparser than its nominal spacing suggests), which is too few destinations
for a transport to express any preference at all; the asymmetry is what buys
the question its subject matter. The null band is what prices it: a permuted
`P^F` spreads exactly as much, so everything the smoothing is worth falls
inside the band and only what the correspondence is worth falls outside.

### The fingerprint, and its cross-fit

A voxel's fingerprint is its 7-vector of `sharedreward` condition betas — only
the seven event-level conditions present in **all 24 runs**
([`DECISION.md`](DECISION.md) risk 7) — centred and unit-normed, so a dot
product between two fingerprints is their correlation. A group node's
fingerprint is the mean of the **other eleven** subjects' fingerprints over the
voxels within 6 mm of it: leave-one-subject-out, so a subject's transport is
never built from that subject's own data even inside the fitting task.

Seven conditions are *read off*, but **all fifteen `sharedreward` trial types
are modelled** — the six block-level levels and the two neutral event types too
sparse to use as conditions are regressors of no interest, not omissions. This
is not fastidiousness. Measured on `sub-104`, fitting the seven conditions
alone instead of all fifteen levels changes the resulting seven-condition
fingerprint beyond recognition: the correlation between the two fingerprints
has a **median of 0.13** across 63,860 voxels, with a quarter of voxels below
−0.20. An event that happens and is not modelled does not vanish; its variance
is absorbed by whatever regressor it correlates with, which here is a regressor
of interest. (`trust` is unaffected — its only levels are the nine conditions
plus `missed_trial`, so nothing is added.)

`provenance$cross_fit` records `"task-sharedreward"`, `plan_population()`
carries it into `subject_index$cross_fit`, and 04 refuses to evaluate η on any
partition named there. §7.2 measures a *circular* transport — the same
estimator, differing only in being fitted on a partition it is also evaluated
on — reporting **3.15×** the honest gain and landing essentially at the oracle
ceiling. In a result object the two are indistinguishable; the provenance
string is the only thing that tells them apart.

The softmax temperature is **fixed at 6, not tuned.** Tuning it against the
held-out `trust` data is precisely that circularity.

---

## The six §7.5 diagnostics

`population-form-v1` §7.4 is normative: **η may not be printed, plotted or
returned without these in the same object.** All six are in
`results/population-slice2-receipts.csv`, per subject where they are
per-subject quantities. Ranges below are across the twelve subjects.

| diagnostic | `P^A` | `P^F` |
|---|---|---|
| **V^C / V^W**, reported separately and not only as their ratio | +475.86 / +1408.75 ± 3167.5 | +442.65 / +350.03 ± 1964.0 |
| **η null band** | — | mean +0.0143, sd 0.1047, q95 +0.1945; η ranks 200 / 200 |
| **displacement** median / p90 / max | 5.30 / 7.39 / 8.92 mm | 6.52 – 6.73 / 9.79 / 11.89 mm |
| — row-mass-weighted mean | 5.22 mm | 6.59 mm |
| **entropy**, row mean, nats | **0.000** exactly — a hard assignment | 0.757 – 0.851 |
| **perplexity**, the row mean of exp(H) | 1.000 | **2.407 – 2.642** |
| — exp(row mean H), *the other summary* | 1.000 | 2.132 – 2.342 |
| **sink territory**, data-free, from the operator alone | **18.2 – 25.2 %** | identical, by construction |
| **sink budget**, per subject, ledger units | in the receipts, per query | in the receipts, per query |
| **all-sink rows**, excluded from displacement and entropy and counted here | 11,319 – 17,603 | identical |
| **group-node subject coverage** | **min 12 of 12**, 0 nodes below the declared floor of 12 | same |

Displacement and entropy follow §7.5 literally: the group part of each row is
renormalized to sum 1, displacement is ‖center(x) − Σ_j p̃_xj center(j)‖,
entropy is −Σ_j p̃_xj log p̃_xj in nats, and rows with zero group mass are
**excluded and counted separately**.

Three notes on reading the table:

- **The two entropy summaries are both reported because they disagree**, and
  §7.5 says an implementation must label which it means. The row mean of exp(H)
  is 2.407–2.642; exp of the row mean of H is 2.132–2.342. Jensen guarantees
  the first is the larger. `perplexity_mean_PF` and `exp_mean_entropy_PF` are
  separate receipt rows for exactly this reason.
- **`P^A`'s displacement is not zero**, unlike slice 1's and unlike the
  contract fixture's. A hard nearest-centre assignment still moves a voxel to
  its group node, and on an 8.92 × 8.92 × 9.66 mm grid that is 5.3 mm at the
  median and 8.92 mm at worst — bounded above by the 9 mm assignment radius, as
  it must be. `P^F`'s maximum, 11.89 mm, is likewise bounded by its 12 mm
  support.
- **Every one of the 1,241 group nodes is reached by all 12 subjects.** That is
  a consequence of drawing the grid from the consensus mask, not a discovery,
  which is why `coverage_nodes_below_floor` is 0. The floor itself (12) is a
  declared choice: `population-form-v1` §14.3 leaves the threshold to the
  maintainer and requires only that the number be reported and the marking
  mechanism exist.

### A seventh diagnostic, and why six were not enough

The six catch a transport that wins `η` by **discarding** territory. None
catches one that wins by **concentrating** it: a `P^F` that piled the brain
onto a handful of group nodes would have full subject coverage, bounded
displacement, respectable entropy and `P^A`'s exact sink, and would pass all
six while making the group grid a fiction. So the arriving mass per group node
is summarized by its inverse participation ratio — `(Σx)² / Σx²`, the effective
number of group nodes actually carrying the territory:

| | `P^A` | `P^F` |
|---|---|---|
| effective group nodes, of 1,241 | **1,156 – 1,166** (0.94 of the grid) | **939 – 1,084** (0.82) |
| heaviest node's arriving mass | 87 | 107 – 159 |
| lightest node's arriving mass | 14 – 26 | 3.0 – 8.8 |
| nodes receiving zero mass | 0 | 0 |
| total arriving mass | identical to `P^F`'s | identical to `P^A`'s |

`P^F` concentrates moderately, which is what choosing a destination *means*,
and it does not collapse: no node is starved and 82 % of the grid is still
effective. The ratchet test holds this at a floor of 0.5, and the identical
totals are the sink control read from the group-node side.

That clears `P^F` of the concentration charge — and `η` was *still*
uninterpretable, for the `V^W` reason in the headline. **Seven diagnostics were
not enough either.** What caught it was the one thing §7.3 requires beyond any
diagnostic: reporting `V^C` and `V^W` separately instead of only their ratio.

---

## The identity acceptances

The same identities slice 1 asserts, at the same tolerance, on an operator that
is no longer trivial. Every number is read out of the fitted objects, never
re-derived, so a change in the package moves them.

| identity | `P^A` | `P^F` | tolerance |
|---|---|---|---|
| budget preservation — the fit's own certificate, relative to the ledger L1 norm | 1.57e-14 | 1.61e-14 | 1e-12 |
| executor plumbing — `transport_values()` of the subject's own `contrast_energy()` ledger equals the group fit's response | **0.00e+00** | **0.00e+00** | 1e-12 |
| sum over group nodes **and the sink** = transported native total, relative | 1.92e-14 | 2.01e-14 | 1e-12 |
| Σ_u Θ — aggregated ledger vs the group coefficient of the summed response, relative | 1.27e-15 | 1.77e-15 | 1e-12 |
| **commutation (claim 3)** | **1.78e-15** | **8.88e-16** | 1e-12 |
| coherent + configuration = total | 4.69e-13 | — | 1e-12 |
| `P^F` sink equals `P^A` sink, per native row | — | **5.55e-16** | 1e-12 |

**The commutation acceptance is the one that matters.**
`estimate_population()` contracts the `choice_friend − choice_computer`
contrast into a *single* packed operator and transports one number per node:
**query, then transport**. `materialize_population()` transports all **45**
packed coordinates of the form and `rdm()` contracts that same edge out of them
afterwards: **transport, then query**. Different arithmetic, and they agree to
1.8e-15 at 1,241 group nodes under a transport with 5–7 mm of displacement and
a 21 % sink. Slice 1 could only demonstrate this on a 3×3 permutation.

The executor-plumbing row is **exactly zero** and is deliberately *not* the
commutation claim: both sides query before transporting, so it is the same
arithmetic twice. It is recorded because an executor that started rescaling or
reordering rows would show up there first.

Two of these are asserted **relatively** where slice 1 asserted them
absolutely. Slice 1's region totals were O(1) and an absolute 1e-12 said the
same thing; here the ledgers run to O(100) summed over 1,242 rows, and an
absolute 1e-12 would be asserting ~1e-14 *relative* — below the accumulation
floor of the summation itself, so it would be a tolerance on summation order
rather than on the transport. The scale used is the one
`estimate_population()`'s own certificate uses,
`relative_to_ledger_l1_norm`.

### What the sink costs a reader who ignores it

Slice 1 recorded "sum over nodes with the sink deleted" as *still zero*,
because its region map had full coverage. Here the same rows are large:

```
node sum without the sink   1.544  (P^A) / 1.587  (P^F)   relative
Θ without the sink          2.949                          relative
```

Neither is asserted, because neither is an identity. They are the size of the
error a reader makes by dropping the sink row, relative to the L1 norm of what
is left — **larger than the surviving ledger itself**, because a signed
contrast ledger cancels heavily across nodes while the sunk part does not. The
sink is not decoration, and these are the receipts that say so.

---

## Two cross-fit axes

[`DECISION.md`](DECISION.md) §5.3 promised two independence levels, and risk 8
says what a disagreement between them would mean: a negative `η` across task
with a positive `η` across run would indicate the fingerprint is capturing
task-specific rather than anatomical idiosyncrasy. Both are computed.

| axis | fitted on | evaluated on | result |
|---|---|---|---|
| **across task** (headline) | `sharedreward` runs 1–2 | `trust` {1,2} vs {3,4,5} | `η = +0.9268`, and see above for what it is made of |
| **across run** | `trust` run 1 | `trust` {2,3} vs {4,5} | `R` and `η` are **`NA`** |

### The across-run axis returns `NA`, and that is the result

```
P^A    V^C = −426.82    V^W = −30826.03 ± 20537.2    R = NA
P^F2   V^C = −144.57    V^W = −17476.55 ± 10284.9    R = NA
```

`V^W` came out **negative** on this split, at **−1.5 standard errors** from
zero (`P^F2`: −1.7). `population-form-v1` §7.1 is explicit about exactly this
case: `V^W` is itself a cross-partition product, it can be zero or negative
when nothing reproduces, and *the share is then `NA`, not a large number*. So
`R` and `η` are `NA`, **no null band is drawn**, and no p-value is
manufactured for a quantity that does not exist. The receipts record
`eta_across_run_shares_formable = 0` and the null CSV has zero rows, so the
absence is itself a recorded fact rather than a gap.

What the axis says is *not* "the functional transport does not help across
runs". It is: **four runs split 2 + 2 do not contain enough reproducible
searchlight geometry for a consensus share to be defined at all.** Both
components are still reported, because they are the entire evidence for the
`NA`, and both structural assertions still run and still pass — `P^F2`'s sink
matches `P^A`'s to 5.6e-16, and rebuilding `P^F2` under the identity
permutation reproduces the sealed operator exactly.

The contrast with the across-task axis is instructive rather than damning: that
one uses all five runs split 2 + 3, and its `V^W` is +0.4 SE from zero rather
than −1.5. Neither is a well-determined denominator. The across-task axis
clears the declared floor and the across-run axis does not, and the difference
between them is one run and a coin flip.

### Why the partition is `{2,3}` vs `{4,5}`

Both held-out sides must support a *cross-validated* conservative form, which
needs two runs each. Five runs, one spent on the fingerprint, leaves exactly
2 + 2. Fitting on runs `{1,2}` as [`DECISION.md`](DECISION.md) §5.3 sketched
would leave `{3,4,5}`, which cannot be halved into two two-run partitions — the
same property that eliminated the AOMIC datasets in
[`DECISION.md`](DECISION.md) §4. So the across-run fingerprint rests on a
single run and is correspondingly weaker; that is the honest price of keeping
both evaluation sides cross-validatable, and it is part of why this axis had
less to work with.

Both axes use the same estimator on different data, not two estimators:
`eta-common.R` holds the held-out form construction, `V^C`/`V^W`, the softmax
transport builder and the null summariser, and both scripts assert that
rebuilding their `P^F` under the **identity** permutation reproduces the
operator `location_transport()` sealed — to 0.00e+00. Without that assertion a
null band could be a band around a subtly different estimator, and `η`'s rank
inside it would mean nothing.

---

## The group form

Five zero-sum contrasts over the nine `trust` conditions, read at every group
node in one pass:

| query | what it is |
|---|---|
| `choice-vs-outcome` | the three choice cells vs the six outcome cells |
| `social-vs-computer-choice` | choosing with friend or stranger vs with the computer |
| `recip-vs-defect` | reciprocated vs defected outcomes |
| `friend-vs-stranger` | friend cells vs stranger cells |
| `friend-vs-computer-choice` | a bare pairwise difference — i.e. exactly an RDM edge, present so the commutation acceptance has an edge to compare against |

Descriptively, at uncalibrated |t| > 3 with 11 residual df, of 1,241 nodes:
`choice-vs-outcome` **393** nodes (max |t| 6.59), `social-vs-computer-choice`
**19**, `friend-vs-computer-choice` **11**, `recip-vs-defect` **1**,
`friend-vs-stranger` **0**. The phase contrast is a large BOLD effect and
behaves like one; the subtle social contrasts are near null at n = 12 and are
reported as near null.

Cross-fitted heterogeneity gives `n_eff` ≈ **2.68** (`P^A`) and **3.03**
(`P^F`) effective between-subject modes. One draw, descriptive, not an estimate
of a between-subject trace.

The two transported component ledgers carry the names
`population-form-v1` §8.1 requires — `native_coherent_ledger` and
`native_configuration_ledger`, never the bare words, which would invite reading
them as properties of a group node's geometry rather than as ledgers of native
node coherence.

---

## Honest boundaries

1. **The headline `η` does not mean what it is named after.** See
   [The headline](#the-headline-η-said-one-thing-and-its-components-said-another).
   `V^C` is at the null median; the entire effect is a 75 % collapse in `V^W`.
   No claim is made here that a functionally-informed transport carries this
   population better than an anatomical one, and none should be read out of
   these files.
2. **This slice's `V^W` floor is a bare sign test, and that is too weak.**
   §7.1 asks for `V^W` "bounded away from zero by a declared criterion"; the
   criterion declared here is `V^W > 0`, and `V^W` clears it by 0.4 of its own
   standard error. A criterion in units of the standard error would be the
   better reading and is a change worth making — it would have marked the
   across-task axis as unformable too, rather than only the across-run one.
   The standard errors are recorded (`V_W_PA_se`, `V_W_PF_se`,
   `V_W_PA_acrossrun_se`, `V_W_PF2_se`) so the stricter reading is available
   from the committed evidence without re-running anything.
3. **The six §7.5 diagnostics are necessary and were not sufficient.** All six
   are clean here — pinned sink, full coverage, bounded displacement, healthy
   entropy — and `η` was still uninterpretable. Slice 2 adds a **seventh**, the
   effective number of group nodes carrying the arriving mass
   (`group_node_effective_A/F`: 1161 of 1241 under `P^A`, 940–1084 under
   `P^F`), which rules out the concentration analogue of the §7.4 sink attack.
   It does not rule out the `V^W`-suppression attack, which only the separated
   components caught.
4. **Twelve subjects is twelve subjects.** Every interval here is
   **uncalibrated**, and `population_uncertainty()` says so in every row it
   emits. The one-sided p of 0.005 is a permutation p against one null, on one
   dataset, under one fingerprint.
5. **The 6 + 6 age split is a design choice, not a sample.** No
   younger-versus-older contrast is computed anywhere in this slice. The older
   group's larger coverage propagates into slightly larger sink territory
   (mean 22.4 % vs 20.3 %) and must not be read as a data-quality difference.
6. **The older group carries substantially more head motion.** With FD > 0.5 mm
   censored one volume at a time, `sub-129` loses **61–112 of 217 volumes per
   run** and its per-run residual df falls to 76–127, against 180–189 for the
   quietest subjects. Per-subject censoring counts, mean FD and minimum
   residual df are all in the receipts. This inflates those subjects' beta
   variance and is a reason to read the *group* fit rather than the subject
   spread.
7. **The functional transport is simple, and deliberately so.** Seven condition
   betas from two runs, a correlation, a fixed-temperature softmax over a 12 mm
   ball. `population-form-v1` §9.2 keeps transport *learning* outside crossform;
   this is an exemplar-built operator and its construction is recorded verbatim
   in `provenance$details`. A hyperalignment-grade fingerprint would very likely
   behave differently, and nothing here is evidence about what such a transport
   would be worth.
8. **The GLM models neither slice timing nor autocorrelation.** The fMRIPrep
   derivatives are not slice-time corrected and no reference slice is declared,
   so volume times are used as-is. OLS betas are unbiased under this; their
   standard errors would not be, which is one more reason the uncertainty here
   is labelled uncalibrated. aCompCor is deliberately not used: its component
   count is run-specific (the confound tables run 154 to 364 columns wide), so
   "all of them" would give different runs different model complexity for no
   stated reason.
9. **The group grid is drawn from the consensus mask**, so full subject
   coverage of every group node is built in rather than discovered. A grid drawn
   from the union would put nodes below the floor and would exercise the
   *marking* mechanism; this one exercises the *reporting*.
10. **`sharedreward` is used for nothing but the fingerprint.** No `trust`
    number anywhere in this slice is computed from it, and no `sharedreward`
    result is reported.

---

## Fetching

```sh
./fetch.sh --dry-run      # list what would be fetched, print the total
./fetch.sh                # fetch all 392 files (7.21 GB) into ./data
./fetch.sh --verify       # re-check sizes and md5s of what is on disk
./fetch.sh --verify --repair   # ...and delete whatever fails, so a plain
                               #    run afterwards re-fetches it
```

Everything lands under `data/`, which is git-ignored, with BIDS relpaths
preserved — so `data/` is a partial BIDS tree with a
`data/derivatives/fmriprep/` subtree beside it. Downloads go to a `.part` file
and are promoted only after both the size **and** the md5 match, so an
interrupted run never leaves a truncated NIfTI in place. The corollary: a plain
run trusts the byte count, so a file that is corrupt *at the right size* is
caught only by `--verify`.

`manifest.csv` has one row per file — `relpath`, `role`, `tier`, `source`,
`bytes`, `md5`, `url`, `s3_uri`, with `sha256` intentionally empty because S3
does not publish it. To change the subset, regenerate the manifest against the
snapshot rather than hand-editing byte counts.

## Running it

```sh
./fetch.sh                 # 392 files, 7.21 GB, md5-verified, resumable
bash run-all.sh            # ~80 min, ~4 GB peak RSS
```

or stage by stage — each is idempotent and skips work already done:

```sh
export SLICE2_DIR="$PWD"
Rscript 01-glm.R            # ~3 min   -> data/derived/<subject>_betas.rds
Rscript 02-transports.R     # ~1 min   -> data/derived/transports.rds
Rscript 03-population.R     # ~35 min  -> results/population-slice2-receipts.csv
Rscript 04-eta-transport.R  # ~11 min  -> results/population-slice2-eta*.csv
Rscript 05-eta-across-run.R # ~12 min  -> results/population-slice2-eta-across-run*.csv
```

**03, 04 and 05 must run in that order**: 03 writes the receipts CSV from
scratch and 04 and 05 merge their rows into it.

| environment variable | effect |
|---|---|
| `SLICE2_DIR` | this directory; inferred from the script path otherwise |
| `SUBJECTS` | comma- or space-separated subset (01 only) |
| `ETA_DRAWS` | null-band draw count (04, 05); `ETA_DRAWS=5` is a smoke run |

If fewer than twelve subjects have been fetched, 02 says so and runs on the
ones present rather than failing or padding.

Everything under `data/` is git-ignored; everything under `results/` is
committed. `tests/testthat/test-population-slice2.R` is the ratchet that holds
it. It recomputes nothing — that would need the 7.21 GB back — and asserts that
the committed record still says what this README claims: twelve subjects, every
identity inside its acceptance, the two sinks identical, η reported signed with
its null band, and all six §7.5 diagnostics present. It skips cleanly when
`results/` is absent. Notably it does **not** assert η > 0; a test that
required a positive η would be the clamping §7.3 forbids, moved into the test
suite.

## Files

| file | what it is |
|---|---|
| `DECISION.md`, `manifest.csv`, `fetch.sh` | E12's dataset decision, the pinned 392-file manifest, and the verifying downloader |
| `00-common.R` | cohort, conditions, query bank, geometry settings, tolerances, path helpers, the receipts recorder |
| `01-glm.R` | per-subject per-run GLM → per-condition betas on the shared lattice |
| `02-transports.R` | group grid, `P^A`, `P^F`, and the data-free §7.5 diagnostics |
| `03-population.R` | two `plan_population()` + `estimate_population()` runs, and the identity acceptances |
| `eta-common.R` | the held-out forms, `V^C`/`V^W`, and the permuted-transport machinery 04 and 05 share |
| `04-eta-transport.R` | η across task, with its 200-draw null band |
| `05-eta-across-run.R` | η across run — the second independence axis |
| `run-all.sh` | all five stages in order |
| `results/*.csv` | committed evidence: receipts, transport diagnostics, both η tables, both null bands |

**Citation.** `dataset_description.json` asks that two papers be cited by
anyone using these data; both belong in any write-up of this slice.

- OpenNeuro **ds003745** v2.1.1, `doi:10.18112/openneuro.ds003745.v2.1.1`,
  **CC0**. Smith DV, Ludwig RM, Dennison JB, Reeck C, Fareri DS. Ethics:
  Temple IRB #24452.
- Fareri DS, et al. *NeuroImage* (2022). `10.1016/j.neuroimage.2022.119267`
- Smith DV, et al. *Scientific Data* (2024). Preprint `10.31234/osf.io/k7d56`
- Esteban O, et al. *fMRIPrep: a robust preprocessing pipeline for functional
  MRI.* Nature Methods 16:111–116 (2019). Pipeline version 21.0.2.
