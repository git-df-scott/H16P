# FASTRA D1 — fold sheets with full compactified counts

**D1 STATUS: OPEN. No four-origin-cycle field was verified.**
The new reproducible outcome is a **3+1 field on the finite-fold unfolding**,
which the earlier fold strike had not constructed. It has four cycles total.
Neither numerical continuation nor the proposed endpoint-sign test establishes
nonexistence of a fourth origin cycle on the entire component.

## Requested report

**FOLD SHEETS CONTINUED:** both positive-K directions and the separate
center-selected negative-K sheet; endpoints below. All radii in this table
are on the **positive horizontal section**, not the older curved maximum-x section.

| Direction | Accepted points | Endpoint c, approximately | Endpoint K, approximately | Origin fold radius |
|---|---:|---:|---:|---:|
| Positive K toward center | 10 | 0.968620633567205 | 1e-10 | 6.75794351533 |
| Positive K toward infinity | 25 | 1.59370290262731 | 7.06938707146707 | 1.73927494152e+18 |
| Negative K, continued in m | 29 | 0.336475146978924 | -2597546766505.68 | 40489446294 |

Negative endpoint actual coefficient: m = 10000000000000;
alpha = -m, and K = m(11c-5)/5-42. Decimal endpoints are numerical approximations.
The exact rational vectors at every endpoint are in its `events_*.jsonl` record.

**MAX ORIGIN COUNT ON ANY SHEET (compactified, reproduced): 3**, for example at

- c = 9688912553490597/10000000000000000,
- K = 1/512,
- beta = -1/10000000.

Maximum on the tested **beta=0 pair-present fields: 2**.
Maximum after the tested Hopf steps: **3**.
An exact semistable fold is not counted by a sign-change-only counter; the
rational vectors at nominal folds are approximate fold locations.

**FOUR-IN-ORIGIN-NEST FIELD: NONE verified.**

**SIGN-MAP VERDICT:** no additional pair-birth region detected on the tested grid;
**nonexistence is not established**. After precision repairs, the stored
outside/terminal comparisons have 198 agreements,
0 disagreements, and
0 unresolved comparisons. A terminal
sample or an integration failure is not a certified return-domain edge.
Binary128 matching terminal signs are explicitly an auxiliary diagnostic where
a complete return has not been established.

**D1 STATUS: OPEN.**

## One exact rational 3+1 field for reproduction

Coefficient convention:
P = a0 + a1*x + a2*y + a3*x^2 + a4*x*y + a5*y^2,
and Q uses the next six entries in the same order.

```
[0, 0, 1, 1, 1, 0,
 0, -2100097656250000000/56578038088396567,
 -1/10000000, -10, 11/5,
 9688912553490597/10000000000000000]
```

These are exact coefficients, not rounded decimal coefficients. Numerical roots
and multipliers from binary128 two-sided shooting are:

| Nest | Section radius, approximately | Stability | Multiplier |
|---|---:|---|---:|
| origin | 0.0949678324657 | U | 1.00000010212 |
| origin | 5.26090310172 | S | 0.999995858247 |
| origin | 8.56090804672 | U | 1.00000688397 |
| remote | 10395387.4971 | S | 0.994434780972 |

Origin radii use the positive horizontal ray. The remote radius is distance
along the negative horizontal ray from its equilibrium. Fable's default
away-from-the-other-equilibrium ray gives different section radii for the
same cycles; its coarse brackets must not be compared digit-for-digit with
horizontal radii. The remote equilibrium was solved again in binary128 before
translation. Tightening shooting tolerance from 2e-28 to 2e-30 preserves the
reported stability and digits. The full compactified engine also reports
U/S/U at the origin and S remotely for this exact coefficient vector.
See `verified_precursor.json` and the `rational_precursor` field record.
This is numerical reproduction, not an interval certificate.

## What was executed

The branch combines Fable's `4cece20` record with the previous main fold strike
`7db8597`. No random/Sobol sweep, descent, or section-theorem rederivation was run.
The new continuation uses the archived augmented two-half equations F=G=0,
with binary128 correction. The negative branch switches to actual coefficient
m before the (c,K) chart pole. An overlarge negative-K predictor was recorded
as unresolved and recovered from an archived nearby seed; it was not interpreted
as an endpoint. The positive branch is followed using log radius after its
rapid radial growth makes fixed-K prediction inappropriate. The first large-radius continuation was
interrupted after its counter failure so the controller and beta scale could
be repaired. `finish_positive.py` contains the final endpoint recovery.

There are **64 accepted fold points**,
**266 appended full-counter field records**, and
**270 distinct named fields** including controls and precision
follow-ups. These are deterministic continuation/unfolding fields. Every newly
accepted fold invokes the inherited `sweep_log.evaluate(coef)`. Each accepted
point also has a pair-present field and at least one beta level (usually two). Additional
small-beta fields preserve the large positive-sheet pair when the initial
levels are too large. Every such field has a rational 12-vector.

The initial pair offset is
`delta_c = -0.02 * exp(M_backward) * G_z / F_c`, evaluated in the appropriate
(K or m) coefficient chart. Positive-K pairs have S/U stability and negative-K
pairs U/S. Reversing the weak-focus stability with beta creates an innermost
U cycle on the positive sheet and S on the negative sheet, giving U/S/U and
S/U/S respectively. The origin itself then has the opposite stability from
its new innermost cycle.

For the large positive folds, beta is additionally scaled using the derivative
of the **bounded matching residual**, with a perturbation one tenth of the
estimated pair-destruction scale. Using a fixed small beta or a beta based on
a numerically ill-conditioned full displacement can destroy the outer pair.
The supplemental grid extends inward when the Hopf radius becomes very small.
`preserved_hopf.jsonl` records these controls and exact vectors.

## The compactified counter is not an exhaustive instrument on these sheets

Four concrete limitations were reproduced:

1. **False positive at large radius.** The nominal fold with log(horizontal
   radius)=20 produced six apparent sign changes in the denser double-precision
   full-return profile. The automatic four-root trigger stopped that continuation.
   Binary128 shooting found positive matching residuals at **both endpoints of
   all six suspect brackets**, and the full binary128 matching profile rejected
   the reported extra crossings. For example, at log radius 19.58 the
   double-precision profile gave D approximately -1.4047e-4, whereas the
   well-conditioned matching residual was +5.29655e-14. Their magnitudes use
   different coordinates, but their signs must agree whenever that full return
   exists. This is an integration/conditioning failure, not a candidate.
   `FOUR_ORIGIN_TRIGGER.json` is retained as a **rejected raw trigger**, with
   `trigger_binary128_recheck.json` and `trigger_resolution.json` explaining it.
2. **False negatives near the center and at large m.** The baseline counter can
   report zero even though binary128 resolves the pair and its Hopf completion.
   Its 1e-10 sign threshold and original-coordinate conditioning do not cover
   those regimes. A zero returned by that counter is a lower bound, not a
   nonexistence statement.
3. **Discarding a small endpoint misses an ordinary crossing.** Five apparent
   outside/edge sign disagreements arose because one endpoint of an ordinary
   root bracket fell below the profile's noise threshold. Binary128 restored
   U/S/U or S/U/S in all five fields. None supplied a fourth origin cycle.
   See `sign_disagreement_checks.log` and `precision_repairs.jsonl`.

4. **Coefficient resolution.** On the far-positive sheet, distinct exact
   nominal-fold and pair vectors can round to the identical float64 coefficient
   array. Even an otherwise exact integrator receiving that array cannot
   distinguish those fields. `float64_field_collapses.json` records these
   collisions. The binary128 matcher reads fraction strings directly; the
   baseline full counter necessarily evaluates a float64 approximation.

The untouched Fable evaluations remain in every field record. Supplemental
profiles scan a whole configured log-radius interval; they do not merely
track the inherited roots. The double grid uses spacing 0.1 plus fold samples;
the binary128 grid uses spacing 0.5 or 1 plus denser fold samples and root
refinement. Exact settings and all sampled signs are stored per profile.

The first two-state matching controller itself hit a sensitivity step guard
at the far-positive endpoint. The final controller retains position tolerance
and uses a separate 1e-24 floor for the sensitivity tolerance. Final endpoints
were re-polished in the actual-m chart; their pair and Hopf profiles were
recomputed with that controller. Tolerance controls and the initial failed
profiles are retained, rather than interpreting guard failures as zero cycles.

The binary128 matching implementation uses y=sqrt(m)*Y and positive time
rescaling to condition the large-m regime. Its origin values were compared
against the archived nine-state shooting engine. Matching zeros with valid
half passages close periodic orbits. The identity relating matching signs to
full displacement applies where the composed full return exists; two half
passages alone do not certify that a full return exists at every nonroot.
An angular-chart failure is recorded as unresolved. These finite grids still
cannot exclude unsampled close pairs or cycles outside the selected chart.

## Why sign agreement cannot close D1

If displacement is continuous on an actual returning interval, opposite signs
at its ends imply an odd number of crossings. Equal signs allow an even number,
including two new simple roots. For example, a local displacement
`D(u;t) = (u-a)^2 - t` creates a pair as t crosses zero while its two exterior
sample signs can remain positive throughout. Thus agreement outside the known
outer root and at an edge does **not** exclude another pair between them.
A sign disagreement is a useful missing-root detector, not a necessary
condition for pair birth. No monotonicity or bound on interior extrema was
proved here that would turn this diagnostic into a component obstruction.

The ledger distinguishes `scan_cap`, `integration_failure`, and
`angular_chart_failure`. `redge` from the inherited counter is its last sampled
successful radius, including its limited failure bisection; it is not necessarily
a dynamical boundary. The full displacement at a failed return is undefined.
Calling these values exact return-domain edges would overstate the computation.

## Remaining mathematical gap

There is no verified fourth origin cycle, no five-cycle field, no rigorous
component-wide upper bound, and no interval isolation of all roots or all
return-domain boundaries. The unsampled parameter intervals, both limiting
organizers, continuation beyond these endpoints, possible additional branches,
and possible extra pairs between grid samples remain open. The observed
2 -> 3 Hopf completion does not establish a universal ceiling of three.

The computational work requested here is recorded and reproducible; the proposed
D1 closure implication is not valid. A further decisive step needs an actual
additional root or a proved/enclosed restriction on the displacement's extrema.

## Files and replay

- `events_*.jsonl`: accepted/rejected corrections, exact nominal-fold/pair vectors,
  and labels of the associated full-count fields.
- `shooting.jsonl`: every new augmented-shooting call, including failures.
- `fields.jsonl`: unchanged baseline counter results, rational vectors, complete
  sampled profiles, root brackets, stability, and numerical edge diagnostics.
- `precision_repairs.jsonl`, `preserved_hopf.jsonl`: precision repairs and
  pair-preserving beta completion.
- `root_ledger.csv`, `sign_map.csv`, `summary.json`: generated review tables.
- `verified_precursor.json`: exact rational 3+1 field and binary128 reproduction.
- `README.md`: build and replay commands; `MANIFEST.json`: source/data hashes.

All results in this report beyond the elementary sign argument and exact
coefficient identities are numerical. Nothing is labeled an interval certificate.
