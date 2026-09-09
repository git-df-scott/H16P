# H16P working report

**Target:** one real planar polynomial vector field of total degree at most
two with at least five distinct isolated periodic orbits.

**Status: not found.** The strongest result obtained is a rigorous proof of
**three** distinct periodic orbits in one field, which is below the target and
below the published record of four. Everything else is floating-point
evidence. Nothing here is a counterexample.

Branch `claude/jolly-shannon-us5b4b`, draft PR #6. Base commit `498b58f`.

---

## 1. Provenance: what was read, and how

- **Read literally.** `H16P_NEXT_STRIKE.md`; `reversible_reseed/verify_control.py`
  and `README.md`; `research_2026_09_08/outputs/H16P-focus-route-correction.md`,
  `H16P-center-locus-compatibility.md`, `H16P-extreme-boundary-search.md`;
  openings of `STAGED_SHI_2026_09_05.md` and `FOUR_CYCLE_FRONTIER.md`.
- **Machine-inspected.** `reversible_reseed/data/verified_control.json`
  (replayed bit-identically); every JSON under `outer_fold_2026_09_09/data/`,
  `certify/data/`, `five_cycle_2026_09_09/data/`; every `.py` I wrote (all
  executed).
- **Listed but not read.** The bulk of the repository: 300 Python files in
  total, of which I wrote 15. Unread trees include `research_2026_09_08/`
  (151 py, 38 md, 180 json), `audit/` (42 py), `fold_surface_2026_09_05/`
  (34 py), `q4/` (22 py), `staged_2026_09_05/` (13 py), `kkl/` (8 py),
  `fold_closure_2026_09_05/` (6 py), plus 52 root-level markdown files of
  which I read four. **No claim below depends on unread material**, and where
  earlier documents in those trees conflict with this one, they were written
  before it and I have not reconciled them.
- **Papers.** Only abstracts, never proofs. Stated explicitly in §8.

---

## 2. The field

Seed family, five coefficients `theta = (a, b, e0, e1, e2)`:

```
x' = (b-2)/4 + (1-b) y + a x^2 + b y^2 + e1 x + e2 x y
y' = e0 - 2 x y
```

Seed point: `a = -7/4`, `b = 1/3`, `tau = 1e-4`, `e0 = -4 tau/11`,
`e1 = -31379 tau/25000`, `e2 = -7517 tau/5000`.

Two finite equilibria: `(-3.6e-05, 0.5)` and `(7e-06, -2.5)`.

Section `{x = 0}`, coordinate `s = log|y|`, `side = sign(y)`. The displacement

```
D(s) = s_forward(s) - s_backward(s)
```

is the mismatch between the forward and backward half returns. `D(s) = 0`
means the two half orbits close up: a periodic orbit.

---

## 3. What is proved

`certify/prove_cycles.py`. **At least three distinct periodic orbits of the
seed field**, by validated interval integration plus the intermediate value
theorem.

| bracket | half-width | `D` at the two ends | enclosure widths |
|---|---|---|---|
| upper inner, `s0 = 0.9802426972516128` | `1e-6` | positive / negative | `5.0e-29` |
| upper middle, `s0 = 1.3899487822846022` | `1e-6` | negative / positive | `6.5e-26`, `9.1e-26` |
| upper outer, `s0 = 2.0452069404076005` | `3e-7` | positive / negative | `7.6e-21`, `7.5e-21` |

The argument per bracket `I = [s0 - d, s0 + d]`:

1. **Segment validation.** Both half returns are computed with `I` itself as
   interval initial data. Success proves that for *every* `s` in `I` the
   forward and backward returns exist, leave transversally, reach the section
   again, and cross transversally. Hence `D` is defined on all of `I`, and
   continuous there by smooth dependence plus transversality.
2. **Endpoint signs.** `D(s0 - d)` and `D(s0 + d)` are enclosed from thin data;
   both exclude zero and have opposite signs.
3. IVT gives `s*` in the open interval with `D(s*) = 0`.
4. **Distinctness.** Each validated half return stops at the *first* crossing,
   so each closed orbit meets the upper section exactly twice: once in its
   bracket, once at its partner. The partner enclosures are
   `-1.4558814768662185`, `-1.5357456235202572`, `-1.6389947020748460`,
   pairwise disjoint and disjoint from all three brackets.

**Explicitly not proved.**

- That any of the three is **isolated** — i.e. a limit cycle. That needs a
  rigorous enclosure of `dD/ds` across each bracket, hence validated
  variational equations, which are not implemented. The target asks for
  *isolated* orbits, so this gap is material.
- Anything about the **lower** bracket near `s = 9.1437`. At `|y| ~ 9355` the
  field has `x' ~ b y^2 ~ 3e7`, which drives the validated step below any
  usable floor. A validated logarithmic reformulation is needed; not written.
- Three is not five. Four is already published (Shi 1980; Chen and Wang 1979;
  Galias and Tucker's interval-arithmetic certification). **This proof adds no
  new mathematics** — its value is a working pipeline.

---

## 4. How the continuity gap was closed

This was the blocker and the one genuine technical turn in the work.

Propagating a bracket as interval initial data first **failed on the wrapping
effect**, measured directly. From a section segment of width `2e-5`
(y-width `1.5e-4`) the enclosure grows to x-width `4.2e-03` at `t = 0.4`,
`3.97e-02` at `t = 0.8`, `1.43e-01` at `t = 1.2`, `3.33e-01` at `t = 1.6`, and
overflows before the crossing at `t ~ 2.7`. Independently, the measured
amplification for a full revolution is about `267x`.

Subdivision does not rescue it at that rate: holding the final width under
`1e-3` needs initial segments below `5e-7`, about 4000 subintervals per
bracket per time direction, roughly 11 hours of compute per bracket.

**The way through was to shrink the bracket rather than subdivide it.** `|D|`
at the endpoints falls only *linearly* in the half-width, while the enclosure
error stays near `1e-20`. So a half-width around `1e-6` keeps the endpoint
signs rigorous by fifteen orders of margin while keeping the segment small
enough to propagate. The half-width is searched per bracket over a ladder
`1e-6, 3e-7, 1e-7, 3e-8, 1e-8`, taking the first that satisfies both
conditions.

### A fix that changed a result

A fat box can straddle the section for more than one step, because its members
cross at different times. The single-step crossing routine cannot handle that
soundly — its sign tracking assumes the step begins strictly on one side. The
routine now **enlarges** the step until one validated step carries the whole
box across (shrinking would make the spread worse), and refuses the return if
no enlarged step does.

Under the earlier, unsound code the outer bracket passed through that path. Its
first "validated" segment did not stand. It now passes honestly at half-width
`3e-7`.

---

## 5. Numerical findings (not proofs)

### 5.1 The seed reproduces

`reversible_reseed/verify_control.py` replays bit-identically: all six moment
rows, all 24 return differences, all four sign brackets.

Repaired refinement gives four roots with audited derivatives. Spreads are
across three finite-difference step sizes and two tolerances, and are
**empirical uncertainties, not error bounds**:

| side | `s` root | `dD/ds` | rel. spread | worst parameter rel. spread |
|---|---|---|---|---|
| upper | 0.9802426972516128 | -1.2481e-07 | 7.1e-02 | 1.5e-02 (`a`) |
| upper | 1.3899487822846022 | +8.5268e-08 | 9.6e-02 | 1.5e-02 (`a`) |
| upper | 2.0452069404076005 | -2.4629e-07 | 9.3e-02 | 1.5e-02 (`a`) |
| lower | 9.1437116433290630 | +1.0228e-04 | 1.1e-05 | ~1e-06 |

Worst position uncertainty `|D|/|dD/ds| = 5.86e-07` against minimum upper
separation `0.4097` — a ratio of `1.43e-06`. Log-coordinate cross-check gaps
are `~1e-12`; that formulation shares the DOP853 stepper and event machinery,
so it is a **coordinate cross-check, not independent verification**.

### 5.2 Two counting errors — the most important numerical findings

**(a) The upper section coordinate is two-to-one on cycles.** The equilibrium
at `(0, 1/2)` lies *on* the section, so a closed orbit around it meets the
section twice. Scanning inward from the tracked roots produced three
apparently new sign changes. They are the same three cycles:

| inner root | its partner crossing | known root |
|---|---|---|
| -1.4558815 | 0.9802426 | 0.9802427 |
| -1.5357456 | 1.3899484 | 1.3899488 |
| -1.6389947 | 2.0452069 | 2.0452069 |

**Seven sign changes, four cycles.** A sign change of `D` is not a distinct
periodic orbit. Candidates are now clustered by their unordered crossing pair
before anything is counted.

**(b) `D` changes sign across the equilibrium itself.** Near `s = log(1/2)`
there is a further sign change whose return period collapses to `6.28318` —
the linearised value `2 pi` — with the section point `1.6e-03` from the
equilibrium. An artefact. The old `EQ_GUARD = 1e-3` measured distance along
the orbit and missed it; a dedicated `SECTION_EQ_GUARD = 5e-2` on the section
point catches it.

Neither error corrupted the previously reported counts, which used only outer
crossings. But the counting *procedure* was unsafe.

### 5.3 Outer-region scans

At the seed, 60 of 60 sampled points in `s` in `[2.046, 11.0]` resolved, with
`D` negative and decreasing at every sampled point. Along a 50-step branch that
scaled `(e0, e1, e2)` to `13.5x` the seed with all four trackers holding, the
same held at every accepted step.

**Sampled monotonicity is not continuous monotonicity.** Nothing was evaluated
between grid points, the windows are finite, and a double zero need not change
sign. No absence statement follows, and the earlier text claiming one is
withdrawn.

---

## 6. Failures, false positives, and defects — full list

Every one of these was found in my own work, most of them after the results
had already been written up.

| # | Defect | How it showed | Status |
|---|---|---|---|
| 1 | `double_zero_found` only tested whether a sampled slope turned positive | Read in source | Replaced by a real 6- then 7-equation solve |
| 2 | Negative branch reversed its own compensation gradient | `delta = step * direction * THETA_SCALE` negates the gradient when `step < 0` | Scheme abandoned |
| 3 | Tracker required `0.35 x gap >= 1e-3`, forcing a stop at gap `0.0028571` | Reported as an "incipient saddle-node at 2.8e-3" | **Withdrawn**: a built-in cutoff, not a bifurcation |
| 4 | Newton accepted on the *pre-update* residual, never checked the returned point | Read in source | Bracket-preserving Newton with a final acceptance test |
| 5 | Only `dD/ds` was audited, though five parameter sensitivities drove prediction | — | All five audited |
| 6 | Identity check compared corrected against *predicted* ordering | — | Moot under anchoring; then guarded explicitly |
| 7 | Failed gates discarded their measurements | — | `ReturnFailure.as_record()` |
| 8 | A rank-2 `2x6` Jacobian was described as making the solve "well posed" | — | **Withdrawn** |
| 9 | Approaching root pairs described as a saddle-node and its unfolding | — | **Withdrawn** |
| 10 | **Fixed-anchor fold solve reported 15 converged folds out of 15 — all spurious** | See below | Normalisation constraint added; regression-tested |
| 11 | **Released anchor lost its identity**, moving `1.288775 -> 1.431203` past the pinned anchor at `1.3899488` in one accepted step | Two roots of one `D` cannot cross | Motion and gap guards; run redone |
| 12 | **Fat-box crossing was unsound**; the outer bracket's first "validated" segment went through it | Found by tightening | Step-enlargement; bracket re-proved at `3e-7` |
| 13 | Broad random sampling had failing positive controls | 0 of 576 samples reproduced even a three-cycle nest | Not resumed |
| 14 | Guarded continuation stalls against its own guard, consuming accepted steps on microscopic moves | `theta` frozen to 6 digits over steps 12-24 | Recorded; step-size policy defect |
| 15 | **The fold solver is attracted to existing roots.** 466 of 634 solves walked inward to the known outermost cycle | Newton on `(D, dD/ds)` reaches the `D = 0` manifold at a simple root before the second equation bites | Recorded; the outer region was therefore never searched |

### Detail on #10, the instructive one

Fixing all four anchors and solving `D(s_i; theta) = 0` gives a scaled `4x5`
Jacobian of **rank 4**, singular values `(7997.7, 4123.0, 25.03, 0.926)`,
condition on range `8.64e+03`, with kernel

```
(-7.3e-08, -9.8e-08, 0.8804, 0.3039, 0.3640)
```

— no measurable `(a, b)` component. `D` is very nearly linear in
`(e0, e1, e2)`, so that kernel is essentially the ray that **shrinks the
perturbation**. Newton slid down it to `|(e0,e1,e2)| = 1.7e-10`, a factor
`1e-6` below the seed, where the field is numerically the unperturbed
reversible one and `D` vanishes identically — satisfying *both* fold equations
everywhere at once. Conditioning `1e+09` to `1e+10` and `|d2D/ds2| ~ 1e-10`
flagged it.

**A small residual is not a candidate.** Neither is a positive sampled slope,
nor a well-conditioned unconstrained Jacobian.

---

## 7. The experiments

### step1 / step2 — first pass (superseded)
Reproduced the seed, scanned outward. Its independence and absence claims were
wrong and are withdrawn (§5.3, §6).

### step3 — slope maximisation (wrong criterion)
`a` driven, `(b, e0, e1, e2)` compensating along the scaled gradient of
`sup dD/ds`. Branches: 11 accepted / 44 rejected and 10 / 44, both stopping at
the halving cap. Best `sup dD/ds` scaled: `-2.966e-01` and `-4.101e-03`. The
ascent flattened `D` around an adjacent pair of *existing* cycles — which, per
defect #3, approached only to the tracker's own cutoff.

### step4 — fixed anchors (false positives)
Rank 4, 15/15 spurious folds. See §6. Data retained in
`data/step4_anchored.json` as the record of the failure mode.

### step5 — normalised anchors with one released (current)
Stated before running:

- **Unknowns (7):** `theta`, released anchor `s_1`, fold coordinate `s*`.
- **Equations (7):** `D(s_1) = 0`; `D(s_i) = 0` for the three pinned anchors;
  `N(theta) = |(e0,e1,e2)|^2 / |(e0,e1,e2)_seed|^2 - 1 = 0`; `D(s*) = 0`;
  `dD/ds(s*) = 0`.
- Dropping the fold pair and `s*` leaves **5 equations in 6 unknowns**: one
  degree of freedom, followed by pseudo-arclength.

`N` removes the scaling degeneracy that produced #10 — the collapsed field is
no longer in the feasible set. Releasing an anchor is what supplies the freedom
to move `a` and `b` at all; under four fixed anchors that freedom does not
exist to first order.

Curve system: **rank 5**, singular values `(7997.7, 4123.0, 25.06, 4.98, 0.901)`,
condition on range `8.87e+03`, tangent dominated by the released anchor
(component `0.9974`).

Acceptance for a fold candidate requires **all** of: residual `< 1e-3` in the
fixed scales; `7x7` condition `<= 1e+06`; `|d2D/ds2| / DS_SCALE >= 1e-3`; every
anchor slope at least `0.1x` its seed value; and `s*` a **distinct** cycle from
all four anchors.

**Result (complete).** 172,980 return evaluations, 10,800 s.
**634 fold solves attempted, 0 converged, 0 passing the gates.**

| branch | accepted | rejected | stop | released anchor | curve condition |
|---|---|---|---|---|---|
| `curve_forward` | 24 | 40 | rejected cap | `0.980 -> 1.369948` | `1.96e+05` |
| `curve_backward` | 19 | 0 | **wall-time cap** | `0.980 -> 0.032454` | `3.15e+03` |

`pert_ratio` stayed exactly `1.0` throughout both, so the collapse mode of
defect #10 never recurred.

*Forward.* All 40 rejections are the identity guard: the released anchor moved
outward and stalled asymptotically at `s_free -> 1.369948`, precisely the
`0.02` guard below the pinned anchor at `1.3899488`. **That is the guard
firing, not a located merge** — the same mistake as defect #3, avoided here
only because the caveat was written before the result. It is also the wrong
direction: a merge of two existing cycles would *reduce* the count. The curve
conditioning degraded to `1.96e+05` as it stalled.

*Backward.* **Zero rejections.** The released anchor moved inward from
`0.980` to `0.032454` (`|y|: 2.665 -> 1.033`), with `a: -1.75 -> -1.76472` and
`b: 0.3333 -> 0.388666` — a far larger coefficient excursion than anything the
earlier schemes achieved — and the curve conditioning stayed healthy at
`3.15e+03`. **It was still running cleanly when the wall-time cap ended it.**
This branch is unfinished, not exhausted, and is the one to resume.

### The fold solver has a nameable failure mode

The 634 outcomes:

| count | status |
|---|---|
| 466 | iterate entered the excluded region |
| 151 | iteration cap without convergence |
| 17 | return failure during solve |

The dominant mode is the iterate walking **inward**: started across the outer
window, `s*` is consistently pulled back toward `1.88`-`2.06`, i.e. toward the
existing outermost root at `2.0452`, where the exclusion margin correctly
rejects it.

This is a defect in the formulation, and it should have been anticipated.
A damped Newton on `(D, dD/ds)` is attracted to the `D = 0` manifold, and the
nearest point of that manifold is an *existing simple root*. The second
equation is supposed to exclude it, but not before the trust-region steps have
already walked there. So the search never effectively explored outward at all.
**The correct reading is that this run did not test the outer region**, not
that the outer region is empty. A usable third attempt needs deflation —
dividing out the known roots — or solving `dD/ds = 0` first and only then
checking `D`.

### five_cycle probe — null, and its own controls failed
A general-quadratic nest counter, validated against the Chen-Wang
visualisation field: three cycles about the origin, outermost `r = 0.33380`
against `0.332839` recorded in `STAGED_SHI_2026_09_05.md`; the two inner radii
differ from that document's and **that discrepancy is unexplained**. The
bounded search evaluated 576 samples in 1658 s and found nothing — and not one
sample reproduced even a three-cycle nest, so the sampler left the thin region
immediately. **A verdict on the sampler; it supports no conclusion about
quadratic fields.** Not resumed.

---

## 8. Literature

- **Huang and Reyn (1995)**, "On the limit cycle distribution over two nests in
  quadratic systems", *Bull. Austral. Math. Soc.* **52**, 461-474,
  [DOI 10.1017/S0004972700014945](https://doi.org/10.1017/S0004972700014945).
  In a quadratic system with two nests, one nest contains exactly one cycle —
  **but the abstract restricts to systems where the sum of the multiplicities
  of the finite critical points equals three.** The seed field has two finite
  equilibria, so the hypothesis is not verified for it and the result must not
  constrain this family. An earlier version of my notes attributed an
  unrestricted form of this to Zhang Pingguang; **that was wrong and is
  withdrawn.** I read the abstract, not the proof.
- **Li Chengzhi (1986)**: no limit cycle surrounds an **exact** third-order
  weak focus of a real quadratic system. This implies **no** positive
  parameter-distance exclusion around that stratum, and earlier text of mine
  that leaned on it that way is corrected. Taken from the repository's own
  `H16P-focus-route-correction.md`; not independently audited.
- **Bautin (1954)**: cyclicity of a focus in quadratic systems is 3. Not
  audited.
- **Shi (1980)**, **Chen and Wang (1979)**: four limit cycles, so `H(2) >= 4`.
  **Galias and Tucker**: interval-arithmetic certification of exactly four for
  Shi's field. Whether `H(2) = 4` is open.
- I earlier used "most experts believe `H(2) = 4`" as a reason to discount the
  search. That is not evidence and should not have gated anything. Withdrawn.

---

## 9. Complete code inventory

3,398 lines of Python across 15 files, plus 4 READMEs and this report.

### `certify/` — validated integration (558 lines)

| file | lines | what it does |
|---|---|---|
| `taylor.py` | 156 | Interval arithmetic on `mpmath.iv` (outward rounded, arbitrary precision); tight interval square; Taylor coefficient recurrences for quadratic fields; Picard-Lindelof rough enclosure (`z0 + [0,h] f(B)` contained in `B`) proving existence on the step; order-K step with remainder bounded by re-running the recurrence over `B` |
| `poincare.py` | 168 | Validated half returns for the seed family on `{x = 0}`; rigorous crossing localisation (checks `x'` holds one sign across the step, then bisects the crossing time with interval evaluations); adaptive step with a floor; step **enlargement** when a fat box straddles the section |
| `certify_seed.py` | 88 | Validated `D` enclosures at bracket endpoints; writes `data/certify_seed.json` |
| `prove_cycles.py` | 146 | The existence proof: per-bracket half-width search, segment validation, endpoint signs, distinctness; writes `data/prove_cycles.json` |

Measured: harmonic oscillator one period, thin data, width `4e-36` (315 steps,
dps 40); Chen-Wang one revolution `7.6e-18`; seed half return `3.8e-21` in
`4.4 s` (673 steps); interval initial data amplifies about `267x` per
revolution.

### `outer_fold_2026_09_09/` — return-map machinery and continuation (2,529 lines)

| file | lines | what it does |
|---|---|---|
| `machinery.py` | 425 | The repaired evaluator. Itinerary gates (launch and terminal transversality, dense-output side preservation, interior crossings, orbit and section equilibrium distance) with structured failure evidence; event-corrected variational derivatives in `s` and all five coefficients; `derivative_audit` and `parameter_sensitivity_audit`; distinctness machinery (`partner_intersection`, `cycle_key`, `distinct_cycles`, `candidate_gate`); log-coordinate cross-check; fixed residual scales |
| `anchored.py` | 273 | `Budget` (evaluation and wall-time caps); bracket-preserving Newton with a final acceptance check and position uncertainty; the anchored `4x5` system with SVD diagnostics; pseudo-arclength corrector; the 6-equation fold system and solver |
| `step5_normalized.py` | 452 | The current experiment: normalisation constraint, released anchor, curve continuation, 7-equation fold solve, and the non-degeneracy acceptance gates |
| `step3_compensated.py` | 379 | Superseded slope-maximisation run (wrong criterion, kept as record) |
| `step4_anchored.py` | 315 | Superseded fixed-anchor run (the 15 false folds, kept as record) |
| `step2_outer_scan.py` | 245 | First-pass outward scan and continuation (superseded) |
| `step1_seed.py` | 225 | First-pass seed reproduction (superseded) |
| `controls_repairs.py` | 215 | 12 focused controls, all passing |

### `five_cycle_2026_09_09/` — general-quadratic probe (311 lines)

| file | lines | what it does |
|---|---|---|
| `counter.py` | 147 | Nest-aware cycle counting for a general quadratic chart: equilibria by resultant, foci by eigenvalues, ray section per focus with crossing direction taken from the rotation sense |
| `search.py` | 132 | The bounded parameter probe (null result) |
| `validate.py` | 32 | Validation against the Chen-Wang field |

### Controls

`controls_repairs.py`, **12/12 passing**: a deliberately failed return carrying
measurements; an unresolved return not read as a zero; `e0 = 0` shown to be the
invariant-line class rather than a centre (`D = 4.6e-03` there, so "centre
limit" was the wrong justification); the seven-sign-changes-are-four-cycles
test; the equilibrium artefact rejected by the section guard; a bad guess
rejected; position uncertainty against separation; a deliberate tracker failure
labelled as such; a positive sampled slope *not* accepted as a fold; all five
sensitivities audited; the anchored rank claim stated with its scales; and the
**collapsed-field regression** that would have accepted defect #10.

### Data

`outer_fold_2026_09_09/data/`: `step1_seed.json`, `step2_outer_scan.json`,
`outer_extension.json`, `step3_compensated.json` (6,195 lines),
`step4_anchored.json` (3,659 lines, the false-fold record),
`step5_firstrun_defect.log` (the identity-loss record), `controls_repairs.json`.
`certify/data/`: `certify_seed.json`, `prove_cycles.json`.
`five_cycle_2026_09_09/data/`: `validate.json`, `search.json`.

Failed runs and rejected steps are retained deliberately, with their reasons.

---

## 10. Reproduction

```bash
pip install numpy scipy sympy mpmath

python3 reversible_reseed/verify_control.py            # replays bit-identically
python3 outer_fold_2026_09_09/controls_repairs.py      # 12/12
python3 certify/prove_cycles.py                        # the 3-orbit proof, ~134 s
python3 certify/certify_seed.py                        # endpoint enclosures
python3 outer_fold_2026_09_09/step5_normalized.py       # the experiment (hours)
python3 five_cycle_2026_09_09/validate.py               # counter vs Chen-Wang
```

`step5_normalized.py` honours `MAX_ACCEPTED`, `MAX_EVALUATIONS`, `MAX_SECONDS`.

---

## 11. Caps and budget

No account usage meter is readable from this session, so no baseline or
percentage can be reported honestly. Explicit computational caps are imposed
instead and recorded in every output: 50 accepted and 40 rejected steps per
branch, plus separate return-evaluation and wall-time caps enforced by the
`Budget` object. Actual use: step3 4,769 evaluations / 296 s; step5 forward
branch and backward branch together 172,980 evaluations / 10,800 s; the proof 134 s; the abandoned
random search 1,658 s.

---

## 12. Open, in priority order

1. **Isolation.** Validated variational equations, then interval Newton on each
   bracket, to upgrade "three periodic orbits" to "three limit cycles". This is
   what the target actually asks for and is the clearest next build.
2. **The lower bracket.** A validated logarithmic reformulation, to bring the
   fourth orbit into the same rigorous framework and match the published
   four-cycle record with my own machinery.
3. **Resume the backward branch.** It took 19 accepted steps with **zero**
   rejections, moved the coefficients further than any earlier scheme
   (`a -> -1.76472`, `b -> 0.388666`), kept conditioning at `3.15e+03`, and
   stopped only because the wall-time cap expired. It is unfinished, not
   exhausted.
4. **A fold solver that does not fall into known roots.** Deflation of the
   existing roots, or solving `dD/ds = 0` first and checking `D` after. As it
   stands, 466 of 634 solves walked back to the outermost existing cycle, so
   the outer region was never actually searched. This is the third fold
   criterion of mine to fail, after slope maximisation and global collapse;
   the next one needs its degeneracies identified *before* it runs.
5. A step-size policy that does not stall against its own guard (defect #14).
6. The unexplained discrepancy between my Chen-Wang inner radii and those in
   `STAGED_SHI_2026_09_05.md`.

**No counterexample has been found, and nothing here bears on whether one
exists.**
