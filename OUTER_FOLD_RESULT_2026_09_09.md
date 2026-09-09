# Anchored continuation for a fifth cycle: repairs, corrections, and results

**No counterexample. No fifth cycle. Everything here is floating-point
evidence: no interval arithmetic, no validated integration, no certified
periodic orbit.** Later sections supersede earlier claims wherever they
conflict.

## What was read, and how

Source commit reviewed: `1a6e8e7d5d0ca07de333cc59add7b5dc66ba721f` (the PR
head at the time of the review).

- **Read literally**: `H16P_NEXT_STRIKE.md`; the previous version of this
  file; `reversible_reseed/verify_control.py` and `README.md`;
  `research_2026_09_08/outputs/H16P-focus-route-correction.md`,
  `H16P-center-locus-compatibility.md`, `H16P-extreme-boundary-search.md`;
  `STAGED_SHI_2026_09_05.md` (opening); `FOUR_CYCLE_FRONTIER.md` (opening).
- **Machine-inspected**: `reversible_reseed/data/verified_control.json`
  (replayed, bit-identical); `outer_fold_2026_09_09/data/*.json`;
  `outer_fold_2026_09_09/*.py` (all executed).
- **Still pending**: the bulk of `research_2026_09_08/outputs/` (roughly 180
  files), the `q4/`, `kkl/`, `staged_2026_09_05/` and `council/` trees, and
  the full text of every cited paper. Claims below that depend on unread
  material are marked.

## The nine review findings: all nine confirmed

| # | Finding | Verification |
|---|---|---|
| 1 | No actual double-zero solve; `double_zero_found` only tested a sampled slope | Read in source. Replaced by a real solve. |
| 2 | Negative branch reverses the compensation gradient; no objective-improvement check | `delta = step * direction * THETA_SCALE` with `step < 0` negates the gradient component too. Confirmed by inspection. |
| 3 | Tracker requires `0.35 x gap >= 0.001`, forcing a stop near gap 0.002857 | `1e-3 / 0.35 = 0.0028571`. Reported "0.0028". **The stop was a built-in cutoff, not a located bifurcation.** |
| 4 | Newton accepted on a pre-update residual, no final acceptance check, no preserved bracket | Read in source: `converged` used `r['D']` from before the update; `final` was computed but never tested. |
| 5 | Only `dD/ds` was audited, not the five parameter sensitivities driving prediction | Confirmed; now audited. |
| 6 | Identity check compared corrected against *predicted* ordering | Confirmed. Made moot by the anchored scheme, which pins identities. |
| 7 | Failed gates discarded their measurements | Confirmed; `ReturnFailure` now carries them. |
| 8 | Rank-2 2x6 Jacobian does not establish a well-posed fold solve | Accepted. The claim is withdrawn. |
| 9 | Approaching root pairs do not establish a saddle-node or its unfolding | Accepted, and finding 3 makes the reported approach partly an artefact. |

### Corrections to earlier text, which no longer stands

- The "incipient saddle-node at separation 2.8e-3" was **at the tracker's own
  cutoff**. It is not evidence of a bifurcation. Withdrawn.
- "The fold isn't there", "a cycle was lost", and "coordinated five-parameter
  motion is necessary" are all withdrawn; none was established.
- The rank-2 Jacobian was described as making the solve "well posed". It does
  not. Withdrawn.

## Two correctness findings that change how cycles are counted

These emerged while testing the repairs and matter more than anything else
here.

**1. The upper section coordinate is two-to-one on cycles.** The equilibrium
at `(0, 1/2)` lies *on* the section `{x = 0}`, so a closed orbit around it
meets the section twice. Scanning inward from the tracked roots produced what
looked like three additional sign changes at `s = -1.4559, -1.5357, -1.6390`.
They are not new cycles: each one's partner crossing is exactly a known root.

| inner root `s` | partner crossing | known root |
|---|---|---|
| -1.4558815 | 0.9802426 | 0.9802427 |
| -1.5357456 | 1.3899484 | 1.3899488 |
| -1.6389947 | 2.0452069 | 2.0452069 |

Seven sign changes, four cycles. **A sign change of `D` is not a distinct
periodic orbit**, and `distinct_cycles()` now clusters candidates by their
unordered crossing pair before anything is counted.

**2. `D` changes sign across the equilibrium itself.** Near `s = log(1/2)`
there is a further sign change whose return period collapses to `6.28318` —
the linearised value `2*pi` — with the section point `1.6e-03` from the
equilibrium. It is an equilibrium artefact. The old `EQ_GUARD = 1e-3` measured
distance along the orbit and did not catch it; a dedicated
`SECTION_EQ_GUARD = 5e-2` on the section point does.

Nothing in the earlier reported counts is invalidated by this — they used only
the outer crossings — but the counting procedure was not safe, and now is.

## Repairs

`machinery.py`: structured failure evidence (`ReturnFailure.as_record()`);
launch/terminal transversality, dense-output side preservation, interior
crossings, orbit and section equilibrium distance; event-corrected variational
derivatives in `s` and in all five coefficients;
`parameter_sensitivity_audit()` covering all five; distinctness machinery.

`anchored.py`: `refine_root()` keeps a sign bracket throughout, bisects when
Newton would leave it, and **re-evaluates at the returned point before
accepting**, rejecting a final residual above `5e-14`. It reports a position
uncertainty `|D| / |dD/ds|`. On the seed the worst is `5.86e-07` against a
minimum upper separation of `0.4097` — a ratio of `1.43e-06`.

Audited seed quantities (spreads are empirical, not error bounds):

| Side | `s` root | `dD/ds` | `dD/ds` rel. spread | worst parameter rel. spread |
|---|---|---|---|---|
| upper | 0.9802427 | -1.2481e-07 | 7.1e-02 | 1.5e-02 (`a`) |
| upper | 1.3899488 | +8.5268e-08 | 9.6e-02 | 1.5e-02 (`a`) |
| upper | 2.0452069 | -2.4629e-07 | 9.3e-02 | 1.5e-02 (`a`) |
| lower | 9.1437116 | +1.0228e-04 | 1.1e-05 | ~1e-06 |

The Cartesian and logarithmic formulations share the DOP853 stepper and the
same event machinery: a **coordinate cross-check**, not independent
verification. Cross-check gaps are `~1e-12`.

**Controls** (`controls_repairs.py`, 12/12 passing): a deliberately failed
return carrying measurements; an unresolved return not read as a zero; a bad
guess rejected; position uncertainty against separation; a deliberate tracker
failure labelled as such; a positive sampled slope *not* accepted as a fold;
all five sensitivities audited; the two-to-one counting test; the equilibrium
artefact rejected; and the false-fold regression below.

## The experiment, and a false positive worth recording

### step4: fixed anchors — and 15 spurious folds

Fixing all four anchors `s_i` and solving `D(s_i; theta) = 0` gives a scaled
`4x5` Jacobian of **rank 4**, singular values
`(7997.7, 4123.0, 25.03, 0.926)`, condition on its range `8.64e+03`. Its
kernel is

    (-7.3e-08, -9.8e-08, 0.8804, 0.3039, 0.3640)

— no measurable `(a, b)` component. Following that kernel, the six-equation
fold solve reported **converged 15 times out of 15**, at every multi-start.

All fifteen are spurious. `D` is very nearly linear in `(e0, e1, e2)`, so the
kernel is essentially the ray that *shrinks the perturbation*. Newton slid
down it to `|(e0,e1,e2)| = 1.7e-10`, a factor `1e-6` below the seed, where the
field is numerically the unperturbed reversible one and `D` vanishes
identically — satisfying both fold equations everywhere at once. Conditioning
`1e+09`–`1e+10` and `|d2D/ds2| ~ 1e-10` flagged it. The raw data is kept in
`data/step4_anchored.json` as the record of the failure mode, and a control
now regression-tests that this field is rejected.

**This is why a small residual is not a candidate.**

### step5: normalisation and a released anchor

Stated before running:

- **Unknowns (7)**: `theta = (a, b, e0, e1, e2)`, the released anchor position
  `s_1`, the fold coordinate `s*`.
- **Equations (7)**: `D(s_1) = 0`; `D(s_i) = 0` for the three pinned anchors;
  `N(theta) = |(e0,e1,e2)|^2 / |(e0,e1,e2)_seed|^2 - 1 = 0`; `D(s*) = 0`;
  `dD/ds(s*) = 0`.
- Dropping the two fold equations and `s*` leaves **5 equations in 6
  unknowns**: one degree of freedom, followed by pseudo-arclength.

`N` removes the scaling degeneracy: the collapsed field is no longer in the
feasible set. The released anchor is the innermost upper root; releasing it is
what supplies the freedom to move `a` and `b` at all — under four fixed
anchors that freedom does not exist to first order.

The curve system has **rank 5**, singular values
`(7997.7, 4123.0, 25.06, 4.98, 0.901)`, condition on range `8.87e+03`, and
tangent dominated by the released anchor (component `0.9974`).

Acceptance thresholds, fixed and documented, all of which must hold:
residual `< 1e-3` in the fixed scales; `7x7` condition `<= 1e+06`;
`|d2D/ds2| / DS_SCALE >= 1e-3`; every anchor slope at least `0.1x` its seed
value; and `s*` a **distinct** cycle from all four anchors by the crossing-pair
test. Caps: 50 accepted and 40 rejected steps per branch, plus separate
return-evaluation and wall-time caps.

Early steps confirm the repair: `pert_ratio` stays exactly `1.0`, `a` and `b`
move (`a: -1.75 -> -1.74757`, `b: 0.3333 -> 0.3262` by step 2), and **no fold
solve converges** — the false positives are gone. Results are in
`data/step5_normalized.json`.

### A third defect, found in the running experiment

The first step5 run moved the released anchor from `s = 1.288775` to
`s = 1.431203` in a single accepted step, passing through the pinned anchor at
`1.3899488`. Two roots of the same `D(.; theta)` cannot cross — they can only
merge and separate — so that step was the corrector jumping to a different
root, not the root moving. Every step after it is unsound. The run is kept as
`data/step5_normalized_firstrun.json` and `data/step5_firstrun_defect.log`,
and the re-run adds two guards: the released anchor may move at most `0.05`
per accepted step, and must stay at least `0.02` from every pinned anchor.

Under those guards the released anchor advances in controlled increments
(`0.980 -> 1.030 -> 1.080 -> 1.130`), moving outward toward the middle anchor,
with `pert_ratio` pinned at `1.0` and no fold solve converging. If the
continuation stops when the gap reaches `0.02`, that is the **identity guard
firing**, and must not be reported as a located merge — that is exactly the
mistake finding 3 identified in the previous run.

## Validated integration

`certify/` now contains a working validated integrator: interval arithmetic on
`mpmath.iv`, a Picard-Lindelof rough enclosure, an order-K Taylor step whose
remainder is bounded over that enclosure, and a validated return map for this
family with rigorous crossing localisation.

Measured: a full harmonic-oscillator period encloses to width `4e-36`; a
Chen-Wang revolution to `7.6e-18`; a seed half return to `3.8e-21` in `4.4 s`;
interval initial data amplifies about `267x` per revolution.

At `s0 +- 1e-3` around each of the three upper candidate roots, the enclosure
of `D` is **strictly signed, with opposite signs at the two endpoints**:

| bracket | `D(s0 - 1e-3)` | `D(s0 + 1e-3)` | enclosure width |
|---|---|---|---|
| upper inner | positive | negative | `4.9e-29` |
| upper middle | negative | positive | `1.7e-25` |
| upper outer | positive | negative | `7.8e-21` |

Those six signs are rigorous. **They are not an existence proof.** Concluding
a root lies between each pair needs `D` defined and continuous across the
whole bracket, i.e. the return map validated for every initial condition in
it, not just at its ends. That is not done. The lower bracket is not attempted
at all: `x' ~ b y^2 ~ 3e7` at `|y| ~ 9355` drives the validated step below any
usable floor, and needs a logarithmic reformulation that is not implemented.

## Literature corrections

- The distribution result previously cited here as an unrestricted theorem of
  Zhang Pingguang is **Huang and Reyn, "On the limit cycle distribution over
  two nests in quadratic systems", Bull. Austral. Math. Soc. 52 (1995)
  461-474**, [DOI 10.1017/S0004972700014945](https://doi.org/10.1017/S0004972700014945).
  Its abstract restricts to systems where *the sum of the multiplicities of
  the finite critical points equals three*. The seed field here has two finite
  equilibria, so that hypothesis is not verified for it, and the result must
  not be used to constrain this family. The earlier "(1, i) distribution, so
  five cycles need four in one nest" claim is withdrawn. I read the abstract,
  not the proof.
- Li Chengzhi's result concerns an **exact** third-order weak focus. It does
  not imply any positive parameter-distance exclusion around that stratum, and
  the earlier text that leaned on it that way is corrected.
- The broad random sampling in `five_cycle_2026_09_09/` had failing positive
  controls (not one of 576 samples reproduced a three-cycle nest). Per the
  research discipline, it is **not resumed**.

## Budget

No account usage meter is readable from this session, so no baseline or
percentage can be reported honestly. Explicit computational caps are imposed
instead and recorded in the data: 50 accepted steps and 40 rejected steps per
branch, a return-evaluation cap, and a wall-time cap, each enforced by the
`Budget` object and reported in the output.

## Status

**Passed**: control replays bit-identically; four distinct cycles reproduced
under the repaired gates with audited derivatives and distinctness; the
repaired tracker follows directions the first pass called lost roots; 12/12
controls; the normalised scheme eliminates the collapse mode.

**Failed**: step3's slope criterion and step4's fold criterion both produced
false positives, now regression-tested.

**Unresolved**: whether a genuine double zero exists beyond the outermost
upper cycle in this family; whether any tracked root is an actual periodic
orbit in the rigorous sense. `certify/` contains a working validated Taylor
integrator (a full revolution of the Chen-Wang field enclosed to width
`7.6e-18`; measured interval wrapping factor `267` per revolution) but it is
not yet wired to these return maps.

**Most valuable next experiment**: complete the step5 curve in both
directions, and if no fold candidate passes the gates, wire the validated
integrator to the seed return map so that the four existing brackets can be
certified rather than merely sampled.
