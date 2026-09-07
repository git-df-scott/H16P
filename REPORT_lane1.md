# REPORT_lane1.md -- Andronov-Hopf curve sweep for a fourth cycle in one nest

Branch `fable/lane1-ahcurve`, branched from `fable/coordination-2026-09-06`.
No emojis anywhere.  No pull request.

## TWO SESSIONS ARE ON THIS BRANCH

A second Lane 1 session (`session_01UqdMSqDz9KHbgPtDkycpJq`) pushed
`9de8f4e` to this same branch with its own engine in `lane1/`, and its
validation in the top-level `VALIDATION.md`.  This session
(`session_01A6xiV4vzHJud4DWhSdYPGH`) works in `lane1_ahcurve/` with its
validation in `lane1_ahcurve/VALIDATION.md`.  The two trees do not overlap and
the merge was clean.  Auditor: read both validation files.

That collision is worth more than it costs.  The two engines were written
independently and locate the section crossing by different means -- `lane1/`
carries the winding angle as a third integrated state and Newton-solves the
step length, `lane1_ahcurve/` accumulates the signed angle across accepted
steps and bisects the last one.  They agree:

| quantity | lane1/ | lane1_ahcurve/ |
|---|---|---|
| Cherkas row 1 cycle x | 1.2809, 2.0070, 4.0193 | 1.28091, 1.00700+1, 3.01932+1 |
| Cherkas row 4 cycle x | 0.5569, 0.7466, 0.8523 | 0.55695, 0.74658, 0.85232 |
| Cherkas row 8 cycle x | 1.3573, 2.3071, 4.1455 | 1.35730, 2.30708, 3.14553+1 |
| KKL origin nest r | 0.6832, 2.1837, 15.9628 | 0.68321, 2.18370, 15.96278 |
| KKL remote cycle from B | 3706 | 3706.05 |
| interior extrema of beta*, all 8 Cherkas rows | 2 | 2 |

Every Cherkas abscissa matches to four decimals across two independent
implementations.  That is a stronger PROTOCOL rule 2 check than either engine
plus scipy on its own.

The rest of this report is this session's (`lane1_ahcurve/`) work.  Engine hash
is printed at the top of `lane1_ahcurve/VALIDATION.md` and stamped into every
ledger row.

## Checkpoint 1 -- engine built and validated

### What ran

1. **Engine.** `lane1_ahcurve/ahcurve.c` (C, OpenMP) plus `engine.py` (ctypes).
   Dormand-Prince 5(4), field re-expanded about the focus so the focus is the
   origin, section a ray from the focus, return detected by the accumulated
   signed angle reaching `2*pi` with a bisection on the last step.
   `beta*(s)` is obtained by a monotone root solve of the family parameter to
   `ptol = 1e-10`, in two families: the uniform rotation
   `X_b = (P cos b - Q sin b, P sin b + Q cos b)` and a linear coefficient
   direction (used for the Cherkas rotating parameter `a11`).
   `REVIEW_engine.md` bugs A1, A2, A3, B1, B2 are addressed structurally, not
   patched; the table is in `lane1_ahcurve/VALIDATION.md`.
   Cost: about 5000 returns per `beta*` curve, ~0.3 s wall for a 300-point
   curve at two tolerances on four cores.

2. **Second integrator.** `refengine.py`: scipy DOP853 in *global* coordinates
   with a two-event half-turn crossing predicate. Different method, different
   coordinates, different section logic (PROTOCOL rule 2).

3. **Validation** (`validate.py`, `extraprobe.py`, rendered by
   `mkvalidation.py` into `lane1_ahcurve/VALIDATION.md`).  Result: **passed, the lane may
   sweep.**
   - all nine fat seeds and the KKL control give exactly **three** displacement
     sign changes in the primary nest, every bracket clearing the two-tolerance
     noise floor;
   - the same count on two further sections (rays rotated by +/- 0.7 rad):
     thirty counts, all three;
   - every bracket endpoint reproduced by scipy, largest disagreement 1.4e-10
     against endpoint magnitudes of 1e-6 and up;
   - `beta*(s)` has **exactly two interior extrema on all ten**;
   - Perko P3 and the KKL control reproduce the published cycle positions to
     5.3e-5 and 1.6e-5;
   - Cherkas row 4's `a11` curve reproduces the published degree-6
     Andronov-Hopf polynomial: two interior extrema at x = 0.62 and 0.80
     (published fit 0.625, 0.805), height range 5.04e-4 vs 5.28e-4, pointwise
     difference at most 5.1e-5.

### Ledger sizes

See `lane1_ahcurve/data/`.  At this checkpoint: `validation.jsonl` (12 rows,
one per seed plus the two extra probes), `validation_summary.json`,
`rounding_box.json`, `remote_probe.json`.

### Max extrema seen

**2.**  No field examined so far has a `beta*` with three interior extrema, and
no nest has shown four displacement sign changes.  No TRIGGER file written.

### Best candidates

None yet -- the sweep has not started.  The nine fat seeds are the starting
points; their `beta*` height ranges (the width of the rotation window in which
three cycles exist) are, smallest first:

| seed | beta* height range |
|---|---|
| cherkas8 | 4.03e-06 |
| cherkas6 | 9.35e-05 |
| cherkas2 | 1.23e-04 |
| cherkas7 | 1.63e-04 |
| cherkas1 | 2.67e-04 |
| kkl_control | 4.25e-04 |
| cherkas3 | 4.78e-04 |
| cherkas5 | 6.36e-04 |
| cherkas4 | 9.96e-04 |
| perko_p3 | 4.21e-03 |

Perko P3 has by far the widest three-cycle rotation window, followed by
Cherkas 4 and 5.  Those are the seeds where a third extremum has the most room
to appear without the existing two colliding, so they lead the sweep order.

### Open problems

1. **Cherkas rows 7 and 8 remote nests: UNRESOLVED, not "no cycle".**  At the
   tabulated coefficients this engine finds no cycle around the second focus B
   on any of twelve rays.  In the paper's own parameter, the `a11`
   Andronov-Hopf curve *of the remote nest* spans [0.334, 0.658] for row 7 and
   [1.870, 2.040] for row 8, while the table prints 2.1502 and 1.51997.  The
   gap is five orders of magnitude larger than the printed rounding, so unlike
   the primary-nest position residuals it is not a rounding effect.  Either the
   paper's remote cycle sits elsewhere in the coefficient box, or it lies
   outside this engine's return domain around B.  This does not gate the lane:
   the mandated "+1" control is KKL, which reproduces (3,1) exactly.
2. The `beta*` root solve returns UNRESOLVED near the outer end of a nest for
   some seeds (cherkas1: 210/300 points, cherkas6: 252/300) because the rotated
   family member loses its return before the displacement changes sign.  That
   truncates the interval on which a third extremum could be detected at the
   outer end.  Widening it needs a compactified chart, not more precision.
3. The engine is double precision.  Anything that would decide a trigger has to
   be re-run in long double / binary128 or mpmath dps 40 (PROTOCOL rule 6);
   that path is written but not yet exercised, because nothing has triggered.

### Next step

Build the sweep: raw 12-coefficient space with the focus translated to a fixed
point, one scale fixed, rotation integrated out by `beta*`.  Random
perturbations of the nine fat seeds at relative sizes 1e-3, 1e-2, 1e-1, then
hill-climb on (number of interior extrema of `beta*`, then the flatness of the
curve at inflection points outside the two-extremum interval).

---

# Session `session_01UqdMSqDz9KHbgPtDkycpJq` (`lane1/`) -- checkpoint 2: the sweep

Everything below concerns the `lane1/` tree.  Its validation is the top-level
`VALIDATION.md`.  **No field with more than two interior extrema of `beta*` has
been found.  No counterexample is claimed and none is implied.**

## What is in `lane1/`

| file | what it is |
|---|---|
| `retmap1.c` | the return-map engine, one source compiled twice: `double` and, with `-DLANE1_QUAD`, binary128 |
| `engine.py` | ctypes driver, focus-local re-expansion, equilibria, two-tolerance noise, `beta*` |
| `refengine.py` | independent SciPy DOP853 return map (different crossing algorithm) |
| `hiprec.py` | trigger-precision recheck: binary128 at two tolerances plus mpmath `odefun` at dps 45 |
| `ahsweep.py` | search space, curve features, objective, ledger |
| `sweep.py` | random-perturbation and (1+1)-ES drivers, trigger discipline |
| `cusp.py`, `cusprun.py` | Newton continuation onto the multiplicity-three cycle |
| `linescan.py` | one-dimensional scans of the extremum count along random lines |
| `validate.py`, `analyze.py`, `seeds.py` | validation, ledger analysis, seed fields |

## Precision ladder (PROTOCOL rules 2 and 6)

`retmap1.c` carries a precision macro layer, so the double and the binary128
engines are the *same algorithm* at two precisions, and `refengine.py` and
mpmath are two *different* algorithms.  On Cherkas row 4, `s = 0.2`, `b = 0`:

| engine | `D` |
|---|---|
| `retmap1.c` double, rtol 1e-12 | `+3.21481573188620e-05` |
| `retmap1.c` binary128, rtol 1e-18 | `+3.214815748046300e-05` |
| mpmath `odefun`, dps 42 | `+3.2148157480475600e-05` |

binary128 and mpmath agree to `1.3e-19`; the double engine's own error is
`1.6e-13`, which is what the two-tolerance noise rule estimates.  Trigger
brackets are re-decided in binary128 and cross-read by the SciPy engine on a
second section; nothing has triggered, so that path has been exercised only on
known points.

## The search space

A base field is the focus-local 10-vector
`L = (p1..p5, q1..q5)` of the field re-expanded about its focus, so the focus
is at the origin of the integration chart by construction (this is the fix for
`REVIEW_engine.md` bug A2 -- an offset chart manufactured ten spurious cycles
in the earlier engine).  Two directions are projected out of every
perturbation: `L` itself (a time rescaling) and `V_rot = (-q, +p)` (the uniform
rotation, which `beta*` already integrates out).  `||L||` is held fixed.  That
leaves **8 live directions**.

## Counting rule for extrema

The height range of `beta*` is dominated by its blow-up at the boundary of the
nest, not by the part of the curve where the cycles live -- on Cherkas row 4
the whole three-cycle wiggle occupies the top 10% of the range.  So an extremum
is counted by its prominence in **absolute** `b`, against the measured noise
floor of `beta*` rather than against that range.  Measured floor (rtol 1e-12 /
btol 1e-10 against rtol 1e-13 / btol 1e-12): median `6e-13`, worst case
`2.3e-11`.  Threshold `PROM_ABS = 1e-9`, about 50x the worst case; a
conservative count at `1e-7` is recorded alongside.  Both counts are 2 on every
Cherkas row, so no narrow third extremum is being discarded by the threshold.

The `b`-height actually spanned by the two extrema -- the window in which three
cycles exist -- is as narrow as the literature suggests:

| seed | `wiggle_range` (in b) | height range of `beta*` |
|---|---|---|
| cherkas1 | 1.04e-04 | 3.03e-04 |
| cherkas2 | 6.14e-06 | 3.56e-04 |
| cherkas3 | 3.47e-04 | 3.47e-04 |
| cherkas4 | 1.18e-04 | 1.23e-03 |
| cherkas5 | 2.37e-04 | 8.71e-04 |
| cherkas6 | 1.76e-05 | 1.08e-04 |
| cherkas7 | 4.21e-06 | 2.10e-03 |
| cherkas8 | **1.24e-07** | 3.48e-05 |

Row 8's three-cycle window is `1.2e-07` wide in the rotation angle and is still
1000x above the noise floor.

## Campaign 1 -- random perturbations and hill climbing

**11,745 fields evaluated** (ledgers `lane1/ledger/perturb_*.jsonl`,
`climb_*.jsonl`), 11,502 resolved, **0 with three or more interior extrema**.

The perturbation statistics are the substantive result.  Counting only
resolved fields, by relative perturbation size in the 8 live directions:

| relative size | 2 extrema (3 cycles) | 1 extremum | 0 extrema |
|---|---|---|---|
| 1e-3 | 235 | 115 | 550 |
| 1e-2 | 12 | 72 | 787 |
| 1e-1 | 0 | 19 | 683 |

**A relative displacement of 0.1% in the raw coefficients already destroys the
three-cycle configuration three times out of four, and 1% destroys it in 98% of
directions.**  The published examples sit at the centre of a region of relative
width of order 1e-3 in the 8 live directions.  That is a quantitative statement
of why forty years of random search found nothing, and it is the reason the
rest of this lane moves to solving rather than sampling.

### The hill climbs stalled -- recorded as a failure, not a result

18 (1+1)-ES runs of 500 iterations each.  Acceptance was 0-31 out of 500, and
in **every** run the step size collapsed to its floor (`1e-6`) and stayed
there, so the second half of each run proposed nothing.  The climbs therefore
do not license any statement about the landscape; only the perturbation
statistics and the line scans below do.  The fix is not a bigger budget -- it
is that the target is a codimension-one variety and a random walk should not be
asked to find one.

### Remote nests (recorded per the brief)

Four seeds have a second focus.  On each, `beta*` of the remote nest is
**monotone** (0 interior extrema), i.e. the remote nest carries at most one
cycle at any rotation angle -- the "free" cycle of Zegeling's Lemma 6.4, and
consistent with every published `(3,1)`.

| seed | second focus B | trace at B | remote `beta*` |
|---|---|---|---|
| cherkas2 | (-1.7309, +2.3682) | +8.088 | monotone |
| cherkas7 | (-4.2032, +1.3122) | -5.817 | monotone |
| cherkas8 | (-2.6690, +1.5991) | +0.8685 | monotone |
| perkoP3 | (0, 1) | +4.495 | monotone |

The KKL control's remote cycle *is* found, around its second focus
`B = (-6.259641, +7.449768)`, at `s = 3706` on the `-x` ray from B, crossing
the plane at `x = -3712` against the repository's recorded `-3711.56`.

## Campaign 2 -- Newton continuation onto the cusp (`cusp.py`)

With `u = log s` and `A(u, lambda) = beta*(s)` along a direction `v`, an
interior extremum of `beta*` is `A_u = 0` (a multiplicity-two cycle) and a
**cusp** -- PROTOCOL section (c), a multiplicity-three cycle -- is

        A_u = 0,  A_uu = 0,

two equations in two unknowns, solvable by Newton.  Derivatives come from a
degree-4 least-squares fit on a 9-point stencil in `u`, which is what makes
`A_uu` usable against `beta*`'s 1e-12 of noise.  Newton is confined to one of
three windows: `inner` (between the Andronov-Hopf plateau and the first
extremum), `outer` (beyond the last extremum, out to the nest boundary), and
`between` (a control -- the cusp there is where the seed's own two extrema
annihilate).  The confinement is essential: unconfined, Newton runs to the
plateau near the focus, where `A_u = A_uu = 0` holds for every field to within
noise; a nondegeneracy gate on `A_uuu` rejects those.

First result, Cherkas row 4, 8 random directions x 7 starts:

| window | outcome |
|---|---|
| `between` | **7 of 8 directions converge**, at `lambda` between `-4.5e-3` and `+7e-4`, at `s ~ 0.23` |
| `inner` | 20 of 24 hit the degeneracy gate (the plateau), 3 fail to converge |
| `outer` | 22 of 24 fail to converge within `abs(lambda) <= 0.8` |

**Along a random line, the cusp that destroys the existing pair of extrema sits
at a relative field displacement of about `5e-4`, while no cusp that would
create a new pair exists in the outer annulus out to a displacement of `0.8`.**
The three-cycle configuration is some three orders of magnitude closer to being
destroyed than to gaining a fourth cycle, in every direction tried so far.
This is a measurement along random lines from one seed; it is not a theorem and
does not exclude a cusp in a direction not yet sampled.

## Campaign 3 -- line scans (running)

`linescan.py` walks `lambda` across a 41-point geometric grid spanning
`+/- 0.6` for each random direction and counts interior extrema at every point,
so nothing is missed between samples and the far end of the segment -- where
the field no longer resembles the seed -- is seen too.  80 directions per seed
over all 9 seeds are in flight; results in the next checkpoint.
