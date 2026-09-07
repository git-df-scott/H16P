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
