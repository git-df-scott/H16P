# Outer-fold strike: repaired machinery, seed audit, and compensated continuation

**Status: executed and bounded. No counterexample. No fifth cycle. Everything
below is floating-point evidence, not certification.**

Execution of [H16P_NEXT_STRIKE.md](H16P_NEXT_STRIKE.md) Steps 1-3, plus
repairs to the evidence machinery that the first pass of this work got wrong.
Code and data: [outer_fold_2026_09_09/](outer_fold_2026_09_09/).

## Corrections to the first pass

The first version of this document and of the pull-request text overstated
several things. They are corrected here and the wording is not reused.

- **"No outer fold on the branches tried" / "the fold isn't there."** Wrong
  form of claim. Sampled profiles cannot establish absence. What the scans
  show is that at every sampled point of the outer windows examined, `dD/ds`
  was negative; nothing was sampled between grid points, and a double zero
  need not change sign, so no sampled profile can exclude one.
- **"`a_sweep` rejected at step 1, root lost."** Wrong attribution. A single
  unadapted step of `0.05` in `a` failed *the tracker*. It establishes
  nothing about whether a cycle continued to exist. With sensitivity-based
  prediction and adaptive step reduction (below), the `a` direction is
  followed successfully.
- **"A fold would need `a` moving with `b, e0, e1, e2` compensating."** Not
  established. That was one hypothesis, not a necessity, and it is stated as
  a chosen search direction below rather than as a requirement.
- **"Two independent solvers" / "isolated roots" / "roots are distinct for
  reasons far larger than the solver disagreement."** The Cartesian and
  logarithmic formulations share the DOP853 stepper and the same event
  machinery, and in the first pass they also shared a Cartesian launch step.
  They differ only in the coordinate the flow is written in. Their agreement
  is a **coordinate cross-check**, not independent numerical verification,
  and it is not evidence of isolation. The shared launch has been removed;
  the shared stepper has not.

No claim of certified existence or isolation is made anywhere in this work.

## Repairs to the machinery

[machinery.py](outer_fold_2026_09_09/machinery.py) replaces the first-pass
evaluator.

**Itinerary gates.** Each half return now checks, and records before it
raises: launch transversality and the launch side of the section; that the
return segment stays on the launch side for its whole length (dense output,
400 samples, not solver steps); that the side `sign(y)` is preserved; that
exactly one qualifying section crossing occurs, with unintended interior
crossings detected on the dense grid; distance to every equilibrium along the
dense trajectory; and transversality at the terminal crossing. The event is
armed by direction so that the `t = 0` departure is excluded structurally
rather than by an artificial lead-in integration. **These are numerical
safeguards. Passing them is not an enclosure.**

**Derivatives.** `dD/ds` and `dD/dtheta` now come from variational equations
integrated alongside the flow, with the event correction for the crossing
time (`dT/du = -v_x / x_dot`). Every seed root is audited against central
differences at three step sizes and two tolerances, and the *spread* is
reported as an empirical uncertainty:

| Side | `s` root | `dD/ds` (variational) | relative spread | log cross-check gap |
|---|---|---|---|---|
| upper | 0.9802427 | -1.2481e-07 | 7.1e-02 | 1.5e-12 |
| upper | 1.3899488 | +8.5268e-08 | 9.6e-02 | -7.9e-13 |
| upper | 2.0452069 | -2.4629e-07 | 9.3e-02 | 2.9e-13 |
| lower | 9.1437116 | +1.0228e-04 | 1.1e-05 | -1.3e-13 |

The three upper derivatives carry a spread of roughly 7-10%. That is an
uncertainty estimate, not an error bound, and **a derivative that is large
compared with an arbitrary threshold does not certify nonvanishing.**

**Fixed residual scales.** `D_SCALE = 1e-6`, `DS_SCALE = 1e-6`,
`THETA_SCALE = (1, 1, 1e-4, 1e-3, 1e-3)`, declared once in `machinery.py` and
never recomputed from the profile being measured. The first pass's
profile-normalised "fold score" was a heuristic and is not used.

**Root tracking.** Each tracked root is now predicted by its own sensitivity,
`s -> s - (dD/dtheta . delta) / (dD/ds)`, then corrected by a Newton iteration
confined to a guard interval of at most 35% of the distance to its nearest
neighbour. Leaving the guard, changing sort order, or two trackers closing to
within `1e-3` are all reported as **tracker failures**, distinct from any
statement about the field. Newton convergence is declared at the integrator
noise floor (`|D| < 2e-14`), since it cannot be driven below it.

## The compensated continuation

[step3_compensated.py](outer_fold_2026_09_09/step3_compensated.py). `a` is the
driver; `(b, e0, e1, e2)` compensate along the unit scaled gradient of

    J(theta) = sup { dD/ds (s; theta) : s in the outer window }

over a window starting `5e-3` past the outermost tracked upper root and 4.0
wide in `s`. The four tracked roots are *tracked, not constrained*: nothing is
imposed as an overdetermined system on their positions, and the free parameter
is the driver step alone. The target is a genuine double zero, `D = 0` and
`dD/ds = 0` together; **no sampled sign change is required to look for one.**

Caps, all predeclared: 50 accepted steps, 60 rejected steps, 10 consecutive
halvings, 40000 return evaluations per run. Actual use: 4769 return
evaluations, 296 s wall.

| Branch | accepted | rejected | stop | best `sup dD/ds` (scaled) |
|---|---|---|---|---|
| `a_up_compensated` | 11 | 44 | halving cap at `a = -1.749994` | -2.966e-01 |
| `a_down_compensated` | 10 | 44 | halving cap at `a = -1.750013` | -4.101e-03 |

Both branches stopped at the halving cap. That is a **tracker limit**, not
evidence that any cycle ceased to exist.

### What the ascent actually found

In both branches the functional rose toward zero, and in both branches it did
so by flattening `D` around an **adjacent pair of already-tracked upper
cycles**, which approached each other:

- `a_up`: upper roots end at `1.161045, 1.163860, 2.091932` — the inner pair
  separated by `2.8e-03`, with `dD/ds` there down to `-8.0e-10` and `+7.8e-10`.
- `a_down`: upper roots end at `0.848744, 1.780242, 1.783025` — the outer pair
  separated by `2.8e-03`, with `dD/ds` down to `+8.8e-10` and `-8.8e-10`.

This is an incipient saddle-node of a pair that already exists. If completed
it would **remove** two cycles, not add any. The objective was gameable: a
window beginning just past the outermost tracked root rewards flattening
caused by a merging pair inside or at the edge of it. That is a defect of the
functional as posed, and it is why the branches stalled where they did.

At the best `a_down` point the fold Jacobian of `F = (D, dD/ds)` in
`(s, a, b, e0, e1, e2)`, scaled, has **rank 2 and condition number 12.2**
(singular values `2.67e+03`, `2.18e+02`), with residual
`(D/D_SCALE, dD_ds/DS_SCALE) = (-1.24e-05, -4.10e-03)` and
`d^2D/ds^2 = -6.50e-07`. So a fold solve there is well posed numerically — but
the fold it is converging on is the annihilation of the existing pair, not a
new double zero beyond the outermost cycle.

### Outer profiles

Full profiles (22 points, `D` and variational `dD/ds` at each) are retained for
every accepted step in
[data/step3_compensated.json](outer_fold_2026_09_09/data/step3_compensated.json),
along with every rejected step and its reason. No return in any outer window
failed to resolve during the continuation (0 failures out of 22 at every
accepted step).

At every sampled point of every outer window examined, `dD/ds` was negative.
**Sampled monotonicity is not continuous monotonicity**: nothing was evaluated
between grid points, the windows are finite, and a double zero need not change
sign. No absence statement follows.

## Earlier scans, restated correctly

The first-pass scans stand as data. Their correct reading: at the seed, 60 of
60 sampled points in `s` in `[2.046, 11.0]` resolved, with `D` negative and
decreasing at every sampled point; along a 50-step branch that scaled
`(e0, e1, e2)` to 13.5x the seed while all four trackers held, the same held at
every accepted step. These are sampled profiles over finite windows on
particular branches. They do not establish that no fold exists.

## What passed, what failed, what is unresolved

**Passed.** The saved control replays bit-identically. The four seed roots
reproduce under the repaired gates with event-corrected variational
derivatives, and the log-coordinate cross-check agrees to `~1e-12`. The
repaired tracker follows the `a` direction that the first pass reported as a
failure, for 10-11 accepted steps in each direction. Fold Jacobians are rank 2
and well conditioned, so the fold solve is well posed where it was attempted.

**Failed.** Both branches hit the halving cap rather than a fold. The slope
functional is gameable and was gamed: the ascent found the annihilation fold of
an existing pair. The first pass's independence and absence claims were wrong
and are withdrawn.

**Unresolved.** Whether a genuine double zero exists beyond the outermost upper
cycle anywhere in this family. Whether a differently posed objective — one that
excludes a neighbourhood of every tracked root, so that a merging pair cannot
be rewarded — would find one. Whether any of the four numerical roots
corresponds to an actual periodic orbit in the rigorous sense: no interval
arithmetic, validated integration, or enclosure of any kind has been
implemented, and none of this bears on H16P.

## Budget

The only budget information available inside this session is the resource use
of the runs themselves: 4769 return evaluations and 296 s wall for the
continuation, a few minutes for the seed audit. The account usage meter
referenced in `H16P_NEXT_STRIKE.md` is not readable from here, so no percentage
of it can be reported honestly.
