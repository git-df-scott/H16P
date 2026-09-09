# Outer-fold strike: seed reproduced, no outer fold on the branches tried

**Status: executed and bounded. No counterexample found. No fifth cycle.**

This records the execution of Step 1 and a bounded Step 2/3 reconnaissance
of [H16P_NEXT_STRIKE.md](H16P_NEXT_STRIKE.md). Code and data:
[outer_fold_2026_09_09/](outer_fold_2026_09_09/).

## Step 1: the four-cycle seed reproduces

`reversible_reseed/verify_control.py` replays bit-identically on this
machine (all six moment rows, all 24 return differences, all four sign
brackets). Independently of it,
[step1_seed.py](outer_fold_2026_09_09/step1_seed.py) rebuilds the same
field from the exact rationals and re-derives the cycles with two solvers:
the Cartesian state `(x, y)` and a logarithmic state `(x, w = log|y|)`,
which is a different ODE rather than a rescaling. Section `{x = 0}`,
coordinate `s = log|y|`, displacement `D(s) = s_forward - s_backward`.

Four disjoint brackets refine to four isolated roots — three upper, one
lower — and every half return passes the itinerary gate (exactly one
transverse crossing, no interior crossing, no approach within `1e-3` of an
equilibrium).

| Side | `s` root | `|y|` root | `dD/ds` | half-period sum |
|---|---|---|---|---|
| upper | 0.9802427 | 2.665103 | -1.28e-07 | 5.81740 |
| upper | 1.3899647 | 4.014708 | +8.94e-08 | 5.64706 |
| upper | 2.0452066 | 7.730756 | -2.44e-07 | 5.37940 |
| lower | 9.1437116 | 9355.42 | +1.02e-04 | 3.29817 |

The two solvers agree on each root to at most `1.6e-05` in `s`, against a
minimum upper-root separation of `0.410`: four orders of margin, so the
roots are distinct for reasons far larger than the solver disagreement.
The derivatives are the weak part of this table — the upper values carry
roughly one to two significant digits (central differences at `h = 1e-4`
and `h = 1e-5` differ by up to 25% on the middle root), because `|D|`
there is `1e-08`–`1e-09` against a `1e-13` noise floor. Their signs and
their nonvanishing are solid; their magnitudes are not.

The control gate passes, so a fold search was allowed to start.

## Step 2/3 reconnaissance: the outer region has no fold

Beyond the outermost tracked upper root, `D(s)` is monotone and single
signed. At the seed, 48 of 48 grid points in `s ∈ [2.046, 6.0]` returned
cleanly, and extending to `s ∈ [2.046, 11.0]` gives 60 of 60 clean returns
— `|y|` out to about `6e4`. Zero sign changes anywhere in that range;
`|D|` grows monotonically from `2.5e-10` at the outermost root to
`5.0e-03` at `s = 11`. This is not a domain failure being read as an
absence: the return map resolved everywhere it was asked.

Two continuation branches, each capped at 50 accepted steps:

- **`perturbation_growth`** (scale `e0, e1, e2` together, holding `a, b`):
  all 50 steps accepted, reaching a perturbation 13.5x the seed's. All
  four roots tracked throughout with no interval merge and no gate
  failure. Every displacement scaled up roughly in proportion — the upper
  derivatives grew from `~1e-07` to `~1e-06`, the lower from `1.3e-04` to
  `1.4e-03` — and the outer region stayed monotone and single signed at
  every step. Growing the perturbation moves the whole picture without
  adding a zero.
- **`a_sweep`** (vary `a` alone, step `0.05`): rejected at step 1, root
  lost, no sign change within `±0.35` of the tracked position. `a` is not
  a soft direction for this configuration; any continuation through it
  needs the other four coefficients moving in compensation.

So the six-equation fold solve of Step 3 was never reached on a live
candidate: there is no fold residual to drive to zero in the region the
plan pointed at, along the directions tried. The best fold score seen
(a normalized `|D| + |dD/ds|` at a common `s`) was `3.1e-03`, and it sits
immediately outside the outermost existing root — that is the tail of the
known cycle, not a new double zero.

## What this does and does not establish

It establishes that the seed is real as a numerical object and reproduces
under an independent formulation, and that the specific outer-fold
mechanism the plan proposed does not appear in the upper annulus out to
`|y| ~ 6e4`, either at the seed or along a 13.5x perturbation growth.

It establishes nothing about H16P. The family is a restricted subset of
quadratic systems, two continuation directions are not a parameter search,
and the evidence class is NUM — floating-point integration and bisection,
no interval arithmetic, no validated integration, no certificate of
existence or isolation for any of the four roots. The 50-step cap and the
two-branch limit were predeclared bounds, and they held.

The documented obstruction is narrow and worth stating plainly: on this
seed the upper displacement outside the third cycle is monotone, and the
one coefficient direction that preserves all four cycles (uniform
perturbation growth) is also the direction that preserves that
monotonicity. A fold, if it exists in this family, needs a direction that
breaks the proportional scaling — `a` alone breaks the cycles instead,
so the next attempt would have to move `a` with `b, e0, e1, e2`
compensating, and that is a five-dimensional search this run did not do.
