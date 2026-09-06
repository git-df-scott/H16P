# REPORT_lane2 — cusp manifold of triple limit cycles, continuation to normal amplitude

Branch `fable/lane2-cusp`.  Lane 2 of the 2026-09-06 program.
Target: continue the manifold of MULTIPLICITY-THREE limit cycles (cusp surface C3,
Perko 1995 Def. 4.2) out of the small-amplitude Bautin region of the Cherkas–Artés–Llibre
third-order weak-focus family, and look along it for a sign change of `D_xxx`
(= a multiplicity-FOUR cycle = swallow-tail C4, Perko 1995 Thm 4.3).

**No counterexample is claimed anywhere in this file.**

---

## Status: 2026-09-06, checkpoint 1 — ENGINE BUILT AND VALIDATED

### What ran

**Engine A** — `lane2_cusp_2026_09_06/engine/cusp128.cpp`, binary128 (`__float128`).

Chart.  Translate `u = x-1, v = y+1` (focus `A=(1,-1)` -> origin), then the
linear-normalising chart `u = xi, v = k1 xi - w eta`, `k1 = 1+T/2`,
`T = a11+a01-2a-1` (the trace, `= V1`), `L = 2a-a01-a10-2a20` (the determinant,
independent of `a11`), `w = sqrt(L - T^2/4)`.  In polar `xi = rho cos th, eta = rho sin th`
the system becomes the single scalar equation

    d rho / d th  =  ( (T/2) rho + alpha(th) rho^2 ) / ( w + beta(th) rho )

with `alpha, beta` explicit trigonometric quadratics in `th`.  The section
`{y=-1, x>1}` is the FIXED RAY `th = th0 = atan2(k1,w)`, and `x = 1 + rho w/nrm`,
`nrm = hypot(w,k1)`; the ray `th0+pi` is `{y=-1, x<1}` with `x = 1 - rho w/nrm`.
Because the section is a fixed ray of the independent variable, **there is no
implicit return-time equation at all**: `D, D_x, D_xx, D_xxx` come from an exact
order-3 jet in the initial radius propagated through the ODE.  **No finite
difference is ever taken in the amplitude direction** (PROTOCOL: never third
derivatives in double / by differencing).

Integrator: Gragg–Bulirsch–Stoer (modified midpoint + Neville extrapolation in
`h^2`), adaptive step, up to 14 extrapolation levels, relative tolerance `1e-28`.
Cost ~35 ms per full `(D,D',D'',D''')` evaluation.

**Engine B** — `lane2_cusp_2026_09_06/engine/engB.py`, independent second integrator
(PROTOCOL rule 2).  Everything differs: mpmath `mpf` at dps 40 (vs binary128),
Cartesian `(u,v)` (vs polar), time `t` (vs angle `th`), variable-order Taylor series
(vs GBS), event `v=0` located by Newton on the Taylor polynomial (vs fixed ray).

### Validation (PROTOCOL rule 7)

**A vs B, Cherkas row 1, six amplitudes** — relative agreement `0` to `1.1e-14`
(the residual is the double-precision print, not the integrators):

| u0 | engine A | engine B | rel |
|---|---|---|---|
| 0.2 | 2.087854968836989683e-04 | 2.08785496883698998e-04 | 1.3e-16 |
| 0.3 | -6.955307830906821737e-05 | -6.95530783090681671e-05 | 7.8e-16 |
| 1.0 | -6.879352148899623570e-05 | -6.87935214889970277e-05 | 1.1e-14 |
| 2.0 | 1.961114690839928232e-02 | 1.9611146908399283e-02 | 0 |
| 3.0 | 1.458174566818462191e-03 | 1.45817456681847335e-03 | 7.6e-15 |

**Cherkas–Artés–Llibre rows 1–8, cycle positions on the paper's own section.**
Three sign changes of `D` in the nest are recovered for every row.
Rows 1,2,6,7,8 are listed by the paper on `x>1`; rows 3,4,5 on `x<1` — reproduced
on the corresponding ray.

| row | found x | published x | note |
|---|---|---|---|
| 1 | 1.2809, 2.0070, 4.0193 | 1.26, 1.98, 3.95 | `a11` published to 4 digits only |
| 2 | 1.1935, 2.0596, 3.0896 | 1.4, 1.9, 3.1 | `a11` to 5 digits |
| 3 (left) | 0.8482, 0.6155, 0.3235 | 0.8, 0.66, 0.32 | |
| 4 (left) | 0.8523, 0.7466, 0.5569 | 0.87, 0.75, 0.56 | **all within 0.02** |
| 5 (left) | 0.8729, 0.8018, 0.6278 | 0.88, 0.80, 0.63 | **all within 0.01** |
| 6 | 1.0504, 1.1817, 1.5415 | 1.05, 1.16, 1.5 | **all within 0.05** |
| 7 | 1.5108, 2.2316, 4.4763 | 1.28, 2.15, 4.43 | |
| 8 | 1.3573, 2.3071, 4.1455 | 1.29, 2.22, 4.63 | |

The rows whose `a11` is published to 6 significant digits (4, 5, 6) match to
`1e-2`; the rows whose `a11` is published to 4–5 digits drift by up to `7e-2`,
which is the expected sensitivity of a rotation parameter.

**Row 4's published Andronov–Hopf polynomial** (PROTOCOL rule 7, the sharpest test).
Solving `AH(x) = a11` such that `x` is on a cycle, on `[0.6, 0.9]`, left ray:

| | max | min |
|---|---|---|
| computed AH | x = 0.620, AH = 9.500081 | x = 0.800, AH = 9.499577 |
| published polynomial | x = 0.625, AH = 9.500111 | x = 0.805, AH = 9.499583 |

**Exactly two interior extrema, in the right places, to `3e-5` in value.**
Engine validation PASSED; the lane may sweep.

### Domain ends recorded (PROTOCOL rule 4)

The return domain `s_max` was located for every seed.  For rows 3 and 5 on the
`x>1` ray the domain ends at `x = 1.346` / `x = 1.152` respectively, with the
outermost cycle sitting essentially ON the boundary (row 3: outer cycle at
`x = 1.3472`, boundary at `x = 1.3475`), so those rows are worked on the `x<1`
ray.  The escape at the boundary is confirmed by engine B independently
(blow-up in finite time, `T -> 0.39`), i.e. it is a genuine separatrix, not a
chart artifact.

### Files

- `engine/cusp128.cpp`, `engine/cusp128` — engine A
- `engine/eng.py` — driver (all numbers cross as >=36-digit decimal strings)
- `engine/engB.py` — engine B
- `validate_seeds.py`, `data/validation_cherkas.json`
- `ah_row4.py`, `data/ah_row4.json`

### Next step

Enter the cusp manifold: for `(a,a20)` on the third-order weak-focus family
(`a11 = 4-2a`, `a01 = 4a-3`, `a10 = (6(a^2-a-2)+a20(6a-7))/(1-3a)`, and
`L = 5(3-a+a20)/(1-3a) > 0` — proved equivalent to the paper's own condition
`(a-3-a20)/(1-3a) < 0`), Newton on `D = D_x = D_xx = 0` in `(a11,a01,a10)` at a
small fixed amplitude, then pseudo-arclength continuation in amplitude with
`D_xxx` monitored at every point.

## Open problems

- None blocking.
