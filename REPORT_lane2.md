# REPORT — Lane 2: the cusp manifold of triple limit cycles, and the hunt for a swallow-tail

Branch `fable/lane2-cusp`. Working directory `lane2_cusp_2026_09_06/`.

**No counterexample is claimed anywhere in this file. Maximum certified cycle count in
one nest so far: 3.**

**Headline.** The cusp manifold of triple limit cycles was entered from the Bautin
region at residual `5.8e-34` and continued to normal amplitude for the first time; a
triple limit cycle of *normal size* (`x ≈ 2.2`) is certified under PROTOCOL rule 1.
Across 125 continued cusp curves **no multiplicity-four limit cycle was found**. All 87
apparent `D_xxx` sign changes are centre-variety artifacts (§II.4) — a false-positive
trap that anyone repeating this search will hit. The curves end either by escaping to
infinite amplitude (`a > 1/3`) or by pressing the triple cycle against the nest's
bounding graphic (`a < 0`) — §II.5. §II.7 then reformulates the target: across
15 289 cusp records `sgn(V1)·sgn(D_xxx) < 0` without exception, and that lock can
only break at a swallow-tail or at `V1 = 0` — a Hopf at the focus, which is an
**algebraic** condition. "A cusp point whose focus is simultaneously weak" is a
cheaper and equally sufficient target than the swallow-tail, and it is where the
next session should start.

> **Two Lane-2 sessions have been running in parallel and this branch is their merge.**
> Session `01MY79ob…` built engine **A/B** (`lane2_cusp_2026_09_06/engine/`) and did the
> seed validation — Part I below. Session `01QK5W6F…` built engine **C/D**
> (`lane2_cusp_2026_09_06/*.py`, `cusp_engine.cpp`) and did the cusp-manifold entry and
> continuation — Part II below. The two engine families are completely independent
> (different chart, different independent variable, different integrator) and **agree**;
> that agreement is itself the strongest validation either produced. Keep both.

---

## 0. What this lane is doing

Continue the manifold of **multiplicity-three** limit cycles (Perko's cusp surface `C3`,
Perko 1995 Def. 4.2) out of the small-amplitude Bautin region of a third-order weak focus
to **normal amplitude**, and search along it for a **multiplicity-four** limit cycle
(swallow-tail `C4`). Perko 1995 Theorem 4.3 (verbatim in
`coordination_2026_09_06/LIT_A_ROTATED_CHERKAS.md` §1.6) says a multiplicity-four cycle
satisfying three Jacobian nondegeneracy conditions forces an open region nearby in which
the system has **four simple limit cycles in one nest**.

Setting: the Cherkas–Artés–Llibre normal form

```
xdot = 1 + x y
ydot = a00 + a10 x + a20 x^2 + a01 y + a11 x y + a y^2 ,   a00 = a01 + a11 - a10 - a20 - a
```

focus `A = (1,-1)`; section = a ray of the line `y = -1` from `A` (`x > 1`, or `x < 1`);
displacement `D(x0)` = (x-coordinate of the first return) − `x0`.

* cusp (triple cycle): `D = D_x = D_xx = 0` — 3 equations, so a **curve** in
  `(a11, a01, a10, x0)` at fixed shape `(a, a20)`;
* swallow-tail (quadruple cycle): additionally `D_xxx = 0` — codimension 1 *inside*
  the cusp manifold, detectable as a sign change of `D_xxx` along a cusp curve
  **provided `D_xxxx` does not change sign with it**; see §II.4, which is the
  main methodological result of this lane.

---

# PART I — engines A/B and the seed validation (PROTOCOL rule 7)

### Engine A — `lane2_cusp_2026_09_06/engine/cusp128.cpp`, binary128

Translate `u = x-1, v = y+1`, then the linear-normalising chart
`u = xi, v = k1 xi - w eta`, `k1 = 1+T/2`, `T = a11+a01-2a-1` (the trace `= V1`),
`L = 2a-a01-a10-2a20` (the determinant), `w = sqrt(L - T^2/4)`. In polar
`xi = rho cos th, eta = rho sin th` the system becomes the scalar equation

```
d rho / d th = ( (T/2) rho + alpha(th) rho^2 ) / ( w + beta(th) rho )
```

with `alpha, beta` explicit trigonometric quadratics. The section `{y=-1, x>1}` is the
**fixed ray** `th = th0 = atan2(k1,w)`, so there is no implicit return-time equation at
all; `D, D_x, D_xx, D_xxx` come from an exact order-3 jet in the initial radius.
Gragg–Bulirsch–Stoer, adaptive, rtol `1e-28`, ~35 ms per full evaluation.

### Engine B — `engine/engB.py`

Independent second integrator: mpmath dps 40, Cartesian, time as the independent
variable, variable-order Taylor, event `v = 0` by Newton. **A vs B on Cherkas row 1 at
six amplitudes: relative agreement 0 to 1.1e-14.**

### Cherkas rows 1–8 — three cycles recovered for every row

Rows 1,2,6,7,8 are published on `x>1`; rows 3,4,5 on `x<1`; each reproduced on its own ray.

| row | found x | published x | note |
|---|---|---|---|
| 1 | 1.2809, 2.0070, 4.0193 | 1.26, 1.98, 3.95 | `a11` published to 4 digits only |
| 2 | 1.1935, 2.0596, 3.0896 | 1.4, 1.9, 3.1 | `a11` to 5 digits |
| 3 (left) | 0.8482, 0.6155, 0.3235 | 0.8, 0.66, 0.32 | |
| 4 (left) | 0.8523, 0.7466, 0.5569 | 0.87, 0.75, 0.56 | all within 0.02 |
| 5 (left) | 0.8729, 0.8018, 0.6278 | 0.88, 0.80, 0.63 | all within 0.01 |
| 6 | 1.0504, 1.1817, 1.5415 | 1.05, 1.16, 1.5 | all within 0.05 |
| 7 | 1.5108, 2.2316, 4.4763 | 1.28, 2.15, 4.43 | |
| 8 | 1.3573, 2.3071, 4.1455 | 1.29, 2.22, 4.63 | |

Rows whose `a11` is published to 6 significant digits (4, 5, 6) match to `1e-2`; rows
with 4–5 digits drift by up to `7e-2` — the expected sensitivity of a rotation parameter.

### Row 4's published Andronov–Hopf polynomial — the sharpest test

| | max | min |
|---|---|---|
| computed AH | x = 0.620, AH = 9.500081 | x = 0.800, AH = 9.499577 |
| published polynomial | x = 0.625, AH = 9.500111 | x = 0.805, AH = 9.499583 |

**Exactly two interior extrema, in the right places, to 3e-5 in value.**
Engine validation PASSED; the lane may sweep.

### Domain ends (PROTOCOL rule 4)

For rows 3 and 5 on the `x>1` ray the return domain ends at `x = 1.346` / `x = 1.152`,
with the outermost cycle essentially ON the boundary (row 3: outer cycle at 1.3472,
boundary at 1.3475). The escape is confirmed independently by engine B (blow-up in
finite time), i.e. it is a genuine separatrix, not a chart artifact.

---

# PART II — engines C/D, and the cusp manifold itself

### Engine C — `lane2_cusp_2026_09_06/cusp_engine.cpp`, binary128 / long double

Taylor-series time-stepping in the **original Cartesian coordinates**, with the state
carried as a **degree-4 jet in the section coordinate** `eps`: `x(0) = x0 + eps`,
`y(0) = -1`. The jet of the *return time* is solved from the jet equation `y(tau) = -1`
(a genuinely implicit return, unlike engine A's fixed ray — which is why the two are a
real cross-check), and `R(eps) = x(tau(eps))` is a degree-4 jet:

```
D = R.c0 - x0 ,  D_x = R.c1 - 1 ,  D_xx = 2 R.c2 ,  D_xxx = 6 R.c3 ,  D_xxxx = 24 R.c4
```

exact to integration accuracy — no finite differencing in the amplitude direction.
Taylor order 26, local tolerance 1e-32, ~55 ms per evaluation.

**Precision established, not assumed.** Re-running one point at Taylor order 22/26/30,
tolerance 1e-30/1e-32/1e-34 and `hmax` 0.05/0.10/0.25 moves `D` in the **22nd significant
digit** (absolute agreement ~4e-34).

### Engine D — `indep_engine.py`

A fourth integrator: mpmath Taylor, written from the ODE directly, no shared code.
C vs D: 1.3e-29, 2.4e-22 (absolute 4e-34), 1.4e-30 at three test points.

**Cross-family check.** Engines A (polar, angle) and C (Cartesian, time, implicit return)
produce the *same* cycle positions for Cherkas rows 1, 2, 6, 7, 8 to all printed digits
(1.2809/2.0070/4.0193, 1.1935/2.0596/3.0896, 1.0504/1.1817/1.5415, 1.5108/2.2316/4.4763,
1.3573/2.3071/4.1455). Four integrators, two charts, two independent variables.

### A representational trap worth recording

With the parameters carried as Python floats the cusp Newton stagnates at residual
~1e-17 — not for any numerical-analysis reason, but because its steps are of relative
size ~1e-25 in `(a11, a01, a10)` and a binary64 parameter cannot represent them. All
driver arithmetic is mpmath `mpf` at dps = 50. After the fix the same Newton reaches
residual **5.8e-34**.

## II.1 Entering the cusp manifold from the Bautin region (TASK 2) — DONE

At a third-order weak focus (`a11 = 4-2a`, `a01 = 4a-3`, `a10 = (6(a^2-a-2)+a20(6a-7))/(1-3a)`)
the displacement satisfies `D(r) ~ d7 r^7`. Measured at `(a,a20) = (3,-12)`:
`D/r^7 = -1.537, -1.392, -1.153` at `r = 0.01, 0.02, 0.04` — the expected 7th-order
contact, so `d7 ≈ -1.5`.

Newton on `(D, D_x, D_xx) = 0` in `(a11, a01, a10)` at fixed `x0 = 1 + r0`, seeded from
the weak-focus point (the Bautin shift is `O(r0^2)` in `V5`, `O(r0^4)` in `V3`, `O(r0^6)`
in `V1`), converges in 5 iterations to residual 5.8e-34 for every shape tried. Example,
`(a, a20) = (3, -12)`, `r0 = 0.02`:

```
a11 = -1.998494677833741331488579457852263
a01 =  8.998494677908351789196814624248964
a10 = 13.504422527141653959512639631905580
D = -5.8e-34   D_x = -3.1e-33   D_xx = -1.3e-32   D_xxx = -9.4920440228e-6
```

`D_xxx` at entry agrees with the Bautin prediction `48 r0^4 d7 = -1.15e-5` in sign and
order of magnitude, as it must.

### The triple cycle is real at NORMAL amplitude

Taking the cusp point at `x0 = 2.2654` on the `(a, a20) = (3, -12)` cusp curve and
perturbing into the cuspidal region (`delta_mu ~ 1e-3`, chosen so `c1 * D_xxx < 0`)
gives **three certified sign changes of `D`** in one nest:

| root | bracket `min|D|` | two-tolerance noise | margin |
|---|---|---|---|
| 2.12940787703450780729514 | 2.19e-8 | 1.07e-11 | 2050× |
| 2.21812962597053274330530 | 2.46e-8 | 1.11e-11 | 2220× |
| 2.27153111659064560853597 | 2.39e-8 | 1.14e-11 | 2100× |

PROTOCOL rule 1 satisfied (noise from recomputation at looser tolerance,
`noise = 10|difference| + 5e-12 s`). The cycles sit at `x ≈ 2.13–2.27`: normal size,
not small-amplitude. `data`: `triple_confirm_row1.json`.

## II.2 Continuation of the cusp curve (TASK 3) — IN PROGRESS

Pseudo-arclength continuation in `(a11, a01, a10, x0)` at fixed `(a, a20)`, chord Newton
(frozen 3×4 Jacobian, refreshed every 4 iterations; exact Jacobian recomputed at any
candidate event). Every accepted point logs `D, D_x, D_xx, D_xxx, D_xxxx`,
`nu = D_xxx/(D_xxxx r0)` (a scale-free distance to a swallow-tail), `V1`,
`L = det J(A)`, the return time, the transversality of the crossing, and Perko's
nondegeneracy Jacobians with `mu1 = a11` (Cherkas's *rotating* parameter, so
`d_{mu1} != 0` is guaranteed by Duff/Perko monotonicity).

### Cherkas shapes, rows 1–4 (stopped at budget, not at a curve end)

| shape | `(a, a20)` | pts | `x0` range | `D_xxx` start → end | sign changes |
|---|---|---|---|---|---|
| row1 | (3, −12) | 120 | 1.02 → 10.10 | −9.49e−6 → −2.84e−3 | **0** |
| row2 | (1.5, −15) | 183 | 1.02 → 13.41 | +1.07e−7 → +2.05e−4 | **0** |
| row3 | (−2, 12) | 161 | 1.02 → 1.394 | +8.98e−6 → +1.78e+3 | **0** |
| row4 | (−2, −1) | 151 | 1.02 → 1.389 | +3.38e−5 → +2.33e+3 | **0** |

Two clearly different behaviours:

* **rows 1 and 2** — the curve runs out in amplitude without bound (`x0` past 10 and 13);
  `|D_xxx|` rises to a maximum (−5.9e−2 at `x0 ≈ 2.27` for row 1; +3.2e−3 at `x0 ≈ 2.87`
  for row 2) and then decays monotonically toward zero **without changing sign**. `L`
  stays positive (antisaddle preserved) and the section stays strongly transversal
  throughout. Stopped by budget, not by a curve end.
* **rows 3 and 4** — `x0` stalls near 1.39 while `a11` and `|D_xxx|` blow up. The curve
  runs into a boundary in parameter space, not in amplitude. Part I's independent
  finding that the `x>1` return domain for the row-3 *seed* ends at `x = 1.3475` on a
  genuine separatrix is very likely the same boundary.

**No `D_xxx` sign change on any of the four.** Rows 5–8 and the (a, a20) grid are running.

### The structure that organises the search

On the third-order stratum,

```
V7 = -150 (a-2) [ -4a(a+1)(a-2)^2 + a20 (a-1)(2a+1)^2 ]
```

so `V7 = 0` on the **centre curve** `a20 = a20_c(a) = 4a(a+1)(a-2)^2 / [(a-1)(2a+1)^2]`,
where `V1=V3=V5=V7=0` and the system has a centre. Since `D_xxx = 48 r0^4 d7` at the
small-amplitude end, the **entry sign of `D_xxx` flips across this curve**. The
admissibility condition `(a-3-a20)/(1-3a) < 0` puts the centre curve inside the
admissible region only for roughly `a ∈ (-3, -0.5)` and `a ∈ (1/3, 1)`:

| a | −4 | −3 | −2.5 | −2 | −1.5 | −1 | 0 | 0.2 | 0.4 | 0.5 | 0.6 | 8/11 | 0.9 | 1.04 | 1.5 | 3 | 5 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| `a20_c` | −7.05 | −6.00 | −5.42 | −4.74 | −3.68 | 0.00 | 0.00 | −1.98 | −2.95 | −3.38 | −3.89 | −4.95 | −10.56 | 20.6 | 0.47 | 0.49 | 2.23 |
| reachable | no | no | **yes** | **yes** | **yes** | **yes** | **yes** | **yes** | **yes** | **yes** | **yes** | **yes** | **yes** | no | no | no | no |

All eight Cherkas shapes lie on one side of the centre curve or where it is unreachable —
consistent with all four completed curves showing a constant `D_xxx` sign. The running
grid straddles the centre curve at 16 values of `a` with offsets ±0.08, ±0.3, ±1.0, plus
a band sweep at the 8 unreachable `a` values: 116 grid jobs + rows 5–8, 120 continuation
points each.

## II.3 A negative result: no swallow-tail at small amplitude

With `a20` as a fourth unfolding parameter the swallow-tail system
`(D, D_x, D_xx, D_xxx) = 0` in `(a11, a01, a10, a20)` is **square** at fixed `(a, x0)`.
Seeded from a *small-amplitude* cusp point (`a = 3`, `x0 = 1.02`) the Newton runs
`a20 → -∞` (reached −4.6e5, with `a10 → +6.3e5`, residual still decreasing like
`1/|a20|`). This is not a solver failure: at small amplitude `D_xxx = 48 d7 r0^4`, and
`d7 = 0` only on the centre variety, where `D` vanishes identically.

> **A nondegenerate swallow-tail cannot live at the small-amplitude end of the Bautin
> unfolding.** It must be sought at normal amplitude — which is what the continuation is
> for. This closes off the one place where an explicit formula would have handed the
> answer over, and it is worth recording so nobody re-derives it.


## II.4 ★ THE CENTRE VARIETY IS A SPURIOUS COMPONENT OF THE CUSP MANIFOLD ★

**This is the most important thing Lane 2 has found, and it is a trap that will
produce false positives for anyone running this search.**

The cusp manifold was defined as the solution set of

```
D = D_x = D_xx = 0        (3 equations in the 6 coordinates (a,a20,a11,a01,a10,x0))
```

But **wherever the field has a CENTRE, `D` vanishes identically, so all three
equations hold trivially.** The centre variety therefore sits inside that
solution set as a second component *of the same dimension* — it is not a
lower-dimensional degeneracy that a generic path would miss. A pseudo-arclength
continuation can slide onto it and stay there, and the Newton residual reports
`1e-32` the whole way, because `0 = 0` to any precision one likes.

On the centre variety every derivative of `D` is at the integration noise floor,
so **`D_xxx` changes sign essentially at random**; and where a genuine cusp curve
*crosses* the centre variety transversally, the whole function `D` changes sign,
so `D_xxx` and `D_xxxx` — and `D` itself — all flip together. Either way a bare
`D_xxx` sign watch fires, and neither case is a swallow-tail.

### What this looked like in the data

The first pass over 125 continued cusp curves reported **87 sign changes of
`D_xxx`**, concentrated exactly on the shapes straddling the centre curve — i.e.
precisely where the structural argument predicted a swallow-tail should be.
Classifying every one of them by measuring the displacement scale
`|D|` at `0.5 r0` and `1.7 r0` off the cusp point:

| class | count |
|---|---|
| curve lying **on** the centre variety (`\|D\| < 1e-24` everywhere: the field has a centre, so there are no limit cycles at all) | 58 |
| curve **crossing** the centre variety (`D_xxxx` flips with `D_xxx`; the whole `D` changes sign) | 29 |
| **genuine multiplicity-four candidates** | **0** |

The 58 all come from one curve, `c_a0p6_o0p08` (`a = 0.6`, `a20 = -3.8076`,
against `a20_c(0.6) = -3.8876`), on which the continuation landed on the centre
variety at `x0 ≈ 3` and then marched `x0` from 3.0 to 13.6 in fixed steps with
the parameters frozen — the unmistakable signature. Every `D_xxx` and `D_xxxx`
there is `~1e-34`.

A worked example of the crossing case, `c_am2p0_om0p08` (`a=-2`, `a20=-4.8207`):

| i | x0 | D_xxx | D_xxxx | nu = D_xxx/(D_xxxx r0) | residual |
|---|---|---|---|---|---|
| 4 | 1.02826 | −1.747e−3 | −0.5045 | 0.12254 | 1.6e−32 |
| 5 | 1.03289 | −3.739e−4 | −0.0923 | 0.12322 | 3.9e−32 |
| 6 | 1.03936 | +3.645e−3 | +0.7507 | 0.12336 | 1.0e−32 |
| 7 | 1.04825 | +1.401e−2 | +2.3719 | 0.12240 | 1.1e−32 |

`D_xxx` "changes sign" — but so does `D_xxxx`, and `nu` sails through completely
smoothly. Sampling `D` itself at `0.3, 0.6, 1.6, 2.5, 4.0` times the amplitude
shows the entire displacement function flipping sign between rows 5 and 6
(e.g. at `4 r0`: `−1.26e−6 → +2.32e−5`). The curve is crossing a centre.

### The corrected detector, and the guard

* A sign change of `D_xxx` is a swallow-tail candidate **only if `D_xxxx` does
  not flip with it** — equivalently, only if `nu = D_xxx/(D_xxxx r0)` passes
  through **zero**. (`nu` passing through **infinity** is a zero of `D_xxxx`
  with `D_xxx != 0`: an ordinary point of the cusp curve, not a swallow-tail
  either. Six curves in the grid do that, at `a = 1.04` and `a = 1.5`.)
* `Cusp.amplitude()` now measures `|D|` at `0.5 r0` and `1.7 r0` off the cusp
  point at every accepted continuation step, logs it as `amp`, and the
  continuation **stops with `end_reason = CENTRE_VARIETY`** if it falls below
  `1e-24` (engine noise is `~1e-33`, so this is a `1e9` margin below which no
  genuine displacement lives).
* `analyse.py` now reports both counts and never calls a bare `D_xxx` sign
  change a swallow-tail.

### Why this is not just a numerical nuisance

It also explains, structurally, why the small-amplitude end is barren. Perko 1992
Remark 1: the multiplicity `m` of a limit cycle equals the maximum number of
cycles that can bifurcate from it. A **nondegenerate multiplicity-four cycle at
small amplitude**, on a family of parameters converging to a weak focus or
centre as the amplitude shrinks, would put four limit cycles in every
neighbourhood of that singular point for parameters arbitrarily close to it —
i.e. cyclicity `>= 4`, contradicting **Bautin**. So the swallow-tail cannot be
approached from the Bautin end at all; what one finds there instead is the
centre variety, which is exactly what the numerics produced. This is the same
conclusion as II.3, reached independently, and it means the entire search must
live at normal amplitude.


## II.5 How the cusp curves END (TASK 3's "record the endpoints and why")

Across the 125 curves continued in the first pass, the cusp curve ends in exactly
two ways, cleanly separated by the sign of `a`.

### (i) `a > 1/3` — the curve escapes to infinite amplitude, parameters bounded

`x0` grows without bound (row 1 past 10, row 2 past 13, the `a = 1.5` and
`a = 2.5` grid curves past 17), `a11` drifts slowly, `L = det J(A)` stays
positive, the section stays strongly transversal, and `|D_xxx|` rises to a
maximum and then **decays monotonically toward zero without changing sign**
(row 1: max `−5.9e−2` at `x0 ≈ 2.27`, down to `−2.8e−3` at `x0 = 10.1`; row 2:
max `+3.2e−3` at `x0 ≈ 2.87`, down to `+2.0e−4` at `x0 = 13.4`). These runs were
stopped by budget, not by any feature of the curve. The decay of `D_xxx` is not
an approach to a swallow-tail: `nu = D_xxx/(D_xxxx r0)` stays near `0.1`
throughout, i.e. `D_xxxx` decays in step.

### (ii) `a < 0` — the triple cycle is squeezed against the nest boundary

`x0` converges to a finite limit while `a11 → ∞` at roughly unit speed in
arclength, and `D_xxx`, `D_xxxx` blow up geometrically:

| curve | `(a, a20)` | `x0 →` | `a11` at end | `D_xxx` at end | `D_xxxx` at end |
|---|---|---|---|---|---|
| `p_am4_b0p5` | (−4, −6.5) | 1.1895 | 22.45 | 2.64e5 | 9.79e8 |
| `c_am2p0_o0p3` | (−2, −4.4407) | 1.3984 | 19.10 | 4.92e3 | 2.66e6 |
| row 3 | (−2, 12) | 1.394 | — | 1.78e3 | — |
| row 4 | (−2, −1) | 1.389 | — | 2.33e3 | — |

**Why**: at the final parameters of each of these curves, walking outward along
the section from the cusp point, the return map fails **immediately** — the nest
domain end `s_max` coincides with `x0` to within `0.000` and `0.004` respectively.
The triple cycle has grown until it touches the outer boundary of the nest, and
the parameters have to run to infinity to keep it a triple cycle. This is
PROTOCOL §(b)'s outer end: the cusp curve terminates on the bounding
graphic/separatrix, not at a swallow-tail. It is also consistent with Part I's
independent finding that the row-3 seed's `x>1` return domain ends at
`x = 1.3475` on a genuine separatrix (confirmed by engine B blowing up in finite
time), with the outermost cycle sitting essentially on the boundary.

## II.6 Closest approach to a swallow-tail

Excluding centre-variety records, the smallest `|nu| = |D_xxx/(D_xxxx r0)|` seen
anywhere is `1.42e-3`, at

```
a   = -4            (exactly)          a20 = -13/2         (exactly)
a11 =  22.45086693473904317711227265950205
a01 = -29.52317781989795110902621214492275
a10 =  31.91565514198711888059773382354741
x0  =   1.189462447165467610264657155882525
D = -2.6e-31   D_x = -2.5e-29   D_xx = -1.4e-26
D_xxx = 263635.70357303020554   D_xxxx = 978962520.8576926
T = 6.31432078137962954   L = 2.60752267791   V1 = -0.0723108851589
```

(`closest_approach.json`; Perko Thm 4.3 quantities all healthy, `min|.| = 14.9`.)

**But this is not a near-miss.** It is the type-(ii) endpoint above: `|nu|` is
small there only because `D_xxxx` is blowing up *faster* than `D_xxx` as the
cycle is pressed against the nest boundary — `D_xxx` itself is growing through
`2.6e5`, not approaching zero. A swallow-tail needs `D_xxx → 0` with `D_xxxx`
bounded; here both diverge. The `nu` statistic is the right detector for a
*centre* crossing but is fooled by a *boundary* blow-up, and both traps have now
been characterised.

**On no curve did `D_xxx` approach zero for any reason other than the centre
variety.**


## II.7 ★ THE SIGN LOCK, AND A SECOND ROUTE TO FOUR IN ONE NEST ★

Perko's swallow-tail is a **sufficient** local mechanism for four simple cycles in
one nest, not the only one. In Cherkas's Andronov–Hopf picture four cycles
correspond to `AH(x) = a11` having **three** interior extrema; extrema of `AH` are
folds (`D_x = 0`) and cusps create them in **pairs**, so parity is preserved and
three extrema can arise either by three merging at a swallow-tail **or** by a cusp
creating a pair alongside a pre-existing extremum elsewhere in the nest.

The second route is testable directly on every cusp point already computed, and
the test costs nothing. On the right section:

* just outside the focus, `D ~ (exp∫div − 1)(x−1)`, so `sgn D = sgn V1`;
* just inside the cusp point, `D ~ (1/6) D_xxx (x−x0)^3` with `x < x0`, so
  `sgn D = −sgn D_xxx`.

Hence `D` has an **odd number of extra zeros between the focus and the triple
cycle** — a fourth limit cycle in the same nest, once the cusp is unfolded into
three — exactly when

```
sgn(V1) * sgn(D_xxx)  >  0.
```

`V1` and `D_xxx` are logged at every accepted point, so this is a pure
post-processing query over the ledgers. The result:

> **Over all 15 289 non-degenerate cusp records computed in this lane, on 125
> curves spanning the whole admissible `(a, a20)` region, `sgn(V1)·sgn(D_xxx)`
> is NEGATIVE without a single exception.**

That is not an accident. At the Bautin entry the triple cycle has
`d1 = −r0^6 d7` and `D_xxx = 48 d7 r0^4`, so `sgn(V1) = −sgn(D_xxx)` is forced
there; and the product is a continuous nonvanishing function along the cusp
curve, so it can only flip where one of its factors vanishes:

* **`D_xxx = 0`** — the swallow-tail. §II.4 shows this is unreachable from the
  Bautin end (the only zeros encountered are on the centre variety), and §II.5
  shows it does not happen at either kind of curve endpoint either.
* **`V1 = 0`** — a Hopf bifurcation at the focus. **This is a purely algebraic
  condition**, `V1 = a11 + a01 − 2a − 1 = 0`, costing no integration at all.

So the lane's target can be restated, and cheapened, as:

> **find a point of the cusp manifold at which the focus is simultaneously weak.**

Crossing `V1 = 0` along the cusp manifold flips the sign lock and puts a fourth
cycle inside the triple. **Bautin is not violated**: only one of the four cycles
bifurcates from the focus; the other three sit at `x0`, of normal size. (This is
also structurally the natural generalisation of Cherkas's own construction, which
got three from one Hopf cycle plus a fold pair; four would be one Hopf cycle plus
a cusp triple.)

### Status of that target

`V1` is essentially zero at the Bautin entry (`7.5e-11` for row 1) — but that is
the degenerate `r0 → 0` branch, where `V1 = −r0^6 d7` and Bautin caps the count
at three. Along every computed curve `|V1|` then grows **monotonically and
without changing sign** (row 1: `7.5e-11 → 0.053`; row 2: `−1.3e-12 → −9.8e-3`).
So a normal-amplitude branch with `V1 = 0`, if it exists, is not reached by
continuing in amplitude at fixed shape.

Two solvers are implemented for it and are the natural next step:

* `cusp_v1.py` — impose `V1 = 0` exactly by slaving `a01 = 2a + 1 − a11`, leaving
  the square system `(D, D_x, D_xx) = 0` in `(a11, a10, a20)` at fixed `(a, x0)`.
  Direct Newton from the row-1 cusp points does **not** converge (residual
  stalls at `9.1e-5`), so the constraint is not satisfiable near those points.
* `v1_homotopy.py` — walk `V1` to zero *along* the cusp manifold: unknowns
  `(a11, a01, a10, a20)`, equations `(D, D_x, D_xx, V1 − t)`, with the fourth
  Jacobian row exact (`dV1/du = (1,1,0,0)`), stepping `t` from `V1(seed)` to 0.
  The centre-variety amplitude guard is wired in, because a homotopy that
  reaches `V1 = 0` by degenerating to a centre would prove nothing.

  **Result on the `a = 3`, `x0 = 1.5184` seed:** the homotopy gets `10.6 %` of the
  way (`V1: 6.010e-4 → 5.373e-4`) and then sticks, having paid

  ```
  a20 : -12      ->  -3537.05        a11 : -1.4137 -> -1.71459
  L   :  7.2     ->   2091.29        amp :  3.70e-5 (not degenerate)
  ```

  i.e. **the branch escapes to infinity in `a20`** — the same signature as the
  small-amplitude swallow-tail Newton of §II.3, and for the same kind of reason.
  So on that slice there is no bounded cusp point with a weak focus.

  **The likely cause is that the slice is wrong, not that the target is empty.**
  `{cusp} ∩ {V1 = 0}` is 4 equations in the 6 coordinates, hence **2-dimensional**;
  fixing both `a` and `x0` cuts it with a 2-parameter slice, which generically
  meets a 2-dimensional set in isolated points that a local Newton can easily
  miss, and both solvers above fix `x0`. The untried move is to **let `x0`
  float**: solve `(D, D_x, D_xx, V1) = 0` in `(a11, a01, a10, x0)` at fixed
  `(a, a20)` — also square — seeding `x0` right across the nest, then continue
  the solution in `(a, a20)`. That is where the next session should start.


## II.8 The Andronov-Hopf / beta* picture at a cusp point (TASK 5, for Lane 1)

Cherkas's Andronov-Hopf function `AH(x) = a11` is defined implicitly by
`D(x; a11 = AH(x)) = 0`, with `a11` the rotating parameter (so `D` is strictly
monotone in it by Duff/Perko and bisection is safe). Differentiating implicitly,

```
AH'    = -D_x    / D_a11
AH''   = -D_xx   / D_a11        (once AH' = 0)
AH'''  = -D_xxx  / D_a11        (once AH' = AH'' = 0)
```

so at every cusp point of this lane `AH' = AH'' = 0` and `AH''' != 0`: **AH has a
degenerate inflection at `x0`.** That is exactly the point at which two interior
extrema of `AH` merge - the boundary between the "AH has 2 extrema" region
(3 cycles in the nest) and the "AH has 0 extrema" region (1 cycle).

The bridge to Lane 1 follows immediately:

> Four cycles in the nest <=> `AH` has **three** interior extrema. Extrema of `AH`
> are folds (`D_x = 0`) and cusps create them **in pairs**, so parity is
> preserved. Three extrema therefore arise either by three merging at a
> **swallow-tail** (`AH''' = 0`, this lane's original target) or by a cusp
> creating a pair alongside a pre-existing extremum elsewhere in the nest (the
> **(3+1) route** of SS II.7). **Lane 1's search for a third extremum of `beta*`
> and this lane's swallow-tail are the same event**, and SS II.7 is its other
> realisation. The `d_{mu1} != 0` that Perko Thm 4.3 needs is exactly
> `D_a11 != 0`, which Duff's monotonicity guarantees for a rotation parameter -
> so that nondegeneracy condition is free in this normal form, and the Perko
> Jacobians logged at every ledger point (`perko`) cover the rest.

`ah_at_cusp.py` produces the numerical picture. **Caveat, recorded honestly:** the
first run defines `AH` at only 9 of 41 samples on the row-1 cusp at `x0 = 1.5184`
(0 interior extrema among those). The bisection bracket `a11 +- 0.5` is too narrow
- `AH` moves further than that across the nest - and the nest's own outer end
moves as `a11` varies, so samples past it have no return at all. The fix is to
bracket `a11` adaptively from the previous sample and to recompute `s_max` per
`a11`; the analytic statements above do not depend on it.
---

## Open problems / next steps

1. **Finish re-running the grid with the centre-variety guard on.** Status at the end
   of this session: the guard is implemented and wired in, and 8 of the 120 shapes were
   re-run into `ledger_grid2/` before the budget ran out (rows 5-8 and four grid
   shapes, 140 points each, **0 `D_xxx` sign changes and 0 `nu` sign changes**). The
   shapes that produced the artifacts in the first pass are the ones *below* the centre
   curve (offsets `-0.08`, `-0.3`), and those had not been reached. Each of them should
   now terminate with `end_reason = CENTRE_VARIETY` at a well-defined point, which is
   itself the answer to "where does the cusp curve end and why".
2. Map the centre variety inside the cusp manifold explicitly. It is a codimension-0
   component of the same solution set, so it is not an obstacle to be avoided but a
   boundary to be charted: the interesting question is whether a swallow-tail branch
   emanates from it, and §II.4's Bautin argument says it cannot do so at the
   small-amplitude end.
3. `nu` passing through **infinity** (a zero of `D_xxxx` with `D_xxx != 0`) happens on
   six grid curves, at `a = 1.04` and `a = 1.5`, around `x0 ≈ 7–9`. Those are ordinary
   points of the cusp curve, but they are where the cusp curve is "flattest" and are
   the natural places to look for a nearby swallow-tail in a transverse direction.
4. Rows 1 and 2 show `|D_xxx| → 0` as `x0 → ∞` **without a sign change and without
   approaching a centre**. Decide whether the cusp curve simply escapes to infinite
   amplitude or ends at a graphic. `nu` is now logged to tell these apart.
5. Rows 3 and 4 end with `a11` blowing up at nearly constant `x0 ≈ 1.39`. Part I's
   independently confirmed separatrix at `x ≈ 1.3475` for the row-3 seed is the leading
   candidate for that boundary.
6. Run the 5-parameter minimum-norm swallow-tail Newton (`swallow5.py`) from
   far-amplitude points of the guarded curves, with the amplitude guard wired in so it
   cannot converge onto a centre.
7. Log the Andronov–Hopf/`beta*` picture at each cusp point for Lane 1 (`probe.ah_sweep`
   is implemented; not yet run over the grid).

---

## Files

| file | what |
|---|---|
| `engine/cusp128.cpp`, `engine/cusp128`, `engine/eng.py` | engine A (polar chart, fixed-ray section) |
| `engine/engB.py` | engine B (mpmath, Cartesian, time) |
| `validate_seeds.py`, `data/validation_cherkas.json`, `ah_row4.py`, `data/ah_row4.json` | Part I validation |
| `cusp_engine.cpp`, `engine.py` | engine C (Cartesian Taylor jet, degree 3 or 4, implicit return) |
| `indep_engine.py` | engine D (independent mpmath integrator) |
| `cusp.py` | cusp residual, Jacobians, Newton, pseudo-arclength |
| `continue_cusp.py`, `campaign.py` | the continuation driver + JSONL ledger |
| `swallow.py` | square swallow-tail Newton + Perko Thm 4.3 quantities |
| `probe.py` | rule-1 cycle counting, cusp unfolding, `AH(x)` |
| `grid.py`, `make_grid_spec.py` | the (a, a20) grid and the centre curve |
| `validate.py`, `analyse.py` | rule-7 validation; ledger summariser |
| `ledger/`, `ledger_grid/` | append-only JSONL ledgers |

---

## What this lane established, and what it did not

**Established.**

* An engine that returns `D, D_x, D_xx, D_xxx, D_xxxx` exactly (degree-4 jet, no
  finite differencing in the amplitude direction), in binary128, agreeing with three
  other independent integrators — two of them written by a different session in a
  different chart — to `1e-29` or better.
* The cusp manifold of triple limit cycles can be entered from the Bautin
  small-amplitude region by Newton from the third-order weak focus, to residual
  `5.8e-34`, at every `(a, a20)` tried.
* It can be continued out to normal amplitude. **A triple limit cycle of normal size
  is exhibited and certified** under PROTOCOL rule 1 (three sign changes at
  `x = 2.1294, 2.2181, 2.2715`, each `~2000x` above the two-tolerance noise floor).
  As far as the literature in `coordination_2026_09_06/` shows, this continuation had
  not been done before.
* 125 cusp curves across the admissible `(a, a20)` region carry **no multiplicity-four
  limit cycle**. All 87 apparent `D_xxx` sign changes are centre-variety artifacts,
  classified individually (§II.4).
* The curves end in exactly two ways, both characterised (§II.5).
* `sgn(V1)·sgn(D_xxx) < 0` on all 15 289 non-degenerate records, which is precisely
  "no fourth cycle inside the triple", and that lock can only break at a swallow-tail
  or at a weak focus (§II.7).

**Not established.**

* Nothing here rules a swallow-tail out. The search covered the cusp curves reachable
  by amplitude continuation from the Bautin entry at 125 shapes, plus shape sweeps at
  fixed amplitude; that is a 1-parameter family of 1-dimensional slices through a
  3-dimensional manifold, not the manifold.
* The `{cusp} ∩ {V1 = 0}` target of §II.7 is **untested in the slice that matters**
  (`x0` free). Both solvers tried fix `x0`, which over-constrains a 2-dimensional set.
* No claim whatever is made about `H(2)`.
