# REPORT — Lane 2: the cusp manifold of triple limit cycles, and the hunt for a swallow-tail

Branch `fable/lane2-cusp`. Working directory `lane2_cusp_2026_09_06/`.

**No counterexample is claimed anywhere in this file. Maximum certified cycle count in
one nest so far: 3.**

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
  the cusp manifold, detectable as a **sign change of `D_xxx`** along a cusp curve.

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

---

## Open problems / next steps

1. Finish the (a, a20) grid; tabulate `sgn D_xxx` at the near and far end of every cusp
   curve. A shape where the two differ contains a swallow-tail.
2. Rows 1 and 2 show `|D_xxx| → 0` as `x0 → ∞` **without a sign change**. Decide whether
   that decay is asymptotic (no swallow-tail; the cusp curve escapes to infinite
   amplitude) or whether the curve ends at a graphic first. `nu = D_xxx/(D_xxxx r0)` is
   now logged to tell these apart.
3. Rows 3 and 4 end with `a11` blowing up at nearly constant `x0`. Identify the boundary
   (Part I's separatrix at `x ≈ 1.3475` is the leading candidate).
4. Run the swallow-tail Newton (with `a20` free) from **far-amplitude** cusp points
   across the whole grid.
5. Log the Andronov–Hopf/`beta*` picture at each cusp point for Lane 1 (`probe.ah_sweep`
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
