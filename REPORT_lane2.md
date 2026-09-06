# REPORT — Lane 2: the cusp manifold of triple limit cycles, and the hunt for a swallow-tail

Branch `fable/lane2-cusp`. Working directory `lane2_cusp_2026_09_06/`.
Session started 2026-09-06. **Status: IN PROGRESS — this file is updated every 30–45 min.**

**No counterexample is claimed. Maximum certified cycle count in one nest so far: 3.**

---

## 0. What this lane is doing

Continue the manifold of **triple** limit cycles (Perko's cusp surface `C3`) out of the
Bautin small-amplitude region of a third-order weak focus to **normal amplitude**, and
search along it for a **multiplicity-four** limit cycle (Perko's swallow-tail `C4`).
Perko 1995 Theorem 4.3 (restated verbatim in `coordination_2026_09_06/LIT_A_ROTATED_CHERKAS.md`
§1.6) says a multiplicity-four cycle satisfying three Jacobian nondegeneracy conditions
forces an open region nearby in which the system has **four simple limit cycles in one nest**.

Setting: the Cherkas–Artés–Llibre normal form

```
xdot = 1 + x y
ydot = a00 + a10 x + a20 x^2 + a01 y + a11 x y + a y^2 ,   a00 = a01 + a11 - a10 - a20 - a
```

focus `A = (1,-1)`; section = a ray of the line `y = -1` from `A` (`x > 1`, or `x < 1`);
displacement `D(x0)` = (x-coordinate of the first return) − `x0`.

* cusp (triple cycle): `D = D_x = D_xx = 0` — 3 equations, so a curve in
  `(a11, a01, a10, x0)` at fixed shape `(a, a20)`;
* swallow-tail (quadruple cycle): additionally `D_xxx = 0` — codimension 1 *inside*
  the cusp manifold, detectable as a **sign change of `D_xxx`** along a cusp curve.

---

## 1. The engine (`cusp_engine.cpp`)

Taylor-series time-stepping (the field is polynomial, so the series recurrences are
exact), with the state carried as a **degree-4 jet in the section coordinate**
`eps`: `x(0) = x0 + eps`, `y(0) = -1`. The jet of the *return time* is solved from the
jet equation `y(tau) = -1`, and `R(eps) = x(tau(eps))` is then a degree-4 jet, so

```
D = R.c0 - x0 ,  D_x = R.c1 - 1 ,  D_xx = 2 R.c2 ,  D_xxx = 6 R.c3 ,  D_xxxx = 24 R.c4
```

come out **exactly** (to integration accuracy) — *not* by finite differencing. This is
what makes an honest `D_xxx` possible. Templated on `long double` and `__float128`;
all production runs are binary128 (113-bit mantissa), Taylor order 26, local tolerance
1e-32.

**Precision established, not assumed.** Re-running the same point at Taylor order 22/26/30,
tolerance 1e-30/1e-32/1e-34 and `hmax` 0.05/0.10/0.25 changes `D` in the **22nd
significant digit** (absolute agreement ~4e-34).

**Second, independent engine** (PROTOCOL rule 2): `indep_engine.py`, an mpmath Taylor
integrator written from the ODE directly, sharing no code with the C++ engine and using a
different event solver. Agreement at three test points:

| point | relative difference |
|---|---|
| `(a,a20,a11,a01,a10)=(3,-12,-1.398,8.4,15.28)`, `x0=1.28` | 1.3e-29 |
| `(3,-12,-2,9,13.5)`, `x0=1.02` (`D ~ -1.78e-12`) | 2.4e-22 (absolute 4e-34) |
| `(5,-50,-5.49995,16.5,76.45)`, `x0=1.3` | 1.4e-30 |

**A representational trap worth recording.** With the parameters carried as Python
floats the cusp Newton stagnates at a residual of ~1e-17 — not for any numerical-analysis
reason, but because the Newton steps are of relative size ~1e-25 in `(a11, a01, a10)` and
a binary64 parameter cannot represent them. All driver arithmetic is mpmath `mpf` at
dps = 50. After the fix the same Newton reaches residual **5.8e-34**.

---

## 2. Entering the cusp manifold from the Bautin region (TASK 2) — DONE

At a third-order weak focus (`a11 = 4-2a`, `a01 = 2a+1-a11`, `a10` as in Cherkas eq.(16))
the displacement satisfies `D(r) ~ d7 r^7`. Measured at `(a,a20) = (3,-12)`:
`D/r^7 = -1.537, -1.392, -1.153` at `r = 0.01, 0.02, 0.04` — the expected 7th-order
contact, so `d7 ~ -1.5`.

Seeding Newton on `(D, D_x, D_xx) = 0` in `(a11, a01, a10)` at fixed `x0 = 1 + r0`
from the weak-focus point (the Bautin shift is `O(r0^2)` in `V5`, `O(r0^4)` in `V3`,
`O(r0^6)` in `V1`) converges in 5 iterations to residual 5.8e-34 for every shape tried.
Example, `(a, a20) = (3, -12)`, `r0 = 0.02`:

```
a11 = -1.998494677833741331488579457852263
a01 =  8.998494677908351789196814624248964
a10 = 13.504422527141653959512639631905580
D = -5.8e-34   D_x = -3.1e-33   D_xx = -1.3e-32   D_xxx = -9.4920440228e-6
```

`D_xxx` at entry agrees with the Bautin prediction `48 r0^4 d7` (`-1.15e-5` with the
measured `d7 ~ -1.5`) in sign and order of magnitude, as it must.

### The triple cycle is real at NORMAL amplitude

Taking the cusp point at `x0 = 2.2654` on the `(a, a20) = (3, -12)` curve and perturbing
into the cuspidal region (`delta_mu ~ 1e-3`, chosen so that `c1 * D_xxx < 0`) gives
**three certified sign changes of `D`** in one nest:

| root | bracket `min|D|` | two-tolerance noise | margin |
|---|---|---|---|
| 2.12940787703450780729514 | 2.19e-8 | 1.07e-11 | 2050× |
| 2.21812962597053274330530 | 2.46e-8 | 1.11e-11 | 2220× |
| 2.27153111659064560853597 | 2.39e-8 | 1.14e-11 | 2100× |

(PROTOCOL rule 1 satisfied; noise estimated by recomputation at looser tolerance,
`noise = 10|difference| + 5e-12 s`.) The cycles sit at `x ≈ 2.13–2.27`, i.e. they are of
normal size, not small-amplitude. `triple_confirm_row1.json`.

---

## 3. Continuation of the cusp curve (TASK 3) — IN PROGRESS

Pseudo-arclength continuation in `(a11, a01, a10, x0)` at fixed `(a, a20)`, with a
chord Newton (frozen 3×4 Jacobian, refreshed every 4 iterations) for speed and an exact
Jacobian recomputed at any candidate event. Every accepted point logs `D, D_x, D_xx,
D_xxx, D_xxxx`, `nu = D_xxx/(D_xxxx r0)` (a scale-free distance to a swallow-tail),
`V1`, `L = det J(A)`, the return time, the transversality of the section crossing, and
Perko's nondegeneracy Jacobians with `mu1 = a11` (Cherkas's *rotating* parameter, so
`d_{mu1} != 0` is guaranteed by Duff/Perko monotonicity).

### Cherkas shapes, rows 1–4 (stopped at budget, not at a curve end)

| shape | `(a, a20)` | pts | `x0` range | `D_xxx` start → end | sign changes |
|---|---|---|---|---|---|
| row1 | (3, −12) | 120 | 1.02 → 10.10 | −9.49e−6 → −2.84e−3 | **0** |
| row2 | (1.5, −15) | 183 | 1.02 → 13.41 | +1.07e−7 → +2.05e−4 | **0** |
| row3 | (−2, 12) | 161 | 1.02 → 1.394 | +8.98e−6 → +1.78e+3 | **0** |
| row4 | (−2, −1) | 151 | 1.02 → 1.389 | +3.38e−5 → +2.33e+3 | **0** |

Two clearly different behaviours:

* **rows 1 and 2** — the curve runs out in amplitude without bound (`x0` past 10 and 13),
  `|D_xxx|` rises to a maximum (−5.9e−2 at `x0 ≈ 2.27` for row 1; +3.2e−3 at `x0 ≈ 2.87`
  for row 2) and then decays monotonically toward zero **without changing sign**. `L`
  stays positive (antisaddle preserved) and the section stays strongly transversal
  throughout. These runs were stopped by budget, not by a curve end.
* **rows 3 and 4** — `x0` stalls near 1.39 while `a11` and `|D_xxx|` blow up. The curve
  is running into a boundary in parameter space rather than in amplitude.

**No `D_xxx` sign change on any of the four.** Rows 5–8 and the (a, a20) grid are running.

### The grid (running)

`grid.py` records the structural fact that organises the search. On the third-order
stratum,

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

All eight Cherkas shapes sit on one side of the centre curve or in a region where it is
unreachable — which is consistent with all four completed curves showing a constant
`D_xxx` sign. The running grid straddles the centre curve at 16 values of `a` with
offsets ±0.08, ±0.3, ±1.0, plus a band sweep at the 8 unreachable `a` values: 116 grid
jobs + rows 5–8.

### A negative result worth recording: no swallow-tail at small amplitude

Seeding the square swallow-tail system `(D, D_x, D_xx, D_xxx) = 0` in
`(a11, a01, a10, a20)` at fixed `(a, x0)` from a *small-amplitude* cusp point
(`a = 3`, `x0 = 1.02`) makes the Newton run `a20 → -oo` (it reached −4.6e5 with `a10`
→ +6.3e5 and the residual still decreasing like 1/|a20|). This is not a solver failure:
at small amplitude `D_xxx = 48 d7 r0^4`, and `d7 = 0` only on the centre variety, where
`D` vanishes identically. **A nondegenerate swallow-tail cannot live at the
small-amplitude end of the Bautin unfolding.** It must be sought at normal amplitude —
which is exactly what this lane is set up to do.

---

## 4. Engine validation (PROTOCOL rule 7) — re-running

First pass (right-hand section only) reproduced **3 cycles in the nest** at the published
parameters for rows 1, 2, 6, 7, 8, with positions ~1.5 % from the published values.
The discrepancy is explained, not waved away: the 3-cycle window in `a11` is narrower
than 1e-3, and the paper quotes `a11` to 4–6 digits. At `(a,a20)=(3,-12)`, `a01=8.4`,
`a10=15.28`:

| `a11` | certified cycles |
|---|---|
| −1.400 | 0 (V1 = 0 exactly — the Hopf) |
| −1.399 | 1 |
| **−1.398 (published)** | **3** at 1.2809, 2.0070, 4.0193 |
| −1.397 | 0 |
| −1.396 | 0 |

So the engine reproduces the delicate three-cycle structure at exactly the published
`a11` and nowhere else — a sharper validation than matching the rounded positions.
Rows 3, 4, 5 have their published crossings at `x < 1`, i.e. on the **other** ray of
`y = -1` from the focus; they are being re-run with the left section, together with
row 4's published Andronov–Hopf polynomial check (2 extrema on [0.6, 0.9]) and the
fold check (`D_xx` from the jet vs a centred difference of `D_x`).

---

## 5. Open problems / next steps

1. Finish the (a, a20) grid and tabulate `sgn D_xxx` at the near and far ends of every
   cusp curve. A shape where the two differ contains a swallow-tail.
2. Rows 1 and 2 show `|D_xxx| → 0` as `x0 → ∞` **without a sign change**. Determine
   whether that decay is asymptotic (no swallow-tail, and the cusp curve simply escapes
   to infinite amplitude) or whether the curve ends at a graphic first. `nu =
   D_xxx/(D_xxxx r0)` is now logged precisely to tell these apart.
3. Rows 3 and 4 end with `a11` blowing up at nearly constant `x0`. Identify the
   boundary: a graphic, a loss of the antisaddle, or a chart failure.
4. Run the swallow-tail Newton (with `a20` free) from far-amplitude cusp points across
   the whole grid, not just from the small-amplitude entry where it provably cannot work.
5. Log the Andronov–Hopf/`beta*` picture at each cusp point for Lane 1 (implemented in
   `probe.ah_sweep`; not yet run over the grid).

---

## 6. Files

| file | what |
|---|---|
| `lane2_cusp_2026_09_06/cusp_engine.cpp` | the Taylor-jet return-map engine (binary128 / long double, jet degree 3 or 4) |
| `lane2_cusp_2026_09_06/engine.py` | mpmath-precision driver |
| `lane2_cusp_2026_09_06/indep_engine.py` | independent mpmath integrator (rule 2 cross-check) |
| `lane2_cusp_2026_09_06/cusp.py` | cusp residual, Jacobians, Newton, pseudo-arclength |
| `lane2_cusp_2026_09_06/continue_cusp.py` | the continuation driver + JSONL ledger |
| `lane2_cusp_2026_09_06/swallow.py` | square swallow-tail Newton + Perko Thm 4.3 quantities |
| `lane2_cusp_2026_09_06/probe.py` | rule-1 cycle counting, cusp unfolding, `AH(x)` |
| `lane2_cusp_2026_09_06/grid.py`, `make_grid_spec.py` | the (a, a20) grid and the centre curve |
| `lane2_cusp_2026_09_06/validate.py` | PROTOCOL rule 7 validation |
| `lane2_cusp_2026_09_06/analyse.py` | ledger summariser |
| `lane2_cusp_2026_09_06/ledger/`, `ledger_grid/` | append-only JSONL ledgers |
