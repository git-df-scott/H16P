# Finite-amplitude cusp geometry and the fifth cycle

**Date 2026-09-06. Branch `opus/rank-repair-cusp-compatibility-2026-09-06`.**

Source commits actually read:

| branch | commit |
|---|---|
| `main` | `45f4ea9b4ab448bdd36702036244faa4f15c9819` |
| `opus/degeneracy-collision-2026-09-06` | `ca4fe9062cb91b2974a99156bbcc88571f17b4d1` |
| `fable/lane2-cusp` | `e3a617418f8e851b8a26a4e88e2db6fc8e91bced` |
| `fable/coordination-2026-09-06` | `4482dfdced185a3a49fb24c702b116606d44476b` |
| `astra/fastra-afternoon-2026-09-05` | `55537aeb9613a2df05c0b26fa1b0d8921f9ab8bb` |

Lane 2's engines and ledgers were copied into this branch and **used as they
are**. Their validation was not re-run and no broad search was restarted. New
engine calls are logged in `lane2_cusp_2026_09_06/ledger_opus/`.

**NO FIVE-CYCLE FIELD. NO QUADRUPLE-CYCLE CANDIDATE.**

## 1. The one `D_xxx` sign change in the saved data is not a swallow-tail

Mining all 15 saved cusp ledgers, exactly one shows a sign change of
`G = D_xxx`: `ledger_grid/cusp_c_am2p0_om0p08.jsonl`, shape
`(a,a20) = (-2, -4.821)`, which sits `0.08` below the centre curve. It was still
running when `REPORT_lane2.md` was written and appears unanalysed. It is
**rejected**, for a reason visible in the ledger itself:

| quantity | interpolated zero at |
|---|---|
| `D_xxx` | `x0 = 1.03349338` |
| `D_xxxx` | `x0 = 1.03359934` (separation `1.1e-4`) |
| `V1` | `x0 = 1.03330651` (separation `1.9e-4`) |

`D_xxxx` vanishes with `D_xxx`, so **Perko's nondegeneracy `D_ssss != 0` fails**;
and `nu = D_xxx/(D_xxxx r0)` runs `0.1225, 0.1232, 0.1234, 0.1224` straight
through the crossing without approaching zero. This is the whole jet changing
sign together — `V1` crosses zero at the same place — not a multiplicity-four
cycle.

**Consequence for the search strategy.** `REPORT_lane2.md` open problem 1 says
"a shape where the two ends differ contains a swallow-tail". That criterion is
**not sufficient**: this is a counterexample to it. The correct detector is a
zero of `nu`, which Lane 2 already logs and which is exactly the invariant that
separates the two cases.

## 2. `nu` on every saved curve, and why it goes to zero

`nu > 0` on every curve that logs it, with no sign change anywhere:

| ledger | `(a, a20)` | `nu` range | `min|nu|` | location of the min |
|---|---|---|---|---|
| `cusp_row5` | `(-4, -1)` | `[0.00144, 0.0944]` | **0.00144** | last logged point |
| `cusp_c_am2p9_o1p0` | `(-2.9, -4.889)` | `[0.00250, 0.1001]` | 0.00250 | last logged point |
| `cusp_c_am2p5_o0p08` | `(-2.5, -5.344)` | `[0.00320, 0.1243]` | 0.00320 | last logged point |
| `cusp_c_am2p9_o0p08` | `(-2.9, -5.809)` | `[0.00958, 0.1325]` | 0.00958 | last logged point |
| `cusp_row6/7/8`, `cusp_c_am2p0_*` | | `[0.0805, 0.366]` | 0.0805 | — |

Every minimum is at the **last logged point**, and there the amplitude has
stalled (`dx0 ~ 3e-4` per step and shrinking) while the coefficient norm grows
linearly (`d|coef| ~ 0.145` per arclength step). Fitting over the last 30 points,

\[
 \nu \sim \lVert(a_{11},a_{01},a_{10})\rVert^{-5.3\ \text{to}\ -7.1}.
\]

**`nu -> 0` is an escape, not an approach to a finite-amplitude swallow-tail.**

### It is a genuine escape, not a coordinate artefact

The Cherkas form `xdot = 1+xy` is preserved by
`(x,y,t) -> (\alpha x,\ y/\alpha,\ \alpha t)`, under which

\[
 a\mapsto a,\quad a_{01}\mapsto\alpha a_{01},\quad a_{11}\mapsto\alpha^2a_{11},
 \quad a_{10}\mapsto\alpha^3a_{10},\quad a_{20}\mapsto\alpha^4a_{20}.
\]

So `a` is a genuine modulus, and `I_1=a_{11}^2/a_{20}`, `I_2=a_{01}^2/a_{11}`,
`I_3=a_{10}/(a_{01}a_{11})` are scale invariants. Since `D\sim\alpha`,
`D_{xxx}\sim\alpha^{-2}`, `D_{xxxx}\sim\alpha^{-3}`, `r_0\sim\alpha`, the
indicator **`nu` is itself scale-invariant** — so a small `nu` cannot be
manufactured by rescaling, and equally the escape cannot be rescaled away.
Measured in invariants along the three lowest-`nu` curves, `I_1` grows by
`3.4x` to `4.8x` while `nu` falls by up to `65x`:

\[
 \nu \sim |I_1|^{-3.4\ \text{to}\ -4.5}.
\]

**`nu -> 0` only as the invariant shape coordinate `I_1` diverges** — on the
boundary of shape space, not at an interior point. This is the sharp form of the
obstruction and it is the same obstruction in both descent directions (below).

## 3. Constrained derivatives across the cusp locus

Implemented as specified, with linear solves and no explicit inverse
(`constrained.py`). With `u=(a_{11},a_{01},a_{10})`, `q=(a,a_{20},x_0)`,
`F=(D,D_x,D_{xx})`, `G=D_{xxx}`:

\[
 \frac{dG}{dq}=G_q-G_uF_u^{-1}F_q,
\]

and the same for `H=D_{xxxx}`, then `nu` by the quotient rule.

**Validated** (`validate_constrained.py`): the constrained `dG/dx_0` reproduces
the slope of `D_xxx` measured along the saved curves to `2.9e-4` to `1.2e-3`
relative at well-resolved points.

**Honest limitation:** `dnu/dx_0` needs `D_{xxxxx}`, which the degree-4 jet does
not supply. That component is reported as unavailable, **not** as zero. The
shape components `dnu/da`, `dnu/da_{20}` are complete and are what is used below.

### The descent direction is `-a`, and it escapes too

Across 21 saved cusp points on 7 shapes, the steepest-descent direction of `nu`
in shape space is `(-1, \epsilon)` with `|\epsilon| < 0.04` at almost every point:
`|dnu/da|` exceeds `|dnu/da_{20}|` by one to four orders of magnitude. Decreasing
`a` decreases `nu`.

Testing that direction directly at fixed amplitude `r_0=0.05`, `a_{20}=-1`
(far from the centre curve), entering the cusp locus at each `a`:

| `a` | `-2` | `-3` | `-4` | `-6` | `-8` | `-12` |
|---|---|---|---|---|---|---|
| `nu` | 0.0949 | 0.0892 | 0.0827 | 0.0670 | 0.0472 | entry fails |
| `|coef|` | 15.0 | 20.1 | 25.1 | 35.4 | 45.9 | (guard bound) |

`nu` falls by only `2x` while the coefficient norm triples, and `I_1 = -(4-2a)^2/a_{20}`
diverges like `a^2`. **The shape direction is the same escape**, and it is much
slower than the amplitude direction.

## 4. The compatibility result — this is the substantive finding

A fifth cycle needs a second nest. Classifying every finite equilibrium
(`y=-1/x`, `a_{20}x^4+a_{10}x^3+(a_{00}-a_{11})x^2-a_{01}x+a=0`) along each cusp
curve:

| shape | second **focus**? | `nu` along the curve |
|---|---|---|
| `row5` `(-4,-1)` | **no** (2 saddles + node) | `0.094 -> 0.0014` |
| `row6` `(5,-50)` | **no** (2 nodes + saddle) | `0.182 -> 0.280` |
| **`row7` `(0.727,-12)`** | **yes**, persists (tr `-4.83 -> -6.18`) | `0.173 -> 0.300` |
| **`row8` `(1.04,-120)`** | **yes**, persists (tr `+0.97 -> +0.92`) | `0.132 -> 0.174` |
| `am2p5`, `am2p9`, `am2p0` grid | **no** (saddles/nodes) | `0.123 -> 0.0032` |

> **The two requirements pull in opposite directions.** The shapes whose `nu`
> falls toward zero have **no second focus at all**, so no fifth cycle is
> available there however good the fourth becomes. The two shapes that do carry
> a persistent second focus have `nu ~ 0.13` to `0.30` and **increasing** along
> the cusp curve — continuation moves them away from a swallow-tail, not toward
> one.

A direct sweep of `a` at fixed amplitude `r_0=0.05` confirms this is not an
artefact of the four particular curves:

| `a20` | `a=2.0` | `a=1.5` | `a=1.04` | `a=0.727` |
|---|---|---|---|---|
| `-12`: `nu` | 0.393 | 0.1083 | 0.1065 | 0.1053 |
| `-12`: second focus | no (node) | yes | yes | yes |
| `-120`: `nu` | 0.388 | 0.1077 | 0.1062 | 0.1052 |
| `-120`: second focus | no (node) | yes | yes | yes |

`nu` is **flat at `0.105`–`0.108` across the entire second-focus window** and
almost independent of `a_{20}` (a factor of ten in `a_{20}` moves `nu` by `0.1%`).

The `a = 2.0` column is degenerate and carries no information: the third-order
stratum has `a_{11} = 4-2a`, which vanishes at `a = 2`, so `I_1 = a_{11}^2/a_{20}`
is `~1e-55` there. Its `nu ~ 0.39` and its lack of a second focus are artefacts
of that degeneracy, not evidence about large `a`.

The admissible window is set by `(a-3-a_{20})/(1-3a) < 0` and **flips side with
`a_{20}`**: at `a_{20}=-12,-120` it is `a > 1/3`, at `a_{20}=-1` it is `a < 1/3`.
(My first attempt to complete the `a_{20}=-1` row swept `a > 1/3` and every entry
failed on a guard; the window is the other side. Recorded rather than quietly
re-run.) Completing it on the correct side:

| `a20 = -1` | `a = 0.3` | `a = 0.0` |
|---|---|---|
| `nu` | 0.1038 | 0.1032 |
| second focus | no (node) | no (node + saddle) |

So the third `a_{20}` value gives the same flat `nu ~ 0.103` and, where entry
succeeds at all, again **no second focus**. It does not change the conclusion.

## 5. What this does and does not establish

Established, on the data examined:

1. The single `D_xxx` sign change in the saved ledgers is not a swallow-tail
   (`D_ssss` vanishes with it), and the sign-change criterion is not sufficient.
2. `nu > 0` on every saved cusp curve; every minimum is at a stopped endpoint.
3. `nu -> 0` only as the scale-invariant shape coordinate `I_1` diverges, in
   both descent directions, with measured power laws.
4. In the region where a second focus exists, `nu ~ 0.105` and increases along
   the cusp curves.

**Not established.** This is not an exclusion of the swallow-tail. Specifically:

- Only `a_{20} \in \{-1,-12,-120\}` and `a` on a coarse grid were swept, at one
  amplitude `r_0 = 0.05`; the invariant shape space `(a, I_1, I_2, I_3)` is
  four-dimensional and was not covered.
- The curves stopped by Lane 2's budget were not continued further; their
  asymptotics are inferred from power-law fits over the last 30 points, not
  proved.
- `dnu/dx_0` is unavailable without `D_{xxxxx}`, so the constrained analysis is
  complete only in the shape directions.
- Nothing here bears on mechanism (b) of `PROTOCOL.md` (a neutral graphic at the
  nest boundary around a **strong** focus), which the protocol lists as open.
- No claim is made about `a \le 1/3` at these `a_{20}`, where the chart fails.

## 6. Next step this points to

The decisive question is now sharp and cheap to state: **is `nu` bounded below
by a positive constant on compact subsets of the invariant shape space
`(a, I_1, I_2, I_3)`?** Everything measured is consistent with yes. A single
counterexample — one interior point with `nu` small and a second focus present —
would be a quadruple-cycle candidate. Extending the `D_{xxxxx}` jet would also
close the one gap in the constrained derivative and let `nu` be continued
directly rather than inferred.

## 7. Replay

```bash
cd lane2_cusp_2026_09_06
python3 mine_ledgers.py          # all 15 ledgers, D_xxx sign changes
python3 inspect_signchange.py    # the one sign change, in full
python3 adjudicate.py            # its rejection: D_xxxx and V1 vanish with D_xxx
python3 nu_scan.py               # nu on every saved curve
python3 nu_profile.py            # minima are at stopped endpoints
python3 escape_check.py          # nu ~ |coef|^(-5..-7)
python3 invariants.py            # scale invariants; the escape is genuine
python3 constrained.py           # dG/dq, dnu/dq  (linear solves)
python3 validate_constrained.py  # dG/dx0 vs the ledger slope
python3 shape_gradient.py        # steepest descent is -a
python3 a_sweep.py               # the -a direction also escapes
python3 remote_nest.py           # second-nest classification (no engine calls)
python3 compatibility.py         # nu vs second focus across shapes
```

New engine calls are appended to `lane2_cusp_2026_09_06/ledger_opus/*.jsonl`
with the engine banner and running call count. The `4096/4096` figure in the
campaign report is the **KKL/Shi** ledger on `main`; Lane 2 operates under
`PROTOCOL.md`, which imposes ledger discipline rather than a numeric cap, and
its own cusp ledgers already hold ~1500 continuation points.
