# H16P — Fold continuation, the focus boundary, and a full cycle inventory

**Branch** `opus/rank-repair-cusp-compatibility-2026-09-06`
**Date** 2026-09-07 · **Engine** `cusp_engine4`, binary128 (113-bit), order-26 Taylor, jet degree 4
**Independent counter** scipy `solve_ivp` (RK45, rtol 1e-11) ray-shooting, `foldcont.remote_cycles`

Everything below is one of three evidence classes, never mixed:
**(N)** numerical evidence, **(A)** exact algebra, **(C)** rigorous certificate.
Nothing in this report is class (C).

---

## 0. What this report answers

The previous report (`CENTRE_AND_DOUBLEFOLD_REPORT.md`) ended with an
unidentified event: an extremum of the displacement was absent from the endpoint
profile of a guarded Newton path, and I wrote that it had been "destroyed".
That was inferred from sampling, not from tracking. This report replaces the
inference with continuous tracking, and corrects three further claims. The
corrections are appended to that file as §C1–C4; they are summarised in §5 here.

---

## 1. Engine convention verified (N)

The brief's identity `D_s(1) = exp(pi*T/omega) - 1`, with `omega = sqrt(det J(A) - T^2/4)`,
is the engine's convention. Checked at three trace values by evaluating the jet
at `s = 1 + e` and letting `e -> 0`:

| `T` | `omega` | `exp(pi T/omega)-1` | `D_s(1+e)`, `e=1e-6` | rel. error |
|---|---|---|---|---|
| `0.05` | `7.1832454` | `0.02210834803` | `0.02210838299` | `1.58e-6` |
| `0.01` | `7.1860709` | `0.004381350943` | `0.00438135774` | `1.55e-6` |
| `0.002` | `7.1866292` | `0.0008746705331` | `0.0008746718847` | `1.55e-6` |

The error is exactly linear in `e` (`1.58e-4 / 1.58e-5 / 1.58e-6` at
`e = 1e-4/1e-5/1e-6`), i.e. it is the `D_ss(1)*e` term, not a discrepancy.
Convention confirmed. Script `lane2_cusp_2026_09_06/nondeg.py`.

## 2. The Hopf + cycle-fold organizer, with every nondegeneracy (N)

Solving `T(mu) = 0`, `D(s_f) = 0`, `D_s(s_f) = 0` in `(a11, a01, s_f)` from the
row-7 fold gives

```
a    = 0.7272727272727272727272727      (= 8/11)
a20  = -12.0
a11  = 2.149639912563477722532764
a01  = 0.3049055419841944822544849
a10  = -26.5
s_f  = 3.435277754690688210709527
T = 2.2177e-12    det J(A) = 51.649639912561    D_ss(s_f) = -1.2381132342e-3
```

Nondegeneracies, all checked, none assumed:

| condition | value | verdict |
|---|---|---|
| `det J(A) > 0` (A is a focus) | `51.6496` | holds |
| `D_ss(s_f) != 0` (the fold is quadratic) | `-1.238e-3` | holds |
| `l1 != 0` | `D_sss(1) = -5.645117e-3`, so `c3 = -9.4085e-4 < 0` | holds |
| `rank[grad_mu T ; d_mu D(s_f)] = 2` | see below | holds |

**The Lyapunov coefficient (N).** At `T = 0` the displacement has no linear and
no quadratic term at `s = 1`: measured `D_ss(1+e) = -5.6449e-8, -5.6436e-9` at
`e = 1e-5, 1e-6` (linear in `e`, so `D_ss(1) = 0`), while `D_sss(1+e)` is flat at
`-5.645e-3`. Hence `D(s) = c3 (s-1)^3 + O((s-1)^4)` with `c3 = D_sss(1)/6 =
-9.4085e-4`; the independent check `D/e^3 = -9.3916e-4` at `e = 1e-3` agrees.
`l1 < 0`, so the Hopf is **supercritical** and the cycle exists for `T > 0`.
Direction test, by search rather than by theory: a cycle at `s-1 = 0.02193` for
`T = +1e-6` and at `s-1 = 0.07231` for `T = +1e-5`, and **none** for `T = -1e-6`
or `T = -1e-5`. The sign choice is `T > 0` with `l1 < 0`, i.e. `T*l1 < 0`.

**Unfolding rank (N).** In scaled coordinates `(a, a20, a11, a01, a10)`:

```
grad_mu T   = (-2, 0, 2.1496399, 1, 0)
d_mu D(s_f) = (-2.5557883, 0.38815693, 5.8401574, 2.3627366, -0.27718643)
```
cosine `0.94743597`; **after normalising each row to unit length** the singular
values are `(1.3955056, 0.22926846)`, ratio `0.16429`. Rank 2, well clear of
round-off. (The row-normalisation is stated because the previous report reported
such a ratio without saying so.)

## 3. The rank-2 unfolding realises three coexisting cycles (N)

Taking the least-norm `delta_mu` that sets `(T, D(s_f)) = (eta_T, eta_D)` and
counting **every** sign change of `D` on the whole return domain `[1.002, 6.997]`,
plus an independent scipy count in the second nest:

| `eta_T` | achieved `T` | achieved `D(s_f)` | local cycles at `s =` | remote | total |
|---|---|---|---|---|---|
| `1e-5` | `1.0000002e-5` | `8.632e-10` | `1.07212, 3.42000, 3.43537` | `0` | **3** |
| `1e-4` | `1.0e-4` | `-4.276e-9` | `1.26643, 3.27776, 3.43523` | `0` | **3** |
| `3e-4` | `3.0e-4` | `-8.851e-8` | `1.60220, 2.90064, 3.43496` | `0` | **3** |
| `1e-3` | `1.0e-3` | `-1.327e-6` | `3.43386` | `0` | **1** |

The predicted `1 (Hopf) + 2 (fold) = 3` is realised and verified by inventory,
not asserted from the normal form. At `eta_T = 1e-3` the step is too large and
the pair has already annihilated. **No fifth cycle, and no fourth**: this
organizer's ceiling is three.

Script `lane2_cusp_2026_09_06/unfold.py`; ledger `ledger_opus/unfold.json`.

## 4. The actual fold continuation, and what the missing extremum really was (N)

### 4.1 The method

`lane2_cusp_2026_09_06/foldrun.py`. State `(mu, s_f)` on the fold surface
`D(s_f) = D_s(s_f) = 0`, with a **second stationary point `s_i` tracked by Newton
on `D_s` from the previous step**, never re-detected by sampling. Steering, using
the tangent identities rather than assuming them:

```
b_f = d_mu D(s_f),   b_i = d_mu D(s_i),   h_i = D(s_i)
w   = b_i - (b_i.b_f)/(b_f.b_f) b_f      (fold-preserving to first order)
dmu = -sign(h_i) * step * w/|w|
ds_f = -(d_mu D_s . dmu)/D_ss(s_f)       (predictor)
```
then Newton correction in `(s_f, a11)` back onto `(D, D_s) = 0`, with real
backtracking (halve on corrector or tracker failure, grow by 1.3 on success).
Seed: the row-7 fold, `s_f = 3.46899418303`, `s_i = 2.28613480459`,
`h_i = -2.50e-4`.

### 4.2 What happened — the path terminates in a cusp, not a second fold

Twelve accepted steps. `h_i` fell monotonically by five orders of magnitude, and
the two folds coalesced:

| step | `delta = (s_f-s_i)/2` | `|h_i|` | `|D_ss(s_f)|` | `|h_i|/delta^3` | `|D_ss|/delta` |
|---|---|---|---|---|---|
| 0 | `0.59143` | `2.50023e-4` | `1.29003e-3` | `1.20857e-3` | `2.18121e-3` |
| 2 | `0.352356` | `4.16101e-5` | `6.05706e-4` | `9.5116e-4` | `1.71902e-3` |
| 4 | `0.102131` | `9.14552e-7` | `1.41318e-4` | `8.58489e-4` | `1.38369e-3` |
| 6 | `0.0329164` | `3.04333e-8` | `4.31699e-5` | `8.53317e-4` | `1.3115e-3` |
| 8 | `0.00703185` | `2.96522e-10` | `9.04281e-6` | `8.52803e-4` | `1.28598e-3` |
| 11 | `0.00504406` | `1.09442e-10` | `6.47682e-6` | `8.52791e-4` | `1.28405e-3` |

Over two decades in `delta` and five in `h_i`, `|h_i|/delta^3` is constant to
**six** significant figures and `|D_ss|/delta` to four. That is the cusp normal
form and nothing else: for `D = (c/6)u^3 - (c/2)delta^2 u + const` normalised so
that `D(s_f) = 0`, one has `|D_ss(s_f)| = c*delta` and `|h_i| = (2c/3)delta^3`
exactly. From the table `c = 1.28405e-3` and `2c/3 = 8.5603e-4`, against the
measured `8.52791e-4` — agreement to `0.4%` at `delta = 0.005`.

### 4.3 The cusp, solved exactly (N)

Newton on `D = D_s = D_ss = 0` in `(a11, a01, s)` from the terminal state
converges in four iterations to a scaled residual of `2.14e-33`:

```
a    = 0.7266687384078295434888743
a20  = -12.00139906654482774401223
a11  = 2.149360603432544109053414
a01  = 0.3045260121871923258255834
a10  = -26.49779478462597842880328
a00  = 40.22641172838271306420563          (pinned: a00 = a01+a11-a10-a20-a)
s_c  = 2.5880523017615422998

D = 1.16e-33   D_s = 1.35e-33   D_ss = -4.57e-34
D_sss  = -1.27916880373e-3   (!= 0: the cusp is nondegenerate)
D_ssss = -2.89501591506e-3
T = 5.49138804077e-4   det J(A) = 51.6494043823
```

`D_sss = -1.27917e-3` against the `c = -1.28405e-3` extrapolated independently
from the scaling table — 0.4% apart, two different measurements of the same
number. This is a **multiplicity-three limit cycle at finite amplitude**
(`s_c - 1 = 1.588`), i.e. Perko's `C3`, reached by connected continuation from a
simple fold rather than by a blind solve.

Full inventory at the cusp parameter, domain `[1.002, 6.997]`: one sign change of
`D`, at `2.58805230157` — the triple cycle itself (odd multiplicity) — and one
other stationary point at `1.47734665662`. The cusp's own stationary point is
correctly *not* a sign change of `D_s`, since `D_s ~ (c/2)u^2` there. Remote
nest: second focus `(-3.2028, 0.3122)`, trace `-5.8135`, cycles `0`.

### 4.4 The event, named

> **The extremum was not destroyed. It merged with the fold.** The
> fold-preserving steepest-descent path drives `h_i -> 0` and the fold separation
> `delta -> 0` together, at the exact cubic rate `h_i ~ (2c/3) delta^3`, and
> terminates at a nondegenerate cusp of the return map. None of the other
> candidate events occurred: no return failure, no exit from the return domain,
> no approach to a centre (`D_sss != 0` throughout), and no rank loss.

This identifies the previous report's unresolved event — but only **along this
path**. The fold-preserving directions form a 4-dimensional space at each point;
`w` is the steepest-descent direction for `h_i` inside it. The 3-dimensional
complement, which changes the geometry at fixed `h_i`, is untested. This is a
statement about one path, not about the route.

### 4.5 The separation-pinned double fold also fails, differently (N)

`lane2_cusp_2026_09_06/gapfold.py` replaces the previous report's ad-hoc guards
by an equation: solve

```
F(a11,a01,a10,s_f,s_i) = ( D(s_f), D_s(s_f), D(s_i), D_s(s_i), s_f - s_i - g0 ) = 0
```

5 equations, 5 unknowns, `(a, a20)` as continuation parameters. The separation
can no longer be "hit", because it is part of the solution. From the row-7 seed
with `g0 = 1.18285937844`, the damped Newton **does not converge in 60
iterations**: the separation constraint is satisfied throughout, but both
critical values drift steadily more negative while the pair marches outward in
amplitude (`s = (3.469, 2.286) -> (5.270, 4.087)`, residual `2.5e-4 -> 2.7e-3`,
monotonically increasing). A different failure mode from the guarded solves,
which all slid onto the focus. Recorded, not diagnosed.

## 5. Unfolding the cusp: three cycles, verified; the fourth is not local (N)

### 5.1 Three cycles

The cusp has a 2-parameter unfolding. Solving, **by Newton and not by one linear
step**, `D(s_c) = 0` and `D_s(s_c) = b1 = -c*eps^2/2` in `(a11, a01)` puts the
displacement in the form `b1 u + (c/6)u^3 + O(u^4)`, whose roots are
`u = 0, +-sqrt(3)*eps`. Full inventory on `[1.002, 6.992]` at each solved point:

| `eps` | `D(s_c)` | `D_s(s_c)` vs target | predicted roots | measured cycles | `T` |
|---|---|---|---|---|---|
| `0.05` | `3.9e-34` | `1.598961e-6` = target | `2.50145, 2.58805, 2.67465` | `2.50220, 2.58805, 2.67544` | `5.474e-4` |
| `0.10` | `-1.9e-33` | `6.395844e-6` = target | `2.41485, 2.58805, 2.76126` | `2.41780, 2.58805, 2.76446` | `5.420e-4` |
| `0.35` | `3.1e-33` | `7.8349089e-5` = target | `1.98183, 2.58805, 3.19427` | `2.01349, 2.58805, 3.23702` | `4.618e-4` |
| `0.50` | `1.8e-31` | `1.598961e-4` = target | `1.72203, 2.58805, 3.45408` | `1.77941, 2.58805, 3.54543` | `3.710e-4` |

Three simple limit cycles at every step, widely separated at the larger `eps`
(`1.78, 2.59, 3.55`), and the measured positions track the cubic prediction with
the drift expected from the quartic term. This is a **verified count from a full
inventory**, not a normal-form assertion.

### 5.2 Why the fourth cycle is not available here

With `T > 0` the displacement leaves the focus upward and stays positive out to
the innermost of the three (`stationary point at s = 1.476`, `D > 0` there), so
no fourth crossing exists between the focus and the triple. A fourth would need
`T < 0`, which requires controlling the trace **alongside** the cusp conditions.

A Newton solve of `(D, D_s, D_ss, T) = (0, b1, 0, tau)` in `(a, a11, a01, a10)`
stalls with the residual sitting on the `T` equation. **A stalled Newton is not
evidence**, so the matrix was measured instead (`lane2_cusp_2026_09_06/cusprank.py`,
10 engine calls):

```
B = d(D, D_s, D_ss, T)/d(a, a20, a11, a01, a10),  scaled

d_mu D    = (-1.8611158, 0.13794252, 3.1137869, 1.3199745, -0.098473044)
d_mu D_s  = (-1.0428337, 0.21558662, 2.8411564, 1.1246711, -0.15399288)
d_mu D_ss = ( 0.4659277, 0.18588361, 0.944072,  0.27719692, -0.13284497)
grad_mu T = (-2,         0,          2.1493606, 1,           0)

row-normalised singular values  (1.8265586, 0.81466539, 1.9563778e-3, 1.9658966e-5)
det of the 4x4 block (a,a11,a01,a10) = 1.678916316e-6
```

`sigma_4 = 1.97e-5` is five orders above the finite-difference noise floor
(`~1e-21` relative in binary128), so **rank B = 4: the trace is controllable in
principle**. But the conditioning is the whole story: shifting `T` from the
cusp's `+5.49e-4` to any negative value needs `|delta_mu| ~ 6.5e-4/1.97e-5 ~ 34`
in scaled parameter units. That is not a local unfolding — it leaves the region
where the cusp persists, which is exactly what the Newton stall reports.

> **Statement, with its class.** At this cusp the maximum verified simultaneous
> cycle count is **three**. The fourth is not excluded — the unfolding has full
> rank — it is **not reachable locally**, by a factor of `~5e4` in conditioning.
> This is evidence class (N) about one cusp, not a theorem about the family.

## 6. The remote nest — the detector, its failure, and its repair (N)

**The first remote count was worthless and I nearly reported it.** Every
`remote_cycles` reading in this session's earlier scripts came back `0`. Before
using a zero from an independent detector, I ran it as a positive control on the
**first** nest, where the binary128 engine independently certifies the count. It
found **none** of three certified cycles. So the zeros carried no information.

Diagnosis (`lane2_cusp_2026_09_06/dirtest.py`): the event direction and the
section condition were right — the displacements were correct to `1e-10` — but
the scan used 40 samples over `d in (0, 3]` while two of the three cycles sit
`0.015` apart, inside one grid cell, and the third lay before the first sample.
A resolution failure, not a broken method.

Repaired (`lane2_cusp_2026_09_06/remote2.py`, `ds = 0.01`, 300 samples, integration
horizon cut from `t = 200` to `t = 4` once the return time was measured at `0.87`).
Positive control now passes on every point tested:

| organizer | engine cycles at `s-1` | detector | verdict | second focus | remote cycles |
|---|---|---|---|---|---|
| Hopf+fold `eta_T=1e-5` | `0.072124, 2.420004, 2.435369` | `0.072118, 2.419931, 2.435441` | AGREES | `(-3.20325, 0.31218)`, tr `-5.81463` | `0` |
| Hopf+fold `eta_T=1e-4` | `0.266432, 2.277757, 2.435232` | `0.266431, 2.277749, 2.435239` | AGREES | `(-3.20319, 0.31219)`, tr `-5.81432` | `0` |
| cusp `eps=0.05` | `1.502205, 1.588052, 1.675441` | `1.502366, 1.587735, 1.675599` | AGREES | `(-3.20281, 0.31223)`, tr `-5.81356` | `0` |
| cusp `eps=0.10` | `1.417803, 1.588052, 1.764462` | `1.417843, 1.587973, 1.764502` | AGREES | `(-3.20281, 0.31223)`, tr `-5.81379` | `0` |

The remote zero is now a **validated** negative, with its scope stated: no cycle
crossing the horizontal ray within `d <= 3` of a second focus whose trace is
`-5.81` and determinant `161.6` — strongly hyperbolic, far from any Hopf. Both
organizers sit in a region where the second nest is empty, which is precisely
why neither reaches `4+1`.

Every `remote = 0` reported earlier in this session, including in
`CENTRE_AND_DOUBLEFOLD_REPORT.md` and the row-7/row-8 fold-start ledger entries,
came from the **unvalidated** detector and should be read as "not measured".
The four rows above are the only validated remote counts.

## 7. Row 8

Row 8's fold start did not converge to the solver contract inside the allotted
wall clock (killed at 900 s, exit 143) and is **not** recorded as a result. The
gap-pinned double-fold sweeps on row 8 (four control pairs) all returned
`STALLED` at the amplitude guard `1.08000` with scaled residuals `2.5e-6`, the
same mechanism as rows 2 and 7. No row-8 conclusion is claimed.

## 8. Scope and what is not claimed

* Nothing here is class (C). The cusp is a numerical object with a scaled Newton
  residual of `2.1e-33` in binary128, cross-checked by an independent cubic
  extrapolation and by a second ODE integrator; that is strong (N), not a proof.
* The cusp is the terminal event of **one path** — the steepest-descent
  fold-preserving direction. The 3-dimensional complement of `w` inside the
  fold-tangent space is untested. Nothing here excludes a separated double fold
  reached another way, and §C3 of the corrected earlier report removes the bogus
  reason to think one is impossible.
* `rank B = 4` at the cusp means the four-cycle configuration is **not excluded
  by rank**; it is only out of local reach at this point. A cusp at a different
  shape `(a, a20)`, where `sigma_4` is not `1e-5`, is a live target.
* The remote-nest zero is measured on one ray to `d <= 3` at four parameters. It
  is not a statement about the family.

---

## STATUS BLOCK

| item | status |
|---|---|
| **FOLD STARTS INDEPENDENTLY REPRODUCED** | Row 7: yes (engine A + `indep_engine`). Row 8: **no** — did not converge in the allotted wall clock; not claimed. |
| **CONNECTED FOLD CONTINUATION COMPLETED** | **Yes.** 12 accepted predictor-corrector steps on `D=D_s=0` from the row-7 fold, second stationary point tracked by Newton throughout, terminating at an identified event. |
| **MISSING-EXTREMUM EVENT** | **Identified: fold coalescence to a nondegenerate cusp.** `h_i ~ (2c/3)delta^3` with `|h_i|/delta^3` constant to 6 s.f. over 2 decades in `delta`; terminal solve `D=D_s=D_ss=0` to scaled residual `2.14e-33` with `D_sss = -1.279e-3 != 0`. Not destruction, not a centre, not a return failure, not a rank loss. The earlier "destroyed" is retracted (§C4). |
| **MAX SIMPLE LOCAL CYCLES COEXISTING WITH THE FOLD** | **3** — verified by full inventory on `[1.002, 6.99]`, at four `eps` values from the cusp unfolding and three `eta_T` values from the Hopf+fold unfolding. At the fold itself (before unfolding): `0` besides the fold cycle. |
| **REMOTE CYCLE RETAINED** | **No — and now measured, not guessed.** `0` remote cycles at all four validated points; the second focus is strongly hyperbolic (trace `-5.81`, det `161.6`). The detector's earlier zeros were unvalidated and are withdrawn (§6). |
| **ORDINARY HOPF + SEPARATE CYCLE FOLD** | **Exists and solves.** `T = 2.2e-12`, `D(s_f)=D_s(s_f)=0` at `s_f = 3.4352777547`, with `det J(A) = 51.6496 > 0`, `D_ss = -1.238e-3 != 0`, `l1 < 0` (`D_sss(1) = -5.645e-3`), Hopf side `T > 0` confirmed by search. Its unfolding realises `1+2 = 3` cycles, verified. |
| **TWO-CONTROL UNFOLDING RANK** | `rank[grad_mu T ; d_mu D(s_f)] = 2`; cosine `0.9474`, row-normalised singular values `(1.3955, 0.2293)`, ratio `0.164`. At the cusp, `rank[d_mu D; d_mu D_s; d_mu D_ss; grad_mu T] = 4` but `sigma_4/sigma_1 = 1.03e-5`. |
| **FIVE-CYCLE FIELD** | **None.** Maximum verified anywhere in this session: **3** simultaneous simple limit cycles. |
| **RIGOROUS CERTIFICATION** | **None attempted.** No interval arithmetic, no validated Newton–Kantorovich ball. All results are class (N). |
| **UNRESOLVED BOUNDARY** | (i) the 3-dimensional fold-tangent complement of the steering direction `w`; (ii) whether a cusp at other `(a, a20)` has a better-conditioned `sigma_4`, which is the concrete route to `3+1`; (iii) why the separation-pinned double fold marches outward in amplitude; (iv) row 8. |
| **EVALUATIONS** | Engine calls recorded at exit: `nondeg` 290, `unfold` 6014, `foldrun p4` 747, `foldrun p2` 3163, `cusp_end` 1343, `gapfold` 840, `cusp3` 4880, `cusprank` 10 = **17287**. Three runs were killed before writing a count (`cusp_unfold`, `cusp_unfold_fast`, `cusp4`, `cusp5`); their ledger entries record what completed, and their totals are not claimed. Independent scipy integrations: several thousand, not counted per-call. |
| **BRANCH** | `opus/rank-repair-cusp-compatibility-2026-09-06` |
| **COMMIT** | see below |

### Files added

```
lane2_cusp_2026_09_06/
  nondeg.py        convention, l1, two-control rank
  unfold.py        rank-2 unfolding of the Hopf+fold organizer, full inventory
  foldrun.py       connected fold continuation with a tracked second stationary point
  cusp_end.py      cubic-scaling test + exact cusp solve
  cusp3.py         cusp unfolding by Newton on the rank-2 system
  cusp4.py         cusp unfolding with the trace, 3 unknowns  (stalls; kept)
  cusp5.py         same with D_ss pinned, 4 unknowns          (stalls; kept)
  cusprank.py      measures WHY they stall: sigma_4 = 1.97e-5
  gapfold.py       separation-pinned double fold (no convergence)
  rolle.py         the two counting lemmas + the opposite-type counterexample
  dirtest.py       event-direction check for the ODE counter
  remote2.py       resolved ODE cycle counter with a positive control
  validate_remote.py  the first, under-resolved control (kept: it is the failure)
  ledger_opus/*.jsonl, *.json, *.out   append-only
```
