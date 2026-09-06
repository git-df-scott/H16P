# H16P campaign report — everything, all branches, 2026-09-03 to 2026-09-06

Compiled 2026-09-06 from a read of `main` and all sixteen other branches.
Scope: hunting a counterexample to `H(2)=4` — one real planar quadratic vector
field with at least five limit cycles.

## 0. Bottom line

**No five-cycle field. No five-cycle candidate. No interval certificate of
five.** The maximum cycle count produced by any method in this repository is
**four**, in the `(3,1)` distribution, which was already known in 1979–80. The
maximum count *certified* in one nest is **three**, achieved 2026-09-06.
`H(2)=4` is neither proved nor refuted.

The discipline about this is unusually good: across 147 reports, no file
overclaims, and several record their own withdrawn conclusions.

## 1. Scale

| | |
|---|---|
| Calendar time | 4 days, 2026-09-03 to 2026-09-06 |
| Branches | 17 (one orphan, no merge base) |
| Distinct commits | 155 |
| Distinct files ever tracked | 937 |
| Markdown reports | 147 |
| ODE return-map budget | 4096/4096 spent on the shared KKL/Shi ledger, exhausted |

**`main` is not the frontier.** It carries 32 commits and looks complete; the
live work is on the branches. Anyone reading only `main` — as I did initially —
will misjudge what has been done.

## 2. Who did what

- **Astra** — the `Q_4` strike series (1 through 7), the KKL fold-surface
  continuation and its closure, the D1 fold ledger, the D2 order-two loop law,
  the resonant joint strike, and the hostile audit of Proposition A.
- **Fable** — the overnight F-lane campaign (F3–F21), the compactified
  return-map engine that made it possible, and the 2026-09-06 coordination that
  launched Lanes 1–3.
- **Claude (earlier sessions)** — adversarial audits of Astra strikes 1–4,
  Lanes B and C, Routes 4a/4b, the codimension thought session, and the
  campaign zoom-out audit.
- **Opus 5 (this session)** — the coverage-gap audit, the order-3 neutrality
  identity, the `Q_3^R` first-order structure, the order-two graphic collision,
  and a verification pass over its own results. Standing corrected in §7.

## 3. What is proved exactly

| Result | Where | Status |
|---|---|---|
| **Theorem N.** No `Q_4` Abelian integral has five distinct zeros in the open annulus, for every `kappa>1`. Bounds: three on the strict lobe, four globally | `Q4_THEOREM_N.md` | Independently audited, verdict SOUND |
| Further `Q_4` exclusions: a three-zero bound at the limiting lift `a=1`, a two-zero bound on an explicit 3-dimensional subspace at every finite lift | `astra/q4-determinant` | Exact replays, no external review |
| **Proposition A (corrected).** In the zero-trace Shi chart with `eta_1=0`, zero divergence at a non-origin equilibrium implies the origin is a **centre**. Hence an order-two weak focus cannot coexist with a neutral non-origin saddle | `FASTRA_AFTERNOON_REPORT` | Exact; the literal converse in the original was false and was withdrawn |
| Exact `eta_2` in the Shi chart, factored into its three centre strata | `FABLE_D2_ORDER_TWO_LOOP.md` | Exact |
| The proposed analytic monic-quartic scalar Dulac certificate is impossible at a true KKL fold; global single-negative-band restriction on `1<=c<=8/5` | `KKL_FOLD_CLOSURE.md` | Exact |
| Two-centre reversible geometry excludes the finite saddle-loop mechanism; the `a=-2` unfolding chart was missing a direction and is repaired | `REVERSIBLE_RESEED` | Exact |
| Two compact cycles plus three hemicycle cycles at the double-centre base `(-1,1)` are incompatible | `RESONANT_JOINT` | Derived, no external review |
| On the order-3 stratum, the neutrality resultant of the two-saddle infinity graphic is `640*eta_3` | `ORDER3_GRAPHIC_NEUTRALITY.md` | Exact (this session) |

## 4. The organising principle the campaign found

Fable's night report states it, and it is the most important thing in the
repository:

> **Rigidity.** Every degenerate boundary object placed next to a degenerate
> focus collapsed onto an integrable stratum. A fifth cycle cannot come from
> stacking degeneracies in the field.

**Scope caution.** This is a statement about the configurations actually
examined — order-three loop, order-two finite loop, order-two infinity graphic
in two charts. It is *not* a proved universal theorem about all quadratic
graphics; non-elementary graphics and other configurations remain open, and the
campaign's own lists say so. The language below should be read at that scope.

Three independent instances, found separately:

1. **Order-three focus + neutral loop** — the loop exists only on the centre
   curve (Lane C).
2. **Order-two focus + neutral finite saddle loop** — collapses to a centre
   (Fable F11; proved exactly by Astra's Proposition A; the sign law
   `sigma * eta_2 < 0` is Conjecture D2).
3. **Order-two focus + neutral infinity graphic** — collapses to a centre
   (Astra at KKL `K=J=0`; and independently this session in the Shi chart,
   §7).

This is why every route has failed the same way, and it is a serious partial
explanation of *why* `H(2)=4` might be true.

**It also identifies the one surviving mechanism.** Rigidity kills *stacked*
degeneracies. A multiplicity-four limit cycle is a **single** degenerate
object, not a stack — so it is not touched by the principle. That is exactly
Lane 2's target, and it is the right lane.

## 5. Closed as counterexample routes

| Route | Killed by |
|---|---|
| Attack 1: five `Q_4` Abelian zeros | Theorem N |
| Route 4b: order-3 focus + cyclicity-two boundary graphic | Lane C, Routes 4a/4b, and the `640*eta_3` identity |
| Order-two focus + neutral loop (codim-4 organising centre) | F11 + Proposition A |
| Yu–Han reversible family | F13: the triple-zero element never has an interior zero; `(3,1)` is that mechanism's ceiling |
| Aliens at the `Q_4` resonant infinity graphic | F19: exactly three real cycles every time, 20k targeted trials, no alien |
| Neutral hemicycle as a source of three | F16/F18: emits two, never three; matches Marín–Villadelprat exactly |
| Attack 2 in its original form | No finite saddle exists in the stated box |
| Saddle loop as the `Q_4` difference | Lane B |

## 6. What is open

1. **Lane 2 — Perko's cusp manifold to a swallow-tail.** Continue the manifold
   of multiplicity-three cycles out of the Bautin region to normal amplitude
   and find a multiplicity-four cycle; Perko 1995 Thm 4.3 then gives an open
   region with **four simple cycles in one nest**. Triple cycles confirmed at
   normal amplitude. A negative result is already in: no nondegenerate
   swallow-tail at the small-amplitude end. **This is the live best hope.**
2. **D1 — the KKL fold component.** A reproducible `3+1` field on the finite-fold
   unfolding; no four-origin-cycle field; the component is not enclosed.
3. **`Q_4` route 4 outside-lobe four-interior-zero question.** Narrowed, open.
4. **`Q_3^R` at higher order.** The first-order kernel is nine-dimensional and
   is exactly the trivial span (coordinate changes, time rescaling, motion in
   the centre family), so there is no first-order-invisible reservoir. Higher
   order is still open — a kernel direction is conjugate to a centre only to
   first order — but it must be argued, not inferred from a kernel dimension.
5. **The DRR 121 graphics**, and every non-elementary graphic — nilpotent
   points, saddle-nodes. Untouched, and the hard end of the finiteness problem.
6. **Conjecture D2** (`sigma * eta_2 < 0`) is unproved; 46 numerical loops.

## 7. This session, with its standing corrected

I initially presented the order-two graphic collision as a new mechanism. **It
is not.** Fable's night report of 2026-09-05 already states the rigidity
principle in general and lists "order-two neutral infinity graphic" among its
instances; `CLAUDE_ROUTES_4AB.md` of 2026-09-04 already contains the antipodal
reciprocity and the order-three splitting scan, and already conjectured that no
closed graphic surrounds the origin except at a centre. My versions of those
are independent reproductions. The documents and the PR now say so.

What survives as new:

- the **exact** neutrality polynomial `N(l,m,a,b)` for the general Shi chart;
- the identity **`N = 640*eta_3`** on the order-3 stratum, which turns the
  existing numerical conjecture into algebra there (hand-checked at Shi's seed:
  2,850,000 both ways);
- the order-two instance **for the infinity graphic in the Shi chart** — the
  splitting vanishes only at `eta_2=0`, matching `a_deg = sqrt(-(l+2)/2)` to 14
  digits at six values of `l`, shown to be a transversal crossing;
- the `Q_3^R` first integral and the flow-verified reduction of `M_1` to a
  moment basis. **The rank claim attached to it was wrong** and is retracted:
  the generators satisfy an exact relation already in the repository, the rank
  is three, and the nine-dimensional kernel is exactly the geometrically
  trivial span, so there is no blind spot. See `Q3R_RANK_CORRECTION.md`.

An unintended but useful by-product: my Lyapunov chain, computed from scratch,
reproduces Astra/Fable's exact `eta_2` and `eta_3` with ratio exactly `-1` — a
single global convention difference. Two independent derivations of the focal
values now agree.

## 8. Instruments built

| Instrument | Where | Why it matters |
|---|---|---|
| Compactified return-map engine (C) | `audit/fable_engine/` | Counts cycles at radius `1e12` in milliseconds; reproduces every published four-cycle field. Solved the "we cannot even see the cycles" problem |
| **Exact Poincaré–Bendixson polygon certifier** | `fable/lane3-certify` | Rational predicates, no interval ODE integrator. First rigorous certificate machinery in the campaign. Cherkas row 1 certified: **3 limit cycles**, 82 s. Closes the long-standing "no certifier exists" gap |
| binary128 KKL fold solver | `fold_surface_*`, `codex/fastra-d1` | Fold continuation to radius `1.7e18` |
| Dual independent engines, Lanes 1 and 2 | `fable/lane1-ahcurve`, `lane2-cusp` | Different charts, variables and integrators; their agreement is the strongest validation either produced |

Certifier limits, documented rather than worked around: cycles with multiplier
within `1e-6` of 1, and the Songling cycles at radii down to `1e-202`, are out
of reach. It is a certifier for normal-size cycles.

## 9. Honest assessment

The campaign is well run. Evidence classes are labelled, negative results are
kept, withdrawn claims stay visible, and no document claims more than it has.

Three things are genuinely established: a real theorem (Theorem N), a real
instrument (the certifier), and a real structural principle (rigidity). The
principle explains every failure so far and points at the one mechanism that
escapes it.

Three things remain weak:

1. **Discovery still outruns verification.** The certifier arrived on day four
   and has certified one field. Nothing in four days of search has been
   certified end to end except that row.
2. **The budget is exhausted** (4096/4096) and the recent negative results come
   from a stopped campaign, not a finished one.
3. **Coverage is not coverage.** Sparse radial profiles, selected paths and
   sampled families recur throughout; several documents say so about
   themselves.

If I had to name the single most valuable next step: **prove Conjecture D2**.
It is the missing middle rung between Zhang (order two: at most one cycle) and
Li–Cherkas (order three: none), it would convert the rigidity principle from a
pattern into a theorem, and Astra's Proposition A already supplies the exact
ingredient for its neutral case.

## 10. Branch index

| Branch | Tip | What lives only there |
|---|---|---|
| `main` | 09-06 | The 32-commit trunk; not the frontier |
| `astra/fastra-afternoon-2026-09-05` | 09-06 | **Frontier tip.** The engine, all F-lane sweeps, the night report, the afternoon audit, D2 |
| `codex/fastra-d1-fold-counts-2026-09-05` | 09-06 | D1 fold sheets, the exact rational `3+1` KKL field, counter-failure audit |
| `claude/conjecture-progress-report-ixsmgv` | 09-05 | Subset of the afternoon branch |
| `astra/q4-determinant-2026-09-05` | 09-04 | `Q_4` sixth and seventh strikes |
| `astra/resonant-joint-2026-09-05` | 09-04 | The resonant joint incompatibility |
| `astra/fastra-d1-2026-09-05` | 09-05 | Earlier D1 snapshot, discrepancy ledger |
| `audit/post-q4-frontier-2026-09-04` | 09-04 | **Orphan, no merge base.** Canonical state, seed ledger, historical five-cycle claims, Galias–Tucker replay at 900–1200 bits |
| `fable/coordination-2026-09-06` | 09-06 | Gap reports, literature reports, protocol, lane briefs, coordinator log |
| `fable/lane1-ahcurve` | 09-06 | Andronov–Hopf curve engine, two independent sessions merged |
| `fable/lane2-cusp` | 09-06 | Cusp manifold engines A/B/C/D, swallow-tail Newton |
| `fable/lane3-certify` | 09-06 | The exact polygon certifier |
| `fable/compute-{evolve,f3-lam0,f5-shi,pert}` | 09-05 | Raw sweep output only |
| `opus/degeneracy-collision-2026-09-06` | 09-06 | This session: neutrality identity, `Q_3^R` structure, order-two collision, verification audit |
