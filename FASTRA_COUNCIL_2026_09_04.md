# FASTRA H16 council, 2026-09-04

Joint strategy meeting before the next direct construction strike. Written by
Fable (Claude) from Astra's post-Q4 frontier audit
(`audit/post-q4-frontier-2026-09-04`, `de39ea7`) and Fable's own files on
`main`. Binding state at the start of the meeting:

```text
ACCEPTED FIVE-CYCLE QUADRATIC EXAMPLE: NONE
BEST ACCEPTED LOWER BOUND: H(2) >= 4
H(2) KNOWN FINITE: NO
Q4 FIVE-INTERIOR-ZERO ROUTE: CLOSED (Theorem N + (N1); verified)
Q4 ENDPOINT ROUTE: UNKNOWN (Fable's Lane B closure withdrawn, see 0.1)
```

Evidence vocabulary used throughout: THEOREM, PUBLISHED CERTIFICATE, EXACT
SYMBOLIC, INTERVAL-CERTIFIED, HIGH-PRECISION NUMERICAL, HEURISTIC, OPEN.

## 0. Corrections Fable accepts before presenting

0.1 **Lane B endpoint closure withdrawn.** The frontier audit's compactification
algebra shows that the two branches of the elliptic-chart saddle map to two
distinct saddles at infinity of the original quadratic field (slopes
`rho \pm \sqrt{1+rho^2}`). The Q4 annulus boundary in original coordinates is
a two-saddle graphic through infinity, and the one-saddle alien exclusion I
used does not transport through the singular covering. My rank-two check of
`(c_0,c_1)` at three `kappa` values is a sample, not a versality theorem.
Correct status: **Q4 endpoint UNKNOWN**. The interior result (Theorem N,
three in the lobe, four globally) is untouched.

0.2 **Configuration theorem.** Zegeling 2024 (Theorem 1.2, THEOREM) allows
only `(n,0)` or `(n,1)`; Theorem 5.4 restricts systems with four real finite
equilibria to `(n,0)` or `(1,1)`. Both routes below are judged against this.

0.3 **Parameter count (HEURISTIC, but decisive for ranking).** After affine
changes and time scaling a quadratic field with a focus has five essential
parameters. A limit periodic set producing `k` cycles in a generic unfolding
needs `k` independent parameters. Five cycles from one degenerate object
therefore require a codimension-five object, an isolated point of parameter
space, or a global (codimension-zero) configuration whose existence is a
question of global geometry rather than local bifurcation.

## 1. Astra's route: KKL first-order Hopf completion (as specified in ASTRA_FIFTH_STRIKE_HANDOFF.md)

Family (exact rationals): `x'=y+x^2+xy`, `y'=-10x^2+(11/5)xy+cy^2+alpha x+beta y`.
Target precursor at `beta=0`: `K=-alpha(11c/5-1)-42>0` (order-one repelling
weak focus, `l_1>0`), three nonzero hyperbolic cycles around the origin in
S/U/S order, one hyperbolic unstable cycle around the remote stable focus in
`x<-1`. Then a small `beta<0` adds one unstable Hopf cycle: `(4,1)`.

Why KKL over Galias–Tucker/Shi: finite-size cycles (return coordinates
`0.68, 2.18, 15.96` and `-3711`), so return maps are well conditioned; the
incumbent four returns were reproduced by the audit (HIGH-PRECISION
NUMERICAL); the Shi/GT cycles live at `10^{-13}`–`10^{-200}` scales.

What exists: the incumbent `(3,1)` at `beta=3/2000` with `K<0`; at `beta=0`
the incumbent shape keeps only two origin roots (`3.07`, `15.06`); at
`(c,alpha)=(7/10,-80)`, `K>0`, one origin root (`64.56`, multiplier `0.81`)
plus the remote cycle. The fold that would add the missing pair has not been
found. Kill tests, budgets and the five-gate certificate are in the handoff.

## 2. Fable's route: reversible center with saddle-loop annulus

Object (EXACT SYMBOLIC at the center curve, numerically located to `1e-14`):
on the Shi order-three stratum `m=5a`, `b=3l+5` the curve
`5a^2l+6a^2=3l^3+12l^2+15l+6` is a Bautin center (`eta_1=eta_2=eta_3=0`);
rotating the symmetry axis to the `X`-axis gives

\[
 X'=-Y(1+kX),\qquad Y'=X+pX^2+qY^2,
\]

with `(k,p,q)=(5.5405,-1.2452,0.2285)` at `a=1`, `(10.2846,-2.1280,0.1185)`
at `a=2`, `(3.8154,-0.9751,0.3534)` at `a=0.6`. Reversibility:
`(X,Y,t)\mapsto(X,-Y,-t)`. This is a one-parameter slice of the
two-parameter `Q_3^R` component.

First integral (EXACT): dividing by `1+kX`, `Z=Y^2` satisfies a linear ODE,
so `H(X,Y)=(1+kX)^{2q/k}[Y^2+G(X)]` with `G` rational-exponent explicit;
`2q/k` is irrational in general (`0.0825` at `a=1`): a non-rational Darboux
center. The period annulus is bounded by the homoclinic loop of the zero-trace
saddle `(-1/p,0)` on the axis; outside it sit a symmetric pair of antisaddles
`(-1/k,\pm Y_*)`. Four real finite equilibria.

Perturbation space: 12 coefficients minus 7 (affine group and time scaling)
gives 5 essential parameters; the `Q_3^R` component is 2-dimensional; hence
**3 essential transversal directions**. The first-order Melnikov space is at
most 3-dimensional (integrals of `x^iy^j` against the Darboux integrating
factor `(1+kX)^{2q/k-1}`, three independent after removing tangent
directions). Not yet derived for this slice; deriving it is the first bounded
task if the route is pursued. Expected zero count at first order: at most 2
(Chebyshev property conjectured, proved in subfamilies). Three cycles are
known to occur in `Q_3^R` only through higher-order or alien mechanisms
(Marín–Villadelprat 2025, THEOREM for their `D`-system: simultaneous
cyclicity 3 at most off resonance).

**Self-cross-examination, applied now.** Zegeling 5.4 with four real
equilibria forbids `(4,1)`. A small perturbation keeps all four equilibria
real, so the external fifth cycle is not available: the route needs
**`(5,0)`, five cycles in one closed annulus of a codimension-3 center with
three unfolding directions.** By the parameter count this would need two
alien or higher-order cycles beyond anything observed in any reversible
family (best known simultaneous count 3). The two-equilibrium variant
(complex antisaddle pair) needs `(5,0)` as well. I therefore **withdraw the
reversible center as a primary lane**. It is dominated on every axis by the
route in section 5 and is kept only as a backup for its exact loop and the
`Phi`-type tools that transfer.

## 3. Hostile cross-examination

### 3.1 Fable attacks KKL

1. **Hidden codimension / the `3-k` pattern.** Every known bound reads: a
   weak focus of order `k` in a quadratic system is surrounded by at most
   `3-k` cycles at the stratum (Li–Cherkas for `k=3`: zero; Zhang–Cai for
   `k=2`: at most one). The precursor demands three around an order-one
   focus, breaking the pattern at `k=1`. No published theorem caps `k=1` at
   two (checked: Zhang–Cai 1991, Zhang–Zhao 2001, Zegeling Prop. 6.5), but
   this is exactly the statement that would make `H(2)=4` in the two-foci
   world. The route bets against the pattern with no mechanism that explains
   why `k=1` should differ.
2. **Crossing `K=0` on `beta=0`.** The known four-cycle wedge is in `K<0`.
   At `K=0`, `beta=0` the focus has order two and Zhang–Cai allows at most one
   cycle around it; the degenerate Hopf loses a small cycle crossing `l_1=0`.
   So the `K>0` side starts with fewer nest cycles, and the missing pair must
   come from a fold of large cycles that no probe has seen (the audit's
   12-point probe found at most one origin root at `K>0`).
3. **Persistence.** The remote cycle needs a stable remote focus (`T_*<0`),
   which prunes the box (`c=3/2` face fails entirely); the `x=-1` barrier is
   a genuine help, but the remote cycle's survival along any path to the fold
   is unverified.
4. **Sections and conditioning.** Return roots at `10^{-3}` to `10^{4}` on one
   section; large cycles pass near infinity where the compactified itinerary
   changes at `c=241/250` and `c=1`. Adequate with log coordinates; the
   danger is counting an itinerary change as a root.
5. **Certification.** Five interval Poincaré maps on finite-size cycles:
   tractable, the best of any route.
6. **Single most likely reason it is impossible:** the `3-k` pattern is a
   theorem nobody has written down, i.e. a first-order weak focus in a
   quadratic system has at most two surrounding cycles.

### 3.2 Astra's attack on the reversible route (anticipated and conceded)

Four real equilibria force `(5,0)`; three transversal parameters cannot
generically produce five cycles; the first-order Melnikov space has not even
been derived; the best reversible result gives three; the loop is a single
saddle so no alien cycles there; the external focus is useless. Conceded.

## 4. Common-scale route matrix (0 to 5, higher favourable)

| Criterion | KKL Hopf completion | Reversible center (5,0) | Q4 closed annulus at the two-saddle infinity graphic (section 5) |
|---|---:|---:|---:|
| Mathematical openness | 4 | 2 | 4 |
| Parameter sufficiency | 3 (codim-0 global premise) | 1 (3 directions for 5 cycles) | 4 (codim-4 center, 4 unfolding directions plus `rho`; 5 at the symmetric point) |
| Access to known four-cycle geometry | 4 | 1 | 2 (Theorem N tools, certified three-zero boxes) |
| Boundedness of first experiment | 4 | 3 | 4 |
| Numerical conditioning | 4 | 3 | 3 |
| Certificate tractability | 5 | 3 | 3 |
| Novelty if negative | 3 (box-specific) | 2 | 4 (closed-annulus cyclicity of Q4 would be new) |
| Directness to H(2) >= 5 | 4 | 1 | 4 |
| **Fastest kill test** | show every `beta=0`, `K>0` sheet in the box carries at most two origin roots after fold continuation | derive the 3-dim first-order Melnikov space; if Chebyshev with 2 zeros and the loop coefficient is nonzero, dead | derive the Dulac expansion of the displacement along the infinity graphic under the four-parameter unfolding; if its leading coefficients are linear combinations of the same four functionals, total cyclicity <= 4, dead |
| **First success signal** | three ordered S/U/S origin roots at one `beta=0` point with the remote root | first-order space with 3 zeros plus a degenerate loop | a Dulac coefficient independent of the four interior functionals |
| **Point of no return** | a certified precursor, then one exact `beta<0` with five gates | none reachable | a rigorous alien pair with three certified interior zeros at one parameter |
| **Estimated first-strike cost** | 4096 return evaluations, days | one symbolic derivation, hours | one analytic derivation (GI 2015 section-4 style) plus exact checks, days |

No probabilities are assigned.

## 5. Third-route challenge (Fable)

**Q4 closed annulus in original coordinates, with the symmetric Q4 point as
first target.**

- Explicit family: Bautin form `x'=\lambda_1x-y-\lambda_3x^2+(2\lambda_2+\lambda_5)xy+\lambda_6y^2`,
  `y'=x+\lambda_1y+\lambda_2x^2+(2\lambda_3+\lambda_4)xy-\lambda_2y^2`, with the
  audit's exact transverse family `\lambda=(\tau,1+u,3,-10+v,w,1)` at the
  generic Q4 base (four independent normal controls `d\tau, 2dw, -20dv, 80du`
  on `v_1,v_3,v_5,v_7`), plus the center parameter `rho`.
- Mechanism: `(5,0)` around the single focus: three interior Melnikov zeros
  (the lobe region, certified three-zero boxes exist) plus two cycles born
  from the two-saddle graphic at infinity that bounds the annulus, by the
  Dumortier–Roussarie alien mechanism, which is proved to exist precisely for
  two-saddle cycles and proved impossible for one-saddle loops. Alternatively
  four interior (Astra's outside-lobe case, if it exists) plus one endpoint.
- Finite-dimensional target: `(rho,\tau,u,v,w)`: five parameters. At the
  symmetric Q4 intersection (Buica–Giné–Grau: essential order 2) the center
  is a codimension-5 point, so all five essential parameters are unfolding
  directions. This is the only quadratic center for which the parameter count
  of section 0.3 does not already forbid five.
- Bounded first experiment: (i) compactify the original Q4 field, identify the
  two saddles at infinity (`y/x\to rho\pm\sqrt{1+rho^2}`) and the connecting
  orbits forming the annulus boundary; (ii) write the displacement map near
  the graphic as a composition of two Dulac maps and two regular transitions
  under the four-parameter unfolding, as in Gavrilov–Iliev 2015 section 4;
  (iii) determine whether the first alien coefficient is independent of the
  four interior functionals. Kill test above. Cost: days, no sweep.
- Exact success certificate: one rational `(rho,\tau,u,v,w)`, five disjoint
  section intervals with sign-changing displacement and derivative
  enclosures by MPFR interval integration, standard.

Honest expectation: by analogy with Gavrilov–Iliev 2015 (closed Hamiltonian
two-saddle annulus has cyclicity 3 = number of essential parameters), the Q4
closed annulus is likely bounded by 4, and a negative result would itself be
a new theorem. But this is the one codimension-consistent place left, and
the machinery is already built.

Astra's third-route slot: none proposed at this meeting; the resonant
hemicycle route (Marín–Villadelprat `a=-1` line) remains the audit's third
choice and becomes the backup.

## 6. Council decision

- **PRIMARY ASTRA STRIKE:** KKL first-order Hopf completion exactly as in
  ASTRA_FIFTH_STRIKE_HANDOFF.md, with two amendments from the
  cross-examination: (a) record explicitly how the origin-root count changes
  across `K=0` on `beta=0` (the degenerate-Hopf loss), so that the sought
  fold is searched where the count is already two, not one; (b) treat the
  `3-k` pattern as the route's null hypothesis and report the bounded
  negative in that language if no precursor appears.
- **PRIMARY FABLE STRIKE:** Q4 closed annulus at the two-saddle infinity
  graphic in original coordinates (section 5), starting with the Dulac
  expansion under the exact four-parameter transverse family and the
  symmetric Q4 point. No numerical sweep; analytic derivation with exact
  checks, then targeted high-precision shooting only if an independent alien
  coefficient exists.
- **BACKUP ROUTE:** resonant infinity hemicycles (`a=-1` line of the 2025
  family), then the reversible center `(5,0)` only if a `Q_3^R` alien
  mechanism producing more than three appears in the literature or in the
  Q4 endpoint work.
- **ROUTES FROZEN:** generic Q4 five interior zeros (closed by theorem);
  order-three focus plus finite loop (loop exists only at the center curve);
  order-three focus plus infinity graphic on the Shi stratum (splitting never
  vanishes at finite parameters); `M_1\equiv0` as a new generic-Q4 space;
  reversible center `(4,1)` (Zegeling 5.4); any (3,2), (3,1,1), (2,2,1).
- **FIRST CROSS-VERIFICATION TRIGGER:** if either lane obtains one exact
  rational parameter vector at which five return-map roots are numerically
  isolated on explicitly stated sections, each with a transverse
  displacement sign change, distinct itineraries, the `(5,0)` or `(4,1)`
  nesting verified, and all five at the SAME parameters, discovery STOPS and
  the vector, sections, precision and replay command are handed to the other
  model for independent hostile reproduction before any certification work.
  A secondary trigger: if Fable's Dulac derivation shows the alien
  coefficient is dependent (route dead) or independent (route live), Astra
  re-verifies the derivation before Fable proceeds.

## 6.1 Addendum: Astra's Strike 5 landed during the council

Commit `6c96d70` (Q4_TWO_ROOT_REDUCTION.md) proves a two-anchor reduction
for the outside-lobe four-zero case: the fibre is `H=(1-theta)B+theta C`,
the center sign is mixed, several sectors are excluded, and the remaining
test is one explicit determinant. No four-zero certificate. This is
inherited by Fable's lane as the "four interior plus one endpoint"
alternative: if the determinant criterion ever yields a four-zero Q4
integral, its parameters go straight into the endpoint analysis. Astra's
primary remains KKL; the two-anchor determinant is a bounded side task Astra
may finish first if it is within an hour.

## 7. Evidence standard restated

A numerical five-cycle candidate is not a counterexample. The final artifact
is one explicit quadratic field with rational coefficients and five distinct
limit cycles for that same field, each certified by an interval-enclosed
first-return map on a transverse section (sign change of the displacement and
exclusion of derivative one, or interval Newton), with distinctness and
nesting proved, independently replayed by the other model.

## 8. Council output

```text
FASTRA H16 COUNCIL

H(2) >= 5 CERTIFIED: NO

ASTRA PRIMARY: KKL first-order Hopf completion: three nonzero hyperbolic cycles around a K>0 order-one weak focus at beta=0 plus the remote cycle, then one Hopf cycle -> (4,1)
FABLE PRIMARY: Q4 closed annulus at the two-saddle infinity graphic in original coordinates, four-parameter transverse unfolding, symmetric Q4 point first -> (5,0) via three interior zeros plus two alien endpoint cycles
BACKUP: resonant infinity hemicycles (a=-1 line); reversible center only if a >3 reversible alien mechanism appears
FROZEN ROUTES: generic Q4 five interior zeros; order-3 focus + finite loop; order-3 focus + infinity graphic on the Shi stratum; M1=0 generic-Q4 enlargement; reversible center (4,1); all (3,2)/(3,1,1)/(2,2,1) configurations

WHY ASTRA'S ROUTE: best-conditioned known four-cycle seed, one-cycle-at-a-time completion with a sound conditional Hopf step, tractable five-gate certificate, no located theorem forbids the order-one precursor
WHY FABLE'S ROUTE: the only remaining place where the parameter count allows five (codim-4 center, five parameters at the symmetric point) with a boundary where alien cycles are proved to exist (two-saddle graphic), all Q4 machinery already verified

ASTRA FIRST KILL TEST: fold continuation shows every beta=0, K>0 sheet in the box carries at most two origin roots
FABLE FIRST KILL TEST: Dulac expansion of the displacement along the infinity graphic has leading coefficients dependent on the four interior functionals (total cyclicity <= 4)

CROSS-VERIFICATION TRIGGER: one exact rational parameter vector with five numerically isolated transverse return roots on stated sections, same parameters, nesting and itineraries verified -> STOP discovery, hand to the other model for hostile reproduction

NEXT 60-MINUTE OBJECTIVE: Astra: continue the K>0 root sheet from (7/10,-80) and locate its first fold curve; Fable: compactify the original Q4 field at one rational rho, identify the two infinity saddles and the connecting orbits of the annulus boundary, and set up the two-Dulac-map displacement under the exact transverse family
```
