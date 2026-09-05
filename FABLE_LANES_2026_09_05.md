# Fable lanes, 2026-09-05: diversified counterexample attempts

Written before any computation. Astra runs the K1 plan (theory kill test on
"three cycles around a first-order weak focus", KKL cutoff diagnosis, Newton
fold continuation, incumbent-preserving KKL pilot with the coefficients -10
and 11/5 fixed, remote cycle, certification). Every lane below either asks a
different question or attacks a shared question with a method Astra does not
use. Nothing here touches the KKL family with -10 and 11/5 fixed.

## Governing facts

- Hyperbolic limit cycles are structurally stable. If H(2) >= 5 is true, it
  holds on an OPEN set of the five-dimensional parameter space. Random or
  Sobol sampling with a return-map counter is therefore a legitimate tool
  here, unlike the codimension-five reasoning for a single degenerate
  object in CLAUDE_THOUGHT_SESSION.md. The open set may be thin (nearly
  semistable cycles), so sampling must be adaptive around near-tangencies.
- Zegeling: exactly five is (5,0) or (4,1); (5,1) = six is permitted.
- Every four-in-one-nest configuration passes through: three cycles around
  a first-order weak focus (K1), or a cyclicity-two object near the
  order-2 or order-3 strata (closed numerically on the order-3 stratum by
  CLAUDE_ROUTES_4AB.md).

## Lanes

| Lane | Question | Method | Kill test |
|---|---|---|---|
| F1 | Literature 2023-2026: K1 theorems, (4,0), multiplicity-3 cycles, Yu-Han, Galias-Tucker, five/seven-cycle claims, Marin-Villadelprat | Web search via agent | none; informs the others |
| F2 | Engine: compiled adaptive return-map counter over a section from a focus, batched over parameters, with cycle stability from the crossing slope | C + OpenMP, validated on the KKL incumbent (roots 0.68321, 2.18370, 15.96278; remote near -3711.56) | fails validation |
| F3 | K1 by sampling: three cycles around the origin with trace zero in the Shi chart (four free parameters l,m,a,b) | Sobol sampling plus refinement where two cycles and a near-tangency coexist | no three-cycle nest after a dense scan |
| F4 | Fold descent from the KKL incumbent in the directions Astra holds fixed (the -10 and 11/5 coefficients) and in all five Shi directions | Track the four cycles, follow the displacement landscape toward a new near-tangency in the large nest | every direction loses a cycle before a new pair appears |
| F5 | Global five-dimensional sweep of the Shi chart for any nest with three or more cycles | Same engine; adaptive refinement | scan completes with maximum three per nest only at Hopf-type unfoldings |
| F6 | Route 4a: reversible center Q3R at the loop point; first-order Melnikov basis over the ovals; maximal zero count | Numerical oval integrals, Wronskian ECT check | count <= 3 or identically zero without higher-order relief |
| F7 | Q4 original-coordinate check: do actual cycles exceed the first-order count near the two-saddle infinity graphic? | Perturb the Q4 center in a lifted direction with three simple zeros of I and integrate the true field | actual count equals first-order count everywhere tested |
| F8 | Yu-Han reversible two-center (3,1) family: reproduce, then search its full neighbourhood for a fifth | Engine sweep around the published field | no fifth in a dense neighbourhood |
| F9 | Marin-Villadelprat resonant hemicycle family on a=-1 (only lower bounds known there) | Engine sweep of (b, eps0, eps1, eps2) | maximum count matches the published lower bound |
| F10 | Multiplicity-3 cycle around a trace-zero focus (D = D_r = D_rr = 0) | Newton on the engine from near-tangency seeds of F3 | no triple cycle found |

## Order and budget

F1 and F2 first (F2 is the enabler). Then F3, F5 as background sweeps on the
four cores while F4, F8, F9 are set up. F7 and F6 after the sweeps report.
F10 only from F3 seeds. Every candidate with three or more cycles in one
nest is re-integrated at tight tolerance from exact rational coefficients
before being reported, and any five-root field goes to Astra for hostile
reproduction before being called anything.

## Reporting

Results are appended below as they arrive, each marked NUMERICAL. Scripts live
in `audit/fable_*`.

## Results (NUMERICAL unless stated), appended 2026-09-05

### Engine (F2)
`audit/fable_engine/`: compiled Dormand-Prince return-map counter. Reproduces the KKL incumbent
(origin roots 0.68321, 2.18370, 15.96278; remote unstable cycle beyond radius 3665 on the ray away
from the origin), the Yu-Zeng field (origin 0.0244 U, 0.0603 S, 0.1043 U; large cycle around (1,0))
and the relaxed Shi field (0.0041 U, 0.0061 S, 0.0643 U; large cycle around (0,1)).
Hostile review (`REVIEW_engine.md`) found that the first production settings missed two of the three
Yu-Zeng cycles and that tolerances scaled with the focus offset. Fixed: the field is re-expanded about
each focus, tolerance is relative to the orbit radius, a noise floor of 5e-12 relative rejects
integration-noise sign changes, the return-domain edge is bisected, rays point away from the nearest
other equilibrium. After the fix all three seeds report their published counts at production settings.
Every sweep result recorded before the fix (local F3, cloud workers) is a lower bound only.

### F6: reversible center at the Shi loop point
First-order Melnikov span on the annulus is three-dimensional (singular values 18.1, 1.83, 0.43, then
1e-13). Maximum zero count over 200k random directions: two.

### F11: order-two weak focus with a homoclinic saddle loop (lam = 0, eta_1 = 0)
Scan over (a, b) with loops located by separatrix splitting. Two populations only:
loops with saddle trace exactly zero, all with eta_2 = O(1e-9) (centers), and focus-type loops with
nonzero trace. In every focus-type loop the trace and eta_2 have opposite signs, giving parity one
(exactly one existing cycle, as Zhang's theorem requires). Along the focus branch at a = -3 the ratio
trace/eta_2 is constant to four digits (-0.02758, -0.02757, -0.02756, -0.02755 at b = 1.22, 1.24,
1.26, 1.29) through the common zero at b = 1.275, where the branch crosses the center stratum
a^2(b+2l+1) = (b+1)(l+1)^2. Same sign structure at a = -2 and a = -1.5.
Conclusion: a neutral saddle loop around a second-order weak focus does not exist as a focus
configuration; it is a center. The codimension-four organizing center "order-two focus + neutral
loop" collapses onto the integrable stratum, exactly as the order-three loop did (Lane C) and as the
neutral infinity graphic did in KKL (Astra, K = J = 0).

### F12: KKL double center (beta = 0, K = 0, J(c) = 0)
Confirmed a center numerically (displacement 1e-13 across the origin annulus); the second antisaddle
also has trace zero. Melnikov analysis pending with the corrected normalisation.

### F13: Yu-Han reversible two-center family on the curve a4 = (a1-5)/3
At the published point (a1, a4) = (-30/7, -65/21): span dimension three, Taylor matrix at the center
singular (ratio 5e-13), so the element with v0 = v1 = v2 = 0 exists, as the Yu-Han construction
requires. Its interior zero count is zero: the reported crossing at r = 0.034 is the least-squares
residual floor, and a direct dense scan of the actual Yu-Zeng field shows D/r positive and increasing
from r = 0.11 to the annulus edge at 0.2333. No fourth origin cycle there. Grid over a1 running for
both annuli.

### F13 (grid, both annuli, a1 from -1.1 to -15 on the curve and off-curve controls)
Every point: span dimension 3, at most two zeros over random directions, Taylor matrix singular on the
curve (so the v0=v1=v2=0 element exists exactly there), and that element has no interior zero in either
annulus. The Yu-Han mechanism gives (3,1) and nothing more at first order plus Bautin. CLOSED.

### F14 / F4: descents from four-cycle seeds (fixed engine)
Random sampling of the widened KKL family (free quadratic coefficients) finds genuine (3,1) fields at a
rate of about four per 22,000 samples; eleven such seeds plus the KKL incumbent were descended for 60
generations. Maximum total four; the origin-nest displacement between the U and S roots is a single hump
with no forming fold. A reported near-miss was a noise wiggle next to a root (detector fixed).

### F15: neutral hemicycle (a = -1 in the two-center chart, J(c) = 0 in KKL) with the compactified engine
New engine `retmap_log.c`: log-polar returns, validated against Astra's KKL table (origin cycle near 7000,
remote near 6.5e9 at c = 0.9683, K = 1/64) and the incumbent at radius up to 1e17. Targeted unfolding
(270 fields): maximum total three, pattern (2,1). Broad Sobol sweeps (20k two-center, 10k KKL near c*)
running overnight.

### F16: Dulac coefficients at the neutral hemicycle (transverse section: height y above the invariant line)
Unperturbed neutral system returns as a center to 2e-14. Reversible directions da, db give zero first-order
displacement. Effective first order: D = pi e0 - pi (2 e1 + e2) y + O(y^2): connection term e0, linear
term the single combination 2e1 + e2; no w log w term at first order (it arises at second order as
da x (2e1+e2)). So at first order the neutral hemicycle emits at most one cycle; two with the ratio
parameter, matching Marin-Villadelprat's lower bound. Nothing at first order suggests three.
At b = 1 the field is holomorphic (z' = -(z^2 + 1/4)) and the upper-annulus span drops to two.

### F17: second order along the null direction at the holomorphic point
The null direction x(1-2y) leaves the upper center a center to all orders (displacement 2e-14 for
eps = 1e-3, 5e-4, 2.5e-4): it is center-preserving, not a second-order route. Lower annulus first order
is nonzero (D/eps = 1.078). CLOSED.

### Status at 05:00 UTC 2026-09-05
Closed tonight: F3/F5 (sweeps, lower bounds only, max one at trace zero), F6, F11, F12, F13, F17.
Open: F15 (overnight sweeps), F16 (second-order coefficient structure: does da x (2e1+e2) plus the
quadratic terms allow a third graphic-born cycle; the honest expectation is no).

### F19: Q4 alien test at the neutral resonant infinity graphic (compactified counting)
Melnikov span dimension four (singular values 20.3, 6.59, 1.51, 0.087). Random directions: at most two
zeros (47 of 300k); targeted construction gives fifteen three-zero directions and zero four-zero
directions in 20k trials. Along every three-zero direction the true field at eps = 1e-4 and 1e-5 has
exactly three cycles around the origin (S/U/S), including zeros within 2e-4 of the boundary realised as
cycles near radius 50. No count ever exceeds the first-order count: no alien. CLOSED numerically at (3,0).

### F18 / F18b: unfolding the hemicycle near the codimension-two point (e0 = 0, 2e1 + e2 = 0)
Neutral case a0 = -1 (b = 0.5, 1, 1.5; 3600 fields) and the proven alien points a0 = -1/2, b0 = 1/2 and 3/2
(2400 fields): maximum total three, always the (2,1) pattern; no nest with three cycles on these grids.
The grids are coarse in the scaling that matters (e0 against (2e1+e2)^2), so a scaled follow-up runs next.

### F18c: scaled unfolding (e0 ~ delta^2, da ~ delta, delta = 2e1 + e2) at a0 = -1 and -1/2, b = 0.5, 1, 1.5
4536 fields: totals 0/1/2/3 = 630/2016/1638/252, no nest with three cycles, no total above three.
The hemicycle emits its two known cycles and nothing more in any scaling tested. CLOSED numerically.

### F15 overnight sweeps
Two-center neutral unfolding (a = -1 + da, da log-uniform, eps log-uniform), 20032 fields, compactified
counting to radius e^40: totals 0/1/2/3 = 11574/6453/1988/17, never four. CLOSED numerically.
KKL near c* (c, K, beta, p, q free): in progress; total four appears at about 0.25 percent (the (3,1)
configurations), no five so far.

### F15 KKL c* sweep and F20 descent (compactified counting)
KKL near c* with (c, K, beta, p, q) free, 10048 fields: totals 0/1/2/3/4 = 3016/4391/1976/641/24, never
five; every four is (3,1) with the remote cycle at radius 1e6 to 1e9 next to the infinity graphic.
Descent from the 24 four-cycle seeds (40 generations, 64 per generation): maximum total four. CLOSED.

### F21: full Shi chart, compactified counting (radius to e^40), 30016 fields
Totals 0/1/2/3 = 21331/8418/266/1, never four. Large cycles near infinity graphics do not hide extra
nests in the generic chart either. CLOSED.

### Correction, 2026-09-05 morning: compactified counter noise floor
Astra's D1 report noted that the compactified counter missed known pairs near the center. Confirmed: on the
Yu-Zeng field the counter reported two origin cycles instead of three. The integrator was correct (profiles
identical at rtol 1e-12 and 1e-13); the rejection floor of 1e-10 in log-displacement was above the genuine
displacements (1e-11 to 5e-10) of near-integrable cycles. Fixed: floor 5e-12, adaptive refinement around
interior minima of the displacement (near-fold pairs), umax 36. Revalidated: Yu-Zeng 3+1, relaxed Shi 3+1,
KKL incumbent 3+1, Astra's table field 1+1. Consequence: the overnight compactified sweeps (F15, F18, F19,
F21) undercounted cycles with displacement below 1e-10; every record with three or more cycles is being
recounted with the fixed counter, and the Q4 alien test is being rerun. Their verdicts are provisional until
the recount finishes.
