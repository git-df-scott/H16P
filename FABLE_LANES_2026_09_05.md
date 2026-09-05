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
