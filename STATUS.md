# Status

## Latest staged work: 2026-09-05

See [STAGED_RUN_2026_09_05.md](STAGED_RUN_2026_09_05.md). There is a numerical
finite KKL fold/two-cycle-pair signal and no certified M1 or five-cycle field.
The pair field's remote equilibrium is exactly an unstable focus. The earlier
KKL pair persists numerically to c=0.9683 beyond the old cutoff, at different
parameters. The two classical-seed trace paths lose their innermost cycle into
the focus. These findings do not prove a family-wide exclusion or K1 false.

Budget: 206 historical + 400 new KKL + 150 new Shi/Chen–Wang = 756 returns;
3340 remain against the conservative common 4096 ceiling. Stage 2's 400-call
allocation is fully spent. Global graphic connection/coefficient work, general
Shi-chart fold continuation and an interval return-map verifier remain open.

Audit cutoff: **2026-09-04**

## Binding mathematical status

1. For a fixed real planar polynomial vector field, the number of limit cycles
   is finite (the Dulac theorem, completed by Écalle and Ilyashenko; Bamon gave
   a quadratic-case proof).
2. This pointwise finiteness does **not** imply a degree-uniform bound.
3. For quadratic fields, the best proved global lower bound is `H(2) >= 4`.
4. No finite uniform upper bound for `H(2)` is known.
5. `H(2)=4` remains a conjecture as of the 2026 sources audited here.
6. No accepted five-limit-cycle quadratic field was found.

It is safest to define

\[
H(2)=\sup\{\#\text{ limit cycles of }X:X\text{ is a real quadratic field}\}
\in \mathbb N\cup\{\infty\}.
\]

With this convention, the exact knowledge is `4 <= H(2) <= infinity`; the
right inequality expresses the absence of a finite uniform bound, not a claim
that any one field can have infinitely many cycles.

## Evidence classes

| Label | Meaning | Counts toward a five-cycle counterexample? |
|---|---|---|
| `NUM` | Floating-point integration, continuation, or plotted trajectories | No |
| `ASYM` | Formal focus/Melnikov expansion with uncontrolled remainder | No |
| `THEOREM-FAMILY` | Exact existence theorem for a stated parameter region | Yes, if it supplies one explicit admissible coefficient point |
| `CAP-EXIST` | Replayable interval/computer-assisted proof of distinct isolated cycles | Yes |
| `CAP-EXACT` | CAP proving existence and absence of any additional cycles | Stronger than needed |

## Counterexample acceptance test

A valid result must provide:

- exact rational coefficients or outward-rounded coefficient intervals;
- five pairwise disjoint isolating regions or transversal-section intervals;
- for each region, a proof of a Poincare-map fixed point and isolation
  (preferably `1 notin P'(I)` or an interval-Newton inclusion);
- a proof that the return map is defined with the stated orientation;
- replay code, dependency versions, precision, and raw interval logs;
- independent replay on a clean machine.

An exact count of all cycles is unnecessary for refuting `H(2)=4`; proving at
least five distinct cycles is sufficient.

## Q4 route closed (2026-09-04)

Theorem N ([Q4_THEOREM_N.md](Q4_THEOREM_N.md), verified in
[CLAUDE_AUDIT_ASTRA_4.md](CLAUDE_AUDIT_ASTRA_4.md)) with the necessary
condition (N1) proves that no Q4 Abelian integral has five distinct zeros in
the open annulus, for every `kappa>1`. The proved bounds are three distinct
zeros on the strict lobe region and four globally. Lane B
([CLAUDE_LANES_B_C.md](CLAUDE_LANES_B_C.md)) shows the saddle loop cannot
supply the difference. Attack 1 of ATTACK_MATRIX.md is therefore closed as a
counterexample route; its remaining content is the outside-lobe four-zero
question, a pure-mathematics item.

## Campaign decision

### Q4 route 4 update, 2026-09-05

**NARROWED, STILL OPEN.** [Strike 6](Q4_SIXTH_BOUNDARY_REDUCTION.md)
proves that any four-interior-zero candidate outside the lobe must have
`r>1-(7/22)^(3/2)` and `kappa>kappa_*`, where
`kappa_*=2.89924108097...` is defined exactly in the proof. Equality is
excluded. The second-anchor determinant has at most one interior zero
unless it vanishes identically, and two nonpositive boundary signs exclude
the entire second-anchor interval. Those two boundary signs are still
unproved. The global Q4 distinct-interior-zero bound remains **four**;
the strict-lobe bound remains **three**. No new limit-cycle example or
global quadratic upper bound is claimed.

**YELLOW — no generic coefficient sweep.** Only the three attacks in
[`ATTACK_MATRIX.md`](ATTACK_MATRIX.md) are authorized by this audit. Any
expansion requires a new mathematical reduction, a new bifurcation identity,
or a verified five-cycle candidate.
