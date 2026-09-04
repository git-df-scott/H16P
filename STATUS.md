# Status

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

## Campaign decision

**YELLOW — no generic coefficient sweep.** Only the three attacks in
[`ATTACK_MATRIX.md`](ATTACK_MATRIX.md) are authorized by this audit. Any
expansion requires a new mathematical reduction, a new bifurcation identity,
or a verified five-cycle candidate.
