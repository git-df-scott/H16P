# H16P: quadratic limit-cycle feasibility audit

Audit date: **2026-09-04**. Scope: real planar polynomial vector fields

\[
\dot x=P(x,y),\qquad \dot y=Q(x,y),\qquad \max(\deg P,\deg Q)\le 2.
\]

## Bottom line

**Astra verdict: YELLOW — attack only three specified narrow families.**

The exact global status is

\[
H(2)\ge 4,
\]

with **no known finite uniform upper bound**. In particular, even the statement
`H(2) < infinity` is open. The frequently repeated `H(2)=4` is a conjecture,
not a theorem.

A counterexample to `H(2)=4` must be one explicit real quadratic field with at
least five **distinct isolated periodic orbits**. A plot, five long-lived
trajectories, five roots of a truncated focus series, or five cycles in a
piecewise/quadratic-perturbation problem does not qualify.

The best certified frontier is the Galias--Tucker (2022) computer-assisted
interval proof that one exact Songling parameter instance has **exactly four**
limit cycles. This shows that candidate-level certification is viable. It does
not make global discovery bounded: after generic affine coordinate changes and
time scaling, quadratic fields still have five essential parameters, multiple
singular charts, unbounded geometric scales, and bifurcation strata requiring
arbitrarily fine resolution.

## Audit conclusions

| Question | Finding |
|---|---|
| Known lower bound | `H(2) >= 4`, from Shi and independently Chen--Wang (1979/1980) |
| Finite upper bound | None known; uniform finiteness of `H(2)` remains open |
| Best explicit proof | Galias--Tucker certify the Songling instance has exactly four |
| Four-cycle configuration | `3+1`: three nested about one focus, one about another |
| Five-cycle historical claim | Shi's preliminary five-cycle argument failed after correcting a sign error in Bautin's fifth focus quantity |
| Two-focus restriction | One nest must contain at most one cycle; a five-cycle two-nest target must be `4+1`, not `3+2` |
| Candidate certification | Strong: interval Poincare maps can prove five hyperbolic fixed points |
| Global searchability | Weak: no known finite-resolution cover or complete normalization atlas |

## Files

- [`STATUS.md`](STATUS.md) — binding status and evidence rules.
- [`LITERATURE_AUDIT.md`](LITERATURE_AUDIT.md) — primary-literature findings and claim audit.
- [`FOUR_CYCLE_FRONTIER.md`](FOUR_CYCLE_FRONTIER.md) — explicit four-cycle systems and a reproduced control.
- [`BIFURCATION_MECHANISMS.md`](BIFURCATION_MECHANISMS.md) — mechanisms that could or could not make a fifth cycle.
- [`RIGOROUS_CERTIFICATION.md`](RIGOROUS_CERTIFICATION.md) — hostile verifier and replayable certificate schema.
- [`SEARCH_SPACE.md`](SEARCH_SPACE.md) — normalizations, dimensions, and scaling reality.
- [`ATTACK_MATRIX.md`](ATTACK_MATRIX.md) — scores and three bounded attacks.
- [`ASTRA_HANDOFF.md`](ASTRA_HANDOFF.md) — execution handoff and stop rules.
- [`SOURCES.md`](SOURCES.md) — source ledger.
- [`controls/`](controls/) — one non-rigorous four-cycle numerical regression.

## Q4 preparation phase

The follow-on Q4 audit identifies the exact four-dimensional elliptic-integral
target, corrects the non-affine Hamiltonian-chart issue, reproduces three-zero,
zero-zero, and double-root controls, and supplies a bounded candidate screen
plus certificate schemas:

- [Q4_THEORY.md](Q4_THEORY.md)
- [Q4_PARAMETERIZATION.md](Q4_PARAMETERIZATION.md)
- [Q4_CONTROLS.md](Q4_CONTROLS.md)
- [ZERO_TO_CYCLE.md](ZERO_TO_CYCLE.md)
- [Q4_SEARCH.md](Q4_SEARCH.md)
- [Q4_CERTIFICATION.md](Q4_CERTIFICATION.md)
- [Q4_COST_MODEL.md](Q4_COST_MODEL.md)
- [q4/](q4/)

The Q4 verdict remains **YELLOW**. The only exploration was a 0.50-CPU-second
smoke test; no production sweep or five-cycle hunt was launched.
