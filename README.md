# H16P: quadratic limit-cycle feasibility audit

## Latest Q4 route-4 reduction, 2026-09-05

[Strike 7](Q4_SEVENTH_LIMITING_FACE.md) proves that the limiting original
family at `a=1` has at most three interior zeros. Both limiting boundary
determinants have exact negative-sign certificates. An explicit
three-dimensional subspace has at most two zeros at every finite lift.
The remaining finite-lift direction is unresolved: **route 4 is still open,
and no four-zero counterexample has been found.**

[Strike 6](Q4_SIXTH_BOUNDARY_REDUCTION.md) narrows the outside-lobe
four-interior-zero problem. Every survivor needs first primitive anchor
`r>1-(7/22)^(3/2)` and `kappa>2.89924108097...` (an exactly defined
algebraic cutoff). A sufficient global determinant exclusion now reduces
to two second-anchor boundary functions. Their universal signs remain
open. Exact proofs, rational interval certificates, and bounded numerical
controls are included. This does not solve `H(2)=4`.

## Current staged strike, 2026-09-05

Read [STAGED_RUN_2026_09_05.md](STAGED_RUN_2026_09_05.md) first. A new augmented
KKL solver located a numerical finite fold and independently reproduced a
nearby two-cycle sign pattern. Both earlier tracked cycles also persist
numerically well past the old radius cutoff. Neither result supplies M1 or
five cycles. The run used 400 KKL and 150 Shi/Chen–Wang returns, with all failures
recorded. Exact theory restrictions, rational coefficients and replay sources
are included. Older sections below are chronological records, not current
verdicts.


## Claude adversarial audit of Astra Strikes 1–3, 2026-09-04

Read [CLAUDE_AUDIT_ASTRA_1_3.md](CLAUDE_AUDIT_ASTRA_1_3.md) and
[FASTRA_H16_HANDOFF.md](FASTRA_H16_HANDOFF.md). All three strikes survive
independent hostile verification with no correction to the canon. A new
tuning-independent necessary condition for five zeros, `Phi(tau_1)>0`, fails
by 60–99.99% at every point examined and tends to equality only in the
double limit where (S1) fails. The proposed `Y=Y'=0` search is not
justified; Strike #4 is re-targeted at an exclusion theorem (`Z_distinct<=4`
for Q4). Independent checkers: [`audit/`](audit/).

## Astra second strike, 2026-09-04

Read [ASTRA_SECOND_STRIKE.md](ASTRA_SECOND_STRIKE.md) for the latest result.
The weighted-lobe region is a bounded analytic cell with an explicit rational
interior box. Exact Green reconstruction excludes that box for every kappa.
Every five-zero target must have its first primitive root after `5/11` and
`kappa>21636/19043`. No five-zero candidate was found; Q4 remains live.

- [Q4_LOBE_REGION.md](Q4_LOBE_REGION.md)
- [Q4_RECONSTRUCTION_GEOMETRY.md](Q4_RECONSTRUCTION_GEOMETRY.md)
- [Second-strike verification](q4/data/second_verification.txt)

## Astra reasoning strike, 2026-09-04

Read [ASTRA_FIRST_STRIKE.md](ASTRA_FIRST_STRIKE.md) for the inherited result.
The new structural reduction makes the auxiliary geometry independent
of kappa and excludes open neighborhoods of its interior cusp from
producing five original zeros. Q4 remains live; no five-zero candidate
has been found. The previous finite kappa cutoff and linear P2 filter
are corrected. Production-search commands are superseded.

- [Q4_STRUCTURE.md](Q4_STRUCTURE.md)
- [Q4_ZERO_GEOMETRY.md](Q4_ZERO_GEOMETRY.md)
- [Q4_CERTIFICATE_PLAN.md](Q4_CERTIFICATE_PLAN.md)

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
