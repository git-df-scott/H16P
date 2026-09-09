# H16P: quadratic limit-cycle feasibility audit

## Structural obstruction, 2026-09-06

[DEGENERACY_COLLISION.md](DEGENERACY_COLLISION.md) proposes a *mechanism* rather
than another negative result. For a nest bounded by an elementary graphic
through two non-antipodal infinite saddles in the Shi chart, the degeneracy
conditions of the focus and of the graphic are not independent: neutrality on
the order-3 stratum is exactly `eta_3 = 0` (the resultant is `640*eta_3`), and
on the order-2 stratum neutrality plus the graphic connection forces `eta_2 = 0`.
The common solution is a centre whose second finite singularity is a node, so
the remote cycle is lost too. The five independent conditions a fifth cycle
needs collapse to the four the known four-cycle systems already spend.

Everything is audited in [VERIFICATION_2026_09_06.md](VERIFICATION_2026_09_06.md),
including two errors found in checking scripts. The one load-bearing claim
still lacking an exact proof is the numerical identity `a* = a_deg`.
Non-elementary graphics, graphics through finite saddles and charts outside Shi
are explicitly not covered. **No counterexample; this is not a proof of H(2)=4.**

## Two coverage gaps closed and sized, 2026-09-06

[ORDER3_GRAPHIC_NEUTRALITY.md](ORDER3_GRAPHIC_NEUTRALITY.md) settles section 4b
of `CLAUDE_THOUGHT_SESSION.md`. On the third-order weak-focus stratum the first
stability coefficient of the boundary graphic through infinity equals one
exactly on the zero set of `eta_3`: the elimination resultant is `640*eta_3`.
The graphic-stability direction is therefore not independent of the focus
unfolding directions, and the proposed codimension-five `3+2` point does not
exist by this mechanism. Where the graphic's saddles are antipodal the
coefficient is identically one for every quadratic field, but a 268-point scan
finds the required connection nowhere on the stratum.

[Q3R_FIRST_ORDER.md](Q3R_FIRST_ORDER.md) gives the reversible two-centre
family an exact integrating factor `y^(a-1)`, a first integral, and a
four-dimensional first-order generating space, checked against direct flow
integration to `1.3e-12`. `M_1` vanishes identically on an **eight-dimensional**
subspace of the twelve quadratic perturbation coefficients, which no
first-order search can see. Whether `M_1` admits four zeros stays open: at 70
digits the Chebyshev determinants sit at the working precision.

**Neither result produces a five-cycle candidate. `H(2)=4` is neither proved
nor refuted here.**

## Fold closure follow-up, 2026-09-05

[KKL_FOLD_CLOSURE.md](KKL_FOLD_CLOSURE.md) proves a global single-negative-band
restriction for the multiplier polynomial on 1<=c<=8/5,K>0 and an impossibility
theorem for the proposed analytic monic-quartic scalar Dulac certificate at
a true fold. It also completes numerical return-sign reproduction at the
largest saved pair. No K1 candidate or component-wide exclusion was obtained.
The failed triple-root Newton attempt is archived as unresolved. Known shared
accounting is now 4096/4096 calls.

## Current fold-component checkpoint, 2026-09-05

Read [KKL_FOLD_SURFACE_STRIKE.md](KKL_FOLD_SURFACE_STRIKE.md) first. The finite
KKL fold was continued to horizontal radius approximately 2.96e17 on the
positive-K sheet, and through the center organizer onto a separate negative-K
sheet. Selected two-cycle pairs have complete-return numerical reproduction.
No three-origin-cycle field, 3+1 precursor, or five-cycle field was found.
The component and endpoint proofs remain incomplete. Work was stopped at the
user's instruction; all numerical failures and replay data are archived.
Inherited KKL/Shi ledger: 4053/4096 calls used, 43 unspent; the parallel
re-seed ledger below is separate. Older entries below
are chronological records.

## Reversible re-seed, 2026-09-05

Read [REVERSIBLE_RESEED_2026_09_05.md](REVERSIBLE_RESEED_2026_09_05.md).
The two-center geometry excludes the old finite saddle-loop mechanism.
An exact unfolding calculation finds and repairs a missing direction at
`a=-2`. A same-parameter moment search over 64 finite shape samples finds
no five-cycle candidate and supplies a four-cycle numerical control,
independently checked by high-precision quadrature and original-field
integration. **The full reversible route remains open; this is not a
proof of H(2)=4.** Replay code and logs are in `reversible_reseed/`.


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
