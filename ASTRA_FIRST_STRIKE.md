# Astra first strike: Q4 structure and an excluded mechanism

Date: 2026-09-04. Base repository commit: 9db4cb3.
Stop condition reached: **C — a structural theorem substantially reduces
the attack.** No five-zero construction is claimed.

## Internal campaign map, corrected from authoritative formulas

| Category | Exact content |
|---|---|
| EXACT TARGET | Original four-term Q4 I_mu; kappa>1; projective nonzero mu; five distinct simple roots in 1<s<kappa. |
| PROVED CONSTRAINTS | Multiplicity bound five; exact coefficient transport; corrected strip (54-23kappa)/31<beta0<1; safe cubic P2 filter. |
| NUMERICAL EVIDENCE | Inherited three-zero kappa=4 control and tiny smoke records. No four/five lead. Old pruning decisions are superseded. |
| CLOSED SUBFAMILIES | beta1=0 and corrected strip complement; safe-concavity region; newly, open coefficient neighborhoods of each fixed interior universal auxiliary cusp. |
| LIVE PARAMETERS | All kappa>1 remain analytically possible. Universal auxiliary (A,B,eta) geometry can be studied before kappa-dependent reconstruction. |
| CERTIFICATE REQUIREMENTS | Exact admissible parameters plus six rigorous alternating signs and the multiplicity bound; original-field realization and explicit-epsilon validation remain separate. |

This map corrects the inherited material instead of treating its search
bounds as authoritative mathematical facts. The old integral family and
normalization are unchanged.

## Main result

Set d=kappa-1, t=(kappa-s)/d and
\[
F(t)={}_2F_1(1/6,5/6;1;t),\quad
M(t)=1-6(1-t)\frac{F'(t)}{F(t)}.
\]
The exact normalized auxiliary quotient is
\[
q(t)=\frac{g(kappa-dt)}{dt}
=A+Bt-1+(t-\eta)M(t).
\]

The first structural theorem is a positive Stieltjes representation of M.
It proves that the space {1,t,M,tM} has strict Wronskian signs +,+,+,-.
Thus its geometry is globally extended Chebyshev and independent of kappa.
Its triple-root curve is explicit, with eta decreasing from 54/31 to 1.
Proofs and exact coefficient transport appear in Q4_STRUCTURE.md.

The stronger exclusion concerns the ORIGINAL target. For every fixed
interior point on that auxiliary triple-root curve, an open neighborhood
of its normalized coefficients produces **at most three original I zeros,
for every kappa>1**. The reason is that
\[
\mathcal F(kappa-dt)=-d^2J_1(kappa)\int_0^t uF(u)q(u)\,du.
\]
At the cubic contact the primitive has already accumulated a strictly
positive area. Nearby three auxiliary crossings cannot undo it; the
primitive has at most one zero, and the remaining operator adds at most
two. The neighborhood may shrink toward either endpoint. This is not a
global exclusion of Q4. See Q4_ZERO_GEOMETRY.md for the full proof.

Any actual five-zero target must instead satisfy three strict weighted-lobe
inequalities for this universal primitive. These compare signed magnitudes,
not just crossing counts, and remove kappa from the first decisive screen.

These are new deductions in this campaign. No claim of priority over every
unexamined publication is made; the requested literature audit was not redone.

## Inconsistencies resolved

1. The inherited beta strip had its lower numerator reversed. The correct
   (54-23kappa)/31 follows both from a direct endpoint derivative and from
   the universal curve. The supposed kappa<85/23 exclusion is withdrawn.
2. The inherited P2 filter used a linear numerator. Polynomial division
   supplies a factor 2 and a safe CUBIC bound. Operational filters and
   exact rational boundary tests are corrected.
3. The prior 24-CPU-hour command is superseded by the user's bounded
   reasoning instruction. No production search ran.
4. Four arbitrary prescribed zeros generically annihilate all coefficient
   vectors. Three anchors determine a direction; subsequent zeros require
   vanishing determinants.
5. Four simple plus one double, or three simple plus one triple, would
   have multiplicity six and cannot occur. Three simple plus one ordinary
   double is a sufficient conditional mechanism, not an exhibited one.
6. Six rigorously alternating signs already force all five roots to be
   simple, using the multiplicity upper bound.
7. The center Taylor series must retain the area-period normalization and
   the factorials when using derivative coefficients from the source.

Historical negative controls and numerical records are preserved. No
historical computation has been relabeled as a rigorous exclusion.

## Work performed and verification scope

- Inherited the exact Q4 family and read the audited handoff, parameter,
  theory, search, controls, certification and realization material.
- Checked only the cited primary Q4 formulas needed to resolve inconsistencies.
- Derived the universal PF system, singularities/exponents, Riccati equation,
  exact endpoint expansions for its period pair, and positive measure.
- Proved auxiliary Wronskian signs, projective quadric geometry, fold/cusp
  formulas, original five-contact rank conditions and the cusp exclusion.
- Derived a necessary strict weighted-lobe criterion for five distinct I zeros.
- Corrected two false-pruning rules; added exact rational regression checks.
- Prepared a simplified scalar certificate and preserved the realization gate.

The bounded replay checks the coefficient-map determinant, coordinate
identity, rational period recurrence and three numerical diagnostic points.
The operational regressions check exact filter boundaries and independent
area/orbit evaluation. Both are small single-thread runs. Recorded results
and environment are in q4/data/astra_verification.txt.

The following were NOT proved or completed, and are not needed to invoke
stop C: original-basis Wronskian signs; a common global scalar disconjugate
ODE for all four original functions; full original third-kind endpoint
expansions with certified remainders; a five-zero candidate; a realized
perturbation arc; validated return maps. No global sweep, optimization or
large symbolic elimination was attempted.

## Requested result report

FIVE Q4 ZEROS CERTIFIED: NO
FIVE-ZERO NUMERICAL CANDIDATE: NO
NEW STRUCTURAL THEOREM: YES
FIVE-ZERO BIFURCATION MECHANISM: NO
Q4 STILL LIVE: YES

The explicit auxiliary cusp is a three-zero mechanism, and its neighborhood
is now excluded as a five-zero original mechanism. Those are different claims.

## Single strongest next Astra task

Characterize the universal coefficient region satisfying all three strict
weighted-lobe inequalities (L), outside the excluded cusp neighborhoods,
and determine whether its kappa-dependent PF/Green reconstruction can
saturate the remaining two-zero allowance. The decisive target is the
sign structure of that reconstruction, not a larger coefficient scan.

## Deliverables

- Q4_STRUCTURE.md — exact PF/Stieltjes reduction and proof.
- Q4_ZERO_GEOMETRY.md — projective geometry, Wronskians, cusp exclusion and lobe test.
- Q4_CERTIFICATE_PLAN.md — precise scalar and original-field proof gates.
- q4/q4_structure_checks.py — bounded exact/diagnostic replay.
- q4/notes_audit.md — independent filter and certificate audit.
- Corrected handoff, parameterization, filter code and regression tests.

This bounded reasoning strike is complete under stop C. The ultimate
five-cycle construction remains unresolved.
