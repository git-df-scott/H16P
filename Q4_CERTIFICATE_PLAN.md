# Q4 certificate plan after the first reasoning strike

Date: 2026-09-04. No five-zero parameter or certificate exists in this
package. This is a precise proof plan for a future candidate.

## Minimal scalar certificate

Input one exact kappa>1 and one exact nonzero vector mu in the original
four-term area-integral convention of Q4_THEORY.md. Prefer rational rho,
kappa=1+rho^2, and rational mu; rational rho also gives rational original
center coefficients. Rational kappa alone need not give rational c.

Choose six ordered rational points
\[
1<s_0<s_1<\cdots<s_5<\kappa.
\]
Enclose all six values of the SAME I in outward-rounded intervals excluding
zero, with alternating signs. The intermediate value theorem then supplies
five roots in the pairwise disjoint open rational intervals
\(J_i=(s_{i-1},s_i)\).

**Simplicity theorem.** Zhao's upper bound counts multiplicity and is five
for this exact nonzero analytic family. Thus those five roots are simple
and are the only interior roots. No derivative estimate or search of the
complement is needed for this logical conclusion.

If a format requires five pairwise disjoint CLOSED rational intervals,
shrink the five open gaps and certify endpoint signs on each resulting
closed interval; record ten endpoint evaluations. Continuity guarantees
such choices exist, but the actual endpoints and enclosures must be
included in a replayable certificate. Do not call the adjacent closed
gaps disjoint, since they share endpoints.

## Candidate admissibility and evidence

Verify exactly:

1. kappa>1 and every endpoint strictly inside (1,kappa);
2. mu is nonzero, and all four coefficients refer to the same basis;
3. coefficient transport includes the two relabelings after Zhao (20);
4. corrected P1 and P2 conditions hold when beta1 is normalized to one;
5. all signs use one fixed parameter point, or uniform enclosures over one
   common parameter box containing an explicitly selected exact point.

P1/P2 are necessary filters, not the definition of a Q4 center. Their
violation by an alleged five-zero candidate signals an inconsistency that
must be resolved. Passing them alone supplies no zero.

Store rational inputs as exact numerator/denominator pairs or decimal
strings parsed exactly, not binary floats. Record code hashes, dependency
versions, precision, quadrature/propagation bounds and raw interval output.
An independent replay must reproduce the inequalities. The existing
zero_certificate.schema.json describes the stronger interval-Newton route;
it is not itself an implementation or an existing certificate.

## Validated evaluation routes

Use either of the already audited approaches:

- validated PF propagation from a center Taylor expansion with a proved
  remainder, transporting the full original basis; or
- rigorous area quadrature on the oval after resolving the cubic roots
  and endpoint coordinate singularities with interval error bounds.

The universal hypergeometric/Stieltjes formulas provide exact auxiliary
structure and initial data. They do not replace the original third-kind
integral contributions. The existing floating evaluators, hypergeometric
diagnostic checks and exact symbolic checks are not interval quadrature.

Six signs may avoid interval derivatives at the scalar stage. Derivative
enclosures and interval Newton remain useful for tighter locations,
conditioning and later quantitative perturbation bounds.

## Discovery gate before expensive certification

A target with five simple I zeros necessarily induces three simple zeros
of the universal weighted primitive. Test the three strict inequalities
(L) in Q4_ZERO_GEOMETRY.md before lifting a shape through kappa-dependent
reconstruction. Entire neighborhoods of the interior auxiliary cusp are
now excluded. No blanket exclusion of endpoint degenerations was proved.

Any numerical five-zero lead stops exploration immediately. Freeze exact
or rationally approximated coefficients and switch to independent replay
and rigorous signs. No such lead was found in this strike.

## From scalar zeros to an explicit quadratic field

The scalar certificate and explicit field certificate remain different:

1. Reconstruct an analytic perturbation arc in the ORIGINAL Q4 coordinates
   with polynomial degree at most two. An arbitrary quadratic perturbation
   in the cubic Hamiltonian chart is not a substitute.
2. Verify lower Melnikov functions vanish and the first nonzero one equals
   the certified I times a nonzero factor.
3. Use its five simple roots on a compact subannulus. The implicit function
   theorem gives five cycles for sufficiently small nonzero epsilon.
4. For an explicit field, provide the twelve original coefficients and one
   explicit rational nonzero epsilon. Quantify the perturbation remainder
   or validate five return maps directly.
5. Verify first-return orientation, transverse sections, five distinct
   simple fixed points and pairwise disjoint flow tubes/annuli.
6. Perform an independent replay and present the human-readable proof.

No perturbation arc, epsilon, return-map certificate or five-cycle claim
is supplied here. The conditional passage is detailed in ZERO_TO_CYCLE.md.

## Evidence classes for this package

- Analytic proof: universal PF/Stieltjes representation, auxiliary Wronskian
  signs, cusp exclusion, weighted-lobe criterion, and sign-certificate lemma.
- Exact computation: rational symbolic identities and coefficient determinant.
- Numerical diagnostic: three specified hypergeometric identities and
  independent evaluator regression. These do not certify zeros.
- Conjectural/unresolved: original I Chebyshev behavior and whether five
  simple zeros can occur.
