# Q4 certification architecture

There are two separate proofs: five simple zeros of the Abelian integral, then
five fixed points of an explicit quadratic field's return map. Neither can be
replaced by plotting.

## Layer 1: interval Abelian zeros

Input exact rational \(\kappa,\mu_1,\ldots,\mu_4\), or narrow rational
outward enclosures. Give five disjoint rational intervals
\(S_i=[a_i,b_i]\Subset(1,\kappa)\).

Two viable validated evaluators are:

1. express the elliptic periods in Arb/Acb complete elliptic integrals and
   transport the basis with the Picard--Fuchs system using a validated Taylor
   ODE solver; or
2. use the vertical cubic roots and the substitution
   \(x=m+d\sin\theta\), split \([-\pi/2,\pi/2]\), and apply
   interval Gauss--Jacobi or tanh--sinh quadrature with analytic tail bounds.

The second is closest to the control code. On each \(S_i\), enclose both
\(I(S_i)\) and \(I'(S_i)\). Prove \(0\notin I'(S_i)\), opposite endpoint
signs, and the interval-Newton inclusion

\[
 N(S_i)=m_i-\frac{I(m_i)}{I'(S_i)}
 \subset\operatorname{int}S_i.
\]

This proves one and only one simple root in each box. Pairwise disjoint boxes
prove five distinct roots. For “at least five,” nothing needs to be proved
about the remaining interval. A claim of exactly five additionally needs
monotonicity subdivision of the complement and endpoint treatment.

The homoclinic endpoint has logarithmic terms; the center endpoint has a
regular Taylor series. Use Zhao's explicit expansions with interval remainder
bounds, never direct cancellation-prone quadrature arbitrarily close to an
endpoint.

The replay format is q4/zero_certificate.schema.json. It records exact inputs,
root boxes, signs, derivative bounds, Newton images, library/version/precision,
and hashes of all generated tables.

## Layer 2: realization and explicit epsilon

Recover an analytic arc of degree-two polynomials in the original Q4
coordinates. The certificate must contain:

- exact \(b,c\) and the identity \(b^2+c^2=4,c\ne0\);
- all twelve coefficients of \(X_2,Y_2\) as functions of \(\varepsilon\);
- symbolic equalities making lower Melnikov functions vanish;
- the exact linear identity between the first nonzero \(M_k\) and \(I_\mu\);
- one rational nonzero \(\varepsilon\).

Next continue the five predicted orbits at high precision in the original
quadratic field. Define transverse Poincaré sections, orientations, first
return conventions, seed intervals, and return-time windows.

## Layer 3: interval Poincaré proof

For each section interval \(X_i\), CAPD or an equivalent validated ODE library
must compute the return map \(P_i\) and derivative. Certify

\[
 N_i=m_i-\frac{P_i(m_i)-m_i}{P_i'(X_i)-1}
 \subset\operatorname{int}X_i,\qquad
 1\notin P_i'(X_i).
\]

This proves a unique hyperbolic fixed point in each interval. Validated flow
tubes or nested isolating annuli must be pairwise disjoint and must verify the
first-return event, excluding iterates and section-start artifacts. Enclose
equilibria to show no tube crosses one unexpectedly.

The interface is q4/poincare_candidate.schema.json. It deliberately leaves
interval_fixed_points empty until an actual candidate exists. Galias--Tucker's
Songling proof is the architectural control: interval Poincaré maps plus
topological separation, not a claim that its exact code applies unchanged.

## Replay levels

- Q4-ZERO: five interval Newton root inclusions.
- Q4-REALIZED: exact original quadratic perturbation and Melnikov identity.
- CAP-CANDIDATE: five validated return-map fixed points and disjoint tubes.
- COUNTEREXAMPLE: independent clean replay plus a human-readable proof.

Only the last level is evidence that \(H(2)\ge5\).
