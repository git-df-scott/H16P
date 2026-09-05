# Independent angular and two-half shooting review

This is a source and mathematical review, with **zero new ODE evaluations**.
It covers `angular_ld.cpp`, `continue_angular.py`, `half_ld.cpp`, and
`continue_half.py` as read on 2026-09-05. It is not interval certification,
an exhaustive return-domain check, or a completed fold-component analysis.

## Differential equations and derivatives

Put x=exp(w)cos(s), y=-exp(w)sin(s), alpha=-m. The code's H is
the physical logarithmic radial velocity and G is the physical angular
velocity. Consequently w_s=-H/G and t_s=-1/G. The divergence integrand
is -(exp(w)/G)[(21/5)cos(s)+(1+2c)(-sin(s))]. The parameter derivatives
alpha_c=55(K+42)/(11c-5)^2 and alpha_K=-5/(11c-5) are correct.

Writing f=-H/G, N=h2*g1-h1*g2, the implemented identities

    f_w=-exp(w)*N/G^2
    f_ww=f_w*(1-2*exp(w)*g2/G)

and their c,K derivatives are correct. The variational states are w,
M=log(∂w/∂z), w_c, w_K, M_z, M_c, M_K, elapsed physical
time, and accumulated physical divergence. Here z is the initial log
radius. In particular M_z integrates f_ww*exp(M), while
M_c integrates f_ww*w_c+f_wc. The final L derivatives follow correctly:
L_z=exp(M)-1, L_zz=exp(M)M_z, L_zc=exp(M)M_c.

The independent horizontal-section flux expression is also correct:

    w_final,z = (r/R) [Q(r,0)/Q(R,0)] exp(integral div dt).

At a closed orbit the physical multiplier equals exp(integral div dt).
Away from a closed orbit that exponential alone is not the section-map
derivative. The original engine appropriately exposes its flux comparison.

## Two-half map and first return

Let A(z) be the negative-ray log radius reached forward from positive-ray
log radius z, and B(z) the negative-ray log radius reached backward from
the same positive ray. Both half-map derivatives are positive. Replacing
s by direction*s and multiplying **all nine** right-hand sides by direction
is correct, including the backward physical-time and divergence integrals.

The residuals are

    F=A-B, G=log(A')-log(B'),
    F_z=exp(M_b)*expm1(G).

Whenever the half passages have a common smooth return domain, the full
return satisfies P=B^{-1} composed with A. Thus F=0 is equivalent to a
closed orbit, sign(F)=sign(P(z)-z), and at a match P'=exp(G). The physical
period is the forward time minus the backward time. The independent
divergence multiplier exp(I_f-I_b) agrees with exp(G) **at a match**.
The code's names ending in `_at_match` express the needed restriction.

Solving F=G=0 is an equivalent fold system. At a fold,

    F_zz=exp(M_b)*G_z.

Therefore G_z nonzero supplies radial nondegeneracy, and a nonzero
parameter derivative of F supplies an unfolding direction. The fixed-z
continuation tangent J_(c,K)^(-1)*[-F_z,-G_z] is correct. The corresponding
full-map formula in `continue_angular.py` is also correct.

Off-root, G=0 detects stationary points of **F** exactly; it does not in
general detect stationary points of P(z)-z. Counting stationary points of
F would suffice for root-count arguments on a proved common domain, but
sparse stationary sign brackets are not such a count.

## Implementation findings communicated to the continuation owner

1. The profile's choice c_pair=c-sign(F_c)*delta assumes a minimum fold,
   G_z>0. The general two-root side is
   c_pair=c-sign(F_c*G_z)*delta, with nonzero F_c,G_z guards. The full-map
   version analogously uses sign(L_c*L_zz). An explicit minimum-fold gate
   is an equally valid solution. Without either, a change of curvature
   can silently direct the profile to the zero-root side.
2. `continue_half.py`'s `refine` breaks on an UNRESOLVED result and then
   accesses its missing G field. It should retain the last valid
   approximation and preserve the unresolved status. Otherwise a single
   failed refinement can abort logging the current continuation event.
3. Small fold residuals alone do not establish nondegeneracy or connected
   component coverage. The code records the (c,K) Jacobian determinant;
   a continuation report must additionally retain G_z (or L_zz), any
   conditioning loss, and failed chart/corrector steps. A failed fixed-z
   chart can require another continuation coordinate.

These are source findings at review time, not a claim that owner fixes
remain pending after subsequent edits. The mathematical half-map repair
itself is sound and directly reduces the large full-return sensitivity.

## Scope of numerical acceptance

The angular engine requires physical angular velocity G<0 and checks
that condition at numerical stages. It therefore covers trajectories
with monotone angle. It neither proves the sign between stages nor
excludes cycles outside that chart. With actual G<0 throughout the
passage, a full turn is a first return to the positive ray; two matching
half turns have the same property. A failed angular chart is unresolved,
not evidence of no cycle.

Modified-midpoint extrapolation and the in-place Richardson recursion
have the expected form. The extrapolation difference is an error
estimate, not a rigorous enclosure. Decimal input preservation avoids
an earlier binary64 input bottleneck but does not turn long-double
arithmetic into arbitrary precision or intervals. Near a fold, the sign
of a residual smaller than actual accumulated error remains unresolved.

Numerical sign brackets at sampled radii give candidates. To apply the
intermediate value theorem rigorously requires validated signs and a
common continuous return domain across each bracket. No profile here
covers the entire radial domain or excludes further stationary points.
No global bound on origin cycles, K1 exclusion, or complete-component
kill follows from this review.
