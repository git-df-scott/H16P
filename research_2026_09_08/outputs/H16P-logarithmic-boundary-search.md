# Large-radius reversible boundary returns

No counterexample was found. A new logarithmic-coordinate calculation found three numerical cycle brackets in one near-resonant quadratic field: two upper and one lower. It reaches trajectories far beyond the earlier Cartesian radius guard. The cycles remain uncertified, and unresolved returns remain explicitly marked.

## Field and result

The field is

x'=(b-2)/4+(1-b)y+a*x²+b*y²+e1*x+e2*x*y,
y'=e0-2*x*y,

with a=-1.01, b=1, e0=-1e-20, k=sqrt(0.99/1.01), e1=(0.01-1e-6)/(4k), and e2=-(0.01+1e-6)/2. The saved decimal e1 is the coefficient actually used in every return calculation; the same coefficients are used for both annuli.

| Annulus | Starting absolute-height bracket |
|---|---:|
| Upper | 288.14 to 1,348.88 |
| Upper | 1.533e13 to 7.175e13 |
| Lower | 1.613e17 to 7.551e17 |

The profile used 61 logarithmically spaced starting heights from 0.6 to 1e40 on each side, totaling 244 half-passages. Thirty-one pairs per side were unresolved. They do not imply absent cycles. The three brackets were replayed at a tenfold tighter relative tolerance, adding twelve half-passages. All retained their signs; the largest change in the logarithmic displacement was 1.89e-12. No fourth or fifth bracket was found.

## Coordinate calculation and checks

Use z=log|y| and v=asinh(x/|y|), separately on each half-plane. A positive time rescaling preserves the trajectory and its orientation. With u=sinh(v), c=cosh(v), d=exp(-z), sigma=sign(y), and C=(b-2)/4, the implemented equations are

dv/ds=[(a+2)u²+b+sigma*e2*u+(sigma*(1-b)+e1*u)d+C*d²-sigma*e0*u*d²]/(c²+d²),
dz/ds=(sigma*e0*d²-2u)c/(c²+d²).

The time factor is ds/dtime=|y|(c²+d²)/c>0. Each half-return must reach v=0 in the correct direction and within the bounded section between the invariant-line height and the corresponding unperturbed center height.

Twelve moderate-size paired returns agreed with independent Cartesian half-returns to 5.93e-13 in log height. Twelve unperturbed large-orbit controls, with starting heights 1e10, 1e30, and 1e60, all resolved. An independent first-integral calculation gave a maximum log-energy mismatch of 1.47e-11. Equal forward/backward heights alone would not have been an adequate control, because reversibility can force numerical symmetry.

The initial time rescaling matched moderate returns but failed ten of twelve large-orbit controls. Its script and failed records are retained. The revised rescaling passed all twelve. Across the two implementations, the controls used 96 logarithmic and 48 Cartesian half-passages; the actual-field profile and refinement used a further 256 logarithmic half-passages.

The remaining limitations are finite precision, unvalidated error, finite sampling, and unresolved returns beyond the valid numerical domain. This calculation reproduces a three-cycle boundary pattern and supplies a tested tool for more extreme scale comparisons; it does not establish five cycles or a global bound.
