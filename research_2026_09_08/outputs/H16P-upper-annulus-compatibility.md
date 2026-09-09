# Upper-annulus compatibility in the boundary-canceling direction

No counterexample was found. For fixed generic base parameters, the proposed leading direction cannot contribute any finite regular cycle from either annulus. The upper conclusion below is a deduction using a classical uniqueness theorem and the checked boundary expansion; it is stronger than the numerical sample alone.

Assume -2<a<0, a≠-1, 0<b<2, and a+b≠0. Set k=sqrt[b(a+2)/(a(b-2))]. Keep a,b fixed and take epsilon0=0, epsilon1=delta, epsilon2=-2k*delta. The invariant line is y=0. Write M(h) for the first energy-change integral on an upper regular oval of energy h, and W(h) for the strictly positive integral obtained by replacing (1-2ky) with 1 in its integrand.

Two literature inputs are used. Quadratic fields with an invariant straight line have at most one limit cycle; this established result is stated and referenced in the introduction of [Number of Limit Cycles for Planar Systems with Invariant Algebraic Curves](https://link.springer.com/article/10.1007/s12346-023-00746-7). The boundary expansion and coefficient coordinates are from equations (13)–(14) of [Marín–Villadelprat](https://arxiv.org/html/2501.16924v1). The original uniqueness proofs have not been independently rederived here.

## Derivation

1. Near the upper center y=1/2, M/W tends to 1-k. This is nonzero because k²-1=2(a+b)/[a(b-2)], and a(b-2)>0.

2. In the boundary coordinates, this direction has epsilon+=0 and epsilon-=4k*delta. The first displacement term therefore has the sign of a+b near the upper hemicycle. To convert displacement to energy change, use H_y<0 on the lower part of the positive-Y section, 0<y<1/2. Indeed, at x=0,

H_y=y^(a-1)*b*(y-1/2)*(y-(b-2)/(2b))<0.

Thus M near the boundary has sign -sign(a+b)=sign(1-k), the same sign as near the center. This sign statement uses the nonzero leading boundary coefficient and is not extended to a+b=0 or a=-1.

3. Suppose M has a zero on a regular upper oval. Normalize by W>0 and, if needed, reverse its sign so that F=M/W is positive near both ends of the open annulus. If F takes a negative value, choose three finite regular energies with signs positive, negative, positive. The corresponding small-delta return map then has at least two distinct isolated zeros, by continuity and analyticity. This would give two limit cycles in a quadratic field preserving y=0, contradicting uniqueness.

4. If F has a zero but never becomes negative, replace the leading coefficient epsilon1/delta=1 with 1+eta. The energy integral becomes M+eta*W exactly at first order. Taking eta of the appropriate sign and sufficiently small makes the normalized integral negative at the assumed zero while preserving positivity at two fixed regular energies on either side. The same two-cycle contradiction follows. This also handles an even-multiplicity zero. Identical vanishing is excluded by the nonzero center limit.

Consequently M has no regular finite zero under the stated assumptions. The lower integral was already proved strictly signed by direct energy monotonicity. This rules out adding two finite regular cycles to the three boundary cycles along this nondegenerate leading direction.

## Independent numerical evidence

Before completing the argument, I evaluated M/W at 61 energies for each of 26 shapes, using quadrature orders 256 and 512. All 1,586 values had the center sign. The smallest ratio (M/W)/(1-k) was 0.0005279; the largest quadrature-order change was 5.69e-15. The four sampled points with a+b=0 were excluded, as was the resonant line a=-1.

Nine representative full-trajectory integrals agreed with the quadrature ratios within 1.24e-12. The maximum relative return-coordinate error was 1.04e-11. These are numerical cross-checks, not interval enclosures. The nine cases were replayed once to record the initial coordinate and relative return error; they are nine distinct test cases, eighteen trajectory runs.

## Limitations

The argument excludes finite regular zeros for this particular fixed-parameter leading direction. It is not a bound on all cycles of nearby quadratic fields. It does not settle regimes approaching a+b=0, a=-1, the reflected exceptional locus, vanishing leading directions, or parameter-dependent cycles approaching a center or boundary. Such regimes need separate scaling and compatibility arguments. No five-cycle field has been verified.
