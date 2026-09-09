# Dividing out the exceptional center locus

No five-cycle counterexample was found. The finite-orbit obstruction extends across the exceptional center locus for the particular boundary-canceling direction studied here. This is a deduction using the classical invariant-line uniqueness theorem, analytic return maps, and the previously checked generic sign argument. It is not a bound on all cycles of nearby quadratic fields.

Take -2<a<0 and 0<b<2, put c=a+b and k=sqrt[b(a+2)/(a(b-2))], and consider the direction epsilon0=0, epsilon1=delta, epsilon2=-2k*delta. On an upper regular oval define

F(a,b,h)=M(a,b,h)/W(a,b,h),

where M is the first energy-change integral and W is the same positive integral with the factor (1-2ky) replaced by 1. An orbit coordinate on a common regular transversal may replace h whenever an energy normalization changes with parameters.

The prior argument established F/c<0 for a≠-1 and c≠0. Here is the extension.

1. On c=0, k=1, and the upper equilibrium is a center for every sufficiently small delta in this direction. This center locus is identified in Lemma 3.3 of [Marín–Villadelprat](https://arxiv.org/html/2501.16924v1). Thus F vanishes identically there. Analytic division gives an analytic coefficient G=F/c on each common compact regular annular region, including c=0.

2. At a=-1, the usual formula for the first integral has a removable parameter singularity after an additive normalization. In particular, the term involving y^(a+1)/(a+1) can be replaced by (y^(a+1)-1)/(a+1), whose limit is log y. Alternatively, analytic dependence of the original field and return map on a avoids this energy-coordinate issue. Therefore G also extends on regular compact orbits across a=-1.

3. By continuity from the dense generic parameter set, G≤0 everywhere in the stated parameter rectangle. It is not identically zero. Near the upper center,

F → 1-k, and (1-k)/c = -2/[a(b-2)(1+k)]<0.

On c=0 the limiting value is -1/[a(b-2)].

4. Suppose G vanishes at a finite regular orbit. Analyticity and non-identical vanishing give negative values on nearby regular orbits on both sides of such a zero. If c≠0, a small change of the leading epsilon1 coefficient adds a nonzero multiple of W to the first integral. It can create a positive value between these two negative values, giving two isolated return zeros for a sufficiently small perturbation.

5. For c=0, use a transverse analytic path c=t and delta=t. The actual local return displacement vanishes both when delta=0 and when c=0, so it factors analytically as delta*c*A. This factorization excludes an overlooked pure delta² term. Its leading term along the path is t²*A at the base point. A change of epsilon1 by eta*t² shifts this leading energy coefficient by eta*W. The same two sign changes can therefore be produced if G had a finite zero.

Both constructions preserve y=0 as an invariant straight line. They contradict the classical theorem that such a quadratic field has at most one limit cycle, stated and referenced in [Number of Limit Cycles for Planar Systems with Invariant Algebraic Curves](https://link.springer.com/article/10.1007/s12346-023-00746-7). Hence G is strictly negative on every finite regular upper oval.

## Numerical checks and remaining scope

I checked G by dividing the computed integral by c at c=±0.001 and ±0.0001, for seven base points on a+b=0, including (-1,1), and 41 energies per case. All 1,148 divided values were negative. These are nonvalidated finite differences, not the proof. The largest change between the two centered estimates was 6.72e-5; moving energy grids and quadrature error are not separated by that diagnostic.

At c=0 the exact direction itself has a center, not isolated upper cycles. The strictly signed divided coefficient shows why merely approaching that locus does not supply the two missing finite cycles in this leading mechanism.

This conclusion is confined to regular compact orbits and the stated coefficient direction or its divided leading term. It does not supply uniform estimates as an orbit approaches a center or hemicycle, nor does it cover every relative scaling of epsilon0, epsilon+, epsilon-, and c. Those parameter-dependent endpoint regimes remain open in this investigation. The original counterexample target remains unmet.
