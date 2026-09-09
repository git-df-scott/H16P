# Boundary and interior compatibility

No counterexample was found. An exact energy calculation rules out adding a lower finite cycle along one leading perturbation direction used to create two upper and one lower boundary cycles. It does not rule out all five-cycle constructions.

The boundary proof uses coupled upper and lower displacement coefficients. Its simultaneous bounds concern cycles approaching the two hemicycles, not all finite cycles inside both annuli. I checked the coefficient relations and the simultaneous-count proof in [Marín–Villadelprat, The cyclicity of hyperbolic hemicycles, arXiv:2501.16924v1](https://arxiv.org/html/2501.16924v1). This source check does not constitute a full audit of the paper.

For fixed a in (-2,0), a≠-1, and b in (0,2), consider

x'=(b-2)/4+(1-b)y+ax²+by²+delta*x*(1-2ky),
y'=-2xy,

where k=sqrt[b(a+2)/(a(b-2))]>0. This direction makes the upper leading boundary coefficient vanish. Put t=|y| and sigma=sign(y). An unperturbed first integral is

H=t^a[x²+b*t²/(a+2)-sigma*(b-1)*t/(a+1)+(b-2)/(4a)].

Direct differentiation gives the exact identity

dH/dtime = 2*delta*t^a*x²*(1-2k*sigma*t).

For y<0 and delta>0, this is nonnegative and is positive whenever x≠0. A nonconstant periodic orbit cannot stay on x=0. Therefore integrating over a supposed lower periodic orbit gives a strictly positive change in H, a contradiction. This excludes a lower periodic orbit entirely in y<0 for this exact one-parameter direction.

It also makes the first Melnikov integral strictly positive on every regular compact lower orbit for an arc with this nonzero leading perturbation direction. Hence a finite lower cycle cannot be added by a small subordinate perturbation on a fixed compact annular region. This argument is not uniform up to the hemicycle, and does not exclude lower boundary cycles.

For the proposed two-upper/one-lower boundary configuration, two additional finite cycles would therefore have to lie in the upper annulus under this leading-direction assumption. The upper energy derivative changes sign, so its integral still needs analysis. Cases where the leading direction degenerates, or where the base parameters reach a special locus, are not covered by this obstruction.

A separate exact algebra check found a missing nu4 factor in the accessed HTML's displayed equation (23), relative to its preceding determinant. Restoring that locally bounded factor preserves the small-variable positivity argument. This is a formula discrepancy, not a refutation of the theorem and not an H16P counterexample.

The accompanying scripts verify both half-plane energy identities and the determinant calculation symbolically. An initial symbolic simplification failed because t had not been declared positive; adding the required t=|y|>0 assumption made both identities pass. All results and their scope are preserved in the JSON records.
