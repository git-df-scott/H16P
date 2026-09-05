# Strike 6: reduce the determinant sign problem to two anchor boundaries

2026-09-05. This is an analytic reduction relative to the audited Q4 canon.
It proves no universal sign on the two remaining boundaries. In particular,
the numerical boundary probes are not used in the proof.

## 1. A determinant with no baseline or center-mixture normalization

Write a primitive as `H_h=sum h_j K_j`, with coefficient order
`(K0,K1,K2,K3)`. Let

\[
k(t)=(K_0(t),K_1(t),K_2(t),K_3(t)),
\quad j_a(t)=(Y_{K_0}(t),\ldots,Y_{K_3}(t)),
\]

where each `Y` uses the exact original-family reconstruction and its
transported center data. These are linear functionals on coefficient
vectors, not freely chosen homogeneous solutions. Fix `a` and `0<r<1`.
For a second anchor `r<s<1`, put

\[
\Delta_{a,r}(s)=
\det\begin{pmatrix}k(r)\\k(s)\\j_a(r)\\j_a'(r)\end{pmatrix},
\qquad L_r(s)=K_0(r)K_1(s)-K_1(r)K_0(s)>0. \tag{B1}
\]

The denominator is positive because `x=K1/K0` strictly increases. Let
`V=K2-b K0-S K1` and `E=K3-e K0-T K1` be the two chord primitives with
anchors `r,s`. Thus `eta_V=-1`, `eta_E=0`. Their columns, together with
the first two coordinate vectors, give the unit-determinant coefficient
basis `(e0,e1,V,E)`. Evaluating (B1) in this basis yields

\[
\boxed{D_a(r,s):=p_a(r)
 [Y_E'(r)Y_V(r)-Y_E(r)Y_V'(r)]
 =p_a(r)\frac{\Delta_{a,r}(s)}{L_r(s)}.} \tag{B2}
\]

The Strike-5 baseline satisfies `E=B+eta_B V`; adding a multiple of `V`
does not change the Wronskian. Consequently its determinant is exactly

\[
\boxed{\mathcal K_a(r,s)=\lambda_c D_a(r,s),\qquad\lambda_c>0.} \tag{B3}
\]

Thus `D` and the earlier determinant always have the same sign. The
condition `P_B(r)>0` is still a separate necessary gate.

The initial determinant has a particularly simple exact value:

\[
D_a(0;r,s)=-\frac{Y_0(E)}{192}<0.
\]

Here the zero denotes evaluation time, with both anchors held fixed.
The strict positivity of `Y0(E)` is the Strike-5 chord bound. This initial
sign alone does not imply the sign at evaluation time `r`.

## 2. One more anchor zero would contradict Theorem N

As a function of `s`, the determinant `Delta(s)` is a single primitive in
`span{K0,K1,K2,K3}`. Its cofactor coefficient vector `h` satisfies

\[
H_h(r)=0,\qquad Y_h(r)=0,\qquad Y_h'(r)=0. \tag{B4}
\]

If all cofactors vanish, `Delta` and `D` are identically zero, and there
is no positive determinant to pursue. Suppose instead that `H_h` is
nonzero.

**Claim.** `Delta` has at most one zero in `(r,1)`, counted with
multiplicity. Any such zero is simple.

If there were two zeros there, including one double zero, `H_h` would
have at least three zeros counted with multiplicity, the first at `r`.
If its `K3` coefficient vanishes, this contradicts the two-zero bound
for `span{K0,K1,K2}`. Otherwise scale that coefficient to one. For three
distinct roots it is precisely a three-anchor primitive in the strict
lobe region. For a double root it is the finite confluent limit of such
primitives. In both cases Theorem N and its strict fixed-first-anchor
bound give

\[
\Phi_h(r)\le-\int_r^1 W_1(t)H_*(t)\,dt<0. \tag{B5}
\]

But (B4) makes `Z_h(r)=P_h(r)=0`, hence `Phi_h(r)=0`, a contradiction.

The confluent passage in this argument needs justification because a
non-strict limit of negative numbers could be zero. Here it cannot:
the right side of (B5) is a fixed strictly negative tail. Interpolation
coefficients converge at finite repeated anchors since the three lower
primitive functions form an extended Chebyshev system; a singular
confluent interpolation matrix would give a nonzero lower primitive
with three zeros counted with multiplicity. Reconstruction is continuous
on the compact interval up to `r`. The same argument handles a third
root at `1`, using endpoint continuity and the endpoint extension of
Theorem N. No singular lift/anchor limits are interchanged.

## 3. Two continuous boundary functions suffice for exclusion

The normalization in (B2) has finite limits at both second-anchor
boundaries. They are

\[
\boxed{D_{\mathrm c}(a,r)=
\frac{p_a(r)}{K_0(r)^2x'(r)}
\det\begin{pmatrix}k(r)\\k'(r)\\j_a(r)\\j_a'(r)\end{pmatrix},} \tag{B6}
\]

\[
\boxed{D_{\mathrm e}(a,r)=
\frac{p_a(r)}{K_0(r)K_1(1)-K_1(r)K_0(1)}
\det\begin{pmatrix}k(r)\\k(1)\\j_a(r)\\j_a'(r)\end{pmatrix}.} \tag{B7}
\]

They correspond respectively to `s -> r+` and `s -> 1-`. The first
formula is ordinary differentiation of the numerator and denominator,
both of which vanish simply in the denominator at `s=r`. The second
denominator stays positive. These formulas require no singular limit of
the baseline `H_{r,s,1}`; the chord pair `E,V` itself has finite limits.

**Boundary reduction theorem.** For every fixed admissible `(a,r)`,

\[
\boxed{D_{\mathrm c}(a,r)\le0\ \text{and}\ D_{\mathrm e}(a,r)\le0
\quad\Longrightarrow\quad
\mathcal K_a(r,s)\le0\ \text{for every }r<s<1.} \tag{B8}
\]

Unless `Delta` is identically zero, the interior conclusion is strict.

For strict negative boundary values, a positive interior value would
force two distinct interior zeros, contradicting the claim above.
Equality at a boundary cannot evade the argument:

- If `D_c=0`, then `Delta(r)=Delta'(r)=0`, so its primitive has a double
  zero at `r`. A later interior zero would supply the forbidden third
  zero counted with multiplicity.
- If `D_e=0`, then the primitive has a zero at `1`; an intervening zero,
  together with the root at `r`, contradicts the endpoint version of
  (B5).
- If both boundary values are zero, the double root at `r` and the
  root at `1` already give the contradiction for a nonzero cofactor
  primitive.

Likewise a zero of a nonpositive interior determinant would be a double
zero there, also impossible for a nonzero cofactor primitive. This proves
(B8) and its strictness claim.

Equivalently, a positive determinant anywhere on a second-anchor
interval forces at least one of the two boundary determinants to be
strictly positive. It cannot occur only in an isolated positive region
between two nonpositive boundary values. This is a sign-exclusion
reduction, not a monotonicity theorem for the determinant.

## 4. Exact remaining task and numerical controls

Combining (B8) with [the exclusion region](notes_exclusion_wedge.md), a
sufficient remaining theorem for closing route 4 is

\[
D_{\mathrm c}(a,r)\le0,\qquad D_{\mathrm e}(a,r)\le0
\]

for

\[
1-(7/22)^{3/2}<r<1,\qquad 1-1/\kappa_*<a<1.
\]

**Those two inequalities have not been proved.** They concern two
functions of two variables in place of the former unrestricted
three-variable determinant. Proving them would exclude every residual
fibre at its first required maximum and, with the inherited exclusions,
give a global bound of three distinct interior zeros for the Q4 family.
That would still not prove `H(2)=4` for all quadratic fields.

The diagnostic evaluator integrates `D` as an additional ODE variable,
using the exact identity

\[
\frac{dD}{dt}=
\frac{H_VY_E-H_EY_V}{1152t^2(1-at)^{3/2}(1-t)^{3/2}},
\]

and logarithmic time `x=-log(1-t)`. This avoids relying on cancellation
between large final Wronskian products. Direct product subtraction is
recorded as a consistency check. One `a=1` value is also checked by
independent high-precision quadrature of `P` and `Phi`. These numerical
checks validate the implementation at their stated points; they do not
certify either two-variable boundary inequality.
