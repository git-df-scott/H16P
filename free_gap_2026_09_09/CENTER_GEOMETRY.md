# Exact affine reversibility behind a nearly flat return map

This calculation identifies a centre degeneracy relevant to the executed search. It produces no five-cycle field and makes no novelty claim. The formula is algebraic; proximity to the searched field is measured separately.

Consider

\[
P=ax^2+by^2+e_2xy+e_1x+(1-b)y+(b-2)/4,
\qquad Q=e_0-2xy.
\]

For −2<a<0, a≠−1 and b>0 define

\[
k=\sqrt{\frac{-a(a+2)}{b(2-a)}},\quad
e_1=-\frac{(a-1)k(1-b)}{a+1},\quad
e_2=\frac{2(a-1)(a+2)}{(2-a)k},\quad
e_0=\frac{k(a+b)(b-a-2)}{4b(a+1)^2}.
\]

Put A=2(a−1), B=e2, C=e1 and L(x,y)=Ax+By+C. Then

\[
R(x,y)=(x,y)-\frac{2L(x,y)}{A+Bk}(1,k)
\]

is an affine involution and reverses the vector field:

\[
F(Rz)=-D R\,F(z).
\]

The denominator is nonzero because A+Bk=4A/(2−a). This is an affine reflection, not necessarily a Euclidean reflection.

## Derivation

The divergence is L. Let u=(1,k) be the direction reversed by R, while the line L=0 is fixed. It suffices to show:

1. F restricted to L=0 is parallel to u.
2. The homogeneous quadratic part H satisfies H(u) parallel to u.
3. The divergence vanishes on the fixed line.

In affine coordinates (X,Y) with X=0 this line and the X direction parallel to u, the first condition says the second vector-field component vanishes identically at X=0. The second removes its X² term. It therefore contains only X and XY. The third then forces the first component's X and XY coefficients to vanish. That component contains only 1, Y, X² and Y². These are exactly the required reversal parities.

Condition 2 is the identity

\[
bk^2+Bk+(a+2)=0.
\]

On L=0 substitute x=−(By+C)/A into Q−kP. Its three coefficients vanish precisely when

\[
\frac{2B}{A}=k\left(\frac{(2-a)B^2}{A^2}+b\right),
\]
\[
\frac{2C}{A}=k\left(\frac{2BC(2-a)}{A^2}+1-b\right),
\]
\[
e_0=k\left(\frac{(2-a)C^2}{A^2}+\frac{b-2}{4}\right).
\]

Substituting Bk=A(a+2)/(2−a), k²=−a(a+2)/(b(2−a)), and C=−Ak(1−b)/(2(a+1)) verifies all four identities. The last identity reduces to the displayed formula for e0. Condition 3 holds because L is the divergence by definition. This proves the reflection identity for the stated family.

The implementation also checks the identity exactly, with rational arithmetic, for three rational (a,k) choices and the corresponding rational b. Six unisolvent points suffice to compare two bivariate quadratic polynomial fields in each check. These tests verify the implementation; the preceding algebra gives the parameter-family argument.

## Centre configurations near the search endpoint

An equilibrium on the reflection axis with positive determinant is a nondegenerate elliptic equilibrium. Nearby trajectories cross the axis and their reflected, time-reversed arcs close, giving a period annulus. Such a point is a centre, and its surrounding ovals are not isolated limit cycles.

On the explicit locus, equilibria on the axis are obtained from

\[
B y^2+C y+(a-1)e_0=0,\qquad x=e_0/(2y).
\]

Directed interval calculations over a∈[−1.908,−1.907] and b∈[0.365,0.366] verify two distinct nonzero real roots and strictly positive determinants for both equilibria. Trace zero follows from the exact symmetry-axis identity. Detailed outward enclosures appear in `center_interval_checks.json`.

The last accepted strong-search field has

```
a  = -1.9075810324974718
b  =  0.3656912644250558
e0 = -0.12283156115176493
e1 = -0.7137537873409198
e2 = -0.39156064615759056
```

Holding b and e0 fixed, the reversible formula has a solution with

```
-1.907584173654 < a_center < -1.907584173652.
```

Opposite strict interval signs of e0(a,b)−e0(target) establish existence inside that bracket. It lies inside the verified two-centre parameter box. Numerical projection gives changes of approximately −3.14116×10^-6 in a, −1.99808×10^-6 in e1, and +2.07571×10^-6 in e2. This is a nearby point at fixed b and e0, not a claim of globally minimal distance to the centre set.

Consequently, shrinking the strong branch's return mismatch cannot by itself be interpreted as convergence to an isolated outer double cycle. A centre has a continuum of return zeros. The numerical observation that several existing multipliers simultaneously approach one is consistent with this degeneracy.

## How the follow-up used this result

The `--divide-center` search fixes e0 away from zero and defines

\[
r=\frac{\|(e_0,e_1,e_2)-(e_{0,c}(a,b),e_{1,c}(a,b),e_{2,c}(a,b))\|}{|e_0|}.
\]

It optimizes D/(|e0|r), enforces r≥0.001, and still requires the original raw return signs for acceptance. This is distance normalization, **not** a proof that the displacement factors analytically by r or a certified Melnikov expansion. An implicit derivative check agrees with independent corrected finite differences. No positive outer witness preserving the four old brackets was found in that run.
