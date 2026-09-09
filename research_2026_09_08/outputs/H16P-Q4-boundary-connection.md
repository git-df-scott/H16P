# Exact original Q4 annulus boundary

This calculation concerns the unperturbed quadratic field
\(P=-Y-3X^2+2XY+Y^2\), \(Q=X+X^2-4XY-Y^2\). It produces no limit-cycle counterexample. The base is the repository council's generic Q4 seed at modulus one. No novelty claim is made.

Use the previously checked quantities
\[
D=1-4Y+2(X+Y)^2,\quad
N=6Y(1-X-Y)-1+2(X+Y)^3,\quad h=N/D^{3/2}.
\]
On \(D>0\), \(\nabla h=6D^{-5/2}(Q,-P)\).

An exact parametrization of the finite boundary connection is
\[
X(s)=\frac{3+8s+4s^4}{12(1-2s^2)},\qquad
Y(s)=\frac{1-8s-24s^2-16s^3-4s^4}{12(1-2s^2)},
\quad -1/\sqrt2<s<1/\sqrt2.
\]
Put \(g(s)=2s^3-3s-2\). Direct substitution gives
\[
(P,Q)(X(s),Y(s))=\frac{g(s)}6(X'(s),Y'(s)),
\]
\[
D=\frac{2g^2}{9(2s^2-1)^2},\qquad
N=-\frac{2g^3}{27(2s^2-1)^3},\qquad h=-1/\sqrt2.
\]
Here \(g<0\): it decreases on the parameter interval, and its left endpoint value is \(\sqrt2-2<0\). Moreover
\[
X'(s)=-\frac{(2s^2+1)g(s)}{3(2s^2-1)^2}>0.
\]
The endpoint limits of \(X\) are respectively minus and plus infinity. Thus the curve is a smooth graph \(Y=b(X)\) on the entire real axis, with no self-intersection. Physical forward time runs from its positive-X end to its negative-X end. Its two asymptotic slopes are
\[
\lim_{X\to-\infty}b(X)/X=-1+\sqrt2,\qquad
\lim_{X\to+\infty}b(X)/X=-1-\sqrt2.
\]
These are the two saddle directions of the original compactification. The third projective direction, slope minus one, is not an endpoint of this connection. On the lower equator arc between the two endpoints there is no further singularity; at the downward vertical direction the angular component is positive. The equator arc therefore runs from the negative-X saddle to the positive-X saddle and closes the directed graphic.

The region below this graph is the center's full period-annulus domain, including its center before removing that point. Here is the global argument.

1. On the curve, \(X+Y-1=2(2s^3+1)/(3(2s^2-1))<0\) and \(D>0\). At fixed \(X\), the quadratic \(D(X,Y)\) decreases with \(Y\) below \(Y=1-X\). Consequently \(D>0\) throughout \(U=\{Y<b(X)\}\).
2. The unique curve crossing of \(X=0\) has parameter between \(-2/5\) and \(-3/8\). Its Y numerator is at least \(18027/20000>0\) on that bracket. Thus the origin lies in \(U\).
3. There is exactly one critical point of \(h\) on \(D>0\). Indeed, the coordinates \(a=X+Y-1,q=\sqrt D>0\) have nonzero Jacobian and give
\[
h=\frac{1-a^3+\tfrac32a(1+q^2)}{q^3}.
\]
The first critical equation gives \(q^2=2a^2-1\); the second then gives \(1+a=0\). Hence \(a=-1,q=1\), exactly the origin. Its level is minus one and its Hessian is \(6I\).
4. On the finite boundary, \(h=-1/\sqrt2\). At infinity within \(U\), the two graph asymptotes imply \(X+Y\le-c\sqrt{X^2+Y^2}\) for some fixed \(c>0\) outside a large disk. The leading terms \(N=2(X+Y)^3+O(R^2)\) and \(D=2(X+Y)^2+O(R)\) then give \(h\to-1/\sqrt2\) uniformly there.
5. An interior value above the boundary value or below minus one would therefore attain an interior extremum, contradicting the critical-point calculation. Equality with the boundary value at an interior point is likewise impossible. Thus \(-1\le h<-1/\sqrt2\) on \(U\), with equality at minus one only at the origin. Every intermediate level is nonempty, compact, and regular.
6. Every component of a compact regular level is a smooth circle. Since \(U\) is simply connected, its bounded interior lies in \(U\); the constant boundary value forces an interior critical point. It must enclose the origin. Two such circles at the same level would be nested and would force another critical point in the annulus between them. Hence there is exactly one circle at each intermediate level. The nonzero tangent vector field traverses it periodically.

Therefore \(U\setminus\{(0,0)\}\) is the annulus, and the displayed graph is its finite boundary connection. Its compactification adds the lower infinity arc and the two distinct saddle points. All of these are statements about the exact integrable base, whose periodic ovals are nonisolated.

The companion script checks the rational identities, endpoint slopes, crossing bracket, and critical-point reduction with exact arithmetic. Its initial run needed a rational expression simplified before symbolic substitution; the corrected final run passed. No ODE evaluations were used. The next unresolved question is the perturbed connection splitting and its compatibility with three compact Melnikov zeros in the same four-control direction, followed by the full boundary return analysis. Geometry alone supplies no extra cycle.

## Positive angular velocity on the entire annulus

Let \(\Omega=XQ-YP\). Substitution of the boundary parametrization gives
\[
\Omega=\frac{g(s)^2(2s^3+6s^2+3s+1)}{54(2s^2-1)^2}>0.
\]
The cubic factor attains its minimum on the parameter interval at \(s=-1+1/\sqrt2\), with value \(2-\sqrt2>0\). This identity and value have been checked exactly.

Because \(P<0\) on the graph and \(Q=P b'(X)\), positivity of \(\Omega\) implies \(b(X)-Xb'(X)>0\). Along any ray \((X,Y)=r(c,d)\), each graph intersection therefore has
\[
\frac{d}{dr}\{b(rc)-rd\}=\frac{Xb'(X)-b(X)}r<0.
\]
The ray starts inside the subgraph since \(b(0)>0\). Every intersection exits it transversely; two exits would require a reentry, which the same derivative sign forbids. Hence the annulus domain is star-shaped about the origin.

In the plane,
\[
\Omega=X^2+Y^2+(X+Y)(X^2-2XY-Y^2),
\quad \dot\theta=1+rA(c,d),
\]
where \(A=(c+d)(c^2-2cd-d^2)\). If a ray meets the boundary at a finite radius, the last expression is positive both at zero and at that radius, so it is positive in between. A ray that never meets the graph lies in the lower cone determined by its two asymptotes. On that cone \(c+d<0\) and \(c^2-2cd-d^2\le0\), whence \(A\ge0\). Its angular velocity is again positive at every finite radius.

Thus the physical polar angle increases strictly throughout the unperturbed open annulus. Since \(\partial_r h=6rD^{-5/2}\dot\theta>0\), every compact oval has one positive radius at each polar angle. This justifies the algebraic-oval quadrature's global coordinate choice; it does not validate its floating-point integration error.
