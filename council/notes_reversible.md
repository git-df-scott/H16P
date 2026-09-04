# Independent Astra review of the reversible proposal

Council review, 2026-09-04. This is an independent Astra subagent's analysis,
not a statement of agreement by Fable. No construction search was performed.
The local inputs were `CLAUDE_ROUTES_4AB.md`, `CLAUDE_THOUGHT_SESSION.md`,
`audit/claude_route4a_normal_form.py`, the post-Q4 frontier packet at
`de39ea7`, and `FASTRA_COUNCIL_2026_09_04.md` at `2136896`.

**Finding:** the proposed finite saddle-loop seed has an explicit integrable
annulus, but its alleged external focus is a saddle. The proposed `4+1`
completion fails at its starting geometry. This does not prove a five-cycle
bound for all reversible centers or for the closed annulus. The council's
parameter-count arguments do not provide such a bound.

## 1. Exact first integral and the finite annulus

Write the proposed normal form as

\[
 P_0=-Y(1+kX),\qquad Q_0=X+pX^2+qY^2.
\]

Put \(u=1+kX\), \(r=2q/k\), and \(b=p/k\), with \(k\ne0\).
On the half-plane \(u>0\), the equation for \(Z=Y^2\) is

\[
 \frac{dZ}{dX}+\frac{2q}{u}Z=-\frac{2X(1+pX)}u.
\]

Consequently the following is an **exact proved identity**:

\[
 H(X,Y)=u^rY^2+V(X),\qquad
 V(X)=2\int_0^X s(1+ps)(1+ks)^{r-1}\,ds.
\]

For \(r\notin\{0,-1,-2\}\),

\[
 V=\frac2{k^2}\left[
 \frac{b(u^{r+2}-1)}{r+2}
 +\frac{(1-2b)(u^{r+1}-1)}{r+1}
 -\frac{(1-b)(u^r-1)}r\right].
\]

The integral definition also covers the exceptional logarithmic cases. Its
positive integrating factor on this half-plane is

\[
 \mu=2u^{r-1},\qquad dH=\mu(Q_0\,dX-P_0\,dY).
\]

In particular, the origin is a nondegenerate center: \(H=X^2+Y^2+O(3)\).
This calculation proves the normal form's integrability; the old rotation
script itself only computes floating-point coefficients at three inputs and
does not certify an exact conjugacy from the Shi slice.

For the sector actually quoted in the proposal, assume

\[
 k>0,\qquad p<0,\qquad q>0,\qquad p+q<0.                 \tag{R1}
\]

All three displayed approximate seeds satisfy these strict inequalities.
The axis equilibrium is at \(X_s=-1/p>0\), and
\(u_s=1-k/p=1-1/b>1\). Let \(h_s=V(X_s)>0\). The left barrier has finite
energy \(V_L=V(-1/k+)\), because \(r>0\). Direct algebra gives

\[
 V_L-h_s=-\frac2{k^2}u_s^{r+1}
       \frac{r+2b}{r(r+1)(r+2)}>0.                       \tag{R2}
\]

Indeed \(V'=2X(1+pX)u^{r-1}\) is negative between \(-1/k\) and zero,
and positive between zero and \(X_s\). Thus there is a unique
\(X_L\in(-1/k,0)\) with \(V(X_L)=h_s\). The two graphs

\[
 Y_\pm(X)=\pm\sqrt{\frac{h_s-V(X)}{(1+kX)^r}},
 \qquad X_L\le X\le X_s,
\]

join into the finite homoclinic loop at the axis saddle. The periodic ovals
inside it are precisely the bounded components of \(H=h\),
\(0<h<h_s\), around the origin. Their closure stays in \(u>0\).
This supplies an annulus proof under (R1), without assuming the proposed
picture or performing a numerical phase portrait. The argument does not
classify the different boundary configurations when \(p+q\ge0\).

## 2. The equilibrium error and the actual topology restriction

Apart from the origin, the possible finite equilibria are

\[
 S=(-1/p,0),\qquad
 E_\pm=\left(-1/k,\ \pm\sqrt{\frac{k-p}{k^2q}}\right).
\]

The pair exists when the radicand is positive. At the axis equilibrium,
\(\operatorname{tr}J=0\) and \(\det J=k/p-1\). At either off-axis
equilibrium the Jacobian is triangular:

\[
 J(E)=\begin{pmatrix}-kY&0\\1-2p/k&2qY\end{pmatrix}.
\]

Its eigenvalues are \(-kY\) and \(2qY\), both real. Its determinant is
\(-2(k-p)/k\). The entire line \(X=-1/k\) is invariant. These equilibria
are **never foci**. In sector (R1), they and the axis equilibrium are three
hyperbolic saddles. The configuration is one center and three saddles, not
one center, one saddle and two antisaddles. Under sufficiently small general
quadratic perturbations these three saddles remain saddles.

There is also a direct check in the original Shi chart. With
\(m=5a,b=3l+5\), the Jacobian at the point advertised as the other focus,
\((0,1)\), is

\[
 \begin{pmatrix}5a&1\\3l+6&0\end{pmatrix},\qquad
 \det=-3(l+2).
\]

At the quoted center-curve value \(l\simeq-1.1835\), this is a saddle.

The primary PDF of [Zegeling, *Nests of limit cycles in quadratic systems*
(2024)](https://doi.org/10.1515/anona-2024-0012) was read directly from its
[published archive copy](https://d-nb.info/1332906729/34). Theorem 5.4,
page 29, permits only \((1,1)\) or \((n,0)\) for the four-real-singularity
case. Section 5's setup and Lemma 2.3 specify four real **finite** distinct
singularities, two being foci; the other two are separate equilibria.
“Generic” here identifies this equilibrium-count case. The statement is not
an unrestricted numerical upper bound on \(n\). Theorem 1.2 gives the wider
\((n,0)\) or \((n,1)\) distribution restriction. These restrictions exclude
a nearby two-nest \((4,1)\) construction retaining four real simple finite
equilibria. In the actual proposed seed, the external-focus error already
prevents that construction before applying this theorem.

Within this normal form an external weak focus cannot be supplied by an
off-axis equilibrium at any parameters. If the other axis equilibrium is a
nondegenerate antisaddle, reversibility makes it a second center, and it is
no longer the claimed axis saddle. A different two-center or infinity-boundary
proposal would require new geometry; it is not an automatic rescue of this
finite-axis-loop seed.

In particular, the council's sentence that a two-real-equilibrium variant
also necessarily needs \((5,0)\) does not follow from Zegeling. The general
distribution theorem allows \((4,1)\) when the four-real-equilibrium
hypothesis is absent. Such a variant would need an independently specified
second center and annulus geometry; none is supplied by the printed seed.

## 3. A derived first-order quotient, and what its dimension does not say

For arbitrary quadratic perturbations \((P_2,Q_2)\), the first energy
displacement, with the base orientation, is

\[
 M_1(h)=\oint_{H=h}\mu(P_2\,dY-Q_2\,dX)
 =\iint_{D_h}\left[\partial_X(\mu P_2)+
                  \partial_Y(\mu Q_2)\right]dX\,dY.      \tag{R3}
\]

The bounded oval domain \(D_h\) is symmetric in \(Y\). Only
\(1,X,X^2,Y^2\) from \(P_2\) and \(Y,XY\) from \(Q_2\) contribute.
After differentiation the resulting functions lie in the span of

\[
 J_f(h)=\iint_{D_h}u^{r-2}f(X,Y)\,dX\,dY,
 \qquad f=1,X,X^2,Y^2.
\]

There is an exact relation, not just a parameter-count suggestion:

\[
 J_X+pJ_{X^2}+(q+k)J_{Y^2}=0.                            \tag{R4}
\]

To prove it, apply the divergence theorem to
\(u^{r-2}Y(P_0,Q_0)\). The flux across an unperturbed oval is zero,
while its divergence is
\(u^{r-2}[X+pX^2+(q+k)Y^2]\).
Thus the first-order space has dimension at most three, and in (R1) one can
use \(J_1,J_X,J_{X^2}\) as spanning functions. This establishes a quotient
without claiming an unproved zero bound or independent basis proof.

For genericity, rotate \(x=Y,y=-X\). The Kapteyn parameters are

\[
 (\lambda_1,\ldots,\lambda_6)
 =(0,0,-q,2q-k,0,p).
\]

The other center-component intersections are detected by
\(p+q=0\), \(2q-k=0\), or simultaneously
\(k+3q+5p=0\) and \(p(q+2p)=0\).
The quoted strict sector excludes the first and third; require additionally
\(2q\ne k\) to exclude the Hamiltonian intersection. The displayed
approximate seeds satisfy this.

[Françoise–Gavrilov–Xiao, §4.1](https://arxiv.org/html/1610.07582v5#S4.SS1)
gives a three-dimensional bifurcation-function space at a smooth reversible
center, generated by first-order essential perturbations. Their scope is
the open period annulus; §1.1 explicitly separates closed-annulus and
polycycle cyclicity. Thus generic higher-order arcs do not automatically
create a new interior function space.

**Dimension is not a zero bound.** Even a one-dimensional analytic space
can contain a function with five simple zeros: take the span of
\(\prod_{j=1}^5(h-j/6)\) on \((0,1)\). A Chebyshev theorem, or another
argument controlling this particular space, is required to obtain two zeros
from three functions. This example refutes the council's general dimension
inference; it makes no claim to construct a quadratic vector field.
Likewise, counting transverse parameters does not prove a cyclicity cap.

Published two-cycle bounds for particular reversible families cannot be
silently extended to all exponents \(2q/k\). For example,
[Iliev–Li–Yu (2010)](https://www.aimsciences.org/article/doi/10.3934/cpaa.2010.9.583)
treat a reversible class with elliptic phase curves; and
[Liu–Li–Llibre (2021)](https://www.cambridge.org/core/journals/proceedings-of-the-royal-society-of-edinburgh-section-a-mathematics/article/abs/cyclicity-of-the-period-annulus-of-a-reversible-quadratic-system/E8F83BADE81E0E896B939BDA649D7EF5)
treat \(\dot x=y+ax^2,\dot y=-x\), \(a\ne0\).
Neither is a blanket cap for (R1). Conversely, FGX's generic first-order
reduction means the council's assertion that three reversible interior
cycles must occur only through higher-order or alien mechanisms is not
justified by a three-dimensional parameter count.

## 4. The symmetric Q4 claim needs a different launch description

The exact symmetric Q4 intersection is

\[
 \lambda_1=\lambda_2=\lambda_5=
 \lambda_4+5\lambda_3-5\lambda_6=
 \lambda_6(\lambda_3-2\lambda_6)=0,
 \qquad\lambda_4\ne0.                                  \tag{R5}
\]

It consists of two lines in the six-parameter Kapteyn chart, with convenient
nonzero representatives
\((0,0,1,-5,0,0)\) and \((0,0,2,-5,0,1)\).
Calling that intersection codimension five in this chart is legitimate.
It does not turn the center union into a smooth codimension-five center
component or imply five independent displacement functions.

[Buică–Giné–Grau, Theorem 6(vii)](https://arxiv.org/html/1406.7612v1)
gives essential order two at the symmetric Darboux center, with four
essential coefficients, choosing \(\lambda_{2,1}=1\):

\[
 M_2(h)=\lambda_{1,2}hB_1(h)+\lambda_{5,2}h^3B_3(h)
       +\lambda_{4,1}h^5B_5(h)+\lambda_{6,1}h^7B_7(h).
\]

[FGX, §4.2.2](https://arxiv.org/html/1610.07582v5#S4.SS2.SSS2)
independently describes the localized Bautin ideal there as

\[
 (\lambda_1,\lambda_5,
 \lambda_2(\lambda_4+5\lambda_3-5\lambda_6),
 \lambda_2(\lambda_3\lambda_6-2\lambda_6^2-\lambda_2^2)).
\]

Its exceptional divisor is \(\mathbb P^3\); these perturbations are in the
closure of the generic Q4 essential family. This supplies four interior
bifurcation functions, not a fifth one. It settles neither the boundary
cyclicity nor whether a particular infinity graphic has alien cycles.

There is also an immediate chart mismatch. The council's proposed fixed
slice

\[
 \lambda=(\tau,1+u,3,-10+v,w,1)
\]

has no symmetric Q4 point: symmetry forces \(u=-1\), whereas the Q4
condition becomes \(3-2-(1+u)^2=0\). At \(u=-1\) its value is one.
An explicit exact base point, boundary itinerary, and transported perturbation
chart are required before calling this a symmetric-Q4 experiment. A theorem
showing alien behavior in some two-saddle graphics does not prove it for
this unspecified symmetric base or prove simultaneous compatibility with
three interior cycles.

## 5. Defensible council position and bounded next derivation

The reversible proposal's strongest legitimate advantages are the exact
first integral, a proved finite annulus in (R1), and a three-function
interior reduction. Its advertised external-focus advantage is false.
Therefore the specific finite-loop-plus-external-fifth launch should be
withdrawn on exact geometric grounds. The single-annulus five-cycle question
has not been closed by this review.

If reversible work is explicitly retained, the bounded next derivation is
to choose one exact center on the intended slice, certify the conjugacy and
smooth-component conditions, and compute the endpoint asymptotics and
Wronskians of the three functions in (R4). This is a defined analytic task.
Its output must distinguish an interior Chebyshev bound from a bound on the
center-plus-annulus-plus-loop closure. No sweep or higher-order expansion is
justified merely by the assertion that “three parameters cannot make five.”

For a symmetric-Q4 reassignment, the first missing input is the exact base
and its original-coordinate graphic, followed by the correct second-order
essential chart and the endpoint return-map coefficients. Independence of
one proposed coefficient would be a useful mechanism test, not a five-cycle
certificate; dependence alone would not prove total cyclicity at most four
without a uniform zero-count argument for the composed return map.

Validation in this review: exact symbolic checks of the off-axis Jacobian,
(R2), (R4), and both representatives in (R5), completed in 0.25 CPU seconds
with one thread and a ten-second CPU fuse. No numerical shooting, parameter
scan, or construction attack was performed.
