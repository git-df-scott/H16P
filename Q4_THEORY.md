# Q4 theory: the exact infinitesimal target

**Audit date:** 2026-09-04. **Scope:** cycles bifurcating from the open period
annulus of a generic codimension-four quadratic center. This is not a bound for
all quadratic fields or for cycles born at the annulus endpoints.

## Origin and status

Żołądek's center classification calls the codimension-four component Q4.
Gavrilov--Iliev (2009) identified its Poincaré--Pontryagin generating functions
as complete elliptic integrals and proved an upper bound of eight zeros. Zhao
(Nonlinearity **24** (2011), 2505--2522) improved the upper bound to five and
proved that three are attainable. Both papers record the conjectural exact
bound as three. We found no primary source through the audit date resolving
that conjecture.

Thus Q4 is logically capable, under current theorems, of producing five cycles,
but no four- or five-zero Q4 integral is known. Zhao's five is an upper bound,
not a construction.

## Original quadratic system

In the normalized complex coordinate \(z=x+iy\), a generic Q4 center is

\[
 \dot z=-iz+4z^2+2|z|^2+(b+ic)\bar z^2,
 \qquad b^2+c^2=4,\quad c\ne0.
\]

Equivalently,

\[
\begin{aligned}
 \dot x&=y+(6+b)x^2+2cxy-(2+b)y^2,\\
 \dot y&=-x+cx^2+(8-2b)xy-cy^2.
\end{aligned}
\]

Set

\[
Y=cx-(2+b)y,\qquad \kappa=\frac4{2+b}>1.
\]

The rational first integral is

\[
 \mathcal H(x,y)=
 \frac{\left[8y(1+Y)-\frac23(1+\kappa Y^3)\right]^2}
      {\left[1-8y+\kappa Y^2\right]^3}.
\]

A permitted perturbation is a perturbation of this **original quadratic
field**:

\[
\begin{aligned}
 \dot x&=P_0(x,y)+\varepsilon X_2(x,y,\varepsilon),\\
 \dot y&=Q_0(x,y)+\varepsilon Y_2(x,y,\varepsilon),
\end{aligned}
\]

where \(X_2,Y_2\) have degree at most two and coefficients analytic in
\(\varepsilon\). This qualification matters: the changes below include a
double ramified cover and inversion. Adding an arbitrary quadratic polynomial
after those changes does not, in general, pull back to a quadratic polynomial
in the original coordinates.

## Elliptic Hamiltonian chart

Gavrilov--Iliev first take

\[
 X^2=1-8y+\kappa Y^2,\quad X>0,
\]

which is a double ramified covering rather than a birational projective
coordinate change. After translation, rescaling, and
\((x,y)\mapsto(x^{-1},yx^{-1})\), the level curves become the cubic elliptic
curves

\[
 H(x,y)=\frac23(\kappa-1)x^3-(\kappa-1)x^2y
       +\frac\kappa3y^3-y=h.
\]

The center is \((1,1)\), the Hamiltonian vector field is

\[
 \dot x=H_y=-1-(\kappa-1)x^2+\kappa y^2,
 \qquad
 \dot y=-H_x=-2(\kappa-1)x(x-y),
\]

and its clockwise ovals fill

\[
 -\frac23<h<-\frac{2}{3\sqrt\kappa}.
\]

The left endpoint is the center and the right endpoint is a homoclinic loop.
It is convenient to use

\[
 h=-\frac23\sqrt{\frac{s}{\kappa}},\qquad 1<s<\kappa.
\]

## The four-dimensional Abelian-integral space

For

\[
 I_{ij}(h)=\iint_{H(x,y)<h}x^iy^j\,dx\,dy,
\]

the first nonzero generating function is

\[
 I_\mu(h)=\mu_1hI_{00}+\mu_2I_{10}+\mu_3I_{01}
 +\mu_4\left(2I_{-1,0}+3\kappa hI_{-1,1}\right). \tag{Q4-I}
\]

All four coefficients are independent. An equivalent pre-reduction formula is

\[
 \iint_{H<h}\left[
 \widehat\mu_1x^3+\widehat\mu_2x^2y+\widehat\mu_3y^3
 +\widehat\mu_4\frac{\kappa^2y^4-1}{x}
 \right]dx\,dy.
\]

The hats are deliberate: the papers perform linear relabelings while reducing
this expression to (Q4-I). Coefficients must be transported through those
maps, not copied by name.

For a fixed \(\kappa\), multiplying \(\mu\ne0\) by a nonzero scalar does not
change its zeros; the coefficient target is therefore
\([\mu]\in\mathbb{RP}^3\).

## What the known bounds say

Zhao proves, counting multiplicity in the open interval,

\[
 Z(I_\mu;(1,\kappa))\le5.
\]

He also constructs an asymptotic coefficient hierarchy giving at least three
zeros tending to the center endpoint. Corollary 9(ii) says “at most three” in
its statement, but its proof, the abstract, and Theorem 1 all say and prove
**at least three**. The statement is a typographical error.

The five-zero Astra target is therefore exact:

> Find rational \(\kappa>1\) and rational or outward-enclosed
> \([\mu]\in\mathbb{RP}^3\) for which (Q4-I) has five distinct simple zeros in
> \(1<s<\kappa\), then construct the corresponding analytic quadratic
> perturbation arc in the original Q4 coordinates.

Five is the absolute maximum for this Q4 annulus problem under Zhao's theorem.
It is not an upper bound on global \(H(2)\).

## Primary sources

- L. Gavrilov and I. D. Iliev,
  [Quadratic perturbations of quadratic codimension-four centers](https://arxiv.org/abs/0811.4602),
  J. Math. Anal. Appl. 357 (2009), 69--76.
- Y. Zhao,
  [On the number of limit cycles in quadratic perturbations of quadratic codimension four centers](https://arxiv.org/abs/1011.2253),
  Nonlinearity 24 (2011), 2505--2522.
- I. D. Iliev, Perturbations of quadratic centers, Bull. Sci. Math. 122
  (1998), 107--161, doi:10.1016/S0007-4497(98)80080-8.
