# Third strike: an exact threshold-crossing path and a certified late point

2026-09-04. Scope: universal primitive-root construction and endpoint
conditioning. No five-original-zero result is claimed in these notes.
The original reconstruction is the next independent gate.

## 1. PROVED: a continuous path crosses the threshold exactly

Use precisely the canonical definitions

\[
 F(t)={}_2F_1(1/6,5/6;1;t),\quad
 q(t)=A+Bt-1+(t-\eta)M(t),\quad
 H(t)=\int_0^t uF(u)q(u)\,du.
\]

The threshold \(5/11\) concerns the **first simple zero of \(H\)**,
not an auxiliary \(q\) zero or an original \(I\) zero. The second-strike
first-Green-maximum theorem makes \(y_1>5/11\) necessary for five
original zeros. It is not sufficient.

Let \(T\) be the globally defined analytic three-anchor coefficient map
proved in `Q4_LOBE_REGION.md`. Define

\[
 \gamma(r)=T\left(r,\frac{1+r}{2},\frac{3+r}{4}\right),
 \qquad0<r<1.
\]

This entire path lies inside the strict weighted-lobe region. Its three
primitive roots are exactly the displayed anchors, so its first root is
identically \(r\). There is a rigorous, transversal threshold crossing at
\(r=5/11\); the derivative of the first-root coordinate along the path is
one. In particular

\[
 r_-=2/5<5/11<3/4=r_+.
\]

Both endpoint coefficient triples are exact, unambiguous finite-integral
expressions through \(T\). The global anchor theorem proves their roots
and lobe membership without numerical root estimation. The rational point
below is an independently certified nearby replacement of \(\gamma(3/4)\),
not an assertion that rounding preserves the exact anchor equations.

## 2. PROVED: closed primitive moments for evaluation near the loop

Put

\[
 J_n(t)=\int_0^t u^nF(u)\,du,\qquad n=0,1,2.
\]

The self-adjoint Gauss equation is
\([t(1-t)F']'=(5/36)F\). Integration by parts gives

\[
 J_0=\frac{36}{5}t(1-t)F',
\]
\[
 \boxed{J_n=
 \frac{t^{n+1}(1-t)F'-nt^n(1-t)F+n^2J_{n-1}}
 {n(n+1)+5/36}},\qquad n=1,2.
\]

For the four primitive basis functions
\(K_j=\int_0^t uF(u)e_j(u)\,du\), with
\((e_0,e_1,e_2,e_3)=(1,u,M,uM)\), the exact identities are

\[
 K_0=J_1,\quad K_1=J_2,
\]
\[
 K_2=6J_0-11J_1-6t(1-t)F,
 \qquad
 K_3=12J_1-17J_2-6t^2(1-t)F.
\]

They follow from \(M=1-6(1-t)F'/F\). Thus the evaluator needs only
\(F,F'\), three short moment recurrences, and one three-by-three solve.
It does not require quadrature or a long Taylor sum near \(t=1\).
At the endpoint use

\[
 J_0(1)=\frac{18}{5\pi},\quad
 J_1(1)=\frac{738}{385\pi},\quad
 J_2(1)=\frac{113202}{85085\pi}.
\]

The implementation in `q4_threshold_path.py` exposes

```
threshold_anchors(r)
primitive_basis_closed(t)
from_primitive_anchors_closed(anchors)
coefficients_from_r(r)
primitive_value_closed(t, coefficients)
```

The functions use the caller's mpmath precision with guard digits. Their
floating outputs are discovery/evaluation quantities, not interval proofs.
The identities themselves are exact. At `t=1/4,1/2,3/4`, independent positive
period-series evaluation agreed with the closed formulas within `5e-57`
at 75-digit working precision; this is a NUMERICAL implementation check.

## 3. RIGOROUS COMPUTATION: explicit rational coefficients beyond the threshold

Freeze

\[
\boxed{
 A=\frac{1210581187245108808}{10^{18}},\quad
 B=-\frac{125731163118386543}{10^{18}},\quad
 \eta=\frac{1212211767298108636}{10^{18}}.}
\]

These came from the single specified path point \(r=3/4\), whose exact
primitive anchors are \((3/4,7/8,15/16)\). The exact corrected strip
\(1<\eta<54/31\) and \(q(0)>0\) hold for the rational point.
The inherited exact-series evaluator, with its proved geometric tail
majorant and \(N=1024\), gives the following outward rational decimal
intervals:

| t | Lower bound for H(t) | Upper bound for H(t) |
|---|---:|---:|
| 23/32 | 0.000017278951342272812294524053 | 0.000017278951342272812294524054 |
| 13/16 | -0.000014889259450679365959381558 | -0.000014889259450679365959381557 |
| 29/32 | 0.000010156415192365481634868459 | 0.000010156415192365481634868460 |
| 31/32 | -0.000093181306082441272143107063 | -0.000093181306082441099333759810 |

Four alternating signs and \(Z(H)\le3\) prove exactly three distinct
simple primitive roots in

\[
 (23/32,13/16),\quad(13/16,29/32),\quad(29/32,31/32).
\]

In particular the first root is strictly greater than
\(23/32>5/11\). There is no additional primitive zero before \(23/32\),
so \(H\) is positive throughout \((0,5/11]\). A separate exact positive
enclosure at \(5/11\) is also stored in the frozen JSON.
The primitive-root bound and anchored Rolle prove all three strict weighted
lobe inequalities. None of these are original-integral sign witnesses.

The exact endpoint identity

\[
 H(1)=\frac{18}{85085\pi}(9061A+6289B-2431\eta-7242)
\]

also gives a strictly negative endpoint value. Its rational numerator and
all evaluator hashes are recorded in
`q4/data/third_threshold_certificate.json`.

### A whole late-root coefficient box is certified

Allow \(|\Delta A|,|\Delta B|,|\Delta\eta|\le10^{-8}\). On all four
witnesses, using \(M\le1\) and \(F\le1/(1-t)\),

\[
 |\Delta H(t)|
 \le10^{-8}\left[\frac{t^2}{1-t}+\frac{t^3}{3(1-t)}\right]
 \le\frac{122047}{3072}\,10^{-8}<10^{-6}.
\]

Every base sign margin exceeds \(10^{-5}\), so the entire closed box
preserves the three root intervals and stays strictly beyond the threshold.
The certificate uses only Fraction arithmetic and directed integer rounding,
with no floating-point error assumptions. One replay took 1.73 CPU seconds,
with a ten-second CPU ceiling and numerical libraries limited to one thread.

## 4. PROVED: matched endpoint asymptotics along the same path

Let

\[
 \varepsilon=1-r,\qquad L=\log(432/\varepsilon),\qquad
 U=A+B-\eta,\quad v=\eta-1,\quad E=2\pi H(1).
\]

For the three anchors, the endpoint distances are
\(\delta=c\varepsilon\), with \(c=1,1/2,1/4\). The exact endpoint
linear functional rewrites as

\[
 B=-\frac{17}{77}+\frac{9061}{2772}U
 +\frac{1105}{462}v-\frac{85085}{99792}E. \tag{E1}
\]

The universal endpoint expansions of \(F\) and its companion period give,
with \(L_\delta=\log(432/\delta)\),

\[
\begin{aligned}
 2\pi H(1-\delta)={}&E
 -\delta\{U(L_\delta+1)+6v\}\\
 &+\delta^2\left\{\frac{B+1}{2}(L_\delta+1/2)-3\right\}
 +\mathcal R,
\end{aligned}
\]

where, coefficientwise in this affine family,

\[
 \mathcal R=O\bigl((|U|+|v|)\delta^2L_\delta
 +(1+|B|)\delta^3L_\delta\bigr).
\]

For a direct check, the exact integrand can be written
\(tFq=(tA+t^2B-\eta)F-(t-\eta)K\), where the companion period is
\(K(t)={}_2F_1(-1/6,1/6;1;t)\). Expanding this expression and integrating from the endpoint
produces the displayed formula and bound.

Introduce scaled unknowns

\[
 U=\varepsilon L Q,\quad
 v=\frac{\varepsilon L}{6}\{V-(L+1)Q\},\quad
 E=\varepsilon^2L e.
\]

Divide the three exact anchor equations by \(\varepsilon^2L\), use
(E1), and pass to the limit. The resulting system is

\[
 e-cV+cQ\log c+\frac{30}{77}c^2=0,
 \qquad c=1,1/2,1/4.
\]

Its coefficient determinant is \((\log2)/8>0\). The exact equations are
linear in the scaled unknowns; their matrices converge to this invertible
matrix. Therefore the scaled unknowns converge, rather than their boundedness
being an additional assumption. Solving the limiting three-by-three system
gives

\[
 e\longrightarrow-\frac{15}{154},\quad
 V\longrightarrow\frac{45}{154},\quad
 Q\longrightarrow-\frac{45}{154\log2}.
\]

Consequently

\[
\boxed{\eta-1\sim\frac{15}{308\log2}\varepsilon L^2},
\qquad
\boxed{A+B-\eta\sim-\frac{45}{154\log2}\varepsilon L},
\]
\[
\boxed{H(1)\sim-\frac{15}{308\pi}\varepsilon^2L},
\qquad
\boxed{(A,B,\eta)\longrightarrow(94/77,-17/77,1)}.
\]

The leading individual coefficient shifts satisfy

\[
 B+17/77\sim\frac{1105}{462}(\eta-1),\qquad
 A-94/77\sim-\frac{643}{462}(\eta-1).
\]

This is a controlled endpoint limit along the same one-dimensional
threshold-crossing path. It is not a return to the excluded interior
auxiliary-cusp construction. Its leading primitive lobes are of size
\(\varepsilon^2L\), even though its coefficient deviations are of size
\(\varepsilon L^2\); that mismatch explains the numerical conditioning.

### NUMERICAL microscopic checks of the asymptotic scales

Two fixed evaluations, at \(\varepsilon=10^{-8}\) and \(10^{-16}\),
gave respectively

| Quantity | epsilon=1e-8 | epsilon=1e-16 | Proved limiting value |
|---|---:|---:|---:|
| \((\eta-1)/(\varepsilon L^2)\) | 0.05031104109 | 0.05924259353 | 0.07026112212 |
| \(U/(\varepsilon L)\) | -0.2817134558 | -0.3417506558 | -0.4215667327 |
| \(2\pi H(1)/(\varepsilon^2L)\) | -0.06325168565 | -0.07791223332 | -0.09740259740 |

The slow logarithmic convergence is expected from the retained terms.
These two evaluations and the three compact evaluator checks together used
about 1.77 CPU seconds. They are implementation diagnostics, not the proof
of the asymptotic theorem.

At very late anchors, converting the universal coefficients to binary64
can erase the primitive lobe margins: at \(\varepsilon=10^{-16}\) the
natural primitive scale is about \(10^{-31}\). High-precision evaluation
of a coefficient vector already rounded to binary64 cannot recover that
lost path information. The near-loop reconstruction must retain sufficient
coefficient precision or use scaled endpoint variables.

## 5. Remaining original-zero work

The threshold crossing is proved and an explicit rational late-root box is
certified. The endpoint scales are explicit. These achievements only open
the reconstruction gate; they do not force a positive first Green maximum
or establish that the remaining two-zero allowance is attained. The root
session is shooting the original PF/Green equation using this same path and
the permitted center parameter. Any original five-sign lead must still be
independently certified in the original normalization.
