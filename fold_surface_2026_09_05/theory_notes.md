# Theory for the connected KKL fold continuation

Independent theory lane, 2026-09-05, continuing the numerical fold recorded at
commit `79001f7`. No ODE integrations were performed. Exact algebra is in
`theory_exact.py`; a separate 45-digit, nonvalidated center quadrature is in
`theory_melnikov.py`. Both machine records are saved beside the scripts.

**The component has not been excluded by theorem.** No at-most-two origin-cycle
bound is asserted here. The new results are (i) an explicit finite-amplitude
center/Melnikov organizer for the decreasing-K end, (ii) exact no-cycle gates
at and above c=1, and (iii) a proof that the two-parameter family is not a
rotated family around an origin cycle. These support continuation; they do
not satisfy the user's complete-component stopping condition.

## 1. The K=0 endpoint is organized by a finite center orbit

Let c0 be the root of

\[
J(c)=305+634c-11c^2-1000c^3
\]

in \((241/250,39/40)\), and put

\[
m=\frac{210}{11c_0-5},\quad
\sigma=\frac{5(2c_0+1)}{21},\quad \Delta=1+m\sigma^2.
\]

All constants in this section are evaluated at this exact algebraic point.
The inherited reversible involution has even/odd coordinates

\[
u=y-m\sigma x,\qquad v=x+\sigma y.
\]

Direct substitution, reducing coefficients modulo J, gives

\[
\dot u=-mv+Auv,\qquad \dot v=u+Bu^2+Cv^2,
\]

where

\[
A=\frac{10c_0\sigma m+5\sigma^3m^2-\sigma^2m-5\sigma m+100\sigma+11}
{5\Delta^2},
\]

\[
B=\frac{\sigma(5c_0-50\sigma^2-6\sigma-5)}{5\Delta^2},\qquad
C=\frac{5c_0\sigma^3m^2+11\sigma^2m+5\sigma m-50\sigma+5}{5\Delta^2}.
\]

In particular \(A+2C=21/5\). Normalize by

\[
X=Au/m,\qquad Y=Av/\sqrt m,\qquad \tau=\sqrt m\,t.
\]

Then

\[
X_\tau=-(1-X)Y,\qquad Y_\tau=X+bX^2+dY^2,
\quad b=Bm/A,\ d=C/A.
\]

Exact rational interval evaluation over an isolating interval for c0 proves
\(A>0\), \(-3/10<b<-29/100\), and \(101/100<d<51/50\).
The line X=1 is invariant at the organizer. This line need not remain
invariant at positive K.

For \(s=1-X>0\), define

\[
a_0=\frac{1+b}{d},\quad
a_1=\frac{2(1+2b)}{2d-1},\quad
a_2=\frac b{d-1},\quad E_\infty=-(a_0-a_1+a_2).
\]

The exact first integral is

\[
H(X,Y)=s^{-2d}(Y^2+a_0-a_1s+a_2s^2)+E_\infty.
\]

It vanishes at the origin, and its potential derivative is
\(2X(1+bX)s^{-2d-1}\). Consequently H has closed origin ovals for
\(0<E<E_\infty\). The potential increases to infinity at X=1 and has the
finite limit E-infinity as X tends to minus infinity; the signs of b,d above
justify these statements. Numerically, \(E_\infty\simeq18.2655313408\).
These assertions concern the exact center foliation, not limit cycles of
the perturbed field.

## 2. The finite Melnikov calculation reduces to two area moments

Keep the preceding linear transformation fixed at the center and perturb
the original coefficients by \((\delta m,\delta c)\). Let

\[
k=\sigma\sqrt m,\quad \nu=2d+1,\quad
a_m=\frac1{m\Delta},\qquad a_c=\frac{\sqrt m}{A\Delta^2}.
\]

The transformed perturbation is exactly

\[
P_1=a_c\delta c(X+kY)^2+a_m\delta m(kX-Y),\qquad Q_1=kP_1.
\]

The center integrating factor is \(\mu=2s^{-\nu}\), with
\(dH=\mu(Q_0dX-P_0dY)\). If \(\Omega_E\) is the interior of its E-oval,
the first energy displacement is

\[
M(E)=\iint_{\Omega_E}\operatorname{div}(\mu(P_1,Q_1))\,dX\,dY.
\]

Odd powers of Y integrate to zero. Define the weighted moments

\[
J_1(E)=\iint_{\Omega_E}2Xs^{-2d-2}\,dX\,dY,\quad
J_2(E)=\iint_{\Omega_E}2X^2s^{-2d-2}\,dX\,dY,
\]

and define Jy analogously with Y squared. Since the center vector field is
tangent to the oval and \(\operatorname{div}(\mu F_0)=0\), integration of
\(\operatorname{div}((Y/s)\mu F_0)\) gives the exact identity

\[
J_1+bJ_2+(d-1)J_y=0.
\]

Thus

\[
M(E)=a_m k\nu\,\delta m\,J_1+
a_c\delta c\,(A_1J_1+A_2J_2),
\]

\[
A_1=2\Delta-\frac{\nu k^2}{d-1},\qquad
A_2=-2\Delta+\nu-\frac{\nu k^2b}{d-1}.
\]

The exact parameter conversion at the center is

\[
\delta m=\frac5{e}\delta K-\frac{11m}{e}\delta c,
\qquad e=11c_0-5.
\]

Writing

\[
B_K=\frac{5a_mk\nu}{e},\quad
C_1=a_cA_1-\frac{11ma_mk\nu}{e},\quad C_2=a_cA_2,
\]

gives the usable expression

\[
M(E)=B_KJ_1\delta K+(C_1J_1+C_2J_2)\delta c.
\]

Where J1 is nonzero, let \(n(E)=J_2/J_1\). Along
\(c=c_0+\lambda K\), a first-order fold satisfies

\[
n'(E)=0,\qquad
\lambda=-\frac{B_K}{C_1+C_2n(E)}.
\]

This is the finite-amplitude condition missing from the earlier local-focus
analysis. A two-dimensional Melnikov space is not automatically Chebyshev:
the stationary ratio is precisely how a double zero can occur.

## 3. Reproducible, nonvalidated prediction for the decreasing-K end

The weighted moments have one-dimensional quadrature form

\[
J_j(E)=4\int_{s_-(E)}^{s_+(E)}
(1-s)^j s^{-2d-2}
\sqrt{-a_0+a_1s-a_2s^2-(E_\infty-E)s^{2d}}\,ds.
\]

The two endpoints are positive roots surrounding s=1. The script uses
bisection for them and a sine-squared substitution to regularize the square
root endpoints. It uses mpmath point arithmetic, not outward-rounded interval
quadrature; the results below are **numerical predictions**.

| Quantity | Numerical value |
|---|---:|
| Stationary center energy | 1.2051408317 |
| Original positive nullcline coordinate | 6.947539083 |
| J2/J1 | 3.586981585081526 |
| Limiting (c-c0)/K | 0.137109610961532 |
| Second derivative of J2/J1 in E | 0.0169942644 |

Centered derivative steps 0.0001 and 0.00005 give agreeing predictions; the
full figures and differences are in `theory_melnikov.json`. The finite-K fold
at \(K=1/512\), \(r\simeq6.94909\) is consistent with this finite center orbit
as K decreases. It is not consistent with simply identifying this pair with
a small cycle collapsing into the origin.

Here is the precise conditional continuation theorem: if interval analysis
establishes \(J_1\ne0\), \(C_2\ne0\), \(C_1+C_2n\ne0\), and a nondegenerate
stationary ratio \(n'=0,n''\ne0\) at an interior energy, then divide the
analytic displacement by K after substituting \(c=c_0+\lambda K\). The
implicit function theorem applied to the two fold equations gives a unique
local analytic fold branch in \((E,\lambda)\) through K=0. The stated
numerical quadrature has not proved those interval hypotheses and has not
proved that this is the only stationary ratio on the full annulus.

Even a proof of uniqueness for this first-order ratio would initially give
a local-in-parameter result on covered compact annuli. It would not alone
bound every positive-K cycle on an arbitrarily long continued component or
control the infinity endpoint. This is where the proposed full at-most-two
argument remains incomplete.

## 4. Exact no-cycle gates at c=1 and for 1<c<8/5

Use the inherited multiplier polynomial N and u=1+x. At c=1,

\[
N(u)=(5K-132/5)u^3-(48/5)u^2+6u+30.
\]

At \(K_N=6292/1125\), it factors exactly:

\[
N(u)=\frac2{225}(11u+15)(4u-15)^2.
\]

For \(K\ge K_N\), add \(5(K-K_N)u^3\). Therefore N is nonnegative on
u>0. Every origin-surrounding cycle would have strictly positive log
multiplier, using the exact positive-weight multiplier integral and N(0 in
the original x coordinate)=5K>0. But the origin repels, so an innermost
surrounding hyperbolic cycle cannot also repel: the displacement is positive
near the focus and cannot have its first zero cross from negative to positive.
Thus **there is no origin cycle at c=1, K>=6292/1125**.

For \(1<c<8/5\), write \(N=\sum_{i=0}^4 n_i u^i\). Here

\[
n_0=5c(c+1)(2c+1)>0,\qquad n_1=-c(c+1)(40c-43),\qquad n_4>0.
\]

Set

\[
K_2(c)=-\frac{8800c^5-58120c^4+92089c^3+39218c^2-53831c-31720}
{500(c-1)(2c+1)^2},
\]

\[
K_3(c)=-\frac{2200c^4-19095c^3+40912c^2-20175c-3050}
{50c(5c-8)}.
\]

The coefficient n2 increases in K and equals \(n_1^2/(4n_0)\) at K2; n3
increases in K and vanishes at K3. Therefore

\[
K\ge\max\{0,K_2(c),K_3(c)\},\quad K>0
\]

implies N>0 for every u>0: its quadratic part is nonnegative by the
discriminant test and its quartic contribution is positive. This gives a
conservative, explicit **no-origin-cycle region**, not merely a no-fold gate.

The restoring-force condition used here also holds throughout this strip.
Since \(m>50/3\), the polynomial W in u is bounded below by a positive cubic
plus \(5u^2-(13/5)u+1\); the quadratic has negative discriminant. Thus no
extra origin-side equilibrium invalidates the innermost-cycle argument.
These sufficient gates do not show that the currently traced fold reaches
them, and a numerical parameter cutoff before them is not a terminal theorem.

At \(c\ge8/5\), any remote equilibrium \(x=-s<-1\) has trace

\[
s\left[\frac{1+2c}{s-1}-\frac{16-10c}{5}\right]>0.
\]

This is an exact remote-stability boundary. It does not by itself exclude
three origin cycles, nor does it rule out a stable cycle around a remote
repelling focus.

## 5. Why rotated-family and raw N-sign arguments do not close the task

For any parameter direction in the beta-zero family,

\[
\delta F=(0,\delta c\,y^2-\delta m\,x),\qquad
\det(F,\delta F)=P(\delta c\,y^2-\delta m\,x).
\]

An origin cycle has a positive x maximum and negative x minimum. At both
transverse extrema P changes sign. For the determinant to keep one sign,
the second factor would have to vanish at both. If delta-c is nonzero, this
requires delta-m/delta-c to equal both a positive and a negative value,
namely \(y_+^2/x_+\) and \(y_-^2/x_-\). If delta-c is zero, delta-m must also
vanish. Thus **no nonzero parameter tangent in this two-parameter slice is a
rotated field on a neighborhood of a whole origin cycle**. A rotated-family
theorem cannot be invoked for the continued curve without a different,
explicitly justified transformation of the family.

Likewise, counting the sign changes of N only constrains the integrand of
each cycle's multiplier. Different cycles supply different positive weights.
Neither one negative band nor two positive roots of N prove an at-most-two
cycle theorem. The actual fold pair already defeats a naive uniqueness
conclusion from the near-center N pattern.

The remaining theorem obligation is a genuine comparison of distinct orbits,
or a complete validated continuation/endpoint analysis. The present exact
identities and no-cycle regions are admissible filters for that work; they
must not be used to terminate the user-authorized component strike early.

## 6. The proposed c=8/5 infinity endpoint

This section responds to the observed large-radius continuation, without
assuming that its endpoint has already been established. At c=8/5 the
positive-x infinity saddle directions are exactly

\[
z_\pm=-1\pm\sqrt{159}/3.
\]

Their radial/tangential eigenvalues in the signed positive-x chart give
Dulac exponents 5/6 and 6/5, respectively, so their product is one. For
1<c<8/5 the origin-side product can be written

\[
\rho(c)=\frac{\sqrt{40c-964/25}+2c-16/5}
{\sqrt{40c-964/25}-2c+16/5}<1,\qquad
\rho'(8/5)=10/\sqrt{159}.
\]

### A nonunit graphic coefficient does not prevent this endpoint

It would be incorrect to demand both rho=1 and C=1 as necessary conditions
for a fold curve to accumulate on an elementary graphic. The local model

\[
D(s)=\delta+C s^\rho-s
\]

has folds satisfying

\[
\log s=-\frac{\log(C\rho)}{\rho-1},\qquad
\delta=\frac{\rho-1}{\rho}s.
\]

Thus s can tend exponentially to zero as rho tends to one while C tends
to a nonunit positive constant. For the present approach from rho<1, the
compatible sign is C<1. The conditions to establish at a generic endpoint
are the actual interior connection (zero splitting) and rho=1; C then
selects the side and exponential scale. If s is asymptotic to a fixed
positive constant times 1/r, the leading prediction is

\[
(8/5-c)\log r\longrightarrow-\frac{\sqrt{159}}{10}\log C.
\]

This calculation assumes one effective splitting, finite nonzero regular
transition coefficient, and the elementary Dulac asymptotic. It is a
conditional mechanism, not a proof that this component reaches the graphic.
The higher-resonance case C=1 would require extra terms and extra analysis;
it must not be imposed on the generic endpoint.

### Exact separatrix seeds for a splitting equation

In v=1/x, z=y/x, the c=8/5 chart is

\[
z_\tau=-10+\tfrac65z+\tfrac35z^2-v(m+z^2),\qquad
v_\tau=-v(1+z)-zv^2.
\]

The radial invariant manifolds have jets

\[
z(v)=z_\pm+a_\pm v+b_\pm v^2+O(v^3),
\]

\[
a_\pm=\frac{5(m+z_\pm^2)}{11(1+z_\pm)},\qquad
b_\pm=\frac{5a_\pm z_\pm-8a_\pm^2}{16(1+z_\pm)}.
\]

The script verifies the invariance residual through degree two exactly modulo
the quadratic equation for z. Equivalently, in original coordinates,
\(y=z_\pm x+a_\pm+b_\pm/x+O(x^{-2})\).

The negative-slope saddle has the forward interior unstable branch; the
positive-slope saddle has the interior stable branch. Integrating the first
forward and the second backward to a common finite transverse section gives
a scalar splitting in m. For a surrounding connection, the negative-x
nullcline minimum is a possible matching section. Existence of these event
hits, their transverse derivatives, the tail remainders, and zero splitting
still require verification; the formal jets alone do not certify an orbit.

### No invariant conic supplies the connection

Any invariant conic with these two infinity directions has homogeneous part
\(y^2+2xy-(50/3)x^2\). Put

\[
F_2=y^2+2xy-(50/3)x^2+a x+b y+g.
\]

The cubic terms in the invariance equation force cofactor
\((16/5)(x+y)+t\). Matching the x-squared minus xy coefficient and the
y-squared coefficient forces \(b=20/11\), \(t=-10/11\). The constant term
then forces g=0, the linear terms force \(a=-200/121,m=-100/121\), and the
remaining x-squared coefficient equals \(-10180/363\ne0\). Hence **no such
invariant conic exists for any m**. This rejects a convenient exact ansatz;
it does not reject a nonalgebraic graphic connection.

## 7. Precise remaining gap in the paired-friction route

Let V be the inherited Lienard potential, E=V(z), and let x-plus(E)>0 and
x-minus(E)<0 be its two inverse branches on a shared energy interval. For
R=f/g define

\[
\mathcal B(E)=R(x_+(E))-R(x_-(E)).
\]

Writing u-plus/minus=1+x-plus/minus, its derivative is exactly

\[
\mathcal B'(E)=
-\frac{u_+^{3c+3}N(x_+)}{5x_+W(x_+)^3}
+\frac{u_-^{3c+3}N(x_-)}{5x_-W(x_-)^3}.
\]

If N is nonnegative on the left branch, then the second term is
nonpositive. Wherever N is also nonnegative on the right branch, the
whole derivative is negative (away from simultaneous isolated zeros).
Therefore stationary points of the paired friction can occur only while
the right branch crosses a negative-N band. This confines the difficult
comparison to a concrete energy interval.

It does **not** prove there is only one hump within that interval: the
positive first term and the positive magnitude of the second term can
cross repeatedly without a further derivative/comparison bound. Even a
proved zero count for this paired friction would need a verified applicable
cycle-count theorem, not only the no-cycle case of the inherited comparison
lemma. These are the two exact missing steps. They remain missing in this
theory lane; no unverified “two-cycle lemma” is being used.

## 8. A stronger exact reduction for a future Dulac certificate

The following identity applies to the entire beta-zero family, not only its
center or first-order perturbations. Assume c is nonzero and c is not one;
the c=1 expression is obtained separately with logarithms. Keep u=1+x and
write d0=16-10c (to distinguish it from the normalized center constant d).
For the inherited Lienard coordinate z, define

\[
\kappa=\frac c{2c+1},\quad \gamma=\frac{c+1}{2c+1},\qquad
F_{\rm raw}=u^{-c-1}\left[
-\frac{d_0u^2}{5(1-c)}+
\frac{(21-2d_0)u}{5c}+
\frac{d_0-21}{5(c+1)}\right],
\]

and h=gamma F-raw. This primitive satisfies h-prime=gamma f, with derivative
taken in z. If v is the Lienard velocity and T=v+h, then

\[
\dot z=T-h(z),\qquad \dot T=-\kappa f(z)T-\mathcal R(z).
\]

The remaining forcing is **exactly linear after removal of its known factor**:

\[
\mathcal R=(u-1)u^{-c}(a u+b),\qquad
a=-\frac{J(c)}{25(c-1)(2c+1)^2},\qquad
b=\frac{5K+(26-11c)a}{11c-5}.
\]

Indeed, the residual g-kappa f h is obtained by adding a cubic product to
W; its constant and linear coefficients cancel identically, and the two
remaining coefficients give the displayed a,b. The exact script checks
this identity without numerical sampling. At the double center, a=b=0.
For 1<c<8/5 and K>0, both a and b are positive, so R has exactly the sign
of x on the whole origin side. This reduction is stronger than merely
counting signs of the original quartic N.

A concrete quartic Dulac ansatz in these coordinates is

\[
\Psi=T^4+C_3T^3+A_2(z)T^2+A_1(z)T+A_0(z).
\]

For the generalized Dulac expression
\(\Phi=F\cdot\nabla\Psi+4\kappa f\Psi\), canceling every positive power
of T gives the recurrence

\[
A_2'=4\mathcal R-\kappa f C_3,
\]

\[
A_1'=hA_2'-2\kappa fA_2+3\mathcal R C_3,
\]

\[
A_0'=hA_1'-3\kappa fA_1+2\mathcal R A_2,
\qquad
\Phi(z)=-\mathcal R A_1-hA_0'+4\kappa fA_0.
\]

There are the constant C3 and three integration constants to choose. The
recurrence is derived directly here. It is the relevant kind of construction
for the generalized Dulac approach in [Cherkas–Artés–Llibre, section 2](https://ibn.idsi.md/sites/default/files/imag_file/Quadratic%20systems%20with%20limit%20cycles%20of%20normal%20size.pdf).
For this weak-focus family, Phi necessarily vanishes at the equilibrium;
strict positivity hypotheses at every point cannot be applied literally.
One would need a sign proof allowing isolated zeros and the requisite
zero-set/topology analysis of Psi.

No choice of these constants giving the required one-sign Phi over a whole
continued component has been proved in this session. Thus the displayed
recurrence is an explicit certificate construction problem, not an already
available at-most-two theorem. The reduction avoids the missing orbit-weight
comparison in section 7 if that certificate can actually be completed.

There is an additional obstruction at the observed fold itself. A strict
one-sign generalized Dulac expression, nonzero along an open part of each
cycle, normally forces the admitted cycles to be hyperbolic. It cannot be
asserted uniformly on the double-cycle parameter without checking its
degeneracy there. A viable route would need an off-fold certificate with
a controlled limiting argument, or a richer expression allowed to vanish
on the fold cycle. The isolated equilibrium zero alone does not resolve
this obstruction.

## 9. Status after the two-half shooting review

The independent source review in `theory_angular_review.md` verifies the
angular variational equations and the equivalent two-half fold system.
If A and B are the forward and backward half passage maps in log radius,
then F=A-B and G=log(A')-log(B') give F=G=0 at a fold. At a match the
full return multiplier is exp(G). Off-root, G detects stationary points
of F, not necessarily stationary points of the composed full-return
displacement. This distinction does not compromise fold continuation.

The finite-amplitude K=0 Melnikov endpoint and its nondegenerate critical
point remain numerical quadrature findings supported by exact reduction;
neither uniqueness nor an interval enclosure has been proved. The
potential upper endpoint at c=8/5 still requires the actual separatrix
connection condition and an analysis of the corresponding Dulac return.
As shown above, a limiting global coefficient C different from one is
compatible with the fold tending to the neutral exponent; C=1 must not
be imposed as an extra endpoint condition.

The exact algebra gives useful no-cycle regions and organizing
identities. It does **not** give an at-most-two theorem on the connected
component being followed. No continuation radius, sparse profile,
unresolved chart, or one-sided asymptotic observation supplies the
missing global exclusion. This theory lane has not completed the
component boundary analysis and has not established a K1 kill.

## 10. Exact half-passage exponents in the positive horizontal radius

The relation between the measured half-map residual G and a compact
return coefficient needs an explicit exponent conversion. Set

\[
c_* =8/5,\quad d=\sqrt{40c-964/25},\quad t=2c-16/5,\qquad
z_\pm=\frac{-6/5\pm d}{2(c-1)}.
\]

Here z is the projective slope y/x, not the initial log radius. In the
positive-x infinity chart v=1/x, the desingularized field has eigenvalues
p'(z_\pm)=\pm d tangent to infinity and -(1+z_\pm) transverse to it.
The forward half passage leaves the negative-slope saddle; the backward
half passage leaves the positive-slope saddle under reversed time. Their
positive saddle-passage exponents are respectively

\[
\nu_f=\frac{2(c-1)d}{d-t},\qquad
\nu_b=\frac{2(c-1)d}{d+t}.
\]

Both equal 6/5 at c=c_*. Their difference is exactly

\[
\Delta\nu=\nu_b-\nu_f
=\frac{10d(c_*-c)}{61-5c},\qquad
\lim_{c\to c_*}\frac{\Delta\nu}{c_*-c}
=\frac{12}{\sqrt{159}}.
\]

Their ratio nu_f/nu_b is the previously computed compact return exponent
rho. All these algebraic identities are checked in `theory_exact.py`.

To use these exponents for the actual family, assume the following
standard local-to-global passage structure is valid uniformly near the
candidate endpoint: both infinity saddles stay hyperbolic; the finite
branches reach the same transverse negative-ray section without further
singularities; and the saddle and regular transition maps admit the
following differentiated leading asymptotics, with positive finite
coefficients a_f,a_b. Writing ell=log r, let

\[
A(\ell;p)=A_\infty(p)-a_f(p)e^{-\nu_f(p)\ell}
 +o(e^{-\nu_f(p)\ell}),
\]

\[
B(\ell;p)=B_\infty(p)-a_b(p)e^{-\nu_b(p)\ell}
 +o(e^{-\nu_b(p)\ell}),\qquad p=(c,K).
\]

The same relative smallness is required after the ell differentiations
being used. These assumptions identify the actual local passages; the
eigenvalue calculation alone does not establish their global applicability.
They give

\[
F=A-B=\delta-a_f r^{-\nu_f}+a_b r^{-\nu_b}
 +o(r^{-\nu_f}+r^{-\nu_b}),\qquad
\delta=A_\infty-B_\infty,
\]

\[
G=\log A'-\log B'
=\log\frac{\nu_fa_f}{\nu_ba_b}+\Delta\nu\log r+o(1).
\]

Consequently a fold sequence tending to r=infinity and p=p_* must have
delta(p_*)=0, the actual separatrix connection. If in addition c tends
to c_*, put

\[
G_\infty=\log\frac{a_f(p_*)}{a_b(p_*)}.
\]

Then the fold equation G=0 implies the conditional asymptotic law

\[
\boxed{\ (c_*-c)\log r\ \longrightarrow\
-\frac{\sqrt{159}}{12}G_\infty.\ }
\]

If the differentiated remainder in G_ell is o(c_*-c), one also gets
G_ell/(c_*-c) tending to 12/sqrt(159). A uniform positive-power remainder
in 1/r, with its differentiated bounds, suffices along the nonzero-limit
scaling above. Merely knowing an o(1) remainder in G does not justify this
derivative ratio. For a branch approaching from c<c_*, a nonzero positive
limit requires G_infinity<0. If G_infinity=0 the displayed limit is zero
and does not determine the next scaling order.

This does not contradict the earlier coefficient 10/sqrt(159) for the
compact return exponent. At an exact connection the leading full return
in q=1/r follows from B(P(r))=A(r):

\[
q_{\rm return}=C q^{\rho}(1+o(1)),\qquad
C=(a_f/a_b)^{1/\nu_b}.
\]

At neutrality, log C=(5/6)G_infinity. Thus
-sqrt(159)log(C)/10 equals -sqrt(159)G_infinity/12. Equating log C directly
with the measured half-map G_infinity would miss this 5/6 factor.
Neither C=1 nor G_infinity=0 is required for the fold sequence to approach
the neutral connected field. These are conditional endpoint asymptotics,
not an interval proof of the connection or of component completeness.
