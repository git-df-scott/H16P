# Exact restrictions for the KKL fold closure

This follows the already recorded continuation. It uses **zero ODE calls**,
does not restart the literature audit, and does not claim a component kill.
The algebraic replay is `theory_exact.py`; its exact output is
`theory_exact.json`. Two results are proved below: a stronger restriction on
where a fold can obtain its multiplier balance, and a precise impossibility
for the proposed one-sign scalar Dulac certificate at a fold.

## 1. One possible negative multiplier band throughout the upper strip

For the beta-zero KKL field, put

\[
u=1+x,\quad e=11c-5,\quad d=16-10c,\quad m=5(K+42)/e,
\]
\[
W=m+(2m+10)x+(m+111/5)x^2+(61/5-c)x^3,
\]
\[
N=\{du+(c+1)(21+dx)\}W-u(21+dx)W'.
\]

**Theorem.** For every \(1\le c\le8/5\) and \(K>0\):

1. \(N(x)>0\) for \(-1\le x\le0\).
2. \(N\) has at most two positive roots, counted with multiplicity.
3. Consequently its negative set on \(x>0\) is empty or one interval.
   For \(1<c<8/5\), a nonempty negative set is bounded by two simple roots.

These assertions concern the **multiplier integrand**, not the number of
limit cycles. They extend the inherited left-side sign restriction to the
whole upper strip and all positive K.

### Left-side proof

At K=0 substitute \(x=-X\), \(c=1+3C/5\) into \(5eN\).
Its bidegree-(4,4) Bernstein coefficient rows are

| X index | C=0 coefficient | C=1 coefficient | C=2 coefficient | C=3 coefficient | C=4 coefficient |
|---:|---:|---:|---:|---:|---:|
| 0 | 0 | 0 | 0 | 0 | 0 |
| 1 | 693 | 24759/20 | 35343/20 | 11403/5 | 27783/10 |
| 2 | 942 | 15807/10 | 114203/50 | 150129/50 | 92169/25 |
| 3 | 945 | 2835/2 | 82971/40 | 1184841/400 | 103194/25 |
| 4 | 900 | 1440 | 4551/2 | 71127/20 | 137592/25 |

Here column labels refer to Bernstein indices, not evaluation points.
All rows except index zero are strictly positive, proving positivity for
\(0<X\le1\), and N(0)=0 at K=0. Furthermore

\[
\partial_K N=\frac{5u^2}{e}
 \{cd\,u+5(2c+1)(c-1)\}\ge0.
\]

On \(0<u\le1\) this derivative is strictly positive throughout the
closed c strip, including its endpoints. Finally N(0)=5K. This proves (1).

### Positive-root proof

Write \(N=\sum_{j=0}^4 n_jx^j\). Here n0=5K>0 and

\[
n_4=\frac25(c-1)(5c-61)(5c-8)\ge0.
\]

The K slopes of n1 and n2 are positive on the whole strip, and
\(\partial_Kn_3=10c(8-5c)/e\ge0\). In its open upper endpoint define
the coefficient thresholds n_i=0 by

\[
k_1=\frac{21(25c^2-396c+305)}{25(5c^2-19c+5)},
\]
\[
k_2=\frac{4(-1375c^3+6529c^2-13445c+7625)}
{25(20c^2-43c+5)},
\]
\[
k_3=\frac{-14465c^3+47192c^2-56089c+22570}{50c(5c-8)}.
\]

The only coefficient pattern that could have four sign variations is
\((+,-,+,-,+)\). It is impossible:

- On \([1,11/10]\), \(k_2>k_3\). Indeed their difference is
  \(eA/[50c(5c-8)(20c^2-43c+5)]\), where
  \(A=21300c^4-100925c^3+160280c^2-102937c+22570\).
  Its Bernstein coefficients on this interval are
  \((288,1446/5,31909/120,34271/160,26451/200)\), all positive.
  Thus n2>0 forces n3>0.
- On \([11/10,8/5]\), \(k_2>k_1\). Their difference is
  \(-eB/[25(5c^2-19c+5)(20c^2-43c+5)]\), where
  \(B=2500c^4-19280c^3+46119c^2-53272c+24095\).
  The Bernstein coefficients of -B are
  \((17541/25,352493/200,559283/200,101556/25,141561/25)\),
  all positive. Thus n2>0 forces n1>0.

The denominators have the indicated positive signs on these intervals;
their quadratic factors are individually negative there. Descartes' rule
now gives at most two positive roots in the interior. At c=1 and c=8/5
the quartic term vanishes, but the same exclusion removes the only possible
three-variation cubic pattern \((+,-,+,-)\). This proves (2), including
zero coefficients by omission. Positivity at x=0 and the positive quartic
tail in the interior prove (3).

The inherited positive-weight multiplier identity is

\[
\log M=\oint\frac{N(x)\dot x^2}{5(1+x)W(x)^2}\,dt.
\]

Its denominator has no restoring-force exception in this strip. Indeed
m>50/3 and, in u coordinates,
\[
W\ge\frac{53}{5}u^3+\frac{79}{15}u^2-\frac{13}{5}u+1>0
\qquad(u>0),
\]
because the displayed quadratic has positive leading coefficient and
negative discriminant. Thus the origin is the only equilibrium on x>-1.

Therefore any attracting origin cycle or multiplier-one origin cycle in
this strip must cross that single right-side negative band. The left half
of an orbit cannot supply the negative contribution. This does **not**
compare the different weights furnished by distinct orbits. In particular,
it does not bound the number of origin cycles by two or exclude a third.

## 2. A scalar one-sign Dulac residual cannot certify across a true fold

Consider the exact reduction from the preceding theory record,

\[
\dot z=T-h(z),\qquad \dot T=-\kappa f(z)T-\mathcal R(z),
\quad h'=\gamma f,\quad \gamma+\kappa=1.
\]

For the KKL parameters, \(\kappa=c/(2c+1)\). Its divergence is exactly
-f. Suppose a C1 function Psi obeys

\[
F\cdot\nabla\Psi+4\kappa f\Psi=\Phi(z).
\]

This includes the proposed quartic recurrence, which the replay checks
symbolically without any restrictions on its four integration constants.

**Theorem.** If this field has a periodic orbit Gamma with multiplier
M=1, then a continuous Phi having one weak sign on the orbit's z-projection
must vanish on that entire projection. If Phi is real analytic on the
connected coordinate interval, it is identically zero there.

**Proof.** Parameterize Gamma by time over a period P and put
\(\mu(t)=\exp(4\kappa\int_0^t f(z(s))\,ds)>0\). The multiplier formula
gives \(\int_0^P f\,dt=-\log M=0\), so mu(P)=mu(0)=1. Along Gamma,

\[
\frac{d}{dt}(\mu\Psi)=\mu\Phi(z(t)).
\]

Integrating yields \(0=\int_0^P\mu\Phi(z(t))dt\). A continuous integrand
with one weak sign must vanish pointwise. A nonconstant periodic orbit's
z-projection is a nondegenerate interval: constant z would force T=h(z)
constant and hence an equilibrium. Analytic continuation proves the last
claim. This proof never divides by Psi and therefore also covers its zeros.

Thus the proposed escape of merely allowing additional isolated zeros of
Phi does not work at the exact fold: Phi must vanish on a whole interval.
For analytic scalar residuals it must vanish identically. In particular a
nontrivial residual positive except at the equilibrium is impossible there.
The observation is stronger than the generic warning that strict Dulac
conditions often imply hyperbolicity, because it identifies the complete
necessary degeneration for this particular scalar recurrence.

Off the fold the exact identity is instead

\[
(M^{-4\kappa}-1)\Psi(p)=\int_0^P\mu(t)\Phi(z(t))\,dt.
\]

For c>0 and a strictly positive integral, Psi is positive on an attracting
cycle and negative on a repelling one. Hence a certificate covering an S/U
pair must control a zero set of Psi between the cycles; its topology cannot
be omitted. If such certificates depend continuously on parameters and
extend with bounded coefficients to a multiplier-one orbit while retaining
a weak sign, their scalar residual must degenerate identically at the fold.

The next result removes even the identically-zero escape for the monic
quartic analytic ansatz. It uses a local argument suggested by the parent
lane and independently checked here.

### The identically-zero escape is impossible for the monic quartic

**Theorem.** Assume \(c>1/2\), \(m>0\), and \(K\ne0\). There is no
real analytic monic polynomial of degree four in T,
\[
\Psi(z,T)=T^4+C_3T^3+A_2(z)T^2+A_1(z)T+A_0(z),
\]
whose generalized Dulac residual is identically zero near the origin
equilibrium. Consequently, if a multiplier-one origin-surrounding cycle
exists, this ansatz cannot have a one-sign scalar residual Phi(z) on the
cycle's projection when its coefficients are analytic on the connected
interval extending to the origin.

**Proof of the local assertion.** Write \(q=4\kappa\). A zero residual
would give the density equation
\[
X\Psi=q(\operatorname{div}X)\Psi.
\]
Under an orientation-preserving analytic coordinate change with Jacobian
J, the transformed density is \(\widetilde\Psi=J^q\Psi\). Indeed
\(\operatorname{div}\widetilde X=\operatorname{div}X+X\log J\), proving
the transformed equation by the product rule. The positive analytic factor
J^q and the coordinate diffeomorphism preserve the vanishing order of Psi.

Make a finite analytic coordinate change putting the vector field into its
cubic Hopf normal form without a time change:
\[
X=\omega(-y\partial_x+x\partial_y)
  +a r^2(x\partial_x+y\partial_y)
  +b r^2(-y\partial_x+x\partial_y)+O(r^4).
\]
Here omega is nonzero and a is nonzero because the focus is of order one.
For completeness, in the Lienard chart at the origin the exact jets are
\[
f_z=-21/5,\quad f_{zz}=-(11+c)/5,\quad
g_z=m,\quad g_{zz}=20+m(1-c),
\]
and their first focal numerator is exactly
\[
f_zg_{zz}-f_{zz}g_z=2K\ne0.
\]
This also agrees with the inherited normalization
\(l_1=K/(8m^{3/2})\). Only its nonvanishing is needed here.

Let n be the lowest degree of the nonzero analytic density. At homogeneous
degree n the density equation says its leading term lies in the kernel of
rotation, so n is even and the term is \(A r^n\), A nonzero. At degree
n+1 the term is again radial (and hence zero when this degree is odd),
because the quadratic vector-field terms have been eliminated. At degree
n+2, average the equation over a circle. The rotation derivative averages
to zero, and the cubic radial part contributes \(naA r^{n+2}\). The
divergence of the cubic normal form is \(4ar^2\), so the right side
contributes \(4qaA r^{n+2}\). Therefore
\[
n=4q=16\kappa=\frac{16c}{2c+1}.
\]
The quartic is monic even when the equilibrium has a nonzero T coordinate;
its fourth T derivative there is 24. Thus its vanishing order n is at most
four. But
\[
\frac{16c}{2c+1}-4=\frac{4(2c-1)}{2c+1}>0,
\]
a contradiction. This proves the local assertion. The preceding
multiplier-one theorem forces a one-sign analytic scalar residual to vanish
identically on the connected interval, giving the claimed consequence.

The proof applies to every analytic monic quartic coefficient choice, not
only particular integration constants in the recurrence. It also covers
negative K as long as m>0 and K is nonzero. K=0 is deliberately excluded.
For c=1 the same statement uses any analytic local primitive h; the
logarithmic version of the primitive poses no local issue.

This closes the **monic-quartic scalar-residual certificate at the fold**
as a proof strategy in the stated regime. It does not close the fold
component or exclude K1. Off-fold certificates that fail to extend
analytically and with bounded coefficients to the fold, a different
two-variable residual, higher degree, or a separate orbit-comparison
argument are not excluded by this theorem. None is supplied here.

## 3. Closure status

No at-most-two origin-cycle theorem has been obtained. No entire connected
fold component has been excluded. No candidate field is supplied by this
zero-integration theory lane. The new exact facts narrow the difficult
comparison to a single right-side multiplier band and exclude the analytic
monic-quartic one-sign scalar-residual certificate at the fold, including
its identically-zero escape.
They are limited theorems and a precise proof-strategy obstruction, not
the user's stopping condition D.
