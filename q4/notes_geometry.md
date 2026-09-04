# Q4 zero geometry — bounded first-strike notes

Date: 2026-09-04. Status labels below distinguish deductions from established
inputs from unresolved questions. No search or large computation was used.

## Inputs and conventions

The original four functions are exactly the repository's
`q4_integrals.py` basis

\[
 B=(hI_{00},I_{10},I_{01},2I_{-1,0}+3\kappa hI_{-1,1}),\qquad
 h=-\tfrac23\sqrt{s/\kappa},\quad1<s<\kappa.
\]

They span a four-dimensional real analytic space for fixed \(\kappa>1\).
Zhao's upper bound is five interior zeros counted with multiplicity.
The analytic changes of independent variable used here have nonzero derivative
inside the interval and therefore preserve multiplicity.

Primary input: [Zhao, arXiv:1011.2253](https://arxiv.org/pdf/1011.2253),
Picard–Fuchs equation (14), Proposition 3, and proof of Theorem 1.
The independent PF notes establish the positive Stieltjes representation used
in the auxiliary-space proofs below. These proofs concern the auxiliary
function, not the original four-integral space.

## The original projective problem: exact rank conditions

**Theorem-level deduction.** Let \(E_m\) be the matrix with rows \(B(s_i)\)
for distinct specified points. Three rows of rank three give a unique
projective null vector. Four rows generically have rank four and give no
nonzero vector. Four prescribed zeros are possible precisely when
\(\det E_4=0\); uniqueness additionally requires rank three.

For three anchors with rank three, set

\[
 D(t)=\det\begin{pmatrix}B(s_1)\\B(s_2)\\B(s_3)\\B(t)\end{pmatrix}.
\]

Up to a nonzero constant, this is the unique member of the original Q4 space
vanishing at the anchors. A fifth zero is not a single extra condition after
four freely chosen zeros: the fourth and fifth determinant conditions must
both hold. In particular, five distinct prescribed zeros are equivalent to
rank \(E_5\le3\), or to vanishing of all its four-row minors.

**Maximal-contact uniqueness lemma.** Any nonzero member satisfying interior
zero conditions of total multiplicity five is unique projectively. Moreover,
those multiplicities are exact and it has no further interior zeros.

Proof: if the common kernel of the five conditions had dimension at least two,
impose one further derivative condition at any prescribed point. A nonzero
member would remain, with total multiplicity at least six, contradicting the
upper bound. Exact multiplicities and absence of other zeros follow from the
same bound. Thus every attainable five-condition matrix has rank exactly
three; every six-condition matrix has rank four. This includes confluent
matrices containing successive derivatives at a repeated point.

## A feasible original-integral birth, and two impossible proposed births

Four simple zeros plus a double zero, and three simple zeros plus a triple
zero, each have total multiplicity six. Neither is possible for a nonzero
original Q4 integral. They cannot serve as starting configurations.

**Exact sufficient mechanism.** Three distinct ordinary zeros and one distinct
double zero have total multiplicity five and are sufficient. More concretely,
find four distinct points \(s_1,s_2,s_3,t\) for which

\[
 \operatorname{rank}
 \begin{pmatrix}B(s_1)\\B(s_2)\\B(s_3)\\B(t)\\B'(t)\end{pmatrix}=3.
\]

The uniqueness lemma proves that the first three zeros are simple and that
\(I''(t)\ne0\), without separate derivative tests. Because
\(B_1(t)=hI_{00}\ne0\) throughout the annulus, varying \(\mu_1\) unfolds
the double zero transversely. On one side of the fold it splits into two
simple zeros; the other three simple zeros persist by the implicit function
theorem. Thus this rank configuration would rigorously force five simple
zeros nearby. No such original-integral configuration has been exhibited.

The double-root locus for fixed \(t\) is the projectivized kernel of
\(B(t),B'(t)\). When those rows have rank two it is a projective line.
Taking its union over \(t\) gives the usual ruled discriminant surface in
coefficient projective space. At its ordinary points the local change in
interior zero count is two. Endpoint events need a separate endpoint chart.

## A shorter original-integral certificate

**Theorem-level deduction.** For exact admissible \(\kappa,\mu\), six rational
points in the annulus with rigorously alternating signs already prove exactly
five distinct simple zeros. The intermediate value theorem supplies one zero
in each adjacent interval. The multiplicity bound forces all five to be
simple and excludes every other interior zero.

Derivative enclosures and interval Newton remain useful optional checks, but
are not logically required for this certificate. The claim in the earlier
certification document that exactness also requires checking the complement
can therefore be shortened when Zhao's global bound is explicitly invoked.
This does not perform the subsequent original-field realization or return-map
certificate.

## The universal auxiliary curve lies on a quadric

Write \(t=(\kappa-s)/(\kappa-1)\in(0,1)\), with center \(t=0\), and use
\(w(t)=J_2/J_1\). The PF equation becomes independent of \(\kappa\).
After imposing the exact coefficient relation \(g(\kappa)=0\), normalizing
\(\beta_1=1\), and setting

\[
 b=\frac{\kappa-\beta_0}{\kappa-1},\qquad
 M(t)=\frac{1-w(t)}t,
\]

the auxiliary quotient has the form

\[
 f(t)=\frac{g(s(t))}{(\kappa-1)t}
     =A+Bt+(t-b)M(t)-1.
\]

Thus its full linear space is
\(V=\operatorname{span}\{1,t,M,tM\}\), independently of \(\kappa\).
The evaluation curve \([1:t:M:tM]\) lies on the explicit quadric
\(X_0X_3=X_1X_2\). Its hyperplane intersections are intersections of
\(M(t)\) with a fractional linear function. This exact quadric statement
applies to the auxiliary curve only: no such equation has been proved for
\([B_1:B_2:B_3:B_4]\).

## All four auxiliary Wronskians have strict signs

**Theorem-level deduction from the PF/Stieltjes representation.** Assume the
positive representation established in the PF notes,

\[
 M(t)=\int_0^1\frac{d\rho(u)}{1-tu},
\]

where the measure is positive and has support containing more than one
positive point. Its continuous positive density on \((0,1)\) is more than
sufficient. Differentiation under the integral is valid on compact
subintervals of \(0<t<1\).

For the ordered basis \((1,t,M,tM)\),

\[
 W_1=W_2=1,\qquad W_3=M''>0,\qquad
 W_4=3(M'')^2-2M'M'''<0.
\]

Indeed, \(M'=\int u/(1-tu)^2\,d\rho\),
\(M''=2\int u^2/(1-tu)^3\,d\rho\), and
\(M'''=6\int u^3/(1-tu)^4\,d\rho\). Strict Cauchy–Schwarz gives

\[
 \left(\int\frac{u^2}{(1-tu)^3}\,d\rho\right)^2
 <\left(\int\frac{u}{(1-tu)^2}\,d\rho\right)
  \left(\int\frac{u^3}{(1-tu)^4}\,d\rho\right),
\]

because \(u/(1-tu)\) is not constant on the support. Therefore the auxiliary
space is globally extended complete Chebyshev. There is no interior
Wronskian degeneracy in this space.

A direct proof avoiding any Chebyshev terminology is useful. Set

\[
 R(t)=t+\frac{2M'(t)}{M''(t)}.
\]

Then
\(R'=W_4/(M'')^2<0\). For
\(F=a+ct+dM+etM\), one has
\(F''/M''=d+eR\), which has at most one zero unless it vanishes
identically. Rolle's theorem gives at most three zeros of \(F\), counting
multiplicity; the degenerate case is affine and has at most one.

Equivalently, the basis \((1,t,w,v)\), \(v=(w-1)/t=-M\), has
\(W_3=w''<0\) and \(W_4>0\). This follows either by the constant basis
change or directly from its determinant.

## Exact fold and cusp parametrization of the auxiliary space

**Theorem-level deduction.** The normalized auxiliary function satisfies

\[
 f''(t)=M''(t)(R(t)-b).
\]

The endpoint expansions proved in the PF notes give

\[
 R(0)=\frac{54}{31},\qquad \lim_{t\to1^-}R(t)=1.
\]

Consequently every \(1<b<54/31\) has exactly one inflection point,
\(t_*=R^{-1}(b)\). The curvature is positive before it and negative after
it. Outside this range the curvature has one strict sign, so three interior
zeros are impossible. In original coordinates this is

\[
 \frac{54-23\kappa}{31}<\beta_0<1,
\]

which is nonempty for every \(\kappa>1\). The sign-reversed lower bound
in the old preparation material was not a valid finite-\(\kappa\) pruning
rule. Zhao's Proposition 17 and its displayed endpoint derivative agree with
the corrected bound; a later corollary contains the sign reversal.

Every double-root parameter lies on the explicitly ruled surface

\[
 B=-M(t)-(t-b)M'(t),\qquad
 A=1-(t-b)M(t)-Bt.
\]

Its triple-root edge has \(b=R(t)\). Along that edge

\[
 f'''(t)=M''(t)R'(t)<0,
\]

so each triple zero is exactly cubic. These formulas eliminate all three
normalized coefficients in favor of the single interior contact location.
For the equivalent expression \(f=A+Bt-w+bv\), they read

\[
 b=\frac{w''}{v''},\quad B=w'-bv',\quad A=w-bv-t(w'-bv').
\]

At the triple-root parameter, \(f'\) has its unique maximum zero at the
contact point; hence \(f\) is strictly decreasing away from that point and
has no other zeros. Perturbing

\[
 (A,B)\longmapsto(A-\lambda t_*,B+\lambda),\qquad\lambda>0,
\]

adds \(\lambda(t-t_*)\). For sufficiently small \(\lambda\), the cubic
zero splits into exactly three distinct simple zeros: the middle root stays
at \(t_*\), and the other two follow from the negative cubic coefficient.
The three-zero upper bound proves their simplicity and excludes extras.
This is an exact auxiliary bifurcation construction, not a five-zero
original-integral construction.

## Tiny exact calculation performed

The universal Riccati equation is

\[
 6t(1-t)w'=(1-t)(2w-1)-w^2,\qquad w(0)=1.
\]

A six-term rational recurrence evaluated with Python's `fractions.Fraction`
(single short process, no numerical scan) gives

\[
 w=1-\frac t6-\frac{25t^2}{432}
 -\frac{775t^3}{23328}-\frac{305675t^4}{13436928}+O(t^5).
\]

It checks \(R(0)=54/31\). It also gives
\(W_4(1,t,w,v)(0)=3705625/1451188224>0\), consistent with the global
integral proof. These are exact rational calculations; they do not constitute
an original-integral zero experiment.

## The first primitive excludes the entire local auxiliary-cusp route

Use distinct notation for the hypergeometric period and Zhao's forcing:

\[
 \mathscr P(t)={}_2F_1(1/6,5/6;1;t)>0,\qquad
 \mathscr H_q(t)=\int_0^t u\mathscr P(u)q(u)\,du,
\]

where \(q=A+Bt+(t-b)M-1\) is the normalized auxiliary quotient above.
Zhao's equation (28) and Lemma 11 say
\(\mathcal F_s=J_1g\) and \(\mathcal F(\kappa)=0\). Since
\(J_1(s(t))=J_1(\kappa)\mathscr P(t)\), direct integration gives the
exact identity

\[
 \mathcal F(s(t))
 =-(\kappa-1)^2J_1(\kappa)\mathscr H_q(t).
\]

The multiplier never vanishes. The period has only an integrable logarithmic
singularity at \(t=1\), and \(M\) extends continuously there. Therefore
\(\mathscr H_q\) extends continuously to \([0,1]\), and depends
continuously in uniform norm on the three coefficients. All of its zero
geometry is independent of \(\kappa\).

**New exclusion theorem.** For every fixed interior triple-contact point
\(t_*\), an open neighborhood of its normalized auxiliary coefficients has

\[
 Z(\mathscr H_q;(0,1))\le1,
 \qquad Z(I;(1,\kappa))\le3
 \quad\text{for every }\kappa>1.
\]

In particular, the explicit unfolding
\(q_\lambda=q_0+\lambda(t-t_*)\) has three simple auxiliary zeros for
small positive \(\lambda\), but cannot give five original-integral zeros.

Proof: at the cubic contact, \(q_0\) is strictly decreasing, positive on
\([0,t_*)\), and negative on \((t_*,1]\). Consequently

\[
 \mathscr H_{q_0}(t_*)
 =\int_0^{t_*}u\mathscr P(u)q_0(u)\,du>0.
\]

Choose a small closed neighborhood \(U=[a,c]\Subset(0,1)\) of \(t_*\)
on which \(\mathscr H_{q_0}>\eta>0\). Uniform continuity in coefficients
preserves \(q>0\) on \([0,a]\), \(q<0\) on \([c,1]\), and
\(\mathscr H_q>0\) on \(U\), for all sufficiently small coefficient
perturbations. Thus \(\mathscr H_q\) is initially strictly increasing
from zero, stays positive through \(U\), and is strictly decreasing after
\(c\). It has at most one interior zero, which is simple if present.
Zhao's inequality \(Z(I)\le Z(\mathcal F)+2\) proves the original
bound of three.

The coefficient neighborhood is independent of \(\kappa\), because
\(q,M,\mathscr P,\mathscr H_q\) are universal. Its size may depend on
\(t_*\); no uniform neighborhood as \(t_*\to0\) or \(t_*\to1\)
is asserted. This is a closed local route, not a proof that the whole
three-auxiliary-zero region is harmless.

## Exact weighted-lobe test for the remaining auxiliary region

**Theorem-level equivalence.** Suppose \(q\) has exactly three distinct
simple zeros \(0<x_1<x_2<x_3<1\), and let \(\sigma\in\{+1,-1\}\)
be its sign before \(x_1\). Then \(\mathscr H_q\) has three distinct
simple interior zeros if and only if

\[
 \sigma\mathscr H_q(x_2)<0,\qquad
 \sigma\mathscr H_q(x_3)>0,\qquad
 \sigma\mathscr H_q(1)<0. \tag{L}
\]

Proof: \(\mathscr H_q(x_1)\) has sign \(\sigma\), and the derivative
has signs \(\sigma,-\sigma,\sigma,-\sigma\) on the four successive
lobes. There is no zero on the first lobe. Each later lobe is strictly
monotone and contains exactly one simple zero precisely when its two
endpoint values have opposite signs. Those three conditions are exactly
(L). If one inequality is equality, a zero can sit at an extremum and be
double, or at the excluded endpoint. Thus the strict equivalence concerns
three distinct simple zeros, not just multiplicity three.

Equivalently, let the four positive weighted lobe areas be

\[
 L_0=\int_0^{x_1}u\mathscr P(u)|q(u)|\,du,\quad
 L_1=\int_{x_1}^{x_2}u\mathscr P(u)|q(u)|\,du,
\]
\[
 L_2=\int_{x_2}^{x_3}u\mathscr P(u)|q(u)|\,du,\quad
 L_3=\int_{x_3}^{1}u\mathscr P(u)|q(u)|\,du.
\]

Then (L) is

\[
 L_1>L_0,\qquad L_2>L_1-L_0,\qquad
 L_3>L_0-L_1+L_2.
\]

These are genuine magnitude requirements; three derivative crossings alone
do not supply them. They are independent of \(\kappa\) and can therefore
be used before selecting a center parameter.

For five distinct original-integral zeros, these strict conditions are
necessary. The center-anchored first derivative produces five distinct zeros
of \(G\). The next strict Rolle step can be verified directly from Zhao's
Proposition 3(ii): enclose those five roots in \([a,b]\Subset(1,\kappa)\),
choose \(s_0\in(1,a)\), and solve \(L_2y=0\) with
\(y(s_0)=0,y'(s_0)=1\). The Chebyshev property of the homogeneous solution
space says that \(y\) has no further zero, hence \(y>0\) on \([a,b]\).
Write \(L_2=\ell_2D^2+\ell_1D+\ell_0\), with \(\ell_2\ne0\), and
\(p=\exp\int\ell_1/\ell_2>0\). Direct differentiation verifies

\[
 L_2G=\frac{\ell_2}{py}
 D\!\left[py^2D\!\left(\frac Gy\right)\right].
\]

Five distinct zeros of \(G/y\) give four distinct zeros of its derivative,
and then three distinct zeros after the second weighted derivative. The
nonvanishing factors in Zhao's equation (24) therefore give three distinct
zeros of \(\mathcal F\). The bound
\(Z(\mathcal F)\le Z(g)\le3\) makes them simple. The anchored first
Rolle step from \(\mathcal F(\kappa)=0\) then gives three distinct,
hence simple, zeros of \(g\). This supplies the hypotheses of (L).
The conditions remain insufficient for five \(I\) zeros, because the final
second-order inverse step still has to attain its maximal extra oscillation.

## What is still open

1. No global Chebyshev property or sign theorem for the original basis
   \(B\) has been established here. Its Wronskians remain unresolved.
2. The local auxiliary-cusp route is excluded, even though its three
   auxiliary zeros exist. The remaining search requires finite separation
   and weighted lobe dominance, or a separately controlled endpoint limit.
3. No original-integral rank-three confluent configuration, five-zero
   numerical candidate, or five-zero certificate has been found.
4. The strongest next step is to determine whether any universal auxiliary
   coefficient satisfies all three weighted-lobe inequalities (L), outside
   the excluded cusp neighborhoods. Only a surviving shape should be lifted
   through the remaining original PF reconstruction.
