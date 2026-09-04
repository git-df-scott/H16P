# Q4 weighted-lobe region: global geometry and a certified interior box

2026-09-04. **PROVED:** the strict universal weighted-lobe region is a
nonempty bounded analytic cell, parametrized by its three primitive roots.
**RIGOROUS COMPUTATION:** an explicit rational point and a closed coefficient
box of radius `10^-7` lie inside it. **PROVED:** that entire box nevertheless
cannot produce five original Q4 zeros, for any center parameter.

The proofs use the first-strike Stieltjes representation, auxiliary ECT
property, and primitive relation. Their inputs and conventions are in
[Q4_STRUCTURE.md](Q4_STRUCTURE.md) and
[Q4_ZERO_GEOMETRY.md](Q4_ZERO_GEOMETRY.md). The rational certificate is
replayed by [q4/q4_lobe_certificate.py](q4/q4_lobe_certificate.py).

## 1. Exact spaces and normalization

Let

\[
 F(t)={}_2F_1(1/6,5/6;1;t)>0,\qquad W(t)=tF(t),\qquad
 M(t)=\int_0^1\frac{\rho(u)}{1-tu}\,du.
\]

Write

\[
 Q=\operatorname{span}\{1,t,M,tM\},\qquad
 K_j(t)=\int_0^t W(u)e_j(u)\,du,
 \quad(e_0,e_1,e_2,e_3)=(1,u,M,uM).
\]

The primitive space \(\mathcal H=\operatorname{span}\{K_0,K_1,K_2,K_3\}\)
is four-dimensional: differentiation maps it injectively onto \(WQ\).
Every member extends continuously to both endpoints. At zero it vanishes to
order at least two; division by \(t^2\) gives an analytic extension with
value \(q(0)/2\). The logarithmic growth of \(F\) is integrable at one,
and \(M(1)=1\).

The normalized projective chart is

\[
 q(t)=A+Bt-1+(t-\eta)M(t),\qquad
 H_q=(A-1)K_0+BK_1-\eta K_2+K_3.
\]

The weighted-lobe region \(\mathcal L\) consists of these coefficients
for which \(q\) has three distinct simple roots \(x_1<x_2<x_3\) and,
writing \(\sigma\) for its sign before \(x_1\),

\[
 \sigma H_q(x_2)<0,\qquad \sigma H_q(x_3)>0,
 \qquad \sigma H_q(1)<0. \tag{L}
\]

Equivalently, \(H_q\) has three distinct simple interior zeros. This is
exactly the first-strike lobe criterion; there is no relaxation of its strict
inequalities.

## 2. A complete global parametrization

**Theorem.** Let

\[
 \Delta=\{(y_1,y_2,y_3):0<y_1<y_2<y_3<1\}.
\]

There is a real-analytic diffeomorphism \(T:\Delta\to\mathcal L\).
It sends the three prescribed primitive roots to their unique normalized
coefficient vector. In particular \(\mathcal L\) is nonempty, open,
connected, path connected, and contractible.

**Proof of the zero bound.** If a nonzero primitive has interior zeros of
total multiplicity \(m\), its endpoint value \(H(0)=0\) and Rolle's
theorem give at least \(m\) interior zeros of \(H'=Wq\), counted with
multiplicity. The auxiliary ECT theorem gives \(Z(q)\le3\), so
\(Z(H)\le3\). The same reasoning for the three-dimensional subspace
\(\operatorname{span}\{K_0,K_1,K_2\}\) gives at most two zeros, since
\((1,t,M)\) has strict Wronskians \(1,1,M''>0\).

**Proof of existence and exact uniqueness.** Three homogeneous evaluation
conditions in the four-dimensional primitive space have a nonzero solution.
Their matrix has rank exactly three: if the common kernel had dimension at
least two, imposing one additional derivative condition at an anchor would
leave a nonzero function with multiplicity at least four. This contradicts
the bound just proved. The three zeros are automatically simple and there
are no others. The coefficient of \(K_3\) is nonzero, because otherwise
we would have three zeros in the lower-dimensional subspace just excluded.
Normalize it to one. This gives exactly the displayed \((A,B,\eta)\)
chart, with no missing or exceptional anchor triples.

**Explicit finite formula.** Define the nonsingular matrix

\[
 E(y)=\begin{pmatrix}
 K_0(y_1)&K_1(y_1)&-K_2(y_1)\\
 K_0(y_2)&K_1(y_2)&-K_2(y_2)\\
 K_0(y_3)&K_1(y_3)&-K_2(y_3)
 \end{pmatrix}.
\]

Then

\[
 \boxed{\begin{pmatrix}A\\B\\\eta\end{pmatrix}
 =E(y)^{-1}
 \begin{pmatrix}K_0(y_1)-K_3(y_1)\\
 K_0(y_2)-K_3(y_2)\\K_0(y_3)-K_3(y_3)\end{pmatrix}.} \tag{A}
\]

Its determinant cannot vanish: a null vector would give three primitive
zeros with no \(K_3\) term. Formula (A) is an exact construction for every
specified anchor triple, not a generic linear-system heuristic.

**Proof of surjectivity and analyticity.** A coefficient in \(\mathcal L\)
has exactly three simple primitive roots, which form its unique preimage.
The matrix inverse in (A) is analytic on \(\Delta\). Conversely, the
implicit function theorem applies independently at each simple root, so
the ordered root map is locally analytic everywhere in \(\mathcal L\).
The inverse function theorem now makes (A) a global real-analytic
diffeomorphism onto the open set \(\mathcal L\).

The anchored Rolle step also proves that the corresponding \(q\) has
exactly three distinct simple roots: one in \((0,y_1)\), one in
\((y_1,y_2)\), and one in \((y_2,y_3)\). Thus the exact interlacing is

\[
 0<x_1<y_1<x_2<y_2<x_3<y_3<1.
\]

The first-strike inflection theorem then forces
\(1<\eta<54/31\). No \(\kappa\) occurs anywhere in this construction.

## 3. Explicit interior point, rigorous signs, and a certified box

**PROVED.** The primitive triple \((1/4,1/2,3/4)\) together with formula
(A) is an exact interior point of \(\mathcal L\). Its coefficients are
defined by specified finite integrals and a nonsingular three-by-three
matrix. Since \(\mathcal L\) is nonempty and open in \(\mathbb R^3\),
rational coefficient points are dense in it.

**RIGOROUS COMPUTATION.** The following frozen rational replacement of that
anchor construction is independently certified:

\[
\boxed{
 A=\frac{1243911778077}{10^{12}},\quad
 B=-\frac{86917392526}{10^{12}},\quad
 \eta=\frac{1460428426173}{10^{12}}.}
\]

The corrected strip holds exactly, and

\[
 q(0)=A-1-\eta/6=\frac{1014080763}{2000000000000}>0.
\]

The verifier uses rational arithmetic and directed integer rounding to
produce these enclosures:

| Primitive witness \(t\) | Lower bound for \(H_q(t)\) | Upper bound for \(H_q(t)\) |
|---|---:|---:|
| \(1/8\) | 0.000001441401383769815937769764 | 0.000001441401383769815937769765 |
| \(3/8\) | -0.000003720918085980812785719090 | -0.000003720918085980812785719089 |
| \(5/8\) | 0.000016780934607604637814181443 | 0.000016780934607604637814181444 |
| \(7/8\) | -0.000364906970505491567967035918 | -0.000364906970505491458511914236 |

These four alternating signs and the proved bound \(Z(H_q)\le3\) certify
exactly three distinct simple primitive zeros, one in each interval

\[
 (1/8,3/8),\qquad(3/8,5/8),\qquad(5/8,7/8).
\]

Anchored Rolle supplies exactly three simple auxiliary zeros, so the
first-strike lobe equivalence proves all three strict inequalities (L).
No approximate auxiliary-root position or quadrature at an unknown extremum
is used as proof.

These are **primitive signs**, not six alternating-sign witnesses for the
original Abelian integral. They establish lobe membership only. In fact the
entire certified box below has subsequently been excluded as an original
five-zero construction.

### Exact series and error control

Write

\[
 F(t)=\sum_{n\ge0}f_nt^n,\qquad f_0=1,\qquad
 \frac{f_{n+1}}{f_n}=\frac{(6n+1)(6n+5)}{36(n+1)^2}<1.
\]

The companion period has coefficients
\(k_n=-f_n/(6n-1)\) for \(n\ge1\). Therefore
\(tF(t)M(t)=F(t)-K(t)=\sum_{n\ge1}d_nt^n\), where
\(d_n=6nf_n/(6n-1)>0\), and

\[
\begin{aligned}
 H_q(t)={}&\sum_{n\ge0}f_n\left[
 \frac{(A-1)t^{n+2}}{n+2}+\frac{Bt^{n+3}}{n+3}\right]\\
 &+\sum_{n\ge1}d_n\left[
 \frac{t^{n+2}}{n+2}-\frac{\eta t^{n+1}}{n+1}\right].
\end{aligned}
\]

If \(S_N\) includes both sums through \(n=N\ge1\), decreasing positive
coefficients and a geometric majorant give

\[
\begin{aligned}
 |H_q-S_N|\le{}&\frac{f_{N+1}t^{N+1}}{1-t}\left[
 \frac{|A-1|t^2}{N+3}+\frac{|B|t^3}{N+4}\right.\\
 &\left.\qquad+
 \frac{6(N+1)}{6(N+1)-1}
 \left(\frac{t^2}{N+3}+\frac{|\eta|t}{N+2}\right)\right].
\end{aligned}
\]

All quantities are exact rationals at rational inputs. The certificate uses
\(N=256\); its largest analytic tail is below \(5.473\cdot10^{-20}\).
The first three tails are below \(10^{-30}\). The frozen evidence is
[q4/data/second_lobe_certificate.json](q4/data/second_lobe_certificate.json),
which records the verifier hash and exact rational inputs. The verifier has
a ten-second CPU ceiling, no third-party dependencies, and an observed replay
cost of about 0.04 CPU seconds.

### The full closed coefficient box

Allow independent perturbations of the displayed rational point with

\[
 |\Delta A|,\ |\Delta B|,\ |\Delta\eta|\le10^{-7}.
\]

Since \(0<M\le1\) and \(F(t)\le1/(1-t)\), for every witness
\(t\le7/8\),

\[
\begin{aligned}
 |\Delta H_q(t)|
 &\le10^{-7}\int_0^t uF(u)(2+u)\,du\\
 &\le10^{-7}\left[\frac{t^2}{1-t}+
 \frac{t^3}{3(1-t)}\right]
 \le\frac{1519}{1920000000}<8\cdot10^{-7}.
\end{aligned}
\]

Every frozen sign has absolute margin greater than \(10^{-6}\), so the
whole closed box preserves all four signs. Every point of this box lies
strictly inside \(\mathcal L\), with one primitive root in each of the
same three certified intervals.

**NUMERICAL discovery, separated from proof.** A small three-anchor matrix
solve selected the nearby decimal coordinates. Its rounding errors were not
used as certificate bounds: the frozen rational point and box were then
proved by the exact series above.

## 4. The region is bounded and has fixed sign orientation

For the normalized coefficient of \(tM\) equal to one, the first-strike
formula is

\[
 q''=M''(R-\eta),\qquad R' <0,
 \quad R(0)=54/31,\quad R(1^-)=1.
\]

Three simple roots require \(q'\) to have two zeros. The curvature changes
from positive to negative exactly once, so \(q'\) increases and then
decreases: its first zero is a minimum of \(q\), its second a maximum.
Consequently \(\sigma=+1\); the signs of \(q\) are \(+,-,+,-\).
In particular,

\[
 q(0)=A-1-\eta/6>0,\quad q(1)=A+B-\eta<0,
 \quad q'(0)=B+1/6-25\eta/432<0.
\]

The endpoint inequalities are strict: an additional endpoint zero, together
with three interior zeros, would force two interior zeros of \(q''\) by
ordinary Rolle, contradicting its single curvature change.

Moreover \(q'\) must be positive somewhere, whereas

\[
 q'(t)=B+M(t)+(t-\eta)M'(t)<B+1,
\]

since \(M<1\), \(M'>0\), and \(t-\eta<0\). Thus \(B>-1\).
These observations give the explicit bounded box

\[
 \boxed{\frac76<A<\frac{85}{31},\qquad
 -1<B<-\frac{49}{744},\qquad 1<\eta<\frac{54}{31}.}
\]

They hold even in the larger region where \(q\) has three simple roots,
without imposing the lobe inequalities. Therefore \(\overline{\mathcal L}\)
is compact in the ordinary normalized coefficient space; no projective
escape to infinity is hidden in the anchor parametrization.

## 5. Boundary mechanisms and exact equations

The closure is described rigorously as limits of (A) with anchor triples
approaching the boundary of the closed ordered simplex. Every finite
boundary point must exhibit at least one of the following events:

1. **A primitive root reaches the center.** Dividing by \(t^2\) before
   taking the limit gives \(q(0)=0\), hence the exact plane
   \(A=1+\eta/6\).
2. **A primitive root reaches the homoclinic endpoint.** Then
   \(H_q(1)=0\). The exact endpoint evaluation is
   \[
   H_q(1)=\frac{18}{85085\pi}
   (9061A+6289B-2431\eta-7242),
   \]
   so this boundary lies on an explicit affine plane. The integration-by-parts
   derivation is preserved in
   [q4/notes_certificate_second.md](q4/notes_certificate_second.md).
3. **Two interior primitive roots collide.** Then
   \(H_q(t)=q(t)=0\) at some \(0<t<1\). A generic such point has
   \(q'(t)\ne0\), and is an ordinary double primitive root. The two
   weighted-lobe faces are \(H_q(x_2)=0\) and \(H_q(x_3)=0\).
4. **Three interior primitive roots collide.** Then
   \(H_q(t)=q(t)=q'(t)=0\); the primitive root is exactly cubic, since
   multiplicity four is excluded. Its auxiliary contact is double, not
   triple.

Intersections of these events are the higher-codimension boundary cases.
These statements give an exhaustive list of possible mechanisms; they do
not assert that every point on any one of the displayed planes belongs to
the boundary. Sign inequalities and root counts select the relevant pieces.
No claim about a smooth compactified chart at the logarithmic endpoint is
needed for this classification.

To prove exhaustiveness, take any coefficient sequence in \(\mathcal L\)
converging to a boundary point and pass to a subsequence of its ordered root
triples in the compact closed simplex. If the limiting roots were still
three distinct interior points, uniform convergence and the multiplicity
bound would make them simple, and the limit coefficient would be in the
open set \(\mathcal L\), a contradiction. Interior collisions give the
successive derivative conditions by Rolle and local analytic convergence;
center escape uses the analytic quotient \(H/t^2\); endpoint escape uses
continuous convergence of \(H\) on the closed interval.

All four named mechanisms are naturally visible in the anchor coordinates:
letting \(y_1\to0\), \(y_3\to1\), \(y_1-y_2\to0\), or
\(y_2-y_3\to0\) stays inside \(\mathcal L\) until the limiting event.
Compactness ensures finite limiting coefficients. Every interior triple
collision can also be constructed exactly by the confluent version of (A),
with conditions \(H(t)=H'(t)=H''(t)=0\).

The generic weighted-lobe boundary is **not** the auxiliary discriminant:
when \(H\) has a double zero, \(q\) usually has a simple zero there.
The auxiliary triple-root cusp and its excluded open neighborhoods remain
disjoint from \(\mathcal L\), as the first strike proved. The primitive
cubic-contact boundary instead has an auxiliary double zero and an earlier
simple auxiliary zero, as anchored Rolle and the multiplicity bound show.

## 6. Useful one- and two-dimensional spines

Restricting primitive anchors to

\[
 (y_1,y_2,y_3)=(c-r,c,c+r),\quad
 0<r<\min(c,1-c),
\]

gives a two-dimensional exact subfamily entirely inside \(\mathcal L\).
The symmetric line

\[
 (y_1,y_2,y_3)=(1/2-r,1/2,1/2+r),\qquad0<r<1/2,
\]

gives a one-dimensional exact spine. Both produce coefficients by (A), so
neither requires preliminary lobe filtering.

The word spine has a precise topological meaning here. Map any anchor triple
to the symmetric triple with \(r=(y_3-y_1)/2\). Straight interpolation
between the original and projected triples stays in the convex set
\(\Delta\), and fixes the symmetric triples. Transporting this homotopy
by \(T\) gives a strong deformation retraction of \(\mathcal L\) onto
the stated one-dimensional spine. This does not imply that a separate
five-original-zero property, if nonempty, must meet the spine.

## 7. A sign consequence for reconstruction initial data

The exact reconstruction derived in
[Q4_RECONSTRUCTION_GEOMETRY.md](Q4_RECONSTRUCTION_GEOMETRY.md) gives the
universal center datum

\[
 Y_0=\frac3{1361360}
 (1326A+864B-2431\eta-102).
\]

**PROVED.** It is strictly
negative throughout \(\mathcal L\), and even throughout the full
three-auxiliary-zero region. Indeed, \(A<\eta-B\), \(B>-1\), and
\(\eta>1\) imply

\[
 1326A+864B-2431\eta-102
 <-1105\eta-462B-102<-745.
\]

Thus \(Y_0<-2235/1361360<0\). This sign is structural and does not come
from sampling the explicit point. Its significance for the final two-zero
allowance depends on the reconstruction equation and the other initial data.

## 8. The certified box is excluded for every original center parameter

**PROVED, independently audited.** The reconstruction analysis establishes
the necessary condition

\[
 \boxed{y_1>5/11}
\]

for five distinct original Q4 zeros, where \(y_1\) is the first universal
primitive root. This holds for every \(\kappa>1\); the first-peak proof
and exact rational estimates are in
[Q4_RECONSTRUCTION_GEOMETRY.md](Q4_RECONSTRUCTION_GEOMETRY.md).

Every point of the rigorously certified box has \(y_1<3/8<5/11\).
Therefore no point of that box can produce five original Q4 zeros for any
\(\kappa>1\). The universal lobe construction succeeded, and the next
reconstruction gate rigorously rejects it. This is an analytic exclusion of
the whole box, not an inference from failed numerical samples.

The surviving universal anchor domain is consequently restricted to

\[
 5/11<y_1<y_2<y_3<1.
\]

It is still nonempty, and formula (A) constructs a strict lobe-region point
for each such triple. Membership in this late-root region remains only a
necessary screen for the original five-zero target; it does not show that
the reconstruction can attain five zeros.

## 9. Scope of the completed lobe result

The weighted-lobe region is rigorously characterized as one bounded analytic
cell, with exact global root coordinates, connectedness, a one-dimensional
spine, rational points, and exhaustive boundary-event types. Its nonemptiness
is both proved analytically and witnessed by a fully certified rational box.

The original five-zero question remains open. The certified box is closed
as a candidate by the uniform first-root exclusion, while the late-root
portion of the cell remains available for the exact reconstruction tests.
No original five-zero witness, original three-simple-plus-double point, or
quadratic five-cycle field is supplied by this artifact.
