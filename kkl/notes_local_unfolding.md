# The KKL K=J=0 organizer is a double center

Independent bounded analytic addition, 2026-09-04. This note concerns only
the fixed-coefficient KKL family

\[
 \dot x=(1+x)y+x^2,\qquad
 \dot y=-10x^2+\frac{11}{5}xy+cy^2-mx,
\]

with \(\beta=0\), \(1/2\le c\le3/2\), \(m=-\alpha>0\),
\(K=m(11c/5-1)-42\). It extends the analytic calculation in
`notes_lienard.md` without changing that checkpointed file. No orbit
integration, numerical discovery, sweep, or CPU algebra experiment was
performed.

**Result.** The zero of the second weak-focus coefficient on \(K=0\) is
an exactly reversible center, not an order-three weak focus. Its unique
other real finite equilibrium is also a center. On the \(\beta=0\)
slice, sufficiently small perturbations near this organizer have at most
one small nonzero origin cycle and no local stable/unstable fold pair.
The stable-remote-focus condition excludes even a cycle shrinking into
the origin within that cone. A separate finite-amplitude pair remains open.

## 1. Exact parameter and its position relative to infinity gates

Use the already derived polynomials

\[
 e=11c-5,\quad J(c)=305+634c-11c^2-1000c^3,
 \quad m_0(c)=210/e.
\]

The function \(J\) is strictly decreasing on \([1/2,3/2]\), and

\[
 J(241/250)=25281/2500>0,
 \qquad J(39/40)=-11333/800<0.
\]

Let \(c_*\) be its unique zero. Then

\[
 \frac{241}{250}<c_*<\frac{39}{40}<1,
 \qquad m_*=\frac{210}{11c_*-5}.                       \tag{U1}
\]

This is an exact algebraic point, not a rational field proposed for
certification. It has \(K=0\), outside the experimental precursor margin
\(K\ge1/64\).

The organizer lies beyond the \(c=241/250\) change in real infinity
directions, and before the vertical degeneration at \(c=1\). In the
inherited vertical chart the extra directions solve
\(10u^2-(6/5)u+(1-c)=0\). They are distinct and real at \(c_*\).
The vertical eigenvalues \((1-c_*,-c_*)\) are nonzero and of opposite
sign. Thus the finite-center condition is not a vertical infinity
degeneration. A path from the last tracked \(c<241/250\) region still
needs its infinity itinerary checked when crossing that boundary.

## 2. Explicit reversal in the original coordinates

At \(K=0\), normalize the linear part by
\(\zeta=x-iy/\sqrt m\) and \(\tau=\sqrt m\,t\). Then

\[
 \frac{d\zeta}{d\tau}=i\zeta+A\zeta^2+B\zeta\bar\zeta+C\bar\zeta^2,
\]

\[
 A=\frac14\left[\frac{16}{5\sqrt m}+i(c+1+10/m)\right],
 \quad B=\frac12\left[\frac1{\sqrt m}+i(10/m-c)\right],
\]

\[
 C=\frac14\left[-\frac6{5\sqrt m}+i(c-1+10/m)\right].
\]

Substitution of \(m=210/e\) gives

\[
 A=\frac85\bar B,\qquad
 \operatorname{Im}(\bar B^3 C)=\frac{J(c)}{35280\sqrt m}. \tag{U2}
\]

For a hand check of the second identity, put
\(a=2c+1\), \(P=21/(5\sqrt m)\), \(S=P+ia\),
and \(T=-6/(5\sqrt m)+2i(16c-13)/21\). Then
\(\bar B=5S/42\), \(C=T/4\), and
\(\operatorname{Im}(S^3T)=2PJ/125\).

At \(c=c_*\), set \(\ell=-B/\bar B\); \(|\ell|=1\).
The reflection \(\zeta\mapsto\ell\bar\zeta\) reverses time. This can
be checked coefficient by coefficient: the required quadratic identities
are

\[
 A\ell=-\bar A,\quad B=-\ell\bar B,\quad C=-\ell^3\bar C,
\]

and all follow from (U2). This also matches the standard reversible
center equations in [Françoise–Gavrilov–Xiao, equation (25)](https://arxiv.org/html/1610.07582v5#S4),
but the displayed verification does not require treating that classification
as a black box.

In original coordinates put \(\sigma=5(2c_*+1)/21\). The involution is

\[
 \boxed{\quad
 M=\frac1{1+m_*\sigma^2}
 \begin{pmatrix}
 m_*\sigma^2-1&-2\sigma\\
 -2m_*\sigma&1-m_*\sigma^2
 \end{pmatrix},\quad F(Mz)=-M F(z).
 \quad}                                                \tag{U3}
\]

Here \(M^2=I\), \(\det M=-1\), and its fixed line is
\(x+\sigma y=0\). The original time orientation in (U3) is genuinely
reversed; the preceding time normalization was positive.

The origin has eigenvalues \(\pm i\sqrt{m_*}\). A reversible planar
singularity with these nonzero imaginary eigenvalues is a center: nearby
arcs between successive crossings of the reflection axis join with their
time-reversed reflections into closed orbits. Thus every origin focal
coefficient vanishes at (U1). The next coefficient is not a nonzero
third-order focus coefficient. This assertion concerns the exact organizer;
it does not say higher return coefficients vanish at nearby noncenter
parameters or impose a bound on general quadratic systems.

## 3. The other real finite equilibrium is a second center

For completeness, finite geometry can be checked without an orbit plot.
The origin-side cubic \(W=-T\) is positive on \(x>-1\) for \(m>18\)
in the present c-box, by the positive-quadratic bound in
`notes_lienard.md`. The organizer satisfies this condition.

Write a remote equilibrium as \(x=-s\), \(s>1\), and
\(y=s^2/(s-1)\). Its equation is \(m=\mathcal F(s,c)\), where

\[
 \mathcal F(s,c)=
 \frac{(61/5-c)s^3-(111/5)s^2+10s}{(s-1)^2}.
\]

Putting \(w=s-1>0\), differentiation gives

\[
 w^3\mathcal F_s=(61/5-c)w^3+(3c-11/5)w+2c>0.          \tag{U4}
\]

For \(w\le1\), the possibly negative linear term is at least
\(-7/10\), while \(2c\ge1\). For \(w\ge1\), the cubic and linear
terms together are at least \(10w\). Also \(\mathcal F\to-\infty\)
as \(s\to1+\), and \(\mathcal F\to\infty\) as \(s\to\infty\).
There is exactly one simple remote real equilibrium. At it,

\[
 \det DF=s(s-1)\mathcal F_s>0.
\]

The reversal maps equilibria to equilibria and fixes the origin. Uniqueness
forces the remote equilibrium onto its fixed line. Hence at (U1),

\[
 x_*=-\frac1{1-\sigma}<-1,\qquad
 y_*=\frac1{\sigma(1-\sigma)}.                         \tag{U5}
\]

Its trace is zero, since
\(\operatorname{tr}DF=(21/5)x+(1+2c)y\) and
\((1+2c_*)/\sigma=21/5\). Positive determinant and the reversal prove
that it too is a center. The two real finite equilibria are simple; the
remaining algebraic pair is nonreal. Small parameter changes preserve
that finite-equilibrium count.

## 4. At most one small origin cycle on the beta=0 slice

Let \(D(r;K,J)=R(r;K,J)-r\) be the positive-section return displacement
near the origin. The parameters \((K,J)\) are valid analytic local
coordinates near (U1), because \(J'(c_*)\ne0\) and
\(\partial K/\partial m=e/5>0\).

The cubic and quintic coefficients previously proved are

\[
 D(r;K,J)=\frac{\pi K}{4m^{3/2}}r^3+O(r^4),
\]

\[
 D(r;0,J)=\frac{\pi J}{150e\,m_0^{3/2}}r^5+O(r^6).
\]

The exact center gives \(D(r;0,0)\equiv0\). Analytic division therefore
gives an exact local factorization

\[
 D(r;K,J)=r^3[K\,a(r,K,J)+Jr^2b(r,J)],                \tag{U6}
\]

where \(a,b\) are positive in a sufficiently small common neighborhood and

\[
 a(0,0,0)=\frac\pi{4m_*^{3/2}},\qquad
 b(0,0)=\frac\pi{150e_*m_*^{3/2}}.
\]

For clarity, obtain the K term by subtracting \(D(r;0,J)\), which
divides by K; the remainder at K=0 divides by \(Jr^5\). The common
small-return domain follows from the uniformly rotating linear part.

For fixed nearby parameters, \(r^2b/a\) is strictly increasing for
small positive r. Thus (U6) has at most one small nonzero root. Any such
root requires \(KJ<0\), and its return derivative minus one has the sign
of J. It is stable for \(J<0\), unstable for \(J>0\), and is never a
multiplier-one fold cycle. At \(K=J=0\), the center's continuum of periodic
orbits consists of nonisolated orbits, not limit cycles.

Consequently this organizer cannot supply a local S/U pair on the
\(\beta=0\) slice. If an old outer S cycle persists, a separate pair
would still require a finite-amplitude mechanism not established here.

## 5. The stable-remote cone blocks an origin cycle shrinking to zero

The remote trace changes sign at

\[
 m_H(c)=\frac{21(1000c^2+1021c+481)}
 {50(2c+1)^2(8-5c)}.
\]

In the simple-remote regime it is negative for \(m>m_H\). This follows
also by evaluating the trace as
\(s[(1+2c)/(s-1)-(16-10c)/5]\) and using (U4).
The exact relation to the local focus coordinates is

\[
 m-m_H=
 \frac{250(2c+1)^2(8-5c)K+441J(c)}
 {50e(2c+1)^2(8-5c)}.                                 \tag{U7}
\]

The denominator is positive here. At the organizer both traces vanish;
on nearby \(K>0\) points satisfying the positive-numerator inequality,
the remote point is an attracting focus, since its determinant is
positive and its discriminant remains negative by continuity.

If \(K>0,J\ge0\), (U6) has no small positive root. If \(K>0,J<0\),
the stable-remote condition implies

\[
 \frac K{|J|}>
 \frac{441}{250(2c+1)^2(8-5c)}>0.                      \tag{U8}
\]

A root tending to zero in (U6) would instead require

\[
 r^2=\frac K{|J|}\frac ab,
 \qquad \frac ab\longrightarrow\frac{75e_*}{2}>0.
\]

This contradicts (U8). Uniform positive bounds on \(a/b\) make this a
common-neighborhood exclusion: no nonzero origin cycle can shrink into
the center while remaining in the \(K>0\), attracting-remote-focus cone.
No explicit amplitude radius is asserted. In particular the large finite
cycles being tracked are not accounted for by a higher-order local Hopf
unfolding at (U1).

## 6. Bounded recommendation and limits

Do not use (U1) as an order-three-focus seed or as a proposed source of
the missing local pair. The mathematically useful nearby control is the
remote trace gate (U7), together with the no-small-pair result, before
spending any further return evaluations.

For one exact rational illustration beyond the infinity threshold,
\(c=97/100\) gives \(J=-30429/10000\) and
\(K_H=69/350\). Thus continuing the old \(K=1/64\) path to this shape
fails the attracting-remote-focus gate. The rational parameter checkpoint
\(K=1/4\), \(\alpha=-21125/567\) clears that trace threshold by
\(37/700\). It is not a known return-root seed: the focus discriminant,
remote cycle, prescribed section cap, shared itinerary, and persistence
from an existing tracked branch remain to be checked. This observation
does not authorize an unseeded search or bypass the exhausted checkpoint's
limits. If a new bounded continuation is authorized, its first obligation
is transport of already known cycles across the infinity-stratum change
with all original gates retained; the algebraic checkpoint alone supplies
no cycle.

There is no global bound on origin cycles here, no exclusion of a separate
finite-amplitude S/U pair, and no five-cycle certificate. The exact organizer
and its local unfolding explain why a proposed local collapse or
order-three-focus completion is unavailable in this fixed KKL slice.
