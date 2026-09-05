# Exact Liénard and energy restrictions for the KKL precursor

Independent analytic lane, 2026-09-04. The field is

\[
 \dot x=(1+x)y+x^2,\qquad
 \dot y=-10x^2+\frac{11}{5}xy+cy^2+\alpha x,
                                                        \tag{L1}
\]

with \(\beta=0\). This note uses the inherited experimental box
\(1/2\le c\le3/2\), \(10\le m=-\alpha\le200\), and

\[
 K=m(11c/5-1)-42>0.
\]

All results below are direct analytic identities or elementary sign proofs.
No shooting, parameter scan, numerical replay, or CPU algebra experiment was
performed. There is no cycle-count or fold-existence theorem in this note.

## 1. Removing the quadratic velocity term

On the origin half-plane put \(u=1+x>0\) and \(v=\dot x\). Eliminating
\(y=(v-x^2)/u\) from (L1) gives

\[
 \ddot x=\frac{c+1}{u}\dot x^2+B(x)\dot x-\frac{xW(x)}u,
                                                        \tag{L2}
\]

where

\[
 B(x)=\frac{x[21+(16-10c)x]}{5u},
\]

\[
 W(x)=m+(2m+10)x+(m+111/5)x^2+(61/5-c)x^3.
                                                        \tag{L3}
\]

Thus \(W=-T\), with \(T\) the inherited nonzero-equilibrium cubic.
Use the increasing coordinate

\[
 z=\phi(x)=\frac{1-u^{-c}}c,\qquad \phi'(x)=u^{-c-1}.
\]

This is a smooth diffeomorphism from \(x\in(-1,\infty)\) to
\(z\in(-\infty,1/c)\), preserving the original time. Equation (L2)
becomes the standard Liénard equation

\[
 \ddot z+f(z)\dot z+g(z)=0,
 \quad f(\phi(x))=-B(x),\quad
 g(\phi(x))=xW(x)u^{-c-2}.                              \tag{L4}
\]

One may equivalently put \(w=\dot z+F(z)\), \(F'=f,F(0)=0\), to obtain
\(\dot z=w-F(z),\dot w=-g(z)\). An explicit primitive is

\[
 F(\phi(x))=-\frac15\int_0^x
     s[21+(16-10c)s](1+s)^{-c-2}\,ds.
\]

The changes are nonsingular on every compact origin-side cycle, including
large cycles. They do not extend through \(x=-1\) or identify distinct
infinity itineraries.

## 2. Exact restoring-force and damping signs throughout the box

Set \(d=16-10c\) and \(h=61/5-c\). The damping can be written

\[
 B(x)=\frac d5x+(1+2c)\frac{x}{1+x},\qquad
 B_x=\frac d5+\frac{1+2c}{u^2}>0.                       \tag{L5}
\]

Also \(21+dx>21-d=5+10c>0\) on \(x>-1\). Consequently
\(f\) is strictly decreasing, positive on the left of the origin and
negative on the right; its unique zero is the origin.

The restoring force has exactly one zero on this entire half-plane, without
using a numerical remote-equilibrium gate. To see this, rewrite (L3) as

\[
 W=h u^3+(m-72/5+3c)u^2+(11/5-3c)u+c.                 \tag{L6}
\]

Since \(K>0\) and \(c\le3/2\),
\(m>420/23>18\). Therefore \(h>0\),
\(m-72/5+3c>51/10\), \(11/5-3c\ge-23/10\), and \(c\ge1/2\).
For \(u>0\),

\[
 W>h u^3+\frac{51}{10}u^2-\frac{23}{10}u+\frac12>0.
\]

The quadratic has positive leading coefficient and discriminant
\(-491/100\). Hence \(g(z)z>0\) for \(z\ne0\), with
\(g'(0)=m>0\). This establishes only the origin-side equilibrium statement;
the separate remote focus and persistence gates remain necessary.

The common zero of \(f\) and \(g\) is material. A uniqueness or
nonexistence theorem that assumes distinct zeros of its damping and
restoring terms does not apply to this weak-focus slice. Nor does strict
monotonicity of \(f\) by itself supply a cycle-count theorem.

## 3. Energy identity and the infinity boundary

Define

\[
 V(\phi(x))=\int_0^x sW(s)(1+s)^{-2c-3}\,ds,
 \qquad E=\frac12\dot z^2+V(z).
\]

Then \(V(0)=0\), \(V>0\) off the origin, \(V'=g\), and

\[
 \dot E=-f(z)\dot z^2=B(x)u^{-2c-2}\dot x^2.          \tag{L7}
\]

Energy is gained during motion on the right and lost on the left. Every
periodic orbit satisfies the necessary balance

\[
 \oint B(x)u^{-2c-2}\dot x^2\,dt=0.                   \tag{L8}
\]

This does not determine the balance, because the orbit itself determines
the positive weight.

There is also an exact displacement interpretation on the actual positive
\(y=0\) section. Its energy is

\[
 \mathscr E(r)=\frac{r^4}{2(1+r)^{2c+2}}+V(\phi(r)),
 \qquad r>0.
\]

It is strictly increasing, since

\[
 \mathscr E'(r)=r(1+r)^{-2c-3}
    [W(r)+2r^2+(1-c)r^3]>0.                            \tag{L8a}
\]

The bracket has positive coefficients: its cubic coefficient is
\(66/5-2c>0\), its quadratic coefficient is \(m+121/5>0\), and the
remaining coefficients are \(2m+10\) and \(m\).
If a selected full return exists with \(r>0,R(r)>0\), then

\[
 \mathscr E(R(r))-\mathscr E(r)
 =\int_0^{T(r)} B(x)u^{-2c-2}\dot x^2\,dt.             \tag{L8b}
\]

Thus the integral in (L8b) has exactly the sign of \(R(r)-r\), including
away from a fixed point. If this energy-gain function is called
\(\mathscr G(r)\), then at a periodic point
\(\mathscr G_r=\mathscr E'(r)(R_r-1)\); at a fold additionally
\(\mathscr G_{rr}=\mathscr E'(r)R_{rr}\).
These are identities on a differentiable full-return itinerary, not bounds
on the sign of the integral or a replacement for event sensitivities.

The potential diverges at the left boundary:

\[
 V(\phi(x))\sim\frac{c}{2c+2}u^{-2c-2},\qquad u\to0+.
\]

At the right boundary its behavior is

\[
 V(\phi(x))\sim\frac{h}{2-2c}x^{2-2c}\quad(c<1),
 \qquad V(\phi(x))\sim h\log x\quad(c=1).
\]

For \(c>1\), \(V(1/c-)\) is finite. These assertions follow directly
from the leading terms of (L6); they concern the potential, not boundedness
of the actual trajectories. Negative damping can increase the energy.
The original compactified itinerary still has to be checked separately,
including the inherited \(c=241/250\) and \(c=1\) boundary strata.

## 4. Exact multiplier and fold moments

The quotient \(R=f/g\) extends analytically through the origin and is

\[
 R(\phi(x))=-\frac{[21+dx]u^{c+1}}{5W(x)}.
\]

Define the quartic polynomial

\[
 N(x)=\{du+(c+1)(21+dx)\}W(x)-u(21+dx)W'(x).           \tag{L9}
\]

Differentiation gives

\[
 R_z(\phi(x))=-\frac{u^{2c+1}N(x)}{5W(x)^2},\qquad
 N(0)=5K,\qquad R_z(0)=-\frac K{m^2}.                 \tag{L10}
\]

The leading coefficient of \(N\) is \(d(c-1)h\), and

\[
 N(-1)=c(c+1)(5+10c)>0.
\]

For any periodic orbit let \(\mathcal M>0\) denote its transverse
multiplier. There is an exact identity

\[
 \boxed{\quad
 \log\mathcal M=
 \oint\frac{N(x)}{5uW(x)^2}\dot x^2\,dt.
 \quad}                                               \tag{L11}
\]

Proof: in (L4), write \(v_z=\dot z\). On a closed orbit,

\[
 \oint f\,dt
 =-\oint R\,dv_z-\oint Rf\,dz
 =\oint R_z v_z^2\,dt.
\]

The integral of \(Rf\,dz\) is zero because it is a function of \(z\)
times \(dz\). The multiplier is \(\exp(-\oint f\,dt)\), so (L10)
and \(v_z=u^{-c-1}\dot x\) prove (L11). Equivalently, in original
coordinates

\[
 \operatorname{div}(L1)=B(x)+(1+2c)\frac{d}{dt}\log u,
\]

which checks that the integral of the original divergence equals
\(-\oint f\,dt\).

Every stable cycle must sample a region where \(N<0\). A multiplier-one
origin cycle must sample both signs: it crosses \(x=0\), where \(N=5K>0\),
with nonzero velocity, so its positive contribution is strict. If \(N\)
is nonnegative on a proposed cycle's entire x-range, that cycle cannot be
stable or a fold. This is a genuine necessary amplitude/itinerary restriction,
not an existence result.

For comparison with the energy condition put

\[
 d\nu=\frac{\dot x^2}{5uW(x)^2}\,dt,\qquad
 A(x)=x(21+dx)W(x)^2u^{-2c-2}.
\]

The necessary fold moments are

\[
 \oint A\,d\nu=0,\qquad \oint N\,d\nu=0.              \tag{L12}
\]

Here \(d\nu\ge0\), and \(A\) has the sign of \(x\). These equalities
apply only at an actual closed orbit. They do not replace return existence,
\(D=0\), or the ordinary-fold conditions \(D_{rr}\ne0\) and a nonzero
transverse parameter derivative.

## 5. The starting shape defeats a naive monotonic-ratio test

At the exact starting shape \((c,m)=(7/10,80)\),

\[
 W(x)=80+170x+\frac{511}{5}x^2+\frac{23}{2}x^3
      =\frac{115u^3+677u^2+u+7}{10},
\]

and direct expansion of (L9) gives

\[
 N(x)=6+\frac{753}{5}x+\frac{7821}{25}x^2
          +\frac{12291}{100}x^3-\frac{621}{20}x^4.      \tag{L13}
\]

The following exact values suffice to classify its roots, without locating
them numerically:

\[
 N(-1)=\frac{357}{25}>0,\quad
 N(-1/2)=-\frac{13431}{1600}<0,\quad N(0)=6>0,
 \quad N(\pm\infty)=-\infty.
\]

The intermediate value theorem gives four distinct real roots, one in
each of

\[
 (-\infty,-1),\quad(-1,-1/2),\quad(-1/2,0),\quad(0,\infty).
\]

Since \(N\) has degree four, these are all its roots and each is simple.
Call the three origin-side roots \(a_1<a_2<0<a_3\). Then

\[
 N<0\text{ on }(a_1,a_2)\cup(a_3,\infty),\qquad
 N>0\text{ on }(-1,a_1)\cup(a_2,a_3).
\]

Thus \(R_z\) changes sign three times in the allowed half-plane. Every
stable or fold cycle at this shape must either enter the left band
\((a_1,a_2)\) or reach beyond \(a_3\) on the right. The available stable
control does not imply which alternative occurs; that is an orbit-level
question for the single numerical lane.

There is also a precise limitation on the simplest use of (L12). No
constant linear combination \(N+\lambda A\) can have one strict sign
throughout \((-1,\infty)\): positivity on the left negative band requires
\(\lambda<0\), while positivity sufficiently far right requires
\(\lambda>0\); negativity is impossible at zero, where \(A=0,N=6\).
Therefore these two moments alone admit no global constant-combination
sign certificate at the starting shape. This does not preclude more refined
orbit-dependent inequalities or another Dulac function.

More generally, for every \(K>0\) point of the box with \(c<1\),
\(N(0)>0\) and its leading coefficient is negative. Hence at least one
positive root of \(N\) is unavoidable. A globally monotone \(f/g\)
argument cannot apply there without an additional bound restricting the
cycle's range.

## 6. Consequence for the bounded construction

The Liénard chart gives a reproducible independent description of the
origin-side flow and two exact diagnostics for actual cycles: energy
balance (L8) and the multiplier integral (L11). Interval bounds showing
\(N\ge0\) throughout a candidate cycle's x-range would rule out its
stability or fold status. Conversely, reaching a negative band is only a
necessary condition.

Neither shape parameter is automatically a rotated-vector-field control.
In the original coordinates,

\[
 \det(F,\partial_\alpha F)=x\dot x,\qquad
 \det(F,\partial_c F)=y^2\dot x,
\]

which have no fixed sign around a general origin cycle. Therefore one may
not infer monotonic return displacement or a one-way fold crossing from
increasing \(\alpha\) or \(c\) alone.

The missing S/U/S precursor, coexistence with the remote cycle, and ordinary
fold transversality remain unproved. The new analytic restrictions support
the inherited bounded continuation; they neither replace it nor justify a
whole-box exclusion after an unsuccessful seeded search.

## 7. The second weak-focus sign on K=0

This continuation of the analytic lane was requested after the charged
pilot, while the parent alone performed orbit computations. The derivation
below used hand algebra only. It supplies a local obstruction to collapse,
not a finite-amplitude fold prediction.

On \(K=0\), write

\[
 e=11c-5,\qquad m_0=210/e,
 \qquad J(c)=305+634c-11c^2-1000c^3,
\]

\[
 \Delta(c)=\frac{4J(c)}{5e}.                            \tag{L14}
\]

Here \(e>0\) on the box. The sign of the second weak-focus term is the
sign of \(\Delta\). More precisely, on the actual positive \(y=0\)
section, at \(m=m_0\),

\[
 R(r)-r=\frac{\pi\Delta(c)}{120m_0^{3/2}}r^5+O(r^6).
                                                        \tag{L15}
\]

This formula presupposes \(\Delta\ne0\) when calling the origin an
order-two weak focus. If \(J=0\), the displayed coefficient vanishes and
this note does not assign the next focus order.

### Derivation of the coefficient

The harmonic-energy coordinate is
\(s=\operatorname{sgn}(z)\sqrt{2V(z)}\). With the positive time change
\(d\tau/dt=g(z)/s\), equations (L4) become

\[
 \frac{ds}{d\tau}=w,\qquad
 \frac{dw}{d\tau}=-s-h(s)w,\qquad h(s)=sR(z(s)).        \tag{L16}
\]

The quotient \(g/s\) extends positively through zero. Expanding the
potential as a function of the original coordinate gives

\[
 s=\sqrt m\,x[1+a x+O(x^2)],\qquad
 a=\frac{10/m-(2c+1)}3.
\]

Let \(N_1,N_2\) denote the coefficients of \(x,x^2\) in (L9).
Direct multiplication gives, for all these parameters,

\[
 N_1=m(-10c^2+38c-10)+210c-4662/5,
\]

\[
 N_2=m(-20c^2+43c-5)-100c^2+(4056/5)c-1430.            \tag{L17}
\]

At \(K=0\), \(R_x(0)=0\). If
\(R(x)=R_0+R_2x^2+R_3x^3+O(x^4)\), differentiating its exact quotient
gives

\[
 R_2=-\frac{N_1}{10m^2},\qquad
 R_3=-\frac{N_2+(c-4-20/m)N_1}{15m^2}.
\]

Consequently the coefficient of \(s^3\) in \(R(z(s))\) is

\[
 \frac{R_3-2aR_2}{m^{3/2}}
 =\frac{(d/7)N_1-N_2}{15m^{7/2}}
 =-\frac{\Delta(c)}{15m^{7/2}},                        \tag{L18}
\]

where \(d=16-10c\), \(30/m+3-3c=d/7\) on \(K=0\), and substitution
of (L17) yields \(N_2-dN_1/7=\Delta(c)\).

The coefficient of \(s^2\) in the damping \(h(s)\) is zero on \(K=0\);
its first possible nonzero even coefficient is the \(s^4\) coefficient
in (L18). To see its dynamical significance, retain all odd terms of
\(h\) as a reversible center. An even term \(a_4s^4\) first contributes
to the radial return at order five. Its contribution over the linear
period is

\[
 -a_4\rho^5\int_0^{2\pi}\cos^4\theta\sin^2\theta\,d\theta
 =-\frac\pi8 a_4\rho^5.
\]

The lower odd damping terms change this contribution only at higher
order; by reversibility their own return is the identity. The actual
section has harmonic radius \(\rho=\sqrt m\,r+O(r^2)\), so converting
the leading coefficient multiplies it by \(m^2\). This proves (L15).
As a normalization check, the analogous quadratic-even-damping calculation
away from \(K=0\) gives the inherited cubic coefficient
\(\pi K/(4m^{3/2})\).

### The shape threshold and the current c=33/40 slice

On \([1/2,3/2]\), \(J'=634-22c-3000c^2<0\). Moreover

\[
 J(19/20)=15999/400>0,\qquad J(39/40)=-11333/800<0.
\]

Thus there is exactly one zero \(c_\dagger\) of \(J\) in the box,
and it lies in \((19/20,39/40)\). This is an exact local-focus threshold,
not a located finite-cycle fold curve.

At the current shape \(c=33/40\),

\[
 m_0=8400/163,\quad J(33/40)=103619/400>0,\quad
 \Delta(33/40)=207238/4075>0.                          \tag{L19}
\]

The order-two weak focus on \(K=0\) is therefore repelling. Analyticity
of the small return map and the positive cubic and quintic coefficients
give numbers \(\epsilon,\rho_0>0\) such that

\[
 R(r;K)-r>0\qquad(0\le K<\epsilon,\ 0<r<\rho_0)
                                                        \tag{L20}
\]

on this fixed c-slice. One way to justify uniformity is to write
\(D(r,K)=D(r,0)+K\int_0^1D_K(r,tK)dt\):
\(D(r,0)/r^5\) and \(D_K(r,K)/r^3\) are both positive near zero.
In particular, a stable cycle cannot shrink into the origin as \(K\)
decreases to zero from the positive side here.

For orientation only, if \(c\ne c_\dagger\) is fixed, the small
nonzero cycle on the opposite-sign side \(K\Delta<0\) satisfies
\(r^2\sim-30K/\Delta\). It is unstable when \(\Delta>0\), and stable
when \(\Delta<0\). This is a local asymptotic consequence of (L15),
not a certificate of any finite-amplitude cycle or coexistence event.

## 8. A common-amplitude obstruction from both exact moments

The same \(\Delta\) has a direct energy/multiplier interpretation. Expand
the function \(A\) in (L12) at zero:

\[
 A(x)=21m^2x+21m^2\left(\frac d{21}+2-2c+\frac{20}m\right)x^2
       +O(x^3).
\]

Define the analytic function

\[
 \widehat N(x)=N(x)-\frac{N_1}{21m^2}A(x).
\]

Its linear term vanishes identically. On \(K=0\), the coefficient in
parentheses above is \(d/7\), and hence

\[
 \widehat N(x)=\Delta(c)x^2+O(x^3).                     \tag{L21}
\]

At \(c=33/40\), for all sufficiently small \(K\ge0\), continuity gives

\[
 \widehat N(x)=5K+\Delta_Kx^2+O(x^3),\qquad
 \Delta_K>\tfrac12\Delta(33/40)>0,
\]

with a remainder uniform in these parameters. Therefore some common
\(\delta>0\) makes \(\widehat N>0\) on
\((-\delta,\delta)\setminus\{0\}\) for the entire small nonnegative-K
interval. For an actual cycle whose entire x-range lies there, (L12)
eliminates the A moment and gives

\[
 \log\mathcal M=\oint\widehat N\,d\nu>0.
\]

This excludes both a stable cycle and a multiplier-one cycle confined to
that common amplitude range. A missing stable/unstable pair produced by
a fold therefore cannot appear entirely in a shrinking neighborhood of
the origin on this approach to \(K=0\). The radius \(\delta\) is proved
to exist, but has not been quantified; no numerical amplitude gate is
claimed.

This local statement does not contradict the global failure of a
constant-combination sign certificate in §5. It also leaves a finite fold,
large-cycle continuation, escape, and changed return itinerary unresolved.
