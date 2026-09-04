# Q4 Strike 4: a global proof of Theorem N

2026-09-04. Based on the audited canon at `084fae3`.

**Theorem N is proved.** For every finite `kappa>1` and every coefficient
point in the strict lobe region, the first primitive root satisfies

\[
\boxed{\Phi_a(\tau_1)<0,\qquad a=1-1/\kappa.}
\]

The proof is analytic and global. Three strict comparisons replace the
proposed interval cover and both endpoint expansions. In particular, no
restriction on endpoint rates or anchor ratios remains. The proof uses
the inherited fixed-ratio corner limit only as a comparison destination;
it does not interchange two singular limits.

The immediate global original-integral consequence is **at most four
distinct interior zeros**. On the strict lobe region the stronger bound is
**at most three distinct interior zeros**. The proposed deduction of a
global three-zero bound has an additional outside-lobe gap, described in
section 7. This document does not claim to settle that conjecture.

## 1. Inherited objects and exactly what is new

Use the audited universal chart

\[
F={}_2F_1(1/6,5/6;1;t),\quad
M=1-6(1-t)F'/F,\quad w=tF>0,
\]
\[
K_0(t)=\int_0^t w,\quad K_1(t)=\int_0^t uw(u)\,du,
\quad K_2(t)=\int_0^t M(u)w(u)\,du,
\quad K_3(t)=\int_0^t uM(u)w(u)\,du.
\]
\[
H=(A-1)K_0+BK_1-\eta K_2+K_3,
\qquad
Y_0=\frac{3(1326A+864B-2431\eta-102)}{1361360}.
\]

The strict lobe region is globally parametrized by the three simple roots
`0<y1<y2<y3<1` of `H`, with signs `+,-,+,-`. Write `y1=tau1`.
The inherited positive Stieltjes representation implies `M'>0` and
`M''>0`. The inherited Green identities are

\[
\Phi_a(t)=Y_0+\int_0^t W_a(u)H(u)\,du,
\quad W_a=\mathcal R_a\Omega_a>0,
\]
\[
p_a=\sqrt{\frac{1-t}{1-at}},\quad
\mathcal R_a(t)=\int_0^t\frac{du}{p_a(u)y_a(u)^2},
\quad
\Omega_a=\frac{y_a}{1152t^2(1-at)^{3/2}(1-t)^{3/2}}.
\]

Here `y_a>0` is the normalized homogeneous solution, `y_a(0)=1`,
vanishing at the loop, and

\[
(1-at)(1-t)y_a''-\frac{1-a}{2}y_a'+\frac{5a}{36}y_a=0.
\]

The audited fixed-ratio corner and its positive primitive are

\[
(A_*,B_*,\eta_*)=(94/77,-17/77,1),\qquad Y_{0,*}=-3/1232,
\]
\[
H_*(t)=\frac{6t(1-t)^2F(t)}{77}(6M(t)-1)>0.
\]

At `a=1`, interpreted as the limiting differential equation,

\[
y_1(t)=(1-t)^{5/6},\quad
W_1(t)=\frac{3}{2304t^2}
\left[(1-t)^{-17/6}-(1-t)^{-13/6}\right],
\]
\[
\boxed{\int_0^1 W_1H_*\,dt=\frac3{1232}=-Y_{0,*}.} \tag{N0}
\]

The notation `y_1(t)` in this display denotes the homogeneous solution at
`a=1`, not an anchor. The new ingredients below are monotonicity in all
three anchors and the global comparison `W_a<W_1`.

## 2. A convex moment curve controls coefficient variations

Put `x=K1/K0` and `m=K2/K0`. They are weighted averages of `u` and
`M(u)` on `(0,t)`, so

\[
x'=(w/K_0)(t-x)>0,\qquad m'=(w/K_0)(M-m)>0.
\]

The slope of the moment curve is

\[
S(t)=\frac{dm}{dx}=\frac{M(t)-m(t)}{t-x(t)}.
\]

It is a positive weighted average, with weights `w(u)(t-u)`, of the
secants `[M(t)-M(u)]/(t-u)`. Strict convexity of `M` gives `S<M'(t)`;
differentiating the displayed quotient cancels the weight terms and gives

\[
S'(t)=\frac{M'(t)-S(t)}{t-x(t)}>0.
\]

Thus `m` is strictly convex as a function of `x`. The exact endpoints are

\[
x(0)=0,\quad m(0)=1/6,\quad
x(1)=6289/9061,\quad m(1)=11/41,
\]
\[
\boxed{S(1-)=\frac{1-11/41}{1-6289/9061}=1105/462.} \tag{N2}
\]

Consider a nonzero variation

\[
D=\alpha K_0+\beta K_1+\gamma K_2
=K_0(\alpha+\beta x+\gamma m)
\]

with two distinct interior roots and `gamma>0`. Strict convexity shows
that its bracket is positive outside those roots and negative between
them. Consequently

\[
\alpha> -\gamma/6,\qquad
\beta>-(1105/462)\gamma.
\]

The induced variation of the center value is exactly

\[
\mathscr Y(D)=\frac9{3080}
\left(\alpha+\frac{144}{221}\beta+\frac{11}{6}\gamma\right).
\]

It follows that

\[
\boxed{\mathscr Y(D)>
\frac9{3080}\gamma
\left(\frac53-\frac{144}{221}\frac{1105}{462}\right)
=\frac9{3080}\gamma\frac{25}{231}>0.} \tag{N3}
\]

For negative `gamma`, reverse every sign. This also proves that the
center functional and a two-root variation have the same initial sign.

## 3. All three anchors increase both the first lobe and the center value

Let `ell_j` be the cardinal interpolant in `span{K0,K1,K2}` satisfying
`ell_j(yi)=delta_ij`. The interpolation matrix has rows proportional to
`(1,x(yi),m(yi))`. Its determinant is positive by strict convexity, and
its `1,x` minors are positive. Hence the coefficient of `K2` in `ell_j`
has sign `(-1)^(j-1)`. The interpolant has exactly the two other anchors
as simple roots. Both its sign before the first anchor and, by (N3), its
center functional have sign `(-1)^(j-1)`.

Differentiate the anchor equations. The `K3` coefficient stays fixed, so

\[
\partial_{y_j}H(t)=-H'(y_j)\ell_j(t).
\]

Since `sign H'(yj)=(-1)^j`, this proves

\[
\boxed{\partial_{y_j}H(t)>0\ (0<t<y_1),\qquad
\partial_{y_j}Y_0>0,\qquad j=1,2,3.} \tag{N4}
\]

Starting from any strict triple, increase its anchors along

\[
y_j(\theta)=1-(1-\theta)(1-y_j),\qquad0\le\theta<1.
\]

The anchors stay strictly ordered and retain fixed strict ratios of
their distances from one. The inherited corner limit therefore applies
to this path. Taking the limit in the strict monotonicity statements gives

\[
\boxed{Y_0<Y_{0,*},\qquad 0<H(t)<H_*(t)\quad(0<t<y_1).} \tag{N5}
\]

Strictness survives the limit: any positive intermediate path parameter
already gives a strictly larger value. At `t=y1`, `H(y1)=0<Hstar(y1)`.
The comparison holds for each original triple, including triples whose
own endpoint approaches have degenerating ratios or only one late anchor.

## 4. A positive residual bounds the Green kernel for every lift

The second homogeneous solution `z_a=y_a Rcal_a` satisfies
`z_a(0)=0,z_a'(0)=1`. Define

\[
v_a(t)=\frac{z_a(t)}{(1-at)^{3/2}},\qquad
W_a(t)=\frac{v_a(t)}{1152t^2(1-t)^{3/2}}.
\]

Exact substitution in the homogeneous equation gives

\[
L_a v_a=0,\qquad v_a(0)=0,\quad v_a'(0)=1,
\]
\[
L_a=(1-at)(1-t)D_t^2-
\frac{1+5a-6at}{2}D_t+\frac{8a}{9}.
\]

At `a=1`, writing `d=1-t`,

\[
v_1=\frac32(d^{-4/3}-d^{-2/3}),
\]
\[
\boxed{L_av_1=(1-a)
\left[\frac{11}{3}d^{-7/3}-\frac76d^{-5/3}\right]
=\frac{(1-a)(22-7d^{2/3})}{6d^{7/3}}>0.} \tag{N6}
\]

For completeness, the causal Green function has a strict positive sign.
The function `f_a=y_a/(1-at)^(3/2)` is a positive solution of `L_a f=0`.
After division by the positive leading coefficient, write the equation as
`v''+b(t)v'+c(t)v=g`. Substitution `v=f_a u` gives

\[
(\rho u')'=\rho g/f_a,\qquad
\rho=f_a^2\exp\!\left(\int b\right)>0.
\]

With zero initial value and derivative, positive forcing therefore gives
a strictly positive solution at every positive time. Apply this on any
compact interval in `[0,1)` to `v1-va`, using (N6) and the identical
initial conditions. We obtain

\[
\boxed{0<W_a(t)<W_1(t)\quad(0\le a<1,\ 0<t<1).} \tag{N7}
\]

There is no assertion of monotonicity between arbitrary finite values of
`a`, and no uniform singular expansion is needed.

The included endpoint `a=0` uses `y_0(t)=sqrt(1-t)` and
`p_0(t)=sqrt(1-t)`, so `z_0(t)=v_0(t)=2(1-sqrt(1-t))`.
These explicit positive homogeneous data justify the same comparison
there. This endpoint extension is optional for Theorem N, whose original
parameter domain is `0<a<1`.

## 5. Completion of Theorem N and all endpoint parts

Combine (N0), (N5), and (N7):

\[
\begin{aligned}
\Phi_a(y_1)
&=Y_0+\int_0^{y_1}W_aH\,dt\\
&<Y_{0,*}+\int_0^{y_1}W_1H_*\,dt\\
&=-\int_{y_1}^1W_1H_*\,dt<0.
\end{aligned} \tag{N8}
\]

All integrals converge. Near zero, `Hstar=O(t^2)` and `W1=O(1/t)`.
Near one their product is `O((1-t)^(-5/6) log(1/(1-t)))`, which is
integrable. This proves Theorem N on its full stated domain, and also at
`a=1` for every fixed strict anchor triple.

The finite-lift corner itself has the exact sign

\[
\Phi_*(1;a)=Y_{0,*}+\int_0^1W_aH_*\,dt<0\quad(a<1).
\]

Moreover the stronger comparison
`Phi_a(y1)<Phi_star(1;a)` follows directly from (N5) and the positive
tail with `W_a`. Dominated convergence using `W1 Hstar` proves continuity
of this corner functional in `a`, including its limit zero at `a=1`.
It gives a uniform negative margin when `a` stays below any fixed
`amax<1`, even as some or all anchors approach one.

For joint approaches to `a=1` and the loop, (N8) remains a strict
pointwise inequality at every admissible parameter tuple. The right side
may approach zero, but never changes sign. This treats every rate and
every degenerating anchor ratio without an asymptotic remainder estimate.

## 6. An explicit compact certificate

The inherited positive series give `F>=1` and
`M>=1/6+(25/432)t`. Since

\[
W_1H_*=
\frac{F(6M-1)}{9856t}
\left[(1-t)^{-5/6}-(1-t)^{-1/6}\right],
\]

(N8) implies, whenever `y1<=1-delta`,

\[
\boxed{\Phi_a(y_1)<-
\frac{25}{118272}\left(\delta^{1/6}
-\frac{\delta^{5/6}}5\right).} \tag{N9}
\]

In particular choose `delta=1/64`. Throughout the requested compact part
`a in [2593/21636,1-delta]`, `y3<=1-delta`, the stronger bound, independent
of `a`, is

\[
\boxed{\Phi_a(y_1)<-\frac{395}{3784704}<0.}
\]

This is an exact analytic certificate, stronger than a covering by
interval boxes. No numerical margin or sampled sign is used.

## 7. Original-integral consequences and the remaining distinction

The audited necessary condition (N1) says that five distinct original
zeros force strict lobe membership and `Phi_a(tau1)>0`. Theorem N
contradicts that condition. Therefore every nonzero original Q4 integral,
for every finite `kappa>1`, has at most **four distinct interior zeros**.

Inside the strict lobe region, the bound improves to **three distinct
interior zeros**. Here `Z(0)=Y0<0`, `Z(1-)>0`, `P(1-)>0`, and
`P'=-Omega H` has successive signs `-,+,-,+`. If `P0<=0`, `Z` initially
decreases and the remaining monotone pieces allow at most three sign
changes. If `P0>0` and there are at most two crossings of `P`, there are
at most three monotone pieces of `Z`. If there are four crossings
`p1<p2<p3<p4`, the first lies before `tau1` and

\[
Z(p_1)=\Phi_a(p_1)<\Phi_a(\tau_1)<0.
\]

Thus `Z` stays negative through its first maximum and the following
minimum; at most three sign changes remain. A zero of `P` at an extremum
is a tangency and does not create another monotone piece, so it cannot
increase this count. Since `Y=yZ` with `y>0`, `Y` also has at most three
sign changes.

Finally `X(0)=0` and `X'=Y/(1-at)^(3/2)`. If `X` has `n` distinct
interior zeros, each of the `n` consecutive intervals beginning at zero
and ending at those zeros contains a sign change of `X'`: an analytic,
nonzero function returning to the same value has an interior extremum
with a derivative sign change. These intervals are disjoint. Thus the
number of distinct zeros of `X`, and hence of `I`, is bounded by the
number of sign changes of `Y`. This handles nonsimple original zeros
as well as simple ones within the lobe region.

The handoff's appended global three-zero claim applies this lobe argument
without an outside-lobe step. The inherited condition places **five**
distinct original zeros in the lobe region; it does not place **four**
there. Outside the lobe region, the inherited estimate
`Z(I)<=Z(H)+2` can still allow four when `H` has two interior zeros.
Theorem N has no hypothesis covering those coefficient directions.
Consequently a global three-zero conclusion, or closure of the conjecture,
does not follow from Theorem N and the stated canon alone.

The five-zero construction route is closed. Whether the remaining
outside-lobe directions can produce four distinct interior zeros is not
resolved by this strike. No four-zero example is asserted.

## 8. Handoff transcription and verification record

The supplied affine formula has a missing `-K0` in its constant term.
The correct constant is

\[
c_0=-306/1361360+\int_0^{y_1}W_a(K_3-K_0)\,dt.
\]

The other three affine coefficients in the handoff agree with the exact
objects. Existing reconstruction scripts use the full `H` correctly.
This is a correction to the new proposed affine decomposition, not to any
of the three audited strikes.

Independent derivations are recorded in
[q4/notes_N_compact.md](q4/notes_N_compact.md),
[q4/notes_N_loop.md](q4/notes_N_loop.md), and
[q4/notes_N_double.md](q4/notes_N_double.md).
The new exact checks verify the rational constants and transformed
operator residual. The analytic convexity, Green positivity, and strict
comparisons above supply the proof; finite numerical tests do not.
No coefficient sweeps, tuned shots, confluent shots, reverse tangencies,
or searches for `Y=Y'=0` were run. Claude's other lanes were untouched.
