# Strike 4: global anchor comparison and the finite-kappa loop boundary

2026-09-04. This note owns the loop-approach component of Theorem N. The
argument gives a stronger global comparison on the full strict lobe region.
Combined with the independently proved lift-kernel comparison, it proves
Theorem N without any parameter scan or numerical sign estimate.

The universal chart, global anchor coordinates, fixed-ratio corner limit,
center data, and beta identities are inherited from the audited canon.
All new inequalities in this note are analytic. Exact arithmetic checks in
`q4_N_loop_checks.py` verify the rational constants only; they are not a
substitute for the sign arguments.

## 1. Exact objects and two corrections to the proposed decomposition

Write `w(t)=tF(t)>0` on `(0,1)`, and use the inherited moments

\[
K_0(t)=\int_0^t w(u)\,du,\qquad
K_1(t)=\int_0^t u w(u)\,du,\qquad
K_2(t)=\int_0^t M(u)w(u)\,du.
\]

The full normalized primitive is

\[
H(t)=(A-1)K_0(t)+BK_1(t)-\eta K_2(t)+K_3(t).
\]

Its three simple roots are the global anchor coordinates
`0<y1<y2<y3<1`; its signs are `+,-,+,-`. The fixed center datum is

\[
Y_0=\frac{3(1326A+864B-2431\eta-102)}{1361360}.
\]

First, the constant term in an affine expansion of
`Phi=Y0+integral Rcal*Omega*H` includes `K3-K0`, not `K3` alone.
Second, sending only the last anchor `y3` to one does not force the other
anchors to one and does not justify replacement by the corner primitive.
The comparison below treats such partial endpoint strata directly.

Let `W_a(t)=Rcal_a(t)*Omega_a(t)>0`, so

\[
\Phi_a(y_1)=Y_0+\int_0^{y_1}W_a(t)H(t)\,dt.
\]

The corner data, inherited for every fixed strict anchor-ratio collapse, are

\[
(A_*,B_*,\eta_*)=(94/77,-17/77,1),\qquad Y_{0,*}=-3/1232,
\]
\[
H_*(t)=\frac{6t(1-t)^2F(t)}{77}(6M(t)-1)>0\quad(0<t<1).
\]

## 2. A strictly convex moment curve

Normalize the three-dimensional coefficient-variation space by setting

\[
x(t)=K_1(t)/K_0(t),\qquad m(t)=K_2(t)/K_0(t).
\]

These are weighted averages of `u` and `M(u)` over `(0,t)`, respectively.
Thus `x(t)<t`, `m(t)<M(t)`, and

\[
x'(t)=\frac{w(t)}{K_0(t)}[t-x(t)]>0,\qquad
m'(t)=\frac{w(t)}{K_0(t)}[M(t)-m(t)]>0.
\]

To distinguish derivatives with respect to the two coordinates, put

\[
S(t):=\frac{dm}{dx}=\frac{M(t)-m(t)}{t-x(t)}.
\]

Strict convexity of `M` follows from its inherited positive Stieltjes
representation. Moreover

\[
S(t)=
\frac{\int_0^t w(u)(t-u)
  \frac{M(t)-M(u)}{t-u}\,du}
 {\int_0^t w(u)(t-u)\,du}<M'(t).
\]

This is a positive weighted average of strict secant slopes of `M`.
Direct differentiation then gives

\[
S'(t)=\frac{M'(t)-S(t)}{t-x(t)}>0. \tag{N-L1}
\]

Consequently `m`, regarded as a function of the strictly increasing
coordinate `x`, is strictly increasing and strictly convex. In particular
`{1,x,m}` is an extended complete Chebyshev family; the ordinary
interpolation determinants used below are strictly positive in increasing
node order. This conclusion also follows from anchored Rolle applied to
the inherited family `{1,t,M}`, but (N-L1) proves the needed sign directly.

At the center, `x(0+)=0` and `m(0+)=M(0)=1/6`. At the loop, the exact
moment identities give

\[
x(1)=\frac{6289}{9061},\qquad m(1)=\frac{11}{41},
\]
\[
\boxed{S(1-)=\frac{1-m(1)}{1-x(1)}=\frac{1105}{462}.} \tag{N-L2}
\]

All these endpoint limits are finite. Thus every strict secant slope of
`m(x)` between interior nodes is less than `1105/462`.

## 3. The center functional is positive on each positive two-root direction

Consider a nonzero variation primitive

\[
D(t)=\alpha K_0(t)+\beta K_1(t)+\gamma K_2(t)
=K_0(t)[\alpha+\beta x(t)+\gamma m(t)]
\]

with two distinct interior roots and `gamma>0`. The bracket is a strictly
convex function of `x`. It is strictly positive before its first root and
after its second root, and strictly negative between them. Its center
value and root secant slope therefore imply

\[
\alpha+\gamma/6>0,\qquad
\beta=-\gamma\frac{m(x_2)-m(x_1)}{x_2-x_1}
>-\frac{1105}{462}\gamma. \tag{N-L3}
\]

Here `x1,x2` denote the two roots in the `x` coordinate; no assumption is
made on how close they are to either endpoint.

A coefficient variation `D` corresponds to
`delta A=alpha`, `delta B=beta`, `delta eta=-gamma`. The induced center
functional is exactly

\[
\mathscr Y(D):=\delta Y_0
=\frac9{3080}\left[\alpha+\frac{144}{221}\beta
                         +\frac{11}{6}\gamma\right].
\]

Using (N-L3) gives the strict sign with an exact rational margin:

\[
\boxed{
\mathscr Y(D)>
\frac9{3080}\gamma
\left[\frac53-\frac{144}{221}\frac{1105}{462}\right]
=\frac9{3080}\gamma\frac{25}{231}>0.} \tag{N-L4}
\]

For `gamma<0`, apply the same argument to `-D`. Therefore a two-root
variation primitive and its center functional have the same sign before
its first root, namely the sign of the coefficient `gamma`.

## 4. Increasing any anchor increases both required quantities

For fixed anchors `y1<y2<y3`, let `ell_j` be the cardinal element of
`span{K0,K1,K2}` satisfying `ell_j(yi)=delta_ij`. Existence and uniqueness
follow from the strictly convex moment curve in section 2.

Write `gamma_j` for the coefficient of `K2` in `ell_j`. Cofactor expansion
of the interpolation matrix, whose rows are
`K0(yi)*(1,x(yi),m(yi))`, gives

\[
\operatorname{sign}(\gamma_j)=(-1)^{j-1}. \tag{N-L5}
\]

Indeed its determinant is positive, and the two-column minors formed by
`1,x` are positive in increasing node order. Each `ell_j` has exactly the
two other anchors as its roots; hence its sign before the first anchor is
also `(-1)^(j-1)`. By (N-L4), the center functional
`mathscr Y(ell_j)` has that same sign.

Differentiate the three anchor equations `H(yi)=0` with respect to `yj`.
The coefficient of `K3` is fixed, so the differentiated primitive lies in
this three-dimensional variation space. Its exact formula is

\[
\frac{\partial H(t)}{\partial y_j}=-H'(y_j)\ell_j(t).
\]

The strict lobe signs give `sign H'(yj)=(-1)^j`. Combining with (N-L5)
and (N-L4) yields

\[
\boxed{\frac{\partial H(t)}{\partial y_j}>0\quad(0<t<y_1),}
\qquad
\boxed{\frac{\partial Y_0}{\partial y_j}>0,}
\quad j=1,2,3. \tag{N-L6}
\]

These are exact global monotonicity statements on the entire strict anchor
simplex, not local statements near the endpoint corner.

## 5. Global comparison with the corner

Starting from any strict anchor triple, move it along

\[
y_j(\theta)=1-(1-\theta)(1-y_j),\qquad0\le\theta<1.
\]

All three anchors strictly increase, remain distinct, and retain fixed
strict ratios of their distances from one. The inherited corner theorem
therefore gives `(A,B,eta)->(A*,B*,eta*)` as `theta->1`.

For any fixed `0<t<y1`, equations (N-L6) show that `H(t)` strictly
increases along this path, and that `Y0` strictly increases as well.
At `t=y1` the starting primitive is zero and `Hstar(y1)>0`. Taking the
corner limit proves

\[
\boxed{0<H(t)<H_*(t)\quad(0<t<y_1),\qquad
Y_0<Y_{0,*}=-3/1232.} \tag{N-L7}
\]

The first inequality extends to `H(y1)=0<Hstar(y1)` at the first root.
The proof uses the fixed-ratio corner only as a comparison destination
reachable from every strict triple. It never asserts that a given partial
endpoint approach has that corner as its own limit.

## 6. Finite-kappa corner negativity and the global theorem

The independent lift analysis in `notes_N_double.md` proves, for every
`0<a<1` and `0<t<1`,

\[
0<W_a(t)<W_1(t),\qquad
W_1(t)=\frac{3}{2}\,
\frac{(1-t)^{-2/3}-1}{1152t^2(1-t)^{13/6}}. \tag{N-L8}
\]

This pointwise comparison does not require monotonicity in `a`. Its proof
uses the positive causal Green kernel and an explicit positive residual;
no numerical kernel maximization enters.

The exact inherited beta moments give

\[
\int_0^1 W_1(t)H_*(t)\,dt
=\frac32\frac{25-1}{14784}=\frac3{1232}=-Y_{0,*}. \tag{N-L9}
\]

This integral converges: near the loop its integrand is a constant times
`(1-t)^(-5/6) log(1/(1-t))`, which is integrable.
Thus the requested finite-kappa corner statement is proved exactly:

\[
\boxed{\Phi_*(1;a):=Y_{0,*}+\int_0^1W_aH_*<0\quad(a<1).} \tag{N-L10}
\]

More strongly, for every strict lobe point, (N-L7)--(N-L10) give

\[
\begin{aligned}
\Phi_a(y_1)
&=Y_0+\int_0^{y_1}W_aH\\
&<Y_{0,*}+\int_0^{y_1}W_aH_*\\
&<Y_{0,*}+\int_0^1W_aH_*\\
&=\Phi_*(1;a)<0.
\end{aligned} \tag{N-L11}
\]

Therefore **Theorem N holds on the full lobe region for every finite
`kappa>1`**. In fact (N-L7) and (N-L9) also give the strictly negative
limiting functional at `a=1` for each fixed strict anchor triple. Equality
is approached only through a boundary comparison, not attained at any
strict triple covered by the theorem.

For `a` in any fixed compact interval `[amin,amax]` inside `(0,1)`, dominated
convergence using `W1*Hstar` makes `Phi_star(1;a)` continuous. Its strict
negativity gives a uniform negative upper bound on that entire lift
interval. Equation (N-L11) transfers the bound to every strict anchor
triple, including partial endpoint strata, without asserting an incorrect
corner expansion for those strata. No numerical value such as `-0.63` is
needed for the proof.

This note proves the target scalar inequality. The root report separately
states its original-integral zero-count consequence and distinguishes
sign changes, simple zeros, distinct zeros, and multiplicities.

## 7. Independent audit and bounded replay

The strict moment-curve convexity and endpoint slope were independently
re-derived by the geometry agent. The cardinal signs, center functional,
rational residual `25/231`, path comparison, and combination with the
kernel bound were independently checked by the audit agent. Both found
that the global comparison closes Theorem N.

`q4_N_loop_checks.py` replays the exact endpoint moments, center-functional
coefficient ratios, and rational positive residual with one numerical
thread and a ten-second CPU fuse. It performs no scan, optimization,
tangency search, or reuse of the previous frozen shooting trials.
