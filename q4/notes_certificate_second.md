# Second strike: an exact rational weighted-lobe certificate

Date: 2026-09-04. **RIGOROUS COMPUTATION:** the universal weighted-lobe
region contains the rational point and the explicit coefficient box below.
This is a certificate of three simple zeros of the universal primitive.
It is not a certificate of five zeros of the original Abelian integral.

The definitions and the analytic ECT theorem are inherited from
`Q4_STRUCTURE.md` and `Q4_ZERO_GEOMETRY.md`:

\[
F(t)={}_2F_1(1/6,5/6;1;t),\qquad
q(t)=A+Bt-1+(t-\eta)M(t),\qquad
H(t)=\int_0^t uF(u)q(u)\,du.
\]

## 1. Frozen rational point and certified signs

\[
\boxed{
A=\frac{1243911778077}{10^{12}},\qquad
B=-\frac{86917392526}{10^{12}},\qquad
\eta=\frac{1460428426173}{10^{12}}.}
\]

The exact corrected strip `1 < eta < 54/31` holds. The exact initial sign is

\[
q(0)=A-1-\eta/6=\frac{1014080763}{2000000000000}>0.
\]

`q4_lobe_certificate.py` evaluates the following outward decimal intervals
using only rational and integer arithmetic:

| t | Lower bound for H(t) | Upper bound for H(t) |
|---|---:|---:|
| 1/8 | 0.000001441401383769815937769764 | 0.000001441401383769815937769765 |
| 3/8 | -0.000003720918085980812785719090 | -0.000003720918085980812785719089 |
| 5/8 | 0.000016780934607604637814181443 | 0.000016780934607604637814181444 |
| 7/8 | -0.000364906970505491567967035918 | -0.000364906970505491458511914236 |

The finite series ends at index `N=256`. The largest absolute analytic
tail bound is less than `5.473e-20`; the first three tails are smaller
than `1e-30`. The displayed endpoints are themselves exact decimal
rationals, rounded down/up after computing the analytic tail bound.

The frozen result is `data/second_lobe_certificate.json`, including the
verifier's SHA-256 hash. One full replay used about 0.04 CPU seconds.
The script enforces a ten-second CPU ceiling and requires no third-party
libraries. This bounded replay is neither a scan nor a numerical fit.

## 2. Why these four signs prove all three strict lobe inequalities

The auxiliary space `span{1,t,M,tM}` has at most three zeros counting
multiplicity. Since `H(0)=0` and `H'(t)=tF(t)q(t)` with positive interior
weight, anchored Rolle implies `Z(H;(0,1)) <= 3`.

The four signs supply one root in each of

\[
(1/8,3/8),\qquad(3/8,5/8),\qquad(5/8,7/8).
\]

Therefore these are exactly three distinct simple primitive roots. The
root at zero and those three roots give, by ordinary Rolle, three distinct
interior roots of `q`; the ECT bound makes those roots simple and excludes
all others. Thus the hypotheses of the first-strike lobe equivalence hold.
The three simple primitive roots prove all three strict inequalities.

Writing the auxiliary roots as `x1 < x2 < x3`, the initial sign is positive,
so the inequalities in this certificate are exactly

\[
H(x_2)<0,\qquad H(x_3)>0,\qquad H(1)<0.
\]

No approximate auxiliary-root location or quadrature at a critical point
is needed for that conclusion. This uses the proved zero bound and the
strict lobe equivalence, not an inference from a plotted graph.

## 3. Positive period coefficients give a rational error bound

Write

\[
F(t)=\sum_{n\ge0}f_nt^n,\qquad f_0=1,\qquad
\frac{f_{n+1}}{f_n}=\frac{(6n+1)(6n+5)}{36(n+1)^2}<1.
\]

The companion period is
`K(t)={}_2F_1(-1/6,1/6;1;t)`. For `n>=1`, its coefficient is

\[
k_n=-\frac{f_n}{6n-1}.
\]

This follows by cancelling the two rising factorials whose arguments differ
by one. Hence

\[
D(t):=F(t)-K(t)=\sum_{n\ge1}d_nt^n,\qquad
d_n=\frac{6n}{6n-1}f_n>0.
\]

The first-strike companion identity gives `t F(t) M(t)=D(t)`. It follows
that the exact primitive is the absolutely convergent series

\[
\begin{aligned}
H(t)={}&\sum_{n\ge0}f_n\left[
 \frac{(A-1)t^{n+2}}{n+2}+\frac{Bt^{n+3}}{n+3}\right]\\
&+\sum_{n\ge1}d_n\left[
 \frac{t^{n+2}}{n+2}-\frac{\eta t^{n+1}}{n+1}\right].
\end{aligned}
\]

Let `S_N` include both sums through index `N>=1`. Because `f_n` decreases
and `6n/(6n-1)` decreases for `n>=1`, every omitted term is bounded by a
geometric majorant. Explicitly,

\[
\begin{aligned}
|H(t)-S_N(t)|\le{}&\frac{f_{N+1}t^{N+1}}{1-t}\left[
\frac{|A-1|t^2}{N+3}+\frac{|B|t^3}{N+4}\right.\\
&\left.\qquad+
\frac{6(N+1)}{6(N+1)-1}
\left(\frac{t^2}{N+3}+\frac{|\eta|t}{N+2}\right)\right].
\end{aligned}
\]

For rational inputs every quantity in this enclosure is rational. No
floating-point tolerance enters the proof. Cancellation in `S_N` is exact.

## 4. An explicit full coefficient box lies inside the lobe region

The certificate actually proves more than a single point. Let the center
be the boxed rational point above and allow

\[
|\Delta A|,\ |\Delta B|,\ |\Delta\eta|\le10^{-7}.
\]

The first-strike positive-measure representation gives `0<M(t)<=1`.
The decreasing positive Taylor coefficients give `F(t)<=1/(1-t)`.
For any `t<=7/8`, linearity in the three coefficients therefore implies

\[
\begin{aligned}
|\Delta H(t)|
&\le10^{-7}\int_0^t uF(u)(2+u)\,du\\
&\le10^{-7}\left[\frac{t^2}{1-t}
 +\frac{t^3}{3(1-t)}\right]\\
&\le\frac{1519}{192}10^{-7}<8\cdot10^{-7}.
\end{aligned}
\]

All four frozen sign enclosures have absolute margin greater than `10^-6`.
Thus every point of this closed rational box preserves the four alternating
signs and belongs to the strict weighted-lobe region. In particular the
center is an explicit interior point. These assertions, including the
rational perturbation bound, are recorded in the frozen JSON.

## 5. Exact endpoint functional: the third lobe boundary is a plane

**PROVED.** The endpoint primitive has the particularly simple form

\[
\boxed{
H(1)=\frac{18}{85085\pi}
\left(9061A+6289B-2431\eta-7242\right).}
\]

For the frozen point this gives

\[
\pi H(1)=-\frac{1908010250631}{132945312500000}<0.
\]

Here is a derivation independent of the finite-series sign calculation.
Let `J_n(t)=integral_0^t u^n F(u) du`. Integration by parts in the exact
formula `M=1-6(1-t)F'/F` yields

\[
H(t)=-6t(t-\eta)(1-t)F(t)-6\eta J_0(t)
 +(A+11+11\eta)J_1(t)+(B-17)J_2(t).
\]

Each `J_n` also has the closed hypergeometric form

\[
J_n(t)=\frac{t^{n+1}}{n+1}
{}_3F_2\left(\frac16,\frac56,n+1;1,n+2;t\right).
\]

The self-adjoint Gauss equation is
`[t(1-t)F']'=(5/36)F`. Its endpoint expansion gives
`lim[t->1]t(1-t)F'=1/(2 pi)`. Multiplying by `t^n` and integrating
twice by parts yields

\[
\left[n(n+1)+\frac5{36}\right]J_n(1)
=\frac1{2\pi}+n^2J_{n-1}(1),
\]

where the `n=0` case omits the final term. Consequently

\[
\pi J_0(1)=\frac{18}{5},\quad
\pi J_1(1)=\frac{738}{385},\quad
\pi J_2(1)=\frac{113202}{85085}.
\]

The boundary term in the formula for `H` vanishes because
`(1-t)F(t)->0`. Substitution gives the boxed rational linear functional.
On the normalized three-root auxiliary chamber the initial sign is positive,
so the final strict lobe inequality is the explicit half-space

\[
9061A+6289B-2431\eta-7242<0.
\]

The other two weighted-lobe conditions remain necessary: this endpoint
half-space alone does not establish three primitive zeros.

## 6. Scope of the result

The coefficient point was selected by the root agent's small numerical
three-anchor solve. Its frozen decimal-rational replacement was then
certified independently by the exact series above. The numerical discovery
values were not used as error bounds.

For any rational `kappa>1`, the first-strike invertible rational coefficient
map transports these rational normalized coefficients to a rational original
Abelian-integral coefficient vector. The resulting primitive has three
simple zeros for every such center parameter. Nothing here proves that the
remaining `kappa`-dependent reconstruction adds two zeros.

**NOT PROVED HERE:** five original Abelian-integral zeros, a three-simple-
plus-double original boundary point, an original quadratic perturbation
arc, or five quadratic limit cycles.

## 7. Final independent hostile audit

The completed `Q4_LOBE_REGION.md` and `Q4_RECONSTRUCTION_GEOMETRY.md`
were independently checked against the first-strike identities and this
certificate. No unresolved correctness issue was found in the analytic
cell/anchor inverse, bounded closure, boundary-event classification,
confluent contact arguments, or the stated reconstruction exclusions.

The boundary classification states necessary event types and explicitly
does not identify entire affine planes with boundary faces. It does not
assume a smooth compactification at the logarithmic endpoint. The explicit
coefficient box remains a valid lobe certificate, but subsequent analysis
excludes every original lift of the box, for every `kappa>1`.

The reconstruction forcing sign was checked directly:
`D_h=(2s/h)D_s` gives `L_h=-(16s/h)L_s`; hence the scalar `Y` equation
has forcing `-H/[1152 t^2(1-t)]`. The opposite sign printed in the
source's equation (24) is incompatible with that substitution.

The shooting conditions (S1)--(S3) have the required endpoint hypotheses.
They are used only in the strict lobe region, where `H(1)<0`, and for a
fixed `0<a<1`. Put `v=1-t` and write `y(t)~c sqrt(v)`, `c>0`. Then

\[
\Omega(t)\sim\frac{c}{1152(1-a)^{3/2}v},\qquad
P(t)\sim\frac{-cH(1)}{1152(1-a)^{3/2}}\log\frac1v>0.
\]

Using `Z'=P/(p y^2)` and `p~sqrt(v)/sqrt(1-a)` gives `Z(t)->+infinity`
and, more precisely,

\[
Y(t)=\frac{-H(1)}{576(1-a)}\log\frac1{1-t}+O(1).
\]

Thus `X(1)` is finite, because the logarithm is integrable and
`1-at` remains positive. No uniform assertion at `a=1` is used.
The strict alternating-extremum tests consequently have their stated
necessary-and-sufficient meanings at each successive stage.

The independent first-maximum estimate strengthens the original box check:
five original distinct zeros require the first primitive root to exceed
`5/11`. Indeed, if it did not, the necessary first zero of `P` would obey
`p1<5/11`, and

\[
\int_0^{p_1}\frac{dt}{p y^2}
<\frac67\left[\left(\frac{11}6\right)^{7/6}-1\right]<\frac89.
\]

Together with `P0<(9/4)|Y0|-eta/192` and
`|Y0|/eta<601/136136<1/216`, this forces the first `Z` maximum to be
negative, contradicting five original roots. The last fractional-power
inequality reduces to the integer check `1948617<2000000`.

The certified box has its first primitive zero below `3/8`, so it fails
this necessary condition uniformly in `kappa`. This closes the constructed
candidate region; it does not exclude the remaining late-root portion of
the lobe region. The intermediate sign-chain bounds count sign changes,
including odd multiple zeros, rather than inferring a global
multiplicity-three theorem.
