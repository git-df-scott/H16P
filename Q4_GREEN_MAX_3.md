# Q4 third strike: the first Green maximum and the joint endpoint limit

2026-09-04. **No five-original-zero candidate or certificate was obtained.**
The exact delayed-anchor path crosses the necessary first-root threshold.
Eight targeted lift shots reached four intermediate derivative crossings
numerically, but every first Green maximum remained negative. Four
primitive-contact shots and four reverse Green-tangency lines likewise
produced no successful construction.

The new proved result is narrower than a global exclusion: for every fixed
strict endpoint anchor-ratio triple and every fixed finite positive ratio
`lambda=1/(kappa*epsilon)`, sufficiently late anchors fail the first
required derivative minimum. The proof is locally uniform for lambda in
compact subsets of `(0,infinity)`. It covers neither escaping lambda nor
degenerating anchor ratios, and it does not settle saturation of the
original inverse operator's `+2` allowance.

The derivation below distinguishes proved identities and asymptotics from
frozen numerical diagnostics. Its independent hostile review is recorded in
[q4/notes_audit_third.md](q4/notes_audit_third.md). The full original
coefficient transport and corrected Picard–Fuchs forcing sign remain in
[Q4_RECONSTRUCTION_GEOMETRY.md](Q4_RECONSTRUCTION_GEOMETRY.md).

## 1. Exact target and escape conditions

The first primitive root means the smallest zero of
`H(t)=integral_0^t uF(u)q(u)du`, not a root of `q`, `P`, or the original
integral. The relevant first maximum is the first local maximum of
`Z=Y/y`, at the first zero of its weighted derivative `P`, where `Y=G/C`
and `y>0` is the homogeneous solution vanishing at the homoclinic endpoint.

The second-strike theorems prove the strict necessary conditions

\[
t^H_1>5/11,\qquad \kappa>21636/19043.
\]

Neither is sufficient. The non-strict complementary intervals, including
equality, were excluded. The exact rational and integer threshold regressions
already exist in `q4_reconstruction.py`; no replacement threshold is inferred
here. The full shooting conditions (S1)–(S3) remain binding.

The controlled path in this subtask is

\[
(t^H_1,t^H_2,t^H_3)=
\left(r,\frac{1+r}{2},\frac{3+r}{4}\right),\quad0<r<1.
\]

The first root is identically `r`; the exact anchor-coordinate theorem puts
every path point in the strict lobe region. Thus the path crosses `5/11`
constructively and continuously. The stable moment evaluator is
`q4_threshold_path.py`; coefficient and sign certification is owned by the
independent geometry/certificate work.

### The exact sequential shooting target

For clarity, use `k=kappa`, `a=1-1/k`, `s=k-(k-1)t`, and
`C=pi/sqrt(k-1)`. The universal and original reconstruction quantities are

\[
F(t)={}_2F_1(1/6,5/6;1;t),\qquad
M(t)=1-6(1-t)F'(t)/F(t),
\]
\[
q(t)=A+Bt-1+(t-\eta)M(t),\qquad
H(t)=\int_0^t uF(u)q(u)\,du.
\]

The original transported center data, with `Y=G/C`, are fixed:

\[
Y_0=\frac{3(1326A+864B-2431\eta-102)}{1361360}<0,
\qquad Y_1=-\frac32(1+a)Y_0-\frac\eta{192}.
\]

The exact reconstruction is

\[
(1-at)(1-t)Y''-\frac{1-a}{2}Y'+\frac{5a}{36}Y
=-\frac{H}{1152t^2(1-t)},
\]
\[
X(t)=\int_0^t\frac{Y(u)}{(1-au)^{3/2}}\,du,
\qquad I(t)=-\frac{aC}{2}\sqrt{1-at}\,X(t).
\]

The forcing sign here is the independently verified correction of the
source's printed sign; the original area evaluation checks below test
this normalization again. No free homogeneous constant is available.

Let `x(s)=asinh(sqrt(s-1))` and

\[
O(x)=\frac3{10}\sinh(5x/3)+\frac32\sinh(x/3),\qquad
E(x)=\frac{5\cosh(x/3)-\cosh(5x/3)}4.
\]
\[
y(t)=\frac{O(x(s(t)))}{O(x(k))}>0,\qquad
p(t)=\sqrt{\frac{1-t}{1-at}},\quad
Z=Y/y,\quad P=py^2Z'.
\]

Writing `r(a)=y'(0)`, the scalar initial datum and forcing are

\[
P_0=Y_1-r(a)Y_0,\qquad P'=-\Omega H,
\quad\Omega=\frac{y}{1152t^2(1-at)^{3/2}(1-t)^{3/2}}>0.
\]

Suppose the universal point lies in the strict lobe region, so `H` has
three simple roots `tau1<tau2<tau3` and signs `+,-,+,-`. Put
`V_H(t)=integral_0^t Omega(u)H(u)du`. For each fixed `a<1`, the loop
limits are `P->+infinity`, `Z->+infinity`, and finite `X(1)`. These
follow from `H(1)<0`, `y` asymptotic to a positive multiple of
`sqrt(1-t)`, and the resulting integrable logarithmic divergence of `Y`.
The exact necessary and sufficient sequential tests are:

\[
\boxed{\max\{0,V_H(\tau_2)\}<P_0<
\min\{V_H(\tau_1),V_H(\tau_3)\}.} \tag{S1}
\]

This is equivalent to four simple roots `p1<p2<p3<p4` of `P`.
In particular, the first local minimum `P(tau1)` must be **negative**.
If (S1) holds, five simple zeros of `Y=yZ` are equivalent to

\[
\boxed{Z(p_1)>0,\quad Z(p_2)<0,\quad Z(p_3)>0,\quad Z(p_4)<0.} \tag{S2}
\]

The first required positive Green maximum is `Z(p1)`. It is a different
quantity from the first minimum `P(tau1)`. A positive first minimum of
`P` is an obstruction to (S1); four crossings of `P` alone do not imply
the positive maximum required by (S2).

If (S1) and (S2) hold and the five roots of `Y` are `v1<...<v5`, the
original integral has five distinct simple interior zeros precisely when

\[
\boxed{X(v_2)>0,\quad X(v_3)<0,\quad X(v_4)>0,\quad
X(v_5)<0,\quad X(1)>0.} \tag{S3}
\]

Each equivalence follows by strict monotonicity on the successive lobes,
with `Z(0)=Y0<0` and `X(0)=0`. The inherited multiplicity-five upper
bound makes these conditions necessary for any five distinct original
zeros as well. All three gates remain binding in this strike.

A useful exact first-maximum form follows by defining

\[
\mathcal R(t)=\int_0^t\frac{du}{p(u)y(u)^2}
=\frac{2O(x(k))^2}{\sqrt{k(k-1)}}
\left[\frac{E(x(s(t)))}{O(x(s(t)))}-
\frac{E(x(k))}{O(x(k))}\right].
\]

At the first `P` root before `tau1`, integration by parts gives

\[
P_0=\int_0^{p_1}\Omega H,\qquad
Z(p_1)=Y_0+\int_0^{p_1}\mathcal R\Omega H.
\]

Thus the first positive derivative lobe must overcome the exact negative
center offset through this weighted moment. The independent audit derives
the formula and uses it for the four fixed reverse-tangency attempts.

An ordinary original fold remains an exact alternative target: under
(S1), (S2), and `X(1)>0`, replacing exactly one of the four strict
extremum signs in (S3) by equality, with the other three strict, gives
three simple roots and one ordinary double. Indeed
`X''(vj)=Y'(vj)/(1-avj)^(3/2)` is nonzero. A transverse variation of
that height splits the double; the original `mu1` direction is transverse
because its basis function `h I00` is nonzero in the open annulus. No
such original fold point was found.

## 2. The endpoint coefficient direction and its original normalization

Put `epsilon=1−r`, `L=log(432/epsilon)`, and `c=(1−t)/epsilon`.
The three anchors correspond to `c=1,1/2,1/4`.

**PROVED from the scaled anchor system** in the geometry notes: the
normalized coefficient limit is

\[
(A_*,B_*,\eta_*)=(94/77,-17/77,1).
\]

The exact coefficient transport gives, for every finite `k=κ>1`,

\[
\begin{aligned}
\mu_{1,*}&=\frac{27(8k-9)}{49280(k-1)},&
\mu_{2,*}&=-\frac{9(8k-27)}{24640k},\\
\mu_{3,*}&=\frac{9(8k-9)}{24640(k-1)},&
\mu_{4,*}&=-\frac{81(2k-3)}{49280k(k-1)}.
\end{aligned}
\]

In particular `mu3=2mu1/3`, and the reconstruction center datum has the
especially simple value

\[
Y_{0,*}=-\frac3{1232}.
\]

The closed PF moment identities simplify the limiting primitive to

\[
\boxed{H_*(t)=\frac{6t(1-t)^2}{77}
\left[5F(t)-36(1-t)F'(t)\right]
=\frac{6t(1-t)^2F(t)}{77}(6M(t)-1)>0.} \tag{G1}
\]

Its positivity follows from the positive Stieltjes representation:
`M(t)>M(0)=1/6` for `0<t<1`. Thus this boundary direction has no interior
primitive zero, although its nearby three-anchor directions do.

The initial weighted derivative is

\[
P_{0,*}(a)=\frac3{1232}\left[C_a-\frac{77}{36}\right],
\qquad C_a=\frac32(1+a)+r(a),\quad a=1-1/k.
\]

## 3. Two exact cancellations at the joint kappa-infinity endpoint

At `a=1`, interpreted as a compact-interval limit, the positive homogeneous
factor and its integrating factor are

\[
y(t)=(1-t)^{5/6},\qquad p(t)=1,
\qquad\Omega_1(t)=\frac{(1-t)^{-13/6}}{1152t^2}.
\]

The center data give `P0*=1/14784`. The following are exact:

\[
\boxed{\int_0^1\Omega_1H_*\,dt=\frac1{14784}},\qquad
\boxed{\int_0^1(1-t)^{-2/3}\Omega_1H_*\,dt=\frac{25}{14784}}. \tag{G2}
\]

Here is a short positive-series proof. If
`F=sum f_n t^n`, its coefficient recurrence gives

\[
5F-36(1-t)F'=
\sum_{n\ge1}\frac{5n}{n+1}f_nt^n.
\]

Termwise beta integration is justified by positivity. With endpoint weight
`(1-t)^(−1/6)`, cancellation of the `(5/6)_n` factor reduces the integral
to `5 sum[n≥1](1/6)_n/(n+1)!=1`; with weight `(1-t)^(−5/6)`, cancellation
of `(1/6)_n` gives `5 sum[n≥1](5/6)_n/(n+1)!=25`. The unmultiplied sums equal the
integrals of `(1-t)^(−1/6)−1` and `(1-t)^(−5/6)−1` respectively.

Consequently

\[
P_*(1)=0,
\]

and the primitive of `1/y²` is
`R_1(t)=(3/2)[(1-t)^(−2/3)−1]`, so

\[
\boxed{Z_*(1)=-\frac3{1232}
+\frac32\frac{25-1}{14784}=0.} \tag{G3}
\]

Moreover `P_*(t)=integral_t^1 Ω1H*>0` for every `t<1`, hence `Z_*` is
strictly increasing from a negative value to zero. Therefore

\[
P_*(t)>0,\quad Z_*(t)<0\quad(0<t<1).
\]

This is a precisely balanced boundary, not a positive first maximum. It
explains why the coupled limit `r→1`, `k→infinity` deserves separate
analysis from a fixed moderate `k` shooting.

## 4. Fixed finite kappa: a distinct endpoint balance

For fixed `k`, the limiting equation has `H_*>0`, hence `P_*` is strictly
decreasing. A cluster of four `P` roots near the three endpoint primitive
anchors can only approach a value satisfying `P_*(1;k)=0`.

**NUMERICAL diagnostic only.** One bounded evaluation of this scalar
endpoint balance at the three specified values gave

| k | P0* | integral_0^1 Omega H* | P*(1) |
|---:|---:|---:|---:|
| 8 | 3.4274121102e−5 | 4.3136198665e−5 | −8.8620775635e−6 |
| 9 | 6.4983521363e−5 | 4.4651735479e−5 | +2.0331785884e−5 |
| 16 | 1.5230831113e−4 | 5.1569191547e−5 | +1.0073911958e−4 |

The adaptive quadrature's reported errors were below `5e−14`, but these
are not rigorous enclosures. The command used one numerical thread and a
ten-second CPU ceiling and completed in about 0.42 seconds. This answers
the narrow question whether the moderate-k shooting branch is distinct
from the joint infinite-k boundary: the diagnostic is consistent with a
finite endpoint balance between eight and nine, while (G2) proves another
limiting balance at infinity. No new original root count is inferred.

## 5. Endpoint-scaled forcing and a useful general sign lemma

To make the cancellation precise, define the exact scaled coefficients by

\[
A+B-\eta=\epsilon LQ_\epsilon,\qquad
\eta-1=\frac{\epsilon L}{6}
       [V_\epsilon-(L+1)Q_\epsilon],\qquad
2\pi H(1)=\epsilon^2Le_\epsilon.
\]

The geometry subtask proves boundedness and convergence of these three
scaled quantities, with limits `Q,V,e`, and gives the scaled forcing limit

\[
\frac{2\pi H(1-\epsilon c)}{\epsilon^2L}\longrightarrow
h(c)=e-Vc+Qc\log c+Dc^2,
\]

uniformly on compact positive-c intervals, where for the specified path

\[
D=\frac{30}{77},\quad e=-\frac{15}{154},\quad
V=\frac{45}{154},\quad Q=-\frac{45}{154\log2}.
\]

The scaled anchor equations are the invertible three-by-three linear
system `h(1)=h(1/2)=h(1/4)=0`.

The same coefficient-limit argument applies to every fixed strict triple
`(1,c2,c3)`: the interpolation functions `1,c,c log c` form an extended
complete Chebyshev family because their full Wronskian is `1/c>0`.
The limiting three-by-three matrix is therefore invertible, and the
coefficientwise endpoint remainder tends to zero after the same scaling.
Matrix inversion proves boundedness and convergence of the scaled
coefficients `e,V,Q`; those properties are not additional assumptions.

**PROVED sign lemma for any fixed anchor ratios.** Suppose instead that
`h` has its three prescribed zeros at `0<c3<c2<1`, with `D>0` and largest
zero normalized to one. Then

\[
-D<e<0,\qquad Q<0,\qquad0<V<D.
\]

Indeed `h(c)/c=e/c−V+Q log c+Dc` has three zeros. Its derivative has two
zeros `xi1,xi2` in `(0,1)` and numerator `Dc²+Qc−e`. Therefore
`−e/D=xi1 xi2∈(0,1)` and `Q=−D(xi1+xi2)<0`. The equation `h(1)=0`
then gives `V=e+D∈(0,D)`. In particular

\[
e-Vc+Qc\log c<0\qquad(c\ge1). \tag{G4}
\]

## 6. Proved matched joint scaling: the first S1 minimum has the wrong sign

The relevant joint parameter is
`lambda=(1−a)/epsilon=1/(k epsilon)`. Let `lambda` tend to a fixed
positive value while `epsilon→0`. Put

\[
C_O=\frac3{5\,2^{1/3}},\quad
x=\operatorname{arsinh}\sqrt{c/\lambda},\quad
\widetilde y_\lambda(c)=\frac{\lambda^{5/6}}{C_O}O(x),
\]
\[
\omega_\lambda(c)=
\frac{\widetilde y_\lambda(c)}{(c+\lambda)^{3/2}c^{3/2}}.
\]

The normalized positive homogeneous factor is asymptotic to
`epsilon^(5/6) ytilde`, and `1152 Omega` to
`epsilon^(−13/6) omega`.

**PROVED; independently audited.** For each fixed positive lambda and each
fixed strict anchor-ratio triple, the matched formula is

\[
\frac{2304\pi P(1-\epsilon c)}{\epsilon^{5/6}L}
\longrightarrow
\mathcal P_\lambda(c)
=D\mathcal B_\lambda(c)
-\int_c^\infty\omega_\lambda(v)
       [e-Vv+Qv\log v]\,dv. \tag{G5}
\]

An explicit primitive for the growing `Dc²` contribution is

\[
\mathcal B_\lambda(c)=\frac{2\lambda^{5/6}}{C_O}
\left[\frac{23}{18}J_O(x)-O(x)\tanh x+\frac12O_x(x)\right],
\]
\[
J_O(x)=\frac9{50}\cosh(5x/3)+\frac92\cosh(x/3).
\]

Exact differentiation, using `O_xx−2 tanh(x)O_x+(5/9)O=0`, gives
`B'_lambda=omega_lambda c²`. Furthermore `B_lambda(0)>0` and its
derivative is positive, so `B_lambda(c)>0`. By (G4), the integral term
in (G5) also has positive sign for `c≥1`. Thus the proved matched
limit satisfies

\[
\mathcal P_\lambda(1)>0
\]

for every fixed positive lambda and every fixed strict anchor-ratio triple.
Consequently, for all sufficiently small positive epsilon,
`P(t1^H)=P(1−epsilon)>0`. The first strict (S1) sign requires that value
to be negative. Thus five original zeros are impossible in this controlled
joint endpoint regime. The conclusion is locally uniform when lambda ranges in a fixed compact
subset of `(0,infinity)`, with the strict anchor ratios fixed.
It does not exclude finite non-asymptotic delayed shapes, lambda tending to
zero or infinity, or degenerating anchor ratios.

### Proof of the matching constant and uniform remainder bounds

The matching constant must not be chosen freely. Let
`b=5/2^(4/3)` and

\[
R(u)=\frac{5F(1-u)-36uF'(1-u)}{1-u},
\]
\[
K_k(u)=\frac{k^{3/2}}{O(x(k))}
\frac{O(x(1+(k-1)u))\sqrt u}{[1+(k-1)u]^{3/2}}.
\]

Then the exact star forcing integral is
`(1/14784) integral_delta^1 K_k(u)R(u)du`, with

\[
K_k(u)=u^{-1/6}+b k^{-2/3}(u^{-5/6}-u^{-1/6})+cdots.
\]

The two exact moments in (G2) say that the full integrals of the two
displayed kernel terms are `1` and `24b k^(−2/3)`. They cancel the
corresponding center-data expansion

\[
P_{0,*}(k)=\frac1{14784}
\left[1+24b k^{-2/3}\right]+O(k^{-1}).
\]

After setting `u=epsilon v`, the renormalized inner kernel
`K_lambda(v)=omega_lambda(v)v²` has tail

\[
K_\lambda(v)=v^{-1/6}+b\lambda^{2/3}v^{-5/6}+O(v^{-7/6}).
\]

Consequently the growing primitive with no free constant is

\[
\frac65c^{5/6}+6b\lambda^{2/3}c^{1/6}
-\int_c^\infty
[K_\lambda(v)-v^{-1/6}-b\lambda^{2/3}v^{-5/6}]\,dv.
\]

It has the same derivative and infinity expansion as the explicit
`B_lambda`, whose hyperbolic expression has no constant term in its large-c
power expansion. The integral remainder has an integrable `v^(−7/6)`
tail. The parameter perturbation contributes the integrable
`v^(−7/6)log v` tail in (G5); its center-data difference is
`O(epsilon L²)=o(epsilon^(5/6)L)`.

Here are explicit bounds that close the overlap-region argument. Write

\[
E_k(u)=K_k(u)-u^{-1/6}
-b k^{-2/3}(u^{-5/6}-u^{-1/6}).
\]

For `u∈[epsilon c,1]`, fixed `c>0`, and lambda in a fixed positive compact
set, the elementary hyperbolic asymptotic expansion gives

\[
\boxed{|E_k(u)|\le C k^{-1}u^{-7/6}}. \tag{G6}
\]

One can verify this bound directly. With `s=ku+1−u`, use
`C_O k^(5/6)/O(x(k))=1−b k^(−2/3)+O(k^(−1))` and
`O(x(s))/(C_O s^(3/2))=s^(−2/3)+b s^(−4/3)+O(s^(−5/3))`.
Mean-value comparison of `s` with `ku` supplies the bound, since
`ku≥c/lambda`; the normalization errors are absorbed using `u≤1`.
Thus `epsilon^(1/6)|E_k(epsilon v)|≤C'v^(−7/6)`. The ratio
`R(epsilon v)/L` is bounded for `v≥c` and tends pointwise to `5/(2pi)`.
Dominated convergence proves the stated finite-part limit and fixes its
constant. The lower tail terms are dominated by
`v^(−1/6)(1+|log v|)` and `v^(−5/6)(1+|log v|)` on `(0,c]`.

For the coefficient perturbation, split the original integral at a fixed
small `u0`. The coefficientwise Frobenius expansion gives the exact leading
scaled expression plus a remainder bounded by

\[
C\epsilon L^2u^2(1+|\log u|).
\]

The leading `L²` terms cancel by the scaled coefficient definition
`eta−1=epsilon L[V_epsilon−(L+1)Q_epsilon]/6`.
On `epsilon c≤u≤u0<1`,
`Omega≤C u^(−13/6)`. Integrating the remainder over `(epsilon c,u0)`
therefore costs only `O(epsilon L²)`, because
`integral_0^u0 u^(−1/6)(1+|log u|)du` is finite. The outer interval costs
the same order by coefficient smoothness and the center cancellation
`H/t²=O(1)`. Both costs divided by `epsilon^(5/6)L` tend to zero.
The remaining inner integrand is dominated by
`v^(−13/6)(1+v+v|log v|)` on `[c,infinity)`, which is integrable.
This proves (G5), with no freely adjustable matching constant.

The algebraic finite-part primitive and its absence of a constant term are
also checked by `q4_green_endpoint_third.py`. The uniform bounds above are
analytic estimates; the script's exact equalities are not substitutes for
them. An independent agent reviewed (G6), both dominated-convergence steps,
the normalization constants, and the finite-nearby sign conclusion.

## 7. Construction implications and scope

The controlled path reaches the required delayed-root region exactly.
The remaining task is genuinely the Green maximum: neither the delay nor
three primitive zeros imply the needed positive first maximum. The coupled
endpoint has two exact cancellations, and its leading inner forcing has
an explicit sign structure. The fixed-ratio, fixed-positive-lambda regime
is now rigorously closed at (S1), while non-asymptotic delayed shapes and
different anchor separation or parameter rates remain open.

No positive first maximum, five-zero numerical candidate, ordinary original
fold, or rigorous five-zero certificate has been obtained in this analytic
subtask. The finite shooting and independent original-integral tests below remain
numerical evidence and do not promote their signs to proofs.

## 8. Frozen construction attempts after the threshold crossing

**NUMERICAL ONLY.** The next table is generated directly from the frozen
[third_tuned_shoot.json](q4/data/third_tuned_shoot.json) and
[third_shape_shoot.json](q4/data/third_shape_shoot.json). Let path A mean
`(r,(1+r)/2,(3+r)/4)` and path B mean
`(r,1-(1-r)^2,1-(1-r)^3)`. Each coefficient triple comes from the exact
primitive-anchor map; the displayed lift parameters and reconstructed
signs use numerical evaluation.

| Path | r | kappa | Numerical S1 | First maximum Z(p1) |
|---|---:|---:|:---:|---:|
| A | 0.5 | 2.1769936903 | passes | -0.003934275289 |
| A | 0.75 | 2.5851372245 | passes | -0.003420052400 |
| A | 0.9 | 3.2734857016 | passes | -0.002994048345 |
| A | 0.99 | 5.5639324237 | passes | -0.002532118321 |
| A | 0.9999 | 8.1071043201 | passes | -0.002362670910 |
| B | 0.6 | 2.4115997810 | passes | -0.003599582185 |
| B | 0.75 | 2.9750772377 | passes | -0.003139178640 |
| B | 0.9 | 4.2310079066 | passes | -0.002719393737 |

All four reconstructed extrema of `Z` were negative in every shot.
Consequently each numerical trial fails (S2), despite its four intermediate
`P` crossings. The original integral was positive at every saved sample;
those samples do not cover all endpoint slivers and are not a root-absence
certificate. The latest path-A point has a `P` extremum gap near
`2.09e-13`, with primitive residuals near `8e-17`; the shooting band is
numerically delicate. No interval claim is attached to it.

The replays impose one numerical thread and a ten-second CPU ceiling.
The saved path-A and path-B records used approximately 2.93 and 4.04 CPU
seconds respectively. They contain only the listed predetermined shapes
and targeted scalar lift shots.

### Four primitive-contact boundary shots

The frozen [third_confluent_shoot.json](q4/data/third_confluent_shoot.json)
records four targeted boundary attempts where the primitive has a triple
contact and the lift is tuned to `P` approximately zero there. These are
primitive contacts, not double zeros of the original integral.

| Primitive triple contact t | kappa | P at contact | Z at contact |
|---:|---:|---:|---:|
| 0.6 | 1.9620342563 | 0.000e+00 | -0.004397010273 |
| 0.9 | 2.7699150452 | 7.364e-19 | -0.003270467969 |
| 0.99 | 4.9082982650 | 1.885e-19 | -0.002607175260 |
| 0.9999 | 8.0098184603 | -4.141e-19 | -0.002367474135 |

The four values of `Z` remain negative. Thus these attempts do not supply
a boundary organizer with a positive first Green maximum or an ordinary
original fold. The saved bounded replay used approximately 0.46 CPU
seconds. Its tuned parameters and signs are numerical diagnostics.

### Reverse tangency and direct original-integral checks

The [independent audit](q4/notes_audit_third.md) derives the exact
first-maximum moment and investigates exactly four prescribed reverse
Green-tangency lines: `(a,t*)=(.75,.8),(.75,.95),(.9,.999),(.99,.999)`.
The first two failed elementary necessary lobe inequalities numerically;
the latter two yielded no detected three-primitive-root interval in the
bounded line investigation. The saved results are
[third_reverse_tangency.json](q4/data/third_reverse_tangency.json).
Finite meshes can miss roots and do not rigorously exclude whole lines or
the full delayed lobe region.

Two targeted evaluations compare the scalar PF reconstruction with the
original area-integral evaluator and the transported original coefficient
vector, as saved in
[third_independent_checks.json](q4/data/third_independent_checks.json):

| Trial | kappa | t | PF original I | Independent original area I | Absolute difference |
|---|---:|---:|---:|---:|---:|
| A, r=.75 | 2.5851372245 | .7 | .0015492881998878836 | .0015492881998878847 | 1.09e-18 |
| B, r=.9 | 4.2310079066 | .95 | .0012121057635914785 | .0012121057635913627 | 1.16e-16 |

Both coefficient points survive the corrected necessary filters. These
checks support the coefficient transport, forcing sign, and original
normalization; they certify neither a root count nor a five-zero lead.

## 9. Verification and remaining construction domain

The bounded exact replay
[q4_green_endpoint_third.py](q4/q4_green_endpoint_third.py) verifies the
limiting primitive, original center datum, beta moments `1` and `25`,
both endpoint cancellations, the finite-part primitive, its zero constant
at infinity, and cancellation of the first nonanalytic kappa correction.
It completes in about half a CPU second. The uniform matching estimates
are proved above and independently audited; a symbolic identity check is
not being used in place of those estimates.

The delayed-root path and rational late-root certificate are documented in
[Q4_THRESHOLD_PATH.md](Q4_THRESHOLD_PATH.md). The full strike assessment is
[ASTRA_THIRD_STRIKE.md](ASTRA_THIRD_STRIKE.md). The present proof closes only
the controlled regime of fixed strict anchor ratios and fixed finite
positive `lambda`, locally uniformly on compact lambda intervals. Finite
delayed shapes, lambda tending to zero or infinity, and anchor ratios
coalescing or separating at additional rates remain outside that proof.

No positive first Green maximum, five-zero numerical candidate,
three-simple-plus-double original point, or rigorous five-zero certificate
was produced. Saturation of the original reconstruction allowance and the
global Q4 five-zero problem remain unknown.
