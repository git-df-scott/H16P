# Third strike: delayed-anchor Green reconstruction analysis

This subtask continues the exact second-strike normalization. It does not
revisit the old coefficient box or auxiliary cusp. No five-zero parameter is
claimed. Status distinctions below separate exact identities, proved signs,
matched asymptotics, and bounded numerical evidence.

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
joint endpoint regime. The conclusion is locally uniform when lambda and
the anchor ratios range in fixed compact subsets of their strict domains.
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
subtask. The root agent owns the finite shooting and independent-original-
integral tests; this note does not promote their numerical signs to proofs.
