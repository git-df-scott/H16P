# Strike 4: a global anchor envelope replaces compact certification

2026-09-04. **PROVED analytically.** The requested compact-domain bound
follows from a strict, global comparison on the primitive anchor cell.
There is no coefficient scan, interval cover, or new numerical shooting in
this argument. The kernel comparison is proved independently in the
Strike-4 audit notes; the argument below supplies its universal forcing
and center-data comparisons.

The inherited definitions and anchor theorem are in
[Q4_LOBE_REGION.md](../Q4_LOBE_REGION.md). The endpoint coefficient limit,
valid for every fixed strict anchor-ratio triple, is inherited from
[Q4_GREEN_MAX_3.md](../Q4_GREEN_MAX_3.md), sections 4–5. Only this
pointwise-in-ratios limit is needed here; no uniform endpoint asymptotic is
assumed.

## 1. The two normalized moments form a strictly convex graph

Write

\[
 W(t)=tF(t)>0,\qquad
 x(t)=\frac{K_1(t)}{K_0(t)},\qquad
 m(t)=\frac{K_2(t)}{K_0(t)}.
\]

Thus `x` and `m` are the weighted averages of `u` and `M(u)` on `(0,t)`
with weight `W`. The positive Stieltjes representation gives `M'>0` and
`M''>0`. Direct differentiation yields

\[
 x'=\frac{W}{K_0}(t-x)>0,\qquad
 \frac{dm}{dx}=S(t):=\frac{M(t)-m(t)}{t-x(t)}.
\]

The quotient `S` is the weighted average of the strict secant slopes
`[M(t)-M(u)]/(t-u)`, now with positive weight `W(u)(t-u)`. Consequently
`0<S(t)<M'(t)`. Differentiating and cancelling the weighted-average
terms gives

\[
 S'(t)=\frac{M'(t)-S(t)}{t-x(t)}>0.
\]

Therefore `m` is a strictly convex function of `x`, without a separate
Wronskian or quadrature assumption. Its endpoint data are

\[
 x(0)=0,\quad m(0)=\frac16,\qquad
 x(1)=\frac{6289}{9061},\quad
 m(1)=\frac{2431}{9061}.
\]

Since `M(1)=1`,

\[
 \boxed{S(1^-)=\frac{1-2431/9061}{1-6289/9061}
                 =\frac{1105}{462}.} \tag{C1}
\]

## 2. A two-zero primitive has a positive center functional

Let

\[
 V=\alpha K_0+\beta K_1+\gamma K_2
\]

have two distinct interior roots `v<w` and be positive before `v`.
The three-function Chebyshev property excludes `gamma=0`. Dividing by
`K0` expresses `V/gamma` as the convex graph `m(x)` minus its chord at
`x(v),x(w)`. Its sign is positive outside the two anchors and negative
between them, so the specified orientation implies `gamma>0`.

Let that chord have intercept `b` and slope `s`. Strict convexity and
(C1) give

\[
 b<\frac16,\qquad 0<s<\frac{1105}{462}.
\]

Since `alpha=-gamma b` and `beta=-gamma s`,

\[
 \alpha>-\frac\gamma6,\qquad
 \beta>-\frac{1105}{462}\gamma.
\]

The center functional on this lower-dimensional primitive space is

\[
 \ell(V)=\frac9{3080}
       \left(\alpha+\frac{144}{221}\beta+\frac{11}{6}\gamma\right).
\]

The constants are exact reductions of the inherited `Y0` formula.
Therefore

\[
 \boxed{\ell(V)>
 \frac9{3080}\gamma
 \left(-\frac16-\frac{144}{221}\frac{1105}{462}+\frac{11}{6}\right)
 =\frac9{3080}\frac{25}{231}\gamma>0.} \tag{C2}
\]

This statement is the decisive additional sign information; positivity of
`V` on an initial interval alone would not imply positivity of an arbitrary
linear functional.

## 3. Every anchor increases the first lobe and the center value

Let `H=H_{r1,r2,r3}` be the normalized primitive with three ordered roots.
Its initial sign is positive and its root derivatives have signs `-,+,-`.
For `j=1,2,3`, differentiating the anchor equations gives

\[
 V_j:=\partial_{r_j}H\in\operatorname{span}\{K_0,K_1,K_2\},\quad
 V_j(r_i)=0\quad(i\ne j),\quad
 V_j(r_j)=-H'(r_j).
\]

The two-zero chord description in section 2 forces the coefficient of
`K2` in `Vj` to be positive in every case: `rj` is respectively left of,
between, or right of the other two roots, exactly matching the signs
`+,-,+` of `-H'(rj)`. Thus

\[
 \boxed{\partial_{r_j}H(t)>0\quad(0<t<r_1),\qquad
        \partial_{r_j}Y_0=\ell(V_j)>0.} \tag{C3}
\]

These are genuine coordinatewise derivative inequalities on the entire
open anchor simplex, not a statement restricted to a chosen spine.

For an arbitrary fixed starting triple, use the path

\[
 r_j(\theta)=1-(1-\theta)(1-r_j),\qquad0\le\theta<1.
\]

It stays in the open anchor simplex and increases every root. Its endpoint
distances have fixed distinct ratios. The inherited endpoint limit
therefore applies and gives

\[
 (A,B,\eta)\longrightarrow(94/77,-17/77,1),\quad
 Y_{0,*}=-\frac3{1232},\quad
 H_*(t)=\frac{6t(1-t)^2}{77}F(t)(6M(t)-1)>0.
\]

Integrating the strict inequalities (C3) along any initial segment of this
path and then taking the endpoint limit proves

\[
 \boxed{Y_0<Y_{0,*},\qquad 0<H(t)<H_*(t)\quad(0<t<r_1).} \tag{C4}
\]

No conclusion about convergence rates, degenerating anchor ratios, or
uniformity over starting triples is needed: each individual triple is
compared to the same endpoint object along its own fixed-ratio path.

## 4. Theorem N from the separately proved kernel comparison

Put `w_a=Rcal_a Omega_a`. The Strike-4 lift-kernel comparison is

\[
 0<w_a(t)<w_1(t)=\frac1{768t^2}
       \left[(1-t)^{-17/6}-(1-t)^{-13/6}\right]
 \quad(0<a<1,\ 0<t<1).
\]

The inherited endpoint cancellation is precisely

\[
 Y_{0,*}+\int_0^1w_1(t)H_*(t)\,dt=0.
\]

All these integrals converge: the integrand is `O(t)` at zero and
`O((1-t)^(-5/6) log(1/(1-t)))` at one. Applying (C4) only where `H>0`
now proves the strict inequality

\[
 \boxed{\begin{aligned}
 \Phi_a(r_1)
 &=Y_0+\int_0^{r_1}w_a H\\
 &<Y_{0,*}+\int_0^{r_1}w_1H_*\\
 &=-\int_{r_1}^1w_1H_*<0.
 \end{aligned}} \tag{C5}
\]

Thus the compact task closes by an analytic majorant which in fact proves
Theorem N on the whole product `L x (0,1)`, once combined with the
independently proved lift-kernel comparison. The first-root and small-kappa
thresholds are not required for (C5).

## 5. An explicit uniform compact margin

The comparison also gives an elementary quantitative result. Exactly,

\[
 w_1H_*=
 \frac{F(t)(6M(t)-1)}{9856t}
 \left[(1-t)^{-5/6}-(1-t)^{-1/6}\right].
\]

The positive hypergeometric coefficients and Stieltjes moments imply
`F>=1` and `M>=1/6+(25/432)t`. Hence for any `0<delta<1`, whenever
`r1<=1-delta` (in particular whenever `r3<=1-delta`),

\[
 \boxed{\Phi_a(r_1)<-\frac{25}{118272}
 \left(\delta^{1/6}-\frac15\delta^{5/6}\right)<0.} \tag{C6}
\]

This holds uniformly for every `0<a<1`, without an upper cutoff on `a`.
For example, `delta=1/64` gives the exact bound
`Phi_a(r1)<-395/3784704`. Only the displayed rational identities were
checked by a tiny standard-library `Fraction` calculation; none of the
analytic inequalities depends on a numerical estimate.

## 6. Handoff constant and consequence scope

The handoff affine constant needs the inherited `-K0` term:

\[
 c_0=-\frac{306}{1361360}+\int_0^{r_1}w_a(K_3-K_0),
\]

because `H=(A-1)K0+B K1-eta K2+K3`. No Claude-owned file is edited here.

The already verified implication (N1) says that five distinct original Q4
zeros would require `Phi_a(r1)>0`. Therefore (C5) rules out five distinct
original zeros. Any assertion of a global bound of three requires an
additional consequence argument, including coefficients outside the strict
lobe region and nonsimple original zeros; (C5) alone does not replace that
argument.

### Consequence audit: the proved scope and the remaining gate

**PROVED on the strict lobe region.** Theorem N bounds original sign
changes by three. When `P0<=0`, this is the already proved sign-chain
exclusion. When `P0>0`, `P` is decreasing on the first `H` lobe and is
positive near the loop because `H(1)<0`. Thus its number of sign changes is
even and at most four. With at most two, `Z'=P/(p y^2)` allows at most
three sign changes of `Z`. With four, let its roots be `p1<...<p4`.
The first maximum satisfies

\[
 Z(p_1)=\Phi(p_1)<\Phi(r_1)<0,
\]

so `Z` has no zero before `p2`; each of its remaining three monotone
pieces contributes at most one sign change. The positive factor `y` and
the center-anchored positive-weight primitive `X` do not increase that
number.

This gives at most three **distinct** original zeros on `L`, as follows.
Four distinct zeros, under the inherited total-multiplicity-five bound,
would be either four simple zeros or three simple zeros and one ordinary
double zero. The former contradict the sign-change bound. In the latter
case choose a small perturbation by the original nonvanishing basis
member `h I00` with the sign that splits the double zero into two simple
zeros. The other three simple roots persist. The perturbation stays in the
open coefficient region `L` and would create five simple roots,
contradicting Theorem N. This local argument uses only that `h` is nonzero
and `I00` is the positive oval area inside the annulus.

**GLOBAL conclusion from the supplied implications:** at most four
distinct original zeros. A stronger global conclusion of three has an
unaddressed case outside `L`; it is not a consequence of the preceding
paragraph. If `H` has two simple crossings with signs `+,-,+`, the sign
chain alone permits

\[
 P:\ +,-,+,-,\qquad Z:\ -,+,-,+,-.
\]

Indeed `P'=-Omega H` makes its three monotone pieces decrease, increase,
decrease, while `Z(0)<0` would allow four crossings if the three extremal
heights alternate. Theorem N assumes three primitive roots and supplies
no inequality on that two-root stratum. This is an abstract allowed sign
pattern, **not** a constructed Q4 example. It identifies exactly why
neither the first-lobe count nor ordinary anchored Rolle proves the claimed
global bound of three. Additional control of the two-root strata (and the
remaining projective chart, as needed) is required before claiming the
Gavrilov–Iliev/Zhao conjectured sharp bound.
