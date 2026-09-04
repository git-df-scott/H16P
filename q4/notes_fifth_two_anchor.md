# Strike 5: complete two-anchor geometry and the mixed center sign

2026-09-04. **PROVED:** the normalized two-simple-root region has exactly
two branches on every two-anchor fibre; its center sign is mixed. A
strictly positive endpoint value and a negative center value identify the
only branch left by the original sign-chain argument. This note does not
claim to exclude that final branch or construct four original zeros.

The inherited objects, ECT bounds, anchor chart, and center functional
are in [Q4_LOBE_REGION.md](../Q4_LOBE_REGION.md) and
[Q4_THEOREM_N.md](../Q4_THEOREM_N.md). No old shots, tangency search, or
corner asymptotic calculation is used here.

## 1. The exact affine chart for two anchors

Fix `0<r<s<1`. As in Theorem N, set

\[
 x=K_1/K_0,\qquad m=K_2/K_0.
\]

The graph `m(x)` is strictly convex and `x` is strictly increasing. Let
`b+S x` be its chord at the two anchors, and put

\[
 \boxed{V_{r,s}=K_2-bK_0-SK_1.} \tag{T1}
\]

Its coefficient of `K2` is one. Strict convexity gives two simple roots
`r,s`, signs `+,-,+`, and `V(1)>0`. The quotient `V/t^2` has positive
center value `(1/6-b)/2`; here `b<1/6` follows by extrapolating the strict
chord to `x=0`.

Define `B=B_{r,s}=H_{r,s,1}` by the three conditions

\[
 B(r)=B(s)=B(1)=0,\qquad [K_3]B=1.
\]

This is an exact nonsingular interpolation problem. Its three coefficient
rows, after division by `K0`, are `(1,x,m)` at `r,s,1`; strict convexity
makes their determinant nonzero, including the endpoint row. Equivalently,
`B` is the limit of the inherited three-anchor primitive `H_{r,s,u}` as
`u` increases to one. Continuity of the nonsingular interpolation matrix
proves coefficient convergence; no singular asymptotic is needed.

The usual anchored Rolle bound includes an endpoint root at one, counted
once: a primitive continuous on `[0,1]`, differentiable in `(0,1)`, and
vanishing at zero and at prescribed later points produces the intervening
zeros of `H'=tFq` by Rolle's theorem. Thus `B` has exactly the two simple
interior roots `r,s` and the endpoint root at one, and no extra contact.
Its interior signs are `+,-,+`, inherited by the coefficient limit and
the exclusion of any extra or multiple root. Its center coefficient is
strictly positive: the first-lobe comparison with any `H_{r,s,u}`, `u<1`,
and passage to the quotient by `t^2` give
`q_B(0)>=q_{r,s,u}(0)>0`.

Write

\[
 B=(A_B-1)K_0+B_BK_1-\eta_BK_2+K_3.
\]

Every primitive with the fixed anchors and `K3` coefficient one is uniquely

\[
 \boxed{H_\lambda=B+\lambda V,\qquad
 \lambda=\frac{H_\lambda(1)}{V(1)}.} \tag{T2}
\]

Indeed the difference of two such primitives belongs to
`span{K0,K1,K2}` and has the two anchors; its one-dimensional kernel is
spanned by `V`. Explicitly,

\[
 A=A_B-\lambda b,\qquad B_{\rm coeff}=B_B-\lambda S,
 \qquad\eta=\eta_B-\lambda.
\]

The notation `B` for the baseline primitive is distinguished here from
its scalar coefficient `B_B` or `B_coeff`.

## 2. A strictly decreasing residual ratio gives the complete classification

The quotient

\[
 R(t)=B(t)/V(t)
\]

has removable singularities at the common simple roots `r,s`. It extends
analytically across both, is strictly positive in `(0,1)`, and has

\[
 R(1)=0,\qquad
 R(0)=\frac{q_B(0)}{1/6-b}>0. \tag{T3}
\]

**Strict derivative proof.** If `R'(t0)=0` away from the anchors, then
`B-R(t0)V` has a double zero at `t0` and the two anchor zeros: at least
four interior zeros counted with multiplicity. If `t0` is one of the
anchors, the analytic extension and `R'(t0)=0` give a zero of order at
least three there, in addition to the other anchor. Both contradict the
inherited primitive bound three. Hence `R'` never vanishes in `(0,1)`.
Analyticity through the anchors makes its sign constant, and positivity
together with `R(1)=0` determines it:

\[
 \boxed{R'(t)<0\quad(0<t<1).} \tag{T4}
\]

Thus `H_lambda=V(R+lambda)` has exactly the following possibilities.

| Fibre parameter | Interior primitive roots and orientation |
|---|---|
| `lambda>0` | Exactly `r,s`, both simple, signs `+,-,+`; `H(1)>0`. |
| `lambda=0` | Exactly `r,s`, both simple, signs `+,-,+`; additional endpoint root `H(1)=0`. |
| `-R(0)<lambda<0`, except the next two values | Three simple interior roots, hence a point of the strict lobe region. |
| `lambda=-R(r)` | An ordinary double primitive root at `r`, and the simple root `s`. |
| `lambda=-R(s)` | The simple root `r`, and an ordinary double primitive root at `s`. |
| `lambda=-R(0)` | Exactly the simple interior roots `r,s`, signs `-,+,-`; the center quotient vanishes. |
| `lambda<-R(0)` | Exactly the simple roots `r,s`, signs `-,+,-`; `H(1)<0`. |

The two double roots are exactly quadratic, since `V` has a simple zero
and `R'` is strictly negative there. There is no third interior root at
the center boundary. In fact the center primitive is exactly cubic:
`q(0)=0`, and a higher center order would give at least four zeros of `q`
counted with multiplicity on `[0,s]`, using its two interior Rolle zeros
and the center multiplicity. The inherited auxiliary ECT property extends
to zero because its Wronskians there are strict Stieltjes moment
determinants. Thus `q` has a simple center zero.

The word “two-root” must distinguish the two simple-root branches from
the two double-contact levels: the latter have two distinct roots but
total interior multiplicity three. They lie in the closure of the strict
lobe region and have only one primitive sign change.

## 3. Strong exact bounds for the center functional

For any full primitive

\[
 U=\alpha K_0+\beta K_1+\gamma K_2+\delta K_3,
\]

linearity of the inherited reconstruction gives

\[
 \boxed{\mathscr Y(U)=\frac9{3080}
 \left(\alpha+X\beta+\frac{11}{6}\gamma+
       \frac{204}{221}\delta\right),\qquad X=\frac{144}{221}.} \tag{T5}
\]

Here

\[
 0<X<x(1)=\frac{6289}{9061},\quad
 m(0)=\frac16,\quad m(1)=\frac{2431}{9061}.
\]

For a strictly convex graph on a closed interval, every chord with two
interior endpoints, evaluated at any fixed interior abscissa, lies
strictly below the chord joining the full endpoints. If the evaluation
point lies between the chord endpoints this follows by interpolation of
the strict endpoint inequalities; if it lies outside them, extrapolation
of the chord lies below the graph itself. Applied to (T1) at `X`, this gives

\[
 b+SX<\frac16+
 \left(\frac{2431}{9061}-\frac16\right)
 \frac{144/221}{6289/9061}=\frac{9889}{37734}.
\]

Consequently the old qualitative (N3) bound strengthens to

\[
 \boxed{v:=\mathscr Y(V)>
 \frac9{3080}\left(\frac{11}{6}-\frac{9889}{37734}\right)
 =\frac{231}{50312}>0.} \tag{T6}
\]

A useful upper bound also follows from convexity. Every such chord at
`X` lies above the smaller of the two endpoint tangents evaluated at `X`.
For example, to the left of a chord its line lies above the tangent at
the right chord endpoint; to the right use the left tangent; between
them use the graph. The values of all intermediate tangents at `X` have
their minimum at an endpoint, since their derivative with respect to the
tangency abscissa has sign `X-x`. Here the two endpoint tangent values are

\[
 \frac16+\frac{25}{432}X=\frac{271}{1326},\qquad
 \frac{2431}{9061}+
 \frac{1105}{462}\left(X-\frac{6289}{9061}\right)=\frac16.
\]

The strict interior chord therefore gives

\[
 \boxed{\frac{231}{50312}<v<\frac3{616}.} \tag{T7}
\]

Only the lower bound is needed for the sign classification.

## 4. The zero-center endpoint still has positive eta

Let `n=K3/K0`, the weighted average of `tM(t)`. The same moment-slope
proof as in Theorem N shows strict convexity of `n` as a function of `x`:
its underlying function has
`(tM)''=2M'+tM''>0`. Its exact endpoints are

\[
 n(0)=0,\qquad n(1)=\frac{1819}{9061}.
\]

The primitive

\[
 E=B+\eta_BV
\]

has the two anchors, coefficient of `K3` one, and coefficient of `K2`
zero. It is therefore `K3` minus the chord of `n(x)` times `K0`.
Applying the same full-endpoint chord bound at `X` yields

\[
 \boxed{\mathscr Y(E)>
 \frac9{3080}\left[
 \frac{204}{221}-\frac{1819}{6289}\frac{144}{221}\right]
 =\frac{27}{12578}>0.} \tag{T8}
\]

The endpoint baseline belongs to the coefficient closure of the lobe
region. Theorem N's center comparison gives
`Y_B=mathscr Y(B)<=-3/1232<0`. Define

\[
 \lambda_c=-Y_B/v>0,\qquad C=B+\lambda_cV.
\]

Since (T8) says `Y_B+eta_B v>0`,

\[
 \boxed{0<\lambda_c<\eta_B,\qquad
 \eta_C=\eta_B-\lambda_c>0.} \tag{T9}
\]

In particular the center-zero endpoint is not a free homogeneous
solution: its other center datum has a fixed strict sign,

\[
 Y_C(0)=0,\qquad P_C(0)=Y_C'(0)=-\eta_C/192<0.
\]

On the first lobe `C>0`, so `P_C'=-Omega C<0`. Hence

\[
 \boxed{P_C(t)<0,\quad Z_C(t)<0\qquad(0<t\le r).} \tag{T10}
\]

This handles the zero-center value in the orientation argument; it would
be incorrect merely to treat a zero initial value as positive.

## 5. Center sign, exclusions, and the bounded remaining fibre

The exact center transport is

\[
 \boxed{Y_0(\lambda)=Y_B+\lambda v,
 \qquad\eta(\lambda)=\eta_B-\lambda.} \tag{T11}
\]

The negative branch `lambda<=-R(0)` has negative center value and negative
initial primitive sign. The positive branch has negative center value for
`0<=lambda<lambda_c`, zero center value at `lambda_c`, and positive center
value for `lambda>lambda_c`. Thus **the normalized two-root region has
MIXED center sign**, with an exact single transition on each positive
fibre.

The two-root sign-chain argument is recorded independently in the
Strike-5 sign-chain note: four original zeros require opposite center and
initial-primitive signs, three crossings of `P`, and the succeeding
alternating extremal heights. Matching signs therefore exclude the whole
negative branch and the part `lambda>lambda_c`. At `lambda=lambda_c`,
(T10) precludes the needed first P crossing. The case `lambda=0` is the
fixed-anchor endpoint limit of Theorem N, with its strict first-root tail
bound; its complete original sign count is treated in that companion
argument. The original `K3=0` projective chart has a two-root primitive
proportional to `V`, whose center and initial primitive signs agree by
(T6), and is excluded by the same orientation count.

The only unresolved normalized branch is consequently

\[
 \boxed{0<\lambda<\lambda_c<\eta_B,\qquad
 H_\lambda(1)>0,\quad H_\lambda:\ +,-,+,\quad Y_0<0.} \tag{T12}
\]

This is a bounded interval on every exact two-anchor fibre. It is not a
uniform positive separation from the endpoint or from `lambda_c`.

For fixed lift `a` and time `t`, all reconstructed quantities are affine
in `lambda`. In particular `Phi_lambda(r)` has positive slope and
`Phi_C(r)=int_0^r W_a C>0`. Theorem N cannot be extended as negativity of
`Phi(r)` over this whole interval. Likewise `Phi_B<0` does not assert
`Z_B<0` while `P_B>0`, because `Z=Phi+P Rcal`.

One further useful necessary reduction follows directly from (T9).
Write `c_a=3(1+a)/2+y_a'(0)`, so

\[
 P_\lambda(0)=-c_aY_0(\lambda)-\eta(\lambda)/192.
\]

It is affine, and at `lambda_c` is strictly negative. If an intermediate
point has the necessary `P_lambda(0)>0`, then `P_B(0)>0`, its fibre slope
is negative, and its admissible range ends strictly before `lambda_c`.
This follows from interpolation between the two endpoint values, without
any assumption that the slope has the same sign on every lift and fibre.

## 6. A further monotonicity of normalized baseline solutions

This section concerns only the endpoint baseline `B=H_{r,s,1}` and a
fixed lift. Put `d=-Y_B>0`. Increasing either finite anchor gives a
variation

\[
 \partial_jB=\gamma_j V_j,\qquad\gamma_j>0,
\]

where `Vj` is the `K2`-normalized chord primitive at the other two anchors,
one of which is the endpoint one. The usual cardinal-interpolation proof
still applies: its normalized moment matrix is strictly convex including
that endpoint. On the old first lobe `Vj>0`, and
`partial_j d=-gamma_j mathscr Y(Vj)<0`. Hence

\[
 \partial_j(B/d)>0.
\]

There is also a strict center-slope inequality. Write `v_j=mathscr Y(Vj)`.
The primitive `B+eta_B Vj` has `eta=0`, coefficient of `K3` one, and the
other two anchors. The positive bound (T8) continues to hold when its
second chord endpoint is one: the first remains strictly interior, and
the chord at `X` is strictly below the full-endpoint chord. Therefore

\[
 \eta_Bv_j-d>0,\qquad
 \partial_j(\eta_B/d)
 =\gamma_j(\eta_Bv_j-d)/d^2>0. \tag{T13}
\]

The normalized reconstruction has initial data

\[
 Z_B(0)/d=-1,\qquad P_B(0)/d=c_a-\eta_B/(192d).
\]

Its normalized initial momentum strictly decreases with either anchor,
while its positive first-lobe forcing `B/d` strictly increases. Integrate
`(P_B/d)'=-Omega B/d` and then `(Z_B/d)'=(P_B/d)/(p y^2)` to obtain

\[
 \boxed{\partial_j(P_B(t)/d)<0,\qquad
 \partial_j(Z_B(t)/d)<0\quad(0<t\le r),} \tag{T14}
\]

where the comparison uses a fixed time in the old first lobe. In
particular, at a fixed first anchor `r`, the normalized value `Z_B(r)/d`
decreases as `s` increases. An upper bound for that baseline value is
therefore reduced to the confluent endpoint interpolation `H_{r,r,1}`
as `s` decreases to `r`. This is a finite, exact boundary limit of the
anchor matrix. It supplies no sign for the boundary value, and it does
not assert monotonicity of the moving-time quantity `Z_B(r)/d` when `r`
itself varies. Nor does it by itself settle the bilinear Green determinant
at the first anchor. Those distinctions prevent promoting an incomplete
first-peak comparison to an exclusion theorem.

## 7. Exact replay and current scope

[check_fifth_two_anchor.py](check_fifth_two_anchor.py) verifies the rational
coefficients in (T5), the endpoint slope, (T6)–(T8), and the tangent values.
It sets one-thread environment variables, reduced priority, and a ten-second
CPU ceiling. It uses standard-library rational arithmetic only; no sampling
or quadrature result is used as proof. The analytic chord, ratio, and
zero-count arguments above are separate from this arithmetic replay.
The replay passed with approximately `0.000064` process CPU seconds;
it reported the exact bounds `231/50312`, `3/616`, and `27/12578`.

The remaining analytic target is the first original Green maximum on
(T12), particularly the baseline case `P_B(r)>0`. No four-original-zero
example or global three-zero theorem is asserted here.
