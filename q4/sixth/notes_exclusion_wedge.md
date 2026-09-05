# Strike 6: an exact exclusion region for the dangerous two-root family

2026-09-05. Base: `3b94f34`. These are new deductions from the audited
reconstruction, Theorem N, and the Strike-5 two-anchor theorem. The
arithmetic replay is [check_exact.py](check_exact.py). No sampled sign
enters the proof. This note does not settle the remaining determinant.

## 1. Statement and normalization

Use the notation of [Q4_TWO_ROOT_REDUCTION.md](../../Q4_TWO_ROOT_REDUCTION.md).
Every potential four-distinct-interior-zero integral can be normalized as

\[
H=B+\lambda V,\qquad B=H_{r,s,1},\qquad
0<r<s<1,\quad 0<\lambda<\lambda_c=-Y_B/v,
\]

where `v=Y0(V)>0`, `Y_B<0`, `eta=eta_B-lambda>0`, and `H>0` on `(0,r)`.
Set

\[
d=-Y_0(H)>0,\qquad q=\frac{\eta}{192d},\qquad a=1-1/\kappa.
\]

Here `q` is a scalar center-data ratio, **not** the earlier auxiliary
function also called `q(t)`.

**New necessary conditions.** A four-zero candidate must satisfy

\[
\boxed{r>1-\left(\frac7{22}\right)^{3/2}
       =0.8205212489\ldots,\qquad q>\frac{19}{10},\qquad
       \kappa>\kappa_*.} \tag{W1}
\]

Define `u_*` as the unique positive root of

\[
1024u^5+2816u^4+3936u^3+2584u^2+440u-135=0,
\]

and `kappa_*=(1+u_*)(1+4u_*)^2`. Exact rational isolation gives

\[
2.899241080973277432530648<\kappa_*
<2.899241080974989225893177. \tag{W2}
\]

The complementary closed parameter regions are excluded, including
equality in either bound in (W1). These are **necessary** restrictions;
points satisfying them are not four-zero candidates without the remaining
determinant and height tests.

## 2. A global lower bound on the normalized center slope

Write `d_B=-Y_B`. Strike 5 proves that `eta_B/d_B` strictly increases with
each anchor, including finite confluent limits. Along the dangerous fibre,

\[
\frac{d}{d\lambda}\frac{\eta_B-\lambda}{d_B-\lambda v}
=\frac{\eta_Bv-d_B}{(d_B-\lambda v)^2}>0. \tag{W3}
\]

The strict numerator is the positive center functional of
`B+eta_B V`, already proved in Strike 5. Thus the infimum over all
dangerous fibres is attained as a boundary limit: first `lambda -> 0`,
then `r,s -> 0`.

For that limit divide `H` by `t^2` before taking confluent anchor rows.
The limiting interpolation conditions are

\[
A-1-\eta/6=0,\quad B_{\rm coeff}+1/6-25\eta/432=0,
\quad9061A+6289B_{\rm coeff}-2431\eta-7242=0.
\]

Their exact solution is

\[
(A,B_{\rm coeff},\eta)
=\left(\frac{11843}{9623},-\frac{833}{9623},\frac{13320}{9623}\right),
\qquad Y_0=-\frac{81}{19246}.
\]

The limiting interpolation matrix is nonsingular (the displayed system
has nonzero determinant), so analytic divided differences justify the
coefficient limit. It follows that every interior dangerous point obeys

\[
\boxed{q>\frac{185}{108}.} \tag{W4}
\]

Strictness is retained by comparison with any intermediate positive
anchor pair before taking the boundary limit.

## 3. A supersolution that bounds every lift at once

The original reconstructed derivative satisfies

\[
(1-at)(1-t)Y''-\frac{1-a}{2}Y'+\frac{5a}{36}Y
=-\frac{H}{1152t^2(1-t)}.
\]

Put `T=Y/[d(1-at)^(3/2)]`. Exact substitution gives

\[
L_aT=-\frac{H}{1152d\,t^2(1-t)(1-at)^{3/2}}<0\quad(0<t<r),
\]

\[
L_a=(1-at)(1-t)D^2-\frac{1+5a-6at}{2}D+\frac{8a}{9},
\qquad T(0)=-1,\quad T'(0)=\frac32-q. \tag{W5}
\]

As proved in Theorem N, this operator has a positive causal Green
function on every compact subinterval of `[0,1)`. For example its positive
homogeneous solution is `y_a/(1-at)^(3/2)`; reduction of order proves the
Green sign. Its homogeneous solution with initial data `(0,1)` is also
positive. Consequently a positive initial derivative difference and a
nonnegative inhomogeneity imply a positive solution difference.

For a comparison number `q0`, define

\[
A_0=\frac{13}{4}-\frac32q_0,\quad B_0=A_0+1,\quad
S_{q_0}(t)=A_0(1-t)^{-4/3}-B_0(1-t)^{-2/3}.
\]

It has `S(0)=-1`, `S'(0)=3/2-q0`, and its exact residual is

\[
\boxed{L_aS_{q_0}=
\frac{1-a}{9(1-t)^{7/3}}
\left[22A_0-7B_0(1-t)^{2/3}\right].} \tag{W6}
\]

For either `q0=185/108` or `q0=167/90`, `A0,B0>0` and
`22A0-7B0>=0`. The bracket is positive for `t>0`. If `q>q0`,
(W5), (W6), and the positive Green function therefore imply

\[
\boxed{T(t)<S_{q_0}(t)\quad(0<t\le r).} \tag{W7}
\]

This comparison uses neither monotonicity of `a` nor an unproved bound on
the loop homogeneous solution. The strict initial derivative difference
already gives strictness, even at an endpoint where `H(r)=0`.

The comparison function is nonpositive up to its unique zero

\[
t(q_0)=1-\left(\frac{A_0}{B_0}\right)^{3/2}.
\]

If `r<=t(q0)`, then `Y<0` throughout the first primitive lobe, including
its endpoint. The first strictly positive Green maximum required by
Theorem T is impossible. Thus any survivor has `r>t(q0)`.

## 4. Bootstrap with two exact, one-point moment certificates

First use (W4): at `q0=185/108`, `A0=49/72`, `B0=121/72`, so

\[
r>1-(7/11)^3=988/1331>7/10.
\]

Anchor monotonicity then compares the baseline to the finite confluent
primitive `H_{7/10,7/10,1}`. The exact rational interval certificate gives

\[
\frac{\eta}{-192Y_0}
\bigg|_{H_{7/10,7/10,1}}
\in(1.85859672526316062756,
     1.85859672526316062770)
>167/90.
\]

Both actual anchors are larger; (W3) increases the ratio further along
the dangerous fibre. Therefore `q>167/90`. Apply (W7) again: now
`A0=7/15`, `B0=22/15`, yielding

\[
r>1-(7/22)^{3/2}>4/5.
\]

The second exact certificate gives

\[
\frac{\eta}{-192Y_0}
\bigg|_{H_{4/5,4/5,1}}
\in(1.90188147601594,1.90188147602308)>19/10.
\]

The same monotonicity proves `q>19/10` for every remaining candidate.
We do **not** substitute `19/10` in (W7): its residual is not nonnegative
near zero. The second certificate is used only for the momentum bound
in the next section.

The certificates enclose the four positive moment series and the
derivative row with `fractions.Fraction`. They include terms through
index 128, bound each positive tail by its first omitted term divided by
`1-t`, and solve the confluent interpolation system by interval Cramer
determinants. No value of pi needs to be approximated: the endpoint row
is the exact positive multiple `(9061,6289,-2431)` with right side `-1819`.
All division intervals avoid zero. Full outward-rounded intervals are in
[exact_checks.json](exact_checks.json).

## 5. The exact lift cutoff

The required positive initial momentum is

\[
P_0=d(C_a-q)>0,\qquad C_a=\frac32(1+a)+y_a'(0).
\]

Hence `C_a>19/10`. Use the exact elementary homogeneous solution from
the reconstruction. With

\[
u=\sinh^2\!\left(\frac13\operatorname{arsinh}\sqrt{\kappa-1}\right)>0,
\quad\kappa=(1+u)(1+4u)^2,
\]

its center slope simplifies to

\[
y_a'(0)=-\frac{5(4u+3)(2u+1)}{6(8u^2+10u+5)}.
\]

Direct rational algebra gives

\[
C_a-\frac{19}{10}=
\frac{1024u^5+2816u^4+3936u^3+2584u^2+440u-135}
{30(1+u)(1+4u)^2(8u^2+10u+5)}. \tag{W8}
\]

The denominator is positive; the numerator is negative at zero, tends
to positive infinity, and has strictly positive derivative for `u>=0`.
It has exactly one positive root. Since `kappa(u)` is strictly
increasing, (W8) proves the cutoff (W1)-(W2), including the excluded
equality case.

## 6. Scope

These inequalities remove whole parameter regions from the *four-interior-
zero, two-root* problem. They do not assert that the remaining region has
four zeros, that every determinant there is nonpositive, or that endpoint
cycles are bounded by an interior-zero theorem. The global H(2) problem
and the outside-lobe Q4 four-zero question remain unresolved.
