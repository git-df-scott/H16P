# Picard–Fuchs and universal auxiliary geometry notes

Independent derivation, 2026-09-04. These notes use exactly the positive-area
integrals and coefficient convention in `q4_integrals.py`. They contain no
five-zero example and no claim about the original integral being Chebyshev.

## 1. Inherited identities and the coefficient convention

Write `k=κ`, `d=κ−1`, and let `J1=I00'(h(s))`, `J2=I11'(h(s))`.
The two-period Picard–Fuchs equation is Zhao equation (14):

\[
6(s-1)(s-k)\binom{J_1'}{J_2'}=
\begin{pmatrix}1-s&d\\1-s&s-1\end{pmatrix}\binom{J_1}{J_2}.
\]

The source is [Zhao, arXiv:1011.2253](https://arxiv.org/html/1011.2253).
Its statements about the four-dimensional original family and the zero chain
are inherited theorems, rather than results proved here. In particular,
`Z(I)≤Z(G)≤Z(ℱ)+2` and `Z(ℱ)≤Z(g)`, with `g=ℱ'/J1`.

Before using the source's alpha/beta map, the original coefficients must be
changed to `m2=−2μ2/3−2d μ3/(3k)` and `m3=−2μ3/(3k)`.
The repository implements this change. An exact determinant computation gives

\[
\det\frac{\partial(\alpha_1,\alpha_2,\beta_0,\beta_1)}
 {\partial(\mu_1,\mu_2,\mu_3,\mu_4)}
=\frac{17843617792000(k-1)^2}{6561k^4}>0.
\]

Thus no projective directions are lost when these four coordinates are used.
The assertion uses the map as explicitly implemented and the source's
relation `α0=β0−kβ1−kα1−k²α2`.

## 2. Exact universal period reduction

**Derived identity.** Put

\[
t=\frac{k-s}{d}\in(0,1),\qquad C=J_1(k)=\frac{\pi}{\sqrt d}.
\]

The center is `t=0`; the homoclinic endpoint is `t=1`. The constant `C`
follows from the Hessian determinant `det Hess H(1,1)=4d` and the small-oval
area formula. The period equations become

\[
6t(1-t)\dot J_1=-(1-t)J_1+J_2,
\qquad 6t\dot J_2=J_2-J_1.
\]

Eliminating `J2` gives

\[
t(1-t)F''+(1-2t)F'-\frac5{36}F=0,\quad F(0)=1,
\]

and regularity at the center selects

\[
\frac{J_1}{C}=F(t)={}_2F_1(1/6,5/6;1;t),\qquad
\frac{J_2}{C}=(1-t)(F+6tF')={}_2F_1(-1/6,1/6;1;t).
\]

The Taylor coefficient recurrence is
`f[n+1]/f[n]=(n+1/6)(n+5/6)/(n+1)²`; thus these formulas also follow
without any special-function identification. This rank-two system is the
smallest useful system for the auxiliary zero analysis. Its scalar equation
has regular singularities `0,1,∞`, with exponents respectively `(0,0)`,
`(0,0)`, and `(1/6,5/6)` in the usual infinity convention.

In particular the ratio is a single universal function, independent of `k`:

\[
w(t)=(1-t)\left(1+6t\frac{F'}F\right),\qquad
6t(1-t)w'=-(1-t)+2(1-t)w-w^2.
\]

An elementary proof that `F` never vanishes in the slit plane is given below.

## 3. Endpoint expansions proved from this system

At the center,

\[
F=1+\frac5{36}t+\frac{385}{5184}t^2+
\frac{85085}{1679616}t^3+O(t^4),
\]
\[
w=1-\frac t6-\frac{25t^2}{432}-\frac{775t^3}{23328}
-\frac{305675t^4}{13436928}+O(t^5).
\]

Define `M=(1−w)/t`, with its removable value at zero. Then

\[
M=\frac16+\frac{25}{432}t+\frac{775}{23328}t^2+
\frac{305675}{13436928}t^3+O(t^4).
\]

For the homoclinic endpoint put `z=1−t` and `L=log(432/z)`. The regular
singular expansion, obtained by the Frobenius recurrence and matching the
Euler integral, is

\[
F(t)=\frac1{2\pi}\left[L+\frac5{36}z(L-26/5)
+O(z^2|\log z|)\right],
\]
\[
\frac{J_2}{C}=\frac3\pi+\frac{z(L-5)}{12\pi}
+O(z^2|\log z|),\qquad w(t)\sim\frac6L,
\qquad M(t)\longrightarrow1.
\]

For example this identifies the constant left unspecified in the leading
`J1` expansion in the source:

\[
J_1(s)=\frac1{2\sqrt d}\log\frac{432d}{s-1}
+O((s-1)|\log(s-1)|).
\]

The first Taylor coefficient vector of the **original four area integrals**
at `s=k`, including their normalization, is

\[
C\left(\frac2{9k},-\frac1{3k},-\frac1{3k},\frac{2d}{3k}\right).
\]

One must not copy Zhao equation (33) literally as a Taylor series: the
quantities in (34) are derivatives, hence the Taylor expansion needs `1/n!`,
and its common period normalization is suppressed. For instance the
quadratic Taylor coefficient of `hI00` is
`C(13k−18)/(324k²d)`, `C` times one half of the normalized second derivative displayed in (34).
The universal period expansions alone are not expansions of the entire
original four-dimensional family; its third-kind part still needs the
additional equations described in section 7.

## 4. Positive Stieltjes representation: a derived theorem

For `z∈C\[1,∞)`,

\[
\boxed{M(z)=\int_0^1\frac{\rho(u)}{1-zu}\,du},\qquad
\boxed{\rho(u)=\frac3{2\pi^2|F(1/u+i0)|^2}>0}\quad(0<u<1).
\]

Here all powers and logarithms use the principal branch, and the positive
side of the cut means approach from the upper half-plane. This is an exact
representation, not a numerical fit or a conjecture.

**Proof.** The Euler integral is

\[
F(z)=\frac1{2\pi}\int_0^1
u^{-1/6}(1-u)^{-5/6}(1-zu)^{-1/6}\,du.
\]

In the upper half-plane each integrand has argument strictly between zero
and `π/6`; thus `F` has positive real part and cannot vanish. Conjugation
handles the lower half-plane, and the integrand is positive for real `z<1`.
Consequently `M=1−6(1−z)F'/F` is analytic on the slit plane.

For `x>1`, write `F(x+i0)=U(x)+iV(x)`. The real Gauss equation and Abel's
identity imply

\[
UV'-U'V=\frac{c}{x(1-x)}.
\]

The endpoint expansion on the upper bank is
`F(x+i0)=[−log(x−1)+log432+iπ]/(2π)+O((x−1)|log(x−1)|)`.
Its derivative therefore gives `c=−1/(4π)`. It follows that

\[
\operatorname{Im}M(x+i0)=
-6(1-x)\frac{UV'-U'V}{|F|^2}
=\frac3{2\pi x|F(x+i0)|^2}>0.
\]

The regular singular expansion at infinity is
`F(z)=A(−z)^(−1/6)(1+O(z^(−2/3)))`, with `A≠0`; differentiating this
convergent Frobenius expansion shows `M(z)=O(z^(−2/3))`.
At `z=1`, `M(z)` is bounded and tends to one. Thus the large circle and
the small circle about one in the slit-plane Cauchy contour contribute zero.
The two banks are conjugate, giving

\[
M(z)=\frac1\pi\int_1^\infty
\frac{\operatorname{Im}M(x+i0)}{x-z}\,dx.
\]

Substitution `u=1/x` proves the stated formula. Convergence is also explicit:
`ρ(u)` is a positive constant times `u^(−1/3)` at zero and a positive
constant times `log(1−u)^(−2)` at one. These estimates justify differentiation
under the integral for `0<t<1`. The boundary value `M(1)=1` is finite.

This proof uses only the just-derived Gauss equation, its Euler integral,
and elementary complex analysis. Its statement concerns the auxiliary
period ratio, not the original Abelian integral.

## 5. Universal auxiliary curve and exact inflection theorem

In the chart `β1=1`, set

\[
b=\frac{k-\beta_0}{d},\quad A=-(\alpha_1+2k\alpha_2),
\quad B=d\alpha_2.
\]

Direct substitution and the relation for `α0` give

\[
g(k-dt)=d\{At+Bt^2+b(w-1)-tw\}
=dt\{A+Bt+q_b(t)\},
\qquad q_b(t)=(t-b)M(t)-1.
\]

Conversely `α2=B/d`, `α1=−A−2kB/d`, `β0=k−db`; thus this is an
invertible coordinate change on the entire `β1≠0` projective chart.
All the auxiliary zero geometry is therefore universal: it can be studied
once in `(A,B,b)` on `(0,1)`, independent of the center parameter `k`.

The Stieltjes formula gives

\[
q_b''(t)=2\int_0^1\frac{u(1-bu)\rho(u)}{(1-tu)^3}\,du.
\]

Let

\[
b_*(t)=\frac{\int_0^1u\rho(u)(1-tu)^{-3}\,du}
 {\int_0^1u^2\rho(u)(1-tu)^{-3}\,du}.
\]

**Derived theorem.** The function `b_*` is strictly decreasing, with
`b_*(0)=54/31` and `lim[t→1] b_*(t)=1`. Hence:

* for `b≤1`, `q_b''>0` everywhere;
* for `b≥54/31`, `q_b''<0` everywhere in the open interval;
* for `1<b<54/31`, `q_b''` changes from positive to negative at exactly
  one point, and the change is simple.

**Proof.** Give `(0,1)` the probability measure proportional to
`uρ(u)/(1-tu)^3 du`. Then `b_*=1/E_t[u]`. Differentiation gives
`dE_t[u]/dt=Cov_t(u,3u/(1-tu))>0`, because both functions are strictly
increasing and the density has support throughout `(0,1)`. The covariance
identity can itself be written as a strictly positive double integral of
the product of their differences. The Taylor moments in section 3 give
`b_*(0)=(25/432)/(775/23328)=54/31`. As `t→1`, the normalizing integral
diverges and its mass concentrates at `u=1`; this follows either by a direct
comparison using the density estimate or by splitting at `u=1−ε`.
Thus `E_t[u]→1`. Finally `q_b''=2V(t)(b_*(t)−b)`, where `V(t)>0`, so its
unique zero is simple.

Consequently the correct three-zero strip in original coordinates is

\[
\frac{54-23k}{31}<\beta_0<1,
\]

which is nonempty for **every** `k>1`. The opposite sign printed in Zhao's
Theorem 14, and inherited by the initial repository preparation, cannot be
correct. It is inconsistent with the paper's Proposition 17 and with this
independent universal reduction. In particular no upper bound `k<85/23`
follows from this argument.

For each `b` in the open strip, let `t0` be the unique inflection. Choose
`A,B` so that `A+Bt+q_b(t)` vanishes to order three at `t0`. Its third
derivative there is strictly negative. Replacing this expression by itself
plus `ε(t−t0)` for sufficiently small `ε>0` creates three distinct simple
zeros near `t0`. This follows by the analytic implicit function theorem
after putting `t−t0=sqrt(ε)y`; the limiting polynomial is
`y+c y³`, `c<0`, which has three simple real zeros.

Thus three simple auxiliary `g` zeros are possible for every fixed `k>1`.
This is not a five-zero bifurcation mechanism for `I`: the intervening
integration and second-order lift need not realize the maximal Rolle losses.

## 6. A Chebyshev statement that really is proved

Let `v=(w−1)/t=−M`. The normalized four-function space is

\[
\mathcal V=\operatorname{span}\{1,t,w,v\}.
\]

Its ordered Wronskians satisfy `W1=W2=1`,
`W3=w''=−2M'−tM''<0`, and

\[
W_4=w''v'''-w'''v''=-(v'')^2b_*'(t)>0.
\]

For this last identity, `v''=−M''<0` and
`w''/v''=t+2M'/M''=b_*`; differentiating the ratio proves the formula.

The standard Wronskian criterion therefore makes this ordered system an
extended complete Chebyshev system on `(0,1)`. In particular it has at most
three zeros counting multiplicity. The important distinction is that this
is the space of `g/t`, not the space of `I` divided by its center factor.
No original-family Chebyshev theorem follows merely by naming this one.

### The infinitesimal auxiliary cusp is excluded as a five-zero route

The root agent identified the following further deduction, independently
checked against the formulas here. Fix any `b∈(1,54/31)` and its inflection
`t*`. Let `q0=A+Bt+q_b(t)` have its triple zero at `t*`, and put
`qλ=q0+λ(t−t*)`, `λ>0`. Define the universal weighted primitive

\[
H_\lambda(t)=\int_0^t F(u)\,u\,q_\lambda(u)\,du.
\]

**Derived exclusion theorem.** For all sufficiently small positive `λ`,
the original Q4 integral corresponding to this auxiliary direction has
at most three interior zeros, for every `k>1`.

Indeed, `q0` has a triple zero and, by the just-proved ECT property, no
other interior zeros. Since `q0'''(t*)<0`, it is positive to the left
of `t*` and negative to the right. Therefore `H0(t*)>0`. The three
simple zeros `x1<x2<x3` of `qλ` tend to `t*`; ECT excludes any others.
Continuity then gives `Hλ(xi)>0` for all three critical points. The sign
pattern of `Hλ'` is `+,−,+,−`, so `Hλ` has no zero before `x3`, and at
most one after it. Under the exact coordinate change,

\[
\mathcal F(k-dt)=-Cd^2H_\lambda(t),
\]

because `ℱ(k)=0` and `ℱ'=J1g`. Thus `Z(ℱ)≤1` and the inherited
zero chain gives `Z(I)≤3`. This is uniform in the sense that the small
`λ` restriction depends only on the fixed universal cusp, not on `k`.
It is not an exclusion of finitely separated three-zero auxiliary
configurations, nor of all possible original-family bifurcations.

## 7. What the period reduction does and does not say about scalar ODEs

For the original integral use `L1=hD_h−1`. Its image `G=L1I` is a
coefficient-dependent combination of the four periods
`(I00', I11', I−1,0', I−1,1')`. Their first-order PF system is closed;
the first two entries form the universal rank-two subsystem. Its other two
entries account for the third-kind part. This core plus
`I(h)=h∫[−2/3,h] G(ξ)/ξ² dξ` is an exact lift back to the original family.

The source's operator is

\[
L_2=5kh-(9kh^2-8)D_h+h(9kh^2-4)D_h^2.
\]

For each fixed `μ`, `R=L2L1I` is a rational-coefficient combination of
`J1,J2`. If `r` is its row in their rational first-order system `J'=AJ`,
put `r1=r'+rA`, `r2=r1'+r1A`. Wherever `det(r,r1)≠0`, expressing `r2`
in that pair constructs a rational second-order operator `M2` with
`M2R=0`. Consequently there is an order-at-most-five scalar equation
`M2∘L2∘L1(I)=0` for each such fixed coefficient choice. This is the
factorization described in
[Gavrilov–Iliev](https://arxiv.org/html/0811.4602).

Crucially `M2` generally depends on `μ`; this is not a common order-five
operator for all four basis functions. A common order-four equation with
analytic coefficients can always be constructed where their full Wronskian
does not vanish, but using that tautology to claim disconjugacy would be
circular. A globally regular common scalar operator and the sign of the
original-family Wronskian have not been established here.

In `s`, the homogeneous `L2` becomes
`s(1−s)D_s²−(1/2)D_s−5/36`. Its singularities are `0,1,∞`, with
exponents `(0,3/2)`, `(0,1/2)`, and `(−1/6,−5/6)`. The original annulus
is `1<s<k`, so this does not introduce an internal singularity. The
coefficient-dependent operator `M2` can additionally acquire apparent
singularities at zeros of `det(r,r1)`; their absence is not assumed.

## 8. Verification and the remaining mathematical task

`q4_structure_checks.py` independently checks the exact coefficient-map
determinant, the universal change of variables, period coefficient
recurrences, ratio moments, corrected threshold, and three specified
hypergeometric identities at `t=1/4,1/2,3/4`. It sets numerical thread
counts to one and enforces a ten-second CPU ceiling. One replay completed
in approximately half a second; no optimization or parameter scan was run.
The Stieltjes and inflection conclusions are analytic proofs above, not
inferences from these floating-point checks.

The strongest reduction is now precise: the auxiliary period and
three-intersection geometry has no `k` dependence. The remaining question
is whether a universal auxiliary direction survives the weighted primitive
and the `k`-dependent third-kind lift with enough zero alternation to give
five simple zeros of the original `I`. No such direction is certified in
these notes.
