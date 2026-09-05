# Resonant joint strike: compact cycles consume the endpoint degeneracy

Date: 2026-09-05. Calculation base: `3b94f34`; integrated after `79001f7`.

**Result:** the proposed mechanism of **two compact cycles plus three
hemicycle cycles at the double-center base `(-1,1)` is incompatible**.
The argument below gives at most one hemicycle cycle along a sequence
carrying two compact cycles. This is a local compatibility result, not a
proof of `H(2)=4` or a sharp bound for pure endpoint cyclicity.

This is a derived analytic argument in the campaign, with symbolic and
numerical cross-checks; it has not received an independent mathematical
review. Numerical return signs below are explicitly not interval certificates.

## 1. Exact scope and source inputs

We use the five-parameter family

\[
 P=\frac{b-2}{4}+\epsilon _1x+(1-b)y+ax^2+\epsilon _2xy+by^2,
 \qquad Q=\epsilon _0-2xy,
 \qquad \mu_* =(-1,1,0,0,0).
\]

A **compact cycle** here tends to a periodic orbit or a center of the
base field in a bounded subset of one half-plane. An **endpoint cycle**
tends in the Poincare compactification to its upper or lower hemicycle.
We consider simultaneous limits in one sequence of parameter vectors;
we do not add cyclicities attained in different perturbations.

Two inputs from [Marín–Villadelprat](https://arxiv.org/html/2501.16924v1)
are used: Proposition 3.2 gives a parameter-uniform leading Dulac expansion
also at resonance; Lemma 3.3 identifies the upper center set. Their
Theorem 2.1 supplies the persistent-connection logarithmic invariant.
Theorems B and C do not supply a resonant sharp upper bound.

The compact normalization is consistent with the Bautin analysis of
[Françoise–Gavrilov](https://arxiv.org/html/2011.08316v1).
We calculate the necessary functions directly below. Their compact
cyclicity theorem does not bound cycles escaping to infinity.

## 2. The base is explicitly integrable

At the base, for `z=x+iy`,

\[
 \dot z=-z^2-\frac14,\qquad
 H(x,y)=\frac{x^2+y^2+1/4}{y}.
\]

The upper level `H=h>1` is the circle

\[
 x^2+(y-h/2)^2=(h^2-1)/4.
\]

Its center height is `c=h/2`, its radius is
`r=sqrt(h^2-1)/2`, and `sqrt(c^2-r^2)=1/2`.
For its interior disc, direct integration gives

\[
 I_1=\iint y^{-1}=\pi(h-1),\quad
 I_2=\iint y^{-2}=2\pi(h-1),\quad
 I_3=\iint y^{-3}=2\pi(h^2-1).
\]

The integrating factor is `y^-2`. On the base, the perturbation's weighted
divergence has its even-in-x part

\[
 \epsilon _1y^{-2}+\epsilon _2y^{-1}-2\epsilon _0y^{-3}.
\]

Thus the upper first normal Melnikov function, with the orientation fixed
by the forward upper half return, is

\[
 M_u(h)=2\pi(h-1)
 \left[\epsilon _1+\frac{\epsilon _2}{2}
                   -2\epsilon _0(h+1)\right].                 \tag{1}
\]

Reflection to the lower half-plane at the base replaces
`epsilon0 -> -epsilon0`, `epsilon2 -> -epsilon2`. In particular the two
linear combinations are

\[
 c_u=\epsilon _1+\epsilon _2/2,\qquad
 c_l=\epsilon _1-\epsilon _2/2.
\]

## 3. An exact first-normal joint Dulac calculation

Use the source's upper sections `(0,1/s)` and `(0,y)` with `0<s<2`.
The base lower intersection is `y=s/4` and

\[
 h=1/s+s/4,\qquad H_y(0,s/4)=1-4/s^2.
\]

Dividing (1) by this derivative gives the full first parameter derivative
of the **difference of the two half maps**, not merely its endpoint jet:

\[
 D^{[1]}_u(s)=\pi\epsilon _0(1-s^2/4)
  -\frac\pi2c_u\frac{s(2-s)}{2+s}.                            \tag{2}
\]

In positive reflected lower section coordinates,

\[
 D^{[1]}_l(s)=-\pi\epsilon _0(1-s^2/4)
  -\frac\pi2c_l\frac{s(2-s)}{2+s}.                            \tag{3}
\]

For example, the upper jet is

\[
 \pi\epsilon _0-\frac\pi2c_us+
 (\tfrac\pi2c_u-\tfrac\pi4\epsilon _0)s^2+O(s^3\|\epsilon\|).
\]

These are first parameter derivatives. They are not the full nonlinear
five-parameter return map. Neither equation is used with an uncontrolled
`O(epsilon^2)` error at an exponentially small section coordinate.

## 4. The missing compact generator is elementary, with a logarithm

Put

\[
 A=a+1,\ B=b-1,\ C=A+B=a+b,\
 c=c_u,\ g=C\epsilon _1.
\]

The upper center set is

\[
 \{\epsilon _0=c=\epsilon _1=0\}
 \ \cup\ \{\epsilon _0=c=C=0\}.
\]

Consequently every compact return displacement can be divided in the
three generators

\[
                 \epsilon _0,\quad c,\quad g.                \tag{4}
\]

Here no radical-ideal inference is needed: subtract the values at
`epsilon0=0`, then at `c=0`. The remainder vanishes at `epsilon1=0`
and at `C=0`, so two applications of the fundamental theorem of calculus
divide it by `epsilon1*C`. This construction is also useful at infinity.

The third base coefficient is the mixed derivative obtained by setting
`epsilon2=-2epsilon1`, `a=-1`, and varying `b`. For the reversible field
at `a=-1`, its orbit equation is

\[
 x^2=\mathcal H y-by^2-(1-b)y\log y+(b-2)/4.
\]

At `b=1`, the derivative of its right side at fixed `mathcal H` is

\[
 R_b=-y^2+y\log y+1/4.
\]

Varying the area integral, including the moving boundary, gives

\[
 J(h)=\int_0^\pi
  \frac{1-2y}{y^2}(-y^2+y\log y+1/4)\,d\theta,
 \quad y=h/2+\sqrt{h^2-1}\cos\theta/2.
\]

Evaluation yields

\[
          J(h)=2\pi\left[h-1-2\log\frac{h+1}{2}\right].      \tag{5}
\]

For checking the integration, the integrand expands to

\[
 -1+2y+\frac{\log y}{y}-2\log y+\frac1{4y^2}-\frac1{2y}.
\]

The two logarithmic integrals are
`integral log(y) = pi log((h+1)/4)` and
`integral log(y)/y = -2pi log(h+1)` over `[0,pi]`.
The nonlogarithmic part integrates to `2pi(h-1)`.
Variation in `a` gives the same coefficient: on `a+b=0`, `epsilon0=c=0`
the upper field remains a center. There is no independent pure
`epsilon1^2` term, because `a=-1,b=1,epsilon0=c=0` is itself an exact
center family.

For any sequence approaching the base, set

\[
 N=\max(|\epsilon _0|,|c|,|g|),\qquad
 (\epsilon _0,c,g)/N\longrightarrow(p,q,k)
\]

after taking a subsequence with `N>0`. The limiting compact displacement,
divided by its positive base factor `2pi(h-1)`, is

\[
 F(u)=q-4p(1+u)+k f(u),\qquad
 f(u)=1-\frac{\log(1+u)}u,\qquad u=\frac{h-1}{2}>0.          \tag{6}
\]

The functions in (6) are linearly independent, so normalization by `N`
cannot yield the identically zero function. Convergence is analytic on
compact transverse intervals. At the center the usual focal
normalization gives the same expansion with `f(0)=0`.

The integral representation

\[
 f(u)=\int_0^1\frac{ut}{1+ut}\,dt
\]

shows `f'>0` and `f''<0` on `u>=0`. Therefore:

* If `p=0`, a nonzero function (6) has at most one zero, with multiplicity.
* If `k=0`, the same assertion holds because (6) is affine.
* **Two compact cycles in the upper annulus require `p*k != 0`.**

This includes coalescing compact roots. For cycles shrinking to the
center, the expansion is

\[
 F(u)=(q-4p)+(-4p+k/2)u-(k/3)u^2+O(u^3).
\]

When `p=0`, a zero at the center has `q=0` and nonzero next coefficient
`k/2`; this is a first weak focus and allows at most one small cycle.
When `k=0,p!=0`, the next coefficient is `-4p`, with the same conclusion.
Thus no extra center multiplicity evades the condition `p*k!=0`.

## 5. Uniform endpoint division: the step that prevents false addition

Fix `L=3/2` and shrink the parameter neighborhood so `lambda<L`, where
`lambda=-(a+2)/a`. Proposition 3.2, including its resonant case, gives

\[
 D_u(s;\mu)=\delta_u(\mu)+\Delta_u(\mu)s^{\lambda(\mu)}
                         +R_u(s;\mu),\qquad R_u\in\mathcal F_L^\infty.
                                                               \tag{7}
\]

Here `delta_u=epsilon0*rho_u`, with `rho_u(mu*)=pi` by (2).
Both `Delta_u` and `R_u` vanish on the upper center set: on that set
`D_u` and `delta_u` vanish, and `L>lambda` forces the leading coefficient
to vanish separately.

Apply precisely the division in (4) to `Delta_u` and `R_u`. The
parameter derivatives required by the two divisions preserve uniform
flatness. It follows that, with constants independent of the sequence,

\[
 |\Delta_u|\le K N,\qquad
 |R_u|\le K N s^L,\qquad |\partial_sR_u|\le K N s^{L-1}.      \tag{8}
\]

In particular, whenever `epsilon0/N -> p != 0`,

\[
 \frac{D_u(s_n;\mu_n)}{N_n}\longrightarrow\pi p\ne0
          \quad\hbox{for every }s_n\to0.                    \tag{9}
\]

This is a statement about a joint limit, with no restriction to analytic
arcs or polynomial relations between `s_n` and the parameters. It
excludes **every upper endpoint cycle** in such a sequence. It uses the
uniform remainder from (7), not the pointwise compact Taylor expansion.

## 6. Two compact cycles leave at most one endpoint cycle

### Case A: both compact cycles are upper

Section 4 forces `p*k!=0`. Equation (9) excludes upper endpoint cycles.
Since `g=C*epsilon1` and `C->0`,

\[
 \frac{N}{|\epsilon _1|}\to0,\qquad
 \epsilon _0=o(\epsilon _1),\quad c_u=o(\epsilon _1),\quad
 c_l=2\epsilon _1-c_u\sim2\epsilon _1.                       \tag{10}
\]

For the reflected lower difference map the leading coefficient is,
using (3) and smoothness,

\[
 \Delta_l/\epsilon _1\longrightarrow-\pi.
\]

The lower flat remainder can be divided by the reversible parameters
`(epsilon0,epsilon1,epsilon2)`, since it vanishes when all three are zero.
By (10) its `s` derivative is `O(|epsilon1| s^{L-1})`. Thus

\[
 \partial_sD_l=s^{\lambda-1}
       [\lambda\Delta_l+O(|\epsilon _1|s^{L-\lambda})]
\]

has a fixed nonzero sign for small `s` and all sufficiently large `n`.
There is at most one lower endpoint zero, by Rolle's theorem. Reflection
gives the corresponding statement when both compact cycles are lower.

### Case B: one compact cycle in each annulus

Normalize instead by
`E=max(|epsilon0|,|epsilon1|,|epsilon2|)` and take a subsequential limit
`(p0,p1,p2)` of the normalized vector. If `p0=0`, equations (1) and its
reflection show that a compact root in each annulus requires

\[
 p_1+p_2/2=0,\qquad p_1-p_2/2=0.
\]

This contradicts the nonzero normalization. These assertions also hold
for a root shrinking to a center, using its leading trace coefficient.
Therefore `p0!=0`. Divide each flat remainder by the three reversible
parameters. The two difference maps divided by `E` tend to `pi*p0`
and `-pi*p0` in every endpoint joint limit. Neither endpoint is available.

### Compatibility theorem

For a sequence `mu_n -> mu*` carrying two compact limit cycles, the
number of simultaneous cycles tending to the two hemicycles is at most
one. If the compact distribution is `(1,1)`, it is zero.

This excludes the advertised **two compact plus three endpoint**
construction at this base, including nonanalytic parameter sequences.
It does not provide a neighborhood-wide total-cycle bound, because the
hypothesis of two compact limit cycles is essential to the proof.

## 7. What the resonant logarithmic invariant actually says

For `epsilon0=0`, in the infinity chart,

\[
 K(u,v)=1+\frac{2}{a+\epsilon _1u+\epsilon _2v+
            (b-2)u^2/4+(1-b)uv+bv^2}.
\]

At `a=-1`, an exact cancellation gives

\[
 (K_uK_v+K_{uv})(0,0)=2(b-1),\qquad
 G_2=\frac{2\pi\epsilon _1}{\sqrt{2-b-\epsilon _1^2}},
\]

and hence

\[
 F_3=-\frac{4\pi(b-1)\epsilon _1}
                   {\sqrt{2-b-\epsilon _1^2}}.              \tag{11}
\]

This is the source's logarithmic **invariant**. Its actual Dulac
coefficient is a positive unit times `F3` modulo the leading coefficient;
it is not legitimate to substitute (11) for that coefficient outright.
At `b=1` the independent logarithmic invariant vanishes. On `a=-1`,
its leading normal monomial is precisely `g=(b-1)*epsilon1`, the third
compact generator above. Resonance does not furnish a free extra knob.

The uniform first-order expansion (7) suffices for the compatibility
theorem. We have not derived the full broken-connection second-order
compensator expansion or its sharp simultaneous endpoint root bound.

## 8. Reproducible positive control

The limiting compact choice `(p,q,k)=(1,3,20)` gives two simple positive
roots of (6):

\[
 u_1\simeq0.208483927710504,\qquad
 u_2\simeq2.007524695383949.
\]

A rational polynomial arc realizing these leading generators is

\[
 a=-1,\quad b=1-20t,\quad \epsilon _0=t^2,\quad
 \epsilon _1=-t,\quad\epsilon _2=2t+6t^2.
\]

Here `c_u=3t^2`, `g=20t^2`. The lower endpoint's leading balance is
`s/t -> 1`. Original-field shooting at three positive values of `t`
and two tolerances produced the expected three sign brackets. This used
36 return-difference evaluations (72 half-flow integrations). At `t=1/10000`:

| Bracket coordinate | Left `D/t^2` | Right `D/t^2` |
|---|---:|---:|
| Upper `h in [1.2,1.8]` | 0.234513198 | -0.294864783 |
| Upper `h in [4,6]` | -0.247007943 | 0.221140008 |
| Lower `s/t in [0.5,1.5]` | -1.56119847 | 1.59286655 |

These are **NUM** controls. They check orientation, the mixed-generator
sign, and the allowed three-cycle scenario. No finite rational member
has been interval-certified here, and no claim of exactly three cycles
for these sampled fields is made.

Replay from repository root:

```bash
python resonant/check_resonant.py
python resonant/shoot_control.py
```

The first script verifies symbolic identities and independently integrates
the circle moments and mixed derivative at 70 decimal digits. The second
uses original quadratic equations and complete forward/backward half
returns. Logs and versions are in `resonant/data/`.

## 9. Remaining mathematical target

The local two-compact mechanism is excluded by the argument above.
The resonant line as a whole remains unresolved. In particular, this
strike does not exclude five cycles all tending to infinity, or a
one-compact/four-endpoint construction, or a different base `(-1,b0)`.

To reopen this exact base, specify a mechanism outside the two-compact
hypothesis and calculate the broken-connection joint compensator
coefficients with their parameter-uniform divided remainder. A claim
that the individual hemicycle bounds add is still insufficient.

No five-cycle field was found. No proof of `H(2)=4` was obtained.
