# Reversible re-seed: exact geometry, a missing unfolding direction, and simultaneous controls

Date: 2026-09-05. Base commit: `79001f7`.

**The full reversible route is NOT closed, and no five-cycle field was
found.** This strike closes the finite-saddle-loop version for every
two-center seed in the stated reversible normal form. It supplies the
replacement geometry and two locally complete perturbation charts, finds
an omitted direction at the nonhyperbolic-infinity boundary, and executes
a bounded first-order construction search with a reproduced four-cycle
positive control. It does not prove `H(2)=4`.

## 1. What is actually closed

Start with the proposed reversible normal form

\[
 \dot X=-Y(1+kX),\qquad \dot Y=X+pX^2+qY^2.
\]

For `k != 0`, apply the exact affine transformation

\[
 x=kY/2,\qquad y=(1+kX)/2,\qquad a=2q/k,\quad b=2p/k.
\]

No time rescaling is needed. The resulting field is

\[
 P_0=(b-2)/4+(1-b)y+ax^2+by^2,\qquad Q_0=-2xy.       \tag{1}
\]

Its axis equilibria and their Jacobians are

\[
 C_u=(0,1/2),\quad J_u=\begin{pmatrix}0&1\\-1&0\end{pmatrix},
 \qquad
 C_l=(0,(b-2)/(2b)),\quad
 J_l=\begin{pmatrix}0&-1\\(2-b)/b&0\end{pmatrix}.
\]

They are nondegenerate centers precisely when `0<b<2`. Reversibility
and the local first integral establish centers, rather than merely
trace-zero linearizations. The other possible equilibria lie on `y=0`
and satisfy `x^2=(2-b)/(4a)`. Therefore:

\[
 \boxed{\text{exactly two real finite equilibria, both centers}
        \iff 0<b<2,\ a\le 0.}                         \tag{2}
\]

The formula for the off-axis pair is used only when `a != 0`; at `a=0`
the constant `(b-2)/4` is nonzero on `y=0`, so there is no additional
finite equilibrium. In the original parameters, (2) is
`0<p/k<1`, `q/k<=0`. If `k=0`, the possible second axis equilibrium
has determinant `-1`, and there is no second center. Off-axis equilibria
in the original chart have triangular Jacobians and cannot supply a center.

**Consequence:** a re-seed with two centers has no finite saddle at all.
The old construction, a finite saddle-loop boundary with a separate
center/focus supplying the fifth cycle, cannot be repaired within this
two-center geometry. This is an exact structural exclusion, not a
failure to locate a loop numerically. It does not exclude boundaries
through infinity or finite-amplitude folds after perturbation.

## 2. The re-seeded class has three boundary regimes

For `a` outside `{0,-1,-2}` and `y>0`, set

\[
 R(y)=-\frac{b}{a+2}y^2-\frac{1-b}{a+1}y
                       +\frac{2-b}{4a},\qquad
 H=y^a[x^2-R(y)].                                    \tag{3}
\]

The integrating factor is `mu=y^(a-1)` and
`H_x=-mu Q_0`, `H_y=mu P_0`. On the negative half-plane use
`r=-y>0` and replace `(1-b)` by `-(1-b)` in the potential.
All fractional powers are of positive quantities.

At the exceptional values, the upper-half-plane first integrals can be
chosen as

\[
\begin{array}{ll}
a=0:&H=x^2+\frac b2 y^2+(1-b)y+\frac{b-2}{4}\log y,\\
a=-1:&H=\frac{x^2}{y}+by+(1-b)\log y+\frac{2-b}{4y},\\
a=-2:&H=\frac{x^2}{y^2}+b\log y-\frac{1-b}{y}
                                      +\frac{2-b}{8y^2}.
\end{array}                                                        \tag{4}
\]

The geometry is:

| Shape | Outer boundary of each center annulus | Scope of an existing hemicycle theorem |
|---|---|---|
| `-2<a<0` | Invariant line `y=0` plus a semicircle at infinity | Hyperbolic hemicycle; resonance at `a=-1` |
| `a<-2` | A branch of `x^2=R(y)` plus an arc at infinity | A different two-saddle boundary, called a bicycle |
| `a=-2` or `a=0` | Nonhyperbolic hemicycle boundary | Hyperbolic theorem does not apply |

For example, in the bicycle sector the potential
`V(y)=-R(y)y^a` decreases up to `y=1/2`, then increases to `0` as
`y->infinity`. The center energy is negative. Its ovals have
`H_center<h<0`, and the outer level is exactly `x^2=R(y)`.
The quadratic `R` has positive leading coefficient and negative constant
term, so it has one positive and one negative root. This gives the upper
and lower branches and confirms the two distinct boundary components.
In the hemicycle sector, the potential tends to positive infinity at both
ends of the positive axis, including the logarithmic cases in (4).

In the infinity chart `u=1/x`, `v=y/x`, after the usual time change,

\[
\begin{aligned}
 u'&=-u[a+bv^2+(1-b)uv+(b-2)u^2/4],\\
 v'&=-v[a+2+bv^2]-(1-b)uv^2-(b-2)u^2v/4.
\end{aligned}                                                       \tag{5}
\]

The horizontal direction has eigenvalues `-a` and `-(a+2)`.
When `a<-2`, the other directions satisfy `v^2=-(a+2)/b` and have
eigenvalues `2` and `2(a+2)`: they are hyperbolic saddles.
At `a=-2`, the angular equation starts with `-bv^3`; at `a=0`, the
radial eigenvalue vanishes. Neither boundary is covered by taking a
hyperbolic result at a nearby parameter and passing to the limit.

A concrete bicycle seed is

\[
 (a,b)=(-7/3,1),\quad
 \dot x=-1/4-(7/3)x^2+y^2,\quad \dot y=-2xy,\quad
 R(y)=3y^2-3/28.                                      \tag{6}
\]

It has exactly the two centers `(0,+/-1/2)` and no finite saddles.
It is a valid explicit re-seed, not a five-cycle construction.

[Marín–Villadelprat](https://arxiv.org/html/2501.16924v1), Theorems B/C,
give individual hemicycle cyclicity two and simultaneous cyclicity at
most three off resonance in `(-2,0)x(0,2)`. Those results count cycles
approaching the specified boundaries; they do not bound additional
compact cycles. Their discussion of the bicycle sector distinguishes
it from these hemicycles.

## 3. A missing perturbation direction at a=-2

The usual five-parameter family is

\[
 P=P_0+\epsilon_1x+\epsilon_2xy,\qquad
 Q=Q_0+\epsilon_0,\qquad (a,b,\epsilon_0,\epsilon_1,\epsilon_2). \tag{7}
\]

To check that it represents arbitrary nearby quadratic perturbations,
form the twelve coefficient columns supplied by six infinitesimal affine
coordinate changes, time rescaling, and the five displayed parameters.
The exact determinant, in the ordering specified by the checker, is

\[
                         -16(a+2).                  \tag{8}
\]

The inverse function theorem therefore makes (7) a locally complete
chart modulo affine transformations and time rescaling whenever `a!=-2`.
At `a=-2` this argument fails for a real algebraic reason.

Replace the `epsilon2 xy` term in `P` with a `gamma x^2` term in `Q`:

\[
 P=P_0+\epsilon_1x,\qquad Q=Q_0+\epsilon_0+\gamma x^2,
 \qquad (a,b,\epsilon_0,\epsilon_1,\gamma).             \tag{9}
\]

The same determinant at `a=-2` is now

\[
                              -48b\ne0.              \tag{10}
\]

Thus (9) is a complete five-parameter local unfolding there for every
`0<b<2`. This **replaces an ineffective direction**; it does not create
a sixth independent system parameter or a fourth independent first-order
Melnikov function.

The loss is visible directly in the moments. For an upper oval domain
`D_h`, write

\[
 I_j=\iint_{D_h}y^j\,dx\,dy,\qquad
 J=\iint_{D_h}x^2y^{a-2}\,dx\,dy.
\]

The zero flux of `y^(a-2) x (P_0,Q_0)` through an oval gives

\[
 (a+2)J+bI_a+(1-b)I_{a-1}+\frac{b-2}{4}I_{a-2}=0.    \tag{11}
\]

At `a=-2`, the three old moments are dependent. Indeed, the old
perturbation direction

\[
 (\epsilon_0,\epsilon_1,\epsilon_2)
       =((2-b)/(12b),(1-b)/b,1)
\]

has identically zero first Melnikov function in **both** annuli. The
`gamma` direction supplies the omitted `J` moment. A first-order search
using only (7) on `a=-2` would miss this direction entirely.

At the explicit boundary seed `a=-2,b=1`, put `u=x/y`, `v=1/y` in the
upper half-plane. Then

\[
 \dot u=1/v-v/4,\qquad \dot v=2u,\qquad
 H=u^2+v^2/8-\log v.                                 \tag{12}
\]

This makes the missing-direction calculation particularly tractable.
The three moments for (9) are area, the area moment of `v`, and the
area moment of `u^2/v` in the logarithmic potential well.

## 4. Same-parameter Melnikov calculation

For (7), weighted divergence has the even-in-x part

\[
 \epsilon_1 |y|^{a-1}+\epsilon_2 y|y|^{a-1}
                   +(a-1)\epsilon_0 |y|^{a-1}/y.
\]

Let `v_u,w_u` be the upper weighted averages of `y,1/y`, and let
`v_l,w_l` be the lower weighted averages of `|y|,1/|y|` with weight
`|y|^(a-1)`. Up to nonzero orientation factors, the two normalized
Melnikov functions are

\[
 F_u=\epsilon_1+\epsilon_2v_u+(a-1)\epsilon_0w_u,
 \qquad
 F_l=\epsilon_1-\epsilon_2v_l-(a-1)\epsilon_0w_l.       \tag{13}
\]

Set `epsilon1=-c`, `epsilon2=-m`, `epsilon0=1/(a-1)`.
Searching for zeros in both annuli then amounts to intersecting one
line with the upper curve `(v_u,w_u)` and the reflected lower curve
`(-v_l,-w_l)`, using the same intercept `c` and slope `m`.

The code visits the projection-order cells determined by every pair of
sampled vertices, then sweeps the common intercept. It looks for `(4,1)`,
`(1,4)`, `(5,0)`, or `(0,5)` sign-bracket patterns. A synthetic oscillating
curve triggers the five-root detector in either annulus; a parabola
provides a negative software control. These artificial curves are not
claimed to be quadratic-system moment curves.

This is a finite **NUM** search. Three functions alone do not give a
two-zero theorem. The finite line search does not exclude unsampled
oscillations, coalescing roots, cycles shrinking to centers, higher-order
degeneracies, or cycles tending to the boundaries.

## 5. Bounded run and independent four-cycle control

The run contains 54 shapes away from the exceptional logarithmic values,
including the published Yu–Zeng shape and points on either side of the
upper/lower cubic-focus lines. Each annulus has 41 energy samples, with
160-point Gauss quadrature. A separate run uses 10 shapes on `a=-2` and
`a=0`, with the repaired chart where necessary, 41 finite energy samples,
and 200-point quadrature. The complete grids and all curves are saved.

**No five-cycle first-order candidate was detected in these 64 finite
shape samples.** This is not an exclusion for any open parameter region.
The regular run took about four seconds in this environment; it used
moment quadrature, not an ODE sweep. Raw polyline counts are diagnostics.

The published Yu–Zeng shape in the repository transforms exactly to

\[
 (a,b)=(-671/450,7/15),\qquad \lambda=229/671.
\]

It lies strictly in the nonresonant hemicycle sector. The harness finds
a same-parameter `(3,1)` Melnikov sign pattern there. This is consistent
with the near-integrable construction in
[Yu–Han](https://arxiv.org/abs/1002.1055) and with the exact numerical
shape recorded in `FOUR_CYCLE_FRONTIER.md`.

For an independently checked control, the search also supplies the
following exact rational arc:

\[
\begin{aligned}
 \dot x&=-5/12+(2/3)y-(7/4)x^2+(1/3)y^2
              -\tau[(31379/25000)x+(7517/5000)xy],\\
 \dot y&=-2xy-4\tau/11.                              \tag{14}
\end{aligned}
\]

At `tau=0`, both centers and both annuli are genuine. Independent
65-decimal-digit quadrature, using a different integration coordinate
from the search, gives the normalized first-order signs:

| Annulus | Base energy h | Normalized M |
|---|---:|---:|
| Upper | 1 | -0.000987368720149 |
| Upper | 1.5 | +0.000051254063109 |
| Upper | 2 | -0.000060717220792 |
| Upper | 3 | +0.011394645644968 |
| Lower | 10 | +0.513430977359841 |
| Lower | 16 | -0.437418481820869 |

Original quadratic-field integration at `tau=1/10000` and `1/20000`,
each at two tolerances, preserves all four brackets. For example, at
`tau=1/10000`, `rtol=2e-13`, the forward-minus-backward half-return
differences divided by `tau` are

| Annulus | Energies in order | D/tau in order |
|---|---|---|
| Upper | 1, 1.5, 2, 3 | +0.0004045873, -0.0000205639, +0.0000230962, -0.0038631761 |
| Lower | 10, 16 | +0.1633986134, -0.0995511139 |

These are **ordinary numerical controls**, not interval certificates.
No exact root count, Floquet bound, or five-cycle result is claimed.
The first-order quadrature is high precision but is not an interval
enclosure. The direct control used 24 return-difference evaluations,
comprising 48 half-flow integrations, in addition to short startup steps.
No prior KKL or Shi/Chen–Wang sweep was repeated.

## 6. Precise route decision

| Claim or proposed mechanism | Result of this strike |
|---|---|
| Two-center re-seed retaining the old finite saddle-loop boundary | Excluded exactly |
| All two-center re-seeds are the previously studied resonant base | False; (2) and (5) display the other regimes |
| Usual five-parameter chart remains complete at `a=-2` | False; repaired by (9), verified by (10) |
| The missing direction supplies a fourth independent first-order compact coefficient | False; it replaces the lost direction in (11) |
| Five first-order compact cycles in the finite sampled profiles | None detected; NUM only |
| Full reversible route closed | **No** |
| Explicit five-cycle counterexample | **No** |

The remaining mathematical task is a simultaneous zero/cyclicity bound
or construction for the **closed** annuli, with parameter-uniform
endpoint control. For the concrete nonhyperbolic seed `(-2,1)`, this must
use the unfolding (9); the old chart is insufficient. The cubic angular
equation in (5) and the logarithmic Hamiltonian (12) specify the object,
but this strike does not calculate its full joint Dulac expansion or
identify its label in the 121-graphic catalogue. No entire infinity
route is declared closed or admitted by that observation alone.

## 7. Replay and evidence

See [`reversible_reseed/README.md`](reversible_reseed/README.md).
The exact identities, full parameter grids, moment curves, line-search
controls, high-precision values, original-field returns, and hashes are
included. The derivations are campaign mathematics with executable
algebraic checks; they have not received independent mathematical review.

This work deliberately leaves the wider route marked open. A numerical
failure, a three-dimensional function space, and a closure of one
geometric mechanism cannot be promoted to a proof of `H(2)=4`.
