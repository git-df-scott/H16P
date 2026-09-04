# Independent red-team of the proposed finite-distance Hopf completion

Audit date: 2026-09-04. This audit did not rerun any return-map probes. It inspected `work/hopf_probe.py` and its existing JSON, independently checked the local algebra and compactification, and checked the relevant primary literature. The proposal remains **UNKNOWN**, with a valid conditional completion step and an unproved global coexistence premise.

## 1. Precisely what would suffice

Consider the exact rational quadratic family

\[
\dot x=y+x^2+xy,\qquad
\dot y=-10x^2+\frac{11}{5}xy+cy^2+\alpha x+\beta y.
\]

At a parameter pair \((c_*,\alpha_*)\), require \(\beta=0\), \(\alpha_*<0\),

\[
K=-\alpha_*\left(\frac{11}{5}c_*-1\right)-42>0,
\]

three distinct hyperbolic periodic orbits surrounding the origin in stable/unstable/stable order from inside out, and one distinct hyperbolic unstable periodic orbit surrounding the other finite focus. Every cycle must have positive distance from the origin, other cycles and infinity.

If this premise is proved, sufficiently small negative \(\beta\) produces one additional unstable Hopf cycle around the now stable origin, and all four old hyperbolic cycles persist. Thus the resulting configuration is four nested cycles plus one remote cycle. This conditional implication is sound; its premise is already a major global coexistence problem. No identified theorem supplies that premise or proves that it occurs in the proposed box.

The Lyapunov coefficient has the stated sign. Put \(\omega=\sqrt{-\alpha}\), \(u=x\), \(v=-y/\omega\) and \(\tau=\omega t\). At \(\beta=0\), the system becomes

\[
u'=-v+u^2/\omega-uv,\qquad
v'=u+10u^2/\omega^2+(11/(5\omega))uv-cv^2.
\]

The standard planar cubic radial coefficient in this normalization is

\[
l_1=\frac{\omega^2((11/5)c-1)-42}{8\omega^3}.
\]

Changing the positive time normalization changes the numerical coefficient, not its sign. The linear eigenvalues have real part \(\beta/2\) near zero, so the crossing is transverse. The normal form \(r'=\mu r+l_1r^3+\cdots\), with \(\mu<0\) and \(l_1>0\), gives the claimed unstable cycle.

## 2. What the known four-cycle construction does and does not provide

Kuznetsov–Kuznetsova–Leonov's published family is exactly this reduced form with variable \(a_2,b_2,c_2,\alpha_2,\beta_2\). Its displayed example is \((-10,2.2,0.7,-72.7778,0.0015)\). Its construction operates on the side

\[
\alpha_2>\frac{a_2(2+b_2)}{b_2c_2-1},\qquad b_2c_2>1,
\]

so here \(K<0\), followed by a positive trace perturbation. At \(c=0.7\), the boundary is \(\alpha=-700/9\). The proposed \(K>0,\beta<0\) target is on the opposite side. Therefore a known four-cycle example cannot simply be assumed to retain its four finite cycles while moving onto the target Hopf surface. A cycle may collapse into the weak focus, collide in a fold, or reach an infinity graphic during continuation. These are exactly the events that require detection.

Source: N. V. Kuznetsov, O. A. Kuznetsova and G. A. Leonov, *Visualization of Four Normal Size Limit Cycles in Two-Dimensional Polynomial Quadratic System*, Differential Equations and Dynamical Systems 21 (2013), 29–34, DOI [10.1007/s12591-012-0118-6](https://doi.org/10.1007/s12591-012-0118-6), especially the reduced system (2.2), coefficient-domain discussion and Fig. 4. The [author-uploaded full text](https://www.researchgate.net/publication/257789723_Visualization_of_Four_Normal_Size_Limit_Cycles_in_Two-Dimensional_Polynomial_Quadratic_System) is indexed. Direct retrieval of the author's spbu PDF URL returned HTML, so `work/kkl_redteam.pdf` is not a usable PDF. Do not cite it as a downloaded primary PDF.

## 3. Accepted local/distribution results do not close this target

Zegeling's 2024 Proposition 6.5 explicitly allows distributions \((n,1)\) or \((n,0)\), with no upper bound on \(n\), for two foci, a complex pair of finite singularities, at least one elementary saddle at infinity, and at least one weak focus. Its global Theorem 1.2 proves the distribution restriction, not \(n\le3\). A first-order weak focus with three surrounding cycles and one remote cycle would satisfy this distribution restriction.

The prohibition of cycles surrounding an order-three weak focus, and the at-most-one result around an order-two weak focus, do not apply when \(K>0\), because the focus has order one. They would apply on appropriately verified higher-order degeneracy strata, and warn against an infinitesimal passage through such a stratum that assumes remote cycles survive.

Source: André Zegeling, *Nests of limit cycles in quadratic systems*, Advances in Nonlinear Analysis 13 (2024), DOI [10.1515/anona-2024-0012](https://www.degruyterbrill.com/document/doi/10.1515/anona-2024-0012/html), Proposition 6.5, Corollary 6.8 and Theorem 1.2. Zhang–Zhao's [2001 weak/strong focus paper](https://doi.org/10.1007/s11766-001-0018-y) likewise states strong restrictions for weak orders two and three; its abstract does not assert a global two-cycle bound around an order-one weak focus.

No applicable accepted finite-graphic, infinity-graphic or local theorem was identified that forbids the stated order-one base. This is an audit finding, not a theorem proving that no obstruction exists. Theorem N and perturbative center-cyclicity bounds are inapplicable without a demonstrated reduction to their exact integrable perturbation setting.

## 4. Exact finite geometry and useful algebraic gates

The vertical line \(x=-1\) satisfies \(\dot x=1\) everywhere. It is a one-way transverse barrier, so a periodic orbit cannot cross it. The origin nest lies in \(x>-1\), and a remote nest must lie in \(x<-1\).

All other equilibria satisfy

\[
y=-\frac{x^2}{1+x},\qquad
(61/5-c)x^3+(111/5+\beta-\alpha)x^2
 +(10+\beta-2\alpha)x-\alpha=0.
\]

Write \(m=-\alpha>0\). At \(\beta=0\), the cubic discriminant is

\[
-\frac{1}{125}\left(3375c^2m^2+500cm^3-19950cm^2+199500cm
-500000c-605m^2+17424m-60500\right).
\]

The intended two-focus regime needs the appropriate finite root count, nondegeneracy, and negative Jacobian focus discriminants; the rectangular box alone is not such a certificate. The Jacobian to use is

\[
J=\begin{pmatrix}2x+y&1+x\\-20x+(11/5)y+\alpha&(11/5)x+2cy+\beta\end{pmatrix}.
\]

In the proposed box \(c\in[1/2,3/2],m\in[10,200]\), the coefficient gate \(K>0\) is possible only for \(c>11/20\); it requires

\[
m>\frac{42}{(11/5)c-1}.
\]

At a remote equilibrium \(x_*<-1\), its trace at \(\beta=0\) is

\[
T_* = \frac{x_*}{1+x_*}\left[\frac{21}{5}+\left(\frac{16}{5}-2c\right)x_*\right].
\]

For the desired unique remote **unstable** cycle, the remote focus must be stable. Otherwise the outward flow near an unstable focus and inward flow just inside an unstable cycle force an additional separating periodic orbit, contradicting the intended one-cycle remote nest (and, with three origin cycles, the distribution theorem). Thus impose \(T_*<0\), in addition to the focus condition.

When the cubic has one real remote root, this trace condition is equivalent to

\[
m>m_T(c)=\frac{21(1000c^2+1021c+481)}{50(2c+1)^2(8-5c)}.
\]

For example, \(m_T(0.7)=27.31458\ldots\), \(m_T(1)=38.92\), \(m_T(1.3)=75.58055\ldots\), and \(m_T(1.5)=223.78125\). Consequently the entire \(c=1.5\) face with \(m\le200\), in the intended one-remote-root regime, fails the remote stability requirement. This is a genuine prune of the box, not a global exclusion of the mechanism.

## 5. Infinity is not constant throughout the proposed box

The homogeneous angular polynomial is

\[
xQ_2-yP_2=x[-10x^2+(6/5)xy+(c-1)y^2].
\]

There is always the vertical projective direction. In its chart the radial/angular eigenvalues are \((-c,1-c)\), so it is a saddle for \(c<1\), a node for \(c>1\), and nonhyperbolic at \(c=1\). The other directions \(v=y/x\) solve

\[
(c-1)v^2+(6/5)v-10=0,
\]

whose discriminant is \(36/25+40(c-1)\). Hence a pair of real directions appears at

\[
c=241/250=0.964.
\]

For a simple finite slope, chart eigenvalues are \(-(1+v)\) and \(6/5+2(c-1)v\). Thus the box contains several distinct infinity configurations and two degeneration surfaces. The KKL sufficient infinity condition \(4a_2(c_2-1)>(b_2-1)^2\), with the fixed coefficients here, is exactly \(c<241/250\). Outside it the known sufficient mechanism cannot be imported, although the algebra does not forbid finite cycles.

Split continuation/validation at \(c=241/250\) and \(c=1\), or explicitly build the needed compactified transitions. Fixed Euclidean cutoffs cannot determine whether a missing return is an escape, a very large cycle, a long saddle passage or numerical failure.

## 6. Limits of the current probe and a reviewable success criterion

The twelve sampled parameter pairs and ten section coordinates per pair are useful exploratory data only. The \(K>0\) samples did not establish the required three origin roots and the probe does not simultaneously certify the remote orbit. Fixed integration horizon 50, sparse sign brackets, swallowed exceptions and unresolved nonreturns prevent an absence conclusion. The existing JSON correctly labels itself nonrigorous and bounded.

The returned quantity \(\exp\int\mathrm{div}F\,dt\) is the Floquet multiplier **at an actual periodic orbit**. Away from a fixed point, the derivative of the section return map \(P\) on \(y=0\) is

\[
P'(x_0)=\exp\!\left(\int\mathrm{div}F\,dt\right)
\frac{Q(x_0,0)}{Q(P(x_0),0)}.
\]

This section-speed correction is needed for interval derivative bounds on brackets. It becomes 1 at a fixed point, so multiplier values evaluated at sufficiently accurate roots are consistent with the intended stability test; their numerical accuracy still needs validation.

A concrete success certificate would consist of: (i) an exact rational parameter pair with strict \(K>0\) and verified finite equilibrium types; (ii) four disjoint section intervals with validated return-map sign changes and multiplier intervals separated from 1, with origin stability order S/U/S and remote U; (iii) orbit enclosures proving the correct focus is enclosed and finite separation from infinity; and (iv) an explicit negative \(\beta\) range preserving those four enclosures and containing the local Hopf cycle, or direct validation of all five at one rational \(\beta<0\).

The restrictions \(a_2=-10,b_2=11/5\) and the finite \((c,\alpha)\) box are experimentally convenient, not known to be complete for this mechanism. Exhausting this box would close only that bounded slice. It would not exclude a first-order Hopf completion in another quadratic slice, and it would not prove \(H(2)\le4\).

## Verdict

**Valid conditional mechanism; unproved and possibly obstructed global premise; no applicable accepted theorem found that closes the complete order-one target.** Retain it as a bounded exploratory strike. Use the algebraic gates above before continuation, preserve the distinction from the established KKL wedge, and define success by simultaneous validated coexistence rather than by adding cycle counts from different parameter values.
