# A first-order check at the symmetric resonant seed

Date: 2026-09-07, America/Edmonton. This is a bounded calculation from the overnight investigation, not a five-cycle counterexample or a closed-annulus cyclicity theorem. No novelty claim is made before completing the repository and literature review.

The repository's reversible moment search omits `a=-1`; its separate boundary search covers `a=-2` and `a=0`. At the omitted symmetric point `a=-1,b=1`, the first-order compact calculation can be done exactly.

Consider

\[
\dot x=-\tfrac14-x^2+y^2+\tau(e_1x+e_2xy),\qquad
\dot y=-2xy+\tau e_0.
\]

The base field has centers at `(0,±1/2)`. In the upper half-plane,

\[
H=\frac{x^2+y^2+1/4}{y},\qquad X_0H=0.
\]

For `h>1`, the oval and its disk are

\[
x^2+(y-h/2)^2=(h^2-1)/4.
\]

The disk lies strictly above `y=0`. Put `k=h/2`, `R=sqrt(h²-1)/2`, so `k²-R²=1/4`. Its moments satisfy

\[
J_j=\iint_D y^{-j}\,dx\,dy,\qquad
J_1=2\pi(k-\sqrt{k^2-R^2}).
\]

To verify the last formula directly, substitute `y=k+R cos θ` into the area integral. Then

\[
J_1=2\int_0^\pi\frac{R^2\sin^2\theta}{k+R\cos\theta}\,d\theta.
\]

The numerator divided by the denominator is `k-R cos θ-(k²-R²)/(k+R cos θ)`. The remaining integral is `π/sqrt(k²-R²)`, obtained by `t=tan(θ/2)`. Differentiating at fixed `R`, over the translated fixed disk, gives `J₂=-∂kJ₁` and `J₃=-(∂kJ₂)/2`. This differentiation is legitimate locally for `k>R`, where all integrands and their derivatives are bounded on the disk. Substitution yields

\[
J_1=\pi(h-1),\quad J_2=2\pi(h-1),\quad
J_3=2\pi(h^2-1).
\]

The positive integrating factor is `|y|^-2`. Up to a nonzero orientation factor, the first energy displacement is the integral of weighted perturbation divergence over the oval disk. Division by the positive `J₂` gives

\[
F_u(h)=e_1+\tfrac12 e_2-2e_0(h+1),\qquad
F_l(h)=e_1-\tfrac12 e_2+2e_0(h+1).
\]

Thus each **nonidentically-zero first-order function** has at most one compact interior root at this fixed base. If `e₀≠0`, any root is simple. For example, the common direction `(e₀,e₁,e₂)=(1,1,12)` gives roots `h_u=5/2`, `h_l=3/2`. This is a two-root first-order control, not an explicit finite-τ cycle certificate.

The calculation does **not** bound higher-order unfoldings when a first-order function vanishes identically, joint limits as the shape varies, or cycles escaping to the hemicycle boundaries. The energy-domain endpoint `h=1` also requires separate small-amplitude analysis. In particular, it does not close the resonant route.

Verification: independent symbolic checks of the first integral, circular levels, and all three moment formulas passed. A separate 60-decimal-digit quadrature at `h=1.001,1.5,2.5,10,100` agreed to relative error below `10^-45`. That quadrature is nonvalidated numerical corroboration; the algebra and integral evaluation above provide the exact argument. No ODE evaluations were used.

Repository context: [reversible re-seed report](https://github.com/git-df-scott/H16P/blob/45f4ea9b4ab448bdd36702036244faa4f15c9819/REVERSIBLE_RESEED_2026_09_05.md), [moment search](https://github.com/git-df-scott/H16P/blob/45f4ea9b4ab448bdd36702036244faa4f15c9819/reversible_reseed/moment_search.py), [boundary search](https://github.com/git-df-scott/H16P/blob/45f4ea9b4ab448bdd36702036244faa4f15c9819/reversible_reseed/boundary_search.py). These links establish the inherited setup and omitted sample, not external verification of this calculation.
