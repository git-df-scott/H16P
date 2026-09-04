# Search-space reality check

## Raw and essential dimensions

Each polynomial of degree at most two has six coefficients, so `(P,Q)` has 12.
For generic fields, phase portraits are unchanged under:

- an affine change of plane coordinates: 6 dimensions;
- a common nonzero time scaling: 1 dimension.

Thus the generic orbit space has **five essential real parameters**. This
dimension count is local. Degenerate fields have stabilizers, affine charts
break when their chosen coefficients vanish, and orientation/time-reversal
choices must be tracked. There is no single safe five-dimensional rectangular
normal form covering all quadratic phase portraits.

Projectivizing the 12 coefficients removes overall scale and gives a compact
coefficient-direction space before affine quotient. It does not make a finite
resolution exhaustive: limit cycles can approach multiple cycles, centers,
separatrix graphics, infinity, or vanishing scales arbitrarily closely.

## Useful charts

### Two-singularity / Shi chart

When the relevant normalizations are nondegenerate, use

\[
\begin{aligned}
\dot x&=\lambda x-y+\ell x^2+mxy+y^2,\\
\dot y&=x+ax^2+bxy.
\end{aligned}
\]

This is a five-parameter chart `(lambda,l,m,a,b)`. The third-order weak-focus
stratum through Shi's seed has

```text
lambda = 0
m      = 5 a
b      = 3 l + 5
```

with additional nonvanishing inequalities. Shi's point is `(l,a)=(-10,1)`.

### Kuznetsov chart

\[
\dot x=y+x^2+xy,\qquad
\dot y=ax^2+bxy+cy^2+\alpha x+\beta y.
\]

This is already five dimensional and contains a visible four-cycle region.
It is excellent for continuation, though it does not represent every affine
class.

### Center charts

The standard complex classification includes Hamiltonian `Q3^H`, reversible
`Q3^R`, generalized Lotka--Volterra `Q3^LV`, codimension-four `Q4`, and the
Hamiltonian triangle. These charts are adapted to Melnikov/Abelian-integral
calculations but cover center strata, not generic global moduli.

For `Q4`, Zhao reduces the center by affine changes and time rescaling to the
Hamiltonian form

\[
\dot x=-1-(\kappa-1)x^2+\kappa y^2,\qquad
\dot y=-2(\kappa-1)x(x-y),\qquad \kappa>1,
\]

with first integral

\[
H=\frac23(\kappa-1)x^3-(\kappa-1)x^2y+\frac\kappa3y^3-y.
\]

**Correction from the Q4 audit:** this cubic Hamiltonian chart is reached
through a double ramified cover and inversion, not solely affine changes and
time rescaling. An arbitrary quadratic perturbation in this chart need not be
a quadratic perturbation of the original Q4 field. The four-function
generating space is valid, but candidate coefficients must pass the
original-coordinate realization gate in
[ZERO_TO_CYCLE.md](ZERO_TO_CYCLE.md).

The period annulus is
`-2/3 < h < -2/(3 sqrt(kappa))`. Under a general small quadratic
perturbation, Zhao reduces the first nonzero generating function to a
four-term linear combination of complete elliptic integrals. After
projectivizing the four coefficients, the numerical root problem has three
direction parameters plus the center modulus: a focused four-dimensional
attack surface.

## Why a grid is not sane

Even granting one five-dimensional chart:

| Resolution per parameter | Boxes |
|---:|---:|
| 100 | `10^10` |
| 1,000 | `10^15` |
| 1,000,000 | `10^30` |

One validated return calculation may take `10^3`--`10^6` integration steps,
and near-degenerate boxes need subdivision and multiprecision. More important,
no fixed resolution is complete. The Songling construction uses parameter
scales from `10^-13` through `10^-200` (or `10^-250` in a reported original
choice), so even a trillion evenly spaced points along one coordinate can miss
the relevant region by almost two hundred orders of magnitude.

## Structured reductions that are sane

1. Solve algebraic focus conditions and compute Lyapunov quantities exactly.
2. Continue zero sets of a separatrix-splitting or Melnikov function rather
   than sample all coefficients.
3. Use Sturm/resultant/interval-polynomial tests to reject boxes before ODE
   integration.
4. Classify equilibria and Poincare-sphere singularities with interval
   discriminants.
5. Evaluate scalar return maps only on topology-compatible boxes.
6. Promote a candidate to arbitrary precision when displacement roots persist
   under precision and integrator changes.
7. Run interval verification only after five separated candidate roots exist.

## Can interval methods prune parameter boxes?

Yes, locally. Parameter-dependent interval integration can reject a box when a
return displacement has fixed sign or a trapping boundary cannot exist.
Bendixson--Dulac criteria and equilibrium topology are cheaper. But wrapping
grows with both flight time and parameter width; boxes near homoclinic loops or
multiple cycles subdivide without a known uniform stopping depth. Therefore
interval pruning makes a narrow family rigorous, not the global ocean finite.

## Computational assessment

- Candidate generation near explicit bifurcations: bounded and parallelizable.
- Verification of one hyperbolic candidate: bounded and usually efficient.
- Exhaustive coverage of a chosen compact parameter box away from degeneracy:
  possible in principle if quantitative margins are imposed.
- Exhaustive coverage of all quadratic systems: no known finite attack.
