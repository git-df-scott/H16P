# Astra feasibility and attack matrix

## Scores

For `literature saturation` and `near-miss risk`, a high score is adverse; for
the other rows, high is favorable.

| Criterion | Score | Reason |
|---|---:|---|
| Finite searchability | 3/10 | Five essential generic parameters, multiple charts, no known finite-resolution global cover |
| Candidate-generation quality | 6/10 | Strong focus, Melnikov, separatrix, and continuation equations; still no established path to a fifth |
| Rigorous certificate quality | 9/10 | Five hyperbolic Poincare fixed points can be independently interval-certified; exact global count is unnecessary |
| Compute scaling | 4/10 | Narrow scalar-function searches scale well, but raw grids scale as `N^5` and validated integration wraps near degeneracy |
| Literature saturation (adverse) | 8/10 | The obvious local and near-integrable mechanisms have been studied for decades |
| Probability failure teaches mathematics | 7/10 | A certified zero bound or separatrix obstruction in one selected family is publishable campaign progress |
| Endless numerical near-miss risk (adverse) | 9/10 | Weak multipliers, unstable cycles, saddle flight times, and 200-digit parameter hierarchies create convincing artifacts |

## Verdict

# YELLOW — only attack specified narrow families

Do not optimize twelve raw coefficients and do not launch a generic grid. The
three attacks below are ordered by mathematical leverage.

## Attack 1 — `Q4` elliptic-integral zero hunt

### Exact family and normalization

Use Zhao's affine/time-normalized codimension-four Hamiltonian center

\[
\dot x=-1-(\kappa-1)x^2+\kappa y^2,\qquad
\dot y=-2(\kappa-1)x(x-y),\qquad \kappa>1,
\]

and add an explicitly constant-in-`eta` general quadratic perturbation

\[
\dot x=P_0+\eta\!\sum_{i+j\le2}p_{ij}x^iy^j,\qquad
\dot y=Q_0+\eta\!\sum_{i+j\le2}q_{ij}x^iy^j.
\]

Use exact rational `kappa`, `p_ij`, `q_ij`, and `eta`. Zhao's reduction of the
first nonzero Melnikov function is

\[
I(h)=\mu_1 hI_{0,0}+\mu_2I_{1,0}+\mu_3I_{0,1}
+\mu_4(2I_{-1,0}+3\kappa hI_{-1,1}).
\]

The integrals range over the ovals of

\[
H=\frac23(\kappa-1)x^3-(\kappa-1)x^2y+\frac\kappa3y^3-y=h,
\quad -\frac23<h<-\frac{2}{3\sqrt\kappa}.
\]

Normalize `||mu||_2=1` with the first nonzero component positive. Search the
projective three-space `[mu]` and eight exact `kappa-1` charts

```text
[2^-8,2^-4], [2^-4,1], [1,2], [2,4],
[4,8], [8,16], [16,32], [32,63]
```

This is four dimensional, rather than the full perturbation-coefficient
space. Recover exact `p_ij,q_ij` by solving the published linear Melnikov map
only after a root pattern exists.

### Candidate generation

- Interval/Chebyshev-evaluate the four elliptic-integral basis functions.
- For ordered rational samples `h0<...<h5`, solve linear inequalities forcing
  alternating signs of `I(hj)`. Five alternations force five intervening roots.
- Continue near vanishing generalized Vandermonde determinants and endpoint
  degeneracies; use arbitrary precision before declaring any sign.
- If first-order five roots are impossible in a subregion, attempt an exact
  extended-Chebyshev or argument-principle bound instead of adding compute.

### Cheap screen

`mpmath`/Arb evaluation at 128--512 bits, adaptive root isolation, and a
condition-number threshold. Reject repeated roots and roots closer to annulus
endpoints than the stated margin.

### Rigorous gate

Interval-isolate five simple zeros of `I`; map `mu` to exact rational
perturbation coefficients; then bound the full displacement
`eta I(h)+O(eta^2)` on five disjoint intervals for an explicit rational `eta`.
Finish with five CAPD Poincare fixed-point enclosures.

### Budget and stop

- Budget: 500 CPU-hours, at most 5,000,000 basis/sign evaluations, 7 wall days.
- Stop on a five-root pattern, or after the eight compact `kappa` charts
  have either interval exclusion or an explicit unresolved boundary list.
- A proof that this Melnikov space has at most three zeros is a successful
  negative result and closes the attack.

## Attack 2 — third-order weak focus plus outer separatrix

### Exact family and normalization

Use the five-parameter Shi chart

\[
\dot x=\lambda x-y+\ell x^2+mxy+y^2,
\qquad
\dot y=x+ax^2+bxy.
\]

The weak-focus-third-order stratum is `lambda=0`, `m=5a`,
`b=3 ell+5`; the historical seed is `(ell,a)=(-10,1)`. Parameterize the
unfolding by exact dyadic variables

```text
m = 5a + delta
b = 3ell + 5 - 9delta + 8epsilon
lambda < 0, epsilon < 0, delta < 0.
```

Work only in `ell in [-12,-8]`, `a in [4/5,6/5]` and logarithmic dyadic
hierarchies `2^-54 <= -delta <= 2^-6`,
`2^-270 <= -epsilon <= 2^-10(-delta)`,
`2^-900 <= -lambda <= 2^-10(-epsilon)`.

### Candidate generation

- Use exact Lyapunov quantities to require three simple local return-map roots.
- Continue the outer separatrix-splitting function in `(ell,a)` and seek a
  loop-bifurcated fourth cycle in that same nest.
- Simultaneously track the known remote cycle. Zhang's theorem demands a
  `4+1` candidate; a numerical `3+2` report is treated as a bug.

### Cheap screen

Multiprecision shooting with blown-up coordinates near the focus and analytic
saddle transition maps. Require five return roots to persist at 128, 256, and
512 bits and under halved step/tolerance.

### Rigorous gate

Exact Lyapunov-sign boxes for the three local cycles, an interval enclosure of
the separatrix splitting and resulting outer fixed point, and CAPD verification
of the persistent remote cycle. All five cycles must receive independent
section intervals and distinctness data.

### Budget and stop

- Budget: 2,000 CPU-hours, at most 200,000 continuation states, 14 wall days.
- Stop if every continued splitting branch exits the stated box, destroys a
  pre-existing cycle before creating the fourth-in-nest cycle, or reduces to a
  known graphic with a rigorous cyclicity obstruction.
- Do not widen the box in the same run.

## Attack 3 — Kuznetsov four-cycle boundary continuation

### Exact family and normalization

\[
\dot x=y+x^2+xy,
\qquad
\dot y=ax^2+bxy+cy^2+\alpha x+\beta y.
\]

Use the exact compact box

```text
a     in [-12,-8]
b     in [2,12/5]
c     in [3/5,4/5]
alpha in [-85,-65]
beta  in [2^-16,2^-6]
```

seeded at `(-10, 11/5, 7/10, -363889/5000, 3/2000)`.

### Candidate generation

- Continue the four return-map roots from the control field.
- Continue three event surfaces: a saddle separatrix splitting of zero, a
  double return-map fixed point (`D=D'=0`), and a cycle/infinity boundary in
  the Poincare compactification.
- Search only certified or high-precision neighborhoods of pairwise
  intersections of these surfaces while enforcing survival of the existing
  four cycles.

### Cheap screen

Adaptive arbitrary-precision Poincare maps on multiple sections; root count,
Floquet multipliers, equilibrium indices, and section-order checks. No raw
trajectory-count objective.

### Rigorous gate

Five CAPD interval-Newton inclusions, plus interval equilibrium classification
and either disjoint flow tubes or a proved `4+1` nesting order.

### Budget and stop

- Budget: 1,000 CPU-hours, at most 100,000 continuation states, 10 wall days.
- Stop after all seeded surface branches leave the compact box or have an
  interval sign separation proving no intersection in the retained cells.
- Do not turn this into another 32-million-point floating grid.

## Priority

Attack 1 first. It has the smallest scalar search space and can return a clean
theorem even on failure. Attack 2 is the most faithful route from the proven
four-cycle mechanism to `4+1`. Attack 3 is the best numerical-control lane but
has the highest risk of multiscale false positives.
