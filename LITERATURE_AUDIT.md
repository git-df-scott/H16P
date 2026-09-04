# Literature audit

## 1. Exact 2026 status

Two current primary sources remove the main ambiguity. Gasull--Santana define
`H(n)` in `N union {infinity}` and state that even `H(2)<infinity` is unknown.
Villanueva--Tucker (version dated 2026-07-02) again state that the global
problem is unresolved for `n=2`, while recording Bautin's local result
`M(2)=3`. Therefore:

- lower bound: **`H(2) >= 4`**;
- finite global upper bound: **none known**;
- `H(2)=4`: **open conjecture**, not consensus theorem;
- pointwise finiteness of each fixed field: true, but logically weaker.

The distinction matters. Bamon's title, “Quadratic vector fields in the plane
have a finite number of limit cycles,” concerns an individual field. It is not
a constant bounding all quadratic fields.

## 2. Four-cycle history

### Shi Songling (1979/1980)

Shi proved four cycles for

\[
\begin{aligned}
\dot x&=\lambda x-y-10x^2+(5+\delta)xy+y^2,\\
\dot y&=x+x^2+(-25+8\varepsilon-9\delta)xy,
\end{aligned}
\]

under the hierarchy `0 < -lambda << -epsilon << -delta << 1`. The construction
uses a third-order weak focus for three small nested cycles and a global
trapping/no-contact argument for one cycle about the second focus. The original
literature and later reproductions use extremely separated explicit scales;
Yu--Zeng report `delta=-10^-13`, `epsilon=-10^-52`, `lambda=-10^-250` for
Shi's proof. Galias--Tucker later certify the nearby exact instance with
`lambda=-10^-200`.

### Chen Lan-sun--Wang Ming-shu (1979)

An independent family is

\[
\begin{aligned}
\dot x&=-\delta_2x-y-3x^2+(1-\delta_1)xy+y^2,\\
\dot y&=x+\frac{2}{9}x^2-3xy,
\end{aligned}
\]

with `0 < delta_2 << delta_1 << 1`. Their proof constructs trapping regions for
a large cycle and a small cycle, then unfolds a second-order weak focus to add
two small cycles. The original qualitative hierarchy did not give convenient
machine-ready values.

### Later fronts

- Yu--Han proved four cycles in a near-integrable quadratic family using local
  focus quantities and a global Melnikov function, again with `3+1` or `1+3`
  distribution.
- Kuznetsov--Kuznetsova--Leonov gave a visually accessible five-coefficient
  family and a decimal example with four observable cycles.
- Yu--Zeng supplied concrete numerical realizations of the historical and
  near-integrable examples, while explicitly distinguishing slow numerical
  convergence from proof.
- Galias--Tucker used adaptive arbitrary precision and interval arithmetic to
  prove that the chosen Songling instance has **exactly four** and to localize
  them. This is the strongest directly relevant rigorous-numerics result found.

## 3. The failed five-cycle episode

Shi's 1990 retrospective is unusually explicit. In 1978 he believed that
combining Sommerfeld's two-cycle mechanism with Bautin's three-cycle mechanism
gave five. On checking the preprint, he found that the published calculation of
Bautin's fifth Poincare--Lyapunov quantity had the wrong sign. Correcting it
made the required inequalities incompatible, leaving four cycles. This is a
real failed five-cycle construction, not merely an unverified internet claim.

There are two other claims that must not be imported into the campaign:

- Gaiko has published claimed proofs that four is the maximum, but current
  specialist primary sources still list `H(2)` and even its uniform finiteness
  as open; the claim is not an accepted resolution.
- A 2024 paper claimed a formula solving Hilbert's second problem (giving four
  for quadratics). Buzzi--Novaes exhibited conflicts with established lower
  growth and counterexamples to the argument. It supplies neither an accepted
  upper bound nor a five-cycle field.

Counts of five or more in **piecewise**, **discontinuous**, **cubic**, or
**higher-degree perturbing** systems do not concern `H(2)`.

## 4. Strong special-family statements

| Scope | Rigorous statement | What it does not say |
|---|---|---|
| One quadratic focus/center, local unfolding | Bautin cyclicity is exactly 3 | It does not bound large cycles in the same nest |
| Two foci | Zhang: at least one focus has at most one surrounding cycle; distributions are `(0,1)` or `(1,i)` up to exchange | It does not bound `i` |
| Neighborhood of third-order weak-focus family without polycycles | Llibre--Schlomiuk: at most four cycles in the stated neighborhood | It is not a global quadratic bound |
| Reversible center `zdot=-iz(1+a zbar)` under quadratic perturbation | Exact period-annulus bound 2 | Other centers/loops/global cycles remain |
| Generic quadratic codimension-four center `Q4` | Zhao: period-annulus cyclicity at most 5, with examples producing at least 3 | The upper value 5 is not known to be attained |
| Quadratic Hamiltonian infinitesimal problem | Gavrilov obtained a global finiteness result for the relevant Abelian-integral problem | Infinitesimal/near-Hamiltonian is not global `H(2)` |
| Fields quantitatively separated from centers/singular fields; cycles separated from equilibria and infinity | Ilyashenko--Llibre give an explicit upper estimate in terms of the separation parameters | The estimate degenerates at exactly the center, singular, polycycle, and infinity regimes where a global search is hardest |
| Quadratic algebraic limit cycles in several infinity classes | At most one algebraic limit cycle in the stated classes | Most limit cycles need not be algebraic |

For a two-focus five-cycle construction, Zhang's theorem is strategically
decisive: `3+2` is excluded. The target must be `4+1` (up to exchange), so a
fourth cycle in one nest must be genuinely global relative to Bautin's local
three.

## 5. Current computational and rigorous-numerics record

| Work | Parameterization / method | Scale and precision | Largest result | Verification |
|---|---|---|---|---|
| Kuznetsov et al. (2013) | Five coefficients in `xdot=y+x^2+xy`; analytic inequalities plus forward/backward integration | Decimal coefficients; standard numerical ODE methods | Four displayed cycles | Analytical-numerical, not an interval replay package |
| Dem'yanovich--Fefelov (2019) | Seven dimensions including two initial data; RK4 with automatic step; MPI | Two supercomputers; one scan used 32,000,000 `(b,c)` pairs | Three reliable attractors; a fourth was rejected because of floating-point rounding | Numerical only; unstable-cycle coverage incomplete |
| Yu--Zeng (2020/2021) | Shi, Chen--Wang, Kuznetsov and near-integrable families; RK4/ODE23, backward integration for unstable cycles | Small parameters and long convergence; no validated rounding | Four visualized | Numerical plus cited analytic existence arguments |
| González Prieto thesis (2021) | Proposed GPU exploration of quadratic systems | GPU-cluster design | No completed general cycle-localization solution reported | Exploratory, non-rigorous |
| Galias--Tucker (2022) | Songling point; interval arithmetic, adaptive multiprecision, global phase-space exclusion | Precision chosen dynamically; MPFR/CAPD-associated tooling | Exactly four, with rigorous localization | Computer-assisted proof |

No audited source supplies a public, standardized corpus of quadratic
coefficients, return-map brackets, and replayable interval certificates for
broad parameter search. The 2019 32-million-pair experiment is the clearest
warning against repeating a large floating-point grid: even its apparent
fourth attractor was discarded by the authors.

## 6. What the literature does and does not license

It licenses a candidate pipeline based on bifurcation equations and an interval
Poincare verifier. It does not license a finite exhaustive search of all
quadratic fields, inference from a failed sweep, or treating cyclicity bounds
near one object as a global limit-cycle bound.
