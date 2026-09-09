# Q4 boundary variation: direct test

No counterexample found. A new variational-flow calculation supports the previously computed boundary splitting and its leading logarithmic correction. These are nonvalidated numerical checks, not a bound on limit cycles.

The base field is

\[
\dot x=-y-3x^2+2xy+y^2,\qquad
\dot y=x+x^2-4xy-y^2.
\]

Perturbations have the form
\(\delta P=\tau x+(2u+w)xy\),
\(\delta Q=\tau y+u(x^2-y^2)+vxy\).
The tested direction is the exact rational interpretation of
\((\tau,u,v,w)=(0.00041138531923618854,0.42436759177461253,-1,0.1029863328970071)\).

For each initial point \((0,-R)\), forward and backward half-passages meet the positive-y section. The differentiated event coordinate is \(s_y-(Q/P)s_x\), where the sensitivity solves \(s'=DFs+\delta F\). This includes the moving event time. Let \(D_1(R)\) be the difference of the two event derivatives.

The boundary integral predicts
\[
S_1=-0.004268103996660816.
\]
The exact unperturbed energy gap is
\[
h_\partial-h(0,-R)=\frac{3}{4\sqrt2 R^2}+O(R^{-3}).
\]
Combining this with the independently derived variation of the saddle exponent product gives the *conditional asymptotic prediction*
\[
D_1(R)=S_1+\frac{C_{\log}\log R+C_0+o(1)}{R^2},\qquad
C_{\log}=\frac{3}{2h_y(0,b_0)}(5u+5v/4+7w/4)
=0.21509807352701105.
\]
Here \(b_0=0.1594645301477301\). This calculation does not prove a parameter-uniform Dulac expansion; the formula is tested below rather than assumed as a theorem.

Writing \(T(R)=R^2(D_1(R)-S_1)\), direct flow yields:

| Radius R | [T(R)−T(R/2)]/log 2 |
|---:|---:|
| 80 | 0.2109387060 |
| 160 | 0.2161935970 |
| 320 | 0.2173714356 |
| 640 | 0.2171245135 |
| 1280 | 0.2165604979 |
| 2560 | 0.2160543073 |

All 36 half-passages succeeded at requested relative tolerances 2e−12 and 3e−14. The table uses the tighter run. Its largest unperturbed section mismatch was below 5.0e−16. Agreement of tolerances and approach toward the predicted coefficient are numerical evidence only; large-radius rescaling amplifies errors.

A separate compatibility probe imposed both first-order boundary splitting zero and first-order saddle exponent-product variation zero. Their numerical nullspace is two-dimensional. Across all 12 projective sectors determined by the saved 12 compact energy samples, at most one adjacent sampled sign change appeared. This does **not** exclude roots between samples, higher-order effects, or another modulus; neither linear constraint guarantees an exact perturbed connection. It does identify an obstacle for obtaining three compact cycles plus two boundary cycles through this particular first-order degeneration.

The next useful mathematical task is to determine whether that loss of compact oscillation holds throughout the modulus family, and to derive the parameter-dependent boundary remainder needed for a genuine coexistence test.
