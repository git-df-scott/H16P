# Reversible Hopf precursor test

No five-cycle precursor was found in the saved first-order profiles.

For the common perturbation ε₀ = τ/(a−1), ε₁ = −cτ, ε₂ = −mτ, the normalized moment function on side σ is σ(w−mv)−c. At a center with absolute ordinate y_c, its center value vanishes when

c = σ(1/y_c − m y_c).

This is the tangent condition for Hopf, not the exact condition at nonzero τ. I imposed it separately at each center for 54 inherited shape profiles and seven new logarithmic profiles: 122 shape/center tests. Sweeping the finite sampled critical slopes found at most two nonzero crossings around the selected center. Some cases also retain one crossing on the other side, but none supplies the required three plus one precursor.

The calculation discards numerically narrow slope cells and uses finite energy grids. Missing crossings, endpoint effects and higher-order effects remain unresolved. No new ODE integration was used for this test.

A concrete next test uses the exact Hopf surface at finite perturbation. Put q = y₀−1/2 and choose rational k. Then

x₀ = −k/(2−a+bk²),  q = −k²/(2−a+bk²),  y₀ = 1/2+q,

ε₀ = 2x₀y₀,  ε₁ = −ε₂y₀−2(a−1)x₀.

Substitution gives equilibrium and trace zero because q+bq²+(2−a)x₀² = 0. The determinant, weak-focus order and coexistence of all required cycles still need checking. These equations define a testable family, not a counterexample.
