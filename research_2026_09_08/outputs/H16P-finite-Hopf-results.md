# Finite quadratic Hopf-family test

No five-cycle candidate was found on the eleven tested rational fields. This experiment integrates actual perturbed fields, extending the earlier first-order moment test.

I fixed a = −7/4, b = 1/3, m = 1539/1000 and used the exact Hopf parameterization recorded in the preceding report. Rational arithmetic verified the selected equilibrium, zero trace and positive determinant at every tested k. The second equilibrium and its focus classification were checked numerically.

For k = 0.0001, 0.001, 0.01, 0.05, 0.1, 0.2, 0.3 and 0.4, the sampled two-sided passage differences bracket two nonzero candidate cycles around the selected focus and one remote candidate. That is insufficient for the desired three-plus-one precursor. At k = 0.5 and 0.75, only one upper crossing and no remote crossing remained in the tested ranges; at k = 1 none were bracketed. These are sampled counts, not global counts.

The section joins the two shifted equilibria. Forward and backward passages start on the same exterior segment and must meet the interior segment between the equilibria with the prescribed crossing orientation. All 1,668 half-flow evaluations completed. Twelve endpoint differences were repeated at tighter tolerance; all signs persisted, with changes below 1.45×10⁻¹³. There is no interval validation or existence proof for all intervening initial conditions.

The exact cubic weak-focus coefficient changes sign twice on this path. Rational root isolation locates the zeros in

- [328148/779725, 6062721/14405863];
- [1292054/1849663, 1455581/2083763].

An independent derivation from the angular return equation reproduces the coefficient exactly: the quadratic radial contribution has zero mean, and the cubic return coefficient is −∫R(θ)T(θ)dθ. I have not calculated the higher focus coefficient, so these zeros are not asserted to be generic Bautin bifurcations.

The observed inner crossing moves toward the focus near the first zero. This does not establish an additional persistent cycle. Remote and very large cycles may lie outside the sampled sections, which extend at most to radius 10⁷. The next useful check is the joint continuation of the nonzero branches near these focus degeneracies, with explicit tracking of losses and remote-cycle survival.

Exact coefficients, all passage records, failed-candidate outcomes, scripts and coefficient checks are preserved in the accompanying reversible_finite_hopf and reversible_hopf_focus files.


## Two-parameter followup

An additional 16 exact Hopf fields varied k over {0.2, 0.4, 0.5, 0.7} and m over {1.4, 1.5, 1.6, 1.7}. All 2,624 numerical half passages completed. No sampled three-plus-one precursor appeared; at most two upper crossings were bracketed. These profiles reach section radius 10⁸, but still provide no exclusion beyond or between sampled points. The complete records are in reversible_hopf_shape_followup.json.

## Different reversible boundary shapes

Eight further fields used base shapes (a,b) = (−3,1) and (−9/4,1/3), with four perturbation sizes each. Their 1,312 half passages completed. One field showed two upper crossings, the other seven showed one, and none bracketed a remote crossing within radius 10⁶. No five-cycle precursor resulted.

An exact scalar reduction also gives

y″ = (1−a/2)(y′)²/y + D(y)y′ − C(y),

where D(y) = ε₁+ε₂y+(a−1)ε₀/y. On the selected Hopf surface this factors as

D(y) = −(a−1)ε₀(y−y₀)(myy₀+1)/(yy₀).

Since y′ = ε₀ on y = 0, a periodic orbit cannot cross that line when ε₀ ≠ 0. This gives a useful separation of the two nests, but the damping factorization alone does not prove a cycle-count bound. The exact reduction and new profiles are saved in reversible_lienard_exact and reversible_bicycle_hopf files.
