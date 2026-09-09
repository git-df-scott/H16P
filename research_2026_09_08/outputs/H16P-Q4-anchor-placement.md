# Moving the three finite-Q4 cycle anchors

No five-cycle counterexample was found. All 24 new fields have three sampled cycle brackets; none has an additional sampled crossing. This is a finite numerical search, not an upper-bound proof.

The test fixes the finite perturbation parameter at epsilon = 4 and uses moduli r = 0.5, 1, 2. Starting from the previously fitted cycles at unperturbed energy fractions (0.25, 0.50, 0.75), I continued the three target fractions by common shifts of ±0.05, ±0.10, ±0.15, and ±0.20. Their actual section radii vary nonlinearly with the energy fractions, so this changes both placement and spacing. At each step, the three normalized controls tau/epsilon, u/epsilon, and w/epsilon were fitted using complete half-return equalities and event-corrected variational derivatives, with v = -epsilon.

All 24 fits passed the return and residual gates. The largest scaled fitting residual was 6.09e-13. The fit calculations used 1,020 sensitivity half-passages. Return profiles used another 2,460 half-passages, with zero unresolved pairs. These profiles include samples around each target radius and, where a finite saddle is present, points approaching its numerically shot stable branch. They show three sign-change brackets each. No assertion is made about unsampled zeros or tangencies.

At shifts +0.15 and +0.20, the other finite equilibrium is a node, whereas it is a saddle at +0.10. This is consistent with the finite equilibrium passing through infinity as the anchors move outward; the exact transition and any connection at infinity have not been located or verified. The six node configurations were sampled to section radius 100,000. The saddle configurations use a numerical stable-manifold boundary estimate, not a validated enclosure.

All eighteen brackets at the six continuation endpoints (shift ±0.20 for each modulus) were replayed at integration tolerance 5e-14, adding 72 half-passages. Every bracket retained its signs; the largest displacement change was 1.18e-14. The 36 stable-manifold branch shots used to guide the profiles are counted separately from the half-passages.

The intended five-cycle mechanism was an additional pair of outer return-map zeros coexisting with the three interior cycles. This run did not find that prerequisite. It also did not solve the double-zero equations for an outer fold or exhaust other anchor spacings.

The next informative continuation is toward the endpoints of this anchor family, where the inner target approaches the center or the outer target approaches the unperturbed annulus boundary. Any near-zero coefficients or vanishing displacements there must be checked for convergence to a center, rather than counted as isolated cycles.

The matching scripts and JSON records preserve every fit call, profile, equilibrium calculation, manifold shot, failure gate, and tighter replay. They use the earlier finite-Q4 continuation scripts and data in the same audit directory. No source-repository files were modified.
