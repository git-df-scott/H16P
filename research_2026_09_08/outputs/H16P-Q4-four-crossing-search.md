# Full-space Q4 compact-crossing search

No counterexample or four-compact-crossing candidate was found. The search used all four perturbation controls, without imposing the boundary constraints from the previous test.

At each of seven moduli (1/8, 1/4, 1/2, 1, 2, 4, 8), the saved Melnikov vectors were evaluated at 14 energy fractions between 0.02 and 0.99. For every choice of five ordered sample points, a linear program sought one control direction producing strictly alternating signs. Reversing all signs is equivalent to reversing the direction, so only one orientation is needed.

All **14,014** linear programs completed. None returned a positive margin; separate three-crossing controls succeeded at every modulus. Conditioning the four-dimensional function space before optimization reduced numerical scale imbalance. This searches sign patterns on the saved samples, not every possible energy or modulus.

An independent exact-arithmetic check found all **7,007 ordered 4×4 minors** positive when each stored decimal table entry is interpreted as a rational number. This certifies the finite-table obstruction without relying on the linear-program solver:

For any five ordered rows, let Δᵢ be the positive determinant obtained by deleting row i. The cofactor identity gives
\[
\sum_{i=0}^{4}(-1)^i\Delta_i M_i=0.
\]
If one control vector d made all five quantities \((-1)^iM_i d\) positive, taking the dot product would express zero as a sum of positive terms, a contradiction. Any four strict sign changes across the full sampled sequence would contain such a five-point subsequence.

**Scope matters:** the exact certificate applies to the frozen rational tables. Those tables approximate the actual integrals; their quadrature errors have not been rigorously enclosed. Neither this certificate nor the floating-point search proves a three-zero theorem for the true Melnikov space, and neither bounds the cycles of a finite perturbed field.

The next useful Q4 checks must target places this test cannot resolve: closely spaced energy roots, endpoint scaling, or moduli beyond the tested range. Simply testing more random control directions on these same tables cannot produce four alternating sampled signs.
