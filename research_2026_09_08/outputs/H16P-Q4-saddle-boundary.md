# Finite Q4 saddle-boundary check

No counterexample was found. The thirteen fitted fields with a finite saddle still show three interior cycle brackets and no additional sampled crossing near the outer return boundary. These are numerical results, not an upper-bound proof.

The exact target remains one real planar polynomial field of total degree at most two with at least five distinct isolated periodic orbits. This check keeps each field's previously fitted coefficients fixed.

## Method and findings

For each saddle, I computed its two real eigenvalues and shot both branches of the stable manifold backward and both branches of the unstable manifold forward. Each shot used a small displacement along the eigenvector. Two relative offsets, 1e-7 and 1e-9, and integration tolerances 3e-12 and 3e-13 test sensitivity to that approximation. Equilibrium residuals and the opposite eigenvalue signs were checked before shooting.

All thirteen fields have one inward stable branch and one inward unstable branch reaching the negative-Y section. The other two branches reach a numerical norm guard; this does not prove escape. There were 104 branch shots in the final run. The preliminary run failed while serializing NumPy's complex-typed eigenvalues, whose imaginary parts were zero; its results were discarded. The final script explicitly checks the imaginary parts and joins geometry records to fields by parameter keys.

The stable and unstable inward branches intersect the section at distinct radii. Their separation is at least 0.8109895 among these fields. Changing the offset and tolerance changes the stable radii by at most 4.96e-12 and the unstable radii by at most 8.82e-11. This is strong numerical evidence against a homoclinic connection at these particular parameters, not a validated manifold enclosure.

| Modulus r | Perturbation epsilon | Stable section radius | Unstable section radius |
|---:|---:|---:|---:|
| 0.5 | 4 | 4.37003605645 | 1.96289007812 |
| 1 | 8 | 2.57873961443 | 1.44215902305 |
| 2 | 8 | 1.74606900510 | 0.935079470247 |

Using the stable radius as a numerical estimate of the boundary, I sampled seventeen radii per field between 1.05 times the outer fitted cycle radius and a point only 1e-8 of that interval's width inside the boundary. All 221 paired half-return calculations resolved: 442 half-passages. No new sign-change bracket appeared. The paired displacement stayed negative at every sample.

The thirteen closest-boundary pairs were replayed at tolerance 5e-14, adding 26 half-passages. All retained their negative sign; the largest displacement change was 7.50e-15. Thus this check resolves a substantial part of the previously unsampled outer region without treating the old escape-guard failures as absent cycles.

## Scope and next test

Finite sampling cannot exclude an additional pair of zeros between sample points, and numerical shooting does not prove a global return-domain boundary. There is no validated cycle count or global exclusion here. A weak-saddle uniqueness result would require zero saddle trace; these saddles have nonzero negative trace. No applicable general bound was established in this turn's literature check.

The next useful construction test is to vary the placement of the three interior anchors and seek an outer fold of the return displacement. A verified new pair of isolated cycles coexisting with three interior cycles would meet the five-cycle target; the current samples supply no evidence that such a fold exists.

## Reproduction

The companion `q4_saddle_separatrices.py`, `q4_saddle_boundary_profile.py`, and `q4_saddle_boundary_refine.py` scripts read the previously saved finite-Q4 continuation records. Run them in that order in the audit environment. Their matching JSON files preserve all final shots, section itineraries, failed guards, coefficients by reference to the continuation records, and refined displacements. The pinned source repository was not modified.
