# Direct outer-fold search in the finite Q4 family

No outer double cycle or five-cycle counterexample was found. Six searches constrained to the numerical return domain completed without meeting the double-zero equations. This is an unsuccessful bounded search, not an exclusion theorem.

Let D(R) be the difference between the forward and backward half-return heights on the upper section, starting at (0,-R). A regular outer fold would require D(R)=0 and its radial derivative D_R(R)=0 at a radius distinct from the three retained interior cycles. To produce five cycles by perturbing such a fold, one would additionally need to verify a nonzero second derivative, a transverse parameter unfolding, complete return paths, and persistence and isolation of the three interior cycles. No candidate reached even the first two requirements here.

I implemented D_R by integrating the trajectory's variation with respect to its initial radius. The initial variation is (0,-1), and the endpoint derivative includes the event-time correction: delta_y - (Q/P) delta_x. Six independent central-difference checks across three moduli showed decreasing error as the difference step shrank from 1e-3 to 1e-4; the largest absolute error at the smaller step was 1.45e-10. These checks used twelve variational half-passages and 48 ordinary half-passages.

At each outer-search evaluation, the three interior return equalities were fitted again. The outer objective was divided by epsilon times a nonzero reference displacement per unit epsilon. This reduces the risk of accepting small residuals caused only by small perturbation amplitude or a large section radius. Complete return gates remained mandatory.

The first six searches stepped outside the return domain and stopped. Their failed evaluations are retained. They used 1,566 parameter-sensitivity half-passages and 300 radial-variation half-passages. These failures say nothing about the existence of a fold.

For the revised search, the trial radius was mapped inside the numerically shot stable-saddle boundary whenever a finite saddle existed. With outer anchor a, proposed positive gap g, and stable radius S, the actual radius was a + g/[1+g/(S-a)]. Without a finite saddle, the radius was a+g. This does not constitute a validated return-domain enclosure, but it prevented the encountered domain failures in the revised attempts.

The revised search used moduli r=0.5,1,2 and two starting radius factors for each. Perturbation epsilon was bounded between 0.1 and 4 for r=0.5, and between 0.1 and 8 for r=1,2. The proposed gap divided by the outer anchor was bounded between 0.1 and 9. All six searches terminated at or near the upper bounds; none passed the 1e-6 normalized candidate threshold.

| Modulus | Best maximum absolute normalized residual |
|---:|---:|
| 0.5 | 0.03403 |
| 1 | 0.009374 |
| 2 | 0.16945 |

The revised attempts used 2,496 parameter-sensitivity half-passages and 660 radial-variation half-passages. Their stable-manifold shots are retained and counted separately in the JSON records. Optimizer convergence is not a solution of the fold equations: the nonzero residuals above are decisive for rejecting these attempts.

The current evidence does not motivate counting an extra pair of cycles in this branch. A further search would need a materially different parameter region or anchor arrangement, with explicit rejection of center limits if negative perturbations are used. The original goal remains unresolved.

Companion scripts and JSON files preserve the fitting history, domain shots, failed initial searches, normalized residuals, and derivative check. The source repository was not modified.
