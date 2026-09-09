# Extreme-scale boundary search

No five-cycle counterexample was found. Two fixed near-resonant quadratic fields each produced three numerical cycle brackets, including one lower-annulus bracket near 1e249. All six brackets retained their signs at tighter tolerance. These remain nonvalidated numerical results.

The family is the same reversible unfolding used in the previous report. Set b=1, a=-1-alpha, epsilon_minus=0.01, epsilon_plus=0.01*alpha*0.001, k=sqrt[(a+2)/(-a)], e1=(epsilon_minus-epsilon_plus)/(4k), and e2=-(epsilon_minus+epsilon_plus)/2. The two fields use:

| alpha | e0 | Maximum sampled starting height |
|---:|---:|---:|
| 0.01 | -1e-100 | 1e110 |
| 0.003 | -1e-250 | 1e275 |

The saved decimal coefficients define the numerical fields actually integrated. Each field was used unchanged on both annuli.

The initial profiles used 41 starting heights per annulus: 328 half-passages. Each profile had four unresolved pairs at its extreme end. The coarse grid detected two upper crossings in each field but missed the lower crossing between its last resolved and first unresolved sample. Three adaptive midpoint pairs, six half-passages total, found the lower brackets:

- alpha=0.01: 9.502e98 to 2.268e100.
- alpha=0.003: 3.005e247 to 1.577e249.

All six upper/lower brackets were replayed at relative tolerance 3e-13, adding 24 half-passages. Their signs were preserved; the largest change in logarithmic displacement was 3.18e-12. Unresolved returns remain in the record and do not imply the absence of other cycles.

## Solver checks

The improved logarithmic solver evaluates scaled combinations through log-cosh, log-add-exp, and logistic functions. An additional positive time factor controls terms involving e0 without changing trajectories. This avoids explicitly forming very large sinh(v) or exp(-2z) values.

It passed the earlier moderate Cartesian comparisons, with maximum log-height difference 5.89e-13. It also passed twelve unperturbed large-return controls. Four additional controls started at height 1e250 on both annuli of two unperturbed fields; independent exact first-integral comparisons gave maximum log-energy error 3.60e-11. Forward/backward symmetry alone was not used as validation.

The revised solver controls used 48 logarithmic half-passages plus 24 Cartesian half-passages, and the new extreme controls used eight more logarithmic half-passages. The actual-field searches, gap checks, and refinement used 358 logarithmic half-passages.

## What remains unresolved

Large radius does not make a result a counterexample. The target is five distinct isolated cycles in one field; this search found only three numerical brackets per field. There are no interval enclosures, and these finite samples do not exclude additional roots between samples or in failed-return regions. The principal gain is resolving the same coupled boundary pattern over a much larger range of scales.
