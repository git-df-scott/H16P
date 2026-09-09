# Implementation findings

These findings concern inherited numerical checks, not counterexamples to their analytic theorems. Original repository files remain unchanged.

## Hypergeometric parameters rounded before precision is raised

`audit/claude_green_tools.py` constructs `one6` and `five6` at module import. Several callers import this helper before setting `mp.mp.dps` to 30–45. Raising precision afterward does not restore the exact fractions 1/6 and 5/6.

A fresh-process reproduction imports at 15 digits and evaluates at 60 digits. The stored parameters differ from fresh high-precision fractions by approximately −9.25e−18 and +3.70e−17. At t=0.99, the helper's F differs from the freshly parameterized hypergeometric function by 2.3984e−17. Thus output with 40 printed digits is not evidence of 40 accurate digits for this helper. The effect on cancellation-sensitive derived quantities requires separate evaluation; the observed discrepancies do not invalidate the analytic proofs.

The appropriate repair in future independent evaluations is to construct rational parameters inside each function at the current precision, or set precision before import and preserve sufficient guard digits. Our direct original-coordinate Q4 probes do not import this helper.

## A claimed wrong-sign control is only a source-text assertion

`audit/claude_check_reconstruction.py` has a comment promising to flip the reconstruction forcing sign and demonstrate failure. Its actual final check merely asserts that the corrected forcing assignment appears in the source text. That script does compare the corrected reconstruction with area evaluations, but does not execute the advertised wrong-sign negative control. Other repository records may contain separate negative controls; this finding applies to this script alone.

## Endpoint quadrature checks truncate singular tails

`claude_check_endpoint_identities.py` and `claude_check_theoremN.py` explicitly replace a small transformed endpoint region by zero and accept finite numerical tolerances. Their moment checks are numerical consistency tests, not rigorous quadrature enclosures. This is consistent with retaining the separate analytic beta-moment proofs as the actual evidence.

## Early separatrix probes must not be treated as connection equations

The early `claude_laneC_splitting.py` measures a dot product with the unstable eigenvector. Stable and unstable eigenvectors need not be orthogonal, so this is not generally the coefficient in the unstable direction. Its successors `splitting2.py` and `splitting3.py` compare with a stable tangent at a finite distance, which also omits stable-manifold curvature. `splitting4.py` improves this by integrating the stable branch backward and comparing the two crossings. Its finite initial offset still needs error control for a certificate. These are successive experimental versions, not independent confirmations of the same connection.

In `claude_route4b_hemicycle.py`, the documented chart is z=1/x and the finite separatrix tangent is u=u0+v_z z+O(z²). The initializer uses z=+1/R for both x=+R and x=−R. On the negative-x branch, the first-order correction should have the opposite sign in that chart. This is a specific initialization defect in that historical probe, not evidence of a hemicycle or a five-cycle construction. Any new use of this route should derive the signed chart expansion afresh and control the initialization error.
