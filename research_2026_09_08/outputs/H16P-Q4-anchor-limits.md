# Finite-Q4 anchor limits

No counterexample was found. Thirty additional fitted fields approach the two endpoints of the common-energy-shift family. The calculations support three numerical crossings, with no additional sampled crossing; they do not prove a global upper bound.

The three base energy fractions (0.25, 0.50, 0.75) were shifted by ±0.225, ±0.24, ±0.2475, ±0.249, and ±0.2499, at epsilon = 4 and moduli r = 0.5, 1, 2. All thirty fits passed, with maximum scaled residual 3.16e-14. Fitting used 954 sensitivity half-passages; profiles used 3,030 ordinary half-passages. All original profile pairs resolved. Thirty stable-manifold branch shots were used separately to guide the negative-shift profiles.

As the shift approaches -0.25, the smallest imposed cycle approaches the origin and the trace approaches zero. The other two crossings remain separated. This is consistent with a Hopf limit, not evidence that an entire period annulus appears. Exact limiting degeneracy has not been established.

As the shift approaches +0.25, the largest imposed radius grows: at +0.2499 it is about 158.09, 133.23, and 97.20 for the three respective moduli. The largest radius is an input to the fitting procedure; its growth alone does not establish a limiting separatrix connection.

The automatic 1e-11 reporting threshold suppressed one crossing in each of two endpoint profiles. A dedicated replay at normal and tighter tolerance recovered opposite signs around all three expected targets in all six endpoint fields. For r = 0.5 at positive shift, the outer displacement is only about +5.98e-12 at 0.8 times the target radius and -1.77e-12 at 1.2 times it. This is numerical evidence of a very small crossing, not a certified periodic orbit. The r = 2 negative-shift inner crossing is also recovered by the wider dedicated bracket.

The replay used 216 half-passages. One deliberately wider outer sample in the r = 0.5 negative-shift field exceeded the estimated stable-manifold boundary and failed in both tolerances. That failure remains in the record. A narrower bracket, adding four half-passages, resolved with displacements +8.54e-5 and -8.14e-5.

For the positive-shift fields, I computed local data at infinity directly from the quadratic homogeneous terms. In the chart x = 1/z, y = s/z, with time rescaled by z, the radial eigenvalue is -p2(s), and the angular eigenvalue is the derivative of q2(s)-s p2(s). Each endpoint field has two saddle directions. The exponent products for the original two-saddle itinerary are approximately 2.766, 3.052, and 2.807. These are conditional local data: the global connection and the limiting parameter values have not been verified. They provide no evidence of a neutral product near these endpoints.

This continuation did not find the extra outer pair required for five cycles. The next direct construction test is to solve the outer double-zero conditions while retaining the three interior anchors, rather than infer the absence of a fold from a sampled profile. Such a solver must reject center limits, tiny residuals caused solely by large radius, and incomplete return paths.

The companion scripts and JSON records retain the fits, profile gates, local eigenvalues, tolerance comparisons, and failed wider bracket. The pinned repository was not modified.
