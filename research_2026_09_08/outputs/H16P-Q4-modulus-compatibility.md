# Q4 modulus compatibility test

No counterexample found. The obstruction observed at modulus 1 persists in a bounded numerical test at moduli 1/8, 1/4, 1/2, 1, 2, 4, and 8. This is evidence about a particular first-order construction, not a theorem excluding five cycles.

For the original Q4 base field
\[
P=-y-(2+r^2)x^2+2rxy+y^2,\qquad
Q=x+rx^2-(1+3r^2)xy-ry^2,
\]
use perturbations \(\delta P=\tau x+(2u+w)xy\), \(\delta Q=\tau y+u(x^2-y^2)+vxy\), with \(r>0\).

An exact symbolic calculation gives the first variation of the logarithm of the physical two-saddle exponent product:
\[
\delta\log\rho=
\frac{20r^2u+5rv+(6r^2+1)w}{(1+r^2)^{3/2}}.
\]
This reduces to the previously checked expression at r=1. The calculation differentiates the saddle locations as well as their eigenvalues. It concerns local saddle data only.

At each tested modulus, I imposed two linear conditions: the boundary energy integral vanishes, and the displayed exponent variation vanishes. The resulting numerical two-dimensional control space was tested at 14 energy fractions from 0.02 through 0.99. For each quadrature run, all projective sectors determined by these sample functionals were inspected. At most **one adjacent sampled sign change** appeared in every case.

All 196 initial compact quadratures completed. The 256/512-node comparison became inadequate at large modulus, with maximum absolute changes up to 0.000297. Refinement at r=2, 4, and 8 used another 84 compact quadratures at 1024/2048 nodes and retained the same crossing count. The largest componentwise change divided by 1+absolute component value was respectively 1.05e−13, 1.15e−13, and 2.30e−10. These comparisons are not interval enclosures.

The result makes this specific first-order route to three compact cycles plus two boundary cycles less promising. It does not rule out unsampled compact roots, other moduli, higher-order unfoldings, or mechanisms without these simultaneous first-order constraints. Exact perturbed connections and exponents have not been constructed.

A useful next test is to release these constraints and directly seek four compact sign alternations across the full perturbation space, rather than relying only on three prescribed anchors.
