# Three-center Q4 compatibility test

**No five-cycle candidate was found.** The same three interior energy anchors were tested at center moduli \(r=1/2,1,2\). Each produced three numerical sign crossings, but its boundary energy integral remained nonzero. This does not exclude other anchors, other parameter directions, or any interval of moduli.

The original quadratic base is
\[
P=-Y-(2+r^2)X^2+2rXY+Y^2,\qquad
Q=X+rX^2-(1+3r^2)XY-rY^2.
\]
With \(d=1+r^2\), the exact first integral is \(h=N/D^{3/2}\), where
\[
D=1-2dY+d(rX+Y)^2,\qquad
N=3dY(1-rX-Y)-1+d(rX+Y)^3.
\]
The gradient identity is \(\nabla h=3r^2dD^{-5/2}(Q,-P)\). The factor \(r^2\) is essential. The center and boundary energies are respectively \(-1\) and \(-1/\sqrt d\).

For each modulus, anchors were placed at energy fractions \(1/4,1/2,3/4\), using \(h=-1+t(1-1/\sqrt d)\). Their common null direction uses the same normal controls as the earlier test, ordered \((\tau,u,v,w)\).

| Modulus | Approximate common direction | Directional boundary integral |
|---:|---|---:|
| 1/2 | (0.000390235, 0.764492, −1, 0.0687110) | −0.00888140549 |
| 1 | (0.000411385, 0.424368, −1, 0.102986) | −0.03131352131 |
| 2 | (0.000213722, 0.296941, −1, 0.110389) | −0.04725238365 |

At energy fractions0.2,0.3,0.6,0.8, all three directions gave successive signs \(+,-,+,-\). These are numerical signs, not validated enclosures or a proof of simple roots. Each test imposes three compact conditions on one common direction; it does not mix cycles from different fields.

Compact integrals were computed by solving the algebraic oval equation along polar rays and applying Gauss-Legendre quadrature. The boundary integral used the exact rational connection parametrization. Increasing both quadratures from256to512nodes changed the anchor matrices by at most \(9.77\times10^{-15}\), the normalized directions by at most \(4.11\times10^{-14}\), and the boundary vectors by at most \(2.05\times10^{-14}\). These are convergence diagnostics, not rigorous error estimates. The modulus-one result also agrees with the separately computed Cartesian orbit integrals and high-precision boundary quadrature.

The script and JSON preserve the base parameters, matrices, directions, sign profiles, and residual diagnostics. No new ODE returns were used. No nonzero perturbation size or finite perturbed field has been certified. The missing work remains a compatible common perturbation and a full two-saddle return analysis with controlled remainders.
