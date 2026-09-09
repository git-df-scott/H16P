# First original-Q4 compatibility test

**No five-cycle candidate was obtained.** A direct numerical calculation gives one common quadratic perturbation direction with three interior Melnikov sign crossings. Its boundary energy integral is nonzero. Thus this direction does not supply the proposed connection-preserving first-order seed for the three-interior-plus-two-boundary construction. This is a result about one seed, not an exclusion of the route.

The base, first integral, and normal controls are those in the companion original-Q4 preflight. Write
\[
h=-1+t(1-1/\sqrt2),\qquad 0<t<1,
\]
so this report's \(t\) linearly interpolates original energy. It is **not** the repository's hypergeometric variable also called \(t\).

Three base-orbit evaluations at \(t=1/4,1/2,3/4\) determine the following approximate normal-control direction, in order \((\tau,u,v,w)\):
\[
(0.000411385319236189,\;0.424367591774613,\;-1,\;0.102986332897007).
\]
The row matrix has rank three numerically: after column normalization its singular values are approximately \(1.95694,0.409972,0.0481470\). A tighter integration changed the normalized direction by at most \(1.45\times10^{-12}\).

The computed directional first-energy displacements include:

| Energy fraction | Numerical displacement coefficient |
|---:|---:|
| 0.2 | +0.0000346144 |
| 0.3 | −0.0000328739 |
| 0.4 | −0.0000603247 |
| 0.6 | +0.000134982 |
| 0.7 | +0.000174017 |
| 0.8 | −0.000501024 |
| 0.9 | −0.00414732 |
| 0.97 | −0.0142692 |

These support three interior crossings near the prescribed anchors. They are floating-point signs, not interval enclosures or a proof of simple roots. No nonzero perturbation size was selected, and no finite perturbed field was certified.

The exact rational boundary parametrization separately gives four convergent first-energy-variation integrals along its oriented finite connection. Numerical quadrature at 60 and 90 decimal digits gives, in the same control order,
\[
(3.60169731946293,\;-0.258340395577996,\;-0.0997651746319194,\;-0.222642064953825).
\]
The integrands were derived symbolically; each is a rational function times \(\sqrt2\), vanishes at both parameter endpoints, and has no denominator zero on the integration interval. Precision changes were below \(10^{-55}\), but this is still point quadrature, not validated quadrature.

Applying this boundary functional to the tighter three-anchor direction gives approximately
\[
-0.0313135213142677.
\]
This is substantially nonzero relative to the observed numerical changes. Deriving the full perturbed saddle-connection splitting theorem and the two-saddle return expansion remains a separate obligation; the finite-connection integral alone is not a boundary cycle count.

The direct orbit experiment used 15 complete unperturbed returns, including three tighter repeats, with no failed returns. Maximum recorded first-integral drift was \(1.11\times10^{-13}\), and maximum section-return error was \(1.48\times10^{-14}\). All calls, parameters, integrals, and diagnostics are preserved. The next test should impose the boundary functional while examining whether three compact roots can coexist, with independent numerical checks before any claim of exclusion or existence.

## Boundary-constrained follow-up and independent compact check

Five additional directions were computed from the frozen data, imposing the boundary integral and two compact anchors at (0.25,0.5), (0.25,0.75), (0.5,0.75), (0.75,0.9), and (0.9,0.97). No additional sampled sign crossing was located. The late-anchor profiles do not bracket every prescribed root on both sides; this is not a proof of exactly two roots or an exclusion of three. No new ODE returns were used.

An independent method solved the algebraic energy equation for the oval radius at fixed polar angles, then applied 128- and256-node Gauss-Legendre quadrature. It used neither an orbit ODE nor its event solver. At seven energy fractions between0.2and0.8, its four integrals agreed with the Cartesian orbit values to at most5.17e-13. The largest128-to256-node change was3.11e-11, at the outermost tested oval. These are independent numerical consistency checks, not rigorous quadrature errors. The angular speed was positive at every quadrature node; this sampled check alone is not a proof of a valid angular chart everywhere.
