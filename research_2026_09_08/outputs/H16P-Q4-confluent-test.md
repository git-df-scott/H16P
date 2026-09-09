# Q4 clustered-root test

No counterexample or fourth-order compact-zero candidate was found. This test addresses a gap in the previous coarse-grid search: several roots could lie close together between its sample points.

I integrated the four columns of the repository's Picard–Fuchs reconstruction simultaneously. In its coordinates, the free coefficients are `(A,B,eta,c)`, with forcing
\[
q(t)=A+Bt-\eta M(t)+c(tM(t)-1).
\]
The state columns are `(H,Y,V,X)`, where V=Y′ and X′=(1−at)^(−3/2)Y. A nonzero combination with a fourth-order zero requires the four-function Wronskian to vanish. Its zeros can be tested without numerically differentiating the integrals:
\[
W_X=f^3 c_H\det[X,Y,V,H],\quad
f=(1-at)^{-3/2},\quad
c_H=-\frac1{1152t^2(1-at)(1-t)^2}.
\]
Here c_H is an ODE coefficient, distinct from the free control c.

The normalized state determinant was negative at all 18 tested energies from 0.001 through 0.999999, for each of eight moduli r=1/32, 1/8, 1/2, 1, 2, 8, 32, 128. All 16 integrations succeeded at two tolerances. No sampled sign reversal was found. At large modulus near the endpoint the matrices become nearly singular, so agreement between tolerances is not a rigorous sign certificate.

There is also an exact local calculation. Expanding the reconstructed columns at the center in the ordered basis `(A,B,eta,c)` gives
\[
W_X(t)=\frac{t^4}{15416885772288}+O(t^5).
\]
The leading coefficient is **independent of a**. It follows that, for each fixed admissible a, this Wronskian is positive sufficiently close to the center. This rules out a fourth-order zero of a nonzero combination in that sufficiently small neighborhood. No uniform neighborhood in a, global Wronskian theorem, or finite-perturbation cycle bound is claimed.

As an independent normalization check, all four reconstructed columns were compared with the original area-integral evaluator at kappa=2 and t=0.2, 0.6, 0.9. The maximum absolute difference was 1.01e−16. The area quadrature is itself nonvalidated and retains the inherited fixed-bisection and endpoint-clipping limitations.

These results make an ordinary compact confluent mechanism less promising in the tested region. The unresolved cases include parameter-dependent endpoint scaling, untested intervals, and higher-order perturbations. A counterexample still requires one exact quadratic field and verified isolated cycles.
