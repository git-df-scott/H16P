# Finite-field degenerate-Hopf precursor search

**Search-direction correction:** Published focus nonexistence/uniqueness results obstruct the proposed exact-focus precursor. The numerical record below is retained, but its suggested continuation is superseded by [the correction](H16P-focus-route-correction.md).

No five-cycle candidate was found. This search tested a concrete coexistence requirement: two outer cycles around a cubic-degenerate Hopf focus, plus a cycle around another focus, would be a possible starting point for adding two small cycles. The required two outer crossings did not appear in the tested fields.

The family is
\[
\dot x=(b-2)/4+(1-b)y+ax^2+by^2+\epsilon_1x+\epsilon_2xy,
\qquad \dot y=\epsilon_0-2xy.
\]
For rational a,b,k, the selected equilibrium and parameters are
\[
x_0=-k/(2-a+bk^2),\quad y_0=1/2+kx_0,
\quad\epsilon_0=2x_0y_0,
\quad\epsilon_2=-(a-1)\epsilon_0m,
\quad\epsilon_1=-\epsilon_2y_0-2(a-1)x_0.
\]
Equilibrium and zero trace are checked exactly. The cubic Hopf coefficient was reduced symbolically to a polynomial equation in m. Real roots with positive determinant were approximated by rational numbers with about 40 decimal digits. Consequently these are **near-degenerate rational fields**, not certified Bautin points. The largest substituted polynomial residual was 9.02e−40; no second Lyapunov coefficient or generic unfolding condition was verified.

The bounded shape set used a=−3/4,−5/4,−7/4,−5/2; b=1/3,1,5/3; and k=1/10,2/5. Of the resulting near-degenerate fields, 27 passed the numerical two-real-equilibrium gate. Another 21 failed that geometry gate and were retained in the records; this does not exclude other mechanisms in those fields.

For each accepted field, forward and backward half-passages from 33 section radii between 0.03 and 10^6 were checked against the between-equilibria section. All **1,782 half-passages** resolved. Six fields had one sampled outer crossing; 21 had none. None had two, so the conditional search around the other equilibrium was not triggered.

All six brackets retained their endpoint signs in **24 tighter half-passages**. The largest change in endpoint separation was 9.03e−14. This supports the numerical brackets, but does not establish rigorous orbit existence, isolation, or completeness between sample points.

This specific bounded precursor search failed. It leaves open tightly clustered cycles, cycles beyond the sampled range, other shapes, and other bifurcation mechanisms. The next useful step is a targeted continuation toward an outer-cycle fold on the cubic-degeneracy surface, rather than declaring the surface free of the needed coexistence.
