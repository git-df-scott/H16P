# From five Abelian zeros to five quadratic limit cycles

## Perturbative theorem

Let \(\Pi_\varepsilon\) be the Poincaré return map of an analytic
one-parameter family of quadratic fields on a compact subannulus \(K\) of an
integrable period annulus. Parameterize the section by the first integral
\(h\), and suppose

\[
 d(h,\varepsilon)=\Pi_\varepsilon(h)-h
 =\varepsilon^k\left(M_k(h)+\varepsilon R(h,\varepsilon)\right),
\]

where \(M_1=\cdots=M_{k-1}=0\), \(M_k\not\equiv0\), and \(R\) has a uniform
\(C^1\) bound on \(K\).

If distinct \(h_1,\ldots,h_m\in\operatorname{int}K\) satisfy

\[
 M_k(h_i)=0,\qquad M_k'(h_i)\ne0,
\]

the implicit-function theorem gives, for every sufficiently small nonzero
\(\varepsilon\), distinct roots \(h_i(\varepsilon)\) of \(d\). Each is a
simple fixed point of the return map, hence an isolated hyperbolic limit
cycle. Disjoint compact neighborhoods preserve distinctness.

Multiplicity is valid in an upper-bound argument, but a multiple zero alone
does not guarantee that many real cycles.

## Application to Q4

Iliev's classification and the Gavrilov--Iliev reduction state that the first
nonzero Q4 generating function has the four-parameter form (Q4-I), with all
four coefficients independent. Consequently:

- five distinct simple zeros of a realizable \(I_\mu\) directly imply the
  existence of a sufficiently small perturbation that is still one quadratic
  vector field and has at least five cycles;
- the implication is existential until “sufficiently small” is quantified;
- an explicit counterexample also needs explicit original-coordinate
  coefficients, one explicit nonzero \(\varepsilon\), and a validated return
  map.

In the requested categories, the mathematical implication is **A**. Producing
an explicit, independently replayable counterexample operationally also
requires **C**. It is not merely B, and it is not D: the accepted Q4 upper
bound is five, not below five.

## The realization gate

The cubic Hamiltonian chart is not an affine normalization of the original
quadratic field. The double cover and inversion mean that

\[
 (\kappa,\mu)\not\mapsto
 \text{“add any quadratic }(f,g)\text{ to the Hamiltonian chart.”}
\]

Before promotion, reconstruct from Iliev's normal form the analytic
original-coordinate quadratic coefficients whose first nonzero Melnikov
function is the target \(I_\mu\). Verify symbolically that lower Melnikov
functions vanish and the first nonzero one is a nonzero scalar multiple of
\(I_\mu\). The literature proves independence of the four parameters for the
cyclicity result, but the audited formulas do not supply a drop-in inverse
routine for arbitrary input \(\mu\). This is a mandatory symbolic gate.

## Non-certificates

Five zeros of a truncated series, five double-precision sign changes, five
multiple or endpoint zeros, five zeros of an arbitrary Hamiltonian-chart
perturbation, or five plotted trajectories are only leads.
