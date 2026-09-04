# Routes 4a and 4b: computations after the thought session

2026-09-04. Scripts: `audit/claude_route4b_infinity.py`,
`audit/claude_route4b_hemicycle.py`, `audit/claude_laneC_loop01.py`,
`audit/claude_route4a_normal_form.py` (final version inline in the log of this
session; the printed normal forms are reproduced below).

## Route 4b (order-3 weak focus plus a boundary graphic of cyclicity two): closed on the Shi stratum

Three independent computations on the stratum `m=5a`, `b=3l+5`:

1. **Infinity.** In the two-foci region (`3a^2<l^2+2l`) there is exactly one
   real invariant direction at infinity, hence one antipodal pair of infinite
   saddles `N` (`x\to+\infty`) and `S` (`x\to-\infty`), with hyperbolicity
   ratios `r` and `1/r` (for a degree-two field the antipode reverses time).
   Any graphic through both is therefore neutral at first order
   automatically. In the three-direction region a second saddle pair and a
   node-type pair appear.
2. **Hemicycle connection.** The unstable finite separatrix of `N` and the
   stable finite separatrix of `S` both cross `x=0` between the two foci
   (`y\approx0.01\ldots0.09`). Their splitting `D=y_N-y_S` is strictly
   positive throughout the two-foci region and decreases monotonically toward
   zero only as `l\to-\infty` (`D=+4.9\cdot10^{-4}` at `l=-30`, `a=1`;
   `+2.1\cdot10^{-2}` at `l=-6`). No connection, hence no closed graphic
   through infinity surrounding the origin, at any finite stratum point
   sampled. (Marín–Villadelprat 2025 give the cyclicity of hyperbolic
   hemicycles, with cyclicity two in the reversible `D`-system and alien
   cycles at three parameter values, but a hemicycle needs an invariant line,
   which the non-integrable stratum does not have.)
3. **Loop through `(0,1)`.** For `l>-2` the second antisaddle `(0,1)` is a
   hyperbolic saddle (trace `5a\ne0`). Both of its unstable branches escape to
   infinity (`|u|\sim10^{11}` at `t=300`) for every sampled `(l,a)`; there is
   no homoclinic loop through `(0,1)` around the origin.

Together with Lane C (finite saddle loop only on the center curve), the
computational picture is uniform: **on the third-order weak-focus stratum
no closed graphic surrounds the origin except when the origin is a center.**
Conjecture (numerical, strengthening Li–Cherkas): no separatrix cycle
surrounds a third-order weak focus of a quadratic system. If true, route 4b
is impossible, and the order-3 stratum can never give more than three cycles
around the weak focus plus one elsewhere.

## Route 4a (reversible centers): the object to study

At the loop point the origin is a reversible center. Rotating so the
symmetry axis is the `X`-axis gives the `Q_3^R` normal form

\[
 X'=-Y(1+kX),\qquad Y'=X+pX^2+qY^2,
\]

with, at `a=1` (`l_c=-1.183503419072`): `k=5.54048179`, `p=-1.24519487`,
`q=0.22849752`; at `a=2`: `k=10.28455322`, `p=-2.12799381`, `q=0.11847195`;
at `a=0.6`: `k=3.8153825`, `p=-0.97508892`, `q=0.35335404`. Finite
equilibria in these coordinates: the center `(0,0)`, the loop saddle on the
axis (`(0.803,0)` at `a=1`), and a **symmetric pair of antisaddles**
(`(-0.180,\pm0.984)` at `a=1`), one of which is the original `(0,1)` focus.

The first integral is elementary: dividing by `1+kX` gives a linear equation
for `Z=Y^2`, so
`H(X,Y)=(1+kX)^{2q/k}\,[Y^2+G(X)]` with `G` explicit. The exponent `2q/k`
(`0.0825` at `a=1`) is irrational in general: this is a reversible center
with a **non-rational Darboux first integral**, the class for which the
literature has partial results (two or three cycles in specific subfamilies)
and no general theorem (Gavrilov: "almost nothing is known about the generic
reversible case").

What a counterexample here would need: the closed period annulus of one such
center (loop included) to have cyclicity **four** under quadratic
perturbations, with simultaneously one cycle around one of the mirror foci.
Program, in the order the `Q_4` campaign has shown to work:

1. Generating functions: the first-order Melnikov function for a general
   quadratic perturbation as an integral over the ovals of `H`; the dimension
   of the space it spans; whether it vanishes identically on the reversible
   subfamily, and then the second- and third-order functions (Iliev's
   framework).
2. A Chebyshev or `\Phi`-type bound for the interior zeros, and the Dulac
   coefficients at the loop (the analogue of Lane B).
3. The mirror focus: whether the same perturbation can make it weak and
   surround it with one cycle while the annulus carries its maximum.

Expected outcome, honestly: three or fewer, as in every reversible subfamily
computed so far. But this is the only integrable family with a loop, a
second focus and no proved bound, so it is where a fifth cycle could still be
hiding, and it is fully computable.

## Status line

```text
COUNTEREXAMPLE FOUND: NO
ROUTE 4b (order-3 focus + graphic): CLOSED NUMERICALLY ON THE STRATUM
ROUTE 4a (reversible center, closed annulus cyclicity 4): OPEN, SET UP
Q4: CLOSED BY THEOREM
```
