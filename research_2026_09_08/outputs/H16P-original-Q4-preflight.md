# Original-coordinate Q4 preflight

This is exact algebra supporting a construction attempt, not a counterexample or a cycle certificate. The base and four normal controls come from `council/notes_third_routes.md` in H16P commit `45f4ea9b4ab448bdd36702036244faa4f15c9819`. No novelty claim is made while repository reading remains incomplete.

Use the real quadratic field

\[
P=-Y-3X^2+2XY+Y^2,\qquad Q=X+X^2-4XY-Y^2.
\]

Define

\[
D=1-4Y+2(X+Y)^2,\qquad
N=6Y(1-X-Y)-1+2(X+Y)^3.
\]

Direct polynomial differentiation gives

\[
\dot D=-4XD,\qquad \dot N=-6XN.
\]

Therefore, on the open set \(D>0\),

\[
h=N/D^{3/2},\qquad
\nabla h=6D^{-5/2}(Q,-P)
\]

is a real analytic first integral, with \(h(0,0)=-1\). The last gradient identity also gives \(h=-1+3(X^2+Y^2)+O(\|(X,Y)\|^3)\), so sufficiently small regular levels are compact center ovals.

The polynomial \(C=2N^2-D^3\) has degree four, despite its degree-six presentation, and satisfies \(\dot C=-12XC\). Its degree-four homogeneous part is

\[
C_4=-12(X+Y)^2(X^2-2XY-Y^2).
\]

Thus its possible infinity directions include both saddle slopes \(-1\pm\sqrt2\) and the repeated node slope \(-1\). The candidate boundary level is \(h=-1/\sqrt2\); selecting it requires \(D>0,N<0,C=0\). Invariance and matching directions alone do not prove that a chosen real branch is the annulus boundary or establish its oriented itinerary.

For a straight perturbation with controls \((\tau,u,v,w)\), put

\[
\delta P=\tau X+(2u+w)XY,\quad
\delta Q=\tau Y+u(X^2-Y^2)+vXY.
\]

On a compact base oval, the first energy displacement is directly expressible in the original coordinates:

\[
M(h)=\oint 6D^{-5/2}(Q\delta P-P\delta Q)\,dt
=\oint 6D^{-5/2}(\delta P\,dY-\delta Q\,dX).
\]

For an oval bounding a domain contained in \(D>0\), traversed counterclockwise, Green's theorem writes this as

\[
M(h)=\iint 3D^{-7/2}
\left[2D\,\operatorname{div}\delta F-5\nabla D\cdot\delta F\right]\,dX\,dY.
\]

The enclosed domain condition matters: a singular integrating factor cannot be passed through silently. These formulas provide four linear perturbation integrals without treating arbitrary perturbations of a transformed Hamiltonian chart as admissible original quadratic fields.

The companion script checks the polynomial identities and records the four area-density numerators using exact SymPy arithmetic. All assertions pass; no ODE evaluations were performed. The next test is to identify the boundary branch and compute the four original compact return derivatives at three anchors, then check whether their common zero direction can satisfy the boundary conditions. No such compatibility or five-cycle conclusion has yet been established.
