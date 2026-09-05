# Independent review of the staged return implementation

2026-09-05. Read-only mathematical and code review of `compact_return.py`,
`run_kkl.py`, `refine_kkl.py`, the saved control/derivative/path/fold records,
and `STAGED_INFINITY_2026_09_05.md`. A final read also checked the added
orientation guards, `cartesian_check.py`, and ledger entries 392–400.
This review performed **zero ODE
evaluations**. It is an architecture review, not an interval certificate or
an independent numerical reproduction.

**Verdict:** the transformation, divergence integral, section-flux derivative,
bounded-coordinate chain rule, and augmented fold equations are correct on a
valid transverse return branch. The saved controls are consistent with those
identities. No mathematical formula bug was found that invalidates the saved
controls. Acceptance remains numerical and has the limitations below. In
particular, no claim of a global KKL exclusion follows from these scripts.

## 1. Log-polar equations and the two clocks

Write \(\rho=e^w\), \(C=\cos\theta\), \(S=\sin\theta\), and

\[
P=\rho S+\rho^2(C^2+CS),\qquad
Q=\rho(\alpha C+\beta S)+\rho^2(-10C^2+\tfrac{11}{5}CS+cS^2).
\]

With \(dt/d\tau=(1+\rho)^{-1}\), the code's `p` and `q` are
\(P/[\rho(1+\rho)]\) and \(Q/[\rho(1+\rho)]\). Therefore
\(w_\tau=Cp+Sq\) and \(\theta_\tau=Cq-Sp\), exactly as implemented.
The stable logistic evaluation avoids directly forming very large \(\rho\)
inside these equations.

The physical divergence is

\[
\operatorname{div}F=\tfrac{21}{5}x+(1+2c)y+\beta.
\]

Its integral variable correctly has derivative
\(L_\tau=[\rho((21/5)C+(1+2c)S)+\beta]/(1+\rho)\).
The fourth state correctly accumulates physical time
\(T_\tau=(1+\rho)^{-1}\). Consequently `period` means physical time and
`rescaled_period` means elapsed numerical desingularized time. Neither can
be substituted for the other when classifying infinity behavior.

The implementation extends finite-radius integration beyond the previous
Cartesian cutoff. It does not integrate through infinity: \(w\) remains an
unbounded coordinate and has a finite guard, while bounded \(q\) is calculated
on the section afterward. A guard termination is correctly `UNRESOLVED`.

## 2. Curved section, orientation, and event sequence

The section is \(\gamma(r)=(r,h(r))\),
\(h(r)=-r^2/(1+r)\), where \(P=0\). The event expression is
\(P/\rho^2=C^2+CS+Se^{-w}\), so it has exactly the correct zero set and
crossing direction at a transverse event. On that set,

\[
\frac{dP}{dt}=(1+x)Q.
\]

For the origin-side maximum, \(x>0,Q<0\), giving the desired negative event
direction. For the remote-side minimum, \(x<-1,Q<0\), giving positive
direction. Searching first for the opposite event and then the desired event
is the appropriate way to skip the starting section point and perform the
full excursion, provided these branches and transversality are established.

The separate transport map is not a full return. The origin transport seeks
the nearby desired event; the remote transport follows the opposite/desired
sequence and can travel most of a turn. This difference agrees with the saved
transport winding records. Transport's `D` and `D_r` are coordinate differences
between different sections, not cycle displacement tests.

The initial code inferred the branch from `r>0` or `r<-1` without explicitly
checking initial \(Q<0\), and did not assert matching nonzero final flux.
For an arbitrary remote input, \(r<-1\) alone does not imply the input is left
of the remote equilibrium or on the chosen minimum branch. I recommended
explicit initial/final \(Q<0\) checks and \(r<x_*\) for that curved remote
branch. **The final source now implements all three checks:** a unique remote
equilibrium with \(r<x_*\), initial \(Q<0\), and final \(Q<0\).
I checked the actual source; the guarded polar evaluator was used successfully
in ledger entry 398. This fixes the unchecked-orientation issue for accepted
inputs, subject to the same floating-point caveats as the rest of the evaluator.

## 3. Flux derivative and coordinate changes

For a transverse section parameterized by \(r\), variation of flow gives

\[
R'(r)=\frac{\det(\gamma'(r),F(\gamma(r)))}
{\det(\gamma'(R),F(\gamma(R)))}\exp L.
\]

The determinant on the curved section is \(Q-h'P=Q\); on the horizontal
section it is also \(Q\). Thus `log_flux` correctly computes
\(\log|Q|=2\log\rho+\log|Q/\rho^2|\), and the derivative formula is correct
when the section orientations match. The absolute value would silently discard
a negative flux ratio; the final explicit initial/final sign guards now check
that assumption. No derivative of \(h\) is missing, since its term is
multiplied by \(P=0\) at the curved endpoint.

For either same-sign branch, \(q(r)=1/(1+|r|)\) gives

\[
\frac{dq(R)}{dq(r)}=R'(r)\left(\frac{q(R)}{q(r)}\right)^2.
\]

Both `D_q_derivative` and `log_displacement_derivative` follow correctly:

\[
\frac{d}{d\log|r|}\log|R/r|=\frac{rR'}R-1.
\]

At a fixed point the flux ratio is one and \(\exp L=R'\), the transverse
Floquet multiplier. **Off a fixed point, \(\exp L\) is only a divergence
factor**, even though the generic result dictionary retains `multiplier`.
The final source also returns `divergence_exponential` and the explicit
interpretation string “Floquet multiplier only at a fixed point”.
The report should call it a multiplier only for accepted periodic-orbit
approximations. The same warning applies to transport rows.

The saved derivative validation compares the analytic logarithmic derivative
with centered finite differences at the ordinary origin control and cutoff
remote control. The differences are approximately \(2.24\times10^{-8}\) and
\(-6.65\times10^{-9}\). This is useful numerical support for the derivative
and branch selection. It is not a rigorous error bound and does not validate
every parameter or large-radius return.

## 4. Augmented fold equations and acceptance

For finite same-sign section points, let \(z=\log|r|\) and
\(F_1=\log|R/r|\), \(F_2=rR'/R-1\). Then

\[
F_1=F_2=0\quad\Longleftrightarrow\quad R=r,\ R'=1.
\]

This also agrees with the bounded-section fixed-point and derivative equations.
The Newton system in \((z,c)\), with \(K\) fixed, therefore targets the correct
event. It need not have the same off-root stationary points as the original
displacement: \(R'=1\) alone does not imply \(F_2=0\) unless \(R=r\).
That distinction does not invalidate a seed, but the seed must not be reported
as an already stationary point for both coordinate choices.

The forward-difference Jacobian and step damping are a numerical discovery
method. There is no proof of conditioning, convergence, or basin coverage.
Four iterations per seed and six \(K\)-values are a bounded probe; they do not
cover a continuous \(K\)-interval or disconnected fold sheets. An ordinary
fold still needs a verified nonzero second radius derivative and an effective
parameter direction, followed by nearby cycle-count checks.

The initial acceptance uses small absolute logarithmic residuals. These can
be small for a nonroot close to a weak focus, where the return is nearly the
identity. `refine_kkl.py` divides both equations by
\((1-q)^2=(|r|/(1+|r|))^2\), removing their automatic small-radius quadratic
vanishing without altering finite roots. This is an appropriate numerical
improvement. It still does not turn a tolerance test into an existence proof.
The saved initial fold attempts all ended at `ITERATION_LIMIT`, and none met
their candidate gates; no false fold was asserted there. Subsequent refinement
did find a numerical fold candidate, reviewed in section 7 below; the initial
negative search outcome is therefore superseded for fold discovery.

`correct()` accepts a root approximation from a residual gate before proving
a bracket or uniqueness. Near \(R'=1\), a small residual can mean a much larger
radius uncertainty. Approximate radii, periods, and multipliers should retain
that evidence level. Tightened-tolerance replay and nearby displacement
brackets can support the root interpretation; interval existence requires a
separate validated integrator and section argument.

## 5. What the guards and sampled paths establish

The radius, inward-radius, physical-time, numerical-time and CPU fuses bound
the run, and failures are explicitly unresolved. The numerical winding check
is useful for rejecting obvious half returns or wrong itineraries, but uses
unwrapped angles on solver nodes. It is not a proof of winding when angular
sampling about the remote focus is insufficiently resolved.

`min_xy`, `max_xy`, and `max_log_radius` use accepted solver nodes. They are
sampled extrema, **not enclosures of the entire orbit**. In particular, these
values cannot certify the new exact theory requirement involving the maximum
of \(x\), separation of annuli, avoidance of equilibria, or noncrossing by
themselves. The exact one-way barrier \(P(-1,y)=1\) is separate analytic
support: an exact remote trajectory returning to \(x<-1\) cannot have crossed
to its right and subsequently crossed back.

The code parses rational requests into binary floats. Discovery rows therefore
are numerical approximations to the specified rational fields. A final
coefficient-alone interval verifier must use the rational coefficients
directly and not treat the float parser as exact arithmetic.

## 6. Infinity report cross-check

The reported chart equations, reversal of physical orientation at antipodes,
vertical separatrix side argument, and distinction between local eigenvalue
neutrality and a connected graphic are mutually consistent. The candidate
\(J(c)=0\) line is explicitly only an eigenvalue-neutrality locus. It does not
compute the global connection or its transition constant. I found no clear
algebraic or orientation flaw in these statements on inspection.

The warning about periods and multipliers is necessary: a physical time factor
vanishing at infinity can leave a finite physical period despite divergent
desingularized time. A nonunit multiplier can also occur near a graphic.
Neither observation alone supplies a terminal-mechanism diagnosis.

This review did not independently classify every possible graphic, especially
at nonhyperbolic infinity strata. The exact exclusion of the specified
vertical connection below \(c=241/250\) should not be expanded into an
unproved global statement about all bifurcations.

## 7. Final numerical evidence and amended verdict

The final ledger contains exactly **400 charged KKL evaluations**. Entry 392
records the refined fold candidate at

\[
 K=1/512,\quad c\simeq0.9688884793906646,\quad
 r\simeq6.949087993605231.
\]

Its polar residuals are \(D\simeq-9.59\times10^{-14}\) and
\(D_r\simeq6.66\times10^{-16}\). Entry 393 uses a separately implemented
original-time Cartesian return map, with \(D\simeq-1.40\times10^{-10}\) and
\(\exp L\simeq0.999999999943145\). I checked that evaluator's section event,
physical divergence integral, and flux derivative; they have the correct
form. Its method remains floating-point DOP853, so it provides independent
equation/coordinate implementation support, not a different interval solver.

At \(c=0.9688894793906646\), the independent Cartesian entries 394–396 give:

| Section coordinate | Numerical displacement |
|---|---:|
| 4 | \(+6.8122806\times10^{-6}\) |
| 6.949087993605231 | \(-1.1352691\times10^{-6}\) |
| 12 | \(+6.2132990\times10^{-5}\) |

This is numerical evidence for the expected stable/unstable pair on the
fold's cycle side. Entry 397 at the opposite parameter displacement gives
positive middle displacement \(1.1349862\times10^{-6}\). An interval sign
certificate and ordinary-fold nondegeneracy proof are still required to make
these existence statements rigorous. The surviving outer-cycle premise
needed for M1 has not been supplied: the guarded large-radius profile in
entry 398 is positive, and one such profile neither finds nor excludes an
additional outer stable cycle.

Entries 399 and 400 independently replay the terminal origin and remote
approximations at \(c=0.9683,K=1/64\), with relative logarithmic displacements
about \(-1.18\times10^{-9}\) and \(3.55\times10^{-9}\). These support the
reported continuation far beyond the old radius cutoff. Their finite residuals
and nonvalidated error estimates must remain visible in the evidence level.

The amended combined verdict is a supported **numerical continuation past the
old artificial radius cutoff**, a **numerical fold candidate with an
independently implemented pair sign pattern**, and an exact candidate
eigenvalue-neutrality line whose connection and coefficient remain open.
M1 and a five-cycle interval certificate remain absent. The Stage 2 global
kill condition is not established, and a statement that no fold was found
would now be incorrect.
