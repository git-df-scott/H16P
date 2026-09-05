# Stage 3: two explicit trace-to-zero paths

## Result and scope

The conditioned Shi field and the published Chen–Wang visualization field each
have three numerically resolved origin-cycle brackets at the initial negative
trace. Along the sampled path toward zero, the innermost cycle shrinks toward
the focus. At trace zero, two numerical cycles remain. Both first Lyapunov
quantities are strictly positive, so the endpoint origin is an unstable
first-order weak focus. This supplies the requested classical-seed diagnostic,
but neither an M1 point nor evidence covering the four-parameter Shi chart.

There is no interval certificate here. We locate the tracked roots and their
nearby signs; we do not exclude extra roots elsewhere or a short-lived fold
between sampled parameter values. No remote cycle was recomputed.

## Exact inputs and provenance

The chart is

\[
\dot x=\lambda x-y+l x^2+mxy+y^2,\qquad
\dot y=x+a x^2+bxy.
\]

The trace at the origin is exactly \(\lambda\), and its determinant is one.

| Path | Exact \((l,m,a,b)\) | Initial trace | Intermediate trace | Endpoint |
|---|---|---|---|---|
| Shi conditioned | \((-10,499/100,1,-3113751/125000)\) | \(-10^{-14}\) | \(-10^{-16}\) | 0 |
| Chen–Wang visualization | \((-3,99/100,2/9,-3)\) | \(-2\cdot10^{-5}\) | \(-2\cdot10^{-7}\) | 0 |

The first row means \(\delta=-1/100,\epsilon=-1/10^6\). It is a newly
conditioned rational point in the Shi family, **not** Shi's original explicit
point and **not** the Galias–Tucker certified point. We did not replay the
\(10^{-13},10^{-52},10^{-200}\) hierarchy.

The equations and the Chen–Wang visualization parameters were checked against
[Yu–Zeng, equations (10), (12), (15)](https://arxiv.org/html/2002.09987v1).
The source distinguishes the classical existence proofs from its numerical
visualizations. The prior local input was `FOUR_CYCLE_FRONTIER.md` and the exact
focus derivation in `audit/claude_laneC_shi_focus.py`.

## Reproduced cycle locations

Section: positive \(x\) axis, with increasing polar angle. Entries below are
rounded numerical coordinates, not interval enclosures.

| Path/state | Innermost section coordinate | Middle | Outermost |
|---|---:|---:|---:|
| Shi, initial | 0.0000710545 | 0.000708766 | 0.0201662 |
| Shi, trace zero | collapsed into focus | 0.000712309 | 0.0201662 |
| Chen–Wang, initial | 0.0662982 | 0.165763 | 0.332839 |
| Chen–Wang, trace zero | collapsed into focus | 0.178026 | 0.332121 |

At the initial fields the stability pattern is unstable, stable, unstable.
At trace zero the surviving pattern is stable, unstable. For example, the
Chen–Wang initial Floquet multipliers are approximately
\(1.00010561,0.999356340,1.02457589\); at zero they are
\(0.999008274,1.02376443\). The Shi multipliers are extremely close to one:
the initial inner multiplier is about \(1+6.2\cdot10^{-14}\). Its closeness to
one is the expected tiny Hopf scale and is not itself evidence of a fold.

At intermediate trace, independent return signs still alternate:

| Path | Section radii, increasing | Signs of \(\log(P(r)/r)\) |
|---|---|---|
| Shi at \(-10^{-16}\) | 0.000004, 0.000015, 0.0015, 0.025 | −, +, −, + |
| Chen–Wang at \(-2\cdot10^{-7}\) | 0.004, 0.015, 0.2, 0.333 | −, +, −, + |

For the tiny Shi intermediate signs, 40-digit checks give
\(-2.1362815\cdot10^{-16}\) and \(+1.0989833\cdot10^{-15}\) at the first two
radii. At trace zero, the independently checked Shi return at \(r=10^{-6}\)
is positive, \(6.28318549\cdot10^{-18}\).

The exact first focus quantity in this Lyapunov normalization is

\[
\eta_1=-\frac{ab+2al-lm-m}{8}.
\]

It equals \(10^{-6}\) for the conditioned Shi field and \(1/400\) for
Chen–Wang. The local Hopf scaling is consequently
\(r^2\sim-\lambda/(2\eta_1)\). This explains which cycle disappears as trace
reaches zero. Numerical continuation is consistent with that local mechanism;
it does not establish a global no-fold theorem.

## Independent numerical methods

`staged_2026_09_05/shi_trace.py` constructs a rational polynomial
\(V=(x^2+y^2)/2+V_3+\cdots+V_8\) by solving the homogeneous Lyapunov equations.
It verifies the identity, exactly in rational arithmetic,

\[
\dot V\big|_{\lambda=0}
=\eta_1r^4+\eta_2r^6+\eta_3r^8+R_9(x,y).
\]

For nonzero trace the extra term is exactly \(\lambda xV_x\). An independent
DOP853 orbit integrates the polar log radius, period, divergence integral, and
\(\Delta V/V(r,0)\). This avoids subtracting two almost equal radii when
checking very small local displacement. The code never treats the truncated
Lyapunov polynomial as a globally positive function.

The independent checker `shi_mp_verify.py` integrates the direct polar
equation with 40-digit arithmetic, modified-midpoint extrapolation, seven
levels, and 32 or 64 angular blocks. It does not use the Lyapunov polynomial or
the DOP853 orbit. All four alternating initial signs in both fields were
checked with both meshes. Mesh changes in the Shi controls are below
\(3\cdot10^{-20}\); the hardest Chen–Wang outer control changes by about
\(9.6\cdot10^{-15}\), against a displacement of \(1.2\cdot10^{-5}\).
These are convergence checks, not rigorous error bounds.

The runner `shi_run_continuation.py` writes roots, return samples, periods,
Floquet diagnostics, and intermediate signs to `shi_continuation.json`.
SciPy and the multiprecision verifier were separate numerical formulations
implemented in this subtask; no claim is made that Fable executed the verifier.

## Budget, failures, and limits

`shi_returns.jsonl` contains **150 evaluations: 139 successful, 11 failed**.
Every attempted return, including controls, root finder samples, failed
probes, tighter replays, and multiprecision checks, is included. Each individual
ODE attempt had a 10 CPU-second timer. The recorded total is about 13.7 CPU
seconds and the longest recorded attempt is about 1.91 seconds. Symbolic
construction and setup are not ODE return evaluations and are outside those
runtime figures. Ten of the allotted 160 returns remain unused.

Four initial failed Shi outer probes used a log-Lyapunov diagnostic whose
denominator approached or crossed zero. The code was corrected to integrate
\(\dot V/V_{\rm initial}\), with direct polar return used when the Lyapunov
section polynomial was unsuitable. Those failed probes remain in the ledger.
The first seven successful rows retain the earlier log-Lyapunov diagnostic;
their direct polar-return values remain usable. Seven Chen–Wang probes failed
the stated polar-chart/radius guard. A guarded failure says nothing about a
cycle's existence outside that chart.

The current ledger is intentionally budget limited: rerunning a script appends
evaluations and stops at 160. To reproduce the complete recorded campaign,
copy the scripts to a fresh output directory and reconstruct the documented
control calls as well as the runner. The runner alone reproduces the trace
paths, not every exploratory failure.

## Corrections to the Stage 3 premise

Trace zero does **not** by itself require a multiple cycle if three surrounding
cycles exist. For an abstract analytic return displacement, the form

\[
D(r)=\varepsilon r^3(r^2-A)(r^2-B)(r^2-C),\qquad 0<A<B<C,
\]

has a first-order weak focus at the origin and three simple positive zeros.
For sufficiently small \(\varepsilon\), its return-map derivative is positive
on a compact section interval. This is a logical illustration for return maps,
not a construction of a quadratic vector field. A double cycle can be a
boundary mechanism producing a pair, but need not occur at the desired M1
point itself. Other boundary mechanisms also require attention.

Likewise, a fold surface having dimension three in a four-dimensional slice
does not imply random sampling finds it quickly. The exact surface has
measure zero; detection requires bracketing, a continuation method, and a
return-map domain containing the relevant cycles. Thin regions and lost return
charts can defeat naive sampling.

The next authorized mathematical work would be a separately budgeted fold
continuation away from these two trace paths, with a validated return map if
M1 is found. This subtask does not claim that systematic four-parameter fold
continuation or the Stage 5 interval verifier has been completed.
