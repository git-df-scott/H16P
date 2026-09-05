# KKL construction: first branch checkpoint, 2026-09-04

**No five-cycle candidate or four-cycle Hopf precursor was found.** The
first directed continuation reached the prescribed remote-section limit
before a detected origin fold. This is a result about the explored path,
not an exclusion of the KKL box. The full 4096-call allowance is not
exhausted; unvisited sheets and additional roots remain open.

The pilot and its directed extension use the council's exact family

\[
 \dot x=y+x^2+xy,\qquad
 \dot y=-10x^2+\frac{11}{5}xy+cy^2+\alpha x+\beta y.
\]

The missing object is still one beta-zero field with `K>0`, three nonzero
hyperbolic origin cycles S/U/S, and a remote U cycle. Only after those four
coexist may negative beta add the Hopf cycle. No final beta was selected,
and no precursor/five-cycle interval certificate is claimed.

## What was actually computed

The 64-evaluation pilot reproduced the four incumbent controls, the two
beta-zero incumbent origin controls, and the K>0 one-origin-plus-remote
starting controls. It then followed the two known cycles at seven rational
parameter points from `c=7/10` to `33/40` along

\[
 K=6/5,\qquad \alpha=-216/(11c-5),\qquad\beta=0.
\]

Both tracked roots remained transverse numerically. No fold was detected.
This is sampled continuation, not a validated path of periodic orbits or
proof that there are only two cycles.

At the pilot endpoint, derivatives along the root sheet gave approximately
`d mu_origin/dc|K=0.33958` and `d mu_origin/dK|c=-0.04729`.
These motivated a specified continuation to smaller K, rather than a
coefficient sweep. At fixed `c=33/40`, four steps reached `K=1/64`.
The origin root moved to about `41.65092`, with multiplier `0.922932`;
the remote U root remained present numerically. The next segment increased
c at fixed `K=1/64`, following both roots together until the remote section
range became limiting.

Selected common-parameter results, all **ordinary double-precision
numerical evidence**, are:

| c | K | Origin section root | Origin multiplier | Remote section root | Remote multiplier |
|---|---|---:|---:|---:|---:|
| 7/10 | 6/5 | 64.55543434 | 0.80969114 | -5391.14116 | 12.168019 |
| 33/40 | 6/5 | 88.46029410 | 0.86349590 | -36339.12365 | 5.228453 |
| 33/40 | 1/64 | 41.65092386 | 0.92293199 | -26821.86125 | 4.992722 |
| 9/10 | 1/64 | 45.06154708 | 0.96034550 | -237627.19148 | 2.671043 |
| 37/40 | 1/64 | 47.88456627 | 0.97355750 | -781402.78925 | 2.028553 |
| 9301/10000 | 1/64 | 48.69483537 | 0.97634279 | -1048286.51152 | 1.902808 |

The last exact field has

\[
 c=9301/10000,\quad\alpha=-8403125/209244,\quad\beta=0.
\]

Its coefficients are explicit rationals, but its two reported periodic
orbits have not been interval certified. Exact rational interval arithmetic
does certify the finite-equilibrium gates at each accepted parameter point:
one simple remote real root, stable-focus trace/determinant inequalities,
and the required K margin. The selected path stays below `c=241/250`.

## What stopped this direction

The remote root reaches the allowed section magnitude `2^20` near
`c=0.93010462269` on `K=1/64`. At the fixed section `r=-2^20`, two rational
parameter values give numerical displacements of opposite signs:

| c | D(-2^20) |
|---|---:|
| 9301046126889/10000000000000 | -0.5653180 |
| 9301046326889/10000000000000 | +0.5653573 |

These locate an experimental boundary numerically. They are not exact or
interval-certified signs, and reaching that boundary is not mathematical
escape or cycle destruction. The initial predictor failure at `149/160`
was therefore investigated rather than counted as a disappearance.

At the last admissible field, the calculated origin-root derivatives are
`mu_c|K=0.54959`, `mu_K|c=-0.07145`; the remote log-radius derivatives are
`59.7285` and `1.87493`, respectively. Locally, increasing c improves the
origin multiplier but grows the remote cycle; decreasing K improves the
origin multiplier but is prohibited by the experimental margin. These are
numerical directions, not a global optimality or absence proof.

A bounded section profile at this same field used ten coordinates between
`1/64` and `1024`. Its observed displacement signs are positive through
32 and negative from 64 onward. The known stable root lies between them.
No extra pair was detected. The intervening unsampled intervals have not
been excluded. A stationary point was located at about `r=28.17411716`,
where `R_r=1` but `D=+0.2426932516`: it is a positive displacement maximum,
**not a cycle fold**. Thus a multiplier-like derivative value of one did
not produce a false candidate.

## Numerical reliability and the repaired failure

The evaluator integrates the orbit, variational equations and divergence,
finds the opposite crossing and then the full downward return, and carries
the derivatives through the intermediate event. It records both crossings,
focus winding, coordinate ranges and the selected section component.
These are numerical itinerary diagnostics, not interval flow enclosures.

The first-derivative implementation was checked at a nonperiodic point:
at `r=60` in the starting field, finite differences give about
`0.8359517011`, the projected variational derivative gives
`0.8359516971`, while the bare divergence exponential is `0.8572672546`.
The latter is not the return derivative off a periodic point. The second
derivative from the full moving-event formula is `-0.005923254813`, agreeing
with the charged finite-difference check.

Near the large remote orbit at `c=0.9`, subtraction of large fixed-time
sensitivity components produced a discrepancy of about `4e-7`. That
continuation attempt was marked unresolved. The engine was then changed to
integrate the transverse determinants directly:

\[
 h_j=\det(F,w_j),\quad
 h_j'=\operatorname{div}(F)h_j+\det(F,F_j),\qquad
 R_j=-h_j(T)/Q_T.
\]

The scalar radial result and the independent logarithmic-divergence
formula then agreed to about `7.3e-11` at the failed control. The old
projection is retained as a cancellation diagnostic. The continuation was
replayed successfully; this numerical issue was not used to reject the
field. [Derivative derivation and independent review](kkl/notes_return_review.md).

All 202 calls, including controls and retries, are charged. They used about
6.68 CPU seconds inside the evaluators and 65.52 wall seconds including
subprocess startup; these timings exclude mathematical analysis and report
writing. Each evaluator ran serially, on one computational thread, with
a ten-CPU-second fuse. Exact totals, continuation-step counts, software
versions and remaining allowance are in
[strike_summary.json](kkl/data/strike_summary.json).

The original intermediate source hashes were not captured. Historical rows
remain unaltered, with the method change documented: calls 1–131 used the
fixed-time projection, calls 132 onward the transverse determinant. The
final reviewed source and data are hashed in [SHA256SUMS](kkl/SHA256SUMS).
This provides a reproducible final implementation, not a claim of bitwise
provenance for earlier uncommitted engine versions.

## Exact mathematical progress

On `x>-1`, a nonsingular Liénard transformation gives an exact energy
balance and a multiplier integral. The restoring force has the sign of x
throughout the experimental K>0 box. A quartic polynomial controls a
necessary amplitude condition for stable or fold cycles. At the starting
shape its sign changes prevent a naive monotonicity exclusion.
[Proofs](kkl/notes_lienard.md).

There is also a new local obstruction, proved analytically and independently
reviewed. On K=0 put `m=210/(11c-5)` and

\[
 J(c)=305+634c-11c^2-1000c^3,\quad
 \Delta(c)=\frac{4J(c)}{5(11c-5)}.
\]

The small positive-section return satisfies

\[
 D(r)=\frac{\pi\Delta(c)}{120m^{3/2}}r^5+O(r^6).
\]

At `c=33/40`, `J=103619/400>0`. Analytic dependence gives a common small
punctured neighborhood with D>0 for all sufficiently small K>=0. Therefore
the stable branch cannot collapse into the origin on this approach to
K=0; a stable/fold pair cannot be confined to that shrinking neighborhood.
This is a **local theorem with an unquantified radius**, not a finite-size
exclusion of the missing pair. Two small K=0 sign replays agree with its
sign and are explicitly controls, not order-one precursor candidates.

The exact algebra replay passed: Liénard conjugacy, divergence identity,
multiplier quartic, local quintic algebra and rational starting equilibrium
gates. [check_exact.py](kkl/check_exact.py) has a ten-second CPU fuse.
[Geometric segment proof](kkl/notes_geometry.md).

## Strongest next construction task

The missing object is a **separate finite-amplitude fold pair** coexisting
with the tracked stable origin cycle and the remote cycle. Continuing the
same multiplier-improving path beyond its agreed bounds would not resolve
that question. The next bounded task is to seed and continue additional
stationary-return branches inside the existing section/parameter limits,
using `D=0, D_r=0` with second derivatives and the Liénard amplitude
restrictions. A positive stationary maximum alone is insufficient; one
needs the S/U/S sign pattern at one shape, including the remote U root.
Disconnected sheets, additional infinity strata and unseen roots remain
unresolved. The remaining call allowance is preserved for that task.

```text
FIVE-CYCLE FIELD CERTIFIED: NO
NUMERICAL FIVE-CYCLE CANDIDATE: NO
KKL FOUR-CYCLE HOPF PRECURSOR: NOT FOUND
ORIGIN FOLD ON EXPLORED PATH: NOT FOUND
SELECTED PATH STATUS: EXPERIMENTAL RADIUS/MARGIN BOUNDARY REACHED
FULL KKL BOX: OPEN
```
