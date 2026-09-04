# Claude adversarial audit of Astra Strike 4 (Theorem N)

Date: 2026-09-04. Audited commit: `46bca95` (merged as `5bcfe11`).
Documents: [ASTRA_FOURTH_STRIKE.md](ASTRA_FOURTH_STRIKE.md),
[Q4_THEOREM_N.md](Q4_THEOREM_N.md), `q4/notes_N_*.md`,
`q4/check_N_kernel.py`, `q4/q4_N_loop_checks.py`.
Independent checker: [`audit/claude_check_theoremN.py`](audit/claude_check_theoremN.py).

## Verdict

```text
THEOREM N PROOF SOUND: YES
FIVE DISTINCT Q4 INTERIOR ZEROS EXCLUDED (all kappa>1): YES, PROVED
Q4 INTERIOR BOUND ON THE STRICT LOBE REGION: 3 distinct (proved)
Q4 INTERIOR BOUND GLOBALLY: 4 distinct (proved); 3 remains open
ATTACK 1 (Q4 FIVE-ZERO ROUTE): CLOSED NEGATIVELY
CLAUDE HANDOFF CORRECTED: YES (two items, below)
```

## The proof, reconstructed

Theorem N states `Phi_a(tau_1)<0` for every `kappa>1` and every point of the
strict lobe region, where `Phi_a(t)=Y_0+\int_0^t W_aH`, `W_a=Rcal_a\Omega_a`.
Astra proves it by three strict global comparisons and one exact identity:

1. **(N0)** `\int_0^1W_1H_*\,dt=3/1232=-Y_{0,*}`, with
   `W_1=\frac{3}{2304t^2}[(1-t)^{-17/6}-(1-t)^{-13/6}]`. This is the second
   strike's pair of beta moments (`1` and `25`) in disguise; I verified it
   by singularity-removing quadrature to the precision of that quadrature
   (`2.435022e-3` versus `2.435065e-3`; the residual is quadrature error on
   the `(1-t)^{-5/6}\log` singularity, and the identity is proved exactly by
   the positive-series argument already audited).
2. **(N5)** `Y_0<Y_{0,*}` and `0<H(t)<H_*(t)` on `(0,tau_1)`. Mechanism:
   with the `K_3` coefficient fixed at one, `\partial_{y_j}H=-H'(y_j)\ell_j`
   where `\ell_j` is the cardinal interpolant in `span\{K_0,K_1,K_2\}`
   (a space with at most two interior zeros, by the ECT property of
   `\{1,u,M\}` and anchored Rolle). Sign bookkeeping gives
   `\partial_{y_j}H(t)>0` for `t<y_1` and all `j`. For `Y_0` the same
   derivative has the sign of the center functional of `\ell_j`, and the new
   convexity lemma (N3) shows that a two-root variation and its center
   functional share their initial sign. Pushing all anchors to one along
   `y_j(\theta)=1-(1-\theta)(1-y_j)` (fixed ratios, so the audited corner
   limit applies) gives the strict comparisons.
3. **(N3)/(N2)** The moment curve `x=K_1/K_0`, `m=K_2/K_0` is strictly convex
   (`S'=(M'-S)/(t-x)>0` because `S` is a weighted average of secant slopes
   of the convex `M`), with `x(1)=6289/9061`, `m(1)=11/41`, terminal slope
   `1105/462`. For a two-root variation `\alpha K_0+\beta K_1+\gamma K_2`,
   `\gamma>0`: `\alpha>-\gamma/6`, `\beta>-(1105/462)\gamma`, hence the
   center functional exceeds `(9/3080)\gamma\cdot25/231>0`.
4. **(N6)/(N7)** `W_a<W_1` for all `a<1`, `t\in(0,1)`: `v_a=z_a/(1-at)^{3/2}`
   solves `L_av=0` with `L_a=(1-at)(1-t)D^2-\frac{1+5a-6at}{2}D+\frac{8a}{9}`,
   and the `a=1` solution `v_1=\frac32[(1-t)^{-4/3}-(1-t)^{-2/3}]` has residual
   `L_av_1=(1-a)(22-7(1-t)^{2/3})/(6(1-t)^{7/3})>0`; identical initial data
   and a positive causal Green function (a positive homogeneous solution
   exists) give `v_a<v_1`.
5. **(N8)** `\Phi_a(y_1)<Y_{0,*}+\int_0^{y_1}W_1H_*=-\int_{y_1}^1W_1H_*<0`.

## What I checked

| Item | Method | Result |
|---|---|---|
| (N0) | quadrature with `t=1-u^6` | `3/1232` to quadrature precision |
| `x(1)`, `m(1)` | closed moments | exact rationals reproduced to 20 digits |
| `S` increasing, `m` convex in `x` | 99-point grid, second differences | holds |
| transformed operator `L_a` | symbolic conjugation from the `y`-equation | matches |
| residual `L_av_1` | symbolic | matches `(1-a)(22-7d^{2/3})/(6d^{7/3})` |
| `W_a<W_1` | sampled `kappa=1.5,4,20,500`, `t=0.2,0.7,0.99` | holds |
| (N8) on Astra's tuned shots | `Phi(tau_1)` versus `-\int_{tau_1}^1W_1H_*` | holds with margin (`-3.93e-3<-2.43e-3` etc.) |
| Astra's replays | `check_N_kernel.py`, `q4_N_loop_checks.py` | exit 0 |
| Cardinal-interpolant signs | re-derived (Cramer, `x` increasing, `det>0`) | `sign(\gamma_j)=(-1)^{j-1}` |
| Anchored ECT for `span\{K_0,K_1,K_2\}` | Rolle from `\{1,u,M\}` | at most two interior zeros |
| Use of the corner limit | fixed-ratio path; limit universal for any strict ratios | inherited theorem applies |

I looked specifically for the gaps that were plausible: a hidden uniform
asymptotic (none: only the pointwise fixed-ratio limit is used, as a
comparison destination); a non-strict inequality surviving the limit
(strictness comes from any positive path parameter plus the positive tail
`\int_{y_1}^1W_1H_*`); the sign of the `K_2` coefficient of `\ell_j`
(re-derived); the direction of the kernel comparison (residual positive, zero
initial data, positive Green function: correct direction); convergence of the
singular integrals (`W_1H_*=O(t^{-1}\cdot t^2)` at zero and
`O((1-t)^{-5/6}\log)` at one). Nothing broke.

## Corrections to Claude's earlier handoff (accepted)

1. **"Theorem N implies at most three distinct zeros" was overstated.** The
   monotone-lobe count that gives three needs `H` to have exactly three
   simple interior zeros (the lobe region). Five distinct zeros force lobe
   membership; four do not. Outside the lobe region the inherited chain
   `Z(I)\le Z(H)+2` still allows four when `H` has two zeros. Correct
   statement: at most **three** distinct interior zeros on the strict lobe
   region, at most **four** globally. The Gavrilov–Iliev/Zhao conjecture
   (three) is therefore not settled by Theorem N; the outside-lobe four-zero
   question is the remaining proof obligation.
2. **The affine constant `c_0` in the handoff omitted `-K_0`.** Correct:
   `c_0=-306/1361360+\int_0^{y_1}W_a(K_3-K_0)`. No script used the wrong form.

## Consequences for the campaign

- The Q4 five-zero target of Attack 1 is excluded by theorem for every
  `kappa>1`. Combined with Lane B (the saddle loop cannot add cycles beyond
  small zeros of members of the same four-term space, and alien cycles do
  not occur at a single saddle loop), Q4 is closed as a counterexample route.
- STATUS.md and ATTACK_MATRIX.md are updated accordingly; the historical
  artifacts are untouched.
- What remains mathematically interesting in Q4 is the outside-lobe
  four-zero question (would complete a proof of the conjectured bound
  three). It is not a counterexample route.
