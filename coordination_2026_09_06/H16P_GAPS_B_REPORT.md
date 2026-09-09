# H16P GAPS-B (code + data) completeness report

Scope: every path listed in `/Users/scottg/Claude_all/H16P_GAPS_B_codedata.txt`, relative to
repo root `/Users/scottg/Claude_all/H16P` (read-only; nothing was modified).
Baseline: `/Users/scottg/Claude_all/H16P_SUMMARY.md` §2 (lines 108–630).
Every file below was opened. **No file is marked UNREAD.**

---

## PART 1 — COVERAGE MANIFEST

### Ignore files (4)
| Path | What it is |
|---|---|
| `.gitignore` | 2 lines: `__pycache__/`, `*.py[cod]` |
| `kkl/.gitignore` | 1 line: `__pycache__/` |
| `q4/.gitignore` | 2 lines: `__pycache__/`, `*.pyc` |
| `reversible_reseed/.gitignore` | 1 line: `__pycache__/` |

### `audit/` scripts (29) — all Claude-lane, all read in full
| Path | What it is |
|---|---|
| `audit/claude_check_boundaries.py` | Four Q4 boundary controls: cusp construction at `b=1.3` (asserts `q0(t*)≈q0'≈q0''≈0`, `q0'''<0`, `H_q0(t*)>0`), the unfolded `q_λ = q0 + λ(t−t*)` at `λ=1e-3` with `assert changes <= 1`; lobe-equality anchor-collision determinant at `δ=1e-2/1e-4/1e-6`; the P2 cubic identity; coefficient-map orientation. Prints `BOUNDARY CONTROLS PASSED`. |
| `audit/claude_check_strike1.py` | Hostile 40-dps re-derivation of Strike-1: Stieltjes density of `M`, companion identity `K=(1−t)(F+6tF')`, `R(0)=54/31`, `W3>0`, `W4<0`, `R` decreasing, β-strip map, closed `H(1)=18/(85085π)(9061A+6289B−2431η−7242)`, closed moments `J0..J2,K2,K3`, `J0(1)=18/(5π)`. |
| `audit/claude_corner_derivatives.py` | Large-κ corner linearization: `P_*(1;k)`, `Φ_*(1;k)` and ∂/∂(A,B,η), then the directional derivative along `dir=(−643/462, 1105/462, 1)`, for k = 8.5 … 10000. Prints the `v` needed for `P(1)=0`. |
| `audit/claude_corner_map.py` | Strike-4 viability control map: `P(τ1),P(τ2),P(τ3)` and `Φ(τ1)` along `γ(1−ε)` for ε∈{1e-2…1e-6} × k∈{5…200}; plus limit-point quantities at `(94/77,−17/77,1)`. |
| `audit/claude_laneB_loop_functionals.py` | Lane B: the two leading Dulac coefficients at the loop. Numerically differentiates `X(1)` w.r.t. (A,B,η) at κ=2,4,9; reports rank of the `(c0,c1)` map; then solves `c0=c1=0` and counts interior sign changes of `H` and `I`. |
| `audit/claude_laneC_focus_numeric.py` | 30-dps `mp.odefun` displacement `d(r)` at Shi's seed `x'=−y−10x²+5xy+y²`, `y'=x+x²−25xy`, r=0.02/0.04/0.08, checking `d(2r)/d(r)` against 32/128/512 for weak-focus order 2/3/4. |
| `audit/claude_laneC_locate_loop.py` | Bisects `l` at fixed `a∈{1.3,1.5,2.0}` on the stratum `m=5a, b=3l+5` using the signed separatrix offset from `claude_laneC_splitting`; reports `l*`, saddle trace, `η3`, and the type of `(0,1)`. |
| `audit/claude_laneC_loop01.py` | Homoclinic loop through the saddle `(0,1)`: crossing angles of the true stable branch (backward) vs the returning unstable branch on a circle `ρ=0.05`, with a winding gate `|wind|≥0.5` and a same-side gate `⟨eu,es⟩≥0.3`; brentq for `l*`. |
| `audit/claude_laneC_loop_probe.py` | At `(l,a)` = (−2,1), (−1,1), (−3,0.7): sympy-solves equilibria, classifies saddles, integrates all four separatrix branches for t∈[0,±60] and reports closest return to the saddle and winding about the origin. |
| `audit/claude_laneC_portrait.py` | Equilibria/eigenvalues of the Shi seed `(l,m,a,b)=(−10,5,1,−25)`, plus the infinity invariant directions from `Q2 − t P2`. |
| `audit/claude_laneC_saddle_region.py` | Symbolic resultant of `(P,Q)` on the stratum → the quadratic in `x`, its discriminant, a sign grid over `l∈[−20,5] × a∈(0,3]`, then equilibrium classification at four sample points with `η3` values. |
| `audit/claude_laneC_shi_focus.py` | **Exact symbolic Lyapunov quantities** for the Shi chart via the standard V-method to degree 8: prints `η1, η2, η3`, solves `η1=0 → m`, `η2=0 → b`, and evaluates `η3` on the stratum `b=3l+5` and at the seed `(l,a)=(−10,1)`. |
| `audit/claude_laneC_splitting.py` | First splitting attempt: unstable branch integrated to a `ρ=0.3` return event; offset = component along `vu`. Scans 9 `l` values at `a∈{0.7,1,1.5}`. |
| `audit/claude_laneC_splitting2.py` | Second rewrite: section Σ through `q = S + 0.15·vs` orthogonal to `vs`; splitting = signed coordinate along the normal. Scans 30 `l` values at `a∈{1.2,1.5,2,3}` and brentq's any sign change. |
| `audit/claude_laneC_splitting3.py` | Third rewrite: picks the *continuous* finite saddle with `x<0`; `σ` = sine of the angle between `vs` and the re-entry point on `|u−S|=ρ`, ρ=0.08; 28 `l` values at `a∈{1,1.2,1.5,2,2.5,3}`. Exports `saddle()` and `field()`. |
| `audit/claude_laneC_splitting4.py` | Fourth rewrite, the precise one: compares the crossing angle of the true stable branch (backward from `S+1e-8·vs`) with the returning unstable branch on `ρ=0.05`; brentq to `xtol=1e-12`; compares `l*` against the exact center curve `5a²l+6a² = 3l³+12l²+15l+6`. |
| `audit/claude_laneC_stratum_saddles.py` | Counts finite equilibria and finite saddles on the stratum over the Attack-2 box (`l∈{−12,−10,−8} × a∈{0.8,1,1.2}`), printing the factored resultant. |
| `audit/claude_laneC_unfolding.py` | First-order unfolding signs at Shi's seed with `m=5a+δ`, `b=3l+5−9δ+8ε`; prints exact `η1`, `η2` to first order in δ and ε, and `η3` at the seed. Loads `claude_laneC_shi_focus.py` by `exec(open(...).read().split(...)[0])`. |
| `audit/claude_q4_crossing_count.py` | Interior zero count of the first-order Q4 integral on both sides of the `X(1)=0` face, 800 sample points log-dense near `t=1`, for four anchor triples. |
| `audit/claude_q4_endpoint_c0.py` | `X(1−1e−6)` table for 10 lobe anchor triples × 7 lift parameters `a`, plus the corner point `(94/77,−17/77,1)`; counts negative `X(1)` values. |
| `audit/claude_q4_graphic_itinerary.py` | Bisects the annulus boundary on the +x axis at ρ=1, then integrates the two infinity-saddle separatrices in the `x=1/z, y=v/z` chart to confirm the graphic itinerary. |
| `audit/claude_q4_original_infinity.py` | Original-coordinate Q4 field at ρ=1: finite equilibria + eigenvalues, the three infinity directions `v = ρ, ρ±√(1+ρ²)` with angular/radial eigenvalues and saddle ratios, and annulus scans on ±x. |
| `audit/claude_q4_trace_boundary.py` | Traces the annulus-boundary level curve `Hcal = 4/(9κ)` from `(0.2272112871755, 0)` by flowing orthogonal to the gradient, both directions, arclength 400; reports asymptotic slopes, min/max radius, level residual, axis crossings. |
| `audit/claude_q4_triple_center_face.py` | Imposes `Y0=e1`, `η=e2`, `X(1)=0`, solves for `(A,B)` at `a∈{0.3,0.6,0.9,0.97}`, counts interior zeros of `I` on an 800-point log-dense grid; flags `n>=3` with ` <==`. |
| `audit/claude_q4_triple_center_face2.py` | Same face with the near-center hierarchy `|Y0| ≪ |η|/384`, `η<0`, at `a∈{0.3,0.6,0.9}` × `η∈{−1e−1,−1e−2,−1e−3}` × `Y0∈{−1e−5,−1e−7,−1e−9}`. |
| `audit/claude_q4_triple_center_point.py` | The unique triple-center point `(A,B,η)=(1,−17/12,0)`; scans `X(1;a)` on 20 `a` values, brentq's the roots `a*`, then unfolds with hierarchical alternating `(Y0,η,q0)` at 2×4×7 combinations; prints only rows with ≥3 interior zeros. |
| `audit/claude_q4_tworoot_scan.py` | Two-root region: solves for `(A−1,B)` from two prescribed primitive roots `r<s`, keeps only primitive count 2, then counts interior zeros of `I` over `r∈{0.15..0.7} × s × η∈{0.3..2.5} × a∈{0.2..0.99}`; prints a histogram. |
| `audit/claude_route4b_hemicycle.py` | Route 4b: splitting `D(l,a)` between the `x→+∞` saddle's unstable finite separatrix and the `x→−∞` saddle's stable one, measured at the `x=0` section; scans `a∈{0.5,1,1.5,2}` × `l∈{−30…−6}`. |
| `audit/claude_route4b_infinity.py` | Poincaré compactification of the order-3 stratum: real roots of `−t³−5a t²+(2l+5)t+a`, eigenvalues `λu = d/du[Q2−uP2]`, `λz = −P2`, type, and hyperbolicity ratio; Shi seed plus a stratum scan. |

### `audit/` logs (3)
| Path | What it is |
|---|---|
| `audit/large_kappa.log` | 14 rows: `P_λ(1)` claimed limits for λ=0.5 and λ=2, with `P(1−ε)`, its scaled value, `Φ(1−ε)` and `Y0`, for ε=1e-2…1e-7. |
| `audit/lobe_scan_phi.log` | 11 rows: `Φ(τ1)/|Y0|` for 10 anchor triples × k∈{1.2,2,4,8.5,20,100,1000}; final `BEST` line. |
| `audit/test_claude_hostile.log` | Six-test unittest transcript, all `ok`, `Ran 6 tests in 70.907s`, `OK`. |

### `controls/`
| Path | What it is |
|---|---|
| `controls/four_cycle_control.png` | Two-panel matplotlib figure titled *"Four-cycle regression (floating point; not proof)"*: left = three nested cycles about (0,0) at `x0 = 0.68321, 2.1837, 15.9628`; right = the remote cycle at `x0 = −3711.56`. |

### `fold_closure_2026_09_05/` (11)
| Path | What it is |
|---|---|
| `closure_budget.py` | 20-line charged-call harness. Docstring: `"""Frozen inherited accounting:756+3297+24=4077,19 calls remain."""`; hard cap `if n>=19: raise RuntimeError('shared4096 evaluation cap reached')`; records `source_sha256` and `cpp_sha256`. |
| `cusp_return.py` | Builds/runs `cusp_return.cpp` (g++ -O2 -std=c++17), mpmath dps 35, 33-digit decimal strings for `(r,c,K,B)`, `RLIMIT_CPU (10,10)`, 15 s subprocess timeout. |
| `cusp.log` | 16 charged calls, all `NUMERICAL_TWO_HALF_PASSAGES`. Purposes: 4× `fixed-focus augmented cusp residual`, 12× `finite-difference augmented cusp Jacobian`. |
| `cusp_attempt.json` | 4 `AUGMENTED_STEP` records of a damped 3×3 Newton on `(log r, c, B)` at `K=0.001953125`, with jacobian, newton_step, damping, base result and two perturbations each. Ends `"next_predictor": [3.911, 0.98975, 1.6951]`, `"component_excluded": false`. |
| `exact_checks.log` | Byte-identical to `theory_exact.json` (md5 `749623433745bac694f0b217cbc95740`). |
| `theory_exact.json` | Exact-rational Bernstein/Descartes output: 5×5 `N_left_bernstein_rows`, `N_K_derivative`, three `positive_root_coefficient_thresholds`, A- and −B-Bernstein rows on `c∈[1,11/10]` and `[11/10,8/5]`, `first_focal_numerator: "f_z*g_zz-f_zz*g_z = 2K"`, `analytic_density_necessary_order: "n = 4q = 16c/(2c+1)"`, `cycle_count_bound_proved: false`, `fold_component_excluded: false`. |
| `generalized_checks.log` | Byte-identical to `generalized_exact.json`. |
| `generalized_exact.json` | `EXACT_IDENTITIES_ONLY`: focus coefficient `(B*c*m−10*B−m−20)/(8*m^{3/2})`, center parity polynomial `B²c²+B²c−2Bc²−Bc−B+40c³−28c−10`, `center_projection_at_c1: 2*(B−1)²`, `K_at_B1_c1: −30`, third derivative of angular velocity, `cycle_count_bound: false`. |
| `outermost_check.py` | 12 lines. Pulls the last ACCEPTED event from `../fold_surface_2026_09_05/events_quad.json`, takes the pair-profile samples at `log_offset ∈ {−0.3, 0, +0.3}`, and re-runs each as a **complete return** through `angular_quad.py` at `tol=2e-26`. Writes `certified: false`. |
| `outermost.log` | The three resulting charged calls (evaluations 17, 18, 19), all `NUMERICAL_ONLY`. |
| `outermost_check.json` | The paired half-claim / full-return records. See NEW FACT B-1. |

### `fold_surface_2026_09_05/` (66 gap entries)
**Scripts / sources (24):**
| Path | What it is |
|---|---|
| `angular_ld.cpp` | Long-double independent **full-turn** angular shooter, 9-state, Gragg midpoint + polynomial extrapolation, `H<=0.08`, step guard 2e6, rejection guard 2e5, `H<1e-20` resolution guard. Header comment: `// No interval arithmetic. A failed angular chart is unresolved.` Chart guard `if(!(G<0))throw std::runtime_error("angular chart lost monotonicity")`. |
| `angular_quad.cpp` | Byte-for-byte the same algorithm in `__float128` (`cosq/sinq/expq/logq/expm1q/fabsq/powq`), reading `(r,c,K,tol)` as decimal strings via `strtoflt128`, emitting `%.36Qg`. |
| `angular_m_quad.cpp` | The `(c,m)`-chart binary128 sibling: the only differences are `al=-K; ac=0; ak=-1;` and emitting `m` / `K = m(11c−5)/5 − 42`, `L_m`, `L_zm`. |
| `half_m.cpp` | Long-double **two-half** matcher in the `(c,m)` chart: `integrate(r,tol,+1)` and `integrate(r,tol,−1)` over `[0,π]`, emitting `F=a[0]−b[0]`, `G=a[1]−b[1]`, `F_z=e^{a1}−e^{b1}`, `F_c`, `F_m`, `G_z`, `G_c`, `G_m`. |
| `angular_ld.py`, `angular_quad.py`, `angular_m_quad.py` | Three near-identical build/run wrappers (hash-named `/tmp` executable, mpmath dps 35/50/50, 33- or 45-digit input strings, `RLIMIT_CPU (10,10)`, 15 s timeout). Note `angular_quad.py`/`angular_m_quad.py` still carry the copy-pasted docstring *"Build/run the archived long-double source"* and the error string `'long-double process exit '`. |
| `angular_return.py` | Pure-Python/scipy DOP853 full-turn origin return with scalar analytic variations, 9 states, `max_step=.02`, plus an independent axis-flux cross-check `fluxA = r/R · Q0/Qf · exp(div)`; `first_derivative_discrepancy = A − fluxA`. Uses `allow_nan=False` so a NaN becomes `UNRESOLVED`. |
| `build_summary.py` | Saved-events-only summarizer: reads the nine `events_*.json`, recomputes the exact-rational cubic-discriminant equilibrium gate per event, counts `root_sign_brackets`/`stationary_sign_brackets`, writes `component_summary.json`, `EVENT_LEDGER.md` and `source_manifest.json`. Sets `global_component_complete: False`, `three_origin_candidate: False`, `five_cycle_candidate: False`. |
| `continue_angular.py` | Fixed-log-amplitude continuation of the full-turn fold `L=L_z=0` in `(c,K)`, acceptance `max|F| < 2e-10`, 13 pair-side offsets `[−12…+12]`, live trigger `if len(root_sign_brackets)>=3: write K1_CANDIDATE_ANGULAR.json`. |
| `continue_fold.py` | The double-precision predictor/corrector for `D=D_r=0` (Engine-2 `compact` evaluator). Acceptance `max|f|/max(|K|,1e−12) < 2e−7` **and** `|first_derivative_discrepancy| < 1e−6`. Also holds `topology()` (exact-rational cubic discriminant, `J = 305+634c−11c²−1000c³`, `K_H = −441J/(125(16−10c)(1+2c)²)`), `root_refine()` (safeguarded log-coordinate secant) and `paired_profile()` (12 offsets `[−5…+16]`). Trigger `>=3 → K1_TRIGGER.json`. |
| `continue_logm.py` | Continues the negative sheet in `log m` at `tol=2e-25`; acceptance `|F|<1e−22 and |G|<1e−18`; step `1.5`; `if used()>3100: break # reserve at least240 calls`. Trigger `→ K1_CANDIDATE_LOGM.json`. |
| `continue_m.py` | Continues through the removable `(c,K)` pole `c=5/11` in actual `(c,m)`; acceptance `|F|<1e−14 and |G|<1e−12 and exp(z)>1`; targets `.46, 5/11, .44, … .05`. Trigger `→ K1_CANDIDATE_M.json`. |
| `continue_negative.py` | Negative-K center-selected sheet; **relative** acceptance `|F|/|K| < 1e−13`, `|G|/|K| < 1e−11`; guard `c<=5/11 or c>2 → None`; 21 K targets `−1e−4 … −40`. Trigger `→ K1_CANDIDATE_NEGATIVE.json`. |
| `continue_positive_last.py` | Final 3 binary128 increments on the positive sheet at `tol=2e-26`; `if used()>3300: break`. |
| `continue_quad.py` | Binary128 two-sided continuation, `tol=2e-28`, acceptance `|F|<1e−26 and |G|<1e−20`; monkey-patches `base.refine` replacing `"1e-17"` by `"1e-28"` via `inspect.getsource` + `exec`. |
| `endpoint_controls.py` | Three binary128 endpoint controls: the neutral infinity connection at `c=1.6` for `r=1e6,1e10,1e14,1e17`; the `c=1` crossing; and the `K→0` center limit at `K=1/8192,1e−5,1e−7,1e−9` with `cstar = findroot(305+634c−11c²−1000c³, .969)`. |
| `final_controls.py` | Last controls: complete-return replay of the large-m pair brackets at `tol=2e-24`, plus two looser-tolerance retries of the `r=1e17` infinity control. |
| `full_return_replay.py` | Independent complete angular returns at the pair brackets of `events_quad.json` and `events_negative.json`. Docstring: *"This verifies two-cycle evidence only, not >=5."* |
| `graphic_coefficient.py` | Two Newton refinements of `K` at `r=1e17, c=8/5`, then `C = exp(5G/6)` and the conditional `(1.6−c)·log r` limit `−√159·G/12`. |
| `outermost_replay.py` | Same idea as `fold_closure/outermost_check.py` but writing `outermost_full_reproduction.json`. |
| `polish_replay_roots.py` | Newton/bisection polish of already-bracketed complete-return roots at fixed field, `tol=2e-26`, stop `|L|<2e−12` or `used()>=3338`; emits exact rational coefficient vectors. Writes `pair_replay_claims.json` with `independent_five_cycle_trigger: False`. |
| `refine_graphic.py` | `brentq(fun, 7, 7.5, xtol=2e-11)` on the infinity separatrix splitting via `graphic_shoot.py`, then three `v0` robustness re-runs. |
| `render_results.py` | Figure + report generator from `component_summary.json` (no ODEs); writes `fold_continuation.png/.svg` and splices the generated block into `../KKL_FOLD_SURFACE_STRIKE.md`. |
| `theory_negative_checks.py` | Exact algebra + saved-data audit of the negative sheet, **zero ODE calls**. Asserts `len(brackets)==2` for every accepted event, `q(left.F)*q(right.F)<0`, `m > m_H(c)`, `G_z<0`, `K<0`; and for the large-m replay `L_left·L_right<0`, `L_z` same sign, exactly one field. |

**Data / logs (42):**
| Path | What it is |
|---|---|
| `events_decreasing.json` | 4 events, all with 2 root brackets; `K = 1/1024 … 1/8192`, `c ≈ 0.96874 → 0.96864`, `r ≈ 6.9483 → 6.9476`. |
| `events_increasing.json` | 18 events, `K = 1/256 … 6`, `c = 0.96916 → 1.53776`, `r = 6.95 → 823.59`. 17 have 2 brackets; event 17 (`K=6`) has **1**. |
| `events_arclength.json` | 27 events, 17 accepted / 10 failed (`CORRECTOR_ITERATION_LIMIT` ×8, `UNRESOLVED_RETURN` ×2); `c = 1.5431 → 1.57398`, `r = 1212 → 214759`. All accepted have 2 brackets. |
| `events_angular_ld.json` | 2 events: 1 `ACCEPTED` (2 brackets, `r=204584.177…`, `c=1.57398`, `K=6.70042`), 1 `CORRECTOR_UNRESOLVED`. |
| `events_half.json` | 8 `ACCEPTED`, all 2 brackets; `r = 5.56e5 → 6.10e8`, `c = 1.57646 → 1.58588`, `K = 6.7475 → 6.9243`. |
| `events_quad.json` | 10 `ACCEPTED`, all 2 brackets; `r = 4.51e9 → 2.96e17`, `c → 1.59340580527813710990835865677884849`, `K → 7.06390700436779910773804298664181037`. |
| `events_negative.json` | 21 events (20 ACCEPTED + 1 `CORRECTOR_UNRESOLVED`), all 2 brackets; `K = −1e−4 → −30`, `c = 0.96861 → 0.47288`. |
| `events_m.json` | 7 events (6 ACCEPTED + 1 `CORRECTOR_UNRESOLVED`), all 2 brackets; `c = 0.46 → 0.35` **through** `c = 5/11`, `m = 352.7 → 188826.7`. |
| `events_logm.json` | 12 events (11 ACCEPTED + 1 `CORRECTOR_UNRESOLVED`), all 2 brackets; `m = 8.46e5 → 2.77e12`, `c → 0.3366860…` (approaching 1/3), `K → −7.17e11`. |
| `decreasing.log` | 121 charged calls, **all** `NUMERICAL_ONLY`. |
| `increasing.log` | 504 calls: 475 OK / 29 `UNRESOLVED`. |
| `arclength.log` | 661 calls: 591 OK / 70 `UNRESOLVED` — the worst failure rate of any lane (10.6 %). |
| `half.log` | 274 lines / 266 JSON: 259 OK, 7 `UNRESOLVED`; 8 `POINT … BRACKETS 2` lines. |
| `quad.log` | 210 lines / 204 JSON: 201 OK, 3 `UNRESOLVED`; 6 `POINT` lines. |
| `positive_last.log` | 139 lines / 135 JSON: 133 OK, 2 `UNRESOLVED`; 4 `POINT` lines. |
| `negative.log` | 697 lines / 677 JSON: 676 OK, 1 `UNRESOLVED`; 20 `POINT` lines. |
| `m.log` | 224 lines / 218 JSON: 217 OK, 1 `UNRESOLVED`; 6 `POINT` lines. |
| `logm.log` | 398 lines / 387 JSON: 382 OK, 5 `UNRESOLVED`; 11 `POINT` lines. |
| `angular_ld.log` | 17 JSON: 13 OK, 4 `UNRESOLVED` (`long-double fixed-amplitude fold continuation`); 1 `POINT` line. |
| `endpoints.log` | 35 calls: 34 `NUMERICAL_TWO_HALF_PASSAGES`, 1 `UNRESOLVED`. |
| `full_replay.log` | 12 calls, all `NUMERICAL_ONLY`. |
| `final_controls.log` | 8 calls, all successful (6 `NUMERICAL_ONLY` + 2 half-passage). |
| `graphic.log` | 18 calls, all `NUMERICAL_SEPARATRIX_PASSAGE`. |
| `graphic_coefficient.log` | 2 calls, both `NUMERICAL_TWO_HALF_PASSAGES`. |
| `polished.log` | 13 calls, all `NUMERICAL_ONLY`. |
| `exact_recheck.log` | Pretty-printed JSON, byte-identical to `theory_exact.json`. |
| `negative_recheck.log` | Pretty-printed stdout of `theory_negative_checks.py` (the `files` key elided) + three `accepted/failed` count lines. |
| `outcome_recheck.log` | Pretty-printed JSON, **byte-identical** to `theory_outcome_audit.json`. |
| `theory_exact.json` | `EXACT_IDENTITIES_PASS`, `orbit_evaluations: 0`; twelve named checks; `c1_no_cycle_K_threshold: "6292/1125"`; `algebraic_center_isolation`; `proved_normalized_bounds b∈(−3/10,−29/100)`, `d∈(101/100,51/50)`, `A` positive; closed forms `K2`, `K3`; `scope: "No global at-most-two theorem or complete fold-component exclusion."` |
| `theory_melnikov.json` | `NUMERICAL_QUADRATURE_ONLY`, 45 digits, 46 moment calls, 0 orbit evaluations; center `c=0.96862063355349…`, `m=37.13641480949…`; `limiting dc/dK = 0.13710961096153199…` reproduced at two FD steps. |
| `theory_negative_checks.json` | Full output of `theory_negative_checks.py`, including per-file sha256, accepted/failure counts and per-row remote-Hopf gaps. |
| `theory_outcome_audit.json` | Per-file saved-data audit: sha256, bracket histograms, min fold curvature, max residual components, max flux discrepancy, min `K_H − K`, exact remote-trace enclosures, stability-label counts. |
| `angular_transfer.json` | 3 double-precision full-turn returns transferring the fold seed to the horizontal section (evaluations 1300–1302). |
| `long_double_control.json` | Same point in `double` (eval 1302) and `long_double` (eval 1303) side by side. |
| `half_control.json` | Single long-double two-half record at `r=204584.177…`, `F=−2.55e−18`, `G=3.34e−15`, evaluation 1321. |
| `quad_control.json` | The same point in `long_double` (eval 1557) and `quad` (eval 1606), including `jacobian_cK_determinant = 0.042637228497748742888…`. |
| `derivative_controls.json` | The Cartesian seed record plus three analytic-vs-finite-difference checks (`r`, `c`, `K`), max relative error ≈ 2.3e−7. |
| `log_derivative_controls.json` | Two (cartesian, logarithmic) pairs at `K=1/512` and `K=6`. |
| `center_binary128.json` | `cstar = 0.96862063355349428616412539953798547325`; four rows `K = 1/8192, 1e−5, 1e−7, 1e−9`, all half-passages, with `secant_dc_dK` converging to `0.137109610975568…`. |
| `crossing_c1.json` | The `c=1` finite fold crossing, 8-step history: `r = 6.95315434526927557…`, `K = 0.224215783967328786…`, `F ≈ −6.4e−32`, `G ≈ 2.2e−31`. |
| `infinity_binary128.json` | Four finite-radius connection rows at `c=8/5`; `r=1e6/1e10/1e14` succeed with `K → 7.184994696406620401…`, `G → −0.2356749595172621…`. `r=1e17` is `UNRESOLVED`. |
| `infinity_tolerance_control.json` | Two retries at `r=1e17`, `tol=2e−24` and `2e−25`, both succeeding. |
| `graphic_coefficient.json` | Refined `K = 7.18499469640662101398737799635205628`, `G = −0.235674959517308217…`, `C ≈ 0.82168694694497326755611…`, conditional `Δlog r` limit `0.24764568047602043…`, `exact_or_interval: false`. |
| `graphic_initial.json` | Separatrix-shoot branches at `K=6.7` and `K=7`, `c=8/5`, with 9-term series coefficients and `endpoint_x` per branch. |
| `graphic_connection.json` | 9 refinement rows; `root_K = 7.184994696941435`, `certified: false`. |
| `full_return_reproduction.json` | Complete-return replays of the pair brackets from `events_quad.json` (`c=1.59162…, K=7.03092`) and `events_negative.json` (`c=0.472839…, K=−30`). |
| `large_m_full_reproduction.json` | The same for the `events_logm.json` pair, at `tol=2e−24`. |
| `pair_replay_claims.json` | `NUMERICAL_TWO_CYCLE_FIELDS_ONLY`, `independent_five_cycle_trigger: false`; three fields with exact rational coefficient vectors and two polished cycles each. |
| `runtime.json` | `python 3.12.13 … [Clang 22.1.3]`, `Linux-6.18.35-x86_64-with-glibc2.39`, numpy 2.3.5, scipy 1.17.0, sympy 1.14.0, mpmath 1.3.0, `g++ (Ubuntu 13.3.0-6ubuntu2~24.04) 13.3.0`. |
| `fold_continuation.png` / `.svg` | Three-panel figure: positive-K sheet in `(K,c)`; log10 section radius vs `c` (curved vs horizontal section); center-selected negative-K sheet in `(log10 m, c)` with the `1/3` organizer and the `5/11` chart pole. |

### `kkl/` (18)
| Path | What it is |
|---|---|
| `kkl/check_new_gate_review.py` | Bounded exact-only reviewer, `RLIMIT_CPU (10,10)`, `os.nice(10)`, **no ODEs**: verifies the original-coordinate reversing involution `M` and `num² = den²·I`, `linear·num + num·linear = 0`, the quadratic-field residual divisible by `J = 305+634c−11c²−1000c³`, `J(241/250)=25281/2500 > 0`, `J(39/40)=−11333/800 < 0`, the remote Hopf threshold `m_H` and `K_H = −441J/(125·d·e²)`, and the exact multiplier quartic at `c=1001/1000, m=196/5` (`K = 32039/6250`). |
| `kkl/follow_segment.py` | Segment-wise seeded continuation of the two known cycles; hard ceilings `assert 1<=steps<=256` and `used_steps()>=256`; predictor trust region `|dlog|>.65 → subdivide`; refuses `c >= 241/250` (*"infinity boundary: separate chart preflight required"*); asserts `origin R_r<1`, `remote R_r>1`, and windings `∓1`. |
| `kkl/pilot.py` | The charged-call harness: `if count>=4096: raise RuntimeError('full-strike evaluation budget exhausted')`; records `evaluator_sha256`; `json.dumps(..., allow_nan=False)`; `used_steps()` excludes `returns.jsonl` and `continuation_events.jsonl`. |
| `kkl/summarize.py` | Static summarizer: asserts contiguous evaluation numbering and `<=256` events; writes `strike_summary.json` and `SHA256SUMS`. Hard-codes `precursor_found: False`, `five_cycle_candidate_found: False`, `five_cycle_certificate: False` and the `not_excluded` list. |
| `kkl/data/pilot_summary.json` | `PILOT_COMPLETE_NUMERICAL_ONLY`, 64 evaluations, 7 accepted steps, `K=6/5`, `c∈[7/10, 33/40]`, `precursor_found: false`, `fold_found: false`. |
| `kkl/data/constant_K_path.jsonl` | 7 rows, all `ACCEPTED_NUMERICAL_POINT`; origin `r` 64.40→88.46, remote `r` −6081→−36339. |
| `kkl/data/reduce_K.jsonl` | 4 rows, all accepted; `c=33/40` fixed, `K` from `1157/1280` down to `1/64`. |
| `kkl/data/K_margin_toward_boundary.jsonl` | 5 rows: 4 accepted, 1 `UNRESOLVED` with `"error": "section range boundary"` at `c=149/160`. |
| `kkl/data/K_margin_to_c09.jsonl` | 6 rows: 5 accepted, 1 `UNRESOLVED` with `"error": "derivative cross-check failed"` at `c=9/10`. |
| `kkl/data/K_margin_c09_repaired.jsonl` | The single repaired row: `c=9/10`, `alpha=−67225/1568`, `K=1/64`, origin `r=45.061547076412836`, remote `r=−237627.19147658575`. |
| `kkl/data/last_inside_boundary.jsonl` | One row, `c=9301/10000`, `alpha=−8403125/209244`, `K=1/64`, origin `r=48.69483536581642`, remote `r=−1048286.5115219664`. |
| `kkl/data/continuation_events.jsonl` | The merged, evaluation-ordered ledger: 24 events, 22 accepted, 2 unresolved. |
| `kkl/data/exact_controls.json` | Exact field coefficients only: incumbent `c=7/10, alpha=−363889/5000, beta=3/2000`; positive-K start `c=7/10, alpha=−80, beta=0, K=6/5`. |
| `kkl/data/fold_direction_at_boundary.json` | `NUMERICAL_LOCAL_DIRECTIONS_NOT_GLOBAL_OBSTRUCTION`: `origin_mu=0.97634`, `mu_c=0.54959`, `mu_K=−0.071451`, remote log-radius sensitivities `59.7285` and `1.87493`. |
| `kkl/data/last_field_profile.json` | `SAMPLED_RETURNS_ONLY_NOT_ZERO_COUNT`: 10 samples on `r∈{1/64 … 1024}`. |
| `kkl/data/stationary_control.json` | `NUMERICAL_POSITIVE_MAXIMUM_NOT_A_FOLD`: three refinements to `R_r = 0.9999999999999931` at `r=28.174117155522246`, with `R_rr = −0.0011220650649…` (negative curvature → a maximum, not a fold). |
| `kkl/data/c_gt_1_remote_controls.json` | Four remote returns at `c=1001/1000, alpha=−196/5` for `r = −8, −512, −32768, −1048576`. |
| `kkl/data/remote_radius_boundary.json` | `NUMERICAL_PARAMETER_BOUNDARY_ONLY`, `K=1/64`, `section_r = −1048576`, 5 records bracketing `c ≈ 0.9301046…`. |
| `kkl/data/next_seed_geometry.json` | `EXACT_GEOMETRY_ONLY` for `c=1001/1000, alpha=−196/5, K=32039/6250`: exact rational cubic discriminant `−145984371027/1562500`, remote-x / trace / determinant / focus-discriminant rational intervals; infinity finite directions `−600±100√37`, both saddles, vertical eigenvalues `−1001/1000` and `−1/1000`, vertical type node; `"infinity_region": "requires separate itinerary audit"`. |

### `q4/data/` (10 gap entries)
| Path | What it is |
|---|---|
| `controls.json` | `NONRIGOROUS_NUMERICAL_CONTROLS`: a positive control (κ=4, three forced roots at s=2,3,3.7 with 60-digit residuals ≈1e−69/1e−70, `root_count: 3`), a negative control (`mu=[1,0,0,0]`, `root_count: 0`), a degenerate control at `s0=2.5`, and a three-way independent-evaluator agreement (`max_abs_gauss_difference 1.31e−14`, `max_abs_orbit_difference 5.54e−13`). |
| `smoke.json` | `NONRIGOROUS_SCREEN_ONLY`, seed 160926, 9 κ values, 2060 projective points, 3 survivors, root-count distribution `{0:0,1:0,2:0,3:3,4:0,5:0}`, `stopped_on_five_zero_lead: false`. |
| `second_green_shoot.json` | One `NUMERICAL_ONLY` shot with exact rational `A,B,η`; four P crossings; `S2_diagnostic: "fails: all four Z extrema are negative"`; `five_zero_candidate: false`. |
| `second_spine_diagnostic.json` | 6 records of fixed compact samples (65 points each), `sampled_sign_changes: 0` in all. |
| `third_initial_shoot.json` | 16 rows; 0 P crossings, 0 sampled original sign changes, `five_original_zero_candidate: false` everywhere. |
| `third_shape_shoot.json` | 3 rows on path `(r,1−(1−r)²,1−(1−r)³)`; 4 P crossings each, `S1_numerical_pattern: true`, but 0 original sign changes. |
| `third_tuned_shoot.json` | 5 rows on `(r,(1+r)/2,(3+r)/4)`; 4 P crossings, `S1_numerical_pattern: true`, 0 original sign changes. |
| `third_confluent_shoot.json` | 4 rows of primitive triple contacts; `"All four shots fail the first positive Z maximum numerically."` |
| `third_reverse_tangency.json` | 4 preselected reverse-tangency lines; `"No three-root H profile was detected."` |
| `third_independent_checks.json` | 2 cross-checks of the scalar Picard–Fuchs `I` against an independent original-area evaluation; `absolute_difference 1.08e−18` and similar. |

### `reversible_reseed/` (11)
| Path | What it is |
|---|---|
| `moment_search.py` | The main NUM scan: `profile(a,b,side)` builds each annulus's sampled moment curve (41 energies, Gauss order 160, logit η range [−9,9]); `arrangement(upper,lower)` sweeps every projection-order slope cell and counts upper/lower zeros of the shared perturbation; explicit `raise ValueError('logarithmic/degenerate shapes are outside this scan')` for `a∈{−1,−2,0}`. |
| `boundary_search.py` | The repaired-chart probes at the two excluded shapes `a=−2` (with `gamma` replacing `eps2`) and `a=0`, 41 energies, Gauss order 200. |
| `check_arrangement.py` | Software controls for `arrangement()`: two synthetic curves that must produce a five-candidate, and a parabola control that must not (`assert negative['maxima']['upper']<=2`, `assert negative['five_candidate'] is None`). |
| `data/moment_search.json` | 54 shape records, 0 errors. |
| `data/boundary_search.json` | 10 records (a=−2 and a=0 × five b). |
| `data/boundary_search.log` | The same 10 records as streamed one-line JSON, with `slope_cells` (1682–3322). |
| `data/control_search.json` | 1 record, the published Yu–Zeng shape `a=−671/450, b=7/15`. |
| `data/arrangement_controls.json` | The three control outputs; both synthetic cases fire a `five_candidate`, the parabola gives `null` and `maxima.upper = 2`. |
| `data/geometry.json` / `data/geometry.log` | Identical content: the exact chart algebra, incl. `full_chart_determinant: "-16*(a + 2)"`, `boundary_replacement_determinant: "-48*b"`, `boundary_invisible_old_direction {eps0:(2−b)/(12b), eps1:(1−b)/b, eps2:1}`, `infinity_horizontal_eigenvalues ["-a","-a - 2"]`, `yu_zeng_shape {a:−671/450, b:7/15, lambda:229/671}`, `bicycle_seed {a:−7/3,b:1}`. |
| `data/verify_control.log` | The four-cycle positive control: 6 high-precision moment rows (65 digits) then 24 original-field integration rows at `τ = 1e−4, 5e−5` × `rtol = 2e−11, 2e−13`; final `seconds 3.168144650000613`. |

### `review_2026_09_05/`
| Path | What it is |
|---|---|
| `check_nullcline_section.py` | Ten-CPU-second exact sympy replay of the curved section `(r, −r²/(1+r))`: `P|section = 0`, `Q|section = r·T/(1+r)²` with `T = (c−61/5)r³+(alpha−111/5−beta)r²+(2alpha−10−beta)r+alpha`, `Ṗ|section = (1+r)Q`, the first and second `σ_y` derivatives, the transversality determinant `det[F, tangent] = −Q_s`, and `det ∂(x,P)/∂(x,y) = 1+x`. No ODEs. |

### `staged_2026_09_05/` (23)
| Path | What it is |
|---|---|
| `cartesian_check.py` | Independent original-time Cartesian KKL section return: DOP853 `rtol=2e-13, atol=2e-15, max_step=.01`, curved nullcline section as the event, initial and final flux gates (`if field(0,state)[1]>=0: raise ValueError('bad initial flux')`, `if Q>=0 or (r>0)!=(R>0): raise ValueError('bad return branch/flux')`), 3 states (x, y, divergence). |
| `infinity_check.py` | Exact sympy infinity algebra, `orbit_evaluations: 0`: the vertical chart `(u,v)`, the finite-slope chart `(z,v)`, the vertical separatrix expansion `u=−v+a v²` with `a=−1/(1+c)`, the one-way `x=−1` crossing, `discriminant(p,z) = 40c − 964/25`, the mixed-stratum neutrality resultant `(c−1)J/25`, `J` sign at `241/250` and `39/40`, `J'(1/2)<0`, `J'' = −6000c−22`, the explicit neutral saddle root `z_n = (11c−5)/(5(1+c−2c²))`, and `c=8/5` as the only `c>1` root-sum solution. |
| `run_kkl.py` | The Stage-2 driver with four blocks. `correct()` accepts at `|log displacement| < 2e−8`, breaks if `|derivative| < 1e−6`. `validation()` asserts FD-vs-analytic `< 2e−5`. `terminal()` walks `c = 0.9301 … 0.968` at `K=1/64`. `profile()` runs 9 c × 6 K × **3 radii `(2., 20., 20000.)`**. `fold()` runs a 4-iteration damped augmented Newton seeded at `log(29.4), c=0.9301` for six K values. |
| `shi_trace.py` | The exact rational degree-8 Lyapunov polynomial builder + polar-angle DOP853 engine (`rtol=2e-11, atol=[2e-14, 2e-22, 2e-14, 2e-14]`, `max_step=.10`), with the cancellation-reduced diagnostic `H = ∫ (dV/dt)/V(initial) dt`. Guards `omega<=.03 or radius>3`, 10-CPU-second `SIGPROF` fuse, hard `CAP=160` return ledger. Families: `shi` = `(−10, 5+δ, 1, −25−9δ+8ε)` with `δ=−1/100, ε=−1/10^6`; `chen` = `(−3, 99/100, 2/9, −3)`. |
| `shi_mp_verify.py` | Independent 40-digit mpmath midpoint-extrapolation polar return (32 blocks × 7 levels), sharing the same 160-return cap; guard `if w<=mp.mpf('.03'): raise ValueError('polar chart lost')`. |
| `shi_run_continuation.py` | Drives two explicit trace paths (`shi` at λ from `−1e−14` to 0, `chen` from `−2e−5` to 0), brentq-locating roots in three (then two) fixed brackets, plus four intermediate-radius sign probes at an intermediate λ. |
| `controls.json` | 4 named controls (`origin`, `remote`, `cutoff_origin`, `cutoff_remote`), each with `old` (horizontal section) and `curved` records. |
| `derivative_validation.json` | Two analytic-vs-FD checks (errors `2.24e−8` and `−6.65e−9`) plus the re-corrected rounded remote seed at `r=−5656.065272434935`. |
| `validation.json` | `STATIC_NUMERICAL_EVIDENCE_CHECKS_PASS_NOT_INTERVAL_CERTIFICATE`: 206 prior + 400 new KKL + 150 new Shi = 756 of 4096; `kkl_historical_source_hashes_verified: true`; `shi_historical_source_hashes_available: false`; 394/400 KKL completed, 139/150 Shi completed; `exact_rational_pair_remote_focus: "UNSTABLE: 0<trace<0.006, determinant>88"`. |
| `theory_exact_checks.json` | `EXACT_CHECKS_PASS`, 0 orbit evaluations, domain `c∈[9/10,1] × K∈[0,6/5]`, 6+8+8 numerator coefficients and 50 left-interval `N` Bernstein coefficients with `min "0"`, plus a degree-9 `nonquadratic_logic_control` whose claim is `"Trace zero alone does not require a double cycle."` |
| `profiles.json` / `profile_run.log` | 162 `origin_displacement_profile` calls: 156 OK, 6 `UNRESOLVED` with `"no full return: A termination event occurred.; guard/time reached"`. |
| `terminal_path.json` | 12 rows `c = 0.9301 … 0.968`, `K=1/64`, both nests `ACCEPTED_NUMERICAL_ROOT` throughout; origin `r` 50.87 → 847.21, remote `r` −1.10e6 → −3.95e8. |
| `terminal_followup.json` / `terminal_followup.log` | 3 further rows `c = 0.9682, 0.9683, 0.9684`; the last has origin `CORRECTOR_NOT_CONVERGED`. Log = 22 charged calls, all `NUMERICAL_ONLY`. |
| `fold_attempts.json` | 6 rows (`K = 1/64, 1/512, 1/16, 1/4, 3/5, 119/100`), **every one `ITERATION_LIMIT` after 4 iterations**. |
| `fold_followup.json` / `fold_followup.log` | 2 further attempts (`K=1/512`, `K=1/64`), both `ITERATION_LIMIT`. Log = 24 calls, all `NUMERICAL_ONLY`. |
| `fold_refined_candidate.json` | 3 Newton refinements converging to `c = 0.9688884793906646`, `log_r = 1.9386104270243927` (`r = 6.949087993605231`), residual `(−1.80e−14, 1.89e−14)`. |
| `fold_pair_coefficients.json` | `NUMERICAL_TWO_CYCLE_PAIR_ONLY` with an exact rational field (see NEW FACT G-2). |
| `final_verification.json` | 8 independent-engine cross-checks (see NEW FACT G-3). |
| `shi_continuation.json` | The Shi/Chen three-cycle continuation results (see NEW FACT F-1). |

---

## PART 2 — NEW FACTS (not in `H16P_SUMMARY.md`)

### A. Scripts whose purpose or result is not captured by the summary

**A-1. `fold_closure_2026_09_05/closure_budget.py` records a *different* budget arithmetic than any figure in the summary.** Verbatim docstring:
`"""Frozen inherited accounting:756+3297+24=4077,19 calls remain."""`
and `PRIOR=4077`, with the hard stop `if n>=19: raise RuntimeError('shared4096 evaluation cap reached')`. So the `fold_closure` lane ran under a 19-call allowance, of which 16 went to the cusp Newton (`cusp.log`) and 3 to the outermost complete-return check (`outermost.log`) — the whole 4096-call campaign budget was exhausted at exactly **4096**.

**A-2. `fold_closure_2026_09_05/cusp_attempt.json` is the record of a *failed* cusp hunt, and the summary describes only that `cusp_test.py` "runs a 3x3 damped Newton".** The four recorded `AUGMENTED_STEP`s at `K = 0.001953125` **diverge in the third residual component**:
- step 0: `x=[1.9109341287192523, 0.9688884793906646, 2.2]`, `residual=[-1.14e-09, 4.26e-09, 0.011459427587586446]`, `damping=0.0652`
- step 1: `x=[2.4109…, 0.9763994442537136, 2.0515325580508392]`, `residual=[-0.006364, -0.022248, 0.030729]`, `damping=0.1107`
- step 2: `x=[2.9109…, 0.9822630517017438, 1.9157386242434187]`, `residual=[-0.012870, -0.046177, 0.047580]`, `damping=0.2850`
- step 3: `x=[3.4109…, 0.986602151575205, 1.7968214411725922]`, `residual=[-0.016056, -0.059223, 0.051657]`, `damping=0.4733`

`"next_predictor": [3.9109341287192523, 0.9897468071626071, 1.6951283547166496]`, `"component_excluded": false`. The `G_z` residual grows 4.5× over four steps; the run ends because the 19-call cap is reached, not because the cusp was located or excluded. **No cusp point exists anywhere in the repository.**

**A-3. `staged_2026_09_05/run_kkl.py::fold()` — the staged fold search never converged once.** `fold_attempts.json` has six rows (`K = 1/64, 1/512, 1/16, 1/4, 3/5, 119/100`) and **all six carry `"status": "ITERATION_LIMIT"` after exactly 4 iterations** (`for iteration in range(4)`). `fold_followup.json` retries `K=1/512` and `K=1/64` and both are `ITERATION_LIMIT` again. The augmented-fold acceptance is `if abs(f[0])<1e-8 and abs(f[1])<1e-7: status='CANDIDATE_REQUIRES_NONDEGENERACY_REPLAY'` — **that status string never appears in any saved artifact.** This is a second, independent instance of the "no augmented fold solver converged" gap the summary attributes only to `kkl/`.

**A-4. `kkl/check_new_gate_review.py` is an exact-only reviewer with zero ODE calls, and it verifies the *reversing involution* — a structure the summary never mentions.** It constructs `num = [[1000c²+769c+355, −10e(11c−5)], [−2100e, −(1000c²+769c+355)]]` with `den = 1000c²+1231c+145`, `M = num/den`, and asserts `num² = den²·I`, `det(num) = −den²`, `linear·num + num·linear = 0`, and that the quadratic-field residual under `z ↦ num·z` has every coefficient divisible by `J = 305+634c−11c²−1000c³`. Printed conclusion: `PASS exact original-coordinate reversing involution modulo J`.

**A-5. `review_2026_09_05/check_nullcline_section.py` proves the curved-section transversality identity exactly.** `Q|section = r·T/(1+r)²` where `T = (c−61/5)r³+(alpha−111/5−beta)r²+(2alpha−10−beta)r+alpha`; `Ṗ|section = (1+r)·Q|section`; `det[F, tangent] = −Q|section`; `det ∂(x,P)/∂(x,y) = 1+x`. This is the exact justification for Engine-2's curved section that the summary describes only prose-wise.

**A-6. `audit/claude_laneC_shi_focus.py` computes the exact Lyapunov quantities of the Shi chart symbolically and derives the third-order stratum from scratch** — the summary lists this only as "Shi focus quantities". It solves `η1=0 → m`, `η2|η1=0 = 0 → b`, substitutes the stratum `b=3l+5`, and prints `η3` at the seed `(l,a)=(−10,1)` and the curve where `η3` vanishes on the stratum. The stratum curve reappears explicitly in three downstream scripts as `5a²l+6a² = 3l³+12l²+15l+6` and as
`eta3 = -25*a*(2*a*a+l+2)*(5*a*a*l+6*a*a-3*l**3-12*l*l-15*l-6)/64`.

**A-7. `audit/claude_laneC_splitting*.py` is four *successive rewrites* of one function, and the summary counts them as generic scripts.** They differ materially in the section used: `splitting.py` uses a `ρ=0.3` re-entry event and the offset along `vu`; `splitting2.py` uses a line through `S + 0.15·vs`; `splitting3.py` uses `ρ=0.08` and the sine to `vs`, and adds a *continuity* selection (`px.real < -1e-6`) so the saddle varies continuously in `(l,a)`; `splitting4.py` compares the backward-integrated true stable branch against the returning unstable branch on `ρ=0.05` with `rtol=1e-12, atol=1e-15`. Only `splitting4.py` is a legitimate splitting function; the first three measure quantities that do not vanish at a loop.

**A-8. `reversible_reseed/check_arrangement.py` is a genuine positive+negative software control for the counting algorithm** — the summary does not list it. It builds synthetic wave/line/parabola curves and asserts the counter fires on the synthetics and returns `five_candidate is None` with `maxima['upper'] <= 2` on the parabola.

**A-9. `fold_surface_2026_09_05/angular_return.py` is a *third* independent full-turn evaluator (pure Python/scipy) that carries an axis-flux cross-check the summary attributes only to Engine 1.** Verbatim:
```python
Q0=-10*r*r+alpha*r;Qf=-10*R*R+alpha*R
fluxA=r/R*Q0/Qf*math.exp(div)
... first_derivative_discrepancy=A-fluxA
```
Same formula is compiled into `angular_ld.cpp`/`angular_quad.cpp`.

**A-10. `fold_surface_2026_09_05/build_summary.py` recomputes an exact-rational equilibrium gate per event, independent of the continuation.** For each accepted event it forms the exact cubic `A=c−61/5, B=alpha−111/5, C=2alpha−10, D=alpha`, its discriminant, and `remote_sign` from the value at `x_h = −21/(16−10c)` — then writes `exhaustive_root_coverage: False`, `global_component_complete: False`, `three_origin_candidate: False`, `five_cycle_candidate: False`.

**A-11. `continue_quad.py` and `continue_logm.py` patch their base module at run time by string substitution on source.** Verbatim:
```python
s=__import__('inspect').getsource(base.refine).replace("1e-17","1e-28");exec(s,base.__dict__)
```
and, in `continue_logm.py`/`continue_m.py`,
```python
s=__import__('inspect').getsource(base.profile).replace("K=m.mpf(a['K'])","K=m.mpf(a['m'])");exec(s,base.__dict__)
```
This is a real reproducibility hazard: the executed refinement tolerance and the executed profile chart are not what `continue_half.py`'s recorded source hash says they are.

**A-12. The binary128 wrappers still identify themselves as long-double.** `angular_quad.py`, `angular_m_quad.py` (and `fold_closure/cusp_return.py`) all carry the docstring `"""Build/run the archived long-double source; one process per charged call."""` and emit the error string `'long-double process exit '` even though they compile with `-lquadmath`. Any log triage keyed on that string will mis-attribute binary128 failures.

**A-13. `staged_2026_09_05/shi_trace.py` uses a cancellation-reduced return diagnostic that is not the displacement.** `H = ∫ (dV/dt)/V(initial) dt` where `V` is the exact rational degree-8 Lyapunov polynomial; the docstring states the caveat verbatim: *"H=0 iff the radius returns, provided V is monotonic on the section locally."* Root finding in `shi_run_continuation.py` brentqs on `H` (or, for the `shi` family, on the independent mpmath `raw_log_return`), never on `R−r`.

### B. Cross-engine disagreement and near-noise sign decisions

**B-1 [most important]. `fold_closure_2026_09_05/outermost_check.json` records a direct contradiction between the two-half matcher and the full-turn integrator at the outermost positive-sheet pair.** Same field (`c = 1.59340580527813710990835865677884849`, `K = 7.06390700436779910773804298664181037`, `alpha = -19.5825378385756324553907630103687574`), same radii, both binary128 at `tol=2e-26`:

| log offset | `r` | half-map `F` | half `F_z` | half `multiplier_at_match` | full-return `L` | full `multiplier` |
|---|---|---|---|---|---|---|
| −0.3 | `219194266668433180.697539869930428019` | `2.766e-25` | `-2.405e-24` | `0.998130066918646105` | `+2.13829603794767822e-04` | `0.999024725357680037` |
| 0 | `295881311432547652.547112249151329666` | `-4.509e-26` | `-2.524e-44` | `0.999999999999999999999971961` | `-5.12233717944015588e-05` | `0.999785400212164869` |
| +0.3 | `399397994234362784.599615340657854157` | `1.552e-25` | `+1.180e-24` | `1.001873436281833104` | `+2.43937901467019632e-04` | `1.002897959574552465` |

The half-map residual `F` says "closed orbit" (`|F| ~ 1e-25/1e-26`, i.e. at or below `continue_quad.py`'s own acceptance floor `abs(F[0])<m.mpf('1e-26')`) at **all three radii**, while the independent full return gives displacements of order `1e-4`–`1e-5`. If `F=0` really closed the orbit, `L` would be ~0. The half-map sensitivity has also collapsed (`F_z ~ 1e-24`, and `2.5e-44` at the centre), so `F` carries no discriminating power at `r ≈ 3e17`. The `+,−,+` pattern in `L` still supports exactly two cycles, but **the two bracket claims recorded for `events_quad.json` event 9 rest on `F` values that are numerically indistinguishable from zero**, and the multipliers disagree in the third decimal (`0.99813` vs `0.99902`; `1.00187` vs `1.00290`). The file itself is marked `"certified": false`. This is the concrete instance of the risk `theory_outcome_review.md` warns about; it is not recorded in the summary.

**B-2. `staged_2026_09_05/shi_continuation.json` infers signs from `H` values of order 1e−16 at an integration tolerance of `rtol=2e-11`.** The `shi` intermediate probes at `λ=−1e−16` are
`(r=4e-06, H=-4.272597137627752e-16)`, `(r=1.5e-05, H=+2.1980320092337277e-15)`, `(r=0.0015, H=-9.742569450079877e-11)`, `(r=0.025, H=-0.00021291187750138718)` — two sign changes, the first pair separated by `~2e-15`. The `origin_endpoint` at `r=1e-06, λ=0` reports `H=1.2566396061644519e-17`. Nothing in `shi_trace.py` or `shi_run_continuation.py` applies a noise floor of any kind.

**B-3. `fold_closure_2026_09_05/cusp.log` evaluation 1 has `F = -2.2266e-12` and `F_z = 9.3007e-13`.** The augmented Newton then divides by that `F_z`-sized Jacobian entry. `cusp_return.py` runs at `tol='2e-17'` in long double, so `F_z` is only ~4 orders above its own noise. This is the likely proximate cause of the divergence in A-2.

### C. Errors, crashes, and failure statistics not in the summary

**C-1. `fold_surface_2026_09_05/infinity_binary128.json` contains a hard process kill.** The `r=1e17` row is
`{"status": "UNRESOLVED", "error": "long-double process exit -9"}` — exit `-9` is `SIGKILL`, the `RLIMIT_CPU` fuse firing, not a chart failure. Two looser-tolerance retries in `infinity_tolerance_control.json` (`tol=2e-24` and `2e-25`, `cpu_seconds 1.79` and `2.26`) then succeed and agree to 30 digits, so the "hit the CPU fuse" phrase in the generated report corresponds to a real killed process, and the published `r=1e17` numbers come from the *retry*, not the original control.

**C-2. Charged-call failure rates per lane (from the logs; none of these appear in the summary):**
| Log | calls | UNRESOLVED | rate |
|---|---:|---:|---:|
| `arclength.log` | 661 | 70 | 10.6 % |
| `increasing.log` | 504 | 29 | 5.8 % |
| `angular_ld.log` | 17 | 4 | 23.5 % |
| `half.log` | 266 | 7 | 2.6 % |
| `logm.log` | 387 | 5 | 1.3 % |
| `quad.log` | 204 | 3 | 1.5 % |
| `positive_last.log` | 135 | 2 | 1.5 % |
| `negative.log` | 677 | 1 | 0.1 % |
| `m.log` | 218 | 1 | 0.5 % |
| `endpoints.log` | 35 | 1 | 2.9 % |
| `decreasing.log` | 121 | 0 | 0 % |
The four `angular_ld.log` failures are all in the `long-double fixed-amplitude fold continuation` purpose, which is why `events_angular_ld.json` accepts exactly one point and then records `CORRECTOR_UNRESOLVED`.

**C-3. `staged_2026_09_05/profile_run.log`: 6 of 162 calls fail with the identical message** `"no full return: A termination event occurred.; guard/time reached"` — all at `c=1.2, r=20000`, i.e. the profile grid's outer radius is outside the return domain for that stratum.

**C-4. Two distinct `kkl` continuation failure modes, both with verbatim error strings:**
- `kkl/data/K_margin_to_c09.jsonl` row 6: `"error": "derivative cross-check failed"` at `c=9/10`. This is the `4e-7` cancellation the summary describes; the *repaired* row lives in the separate file `kkl/data/K_margin_c09_repaired.jsonl` — so a naive read of `K_margin_to_c09.jsonl` alone sees only a failure.
- `kkl/data/K_margin_toward_boundary.jsonl` row 5: `"error": "section range boundary"` at `c=149/160`.

**C-5. `staged_2026_09_05/terminal_followup.json` shows exactly where the terminal path dies:** origin roots are accepted at `c=0.9682` (`r=2538.12`) and `c=0.9683` (`r=7036.91`), then `c=0.9684` gives `"status": "CORRECTOR_NOT_CONVERGED"` while the remote nest still converges (`r=-6.54e9`). The origin cycle radius grows by ~2.8× per `Δc=1e-4` in the final two steps.

### D. Coefficient vectors for 3-cycle fields (new — the summary lists no 3-cycle vector on `main`)

**D-1. `staged_2026_09_05/shi_continuation.json` — two exact rational three-cycle fields, both with three located origin roots.**
Chart: `x' = λx − y + l x² + m xy + y²`, `y' = x + a x² + b xy`.

*Shi conditioned*, `(l, m, a, b) = (-10, 499/100, 1, -3113751/125000)` (i.e. `δ = -1/100`, `ε = -1/10^6`), exact focus quantities
`η1 = 1/1000000`, `η2 = -738648442471027523/375000000000000000`, `η3 = 393801202719981275366738799578213/37500000000000000000000000000`:
- at `λ = -1e-14`: **three roots** at `r = 7.105453741451317e-05`, `7.087663979210806e-04`, `2.0166175831359916e-02`
- at `λ = 0` (trace zero): only **two** roots survive, `r = 7.123090765406943e-04` and `2.0166175706286565e-02` — the innermost cycle is lost.

*Chen–Wang visualization*, `(l, m, a, b) = (-3, 99/100, 2/9, -3)`, exact
`η1 = 1/400`, `η2 = -2097205519/7776000000`, `η3 = 19274211731911232971/604661760000000000`:
- at `λ = -2e-5`: **three roots** at `r = 0.0662982152603612`, `0.16576344873307564`, `0.33283851437480927`
- at `λ = 0`: two roots, `0.17802591969203013` and `0.3321212831537436`.

File-level caveats, verbatim: `"validated": false`, and
`"notes": ["Rational conditioned Shi point differs from Shi/Galias certified hierarchy.", "No interval integration or exhaustive cycle count."]`

**D-2. `reversible_reseed/data/verify_control.log` gives the four-cycle Melnikov sign pattern numerically, at 65 digits.** Upper annulus `normalized_M` at `h = 1, 1.5, 2, 3`:
`-0.00098736872…`, `+0.000051254063…`, `-0.000060717220…`, `+0.011394645644…` → `−,+,−,+` = **3 upper sign changes**.
Lower annulus at `h = 10, 16`: `+0.51343097735…`, `-0.43741848182…` → **1 lower sign change**. Total 4.
The original-field integration rows confirm the same pattern at `τ = 1e-4` and `5e-5` and at both `rtol = 2e-11` and `2e-13`; `D/τ` is stable to 8 digits between the two τ values (e.g. `0.0004045872659119709` vs `0.00040459536054804346`).

### E. Exact rational coefficient vectors for two-cycle fields (new — three of them, never quoted in the summary)

**E-1. `fold_surface_2026_09_05/pair_replay_claims.json`** — status `NUMERICAL_TWO_CYCLE_FIELDS_ONLY`, `independent_five_cycle_trigger: false`. Basis `['1','x','y','x^2','xy','y^2']`, `P = [0,0,1,1,1,0]`, section `y=0, x>0; clockwise full return`. Three fields, each with two polished complete-return roots:

| Source | `Q` (the `x` and `y²` entries; the rest are `0, 0, -10, 11/5`) | cycle 1 `r` | cycle 2 `r` | multipliers |
|---|---|---|---|---|
| `events_quad.json` | `x: -817181931567841384795366196973047495500/41692728942126575817276653714903670687`, `y²: 15916198802398157041075451013155546551/10000000000000000000000000000000000000` | `87981935401979.0568440822798367144836` | `113407241573834.922758790618347885218` | `0.999045918876561282691644902214499217`, `1.00105554343903999348616293643434359` |
| `events_negative.json` | `x: -6000000000000000000000000000000000000000/20123311963189102087151310255630448589`, `y²: 47283937451199009280650119114148222599/100000000000000000000000000000000000000` | `14.3542298817793356349251154370540839` | `18.4890091237174209870509476821525999` | `1.00563685295495388227003662498242614`, `0.99267545779659245107136032271250046` |
| `events_logm.json` | `x: -8645145685991851026666419025163957/3125000000000000000000`, `y²: 168343696162308681404693782077301251/500000000000000000000000000000000000` | `10524664759.6194216618298398678830992` | `13555129779.1358921811983386620604394` | `1.00108409791102483469108597236383389`, `0.998672869213268943582140987806142115` |

The `events_negative.json` field is a **small-amplitude** two-cycle field (`r ≈ 14.35` and `18.49`, periods `0.3278` and `0.3233`) — much more tractable than the 1e14/1e17 positive-sheet objects, and the file explicitly warns: *"Each row group is a separate two-cycle field; their cycles must not be added across fields."*

**E-2. `staged_2026_09_05/fold_pair_coefficients.json`** — a fourth exact rational two-cycle pair, on the curved section:
```
"field": {"c":"4844447396953323/5000000000000000",
          "alpha":"-1050048828125000000/28288921366486553", "beta":"0"},
"K":"1/512", "section":"(r,-r^2/(1+r))", "radii":[4, 6.949087993605231, 12],
"numerical_log_return_signs":["+","-","+"],
"remote_stability":"positive trace, exactly checked; not precursor"
```

**E-3. `fold_surface_2026_09_05/theory_negative_checks.json` gives the large-m two-cycle field as an exact rational pair**, cross-verified by the independent complete-return engine:
`c = 168343696162308681404693782077301251/500000000000000000000000000000000000`,
`m = 8645145685991851026666419025163957/3125000000000000000000`,
with `two_brackets_verified: true`, `left_L/right_L` of opposite sign in both brackets, `interval_certified: false`, `third_cycle_trigger: false`.

### F. Candidate / trigger files — the definitive negative result

**F-1. Six distinct live `>=3 root bracket` triggers exist in the code, and *not one* of the six output files exists.** Grepped and confirmed absent:
`K1_CANDIDATE_ANGULAR.json` (`continue_angular.py`), `K1_CANDIDATE_QUAD.json` (`continue_quad.py`, `continue_positive_last.py`), `K1_CANDIDATE_NEGATIVE.json` (`continue_negative.py`), `K1_CANDIDATE_M.json` (`continue_m.py`), `K1_CANDIDATE_LOGM.json` (`continue_logm.py`), `K1_TRIGGER.json` (`continue_fold.py`). The summary mentions only `K1_CANDIDATE_HALF.json`.

**F-2. Exhaustive bracket census across all nine `events_*.json` (109 events, 85 accepted):** bracket count is **2 for 84 of them and 1 for exactly one** (`events_increasing.json` event 17, `K=6`, `c=1.53776`, `r=823.586`). All 56 `POINT … BRACKETS n` lines across the seven continuation logs read `BRACKETS 2`. Maximum stationary-bracket count is likewise 2. `theory_outcome_audit.json` independently confirms via `saved_root_bracket_count_histogram` and a recomputed `adjacent_profile_sign_change_histogram` (identical in every file), with `bracket_or_label_problems: []` everywhere.

**F-3. `theory_negative_checks.py` hard-asserts the two-bracket count.** Verbatim: `assert len(brackets)==2`. The script would *crash* on a three-bracket event. Any future three-cycle discovery on the negative sheet will break this audit script rather than be reported by it.

**F-4. `q4/data`: the maximum interior sign-change count over every gap-list shoot file is 0.**
| File | rows | max `sampled_original_sign_changes` | max `P_crossings` |
|---|---:|---:|---:|
| `third_initial_shoot.json` | 16 | 0 | 0 |
| `third_shape_shoot.json` | 3 | 0 | 4 |
| `third_tuned_shoot.json` | 5 | 0 | 4 |
| `third_confluent_shoot.json` | 4 | 0 | 0 |
| `second_green_shoot.json` | 1 | 0 | 4 |
| `second_spine_diagnostic.json` | 6 | 0 | 0 |
| `third_reverse_tangency.json` | 4 | 0 | 0 |
`five_original_zero_candidate` is `false` in every row that has the field. `q4/data/smoke.json` gives `root_count_distribution_among_survivors: {"0":0,"1":0,"2":0,"3":3,"4":0,"5":0}` and `stopped_on_five_zero_lead: false`. The `S1_numerical_pattern: true` rows (shape/tuned shoots, 4 P crossings) are the closest approach, and their `P_extremum_gap` values are `5.61e-08` and `2.79e-08` — i.e. the fourth crossing is a near-tangency that the file's own warning calls *"floating diagnostics, not certificates."*

**F-5. `reversible_reseed`: 65 shapes scanned across three files, maximum total cycle count 4, no `five_candidate` anywhere.**
- `moment_search.json`: 54 shapes, 0 errors, best is `(upper 3, lower 1)` = 4 at the published Yu–Zeng shape `a=−671/450, b=7/15` (`m=3.80465698501488, c=0.09768891358101983`) and at `a=−1.25, −1.75, −1.8078` etc.
- `boundary_search.json`: 10 records at the two *excluded* shapes. The repaired `a=−2` chart also reaches `(3,1)` = 4 at `b=1/3` (`m=0.6671023222598771, c=1.9999999552160357`) and `b=2/3`. The `a=0` chart reaches only 2.
- `control_search.json`: 1 record, identical to the Yu–Zeng row.
The counter is validated by `arrangement_controls.json`, so the `maxima.total = 4` ceiling is a real result of the algorithm, not a detection failure. **This 65-shape negative sweep is entirely absent from the summary.**

**F-6. `q4/data/controls.json` documents a validated three-root positive control with 60-digit forced-root residuals** (`-1.70e-69`, `-7.81e-70`, `-2.24e-70` at `s = 2, 3, 3.7`, κ=4) and a three-way independent-evaluator agreement (`max_abs_gauss_difference 1.31e-14`, `max_abs_orbit_difference 5.54e-13`). The three-root control is the strongest positive control in the Q4 lane and the summary does not mention it.

### G. Independent-engine agreement results the summary does not record

**G-1. `fold_surface_2026_09_05/quad_control.json` — the long-double/binary128 agreement at `r = 609856836.87…`.** `F` differs in sign (`3.198e-18` vs `-2.073e-18`) but both are ~`1e-18`; `F_c` agrees to 16 digits (`-0.048316780804242417` vs `-0.0483167808042424208`); `G_c` agrees to only **8** digits (`-16.640545247931040881` vs `-16.640545253266223017`), and `G_z` to 11 (`0.013268726459541504` vs `0.0132687264592902187`). `jacobian_cK_determinant = 0.04263722849774874288863833383521200253` is recorded only on the long-double side.

**G-2. `fold_surface_2026_09_05/derivative_controls.json` — analytic-vs-finite-difference at the seed:**
`R_rr`: `3.180065322183623e-06` (FD) vs `3.179349684689914e-06` (analytic) — agreement to only **3 significant figures**, relative error `2.25e-4`. `R_rc` and `R_rK` agree to 6–7 figures. The `R_rr` disagreement is worth noting because `R_rr` is exactly the quantity that a fold detection depends on.

**G-3. `staged_2026_09_05/final_verification.json` — eight cross-checks, all `NUMERICAL_ONLY`, all passing:**
| Purpose | `r` | `D` |
|---|---|---|
| `fold_independent_cartesian` | 6.949087993605231 | `-1.4023981975697097e-10` |
| `fold_pair_lower_sign` | 4.0 | `+6.812280587631392e-06` |
| `fold_pair_middle_sign` | 6.949087993605231 | `-1.1352691329236109e-06` |
| `fold_pair_upper_sign` | 12.0 | `+6.21329898127243e-05` |
| `fold_absent_side_middle` | 6.949087993605231 | `+1.1349861575027376e-06` |
| `fold_pair_outer_profile` | 20000.0 | `+195.28978863163866` |
| `terminal_independent_origin` | — | `-8.2882788774441e-06` |
| `terminal_independent_remote` | — | `-23.180771827697754` |
The `fold_pair` triple gives the clean `+,−,+` pattern at `c = 0.9688894793906646` and the `fold_absent_side` gives `+` at the same middle radius at `c = 0.9688874793906646` — a textbook two-sided fold witness, produced by a *different* engine (`cartesian_check.py`) than the one that found it.

**G-4. `fold_surface_2026_09_05/theory_outcome_audit.json` reports maximum flux discrepancies per lane** (the summary quotes only `7.3e-11` and `6.9e-24`): `events_arclength.json` has `maximum_recorded_flux_discrepancy: 8.835876080470761e-06` — **eighteen orders of magnitude worse than the quad lane**, and its `max_abs_residual_components` are `[1.01e-08, 1.49e-08]`. The arclength lane is by far the least reliable accepted lane, and it is also the one with the 10.6 % failure rate (C-2).

**G-5. `staged_2026_09_05/validation.json` gives numbers not in the summary:** `new_kkl_evaluations: 400` of which `new_kkl_completed: 394` / `new_kkl_unresolved: 6`; `new_shi_evaluations: 150` of which `new_shi_completed: 139` / `new_shi_failed: 11` (7.3 % Shi failure rate); `kkl_recorded_success_cpu_seconds: 8.767690056`; and the exact-rational statement `"exact_rational_pair_remote_focus": "UNSTABLE: 0<trace<0.006, determinant>88"`.

### H. Correctness risks in the code

**H-1 [section detection]. `staged_2026_09_05/run_kkl.py::profile()` counts cycles on a three-point radial grid.** Verbatim: `for r in (2.,20.,20000.):`. That is a **10×** ratio and then a **1000×** ratio per interval. The `audit/fable_engine/REVIEW_engine.md` bug **B1** rejects a 1.50× ratio as too coarse. Across all 54 `(c,K)` cells in `profiles.json` the maximum sign-change count is **1**; two cycles closer than a factor 10 in radius — which is exactly the separation of the KKL origin pair at `r ≈ 4` and `12` (a factor 3) — are invisible to this grid by construction. The script's own comment admits scope but not resolution: *"Written finite design, not parameter-space coverage."*

**H-2 [section detection]. `kkl/data/last_field_profile.json` has the same defect at the campaign's terminal field.** Samples are `r = 1/64, 1/4, 1, 4, 16, 32, 64, 128, 512, 1024` — ratios of 16, 4, 4, 4, 2, 2, 2, 4, 2. Exactly one `D` sign change is found (between `r=32`, `D=+0.2344` and `r=64`, `D=-0.4921`). The file's own status is honest — `"SAMPLED_RETURNS_ONLY_NOT_ZERO_COUNT"` — but this is the only origin-side cycle census recorded at the last accepted KKL field.

**H-3 [tolerance floor]. Not one of the fold-surface continuation scripts applies a noise floor to a sign change.** `continue_fold.py::paired_profile` accepts a bracket on `if fa*fb<0:` alone; `continue_half.profile` (patched into `continue_m/logm/negative/quad`) on `if m.mpf(b['L'])*m.mpf(d['L'])<0:` alone; `continue_angular.py` likewise. The acceptance floors that *do* exist (`1e-26`/`1e-20` in `continue_quad.py`, `1e-22`/`1e-18` in `continue_logm.py`, `2e-10` in `continue_angular.py`, `1e-14`/`1e-12` in `continue_m.py`) apply only to the **corrector residual**, never to the bracket decision. B-1 above is what this costs at `r ≈ 3e17`.

**H-4 [tolerance floor, relative vs absolute]. `continue_negative.py` uses a *relative* acceptance and the others use *absolute* ones, on the same mathematical object.** Verbatim: `if abs(F[0])/abs(K)<m.mpf('1e-13') and abs(F[1])/abs(K)<m.mpf('1e-11')`. Since `events_negative.json` runs `K` from `-1e-4` to `-30`, the *effective* absolute tolerance varies by **five and a half orders of magnitude** along one continuation branch (`1e-17` at `K=-1e-4`, `3e-12` at `K=-30`). The same object in `continue_quad.py` is accepted at a fixed `1e-26`.

**H-5 [sign handling]. `half_m.cpp` computes `F_z` as a difference of exponentials, not `expm1`.** Verbatim: `emit("F_z",ea-eb);` with `ea=expl(a[1]), eb=expl(b[1])`. Near a fold both half-sensitivities are ~1, so `ea-eb` cancels catastrophically — and `F_z` is precisely the fold-detection quantity. The *full-turn* siblings do it correctly: `angular_ld.cpp` and `angular_quad.cpp` both emit `emit("L_z",expm1l(y[1]))` / `expm1q`. In `cusp.log` evaluation 1 this shows up directly: `F_z = 9.3007493931962156863413e-13` while `forward_log_sensitivity = -2.1915437098749747687858` and `backward_log_sensitivity = -2.1915437098832980253663` — 11 digits of cancellation.

**H-6 [sign handling]. The angular chart guard is `if(!(G<0))`, which conflates three distinct conditions.** In every `*.cpp` here it fires for `G > 0` (chart reversal), `G == 0` (tangency), and `G == NaN` — all reported as the single string `"angular chart lost monotonicity"`. NaN-safety is incidental, not designed; the log triage in C-2 cannot separate a genuine chart reversal from a numerical blow-up.

**H-7 [section detection]. `audit/claude_laneC_loop01.py` and `claude_laneC_splitting3.py` silently discard returns from the wrong side.** Verbatim in `loop01.py`: `if np.dot(eu, es) < 0.3: continue`; in `splitting3.py`: `if np.dot(e, vs) < 0: continue  # returned from the wrong side`. A homoclinic loop whose return grazes the stable direction at more than ~72° is dropped and reported as `None` (a `.` in the printed table), indistinguishable from "no return at all". Combined with the winding gate `if abs(wind) < 0.5: continue`, a genuine loop that winds less than half a turn is also invisible.

**H-8 [section detection]. `audit/claude_laneC_splitting3.py::saddle()` selects a saddle by `px.real < -1e-6` and takes `cands[0]`.** If two finite saddles with `x<0` coexist the choice is whatever order `S.solve` happens to return, which is not guaranteed stable in `(l,a)` — so the "continuous in `(l,a)`" claim in its own docstring is not enforced. `claude_laneC_splitting2.py` has a worse version: it iterates over all saddles and `return`s from inside the loop on the *first* one that produces any admissible crossing.

**H-9 [tolerance floor]. `run_kkl.py::correct()` accepts a root at `abs(d)<2e-8` while integrating at `rtol` defaulting to `2e-13`, and abandons at `abs(der)<1e-6`.** The `1e-6` derivative bail-out means the corrector *refuses to look* precisely in the neighbourhood where a fold lives — which is a mechanical explanation for A-3 (all six fold attempts hitting `ITERATION_LIMIT`).

**H-10 [reproducibility]. The three `theory_*` / `*_recheck` pairs are byte-identical duplicates, but one pair is not.** `theory_outcome_audit.json` ≡ `outcome_recheck.log` (md5 `97c167d099ae2e9adb183612de58bd06`) and `theory_exact.json` ≡ `exact_recheck.log` (md5 `749623433745bac694f0b217cbc95740`) — in `fold_closure`, `theory_exact.json` ≡ `exact_checks.log` and `generalized_exact.json` ≡ `generalized_checks.log` as well. But `theory_negative_checks.json` (a0fc20e…) **differs** from `negative_recheck.log` (37bf4aa…): the log is the abridged stdout with the `files` key stripped. Any integrity check that assumes the `.log`/`.json` pairs are interchangeable will produce a spurious mismatch on exactly one of the four.

**H-11 [tolerance floor]. `polish_replay_roots.py` stops polishing at `abs(m.mpf(a['L']))<m.mpf('2e-12')` even though the calls run at `tol='2e-26'`.** So the polished roots in `pair_replay_claims.json` (E-1) carry ~12 digits of displacement accuracy, not the 26 the tolerance suggests — and one of them, `events_negative.json` cycle 1, actually reports `L = -1.1245921497091985119905901310006487e-20`, far past the stopping rule, while `events_quad.json` cycle 0 stops at `L = +7.70554088054345293309001181247217091e-13`. The recorded root precisions are inhomogeneous by eight orders of magnitude across one file.

**H-12 [budget interaction]. Three continuation scripts abandon work based on a global counter rather than on convergence.** `continue_logm.py`: `if used()>3100:break # reserve at least240 calls for complete-return controls and handoff`; `continue_positive_last.py`: `if used()>3300:break`; `polish_replay_roots.py`: `if used()>=3338:break`. Every "the continuation stopped here" statement downstream of these is a budget artifact, not a mathematical boundary — and `closure_budget.py` (A-1) shows the campaign then ran to exactly 4096.

**H-13. `audit/claude_laneC_unfolding.py` loads its dependency by textual `exec` of a sibling file.** Verbatim:
`exec(open(__file__.replace('unfolding','shi_focus')).read().split("for i, e in enumerate(etas, 1):")[0])`
This silently breaks if that exact loop header ever changes in `claude_laneC_shi_focus.py`; there is no import, no hash, and no version guard.

**H-14. `audit/claude_laneC_locate_loop.py` computes the `(0,1)` Jacobian by hand with the derivation in a trailing comment, and the comment does not match the code.** Verbatim:
```python
J01 = np.array([[m*1, -1+2*1], [1+b*1, b*0]])  # dP/dx = 2l x + m y = m ; dP/dy = -1 + m x + 2y = 1 ; ...
```
`dP/dx = 2l·x + m·y` evaluated at `(0,1)` is `m` — correct — but the literal `m*1` hard-codes `y=1` into a slot whose comment says the `l` term vanishes; and `dQ/dy = b·x = 0` at `x=0` is written `b*0`. It is right at this one point only and will silently produce nonsense if reused. The same matrix is duplicated in `claude_laneC_loop01.py::sigma` as `Jn = np.array([[m, 1.0], [1.0+b, 0.0]])`, so the two files agree — but neither derives it symbolically the way `claude_laneC_splitting*.py` do.

---

## Bottom line for the five-cycle question

Nothing in the GAPS-B set contains a record of `>=5` cycles, `>=3` origin cycles, or any fired candidate trigger. The maxima observed anywhere in these files are:
- **4** — reversible Melnikov (3 upper + 1 lower), at the published Yu–Zeng shape and, newly, in the repaired `a=−2` chart at `b=1/3` and `b=2/3` (F-5, D-2);
- **3** — Shi conditioned and Chen–Wang origin nests at small negative trace, with exact rational coefficient vectors (D-1); and the Q4 positive control (F-6);
- **2** — every accepted fold event across all nine `events_*.json` (F-2), with four new exact rational two-cycle fields (E-1 … E-3).

The two things worth acting on: **B-1** (the outermost binary128 pair's bracket evidence is at the noise floor and its two engines disagree on the multiplier in the third decimal) and **H-5** (the two-half matcher's fold-detection derivative `F_z` is computed by an `exp − exp` cancellation where the full-turn engines correctly use `expm1`).
