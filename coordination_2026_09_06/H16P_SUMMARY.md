# H16P repository summary

Prepared 2026-09-05 from a read of the full working tree at
`/Users/scottg/Claude_all/H16P` **and of all seven unmerged remote branches**.
Subject: hunting a counterexample to `H(2)=4` — one real planar quadratic
vector field with at least five limit cycles.

---

## 0. READ THIS FIRST: the checked-out tree is not the frontier

`git log` on the local clone ends at `45f4ea9` (2026-09-06 01:44 JST,
"Record exact fold certificate obstruction and complete outer pair replay").
`main` has 32 commits and looks complete. **It is not.** Seven remote branches
carry work that never merged, including the most current results, the entire
production return-map engine, the overnight cloud campaign, and two further Q4
strikes. Nothing in the working tree references most of them.

| Remote branch | Tip date | What lives only there |
|---|---|---|
| `origin/astra/fastra-afternoon-2026-09-05` | 2026-09-06 02:57 JST | **The frontier tip.** Superset of `claude/conjecture-progress-report-ixsmgv`. Contains `audit/fable_engine/` (the C return-map engine, all F-lane sweeps, `CLUES.md` night watch, `REVIEW_engine.md`), `FABLE_NIGHT_REPORT_2026_09_05.md`, `FABLE_LANES_2026_09_05.md`, `FABLE_D2_ORDER_TWO_LOOP.md`, `FASTRA_AFTERNOON_REPORT_2026_09_05.md`, `astra_afternoon_2026_09_05/` (the binary128 `full_return128` interface) |
| `origin/claude/conjecture-progress-report-ixsmgv` | 2026-09-05 17:38 UTC | Subset of the above |
| `origin/codex/fastra-d1-fold-counts-2026-09-05` | 2026-09-06 02:24 JST | `fastra_d1_2026_09_05/` — the D1 fold-sheet continuation, the **exact rational 3+1 four-cycle field in the KKL chart**, the counter-failure audit, `FASTRA_D1_REPORT_2026_09_05.md` |
| `origin/astra/fastra-d1-2026-09-05` | 2026-09-05 11:23 MDT | Earlier, smaller D1 snapshot with `counter_check/` (`discrepancy_ledger.json`) |
| `origin/astra/q4-determinant-2026-09-05` | 2026-09-04 20:13 PDT | `q4/sixth/`, `q4/seventh/`, `Q4_SIXTH_BOUNDARY_REDUCTION.md`, `Q4_SEVENTH_LIMITING_FACE.md` — the determinant reduction and the `a=1` three-zero theorem |
| `origin/astra/resonant-joint-2026-09-05` | 2026-09-04 18:59 PDT | `resonant/`, `RESONANT_JOINT_2026_09_05.md` |
| `origin/audit/post-q4-frontier-2026-09-04` | 2026-09-04 17:17 MDT | **Orphan branch, no merge base.** `CANONICAL_STATE.md`, `FOUR_CYCLE_SEED_LEDGER.md`, `HISTORICAL_FIVE_CYCLE_CLAIMS.md`, `FIVE_CYCLE_MECHANISMS.md`, `H16_CERTIFICATION_PLAN.md`, `frontier_2026_09_04/` (Galias–Tucker replay at 900–1200 bits) |
| `origin/fable/compute-{evolve,f3-lam0,f5-shi,pert}` | 2026-09-05 | Worker branches: raw sweep `.jsonl` output only |

Origin is `https://github.com/git-df-scott/H16P.git`.

**Headline result across everything: no five-cycle field, no five-cycle
candidate, no interval certificate. Maximum cycle count ever produced by any
method in this repository is four.** Every document says so explicitly; the
discipline about this is unusually good and no file overclaims.

---

## 1. REPO MAP

### 1.1 Top-level files on `main` (all prose unless noted)

| File | What it is | Status |
|---|---|---|
| `README.md` | Reverse-chronological index; newest section first | Current index, but blind to all seven branches |
| `STATUS.md` | Binding status, evidence classes, acceptance test for a counterexample | Current and authoritative on `main` |
| `LITERATURE_AUDIT.md` | Primary-literature findings, four-cycle history, Shi's failed five-cycle episode, special-family bounds | Current |
| `SOURCES.md` | 29-entry source ledger with DOIs/links | Current |
| `RIGOROUS_CERTIFICATION.md` | The CAPD/interval certification protocol a counterexample must satisfy; certificate directory layout | Current, **never executed** |
| `SEARCH_SPACE.md` | Five essential parameters, Shi/Kuznetsov/center charts, why a grid is insane | Current |
| `BIFURCATION_MECHANISMS.md` | Mechanism ledger: what can and cannot emit a fifth cycle | Current |
| `FOUR_CYCLE_FRONTIER.md` | The four explicit four-cycle families + the reproduced Kuznetsov control | Current |
| `ATTACK_MATRIX.md` | The three authorized attacks with scores/budgets | **Largely superseded**; Attack 1 marked CLOSED, Attack 2 marked topologically defective |
| `ASTRA_HANDOFF.md` | Second-strike handoff | Superseded (banner) |
| `ZERO_TO_CYCLE.md` | The perturbative theorem taking five Abelian zeros to five cycles, and the realization gate | Current |
| `Q4_THEORY.md`, `Q4_PARAMETERIZATION.md`, `Q4_STRUCTURE.md`, `Q4_ZERO_GEOMETRY.md` | Q4 setup: normal form, coefficient transport, universal chart, Stieltjes/ECT structure, cusp exclusion | Current foundations |
| `Q4_CONTROLS.md`, `Q4_SEARCH.md`, `Q4_COST_MODEL.md` | Numerical controls, screening design, cost estimates | Historical; old pruning "has no mathematical force" |
| `Q4_CERTIFICATION.md`, `Q4_CERTIFICATE_PLAN.md` | The two-layer Q4 certification plan (interval Abelian zeros, then interval Poincaré) | Plan only, never executed |
| `Q4_LOBE_REGION.md`, `Q4_RECONSTRUCTION_GEOMETRY.md`, `Q4_THRESHOLD_PATH.md`, `Q4_GREEN_MAX_3.md` | Strikes 2–3: lobe cell, Green reconstruction, `5/11` threshold, first-Green-maximum obstruction | Current |
| `Q4_THEOREM_N.md` | **The load-bearing theorem.** Full proof of `Phi_a(tau_1)<0` | Current, independently audited |
| `Q4_TWO_ROOT_REDUCTION.md` | Strike 5: two-anchor reduction, mixed center sign, the residual determinant `K(r)` | Current |
| `ASTRA_FIRST_..._FIFTH_STRIKE.md`, `ASTRA_STRIKE5_CHECKPOINT.md` | Astra's five Q4 strike reports | Chronological record; all superseded by their successors |
| `CLAUDE_AUDIT_ASTRA_1_3.md`, `CLAUDE_AUDIT_ASTRA_4.md` | Adversarial audits; both verdicts SOUND | Current |
| `CLAUDE_LANES_B_C.md` | Lane B (Q4 saddle loop) and Lane C (Shi order-3 stratum) | Lane B's endpoint closure **later withdrawn** |
| `CLAUDE_ROUTES_4AB.md` | Route 4a (reversible `Q3R`) and 4b (order-3 focus + infinity graphic) | 4b closed numerically; 4a's premise **refuted by the council** |
| `CLAUDE_Q4_ENDPOINT_LANE.md` | Fable's Q4 two-saddle-infinity endpoint lane | Conclusion ruled **premature** by `FRONTIER_AUDIT.md` |
| `CLAUDE_THOUGHT_SESSION.md` | The codimension inventory of five-cycle mechanisms; the best strategic document in the repo | Current |
| `FASTRA_COUNCIL_2026_09_04.md` | The two-model council: route matrix, cross-examination, §9 Astra corrections | §9 governs |
| `FASTRA_H16_HANDOFF.md`, `FASTRA_H16_HANDOFF_5.md` | Handoffs to strikes 4 and 5 | Historical |
| `FASTRA_ZOOM_OUT_2026_09_05.md` | The campaign audit: what was actually tested, what was overstated, the next experiment | **Most useful single document for a newcomer** |
| `STRIKE5_PRECURSOR.md` | KKL pilot checkpoint, 206 calls, no precursor | Current for that pilot |
| `STAGED_RUN_2026_09_05.md` | The staged strike: finite KKL fold signal, 550 new calls | Current |
| `STAGED_K1_THEORY_2026_09_05.md` | K1 theory attempt; Bernstein amplitude restriction | Current |
| `STAGED_SHI_2026_09_05.md` | Two trace-to-zero paths on Shi and Chen–Wang | Current |
| `STAGED_INFINITY_2026_09_05.md` | Exact infinity audit of the KKL family | Current |
| `REVERSIBLE_RESEED_2026_09_05.md` | Reversible two-center reseed; the `a=-2` chart repair | Current |
| `KKL_FOLD_SURFACE_STRIKE.md` | The 3297-call fold continuation to `r~3e17` | Current |
| `KKL_FOLD_CLOSURE.md` | The Dulac-certificate impossibility theorem; budget exhausted at 4096 | Current, newest on `main` |

### 1.2 Directories on `main`

| Directory | Contents | Status |
|---|---|---|
| `audit/` | 40 Python checkers (Claude's independent re-implementations), 5 `.log` files, `test_claude_hostile.py` | Current; only 6 of 40 scripts have saved output |
| `controls/` | The Kuznetsov four-cycle numerical control: `reproduce_four_cycle.py`, `four_cycle_control.json/.png` | Current, non-rigorous |
| `council/` | `check_council_algebra.py` + three hostile note files | Current; contains the sharpest refutations |
| `q4/` | 22 scripts, 3 test suites, 15 note files, 2 JSON schemas, 17 data files | Current; the most mathematically serious code |
| `kkl/` | The KKL return-map pilot: 8 scripts, 6 notes, 17 data files, `SHA256SUMS` (34 entries, all verify) | Current |
| `staged_2026_09_05/` | The staged strike: augmented fold solver, Shi trace paths, exact Bernstein checks, 400+150 call ledgers, `MANIFEST.json` (37 hashes, all verify) | Current |
| `reversible_reseed/` | Reversible reseed: moment search, boundary search, 65-digit control, `MANIFEST.json` (19 hashes, all verify) | Current |
| `fold_surface_2026_09_05/` | 102 files. The long-double/binary128 two-half shooting engine, 9 continuation drivers, 3297-record ledger, 9 `events_*.json`, 4 theory notes | Current |
| `fold_closure_2026_09_05/` | 21 files. The cusp Newton attempt, the Dulac obstruction proof, 19-record ledger, budget exhausted | Newest on `main` |
| `review_2026_09_05/` | `COVERAGE_AUDIT.md`, `FRONTIER_AUDIT.md`, `KKL_NEXT_CONSTRUCTION.md`, `KKL_SECTION_REPAIR.md`, `EVIDENCE_MANIFEST.json`, transcript exporter | Current; the adverse self-review |

### 1.3 Directories that exist only on branches

| Directory | Branch | Contents |
|---|---|---|
| `audit/fable_engine/` | `astra/fastra-afternoon`, `claude/conjecture-progress-report` | **The production return-map engine** (`retmap.c`, `retmap_log.c`, `retmap.py`, `sweep_*.py`, `evolve*.py`, `f1x_*.py`), `REVIEW_engine.md`, `queue*.sh`, and `data/` with ~60 F-lane result files + `CLUES.md` |
| `astra_afternoon_2026_09_05/` | `astra/fastra-afternoon` | `full_return128.cpp/.py` — a standalone binary128 full-return interface with a documented API, 6 analytic controls, and the 24-K centerward sign map |
| `fastra_d1_2026_09_05/` | `codex/fastra-d1-fold-counts` (A) and `astra/fastra-d1` (B) | **Two disjoint D1 campaigns sharing a directory name.** A: `engine.py`, `matching_quad*.cpp`, `verified_precursor.json`, 64 folds. B: `run_d1.py`, `*_beta_quad.cpp`, `counter_check/`, `rational_3_plus_1.json`, 44 folds |
| `q4/sixth/`, `q4/seventh/` | `astra/q4-determinant` | The boundary-reduction determinant and the `a=1` three-zero theorem |
| `resonant/` | `astra/resonant-joint` | The resonant hemicycle compatibility obstruction |
| `frontier_2026_09_04/` | `audit/post-q4-frontier` (orphan) | Galias–Tucker/Shi replay at 900–1200 bits, Taylor orders 112/128 |

---

## 2. CODE INVENTORY

Everything is Python 3 + NumPy/SciPy/SymPy/mpmath, plus hand-written C and C++.
Declared pins: `mpmath==1.3.0`, `numpy==2.3.5`, `scipy==1.17.0`, `sympy>=1.13`,
`matplotlib==3.10.8`. Recorded runs actually used Python 3.12/3.14, sympy 1.14,
numpy 2.5.2, scipy 1.18.1 — **the pinned requirements do not match the versions
in the verification logs**. Nearly every script pins BLAS threads to 1, calls
`os.nice(10)`, and installs `resource.setrlimit(RLIMIT_CPU, (10,10))`.

### 2(a) The limit-cycle counting / return-map `D(r)` engines and their bug-fix history

There are **five distinct engines**, written by two different agents, and they
disagree with each other in ways that were only partly resolved.

#### Engine 1 — `kkl/return_map.py` (the original KKL pilot evaluator)

Cartesian DOP853, `rtol=1e-11`, `atol=tol*0.01`, `max_step=0.025`, **double
precision only**. Section: horizontal `y=0`, downward crossing. Integrates to
the *opposite* (upward) crossing, then to the *desired* (downward) crossing —
a full turn. `D = R - r`. Carries 12 states (2 orbit + 6 sensitivity +
1 divergence + 3 transverse determinant), 15 with second order.

Three independent first-derivative routes are computed and cross-checked:
the moving-event projection `V[0] - (P/Q)V[1]`, the transverse determinant
`-h_j(T)/Q_T`, and the section-speed × divergence exponential. Acceptance gate:
`abs(derivative_discrepancy) > 1e-7*max(1,|R_r|)` raises.

**Bug-fix history:** calls 1–131 used the fixed-time projection; at large remote
amplitude (`c=0.9`) subtraction of two large sensitivity components produced a
`4e-7` discrepancy and the continuation was marked unresolved. Calls 132+ switched
to integrating the transverse determinants directly,
`h_j' = div(F)h_j + det(F,F_j)`, `R_j = -h_j(T)/Q_T`. Agreement then improved to
`7.3e-11`. The old projection is retained as a cancellation diagnostic.
**Historical rows 1–131 carry no source hash and are not reproducible bitwise**
(`strike_summary.json` says so).

NaN handling: none explicit. `json.dumps(..., allow_nan=False)` *raises* on NaN,
and the outer `except` converts that to `{"status":"UNRESOLVED"}`. That is the
entire NaN policy — a NaN becomes an unresolved call, never a result. This is
actually sound.

Guards: coordinate cap `1e7`, horizon `t=10`, barrier `x=-1` for remote runs,
section-component check. A stopped integration is `UNRESOLVED`, never an escape
or a nonexistence proof. This discipline is uniform and correct.

**Fold detection: deliberately absent.** `kkl/continue_path.py` raises
`'near fold: switch to augmented fold equations'` when `|slope| < 1e-5`, **and
no augmented solver exists in `kkl/`.** `review_2026_09_05/COVERAGE_AUDIT.md`
identifies this as the single most consequential engine gap of the 206-call pilot.

#### Engine 2 — `staged_2026_09_05/compact_return.py` (log-polar, curved section)

Log-polar atlas `w = log hypot(x,y)`, angle, desingularized time
`dt/dtau = 1/(1+rho)`, on the **curved nullcline section** `(r, -r^2/(1+r))`.
Removes the old Cartesian radius cutoff via bounded coordinate `q = 1/(1+|r|)`.
Guards: log radius < 32, physical time < 10, desingularized time < 20000.

**`compact_return_v1.py` vs `compact_return.py` is a real bug fix**, and the diff
is exactly four changes:
- added a remote-branch guard (`remote point is not left of its equilibrium`);
- added an initial flux-orientation guard;
- added a final flux-orientation guard;
- renamed the output field and added
  `multiplier_interpretation='Floquet multiplier only at a fixed point'`.

Why it matters: on the curved section, `r < -1` does **not** imply the point is
left of the remote equilibrium or on the intended branch, and `log_flux` takes
`abs()` inside the log, so a **flux sign flip would be silently swallowed** and
produce a wrong-signed derivative with no complaint. The multiplier rename fixes
the exact conflation that `council/notes_kkl_hostile.md` §4 flagged: away from a
fixed point, `exp(∫div)` is a divergence factor, not a Floquet multiplier, and
using it as `P'` "corrupts fold location and derivative certification."
`compact_return_v1.py`'s hash matches ledger rows 1–392; rows 393–400 used the
guarded version.

**Fold detection here is real**: a damped 2x2 forward-difference Newton on
`(log r, c)` at fixed `K`, solving `log(R/r) = 0` and `d log(R/r)/d log r = 0`.
`refine_kkl.py` adds the essential rescaling `scale = (1-q)^2` that removes the
automatic quadratic vanishing of both equations near the focus — without it the
Newton system is degenerate at small radius.

#### Engine 3 — `fold_surface_2026_09_05/half_*.cpp` (the two-half matching engine)

This is the engine that carries the deepest continuation. Hand-rolled Gragg
modified-midpoint + polynomial (Bulirsch–Stoer) extrapolation on a 9-vector,
with **angle** as the independent variable. Precision ladder:
`double` (scipy DOP853, tol 2e-12) → `long double` (`half_ld.cpp`, tol 2e-17/2e-18)
→ **`__float128` / binary128** (`half_quad.cpp`, `half_m_quad.cpp`, with
`quadmath.h`, `strtoflt128`, `%.36Qg`, tol 2e-25 … 2e-28). Parameters are passed
as 33- or 45-digit decimal strings from mpmath specifically so binary64 is never
an input bottleneck.

The method: from `(r,0)` integrate clockwise to the negative ray, log radius
`A(z)`; from the same point integrate backward, log radius `B(z)`; solve
`F = A - B = 0`, `G = M_f - M_b = 0`. Then `F = 0` closes one full orbit and the
multiplier is `exp(G)`; `F_z = e^{M_b} expm1(G)` and, at a fold, `F_zz = e^{M_b} G_z`.
A fold is `F = G = 0`, solved by a 2x2 mpmath Newton with analytic derivatives.

Step control: `H *= min(1.5, max(0.5, 0.9*err^{-1/12}))`, `H <= 0.08`; rejection
guard at 2e5, step guard at 2e6, step-resolution guard at `H < 1e-20`.

**NaN handling:** no `isnan`/`isfinite` anywhere. Chart failure is caught by
`if(!(G<0)) throw std::runtime_error("angular chart lost monotonicity")` — the
`!(G<0)` form catches NaN incidentally. This is fragile but works.

**Noise floor:** hardcoded per lane, not adaptive. Acceptance `1e-16/1e-12` (long
double), `1e-26/1e-20` (binary128), `1e-22/1e-18` (log-m). The reviews are honest:
*"the sign of a residual smaller than actual accumulated error remains unresolved"*
(`theory_angular_review.md`), and `theory_outcome_review.md` warns that at very
large `r` a small `|F|` is weak evidence because sensitivity collapses too
(`F_z ~ -4.7e-22` at the final infinity control).

**Cycle counting: there is none.** `continue_half.profile()` samples ~14 fixed
log-radius offsets on the pair side and calls a sign change a bracket. There is a
live trigger `if len(root_sign_brackets) >= 3: write K1_CANDIDATE_HALF.json` —
**that file does not exist.** Every accepted event has 2 brackets (one has 1).

**Bug fixes caught by the independent reviewer** (`theory_angular_review.md`):
(1) the pair-side offset originally assumed positive curvature, `c_pair = c - sign(F_c)*delta`;
corrected to `sign(F_c * curvature)` before the negative-curvature sheet was explored.
(2) an unresolved-refinement handling defect. Both were fixed before publication.

#### Engine 4 — `audit/fable_engine/retmap.c` + `retmap_log.c` (branch only)

**This is the engine with the documented bug-fix history the question asks about.**
Dormand–Prince 5(4) in C with OpenMP, batched over parameter sets, plus a
compactified log-polar variant `retmap_log.c` in coordinates
`x = e^u cos θ`, `y = e^u sin θ` with rates multiplied by `e^u/(1+e^u)` so the
right-hand side stays bounded as `u → +∞`. Counts limit cycles as **sign changes
of `D = R - r`** on a radial grid from a ray leaving the focus, with a winding gate.

`audit/fable_engine/REVIEW_engine.md` is a 21 KB hostile correctness review that
found and quantified the following. Severity C = "wrong published counts".

| ID | Bug | Fix | Verified effect |
|---|---|---|---|
| **A1 [C]** | `atol = 1e-14*(1+r)` is an *absolute* floor, so any sign change with `|D| < ~1e-13` is integration noise. In `--lam0` mode (trace exactly 0) `D ~ V3 r^3` is tiny by construction | `atol = 1e-16*r` and `rtol <= 1e-12` | 3 of 4 re-checked `count=2` records in `F3_lam0_L4.jsonl` were actually `count=1`; one record flipped sign from `-2.06e-14` to `+1.27e-13` on tightening |
| **A2 [C]** | Step tolerance scaled by `\|x\|` (the *focus offset*), not the orbit radius. With the focus at `X0=1e4` a translated KKL field produced **10 spurious limit cycles** where the untranslated one produced 0, and tightening `rtol` made it *worse* (5 → 10) because the error is roundoff in the global coordinate | Re-expand the quadratic field about the focus and integrate in local coordinates | 0 spurious sign changes at every offset up to `1e8`; KKL and Yu-Zeng unchanged |
| **A3 [M]** | A dead third disjunct in the ray-crossing filter subsumes the two direction checks, so wrong-direction crossings are accepted. 22 of 10699 returns reported a different `R`, some negative | Delete the clause; reject `R <= 0` after bisection | — |
| **B1 [C]** | Radius grid `NR=40` over 6.5 decades is a **1.50x ratio per interval**. Two cycles closer than 1.5x cancel and count as zero. Yu-Zeng's three inner cycles are 1.3–2.3 intervals apart | `NR >= 160`, or coarse+refine at local extrema | Yu-Zeng under production settings reported **2 cycles, not 4**; at `NR=200, rtol=1e-12` it reports 4 |
| **B2 [H]** | Any single failed radius truncates the whole nest (`argmin(ok)`), and statuses 2/3 mean "gave up", not "boundary". Only 24% of 44240 radii returned successfully in a 512-set batch | Mark failures NaN, bracket only between adjacent successes | Yu-Zeng's origin nest truncated at grid index 15 of 40 |
| **B3 [H]** | `Tmax=2e3` is 50x below `count_nest`'s own default; `RMAX=3e3` is *smaller than the KKL remote cycle at r≈3711* the engine is validated against | Raise all four | — |
| **D1 [H]** | `c5 == c11 == 0` makes the resultant identically zero and **all equilibria are silently lost**. The Shi chart has `c11 ≡ 0` always | Handle the degenerate case explicitly | Any evolve candidate landing on that wall scores `total = 0` and is discarded |
| **E4 [L]** | `evolve.py`'s near-miss gate `rad > 3e-4*rad[-1]/1e8*1e8` reduces to `rad > 3*scale`, i.e. it **ignores every radius inside 3x the nearest equilibrium** — exactly where small nested cycles live | — | Mis-steered the whole evolutionary search |

The review also proposes a two-part noise-floor rule (roundoff floor
`10 u (|focus|+r) sqrt(N)` plus a differential truncation floor from recomputing
at `rtol/100`), with accept/reject/unresolved bands at `8z` and `2z`, validated
on six known points.

**What was actually implemented afterwards.** `retmap.c` now has
`atol = 1e-16*r + 1e-300` (A1 fixed). `retmap.py` has a module-level
`NOISE = 5e-12` and

```python
def count_signs(rad, D, noise=NOISE):
    for i in range(len(D)-1):
        if D[i]*D[i+1] < 0 and min(abs(D[i]), abs(D[i+1])) > noise*rad[i]:
```

`retmap_log.c` line 55 has the NaN fix:
`if (!(err <= 1.0)){ h *= isfinite(err) ? fmax(0.2, 0.9*pow(err,-0.2)) : 0.1; ... }`
— note `!(err <= 1.0)` is NaN-safe where the original `err > 1` was not.

The **"two-tolerance adaptive noise floor"** lives in
`audit/fable_engine/sweep_log.py::count_nest_log`:

```python
NOISE = 5e-12
def ret(us):
    """two-tolerance return: value at rtol 1e-12 and a per-point noise estimate from the 1e-11 run"""
    u1, S, st  = rm.returns_log(..., rtol=1e-12, **RT)
    u2, S2, st2 = rm.returns_log(..., rtol=1e-11, **RT)
    noise = np.where((st[0]==0)&(st2[0]==0), 10*np.abs(u1[0]-u2[0]) + NOISE, np.inf)
    return u1[0], st[0], noise
...
# adaptive refinement: interior local minima of |D| without a sign change (near-fold pairs)
idx = [i for i in range(len(D)-1)
       if D[i]*D[i+1] < 0 and min(abs(D[i]),abs(D[i+1])) > max(NZ[i], NZ[i+1])]
```

So: a per-point noise estimate from a *differential* tolerance comparison, plus
two rounds of refinement targeted at interior minima of `|D|` (the near-fold
pairs the coarse grid would miss), plus edge bisection in `u`. This is the
correct design and it is what the F15/F18/F19/F21 recounts used.

**Did the fixes change any verdict? No.**
`FABLE_LANES_2026_09_05.md`: *"All 24 four-cycle KKL c* fields recount to exactly
four. F19 rerun: the fifteen three-zero Q4 directions give at most three real
origin cycles at every amplitude, no alien, no four-zero direction. The KKL c*
and Q4 verdicts survive the noise-floor correction."*

The response to Astra's counter audit identifies three failure classes and
**one genuine bug**: *"Genuine bug found and fixed by the audit: NaN error
estimates were accepted by the step control."* The three classes are
(i) displacement below `1e-12` — precision limit, needs binary128;
(ii) radius above `e^36` — outside the domain cap by design;
(iii) coefficients above `1e6` — double-precision coefficient-magnitude limit.

#### Engine 5 — `fastra_d1_.../matching_quad*.cpp` (A) and `*_beta_quad.cpp` (B)

Binary128 two-state (A) or nine-state (B) angular matchers exposed to Python.
A's uses `y = sqrt(m) Y` and positive time rescaling to condition large `m`, and
is loaded as a shared library so profiles run in-process. A's controller fix that
matters: **the sensitivity tolerance is floored at `1e-24`**
(`(v==0 ? tol : std::max(tol, 1e-24Q))`), which is what unblocked the far-positive
endpoint after the first controller hit a step guard.

`astra_afternoon_2026_09_05/full_return128.cpp` is a cleaner sibling: a
documented one-call binary128 **full-turn** interface (not half-map matching)
that rejects float inputs, translates the focus with exact `Fraction`s before
rounding to binary128, rejects nonfinite intermediates *and error estimates*,
and returns 36-digit decimal strings with explicit statuses
(`OK_NUMERICAL`, `ANGULAR_CHART_UNRESOLVED`, `EVALUATION_LIMIT`,
`NONFINITE_OR_RANGE`). Six analytic controls pass, covering a 1e-17 linear-focus
displacement, log radii above 36, coefficients of magnitude 1e14, an anisotropic
center on a nonzero ray, a nonlinear rational reversible center, and a
deliberate chart failure. It is the best-engineered numerical artifact in the
repository — and its own README says *"They are numerical approximations, not
36 certified digits… No outward rounding or interval certification is performed."*

#### The most important engine finding of the whole campaign

On branch A, `FOUR_ORIGIN_TRIGGER.json`: the double-precision profile at
`positive_infinity_019_fold`, log-radius 20, reported **six sign changes,
stability `SUSUSU`** — which would have been a six-cycle origin nest. The
automatic four-root trigger halted the continuation.
`trigger_binary128_recheck.json` re-evaluated all twelve bracket endpoints in
binary128: **every one returns `F > 0`.** At log radius 19.58, double precision
gave `D ≈ -1.4047e-4` while the well-conditioned matching residual was
`+5.29655e-14`. Verdict: `"REJECTED_NUMERICAL_FALSE_POSITIVE"`.

That is the closest this campaign ever came to a five-cycle claim, and it was a
float64 artifact.

### 2(b) Interval arithmetic / Arb / rigorous certification code

**There is none. Not one validated ODE integration exists anywhere in this
repository or on any branch.** No Arb, no MPFI, no CAPD, no VNODE, no INTLAB,
no `mpmath.iv`, no Krawczyk or Kantorovich test. This is the single most
important fact about the codebase.

What does exist, and is genuinely rigorous within its scope:

1. **`kkl/geometry.py`** — a hand-rolled `Interval` class over
   `fractions.Fraction` with `+ - * /` and a division-by-zero-containing-interval
   guard. Used **exclusively** for exact rational finite-equilibrium gates:
   cubic discriminant, bisection isolation of the remote root, trace/determinant/
   focus-discriminant signs. Sound. **Never touches an orbit, a return map, or a
   cycle.**
2. **`q4/q4_lobe_certificate.py`** — pure `Fraction` + directed integer rounding,
   **zero third-party imports**. Certifies the signs of the primitive `H` at
   `t = 1/8, 3/8, 5/8, 7/8` for a frozen rational `(A,B,eta)`, with an
   **analytically proved geometric tail majorant** (N=256, largest tail
   `< 5.473e-20`), plus a full `1e-7` L-infinity box. This is arguably stronger
   than naive interval arithmetic. Output `q4/data/second_lobe_certificate.json`.
3. **`q4/q4_threshold_path.py::certify_frozen`** — the same method, N=1024, box
   radius `1e-8`, certifying a first primitive root beyond `23/32`.
   Output `q4/data/third_threshold_certificate.json`.
4. **`q4/sixth/check_exact.py`** (branch) — `Fraction` interval arithmetic on the
   four positive moment series through index 128, each tail bounded by its first
   omitted term over `1-t`, with interval Cramer determinants and the guard
   `assert other.lo*other.hi > 0, 'division interval meets zero'`. Produces the
   two confluent certificates
   `eta/(-192 Y0)|_{H_{7/10,7/10,1}} > 167/90` and `> 19/10` at `4/5`.
5. **`staged_2026_09_05/theory_exact_checks.py`** — exact rational **Bernstein
   coefficient enclosures** of the multiplier polynomial `N` on the closed cube
   `c in [9/10,1] x K in [0,6/5]`, 50 exact coefficients. A real, sound
   polynomial enclosure. Certifies the multiplier *density*, not any cycle.
6. **`fold_closure_2026_09_05/theory_exact.py`** — Bernstein positivity plus
   Descartes sign exclusions for `N > 0` on `-1 <= x <= 0` over `1 <= c <= 8/5,
   K > 0`.
7. **`fold_surface_2026_09_05/theory_exact.py`** — exact rational interval
   arithmetic isolating the algebraic center root `c_0`.
8. **`theory_outcome_audit.py::remote_enclosure`** — 64 exact `Fraction`
   bisections of the remote equilibrium and its trace, at finitely many stored
   decimal parameter points.

Everything else labelled "high precision" is **mpmath point arithmetic**
(dps 25–75) or **`__float128` shooting**, and every file says so. Representative
self-descriptions, verbatim from the code:
`q4/q4_integrals.py:6` — *"Neither routine is interval-rigorous."*;
`kkl/return_map.py` docstring — *"NOT a certificate."*;
`reversible_reseed/verify_control.py:3` — *"No interval claim: this is a
numerical positive control for four cycles."*;
`fold_surface_2026_09_05/theory_melnikov.py` — *"It uses mpmath point
arithmetic, not outward-rounded interval quadrature"*;
every `*_quad.cpp` line 2 — *"No interval arithmetic. A failed angular chart is
unresolved."*

Machine-readable status flags carried by essentially every JSON artifact:
`certified: false`, `interval_certified: false`, `exact_or_interval: false`,
`exhaustive_root_coverage: false`, `cycle_count_bound_proved: false`,
`five_cycle_candidate: false`, `five_cycle_certificate: false`.

The certification protocol that *would* be needed is written down three times —
`RIGOROUS_CERTIFICATION.md`, `Q4_CERTIFICATION.md` + `Q4_CERTIFICATE_PLAN.md`,
and (best) `H16_CERTIFICATION_PLAN.md` on the orphan branch. The last names
**CAPD::DynSys with an MPFR interval backend**, specifies seven algorithm steps,
gives precision guidance (128–256 bits for KKL, 512 near multiplier one,
1024–2048 for the Galias–Tucker stress control), lists a required per-orbit
manifest, and five controls including a deliberate **negative** control (a
quadratic center's period annulus must *fail* derivative isolation). It states
plainly: *"The current repository contains numerical replay scripts, not a
finished CAPD verifier."* Its four commands
(`certify-field`, `certify-return`, `certify-distinctness`, `verify-certificate`)
are a required contract, not installed executables.

There are two JSON schemas defining the promotion interface, and **nothing in
the repository validates against either**:
- `q4/zero_certificate.schema.json` — requires `format` const
  `"H16P-Q4-ZERO-v1"`, `kappa.rational`, `mu` (exactly 4 rational strings),
  `root_boxes` (**minItems 5**, each requiring `s`, `I_left`, `I_right`,
  `I_prime`, `interval_newton_image`, all as 2-element string interval
  endpoints), `endpoint_claim` in `{none, excluded_with_expansions}`,
  `evaluator{library, version, precision_bits, method}`, `artifact_hashes`.
- `q4/poincare_candidate.schema.json` — requires `format` const
  `"H16P-Q4-POINCARE-v1"`, `original_q4{rho,b,c}`,
  `perturbation{x_coefficients, y_coefficients, melnikov_realization_proof}`
  (note: a *proof artifact* is required, not just numbers), `epsilon.rational`,
  `sections` (minItems 1), `predicted_roots` (**minItems 5**), and optional
  `interval_fixed_points`, `disjoint_tubes`, `software`, `artifact_hashes`.

### 2(c) The certified four-cycle seeds and the "persistence tube"

**There is no certified four-cycle seed in this repository.** The only rigorous
four-cycle result in existence is the *published* Galias–Tucker (2022)
computer-assisted proof for the exact Songling point
`delta=-10^-13, epsilon=-10^-52, lambda=-10^-200`, and this repository
**reproduces its brackets numerically but does not reproduce its certificate**.
`FOUR_CYCLE_FRONTIER.md` notes that no maintained public one-command replay of
that proof was located.

The repository's own four-cycle objects, all **NUM only**:

| Seed | Where | Parameters | Result |
|---|---|---|---|
| **Kuznetsov (KKL) decimal control** | `controls/reproduce_four_cycle.py`, `controls/four_cycle_control.json` | `a=-10, b=2.2, c=0.7, alpha=-72.7778, beta=0.0015` — **binary64 decimal literals, not rationals**; `72.7778` is a truncation, not `655/9` | Section `y=0` downward: `0.6832102187597299` (S, mult `0.999226902940759`), `2.183699825305492` (U, `1.0024200551375815`), `15.962783983169892` (S, `0.9620208098988815`), `-3711.560806385011` (U, `11.46226771347891`). DOP853 `rtol=2e-11`, Brent on four **hard-coded brackets**. Closure error on cycle 4 is `1.07e-5`, residual `5.4e-9` — six orders worse than the others |
| **Exact rational KKL 3+1 (branch A)** | `<codex/fastra-d1>:fastra_d1_2026_09_05/verified_precursor.json` | `[0,0,1,1,1,0, 0, -2100097656250000000/56578038088396567, -1/10000000, -10, 11/5, 9688912553490597/10000000000000000]` i.e. `c=9688912553490597/10^16`, `K=1/512`, `beta=-10^-7` | Origin U/S/U at `0.0949678`, `5.26090`, `8.56091`; remote S at `10395387.5`. Multipliers `1.0000001021`, `0.9999958582`, `1.0000068840`, `0.9944348`. **Binary128** two-sided angular matching; remote equilibrium re-solved in binary128; tolerance tightened `2e-28 -> 2e-30` preserving digits |
| **Exact rational KKL 3+1 (branch B)** | `<astra/fastra-d1>:fastra_d1_2026_09_05/rational_3_plus_1.json` | `c=242288563571/250000000000`, `beta=-1/25600000`, `K=4660642062301237256681/1.25e24` — **a different field** | Origin U/S/U at `0.06955`, `7.88025`, `18.0424`; remote S at `19913343.5`. Cross-verified by two independent engines (log-polar and Cartesian). Status string: `"NUMERICAL_3_PLUS_1_NOT_A_COUNTEREXAMPLE"` |
| **Reversible rational four-cycle arc** | `reversible_reseed/data/verified_control.json` | `a=-7/4, b=1/3`, `eps0=-4τ/11, eps1=-31379τ/25000, eps2=-7517τ/5000` | Four Melnikov sign brackets (upper `-,+,-,+`; lower `+,-`) at 65 digits, persisting under original-field integration at `τ=1e-4, 5e-5` and two tolerances |
| **Visual Shi, Chen–Wang, Yu–Zeng, GT** | `<audit/post-q4-frontier>:FOUR_CYCLE_SEED_LEDGER.md` + `frontier_2026_09_04/data/exact_seed_coefficients.json` | Exact rationals given for all four | All four reproduced numerically. GT's three tiny cycles are only *bracketed*, at 900-bit MPFR Taylor orders 112 and 128, agreeing to 35 printed digits — no remainder enclosure |

**The exact rational 3+1 fields on branch A and branch B are the campaign's
single most valuable concrete asset**, and they are on unmerged branches. They
are what `review_2026_09_05/COVERAGE_AUDIT.md` §1 lists as never obtained
("Establish S/U/S plus remote U at one parameter — Never obtained. Central
construction gate remains unmet"). D1 closed that gap. `main` does not know.

**The "persistence tube" does not exist, in code or in prose.** I grepped the
entire tree and all seven branches. What exists is two adjacent and distinct
concepts that should not be conflated:
- **flow tube / isolating annulus** — a *phase-space* device for proving five
  cycles are pairwise **distinct**. Appears as a requirement in
  `RIGOROUS_CERTIFICATION.md:56`, `Q4_CERTIFICATION.md:78`,
  `H16_CERTIFICATION_PLAN.md:52,73`, `ATTACK_MATRIX.md:224`, and as the
  unconstrained optional field `disjoint_tubes` in
  `q4/poincare_candidate.schema.json`.
- **persistence box** — a *parameter-space* device.
  `H16_CERTIFICATION_PLAN.md:59` says "Taylor models help retain parameter
  dependence in persistence boxes", and the Hopf strike requires "a common
  parameter box proving persistence of J1–J4".

Neither has ever been computed. The nearest prose is one aspirational line in
`review_2026_09_05/COVERAGE_AUDIT.md:374` and one row of
`FASTRA_ZOOM_OUT_2026_09_05.md` §7 proposing to "interval-certify four relaxed
controls and a preserving neighborhood" for the GT relaxed field. Since all four
KKL multipliers are hyperbolic, a persistence neighbourhood **exists** by the
implicit function theorem — but nobody computed, enclosed, or even estimated one.

### 2(d) The five-parameter chart and the local-completeness determinant

The phrase "(A, m, beta, B, c)" matches the **generalized KKL chart**, which is
the family with the two normally-frozen quadratic coefficients released:

```
xdot = (1+x)y + x^2
ydot = A x^2 + B xy + c y^2 - m x + beta y
```

with `A = -10` and `B = 11/5` frozen in almost all work. `K = m(Bc-1) - A(B+2)`
is the first focus quantity, `l1 = K/(8 m^{3/2})`.
`fold_closure_2026_09_05/generalized_exact.py` is the one place `B` is released:
it derives `K = m(Bc-1) - 10(B+2)` and the exact reversible-center parity
polynomial `H(B,c) = B^2c^2 + B^2c - 2Bc^2 - Bc - B + 40c^3 - 28c - 10`
(which equals `2(B-1)^2` at `c=1`). `fold_closure_2026_09_05/cusp_return.cpp`
promotes `2.2` to a variable `Beta` and `cusp_test.py` runs a 3x3 damped Newton
on `(log r, c, B)` for `F = G = G_z = 0`.

**No local-completeness determinant is ever computed for this chart.** The
Jacobians that are computed in the KKL lane are all 2x2 or 3x3 *continuation*
Jacobians (`jacobian_cK_determinant`, `jacobian_zm_determinant`,
`jacobian_zc_determinant`, and `arc_continue.py`'s SVD `rank_singular_values`).

The genuine **local-completeness determinant computations in the repository are
two, and both are elsewhere**:

**(i) The reversible chart — `reversible_reseed/check_geometry.py`.** This is the
real one. Twelve columns spanning the 12-dimensional space of planar quadratic
fields: six infinitesimal affine changes, one time rescaling, `∂_a`, `∂_b`, and
the three perturbation directions of `(a, b, eps0, eps1, eps2)` with
`P = P0 + eps1 x + eps2 xy`, `Q = Q0 + eps0`:

```python
columns=[]
for g in [Matrix([1,0]),Matrix([0,1]),Matrix([x,0]),Matrix([y,0]),Matrix([0,x]),Matrix([0,y])]:
    columns.append(coeffs(J*g - g.jacobian(variables)*V))
columns.append(coeffs(V))
columns += [coeffs(V.diff(a)), coeffs(V.diff(b)), coeffs(Matrix([0,1])),
            coeffs(Matrix([x,0])), coeffs(Matrix([x*y,0]))]
chartdet = factor(Matrix.hstack(*columns).det())
```

**Result: `det = -16(a+2)`.** By the inverse function theorem the chart is a
locally complete unfolding modulo affine changes and time rescaling **iff
`a != -2`**. At `a = -2` the fix is to replace `eps2 xy` in `P` by `gamma x^2`
in `Q`, and then

```python
assert checks['boundary_replacement_determinant'] == '-48*b'
```

**`det = -48b != 0` for `0 < b < 2`.** The mechanism is the exact area-flux
identity `(a+2)J + b I_a + (1-b) I_{a-1} + ((b-2)/4) I_{a-2} = 0`: at `a=-2` the
three old moments become dependent, and the direction
`(eps0,eps1,eps2) = ((2-b)/(12b), (1-b)/b, 1)` has **identically zero first
Melnikov function in both annuli**. `gamma` restores the omitted moment
`J = ∬ x^2 y^{a-2}`. Correctly flagged as a *replacement*, not a fourth
independent Melnikov coefficient. This is the sharpest technical result in
`reversible_reseed/`.

**(ii) The Q4 coefficient-transport determinant — `q4/q4_structure_checks.py`.**

```
det ∂(alpha_1, alpha_2, beta_0, beta_1) / ∂(mu_1, mu_2, mu_3, mu_4)
  = 17843617792000 (kappa-1)^2 / (6561 kappa^4) > 0
```

Strictly positive for every `kappa > 1`, hence full rank 4. What it excludes: it
certifies that the transport into Zhao's `(alpha, beta)` coordinates loses no
projective direction, so filters proved there apply to the whole `RP^3` of
original coefficient directions. Exact/symbolic. A statement about coordinates,
not about zeros.

**Do not confuse either of these with the `q4/sixth`, `q4/seventh` determinants.**
Those are Wronskian/boundary determinants of the *reconstructed* function space
and are the subject of §3 below.

The other five-parameter charts, for completeness: the Shi chart
`(lambda, l, m, a, b)` (`SEARCH_SPACE.md`), the Kuznetsov chart
`(a, b, c, alpha, beta)`, the Marín–Villadelprat resonant chart
`(a, b, eps0, eps1, eps2)`, and the Q4 universal coefficient chart
`(A, B, eta)` with modulus `kappa` and universal moment function `M`
(`Q4_PARAMETERIZATION.md`: `eta = (kappa - beta_0)/d`,
`A = -(alpha_1 + 2 kappa alpha_2)`, `B = d alpha_2`, `1 < eta < 54/31`).
**No local-completeness determinant is computed for the Shi or Kuznetsov charts.**

### 2(e) Remaining runnable scripts, by directory

**`q4/` (22 scripts).** `q4_integrals.py` (the Abelian-integral primitives:
`basis_mp` at dps 60–75 via `mp.quad` with trigonometric cubic-root selection,
`basis_float` at Gauss order 96, and `basis_orbit_float`, an *independent*
evaluator integrating the Hamiltonian orbit with DOP853 and converting areas to
Green line integrals; plus `alpha_beta_from_mu` implementing Zhao's (20)→(21)
relabeling and `q4_coefficients(rho)`); `q4_search.py` (the bounded screen, with
the **corrected** strip and cubic filters, SVD triple-zero directions,
`robust_crossings` at relative floor `1e-10`, CPU fuse); `q4_controls.py`;
`q4_structure_checks.py`; `q4_lobe_certificate.py`; `q4_reconstruction.py` (the
Green/PF lift, DOP853 `rtol=3e-12`); `q4_lobe_anchors.py`; `q4_green_shoot.py`;
`q4_threshold_path.py`; `q4_green_max_3.py` (16 initial + 5 affine-tuned + 3
power-path shots, with a **live five-sign trigger that raises and freezes
`third_candidate_trigger.json` — never fired, that file does not exist**);
`q4_green_boundary_3.py`; `q4_reverse_tangency_third.py`;
`q4_third_independent.py`; `q4_green_endpoint_third.py`; `check_N_kernel.py`;
`q4_N_loop_checks.py`; `q4_fifth_green_checks.py`; `check_fifth_two_anchor.py`;
`check_fifth_loop_gate.py`; and three unittest suites `test_q4.py` (7),
`test_q4_second.py` (4), `test_q4_third.py` (4), all with SHA-256 equality
against frozen JSON. **Note:** these modules install a 10-CPU-second fuse at
import, so they must be run as separate processes; collecting them into one
pytest session kills the session and reports nothing.

**`audit/` (40 scripts).** Claude's independent re-implementations — the value
is that `audit/claude_green_tools.py` is a deliberate *second* implementation of
the Q4 engine, quadrature-only, which is what gives the strike checks their
force. Groups: six strike checkers (`claude_check_strike1/2`,
`claude_check_theoremN`, `claude_check_threshold`, `claude_check_boundaries`,
`claude_check_endpoint_identities`, `claude_check_reconstruction`,
`claude_check_large_kappa`); corner map/derivatives; Lane B loop functionals;
eleven Lane C scripts (Shi focus quantities, saddle regions, four successive
rewrites of the separatrix-splitting function, unfolding); ten Q4-endpoint
scripts; three route-4a/4b scripts; `claude_center_identify.py`;
`claude_lobe_scan_phi.py`; `claude_shots_phi.py`. **Only five have saved logs.**

`audit/test_claude_hostile.py` is **not** an adversarial numerics test — it is a
40-line unittest wrapper that subprocesses six checkers and asserts
`returncode == 0` plus a hard-coded success marker string. All six green,
`Ran 6 tests in 70.907s`. It does **not** cover `claude_check_theoremN.py`, the
one checker for the load-bearing result. The genuinely hostile content is *inside*
the checkers: `claude_check_strike2.py` brute-forces an attempt to violate
`Z(p_1) < 0` (fails to violate it), `claude_check_reconstruction.py`
cross-checks random `(A,B,eta)` at four `kappa` against independent 40-digit
quadrature, and `claude_check_threshold.py` re-derives the certificate from
scratch with its own rational series and tail bound.

**`kkl/` (8 scripts), `staged_2026_09_05/` (13), `reversible_reseed/` (5),
`fold_surface_2026_09_05/` (~30), `fold_closure_2026_09_05/` (6),
`review_2026_09_05/` (2), `council/` (1)** — described in §2(a) and §4.

**Integrity.** `kkl/SHA256SUMS` (34 entries) verifies; `staged_2026_09_05/MANIFEST.json`
(37 hashes) verifies; `reversible_reseed/data/MANIFEST.json` (19 hashes) verifies;
`fold_surface_2026_09_05/source_manifest.json` (41) and
`fold_closure_2026_09_05/source_manifest.json` (7) are present. **But every one
of these is a post-hoc snapshot.** `kkl/data/strike_summary.json` says it:
*"Original intermediate source hashes were not captured; this manifest records
the final reviewed implementation, not bitwise provenance for those earlier
rows."* `staged_2026_09_05/MANIFEST.json` is labelled
`"FINAL_FILE_HASHES_NOT_HISTORICAL_EXECUTION_HASHES"`, and
`validate_artifacts.py` records `shi_historical_source_hashes_available: false`.
KKL rows 1–131 and all 150 Shi rows have **no evaluator hash at all**.

---

## 3. MATHEMATICAL STATE

Graded bluntly: **E** = exact/symbolic or exact-rational-with-proved-tail;
**A** = analytic proof, not machine-checkable end to end; **N** = numerical only;
**W** = the argument as written does not carry the weight put on it.

### 3.1 The one result that closed a route

**Theorem N** — `Q4_THEOREM_N.md`, proved by Astra, independently audited in
`CLAUDE_AUDIT_ASTRA_4.md` with verdict `THEOREM N PROOF SOUND: YES`.

> For every finite `kappa > 1` and every coefficient point in the strict lobe
> region, `Phi_a(tau_1) < 0`, where `a = 1 - 1/kappa`,
> `Phi_a(t) = Y_0 + ∫_0^t W_a H` and `tau_1` is the first primitive root.

Grade **A**, resting on **E** ingredients. The proof is three global comparisons
plus one exact identity:
1. **(N0)** `∫_0^1 W_1 H_* dt = 3/1232 = -Y_{0,*}` at the corner
   `(A_*,B_*,eta_*) = (94/77, -17/77, 1)` with `H_* = (6t(1-t)^2 F/77)(6M-1) > 0`.
2. **(N2)/(N3)** The weighted moment curve `x = K1/K0`, `m = K2/K0` is strictly
   convex, with exact endpoints `x(1) = 6289/9061`, `m(1) = 11/41`, terminal
   slope `S(1-) = 1105/462`; hence for a two-root variation with `gamma > 0`,
   the center functional exceeds `(9/3080) gamma (25/231) > 0`.
3. **(N4)/(N5)** All three anchors strictly increase both `Y_0` and `H` before
   the first root (cardinal-interpolant sign `(-1)^{j-1}`); pushing all anchors
   to one along the fixed-ratio path gives `Y_0 < Y_{0,*}` and `0 < H < H_*`.
4. **(N6)/(N7)** `W_a < W_1` for all `a < 1`: the conjugation `v_a =
   z_a/(1-at)^{3/2}` gives `L_a v_1 = (1-a)(22 - 7d^{2/3})/(6 d^{7/3}) > 0`, and
   identical initial data plus a positive causal Green function gives `v_a < v_1`.
5. **(N8)** `Phi_a(y_1) < Y_{0,*} + ∫_0^{y_1} W_1 H_* = -∫_{y_1}^1 W_1 H_* < 0`,
   with an explicit compact certificate `Phi_a(y_1) < -395/3784704` on
   `a in [2593/21636, 1-1/64]`.

**Consequence, and this is the campaign's one genuine closure.** Combined with
the inherited necessary condition (N1) — five distinct original zeros force
strict lobe membership *and* `Phi_a(tau_1) > 0` — Theorem N gives: **every
nonzero Q4 Abelian integral, for every finite `kappa > 1`, has at most FOUR
distinct interior zeros; on the strict lobe region, at most THREE.** The Q4
five-zero construction route is closed by theorem.

**What it does NOT give, and the repository is scrupulous about this.** It does
*not* prove the Gavrilov–Iliev/Zhao conjecture (sharp bound three): outside the
lobe region the inherited chain `Z(I) <= Z(H)+2` still permits four when `H` has
two interior zeros, and Theorem N has no hypothesis there. Claude's earlier
handoff claimed the global three and **was corrected**; the affine constant `c_0`
in that handoff was also missing a `-K_0` term. Both corrections are recorded.

Audit findings worth carrying: (N0) was verified only to *quadrature precision*
(`2.435022e-3` vs `2.435065e-3`, residual attributed to the
`(1-t)^{-5/6} log` singularity) — the exact proof is the positive-series
argument, not the quadrature. `claude_check_endpoint_identities.py` admits in
comments that its closed-moment evaluator "loses ~20 digits to cancellation" and
that the second beta moment is checked only to `2e-3`. These are honest
precision compromises baked into "passing" tests.

### 3.2 Exact results that constrain but close nothing

| Claim | Where | Grade | Content |
|---|---|---|---|
| Corrected Zhao strip `(54-23κ)/31 < β0 < 1`; **no `κ < 85/23` cutoff** | `q4/notes_audit.md` | **E** | Zhao's printed Theorem 14/Corollary 18 has the **wrong sign**; direct differentiation gives `g'''(κ) = -25(23κ+31b-54)/(3888L³)`, matching his own Proposition 17. Every earlier search sampled an arbitrarily truncated domain, so **all old pruning is unsound** |
| Corrected cubic filter `P2(b) > 25(1-b)³/(432(κ-1)²)` | `q4/notes_audit.md` | **E** | The correct identity is `f''(s) = 2P/(s-b)³ + w''(s)`; both the factor 2 and the cubic denominator matter. The old linear bound excluded parameters the valid bound leaves undecided |
| **Reconstruction forcing sign: Zhao eq. (24) is wrong** in the repository's conventions | `q4/notes_reconstruction_second.md`, `audit/claude_check_reconstruction.py` | **E + N signature** | With the source's sign, the independent-area cross-check at `κ=2, s=1.5` differed by `4.54e-9`; with the corrected sign, `2.17e-19`. A real caught bug with a hard numerical fingerprint |
| Positive Stieltjes representation of `M`; auxiliary space is ECT; `Z(q) <= 3`, `Z(H) <= 3` | `Q4_STRUCTURE.md` | **E/A** | `M(z) = ∫_0^1 ρ(u)/(1-zu) du` with `ρ = 3/(2π²|F(1/u+i0)|²) > 0`; complex-analytic proof from the Euler integral. Gives `M' > 0`, `M'' > 0`, which everything downstream uses |
| Lobe region is a bounded contractible analytic cell, globally anchor-parametrized | `Q4_LOBE_REGION.md` | **E/A** | Explicit box `7/6 < A < 85/31`, `-1 < B < -49/744`, `1 < η < 54/31`; explicit 3x3 inverse; exhaustive boundary-event classification |
| Certified lobe box (3 primitive zeros, radius `1e-7`) — **and its own exclusion** | `q4/q4_lobe_certificate.py`, `Q4_LOBE_REGION.md` | **E** | Proved rigorously, then killed by the same strike: its first primitive root is `< 3/8 < 5/11` |
| `tau_1 > 5/11` and `kappa > 21636/19043` are necessary for five zeros | `Q4_RECONSTRUCTION_GEOMETRY.md` | **E/A** | The only nontrivial comparison reduces to the integer inequality `11·3^11 = 1948617 < 2000000 = 2^7·5^6` |
| The exact threshold path `gamma(r) = T(r,(1+r)/2,(3+r)/4)` crosses `5/11` transversally, with a second rigorous certificate (first root `> 23/32`) | `Q4_THRESHOLD_PATH.md` | **E** | Opened the gate; every Green shot through it then failed |
| Local auxiliary-cusp neighbourhoods excluded | `Q4_ZERO_GEOMETRY.md` | **E/A** | Scoped as local |
| **Theorem T** (two-anchor reduction) | `Q4_TWO_ROOT_REDUCTION.md` | **E/A** | Four distinct zeros force `H = (1-θ)B + θC`, `0 < θ < 1`, `H(1) > 0`, `Y_0 < 0`, and `P_B(r) > 0` **and** `K(r) > 0`. Excludes the whole `β1 = 0` sector, `Y_0 = 0`, all multiple-root cases, `H(1)=0`, the negative-orientation branch, and everything at or above the center-sign transition. Exact bounds `231/50312 < v < 3/616`, `Y_0(B+η_B V) > 27/12578` give `0 < λ_c < η_B` hence `η_C > 0` hence `P_C(0) = -η_C/192 < 0` |
| **Strike 6 exclusion wedge** (branch) | `<astra/q4-determinant>:q4/sixth/notes_exclusion_wedge.md` | **E** | Every four-zero candidate needs `r > 1-(7/22)^{3/2} = 0.82052...`, `η/(-192Y_0) > 19/10`, and `κ > κ_* ` with `2.899241080973277 < κ_* < 2.899241080974990`. Uses a supersolution `S_{q0} = A0(1-t)^{-4/3} - B0(1-t)^{-2/3}`, a positive Green comparison, and **two exact `Fraction` interval certificates with analytically bounded series tails** |
| **Strike 6 boundary reduction** (branch) | `<astra/q4-determinant>:q4/sixth/notes_boundary_reduction.md` | **A** | It suffices to prove nonpositivity at the two boundary choices `s -> r+` and `s -> 1-`; the second anchor can be eliminated. Proof: the cofactor primitive would otherwise land in the lobe closure where Theorem N's **strict tail margin** contradicts its forced `Phi(r)=0`. All three equality cases handled |
| **Strike 7: `a=1` three-zero theorem** (branch) | `<astra/q4-determinant>:Q4_SEVENTH_LIMITING_FACE.md` | **E** | Both boundary determinants strictly negative at `a=1`. Endpoint: `det = (27/1600)F³t²(M-1)(t-1)[t(M²+4M)-3M-2] < 0`. Confluent: `det = -(81F⁴t³/2620618000)Q(t,M)` with `Q > 0` proved by a Jensen strip `1/6 < M < 1/(6-5t)`, a substitution `t=z/(1+z)`, and the fact that **all 29 nonzero coefficients of the resulting `P(z,v)` are positive integers**. Hence the four-dimensional reconstructed space at `a=1` is an extended complete Chebyshev space, so at most three interior zeros. **Caveat the document states itself: `a=1` is the `kappa -> infinity` comparison limit, not a finite Q4 center** |
| **Strike 7: finite-lift subspace exclusion** (branch) | same | **E** | The 3-D subspace spanned by the reconstructions of `F`, `tF`, `t(1-t)G` is ECT at every `0<a<1` (Wronskians `F>0`, `F²>0`, `-(F³/6)(2M'+tM'')<0`), so at most two interior zeros there. Does not touch the fourth direction |
| KKL: `x=-1` is a one-way barrier (`P(-1,y) ≡ 1`) | `STAGED_INFINITY_2026_09_05.md` | **E** | Therefore the vertical infinity graphic is **impossible** for `c < 241/250` and every finite `m` — the north stable separatrix lies in `x<-1` (expansion `u = -v - v²/(1+c) + O(v³)`), the south unstable one in `x>-1`. The `2^20` radius cutoff was an artifact, not a mechanism |
| KKL: candidate eigenvalue-neutrality line `J(c)=0` at `c_* = 0.968620633553494...` for `241/250<c<1`, independent of `K` | same | **E** | Via `Res_z(p, (1-c)(1+z)-cp') = (c-1)J(c)/25`, `J = 305+634c-11c²-1000c³`. For `1<c<=3/2`, neutrality would need `c=8/5`, outside the box — excluded. **The connection, its splitting, and the transition constant `C` are all uncomputed** |
| KKL: at `K=J=0` the origin is an exactly reversible **double center**, not an order-three focus | `kkl/notes_local_unfolding.md` | **E** | Explicit involution `M`, `F(Mz) = -MF(z)`, `M²=I`, `det M = -1`. Local return divides as `D(r) = r³[K a + J r² b]` with `a,b>0`, so **at most one small nonzero cycle**, requiring `KJ<0`, never multiplier one. Kills the local missing-pair shortcut |
| KKL: remote Hopf is supercritical for `J<0`, giving a **stable** small cycle on the unstable-focus side | `kkl/notes_other_strata.md` | **E** | `l1_remote = 441J/(5000 e³ ρ² D^{3/2})`. **The ordinary remote Hopf cannot seed the required remote unstable cycle.** `K_H(c) = -441J/(125(16-10c)(1+2c)²)` |
| KKL Liénard: exact multiplier integral `log M = ∮ N ẋ²/(5(1+x)W²) dt`; every stable or multiplier-one origin cycle must reach `x > a(c,K)` | `kkl/notes_lienard.md`, `STAGED_K1_THEORY_2026_09_05.md` | **E** | The amplitude restriction is proved on the whole rational rectangle `9/10<=c<=1, 0<K<=6/5` by **exact rational Bernstein coefficients** (50 of them), not a sample grid |
| KKL: **no constant combination `N + λA` has one sign** on `(-1,∞)`; and neither `(c,m)` tangent is a rotated-field direction (`det(F,∂_α F) = xẋ`, `det(F,∂_c F) = y²ẋ`) | `kkl/notes_lienard.md` | **E** | Kills the naive monotone-ratio/Dulac exclusion and all monotone-displacement reasoning |
| Reversible two-center: exactly two finite equilibria, both centers, iff `0<b<2, a<=0` — **no finite saddle at all** | `REVERSIBLE_RESEED_2026_09_05.md` | **E** | Exactly excludes the old finite-saddle-loop mechanism for the entire two-center sector |
| Reversible chart completeness `det = -16(a+2)`, repaired at `a=-2` to `det = -48b` | `reversible_reseed/check_geometry.py` | **E** | See §2(d) |
| **Fold closure: the Dulac-certificate impossibility theorem** | `KKL_FOLD_CLOSURE.md`, `fold_closure_2026_09_05/theory_obstructions.md` | **E/A** | Two parts. (i) At a multiplier-one periodic orbit, `mu(P) = mu(0) = 1`, so `0 = ∫_0^P mu Φ dt` forces any one-sign continuous `Φ` to vanish on the whole projection; analytic ⇒ identically zero. (ii) The identically-zero escape fails: at a first-order weak focus a nonzero analytic solution of `X(Ψ) - 4κ(div X)Ψ = 0` has vanishing order `n = 16κ`, and a monic quartic in `T` has order at most 4; for `c > 1/2`, `16c/(2c+1) > 4`. **No certificate of this specified analytic monic-quartic scalar-residual form can have one sign across a true fold.** The one soft spot: the cubic-Hopf-normal-form step is invoked as standard rather than constructed |
| Fold closure: `N > 0` on `-1<=x<=0` and at most two positive roots, for all `1<=c<=8/5, K>0` | same | **E** | 5x5 Bernstein certificate on the `K=0` slice, `∂_K N >= 0`, Descartes with two Bernstein-certified threshold orderings. Consequence: the negative multiplier band is empty or a single interval. **Correctly disclaimed as not a cycle bound** — different cycles supply different positive weights in the multiplier integral, and that comparison is the missing proof obligation |
| **New complete section theorem** for KKL | `review_2026_09_05/KKL_SECTION_REPAIR.md` | **A** | Under stated gates, the rational nullcline `sigma(r) = (r, -r²/(1+r))` meets **every** origin cycle exactly once on `r>0` (at its x-maximum) and every remote cycle exactly once on `r < x_*`. Proof by the `(x,v=P)` diffeomorphism (Jacobian `1+x`, never zero), the sign of the cubic `T_β`, and a Jordan-curve ray intersection-number argument. Exact rational identities replayed. **Its numerical evaluator was never implemented and no control was transported** |
| Resonant joint compatibility obstruction (branch) | `<astra/resonant-joint>:RESONANT_JOINT_2026_09_05.md` | **A + E checks** | For a sequence carrying two compact cycles, at most **one** endpoint cycle can accompany them; if the compact distribution is (1,1), **zero**. Hence 2 compact + 3 hemicycle = 5 is impossible at the base `(-1,1)`, including for nonanalytic sequences. Symbolic assertions plus 70-digit quadrature cross-checks with `err < 1e-55` |
| Shi order-3 stratum: having the second focus **forces the extra equilibrium pair to be nonreal** | `<audit/post-q4-frontier>:frontier_2026_09_04/SHI_TOPOLOGY_AUDIT.md` | **E** | `N=(0,1)` is a focus iff `25a²+12(l+2)<0`, hence `l < -2-25a²/12`, hence `3a²-l(l+2) < -(7/6)a² - (625/144)a⁴ < 0`. So there is **no finite saddle** on the two-focus part of the stratum. This is the exact form of the "Attack 2 is topologically defective" finding |
| Q4 boundary graphic in original coordinates has **two** infinity saddles | `<audit/post-q4-frontier>:frontier_2026_09_04/Q4_GRAPHICS_AUDIT.md`, `audit/claude_q4_boundary_level.py` | **E** | The elliptic-chart saddle's two nodal branches map to `lim y/x = ρ-√κ` and `ρ+√κ`. Infinity directions: node at `v=ρ` (eigenvalues `-8,-4`), saddles at `v = ρ ± √(1+ρ²)` (`-4, 8`), 1:2 resonant, hyperbolicity ratios 2 and 1/2, product one. **This is what invalidated Lane B's endpoint closure** |
| Proposition A (order-two loop law), corrected form (branch) | `<astra/fastra-afternoon>:FASTRA_AFTERNOON_REPORT_2026_09_05.md` | **E** | In the zero-trace Shi chart with `eta_1 = 0`: if a non-origin equilibrium has zero divergence then the origin is a center — including the chart degeneracies `b=0` and `l=-1`, and the omitted case `H = a²-b(l+1) = 0`. **The literal converse in the original D2 note is FALSE**, with the exact counterexample `ẋ = -y+xy+y²`, `ẏ = x+x²+xy` (a reversible center with no neutral non-origin equilibrium) |

### 3.3 Claims graded W — where the repository over-reached, and caught itself

The repository's internal audits found these; I am recording them, not
discovering them.

1. **"Theorem N implies at most three distinct zeros globally"** — overstated in
   `FASTRA_H16_HANDOFF.md`, corrected in `CLAUDE_AUDIT_ASTRA_4.md`. The correct
   statement is three on the strict lobe region, four globally.
2. **Lane B's blanket Q4 endpoint closure** (`CLAUDE_LANES_B_C.md`) — **withdrawn**
   by `FASTRA_COUNCIL_2026_09_04.md` §0.1 and `CANONICAL_STATE.md`. The
   one-saddle alien-exclusion argument does not transport through the singular
   double cover; the rank-two check of `(c_0,c_1)` at three `kappa` values is a
   sample, not a versality theorem. `STATUS.md` on `main` still carries the old
   "Attack 1 is therefore closed" language for the *interior* route (correct) but
   the endpoint status is **UNKNOWN**.
3. **"Only two alien cycles remain" (Fable's endpoint lane)** — ruled
   **premature** by `review_2026_09_05/FRONTIER_AUDIT.md` §3, in five itemized
   findings: sampling `a` on `[0.05,0.995]` does not prove `X(1;a)>0` on the full
   admissible interval; 546 two-root samples do not cover the region or its
   discriminant strata; Theorem N bounds *distinct interior* zeros and does not
   give "<=4 after adding a boundary zero"; a boundary zero is not itself an
   isolated periodic orbit; and "first-order count + alien count = 5" is not a
   counting law until the arc and asymptotics are specified.
4. **Route 4b "CLOSED NUMERICALLY on the stratum"** (`CLAUDE_ROUTES_4AB.md`) —
   rests on a stated **conjecture** ("no separatrix cycle surrounds a third-order
   weak focus of a quadratic system") plus sampled `(l,a)` values. It was later
   *upgraded* to a genuine literature theorem on the afternoon branch: 
   Llibre–Schlomiuk 2004 **Theorem 16(III)** shows the only graphic portraits in
   QW3 are W13, W15, W18, each with one graphic surrounding the **strong** focus,
   all vertices at infinity. So statement (C) is verified from the classification
   — but `FASTRA_AFTERNOON_REPORT_2026_09_05.md` also warns that the same paper's
   Theorem 12 explicitly qualifies its connection set `G3` as *numerically
   observed*, so the full bifurcation diagram must not be cited as an exact
   component enumeration.
5. **Route 4a's premise is false.** `council/notes_reversible.md` proves exactly
   that the advertised "external foci" of the reversible seed are **saddles**:
   the off-axis equilibria `E_± = (-1/k, ±sqrt((k-p)/(k²q)))` have a triangular
   Jacobian with `det = -2(k-p)/k < 0`. In the original Shi chart the point
   advertised as the second focus, `(0,1)`, has `det = -3(l+2)`, a saddle at
   `l ≈ -1.1835`. The configuration is **one center and three saddles**, directly
   contradicting `audit/claude_route4a_normal_form.py`'s reported "symmetric pair
   of antisaddles".
6. **"Three parameters cannot make five cycles" is not an argument.**
   `council/notes_reversible.md` supplies the counterexample: a one-dimensional
   analytic span containing `∏_{j=1}^5 (h - j/6)` has five simple zeros.
   **Dimension is not a zero bound.** This kills a whole class of reasoning both
   councils were leaning on, including parts of §0.3 of the council itself.
7. **The divergence exponential is not `P'`.** `council/notes_kkl_hostile.md` §4:
   `P'(r) = Q(sigma(r))/Q(sigma(R(r))) · exp ∫ div F dt`. Using the bare
   exponential away from a fixed point "corrupts fold location and derivative
   certification". This was then fixed in `compact_return.py` (§2(a)).
8. **`corner_map.log` and `corner_derivatives.log` are partial runs presented as
   controls.** The corner map's scan **ends mid-file** ("scan stopped here;
   remaining eps rows not computed"), and the derivatives log contains **exactly
   one line, `kappa = 8.5`**, out of nine intended values.
9. **`claude_check_large_kappa.py` shows only logarithmic convergence.** At
   `lambda = 0.5` the claimed limit is 23.755 while scaled values climb
   7.33 → 15.50 over `eps = 1e-2 … 1e-7`. Sign is right everywhere; magnitude has
   not converged. `CLAUDE_AUDIT_ASTRA_1_3.md` flags "sufficiently small eps is
   not quantified".
10. **The KKL 206-call campaign optimized the wrong object.**
    `review_2026_09_05/COVERAGE_AUDIT.md` §2 is unsparing: the continuation drove
    the *existing* origin cycle's multiplier toward one instead of seeking a new
    coexisting pair; `continue_path.py` literally raises
    `'near fold: switch to augmented fold equations'` and **no such solver
    exists**; the single stationary point found (`r ≈ 28.174`, `D ≈ +0.2427`,
    `D_rr < 0`) is a **positive maximum inside** the known stable root, not a fold.
11. **The remote-section completeness assumption was unproved.** The remote focus
    lies above `y=0`; four inward returns from that line cannot exclude a
    surrounding cycle that never meets it. Repaired analytically by the nullcline
    section theorem (§3.2), but the repair was never implemented.
12. **A "six-cycle" result was a float64 artifact.** See §2(a), last paragraph.
13. **Several exact/asymptotic results are non-effective.** `kkl/notes_lienard.md`
    §8 and `notes_local_unfolding.md` prove existence of a neighbourhood radius
    `delta`, `eps`, `rho_0` and never quantify any of them; `theory_notes.md`
    §3/§6/§10 and the `c=1/3` law are explicitly *conditional* on unproved
    uniform joint-limit expansions.

### 3.4 Binding status, as the repository itself states it

From `STATUS.md`, and I found nothing anywhere that contradicts it:

- `H(2) >= 4` (Shi 1979/80, Chen–Wang 1979).
- **No finite uniform upper bound for `H(2)` is known. Even `H(2) < infinity` is
  open** (Gasull–Santana 2025). Bamón's theorem is *pointwise* finiteness for one
  fixed field, which is strictly weaker.
- `H(2) = 4` is a **conjecture**, not a theorem.
- The best certified frontier is Galias–Tucker (2022): the Songling instance at
  `lambda = -10^-200` has **exactly four**.
- Distribution (Zegeling 2024, Thm 1.2, filling gaps in Zhang 2002): only `(n,0)`
  or `(n,1)`. So five must be **(5,0) or (4,1)**; `(3,2)`, `(3,3)`, `(2,2)`,
  `(3,1,1)`, `(2,2,1)` and all three-nest arrangements are **excluded**. Thm 5.4:
  with four real finite equilibria, only `(n,0)` or `(1,1)`.
- **`(5,1)` — six cycles — is permitted and would also refute `H(2)=4`.** Two
  documents flag that discarding a six-cycle configuration would be a mistake.
- Bautin: at most three small cycles from one quadratic focus or center. So a
  fourth in the same nest must be **global**.

---

## 4. ROUTES

Every construction route attempted, in the order it was tried. "Repo verdict" is
what the repository concluded; "assessment" is mine.

### R1. Q4 — five distinct interior zeros of the Abelian integral (Attack 1)
**What:** Żołądek's codimension-four quadratic center. Zhao (2011) proves the
first nonzero Melnikov function has at most five zeros counted with multiplicity
and constructs three; the conjectured sharp bound is three. Five simple zeros
plus the realization gate would give five cycles. Five strikes were run.
**What happened:** Strikes 1–3 built the universal chart, the Stieltjes/ECT
structure, the lobe cell, the exact Green/PF reconstruction, and the necessary
thresholds `tau_1 > 5/11`, `kappa > 21636/19043`; the certified lobe box was
proved and then excluded by its own threshold. Eight tuned shots reached the
four-crossing S1 pattern numerically and **all eight failed S2** (first Green
maximum negative, `-0.00236` to `-0.00394`). Strike 4 proved Theorem N.
**Why it stopped:** closed by theorem. `Phi_a(tau_1) < 0` contradicts the
necessary condition (N1). At most four distinct interior zeros globally, three
on the strict lobe region.
**Repo verdict:** CLOSED. **Assessment: justified.** The proof is analytic,
independently audited, and the audit specifically hunted for a hidden uniform
asymptotic, a non-strict inequality surviving a limit, a wrong cardinal-interpolant
sign, and a reversed kernel comparison — and found none.

### R2. Q4 — the outside-lobe four-zero question (Strikes 5, 6, 7)
**What:** not a counterexample route. It is the pure-mathematics remainder:
does a Q4 integral with two interior primitive zeros admit four distinct
interior zeros? Settling it negatively would prove the Gavrilov–Iliev/Zhao
conjecture.
**What happened:** Theorem T reduced the problem to a one-parameter interval
`0 < lambda < lambda_c` on each two-anchor fibre plus one determinant condition
`K(r) > 0`. Strike 6 (branch) added the exclusion wedge
`r > 0.82052`, `kappa > 2.899241...` with exact rational interval certificates,
and reduced the second anchor away: nonpositivity at the two boundary choices
`s -> r+` and `s -> 1-` suffices. Strike 7 (branch) proved both boundary
determinants strictly negative **at `a=1`**, giving an ECT structure and a
three-zero bound on the limiting face, plus a two-zero bound on an explicit
3-dimensional subspace at every finite `a`.
**Why it stopped:** the two finite-lift boundary inequalities
`D_c(a,r) <= 0` and `D_e(a,r) <= 0` on
`1-(7/22)^{3/2} < r < 1`, `1-1/kappa_* < a < 1` are **unproved**. The obvious
repair — a Duhamel comparison dropping the negative double-integral term —
**fails**: the upper bound goes positive near the joint corner while the actual
determinant stays negative (frozen in `q4/seventh/explore_upper.json`:
`D = -1.336e-5` vs `upper = +3.543e-5` at `r=.99999, a=1.0`).
**Repo verdict:** NARROWED, STILL OPEN. **Assessment: justified, and this is the
cleanest remaining piece of publishable mathematics in the repository.** 100
diagnostic samples all have the right sign; nobody can prove it.

### R3. Q4 endpoint / two-saddle infinity graphic (Fable's lane)
**What:** the closed annulus. Three interior Melnikov zeros plus two cycles born
at the boundary graphic (the Dumortier–Roussarie alien mechanism, which is proved
to exist precisely for two-saddle cycles and proved impossible for one-saddle
loops) would give `(5,0)`.
**What happened:** the graphic was identified exactly (`claude_q4_boundary_level.py`):
in original coordinates the Q4 annulus boundary is a **finite heteroclinic orbit
between two infinity saddles**, closed by the equator arc avoiding the node.
Slopes `1 +- sqrt 2` at `rho=1`, eigenvalues `(8,-4)`, hyperbolicity ratios 2 and
1/2, product one, each 1:2 resonant. The first-order compatibility face `X(1)=0`
was found to be crossed for nine of ten lobe triples at `a* in (0.87,0.99)`, with
the interior count going 0 -> 1 across it. 546 two-root samples: counts 0 (315)
and 1 (231), **never 2, 3 or 4**. The near-center hierarchy gives exactly two
small interior zeros plus the boundary zero, closed count three. The unique
triple-center point `(A,B,eta) = (1,-17/12,0)` has `X(1;a) > 0` for all `a`, so
**no lift makes a triple center coexist with a boundary zero**.
**Why it stopped:** the lane concluded "the route survives only as a two-alien
hypothesis". That conclusion was then **ruled premature** by
`review_2026_09_05/FRONTIER_AUDIT.md` §3 on five separate grounds (see §3.3
item 3). Independently, `Q4_GRAPHICS_AUDIT.md` shows the Gavrilov–Iliev 2015
closed-annulus bound for a nondegenerate Hamiltonian two-saddle annulus is **3**,
equal to the parameter count, and the additive "two interior plus three endpoint"
loophole is obstructed within that exact family. Night lane F19 then tested it
directly: fifteen constructed three-zero Q4 directions, real cycles counted to
radius `e^40`, **exactly three every time, no alien, no four-zero direction in
20k targeted trials**.
**Repo verdict:** UNKNOWN, surviving only as a two-alien hypothesis.
**Assessment: the repo is right to call it unknown rather than closed**, but the
accumulated evidence against it is now substantial and the honest reading is that
the route is dying rather than open.

### R4. Attack 2 / Lane C — Shi third-order weak focus plus outer separatrix (`4+1`)
**What:** unfold the order-3 focus for three small cycles, add a fourth from a
homoclinic loop in the same nest, keep the remote cycle.
**What happened:** **topologically defective at its stated seed.** On the stratum
`m=5a, b=3l+5` the extra equilibria have discriminant `12(3l+5)^2(3a^2-l^2-2l)`,
so finite saddles exist only for `-1-sqrt(1+3a^2) < l < -1+sqrt(1+3a^2)`; the
entire Attack-2 box `l in [-12,-8]`, `a in [4/5,6/5]` **has none**. Independently
and exactly: `N=(0,1)` being a focus forces `25a^2+12(l+2)<0`, hence
`3a^2-l(l+2) < 0`, so **having the second focus forces the extra pair to be
nonreal**. Where finite saddles do exist, the homoclinic loop around the origin
occurs **precisely on the center curve** `eta_3=0`, verified to `4e-15` at four
values of `a` — i.e. the only way to have both an order-3 focus and a surrounding
finite loop is at a center.
**Why it stopped:** the configuration does not exist in quadratic systems.
**Repo verdict:** NEEDS REFORMULATION / closed as stated.
**Assessment: justified, and the exact algebraic version is decisive.** The
numerical loop-at-center coincidence is only to machine precision at four points,
but the finite-equilibrium obstruction is exact and independent.

### R5. Route 4b — order-3 focus plus an infinity graphic
**What:** if no finite loop, use a boundary graphic through infinity of
cyclicity two.
**What happened:** three computations. In the two-foci region there is exactly one
real invariant direction at infinity, so any graphic through both antipodal
saddles is automatically neutral at first order. The N/S separatrix splitting
`D = y_N - y_S` is **strictly positive throughout**, decaying to zero only as
`l -> -infinity` (`+4.9e-4` at `l=-30`, `+2.1e-2` at `l=-6`) — no connection. Both
unstable branches of `(0,1)` escape to infinity, so no loop there either.
**Why it stopped:** first closed numerically on a conjecture; then **upgraded to a
literature theorem** on the afternoon branch — Llibre–Schlomiuk 2004 Theorem
16(III): the only QW3 graphic portraits are W13, W15, W18, each with one graphic
surrounding the **strong** focus, all vertices at infinity.
**Repo verdict:** CLOSED. **Assessment: justified now**, though the afternoon
audit correctly warns that the same paper's Theorem 12 flags its connection set
`G3` as *numerically observed*, so the bifurcation diagram must not be cited as
an exact component enumeration.

### R6. Route 4a — reversible center `Q3^R` with a loop and a second focus
**What:** the loop point on the Shi stratum is a reversible center with a
non-rational Darboux first integral `H = (1+kX)^{2q/k}[Y^2+G(X)]`, exponent
`2q/k` irrational (0.0825 at `a=1`). Gavrilov: "almost nothing is known about the
generic reversible case." The plan was to run the Q4 playbook on it.
**What happened:** **the premise is false.** `council/notes_reversible.md` proves
exactly that the advertised external foci are **saddles**: `E_± = (-1/k,
±sqrt((k-p)/(k^2 q)))` has a triangular Jacobian with `det = -2(k-p)/k < 0`; in
the Shi chart the advertised second focus `(0,1)` has `det = -3(l+2)`, a saddle at
`l ≈ -1.1835`. The configuration is **one center and three saddles**. Since
Zegeling 5.4 forbids `(4,1)` with four real equilibria, the route would need
`(5,0)` — five cycles in one closed annulus of a codimension-3 center with three
unfolding directions.
**Why it stopped:** Fable **withdrew it** in the council. Night lane F6 then
measured the first-order Melnikov span at the loop point: **dimension three,
maximum two zeros over 200k random directions**.
**Repo verdict:** withdrawn as a primary lane, kept as a backup for its exact loop
and transferable tools. **Assessment: justified — but note the council also
records the counter-argument that "three parameters cannot make five" is *not* a
theorem** (a one-dimensional span containing `∏(h-j/6)` has five zeros), so the
route is dead on its geometry, not on its dimension count.

### R7. Reversible two-center reseed
**What:** re-seed from `ẋ = -Y(1+kX)`, `Ẏ = X+pX^2+qY^2`, affinely reduced to
`P_0 = (b-2)/4 + (1-b)y + ax^2 + by^2`, `Q_0 = -2xy`, and search for a shared
perturbation with five Melnikov zeros across the two annuli.
**What happened:** an exact structural exclusion — with `0<b<2, a<=0` there are
exactly two finite equilibria, both centers, **and no finite saddle at all** — so
the old saddle-loop mechanism cannot be repaired in this geometry. The chart's
local-completeness determinant `-16(a+2)` revealed a genuine degeneracy at
`a=-2`, repaired by swapping `eps2 xy` for `gamma x^2` (`det = -48b`); the lost
direction has **identically zero first Melnikov function in both annuli**. A
64-shape moment search found **no five-cycle candidate**; the best anywhere is
`(3,1)`. A rational four-cycle arc was produced and independently checked at 65
digits and by original-field integration.
**Why it stopped:** ran out of the bounded allowance without a candidate.
**Repo verdict:** *"The full reversible route remains open; this is not a proof
of H(2)=4."* **Assessment: correct and appropriately modest.** The remaining
task is a simultaneous zero/cyclicity bound for the **closed** annuli with
parameter-uniform endpoint control, using the repaired chart (9) at `a=-2`.

### R8. KKL Hopf completion (the primary lane, three campaigns)
**What:** in `ẋ = y+x^2+xy`, `ẏ = -10x^2+(11/5)xy+cy^2+alpha x+beta y`, find a
`beta=0`, `K>0` field with three hyperbolic origin cycles in S/U/S order **plus**
a remote unstable cycle; then a small `beta<0` adds an inner Hopf cycle, giving
`(4,1)` = five. The conditional Hopf step is sound (`l1 = K/(8 omega^3) > 0`,
`∂_beta Re λ = 1/2`); **the unproved premise is the four-cycle precursor.**

*Campaign 1 — the 206-call pilot (`kkl/`, `STRIKE5_PRECURSOR.md`).* Followed two
known cycles from `c=7/10` to `c=0.9301` at `K=1/64`; hit the artificial remote
section cap `2^20`; found one **positive maximum** at `r ≈ 28.174`, `D ≈ +0.2427`,
which is not a fold. Exact by-products: the `K=J=0` double-center obstruction, the
remote-Hopf-supercriticality obstruction, the Liénard multiplier integral, the
`N+λA` non-existence. **`review_2026_09_05/COVERAGE_AUDIT.md` is scathing**: the
campaign optimized the *old* cycle's multiplier instead of seeking a new pair,
and the augmented fold solver its own code calls for was never written.

*Campaign 2 — the staged strike (`staged_2026_09_05/`, 550 calls).* Built the
augmented Newton and **found a numerical finite fold**: `K=1/512`,
`c ≈ 0.9688884793906646`, `r ≈ 6.949087993605231`, normalized residuals
`(-1.8, 1.9)e-14`, with an independent Cartesian replay giving
`log(R/r) = -2.02e-11` and the pair-side sign pattern `+,-,+` at exact rational
coefficients. It also showed the old cutoff was **an artifact**: both cycles
persist to `c=0.9683` with the remote x-extremum at `-6.54e9`.
**But the fold fails its own precursor gate**: at that field `K_H ≈ 0.03836 > K =
0.001953125`, so the remote equilibrium has **positive** trace — an unstable
focus where a stable one is required. Exactly checked in `validate_artifacts.py`.

*Campaign 3 — the fold surface and closure (`fold_surface_`, `fold_closure_`,
3297+19 calls).* Continued the fold on the positive-K sheet to horizontal radius
**2.96e17** (`c=1.5934`, `K=7.0639`) toward a conjectural `c -> 8/5` infinity
graphic, and through the center organizer onto a separate negative-K sheet to
`m ≈ -7.17e11`, `c -> 0.3367`, toward a separate large-m organizer at `c=1/3`.
109 events, **every accepted one with exactly two root brackets** — the live
three-bracket trigger `K1_CANDIDATE_HALF.json` never fired. Then the closure
proved the **Dulac-certificate impossibility theorem** (§3.2) and exhausted the
4096-call budget.

*Campaign 4 — D1 (branches, the current frontier).* Continued three fold sheets
with binary128 (positive-center to `K=1e-10`, positive-infinity to `r ≈ 1.74e18`,
negative to `m=1e13`), 64 accepted folds, 266 field records, and — the genuine
new result — **produced exact rational 3+1 fields**: three origin cycles U/S/U
plus a remote S, at `beta = -10^-7` (branch A) and `beta = -1/25600000`
(branch B). This is precisely the object `COVERAGE_AUDIT.md` listed as never
obtained. It also caught a **six-cycle float64 false positive** and rejected it in
binary128.
**Why it stopped:** budget, and the absence of a fourth origin cycle. Max origin
count is **3**; max total is **4**.
**Repo verdict:** D1 OPEN. **Assessment: justified, and this is the live lane.**
The gap to a counterexample is exactly one cycle in the origin nest. The
repository's own refutation of the proposed shortcut is correct and important:
endpoint-sign *agreement* permits an **even** number of interior roots
(`D(u;t) = (u-a)^2 - t` births a pair with both exterior signs positive), so the
198-agreement sign map is a missing-root detector, not an obstruction.

### R9. Staged Shi / Chen–Wang trace paths
**What:** take the classical four-cycle seeds and drive the trace to zero, looking
for a double cycle (M1).
**What happened:** both fields start with three numerically resolved origin cycles
(U,S,U); the innermost collapses into the focus as `lambda -> 0`, leaving two
(S,U). `eta_1 > 0` at both endpoints (`1e-6` and `1/400`). Method quality is the
best in the tree: an **exact rational degree-8 Lyapunov polynomial** with its
identity verified symbolically, plus a fully independent 40-digit mpmath
midpoint-extrapolation integrator agreeing at the `1e-16`–`1e-18` level.
**Why it stopped:** negative for its purpose, and the strike **retracts its own
motivating premise** — trace zero plus three surrounding cycles does *not* force a
double cycle, shown by `D(r) = eps r^3(r^2-A)(r^2-B)(r^2-C)`.
**Repo verdict:** two paths do not cover the four-parameter chart; no M1.
**Assessment: justified.** The self-retraction is the right call.

### R10. Staged infinity (the KKL compactification audit)
**What:** is the fold curve terminating at an infinity graphic?
**What happened:** two exact results. Below `c=241/250` the vertical graphic is
**impossible** (the `x=-1` barrier), so the `2^20` boundary was an artifact. For
`241/250<c<1` there is a genuine candidate eigenvalue-neutrality line
`J(c)=0` at `c_* = 0.968620633553494`, independent of `K`. For `1<c<=3/2`,
neutrality would need `c=8/5`, outside the box.
**Why it stopped:** the connection, its splitting parameter, and the transition
constant `C` are **all uncomputed**; the first graphic coefficient `C-1` is
therefore unknown.
**Repo verdict:** the proposed Stage-2 kill test remains unmet.
**Assessment: justified.** Two diagnostic corrections here are worth keeping: a
cycle approaching a hyperbolic graphic **need not** have multiplier tending to
one, and a logarithmically divergent *desingularized* residence time can have
**bounded physical duration**.

### R11. K1 theory (can a first-order weak focus be surrounded by three cycles?)
**What:** a theory kill test for the whole KKL construction.
**What happened:** the source attribution was **corrected** — Zhang–Cai 1991 is
about the *strong* focus's nest; the general own-nest result is Pingguang Zhang
1999, Acta Math. Sinica 42(2), 289–304, whose **full proof could not be
retrieved** (landing-page PDF control resolved to an icon; direct retrieval HTTP
403). The accessible Zhang–Cai monotone-ratio argument was tested and **fails
exactly**: its required one-sign derivative hypothesis is false, because
`(f/g)_z = -u^{2c+1}N/(5W^2)` with `N(0)=5K>0` and leading coefficient
`d(c-1)h < 0` for `1/2<=c<1`. A new exact Bernstein amplitude restriction was
proved instead.
**Why it stopped:** K1 unresolved.
**Repo verdict:** K1 UNRESOLVED. **Assessment: justified.** The failed-hypothesis
result is real and is not evidence for three cycles.

### R12. D2 — the order-two loop law (branch)
**What:** prove `sigma * eta_2 < 0` for a homoclinic loop through a hyperbolic
saddle around an order-two weak focus. This is the missing middle rung of the
"3-k" ladder between Zhang 1999 (order two: at most one) and Li–Cherkas (order
three: none), and it would explain rigorously why every attempt to stack a
degenerate boundary on a degenerate focus collapsed onto a center.
**What happened:** Proposition A was claimed as an "if and only if" and is
**FALSE as stated** — Astra's audit supplies the exact counterexample
`ẋ = -y+xy+y^2`, `ẏ = x+x^2+xy` (a reversible center whose only other equilibria
have traces `-1` and `+1`, neither neutral) and an omitted case
`H = a^2-b(l+1) = 0`. **The corrected one-way implication is VERIFIED exactly**,
including the chart degeneracies `b=0` and `l=-1`, via explicit reflection
certificates `X(Rz)+RX(z)=0` on `C3=0`. Statement (C) is verified from
Llibre–Schlomiuk 2004. The Dulac route was closed negatively: `(1+ax+by)^k` is
sign-indefinite on every loop interior for `k in {m/a, -(2l+b)/b, 1, -1, 2}`, and
an order-two focus forces `div(BX)` to vanish to fourth order, so no
polynomial-times-power `B` of degree below four can work.
**Why it stopped:** statement (D), the component enumeration, needs
Artés–Llibre–Schlomiuk 2006 (IJBC 16, 3127–3194) — **publisher PDF returns HTTP
403**, only p. 3127 accessible. The auditor explicitly refused to invent a theorem
or portrait number. 104 sampled loops, all negative, no completeness argument.
**Repo verdict:** D2 OPEN. **Assessment: justified.** This is a rigidity result,
not a counterexample route, but it is the best-argued piece of new mathematics on
the branches.

### R13. Resonant infinity hemicycles (branch)
**What:** the Marín–Villadelprat 2025 resonant line `a=-1`, where the published
theorem gives lower bounds but **omits the sharp upper bound**. Proposed
mechanism: two compact cycles plus three hemicycle cycles.
**What happened:** **the specific mechanism is excluded.** For a sequence carrying
two compact cycles, at most **one** endpoint cycle can accompany them; if the
compact distribution is `(1,1)`, **zero**. Proof by uniform endpoint division
(`D_u = delta_u + Delta_u s^lambda + R_u` with `L=3/2 > lambda`) giving
`D_u(s_n)/N_n -> pi p != 0` for *every* `s_n -> 0`, plus Rolle on the lower side.
The third compact generator on `a=-1` is logarithmic,
`J(h) = 2 pi[h-1-2 log((h+1)/2)]`, and the resonant logarithmic invariant reduces
to the same generator `g=(b-1)eps_1`, vanishing at `b=1`. Verified symbolically
plus 70-digit quadrature with `err < 1e-55`, and a positive control
`(p,q,k)=(1,3,20)` shot in the **original quadratic field** at three `t` and two
tolerances.
**Why it stopped:** the *proposed* construction is dead. The route is not.
**Repo verdict:** does not exclude five cycles all tending to infinity, or
one-compact/four-endpoint, or a different base `(-1,b_0)`.
**Assessment: justified.** Night lanes F16/F17/F18/F18b/F18c independently
confirmed by sweep: the hemicycle emits **two cycles and never a third** in any
scaling tested (4536 fields in the physical scaling alone), and the second-order
null direction at the holomorphic point `(a,b)=(-1,1)` is **center-preserving**.

### R14. The overnight sweeps (F3–F21) — the only genuinely broad search
**What:** ten planned lanes plus eleven opened during the night, using the compiled
C engine, over the Shi chart, the KKL family, the two-center chart, the Yu–Han
reversible family, and Q4 directions.
**What happened, lane by lane:** F3/F5 Shi Sobol — max **1** at trace zero, 3
otherwise (`300032 sets, hist {0:490116, 1:8861, 2:33}`, max 2 post-fix).
F4/F14/F20 descents from 40+ four-cycle seeds, 60 generations — max **4**, "the
displacement between the outer roots is a single hump, no forming fold".
F6 reversible Melnikov span 3, max **2** zeros in 200k directions.
F11 order-two neutral loop — does not exist as a focus; 46/46 focus-type loops
have `sigma*eta_2 < 0`. F12 KKL double center — span 3, max 2. F13 Yu–Han
degenerate curve — the triple-zero element **never has an interior zero**;
`(3,1)` is the ceiling; CLOSED. F15 neutral hemicycle — mvneutral max 3
(20032 fields), kklstar max **4** (10048 fields, 24 fours, every one a `(3,1)`
with the remote cycle at radius 1e6–1e9). F16 Dulac coefficients — first order
gives one cycle, two with the ratio parameter, exactly matching the published
lower bound. F17 second-order null direction — center-preserving; CLOSED.
F18/F18b/F18c — max total **3**, and **no single nest ever held three**.
F19 Q4 alien test — span dimension 4, fifteen three-zero directions, **exactly
three real cycles every time, zero aliens, zero four-zero directions in 20k**.
F21 full Shi chart compactified, 30016 fields — max **3**, and exactly one field
ever reached it.
**Why it stopped:** the night ended. **Repo verdict:** *"A counterexample, if it
exists, is not near any integrable stratum, any degenerate graphic, or any
published four-cycle field, and it is not reachable by descent from them. The
honest prior after tonight is well under one in twenty."*
**Assessment: the *observations* are sound and the honesty is exemplary, but the
"CLOSED numerically" labels are weaker than they look** — see §8.

---

## 5. OPEN THREADS

Everything the repository leaves explicitly open, with the file that says so.

### 5.1 Mathematical

| # | Thread | Where written |
|---|---|---|
| 1 | **The two finite-lift boundary determinant signs** `D_c(a,r) <= 0` and `D_e(a,r) <= 0` on `1-(7/22)^{3/2} < r < 1`, `1-1/kappa_* < a < 1`. Sufficient (with the inherited exclusions) for a global Q4 three-zero theorem | `<astra/q4-determinant>:Q4_SIXTH_BOUNDARY_REDUCTION.md`, `q4/sixth/notes_boundary_reduction.md` |
| 2 | Extend the `a=1` ECT theorem to finite `a<1`. Explicitly does **not** follow from compact-interval continuity ("no uniform control as the anchor tends to 1 simultaneously with `a`") | `<astra/q4-determinant>:Q4_SEVENTH_LIMITING_FACE.md` §5 |
| 3 | The fourth coefficient direction `J = Y_{K_0}`: the exact `a`-connection `J_a = c(a)J - (1-t)J_t/[a(1-a)] + U` exists as an identity but is not a sign theorem; the first Duhamel upper bound is **dead** | same, `q4/seventh/parameter_connection.py`, `explore_upper.json` |
| 4 | Sign of `K(r) = P_B(r)Phi_C(r) - Phi_B(r)P_C(r)` on the remaining baselines with `P_B(r)>0` | `Q4_TWO_ROOT_REDUCTION.md` §5, `ASTRA_FIFTH_STRIKE.md` |
| 5 | Sign of the exact loop functional `I_loop(k,mu) = 9k^{-3/2} ∫ (d-km^2)^2 N_mu(m)/D(m)^3 dm` on the residual fibre | `q4/notes_fifth_loop_gate.md` |
| 6 | **D2 statement (D)**: an exhaustive component list of the loop locus `L` minus the center strata, with a justified sign on each, including chart boundaries. Blocked on an unreadable paper | `<astra/fastra-afternoon>:FABLE_D2_ORDER_TWO_LOOP.md` §13, `FASTRA_AFTERNOON_REPORT_2026_09_05.md` §2 |
| 7 | D2 Corollary A2's nongeneric crossings; the square-factor argument is generic-only | same |
| 8 | **D1**: an exhaustive origin-root count on any accepted field. Endpoint-sign agreement permits an **even** number of extra roots; no derivative sign certificate, no enclosure, no proof this is the only finite-radius birth region | `<codex/fastra-d1>:FASTRA_D1_REPORT_2026_09_05.md` "Why sign agreement cannot close D1" |
| 9 | Resonant: the broken-connection joint compensator coefficients with parameter-uniform divided remainder; pure endpoint cyclicity on `a=-1` (the source omits its own sharp upper bound); one-compact/four-endpoint; all-five-at-infinity; a different base `(-1,b_0)` | `<astra/resonant-joint>:RESONANT_JOINT_2026_09_05.md` §9 |
| 10 | KKL: the **weighted-orbit comparison** — one sign change of the multiplier integrand does not compare the orbit-dependent positive weights of distinct cycles. This is named as "a real proof obligation" | `KKL_FOLD_CLOSURE.md`, `fold_closure_2026_09_05/theory_obstructions.md` |
| 11 | KKL: the exhaustive connected fold component; every possible origin root; endpoint certificates; the graphic connection and its transition constant `C` at `c_*` and at `c=8/5` | `KKL_FOLD_SURFACE_STRIKE.md`, `STAGED_INFINITY_2026_09_05.md` §4 |
| 12 | KKL: quantify the neighbourhood radii `delta`, `eps`, `rho_0` in the local no-collapse results — all three are proved to exist and never bounded | `kkl/notes_lienard.md` §8, `kkl/notes_local_unfolding.md` |
| 13 | KKL: whether a remote unstable cycle can exist in the attracting-remote-focus region with `J<=0`. The Hopf sign blocks only its *local* birth; a Dulac or return-map comparison is needed for a finite-amplitude exclusion | `STRIKE5_PRECURSOR.md` |
| 14 | Reversible: a simultaneous zero/cyclicity bound for the **closed** annuli with parameter-uniform endpoint control, using the repaired chart (9) at `a=-2`; and the label of the `(-2,1)` object in the 121-graphic catalogue | `REVERSIBLE_RESEED_2026_09_05.md` §6 |
| 15 | K1: the general order-one/order-two own-nest theorem. Zhang 1999's full proof was never retrieved | `STAGED_K1_THEORY_2026_09_05.md` §1 |

### 5.2 Engineering / verification

| # | Thread | Where |
|---|---|---|
| 16 | **Implement an interval ODE return-map verifier.** Recommended as the next step in at least four separate documents and never done. `H16_CERTIFICATION_PLAN.md` names CAPD::DynSys + MPFR, gives seven algorithm steps, precision guidance and five controls including a deliberate negative control | `H16_CERTIFICATION_PLAN.md`, `COVERAGE_AUDIT.md` §7, `FASTRA_ZOOM_OUT_2026_09_05.md` §5.6, `STAGED_RUN_2026_09_05.md` |
| 17 | **Certify one known control first** — the KKL four-cycle, or the conditioned GT field `ẋ = -10^{-24}x-y-10x^2+(499/100)xy+y^2`, `ẏ = x+x^2-(311375001/12500000)xy` — "so discovery does not outrun the method needed to prove a candidate" | `FASTRA_ZOOM_OUT_2026_09_05.md` §5.6, `FRONTIER_AUDIT.md` §6C |
| 18 | **Implement and benchmark the nullcline section theorem.** It is proved and exactly checked; its evaluator was never written, no control was transported, and the old radius caps were never translated into the new coordinate | `review_2026_09_05/KKL_SECTION_REPAIR.md` |
| 19 | Run the specified 64-call pilot: at `c=33/40`, bracket `K` so the remote U has old section coordinate `-2^15`, then continue that constraint jointly with `D_r=0`. **Never executed** | `review_2026_09_05/KKL_NEXT_CONSTRUCTION.md` |
| 20 | The **finite-beta incumbent scope amendment**: keep all four incumbent cycles at `beta>0, K<0` and hunt an additional finite pair, giving `(5,1)` = six, which also refutes `H(2)=4` and is permitted by Zegeling. Suggested exact box `c in [69/100,71/100]`, `alpha in [-74,-72]`, `beta in [1/1000,1/500]`. **Requires a scope amendment and was never authorized** | `FASTRA_ZOOM_OUT_2026_09_05.md` §8, `KKL_NEXT_CONSTRUCTION.md` §5 |
| 21 | **Eight of twenty engine review items are still live**, including one M-rated correctness bug in the C crossing logic (A3 dead-code disjunct, E1 bisection with no bracket-validity check and `gnew==0` collapse) and E4, whose broken near-miss gate `rad > 3*scale` **mis-steered every F4/F14/F20 evolutionary descent**. Also unfixed: D2 (a descending radius grid inverts every S/U label), C1 (arithmetic midpoint of a geometric bracket), A4, A5, E2, E3, C2, D3, D4 | `<astra/fastra-afternoon>:audit/fable_engine/REVIEW_engine.md` §G items 6, 7, 8 |
| 22 | The review's §F two-source noise rule was only half-implemented, and **only in the compactified path**. Cartesian `count = 2` records (33 in F3, 445 in `W_kklx_L3`) were **never re-audited** — `recount_hits.py` skips anything with `total < 3` | same, `recount_hits.py:9` |
| 23 | **B4: "This validation claim should be re-derived or retired."** The reviewer could not reproduce the "remote cycle near 3711" on the default ray (`D` smooth and strictly negative for `r in [16,1e5]`, verified against an independent scipy integration). `FABLE_LANES` later restates a *different*, rotated-ray version. Never explicitly retired | same |
| 24 | Ten queued sweep lanes have **no output at all** and no report acknowledges them: `F5_shi_full_L4` (150k sets), `F5_kklx` (60k), `F7_q4pert`, `F8_yz_evolve`, `F8_shi_evolve`, `F9_mv` (40k), `F9_mvpert`, `F6_q3rpert`, `F3b_lam0_store1`, `F3c_lam0_evolve`. `data/QUEUE_DONE` was never written | `audit/fable_engine/queue*.sh` vs `data/` |
| 25 | `data/push_fold_result.json` has `"min_q": null` — the fold descent's objective returned `None` and terminated without crossing zero. **Reported nowhere** | same |
| 26 | The `D2_cross_a-1.5` and `D2_cross_a-2` crossing bisections **both crashed** (`TypeError: 'NoneType' object is not subscriptable`); §10 substitutes hand interpolation | `audit/fable_engine/data/D2_cross_*.log` |
| 27 | Astra's binary128 handoff `full_return128` is built, validated on six analytic controls, and **has never been used for a sweep** | `FASTRA_AFTERNOON_REPORT_2026_09_05.md` §4 |
| 28 | Merge debt: branch A (`codex/fastra-d1`) and branch B (`astra/fastra-d1`) are two **disjoint** D1 campaigns sharing a directory name; A is later and larger but **B's `counter_check/` is not superseded** and contains the only diagnosis of the NaN-acceptance bug that A's own campaign ran on. `astra/q4-determinant` and `astra/resonant-joint` conflict on three files (`README.md`, `STATUS.md`, `Q4_TWO_ROOT_REDUCTION.md`) | — |

### 5.3 Budget state

The KKL/Shi shared ledger is **exhausted: 4096/4096** (`756 + 3297 + 24 + 19`),
per `fold_closure_2026_09_05/accounting.json` and `KKL_FOLD_CLOSURE.md`. The D1
branch work is a separate, unreconciled ledger (64 folds, 266/176 field records,
274 shooting calls). Any resumed campaign needs a new budget decision.

---

## 6. LITERATURE

`LITERATURE_AUDIT.md` and `SOURCES.md` (29 entries) on `main`, plus
`HISTORICAL_FIVE_CYCLE_CLAIMS.md`, `FOUR_CYCLE_SEED_LEDGER.md` and
`STATUS_SOURCE_AUDIT.md` on the orphan branch.

### 6.1 Theorems the campaign relies on

| Result | Source | Used for |
|---|---|---|
| `H(2) >= 4` | Shi, Scientia Sinica 23 (1980) 153–158; Chen–Wang, Acta Math. Sinica 22 (1979) 751–758 | The lower bound |
| **Even `H(2) < infinity` is unknown** | Gasull–Santana, PAMS 153 (2025) 669–677, arXiv:2407.13465; Villanueva–Tucker arXiv:2602.22558v2 | The framing: `H(2)=4` is a conjecture |
| Pointwise finiteness for a fixed quadratic field | Bamón, Publ. IHÉS 64 (1986) 111–142 | Explicitly distinguished from a uniform bound. **Independent of the challenged Ilyashenko proof** |
| Bautin cyclicity of one quadratic focus/center is exactly 3 | Bautin, Mat. Sb. 30 (1952) | The fourth cycle in a nest must be global |
| **Distributions are `(n,0)` or `(n,1)` only** | **Zegeling, Adv. Nonlinear Anal. 13 (2024) 20240012, Thm 1.2** (filling gaps in Zhang, QTDS 3 (2002) 437–463) | Five must be `(5,0)` or `(4,1)`. `(3,2)`, `(3,3)`, `(2,2)`, `(3,1,1)`, `(2,2,1)`, three nests all **excluded**. Thm 5.4: with four real finite equilibria, `(n,0)` or `(1,1)` |
| Songling has exactly four, with rigorous localization | **Galias–Tucker, Appl. Math. Comput. 415 (2022) 126691** | The only accepted computer-assisted proof in the area; the certification benchmark |
| QW3 has 18 portraits; graphics only in W13, W15, W18, each surrounding the **strong** focus | **Llibre–Schlomiuk, Canad. J. Math. 56 (2004) 310–343, Thm 16(III)** | Closes route 4b / statement (C) |
| Q4 Melnikov: upper bound 8 | Gavrilov–Iliev, JMAA 357 (2009) 69–76, arXiv:0811.4602 | The Q4 normal form, double cover, four independent coefficients |
| Q4: upper bound 5, lower bound 3, conjectured sharp 3 | **Zhao, Nonlinearity 24 (2011) 2505–2522, arXiv:1011.2253** | The whole Q4 lane. **Two sign errors found and corrected**: Theorem 14/Corollary 18's strip endpoint, and equation (24)'s reconstruction forcing |
| Classification of higher-order first nonzero Melnikov functions | Iliev, Bull. Sci. Math. 122 (1998) 107–161 | Invoked by both Q4 papers |
| Nondegenerate Hamiltonian two-saddle closed annulus has cyclicity 3 | Gavrilov–Iliev 2015 | Obstructs the additive "2 interior + 3 endpoint" loophole |
| Alien cycles require a polycycle with at least two saddles; they do not occur at a single saddle loop | Caubergh–Dumortier–Roussarie; Gavrilov | Lane B's original argument (later withdrawn for a different reason) |
| Homoclinic loop cyclicity of a quadratic integrable non-Hamiltonian system is 2, with one exceptional case | Han, Sci. China A 40 (1997) | Lane B. **The exceptional case was never identified from the abstract** — recorded as an open caveat |
| Leontovich–Roussarie loop order; first four Dulac coefficients | Han–Yang–Tarta–Gao, JDDE 2008 | Lane B's `c_0 + c_1 w log w + ...` expansion |
| Hemicycle cyclicity 2; simultaneous cyclicity at most 3 off resonance; **upper bound omitted at `a=-1`** | **Marín–Villadelprat, arXiv:2501.16924v1, Thms B–D** | The resonant route. The omitted bound *is* the gap being attacked |
| Reversible center `zdot=-iz(1+a zbar)` under quadratic perturbation: annulus bound 2 | Li, QTDS 6 (2005) 205–215 | Mechanism ledger |
| Global finiteness for the quadratic infinitesimal problem | Gavrilov, Invent. Math. 143 (2001) 449–497 | Mechanism ledger |
| Explicit bound away from centers/singular fields/infinity | Ilyashenko–Llibre, arXiv:0910.3443 | Mechanism ledger; degenerates exactly where a search is hardest |
| Four cycles in a near-integrable quadratic family | Yu–Han, IJBC 22 (2012) 1250254, arXiv:1002.1055 | Seed YH |
| Concrete numerical realizations | Yu–Zeng, IJBC 31 (2021), arXiv:2002.09987 | Seeds; explicitly distinguishes visualization from proof |
| Five-coefficient family with four visible cycles | Kuznetsov–Kuznetsova–Leonov, DEDS 21 (2013) 29–34 | The KKL seed and the primary lane's family |
| CAPD::DynSys | Kapela et al., CNSNS 101 (2021) 105578; Chaos Solitons Fractals 157 (2022) 112125 | The (never-installed) certification toolchain |

### 6.2 What the repository explicitly could NOT verify

This list matters more than the previous one.

1. **Pingguang Zhang 1999**, *On the Uniqueness of the Limit Cycle of the
   Quadratic System with a 2nd-Order Weak Focus*, Acta Math. Sinica 42(2)
   289–304 — **full proof not retrieved.** The landing-page PDF control resolved
   to an icon; direct retrieval returned **HTTP 403**. Only the indexed abstract
   was read. The reviewer explicitly refused to reconstruct the proof from the
   abstract. This is the theorem the whole "3-k ladder" rests on.
   (`STAGED_K1_THEORY_2026_09_05.md` §1)
2. **Artés–Llibre–Schlomiuk 2006**, IJBC 16, 3127–3194 — **HTTP 403**; only
   p. 3127 accessible; the author's supplementary page supplies no article. Its
   advertised 373 parameter subsets and 126 portraits (95 with an order-two
   focus) are exactly what statement (D) of the D2 conjecture needs. The auditor
   refused to cite a theorem or portrait number from an unread catalogue.
   (`FASTRA_AFTERNOON_REPORT_2026_09_05.md` §2)
3. **Han's exceptional quadratic integrable case** — never identified from the
   abstract. If Q4's loop were that case, Lane B's loop-cyclicity bound changes.
   (`CLAUDE_LANES_B_C.md`)
4. **Shi's original 1978 five-cycle preprint** — never located. Only his 1990
   retrospective and the corrected family were recovered.
5. **Shi's 1988 rebuttal to Qin Yuanxun**, BLMS 20, 597–599 — abstract retrieved,
   the three pages of coefficients not.
6. **Chen–Wang's original Chinese full text** — not retrieved; the coefficients
   used are attributed to Yu–Zeng's later reconstruction.
7. **Galias–Tucker** — no maintained public one-command replay archive was
   located. A future campaign must contact the authors or reimplement.
8. **Hernandez Rosales, August 2026, "H(2)=7"** — ResearchGate preprint, "no file
   available", zero citations, no coefficients, no DOI/arXiv. Deliberately
   recorded as **claimed/unverified, not refuted**.
9. **Pedregal arXiv:2103.07193v3** (active, claims `H(2)=4` as Corollary 1.5) —
   the central variational theorem was not verified. The audit records a concrete
   defect: Theorem 1.4 as printed gives an erroneous component count already for
   a line and a hyperbola. Called "a concrete defect as written, not a completed
   refutation".
10. **Gaiko arXiv:math/0611142** — no dedicated primary rebuttal located. The
    audit flags two unproved bridges (the four-distinct-to-multiplicity-four
    merger, and invoking a focus-local Bautin bound where the termination
    alternative is a separatrix cycle) but does not claim a formal disproof.
11. The 2023 six-cycle antecedent and the interior coefficients of the 2024 Riga
    proceedings.

### 6.3 Claims the repository actively refuted with its own algebra

- **Malyarets et al., Mathematica Montisnigri 65 (2026) 23–35, six cycles** —
  published, and the audit shows the proof is invalid as written: the two
  trace-zero conditions are `lambda=0` and `lambda=-5a/n`, which are **different
  vector fields** unless `a=0`; at `a=0` the system is reversible under
  `(x,y,t) -> (-x,y,-t)` so both monodromic equilibria are **centers**, and the
  paper's own `V_3` expression carries a factor `a`. Second independent error:
  a nonzero third Lyapunov quantity gives cyclicity under an unfolding, not
  three present orbits at the degenerate parameter. **No independent published
  rebuttal exists** — this is the repository's own finding.
- **Voronin–Lebedev 2022, six cycles (3:3)** — the argument adds local counts
  from two incompatible values of one coefficient. Also excluded a priori by
  Zegeling: `(3,3)` is not an admissible distribution.
- **Da Silva–Vieira–Leonel, Entropy 26 (2024) 745**, `H(n)=2(n-1)(4(n-1)-2)` —
  refuted by Buzzi–Novaes arXiv:2411.09594; no retraction located. Neither
  establishes `H(2)>4`.
- **Zhao's Theorem 14 / Corollary 18 sign**, and **Zhao equation (24)'s forcing
  sign** — both corrected, the latter with a `4.54e-9 -> 2.17e-19` numerical
  signature.
- **Galias–Tucker Lemma 2's sign labels** — all four endpoint pairs match
  `y - P(y)` rather than the stated `P(y) - y`, checked against Lemma 3's
  multipliers and Figure 1. Verdict: probable transcription error; **the
  bracketing proof is unaffected**. Take stability from Lemma 3.
- The README's own earlier claim that Shi corrected the "**fifth** focus
  quantity" — the source says the **third** (`V_3`).

### 6.4 Scope exclusions the repository enforces (correctly)

Piecewise/discontinuous quadratic systems, 3-D quadratic systems, cubic or
higher-degree fields, and perturbations of degree greater than two are all
excluded as evidence about `H(2)`. Several recurring search hits are logged as
false positives on exactly these grounds (Yang–Liang–Zhang 2010 and Tian–Yu 2014
are 3-D; Lu arXiv:2607.27464 is piecewise-smooth slow-fast; Eshkobilov et al.
arXiv:2604.12883 raises degree to 14/29/31/39).

---

## 7. INFRASTRUCTURE

### 7.1 Two execution environments, neither of them this Mac

The repository was written by two agents ("Astra" = Codex, "Fable" = Claude)
working in parallel on different machines, plus cloud workers.

- **Fable's cloud sandbox: `/home/user/H16P`.** Hard-coded in
  `audit/fable_engine/queue*.sh` (`cd /home/user/H16P/audit/fable_engine`), in
  `audit/fable_d2_center_crossings.py:5` and `fable_d2_crossing_refine.py:4`
  (`spec_from_file_location('f11', '/home/user/H16P/audit/fable_f11_neutral_loop.py')`),
  and in `audit/fable_q3r_melnikov.py:65`
  (`np.save('/home/user/H16P/audit/fable_q3r_basis.npy', ...)`).
  **These paths will not run on this machine.**
- **Astra's sandbox: `/workspace/scratch/97757d9f13f6/H16P`**, leaked once in
  `reversible_reseed/data/moment_search.log:55`.
- The paused Strike-5 checkpoint records a third path,
  `/Users/scottg/Documents/Codex/2026-09-04/he/work/H16P`.
- `main`'s own directories contain **zero** references to `/home/user`, cron,
  night watch, `fable_engine`, or `CLUES.md`. All of that lives only on the
  `astra/fastra-afternoon` and `claude/conjecture-progress-report` branches.

### 7.2 There is no cron. The "night watch" is bash plus a polling loop.

Seven scripts `audit/fable_engine/queue.sh` … `queue7.sh`. Each `cd`s to the
sandbox path and runs a fixed sequence of `evolve.py` / `sweep_family.py` /
`sweep_shi.py` / `sweep_log.py` / `evolve_log.py` / `recount_hits.py` jobs, each
redirecting into a per-lane log. Dependency edges are two `pgrep` spin-waits:

```bash
while pgrep -f "sweep_shi.py data/F3_lam0_L4" > /dev/null; do sleep 30; done   # queue.sh
while pgrep -f "^python3 sweep_log.py kklstar" > /dev/null; do sleep 30; done  # queue6.sh
```

Completion is signalled by sentinel files. `data/F15_QUEUE_DONE`,
`data/F20_DONE` and `data/QUEUE7_DONE` exist and contain `DONE`;
**`data/QUEUE_DONE` does not exist** — queues 1–4 never finished. All six
`data/queue*.log` files are **0 bytes** by construction. Queues 2/3/4 are
near-duplicates of queue 1 with the spin-wait removed and jobs reordered:
evidence of repeated restarts.

Separately, four `origin/fable/compute-*` branches are **cloud worker outputs**:
`compute-f3-lam0` (a 23142-line `F3_lam0_L8_worker.jsonl`), `compute-f5-shi`,
`compute-evolve`, `compute-pert`. Each commits "worker: … results (partial)"
repeatedly — a poor-man's checkpointing scheme, committing data snapshots to git.

### 7.3 The night watch ledger — `audit/fable_engine/data/CLUES.md`

291 lines, titled *"Night Watch Log - Planar Quadratic Vector Field Search"*.
It ran 2026-09-05 04:42–08:24 UTC, nominally 12 iterations at 20-minute
intervals (~220 min elapsed). Each iteration snapshots the last histogram lines
of the live `.log` files and counts records with `total >= 4`.

Its own conclusion, verbatim:

> **Max total seen: 4 — FIVE CANDIDATE found: NO.** … "Max Total Cycles: 4
> (search never reached 5; 13 of 12 iterations found max_total=4) … Anomalies
> Detected: Iterations 11-12 showed max_total=0 …, suggesting potential data
> purge or search transition … FIVE CANDIDATE Result: NOT FOUND - despite 12 full
> iterations with automated scanning of all *.jsonl files, no configuration with
> 5 or more limit cycles was discovered."

**Be blunt about this file: it is low-quality instrumentation, not a result.**
Two interleaved iteration counters write into the same file, so iteration numbers
repeat and go backwards ("13 of 12 iterations"); `max_total` oscillates
`4 -> 0 -> 4 -> 0` because the scanner was reading files other queues were
truncating; the "data purge" is guessed, not diagnosed; and the total-4 record
count dropping `1387 -> 77 -> 25` is never explained. **The actual results are in
the per-lane logs**, and those are sound. The best records the watch saw were
`F4_kkl_evolve2.jsonl` at `score = 4.0` with nests `[roots:1 stab:U]` and
`[roots:3 stab:SUS]` — i.e. the same 3+1 shape that D1 later pinned down as an
exact rational field.

### 7.4 What the ledgers contain

| Ledger | Records | Schema / content |
|---|---|---|
| `kkl/data/returns.jsonl` | 206 | `evaluation, purpose, request{r,c,alpha,beta}, result{status, return_coordinate, D, period, R_r, R_c, R_alpha, speed_derivative, derivative_discrepancy, divergence_exponential, opposite_coordinate, winding_about_focus, min_xy, max_xy, nfev, tol, cpu_seconds}`. All `NUMERICAL_ONLY`. ~6.86 evaluator CPU-seconds total |
| `staged_2026_09_05/kkl_returns.jsonl` | 400 | Adds `stage_evaluation, prior_completed_calls, section, transport, log_displacement, log_displacement_derivative, q, multiplier`. 394 completed, 6 unresolved (all `c=1.2, r=20000` guard stops) |
| `staged_2026_09_05/shi_returns.jsonl` | 150 | `index, family, l_m_a_b_exact, r, lambda, rtol, tag, validated:false, H, raw_log_return, period, divergence_integral, multiplier_at_cycle, return_derivative, nfev, parameter_source`. 139 ok, 11 failed |
| `fold_surface_2026_09_05/returns.jsonl` | **3297** | `evaluation, campaign_evaluation, purpose, request, result, source, source_sha256, recorded_at, wall_seconds` (+`dependency_sha256` on 1692 rows). Statuses: 1908 `NUMERICAL_TWO_HALF_PASSAGES`, 1240 `NUMERICAL_ONLY`, **123 `UNRESOLVED`**, 26 separatrix. Engines: 1075 `compact_return.py` (source **not in that directory**), 944 `half_ld.py`, 401 `half_m.py`, 379 `half_quad.py`, 204 `half_m_quad.py`, 164 `log_variational.py`, 52 `variational_return.py`, 26 `graphic_shoot.py`, 49 angular. Failures: 97 guard/time, 23 SIGKILL (CPU fuse), 2 horizon, 1 section gate |
| `fold_closure_2026_09_05/returns.jsonl` | 19 (the hard cap) | 16 two-half + 3 full. 12 cusp Jacobian, 4 cusp residual, 3 outermost sign check |
| `fold_surface_2026_09_05/events_*.json` | 9 files, 21 KB–1.2 MB | Corrector histories, **rejected** predictors, section data, derivatives, root brackets, stationary brackets, sampled displacement signs |
| `fold_surface_2026_09_05/EVENT_LEDGER.md` + `component_summary.json` | 109 events | 91 `ACCEPTED_NUMERICAL_FOLD`, 18 `UNRESOLVED_CORRECTION`. Positive sheet `K=1.22e-4, r≈6.948` to `K≈7.064, c≈1.593, r≈2.96e17`; negative sheet to `c≈0.3367, m≈2.77e12, r≈1.20e10`. **Every accepted row has 2 brackets except one with 1.** `global_component_complete: false`, `three_origin_candidate: false`, `five_cycle_candidate: false` |
| `<codex/fastra-d1>:root_ledger.csv` | 270 rows | Count distribution `{3: 111, 2: 99, 0: 58, 1: 2}` — **zero rows with count >= 4** |
| `<codex/fastra-d1>:sign_map.csv` | 263 rows | Counts `{3:108, 2:96, 0:58, 1:1}`; methods `{double_full_return:167, binary128_matching:96}`; edge kinds `{scan_cap:201, integration_failure:52, angular_chart_failure:8}` |
| `<astra/fastra-d1>:field_ledger.csv` | 176 rows | Origin counts `{0:60, 3:46, 2:40, 1:30}`; stabilities `{'':60, S:28, SU:26, USU:24, SUS:22, US:14, U:2}`. Plus a list of **44** (label,kind) pairs the stock counter under-reported |
| `<astra/fastra-d1>:counter_check/discrepancy_ledger.json` | 5 fields, 10 roots | The forensic counter audit. Per root: `index, stability, section_theta, approx_r, root_residual, numerical_sign_bracket, grid[2]{native, quad, inside_default_grid, visited_in_original_run}, classification` |
| `audit/fable_engine/data/*.{log,jsonl}` | ~60 files | The F-lane sweeps. Notable: `F3_lam0_L4.log` `DONE 300032 sets 6077s hist {0:490116, 1:8861, 2:33}`; `F15_kklstar_sweep.log` `10048 sets {0:3016,1:4391,2:1976,3:641,4:24}`; `F21_shi_compact.log` `30016 sets {0:21331,1:8418,2:266,3:1}`; `RECOUNT_fixed_counter.log` `DONE {'4->4': 24}` |
| `review_2026_09_05/EVIDENCE_MANIFEST.json` | 1 | `review_additional_ode_evaluations: 0`, `saved_return_records: 206`, `allowance: 4096`, `remaining: 3890`, `unchanged_ledger_sha256: 6512aee9…`, **`five_cycle_candidate: false`, `five_cycle_certificate: false`** |

### 7.5 Discipline that is worth preserving

- Every ledger is **append-only** and records **failures**, with an explicit
  distinction between `scan_cap`, `integration_failure` and
  `angular_chart_failure`. A stopped integration is `UNRESOLVED`, never an escape
  or a nonexistence proof. This is stated in a dozen places and never violated.
- Every evaluator runs single-threaded, at `nice 10`, with a `RLIMIT_CPU` fuse
  (usually 10 s) and a wall timeout. Budgets are enforced in code: `budget.py`
  raises at `PRIOR+n >= 4096`; `finalize_kkl.py` now **refuses new calls**.
- Source hashes are pinned in five manifests, all of which verify — **and all of
  which are honestly labelled as post-hoc**
  (`FINAL_FILE_HASHES_NOT_HISTORICAL_EXECUTION_HASHES`).
- The D1 branches verify that the shared production counter is **byte-identical
  to git** before using it:
  `subprocess.check_output(['git','show','fable/current:'+path]) == (ROOT/path).read_bytes()`.
- `review_2026_09_05/export_visible_transcript.py` redacts GitHub PAT patterns,
  then **re-reads the written files and asserts no pattern survives**, and emits
  a SHA-256 manifest. No credential is stored in the repository.
- Live triggers exist and are honest: `q4_green_max_3.py` raises on five
  alternating signs; `continue_half.py` writes `K1_CANDIDATE_HALF.json` on three
  brackets; `run_d1.py` raises `FOUR_ORIGIN_TRIGGER_REQUIRES_REVIEW`;
  `engine.field_record` raises `'STOP: four origin roots require hostile
  reproduction'`. **Three of the four never fired. The fourth fired once and the
  result was rejected in binary128 as a float64 artifact.**

---

## 8. GAPS — what a reader needs that is NOT in the repository

Ordered by how much they would change a decision.

1. **The seven unmerged branches are the state of the art, and nothing in the
   working tree tells you that.** `main`'s `README.md` and `STATUS.md` are
   confidently written and materially out of date: they say no 3+1 precursor was
   found (D1 found two, at exact rationals), they do not mention the engine, the
   night watch, the D2 lane, the determinant reduction, or the resonant
   obstruction. **First action for anyone resuming: fetch and read
   `origin/astra/fastra-afternoon-2026-09-05` and
   `origin/codex/fastra-d1-fold-counts-2026-09-05`.** Merging is not trivial:
   the two D1 branches are disjoint campaigns sharing a directory name, and
   `astra/q4-determinant` conflicts with `astra/resonant-joint` on three files.

2. **There is no chart-completeness argument for any of the search charts, and
   this quietly undermines every "closed numerically" verdict.** The reversible
   chart has a proved local-completeness determinant (`-16(a+2)`, repaired at
   `a=-2`); the Q4 coefficient transport has one (`> 0` for all `kappa`). **The
   Shi chart, the KKL/Kuznetsov chart, and the generalized `(A,m,beta,B,c)` chart
   have none.** `SEARCH_SPACE.md` even concedes for the Kuznetsov chart that "it
   does not represent every affine class" — and then every sweep proceeds anyway.
   F3/F5/F15/F21 report "max 3 in the Shi chart"; that is a statement about a
   coordinate patch, not about quadratic vector fields. **This is the single
   largest unstated gap in the search.**

3. **No validated numerics exist, and by the repository's own acceptance table
   nothing produced counts as evidence.** `STATUS.md` grades evidence
   `NUM / ASYM / THEOREM-FAMILY / CAP-EXIST / CAP-EXACT` and says only the last
   two count toward a counterexample. Every dynamical result in the repository is
   `NUM`. The certification protocol is written down three times and implemented
   zero times. **A resumed campaign should build and benchmark the interval
   verifier on a known control before generating one more candidate** — this is
   recommended in four separate documents and never done.

4. **Eight of twenty engine review items are still live in the shipped code, and
   one of them corrupted the searches that were then reported as closed.** E4's
   broken near-miss gate (`rad > 3e-4*rad[-1]/1e8*1e8`, which reduces to
   `rad > 3*scale`) makes the evolutionary fitness **ignore every radius inside
   3x the nearest equilibrium — exactly where small nested cycles live**. Every
   F4/F14/F20 descent ran with that term. A3/E1 (wrong-direction crossings
   accepted, bisection with no bracket-validity check) and D2 (a descending
   radius grid silently inverts every S/U label) are also unfixed. **The
   "max 4 by descent" conclusion should be re-run with a corrected fitness before
   anyone treats it as informative.**

5. **The Cartesian `count = 2` records were never re-audited under the noise
   rule.** `recount_hits.py` skips anything with `total < 3`, and the two-tolerance
   floor was implemented only in the compactified path. 33 records in
   `F3_lam0_L4.jsonl` and 445 in `W_kklx_L3_worker.jsonl` still carry pre-rule
   counts. The review's own §A1 demonstration was that 3 of 4 such records were
   actually `count = 1`.

6. **Ten queued sweep lanes produced no output and are acknowledged nowhere**
   (`F5_shi_full_L4` at 150k sets, `F5_kklx` at 60k, `F9_mv` at 40k, `F7_q4pert`,
   `F6_q3rpert`, `F9_mvpert`, `F8_yz_evolve`, `F8_shi_evolve`, `F3b_lam0_store1`,
   `F3c_lam0_evolve`). Anyone reading `FABLE_NIGHT_REPORT_2026_09_05.md` would
   believe the search was broader than it was.

7. **Two load-bearing papers could not be read** (§6.2 items 1 and 2). Zhang 1999
   underpins the "3-k ladder" and the K1 question; Artés–Llibre–Schlomiuk 2006 is
   exactly what D2's statement (D) needs. **Institutional access to these two
   papers is probably the highest-value non-computational action available.**

8. **No independent mathematical review of any strike.** Every branch document
   says so in its own words ("has not received an independent mathematical
   review", "no new sub-agent or outside referee audit was performed"). The
   adversarial audits are genuine and caught real errors — including a false
   biconditional in D2's Proposition A, two sign errors in Zhao, an `-K_0`
   omission, a wrong-root selection in the Lyapunov computation, and a six-cycle
   float64 artifact — but every one of them was performed by the same two agents.
   Theorem N in particular is a genuinely nontrivial ODE-comparison proof whose
   only external check was one of those agents.

9. **No persistence box or tube was ever computed** for any four-cycle field, so
   there is no parameter neighbourhood in which four cycles are known to survive.
   Every four-cycle object in the repository is a **single point** with sign
   brackets in the radial variable only. The multipliers are hyperbolic, so a
   neighbourhood exists by the implicit function theorem — but a construction
   that adds a fifth cycle needs an explicit one, and building it is a prerequisite
   for the whole `(4,1)` and `(5,1)` programme.

10. **Provenance for the early numerics is unrecoverable.** KKL rows 1–131 (which
    used the *superseded* projected-derivative engine) and all 150 Shi rows carry
    no evaluator hash. The manifests are post-hoc snapshots and say so. If a
    number from those rows ever becomes load-bearing, it must be recomputed.

11. **A standing-policy contradiction nobody reconciles.** `STATUS.md`'s campaign
    decision is **"YELLOW — no generic coefficient sweep"**, with only the three
    `ATTACK_MATRIX.md` attacks authorized and expansion requiring "a new
    mathematical reduction, a new bifurcation identity, or a verified five-cycle
    candidate". The overnight F3/F5/F9/F15/F21 Sobol sweeps are, on their face,
    generic coefficient sweeps. They may well have been the right thing to do —
    hyperbolic cycles are structurally stable, so `H(2)>=5` would hold on an open
    set and sampling is legitimate, as `FABLE_LANES_2026_09_05.md` argues — but
    the governing document was never amended.

12. **The budget is exhausted and there is no successor plan.** 4096/4096 on the
    shared KKL/Shi ledger; the D1 branch work runs on a separate, unreconciled
    count. Any resumption needs a new budget decision and a decision about which
    ledger it charges.

13. **Things a reader will look for and not find:** any five-cycle candidate; any
    interval certificate; any executable CAPD/Arb code; a merged branch; a single
    current top-level status document that reflects the branches; a computed
    graphic connection or transition coefficient for the KKL infinity organizers;
    a quantified radius for any of the "sufficiently small neighbourhood" results;
    and a proof that the KKL fold component is connected or exhaustively covered.

---

## Bottom line

`H(2) >= 4`. No five-cycle field, no five-cycle candidate, no interval
certificate. The maximum count produced by any method in this repository, across
roughly 4100 charged ODE evaluations plus several hundred thousand swept fields,
is **four** — and the one apparent six was a float64 artifact that the campaign's
own binary128 recheck correctly rejected.

Three things were genuinely accomplished. **Theorem N** closes the Q4
five-interior-zero route by an exact analytic argument, and is the one
publishable closure. **The exact rational 3+1 fields on the D1 branches** are the
first time three origin cycles and a remote cycle have coexisted at one explicit
parameter vector in this family — the object the campaign's own coverage audit
had listed as never obtained, and the gap to a counterexample is now exactly one
cycle in the origin nest. And the repository accumulated a large, honest ledger of
**exact structural obstructions** (the `x=-1` barrier, the Shi two-focus
finite-equilibrium exclusion, the reversible no-finite-saddle theorem, the
`K=J=0` double center, the remote-Hopf sign, the Dulac quartic impossibility, the
D2 corrected Proposition A, the resonant compatibility bound) that together
explain *why* every attempt to stack a degenerate boundary on a degenerate focus
collapsed onto an integrable stratum.

What is missing is not more searching. It is (a) a validated return-map verifier,
(b) a chart-completeness argument that would let a negative sweep mean something,
and (c) two papers behind paywalls.
