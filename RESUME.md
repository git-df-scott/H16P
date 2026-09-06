# RESUME — Lane 2 (cusp manifold / swallow-tail hunt)

Branch `fable/lane2-cusp`. All commands are run from `lane2_cusp_2026_09_06/`
unless stated otherwise.

```bash
git fetch origin && git checkout fable/lane2-cusp && git pull
cd lane2_cusp_2026_09_06
pip install mpmath                                  # not present in a fresh container
g++ -O2 -std=c++17 -DJETDEG=4 cusp_engine.cpp -lquadmath -o cusp_engine4   # engine C (default)
g++ -O2 -std=c++17            cusp_engine.cpp -lquadmath -o cusp_engine    # degree-3 variant
g++ -O2 -std=c++17 engine/cusp128.cpp -lquadmath -o engine/cusp128         # engine A
```

Both compiled binaries are gitignored (`cusp_engine`, `cusp_engine4`,
`cusp_engine_d3`); engine A's `engine/cusp128` is committed. Rebuild before use.

> **Import gotcha.** `lane2_cusp_2026_09_06/engine.py` (my driver, engine C/D) and
> `lane2_cusp_2026_09_06/engine/` (the other session's engine A/B) coexist. Python
> resolves `import engine` to the regular module `engine.py`, which is what the
> Part II code expects. Do **not** add `engine/__init__.py` — that would turn the
> directory into a regular package and shadow the driver. Engine A/B code is imported
> as `engine.eng` / `engine.engB` by its own scripts via explicit paths.

## Smoke test (5 s)

```bash
printf 'V\nD 3 -12 -1.398 8.4 15.28 1.26\nQ\n' | ./cusp_engine4
# -> OK <D> <Dx> <Dxx> <Dxxx> <T> <transversality> <Dxxxx> <nsteps>
# D should be  6.8316004494028633816e-05
```

`D <a> <a20> <a11> <a01> <a10> <x0>` uses the right section (`y=-1, x>1`);
`L ...` uses the left section (`y=-1, x<1`). Numbers are decimal strings parsed at
full working precision. `--ld` selects long double. Env overrides for convergence
testing: `CUSP_ORDER`, `CUSP_TOL`, `CUSP_HMAX`. `CUSP_EXE` picks the binary.

## Re-run the pieces

```bash
# PROTOCOL rule 7 validation (~30 min): Cherkas rows 1-8 both sections,
# row-4 Andronov-Hopf, fold check, rule-1 noise certification
python3 -u validate.py                      # -> validation.json

# enter the cusp manifold at the 8 Cherkas shapes (~5 min)
python3 -u enter_cusp.py                    # -> entry.json

# one cusp curve
python3 -u continue_cusp.py --a 3 --a20 -12 --maxpts 200 --tag row1 --out row1.json

# the (a,a20) grid: regenerate specs, then 4 workers
python3 make_grid_spec.py 4 120 grid
python3 - <<'PY'
import json
gs=[json.load(open("grid_g%d.json"%i)) for i in range(4)]
for i,g in enumerate(gs): json.dump(g, open("run_g%d.json"%i,"w"), indent=1)
PY
for i in 0 1 2 3; do setsid nohup python3 -u campaign.py run_g$i.json camp_r$i.json \
    > log_r$i.txt 2>&1 < /dev/null & done

# summarise every ledger: endpoints, D_xxx sign, sign changes, closest approach
python3 analyse.py 'ledger/cusp_*.jsonl' 'ledger_grid/cusp_*.jsonl'   # -> analysis.json
```

## Where the state lives

| path | what |
|---|---|
| `ledger/cusp_<tag>.jsonl` | append-only, one JSON record per accepted cusp point |
| `ledger_grid/cusp_<tag>.jsonl` | same, for the (a, a20) grid |
| `camp_r*.json`, `log_r*.txt` | per-worker campaign summaries and logs |
| `analysis.json` | `analyse.py` output: the D_xxx sign table |
| `triple_confirm_row1.json` | the rule-1 certified triple cycle at normal amplitude |
| `validation.json` | rule-7 validation output |

Each ledger record carries the exact `(a, a20, a11, a01, a10, x0)` as 34-digit decimal
strings, `D, D_x, D_xx, D_xxx, D_xxxx`, `nu = D_xxx/(D_xxxx r0)`, the weighted residual,
return time, transversality, `V1`, `L = det J(A)`, and Perko's Thm 4.3 nondegeneracy
Jacobians (`perko`).

## The next step, precisely

1. `python3 analyse.py 'ledger*/cusp_*.jsonl'` and read the `sc` (sign-change) column.
   **Any nonzero entry is a swallow-tail bracket** — go to step 3.
2. If all zero: compare `sign_start` and `sign_end` in `analysis.json` across the grid.
   Two grid neighbours with different `sign_end` bracket a swallow-tail in `(a, a20)`;
   bisect between them.
3. Refine a swallow-tail with the square Newton (4 equations, 4 unknowns; `a20` is the
   fourth unfolding parameter, `a` and `x0` held):

```python
import mpmath as mp
from engine import Engine
from swallow import try_from_cusp_point
mp.mp.dps = 50
e = Engine(quad=True)
out = try_from_cusp_point(e, a, a20, [a11, a01, a10], x0, verbose=True)
# out['perko43'] carries every Perko Thm 4.3 nondegeneracy quantity
```

   Do **not** seed this from a small-amplitude cusp point: `D_xxx = 48 d7 r0^4` there
   and `d7 = 0` only on the centre variety where `D` vanishes identically, so the Newton
   provably runs `a20 -> -inf`. Seed from `x0` of normal size. (Recorded in
   REPORT_lane2.md §II.3.)

4. If a swallow-tail is found and Perko's conditions hold, unfold and count:

```python
from probe import Noise, sign_changes, triple_confirm
```

   `probe.unfold_cusp(c, mu, x0, target)` sets `(D, D_x, D_xx)` to a prescribed target
   through the 3x3 parameter Jacobian; for a swallow-tail the analogous 4-parameter
   solve (including `a20`) puts the field in the swallow-tail region. Count sign changes
   of `D` with `probe.sign_changes(..., noise=Noise())`, which applies PROTOCOL rule 1's
   two-tolerance certification.

5. **Four certified sign changes in one nest is a TRIGGER (PROTOCOL rule 3):** stop,
   recheck every bracket endpoint in binary128 *and* with engine D (`indep_engine.py`)
   and engine A (`engine/cusp128`), write
   `TRIGGER_lane2_<UTC timestamp>.json` with the exact rational coefficient vector,
   section, brackets and all engines' values, push immediately, and **do not announce a
   counterexample** — the auditor decides.

## Budget notes

* One continuation point costs ~1.7 s wall (binary128, ~9 engine calls with the chord
  Newton). 120 points ≈ 3.5 min. Four workers on 4 cores.
* Never take `D_xxx` by finite differencing, and never in double: the jet gives it
  exactly. `long double` (`--ld`) is 20x faster and fine for reconnaissance, but every
  number that decides anything is binary128.
* Do not spend time on Q4, the reversible reseed, or KKL fold continuation at radius
  > 1e10 (PROTOCOL rule 9).
