# Exact centre certificate, solver contract, and the separated double-fold strike

**Date 2026-09-07. Branch `opus/rank-repair-cusp-compatibility-2026-09-06`,
continuing from `5a29d0f`. NO FIVE-CYCLE FIELD.**

## 1. Exact centre certificate (verified)

The family supplied in review, at `a = -2`:

\[ a_{11}=k,\ a_{01}=-k-3,\ a_{20}=-\tfrac{2k^2}{27},\
   a_{10}=\tfrac{5k(k+3)}{27},\ a_{00}=-\tfrac{k^2+5k+9}{9} \]

is consistent with the pinning `a00 = a01+a11-a10-a20-a` (checked), and with

\[ U=y+\tfrac{k+3}{6}-\tfrac{kx}{9},\quad
   V=-1+xy+\tfrac{(k+3)x}{3}-\tfrac{2kx^2}{9},\quad
   C=\tfrac{kx}{3}-y-\tfrac{k+3}{3} \]

**every identity holds exactly** (sympy, difference identically zero):

| claim | result |
|---|---|
| `X(U) - 2CU` | `0` |
| `X(V) - CV` | `0` |
| `X(U/V^2)` | `0` — `H=U/V^2` is a first integral |
| `V(A)` | `(k-9)/9` |
| `trace J(A)` | `0` |
| `det J(A)` | `(k-3)(9-k)/27` |
| `grad H(A)` | `(0,0)` |
| `H_yy(A)` | `729/(9-k)^3` |
| `det Hess H(A)` | `19683(k-3)/(9-k)^5` |

so `3<k<9` gives a nondegenerate centre. **Attribution: this family and its
Darboux factors were supplied in review, not found here; this section is a
verification.**

### It is exactly the Newton limit

With `k=\sqrt{1627}/5`:

| coeff | exact from `k` | Newton limit (`5a29d0f`) | difference |
|---|---|---|---|
| `a11` | `8.067217612039481856486` | `8.067217612039481856440` | `4.6e-20` |
| `a01` | `-11.06721761203948185649` | `-11.06721761203948185644` | `-4.6e-20` |
| `a20` | `-3254/675` exactly | `-3254/675` | `0` |
| `a10` | `1627/135+\sqrt{1627}/9` | `16.53363941409600843877` | `1.4e-20` |

Agreement at the Newton's own convergence level, and `k=8.0672\in(3,9)`, so the
numerical limit **is** this nondegenerate centre. It is used below as a
regression case and as an identified unwanted solution set — not re-discovered
at neighbouring `a_{20}`.

## 2. Solver success contract (fixed)

`solver.py` returns `CONVERGED / STALLED / SINGULAR / RETURN_FAILED /
GUARD_EXIT / ITERATION_LIMIT`, the final residual vector **componentwise**, its
documented scaling (component `i` of `(D,D_x,D_{xx},D_{xxx})` weighted by
`r_0^i/i!`, the Taylor term it controls), and the best diagnostic iterate
**separately** from the root. `Result.root` is `None` unless `CONVERGED`, so no
consumer can read a parameter vector as convergence. Every component must pass,
not just a norm.

Regression cases (`solver_regression.py`):

* the exact centre — returns `STALLED` at the noise floor
  (`2e-24 … 5.7e-24`), as it must: the system is degenerate there
  (`det Jac = det(F_u)\cdot D_{ssss}` with `D_{ssss}=0`), so Newton stagnates
  rather than converging. **This also corrects `5a29d0f`, which described that
  run as "converged with residual 5.7e-24"; it was stalling at the floor.**
* the previously mis-reported offsets — `STALLED`, `root is None`, worst
  components exposed as `1.5e-3` and `4.8e-4`.

## 3. The separated double-fold system

\[ F(\mu,s_1,s_2)=(D(s_1),D_s(s_1),D(s_2),D_s(s_2))=0 \]

square in `(p_1,p_2,s_1,s_2)` for a chosen control pair; `s`-columns free from
the jet, parameter columns by central differences. Guards:

* **separation** `|s_1-s_2| \ge 0.05` — `s_1\to s_2` is a cusp, not two folds;
* **amplitude** `\min|s_i-1| \ge 0.08`, with damped steps. This is essential:
  the system has a **trivial branch at the focus**, since `D(1)=0` identically
  and `D_s(1)=V_1`, so `s\to1` with `V_1\to0` solves it. That is the
  Bautin/Hopf boundary (PROTOCOL mechanism (a)), not a separated double fold.

## 4. Entry census — the promising part

Displacement census on the `(3,1)` controls with a second focus, and on row 2:

| row | `D` zeros | `D_s` zeros | stationary structure |
|---|---|---|---|
| 7 `(3,1)` | 3 | **3** | max `D=+1.5e-5`, min `-3.5e-5`, max `+7.4e-4` |
| 8 `(3,1)` | 3 | **3** | min `-6.8e-7`, max `+5.0e-6`, min `-4.0e-5` |
| 2 `(3,0)` | 3 | **3** | min `-3.4e-6`, max `+7.6e-5`, min `-2.7e-4` |

**Three** interior stationary points, not two — and the outer two are the
**same type with the same sign**. Driving both to zero is exactly two separated
double folds, unfolding to four cycles in the nest. So the target is not
excluded by the seed structure.

## 5. What actually happened

**21 solves** (3 rows x 7 control pairs, including the strong shape controls
`a` and `a20` chosen by conditioning, not frozen): **every one drives `s_1` onto
the focus.** Without the amplitude guard `s_1\to1` with the residual falling to
`~1e-9`; with it, `s_1` pins at the guard `1.08000` and stalls at
`~1e-6`–`1e-5`. No slice reaches a finite-amplitude solution.

**Why — the unfolding rank.** `B_{full}=[\partial_\mu D(s_1);\partial_\mu D(s_2)]`
over `\mu=(a,a_{20},a_{11},a_{01},a_{10})` has **rank 2**, but badly conditioned:

| row | `cos` angle between rows | `|row_1|`,`|row_2|` | `\sigma_{min}/\sigma_{max}` |
|---|---|---|---|
| 8 | `0.9693` | `0.100`, `2.233` | `0.125` |
| 7 | `0.9298` | `0.355`, `5.006` | `0.191` |
| 2 | `0.9892` | `0.231`, `7.859` | `0.0737` |

The rows are nearly parallel and the outer one is `14`–`34x` longer: every
control moves the outer fold far more than the inner, in almost the same
direction. Rank is not lost — the path is.

**The structural obstruction.** Continuing the *single* outer fold instead
(instructed step 4) converges cleanly — genuine finite-amplitude double cycles:

| row | fold `s_2` | residual | `D_{ss}(s_2)` | other stationary point |
|---|---|---|---|---|
| 7 | `3.46899418303` | `1.2e-33` | `-1.290e-3` (max) | `s_1=2.2861`, `D=-2.50e-4`, `D_{ss}=+7.73e-4` (min) |
| 8 | `3.2755755769` | `8.1e-29` | `+9.776e-5` (min) | `s_1=2.1770`, `D=+1.66e-5`, `D_{ss}=-6.10e-5` (max) |

On the fold's constraint surface the surviving neighbour is of the **opposite
type**, and the inner **same-type** extremum present at the seed (row 7: the max
near `1.28`) is **gone**. Two opposite-type folds cannot both sit at zero — `D`
would have to vanish identically between a maximum at `0` and a minimum at `0`,
i.e. be a centre.

> **Failure mode, named precisely: the inner same-type extremum is destroyed
> while the outer fold is driven to zero.** Not a coalescence of the two folds,
> not a rank loss, not a return failure, and not a centre approach along the
> guarded path.

## 6. Scope

This is a **finite-path result, not a route exclusion**. It covers three seeds,
seven control pairs, one amplitude guard value, and a single-slice Newton — not
pseudo-arclength continuation on the full constraint system, not other `(3,1)`
controls or their saved descendants, and not shapes reached by first moving
`(a,a_{20})` far from these rows. The guards define the search scope; they are
not mathematical exclusions. No four-cycle unfolding was constructed, so
sections 5 of the brief (positive construction, remote cycle, certification)
were not reached.

## 7. Replay

```bash
python3 centre_2026_09_07/exact_centre.py       # all Darboux/centre identities
python3 centre_2026_09_07/match_newton.py       # k = sqrt(1627)/5 vs the Newton limit
cd lane2_cusp_2026_09_06
python3 regression_tests.py                     # Bautin baseline + coordinate rank
python3 solver_regression.py                    # solver contract cases
python3 extrema_census.py                       # three stationary points at rows 2,7,8
python3 run_doublefold.py                       # the 21 double-fold solves
python3 rank_Bfull.py                           # unfolding rank and conditioning
python3 singlefold_track.py                     # single fold + the surviving neighbour
```

New evaluations logged under `lane2_cusp_2026_09_06/ledger_opus/`.
