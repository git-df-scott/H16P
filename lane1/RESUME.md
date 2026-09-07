# RESUME -- Lane 1 (`lane1/`), Andronov-Hopf curve sweep

Session `session_01UqdMSqDz9KHbgPtDkycpJq`, branch `fable/lane1-ahcurve`.
A second Lane 1 session works in `lane1_ahcurve/`; the trees do not overlap.
**Nothing has triggered.  No counterexample is claimed.**

## Setup from a cold container

```sh
git fetch origin && git checkout fable/lane1-ahcurve
pip install numpy scipy mpmath matplotlib
cd lane1
python3 -c "import engine"        # builds build/libretmap1*_<hash>.so (double + binary128)
```

`engine.py` compiles `retmap1.c` twice on import (`cc -O3 -march=native
-fno-fast-math -fopenmp`, and again with `-DLANE1_QUAD -lquadmath`) into
`lane1/build/`, stamped with the source hash, so an old binary is never reused
for a changed source and a ledger row's `engine_sha256` always names the binary
that produced it.

## Re-run the validation (PROTOCOL rule 7) -- do this first, always

```sh
cd lane1 && python3 validate.py          # ~80 s, writes ../VALIDATION.md
```

Must reproduce: three sign changes of `D(s,0)` on all 8 Cherkas rows; the KKL
control at `r = 0.6832, 2.1837, 15.9628` to four decimals plus its remote cycle
at `s = 3706` from `B = (-6.259641, 7.449768)`; exactly two interior extrema of
`beta*` on every row; row 4's `beta*` extrema at `x = 0.6207, 0.7981` against
the published polynomial's `0.6238, 0.8056`.  If any of that moves, stop and
fix the engine before sweeping.

## Re-run the campaigns

```sh
cd lane1
# random perturbations + (1+1)-ES  (the ES stalls -- see below)
OMP_NUM_THREADS=1 ./run_lane1.sh cherkas1,cherkas2 p1     # and p2/p3/p4 splits

# one-dimensional scans of the extremum count along random lines  (the good one)
OMP_NUM_THREADS=1 python3 linescan.py --seed cherkas1,cherkas2,cherkas3 \
    --tag s1 --dirs 200 --n 200 --nlam 16 --rngseed 41001
OMP_NUM_THREADS=1 python3 linescan.py --seed cherkas4,cherkas5,cherkas6 \
    --tag s2 --dirs 200 --n 200 --nlam 16 --rngseed 41002
OMP_NUM_THREADS=1 python3 linescan.py --seed cherkas7,cherkas8,perkoP3 \
    --tag s3 --dirs 200 --n 200 --nlam 16 --rngseed 41003

# Newton continuation onto the multiplicity-three cycle
OMP_NUM_THREADS=1 python3 cusprun.py --seed cherkas1,cherkas2,cherkas3,cherkas4,\
cherkas5,cherkas6,cherkas7,cherkas8,perkoP3 --tag c1 --dirs 30 --n 200 --rngseed 99001

python3 analyze.py                       # ledger summary
```

Four single-threaded processes beat one 4-thread process: the OpenMP loop is
over section points inside a single field, and independent fields parallelise
with no synchronisation at all.

`Tmax = 100` in `linescan.SWEEP` is a speed choice, not an accuracy one.  The
longest return over every seed's whole AH domain is `6.58` time units, and the
status arrays at `Tmax = 400` and `Tmax = 100` are bitwise identical on all
nine seeds (the check is three lines; re-run it if you change the seed set).
Leaving it at 400 makes a failed return cost 4x more and slowed the first
launch of the line scan by roughly 30x.

## What is known

* No field with three interior extrema of `beta*` has been produced by any
  campaign.  Everything the sweep has found has 2, 1, or 0.
* The three-cycle region is of relative width `~1e-3` in the 8 live
  coefficient directions: 1e-3 perturbations keep it 26% of the time, 1e-2
  keep it 1.4%, 1e-1 never.
* The cusp that *annihilates* the seed's two extrema is at relative
  displacement `~5e-4`; no cusp that would *create* a pair was found in the
  outer annulus out to `0.8`.
* Every seed with a second focus has a **monotone** remote `beta*` -- one
  cycle there at any rotation angle, as Zegeling's Lemma 6.4 predicts.

## Known defects, in the order they should be fixed

1. **The (1+1)-ES in `sweep.py` is broken as a search.**  All 18 runs collapsed
   the step size to its `1e-6` floor by roughly iteration 100 and proposed
   nothing afterwards (acceptance 0-31 of 500).  Do not read anything into
   those ledger rows beyond the extremum counts.  Either delete the ES or
   replace it: restart on collapse, floor `sigma` at `1e-3`, and rank on the
   discrete count first and the continuous surrogate only within a tie.
2. **The `outer` window of the cusp solve almost never converges** (22 of 24
   starts on row 4).  The residual there is large -- `beta*` plunges at the
   nest boundary -- so a plain Newton with a `+/-0.05` trust region in
   `lambda` cannot get anywhere.  It needs pseudo-arclength continuation in
   `lambda` with the fold branch as the predictor, not a cold Newton.
3. **`beta*` is UNRESOLVED near the outer end for cherkas1 (183/322 resolved)
   and cherkas6 (228/322).**  The rotated member loses its return before `D`
   changes sign, so the interval on which a third extremum could be seen is
   truncated exactly where Perko's route (b) would put one.  This needs a
   compactified chart (`retmap_log.c` on `astra/fastra-afternoon-2026-09-05`
   is the model), not more precision.
4. `Rmax = 1e4` and the `bmax = 1.5` cap on the rotation search are inherited
   defaults, not measured ones.  Neither has been shown to bind, but neither
   has been checked.

## The next step, concretely

The line scans and the cusp solve agree that a third extremum is not reachable
by moving a *published seed* along a random line.  The instrument that has not
been tried is the one PROTOCOL section (c) names: **start from the cusp, not
from the seed.**

Triple cycles are known to exist at small amplitude in the Bautin unfolding of
a third-order weak focus,
`D(r) ~ V1 r + V3 r^3 + V5 r^5 + V7 r^7` with
`V5 = -3 r0^2 V7`, `V3 = 3 r0^4 V7`, `V1 = -r0^6 V7` giving a triple cycle at
`r0`.  Cherkas--Artes--Llibre's own two-parameter family of third-order weak
foci is in `SEEDS.json`
(`a11 = 4 - 2a`, `a01 = 2a + 1 - a11`,
`a10 = (6(a^2 - a - 2) + a20(6a - 7))/(1 - 3a)`, parameters `a` and `a20`),
and `V3, V5, V7` are given in closed form in `LIT_A` section 3.1.  So:

1. Solve the three focal-value conditions for a triple cycle at small `r0` on
   that two-parameter family -- this is algebra, not a search, and it lands
   exactly on the cusp manifold.
2. Continue that cusp in `r0` up to normal amplitude with `cusp.py`'s Newton,
   using the small-`r0` solution as the first predictor.  **Nobody has done
   this** (PROTOCOL section (c), last sentence).
3. Watch for the swallowtail `D = D_s = D_ss = D_sss = 0`.  By Perko 1995
   Theorem 4.3 a multiplicity-four cycle with the nondegeneracy conditions
   *forces* four simple cycles in the nest nearby -- it does not merely permit
   them.  That is the only route in this literature that converts a solved
   equation directly into four cycles.

The apparatus for step 3 already exists here: extend `cusp.py`'s stencil fit
from degree 4 to degree 6 and solve `A_u = A_uu = A_uuu = 0` in
`(u, lambda1, lambda2)`.

## Trigger discipline (PROTOCOL rule 3) -- wired, never fired

`sweep.check_trigger` fires on `n_extrema >= 3`.  It writes
`TRIGGER_lane1_<UTC timestamp>.json` at the repository root containing the
exact `local10` coefficients (both `repr` decimal and `float.hex`, each
round-tripping the double exactly), the section angle, the overlap window of
the three extrema, the sign-change brackets of `D(.,b)` at the mid-window
rotation, the same count on a **second section** rotated by `pi/2`, the SciPy
engine's `D` at every bracket endpoint, and `hiprec.recheck`'s binary128 values
at two tolerances with their noise estimates.  Then it commits and pushes.
**Do not announce a counterexample; the auditor decides.**
