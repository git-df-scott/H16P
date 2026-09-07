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
been made to work yet is the one PROTOCOL section (c) names: **start from the
cusp, not from the seed.**  `bautin.py` and `bautinrun.py` build that, and they
are one calibration away from being correct.

### What is already verified

`bautin.focal_values` is a faithful transcription of Cherkas-Artes-Llibre
eq. (15): on the paper's own third-order weak focus family it returns
`V1 = V3 = V5 = 0`, `V7 != 0`, and the engine confirms the focus really is
third-order -- in binary128 with the two-tolerance gate, `log|D|` against
`log s` has slope **7.10** and `D/s^7` is constant to four figures over a
decade in `s`.  The three cusp conditions are solved exactly (residuals
`1e-20` to `1e-13`), the unfolding in `(e1, e3)` is closed-form, and the
`(a, a20, r0)` sweep runs.

### What was wrong, and the fix -- now built in `bautin2.py`

The printed `V1, V3, V5, V7` are focal values **up to per-field normalisation
constants**.  Measured against the true leading Taylor coefficient of `D`:

    c7 / V7 = -6.7e-05   (a = -2,  a20 = -1)
    c7 / V7 = -9.0e-06   (a =  3,  a20 = -12)

Not constant across fields, and negative -- the sign convention differs too.
So `V5 : V3 : V1 = -3 r0^2 : 3 r0^4 : -r0^6` imposed on the *printed* values is
not the triple-root condition on the true series, which is why the 150 fields
in `ledger/bautin_*.jsonl` have no triple cycle: `D` has no sign change near
the focus on any of them.

`bautin2.py` measures the true coefficients instead, and no normalisation
constant survives that.  It is written and validated; what remains is to run
the continuation with it.

    bautin2.taylor_D(L, phi, s_lo, s_hi)      # c1, c3, c5, c7 of D(s)
    bautin2.fit_window(L, phi, s_ref)         # picks the window, ~8 s
    bautin2.third_order_seed(a, a20)          # the exact r0 = 0 point
    bautin2.newton(a, a20, x0, r0)            # the three calibrated conditions

Validated on the third-order family, where the direct measurement gives
`c7 = 4.0758` and `c5 = 0` exactly:

| a, a20 | window in s | rel resid | c1 | c3 | c5 | c7 |
|---|---|---|---|---|---|---|
| -2, -1 | [4.8e-3, 1.0e-2] | 5.3e-07 | -7.7e-17 | 1.3e-11 | -9.3e-07 | 4.09999 |
| 3, -12 | [5.5e-3, 1.4e-2] | 1.3e-05 | -5.5e-16 | 6.2e-11 | -3.1e-06 | -1.57523 |
| 1.5, -15 | [1.2e-2, 2.7e-2] | 7.4e-06 | 4.8e-16 | -1.3e-11 | 1.5e-07 | 0.01526 |

`c7` to 0.6% and the spurious `c5` down to 1-7% of the `c7` term.  Three
things had to be right at once, and each of them silently ruined the fit on
its own:

1. **The two-tolerance gate.**  `|D|` runs down to `1e-28` while the binary128
   floor is around `1e-24`.  Ungated, the first pass reported the weak focus as
   second-order (slope 5.03 instead of 7.10).
2. **Relative weighting.**  `D/s` spans seven decades across a fit window, so
   an unweighted least squares is decided entirely by the top of the window --
   and the small-`s` points are exactly the ones that pin `c5`.  Unweighted,
   `c5` came out `-1.2e-04` where the truth is 0.
3. **Fitting the tail.**  A bare cubic in `w = s^2` makes `c5` absorb the
   `O(s^9)` term.  Degree 5 with a narrow window: `c5` improves by a factor of
   40 and `c7` from 9% off to 0.6% off.

Cost is 8 s per fit at rtol 1e-18 / 1e-16.  Newton needs four fits per
iteration, so a solve is a few minutes; a `(a, a20)` grid crossed with an `r0`
ladder is an overnight run.  If that is too slow, the lever is the integrator,
not the tolerance: the binary128 build is Dormand-Prince 5(4), so `rtol 1e-20`
costs 4.4 s per return against 0.12 s at 1e-18.  An order-8 scheme in the quad
build would cut that by about two orders of magnitude.

### Then the actual question

Three small cycles out of the focus is Bautin's cap and is not new.  What is
new is whether a field **on** the cusp manifold also carries an extremum of
`beta*` further out: an extremum beyond the cusp is a third extremum, and a
third extremum is four cycles in the nest.  Bautin caps what comes out of the
focus, not what the rest of the nest does.  `bautinrun.py` already records the
full `beta*` for every constructed field, so once the calibration is in, that
question is answered by re-reading the same ledger column.

Watch also for the swallowtail `D = D_s = D_ss = D_sss = 0`: by Perko 1995
Theorem 4.3 a multiplicity-four cycle with the nondegeneracy conditions
*forces* four simple cycles in the nest nearby -- it does not merely permit
them.  With the calibration in place that is `c1, c3, c5, c7` all tied to a
fourth-order root, one more condition than above and one more free parameter
(`a20` joins `a11`, `a10`, `a01`).

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
