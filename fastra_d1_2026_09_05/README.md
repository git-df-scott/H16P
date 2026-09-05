# FASTRA D1 reproduction

The authoritative report is `../FASTRA_D1_REPORT_2026_09_05.md`.
This directory contains a deterministic continuation and its full sampled data,
not a new random parameter sweep. Each vector is a list of exact fraction strings
in the six-monomial coefficient order for P followed by Q.

Requirements: Python, NumPy, SciPy, mpmath; GCC/G++ with OpenMP and libquadmath.
Run from the repository root:

```bash
gcc -O3 -shared -fPIC -fopenmp audit/fable_engine/retmap.c -lm -o audit/fable_engine/libretmap.so
gcc -O3 -shared -fPIC -fopenmp audit/fable_engine/retmap_log.c -lm -o audit/fable_engine/libretmap_log.so
g++ -O2 -std=c++17 -fext-numeric-literals -shared -fPIC fastra_d1_2026_09_05/matching_quad.cpp -lquadmath -o fastra_d1_2026_09_05/.matching_quad.so
g++ -O2 -std=c++17 -fext-numeric-literals -shared -fPIC fastra_d1_2026_09_05/matching_quad.cpp -lquadmath -o fastra_d1_2026_09_05/.matching_quad_v2.so
g++ -O2 -std=c++17 -fext-numeric-literals -shared -fPIC fastra_d1_2026_09_05/matching_quad.cpp -lquadmath -o fastra_d1_2026_09_05/.matching_quad_v3.so
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 python fastra_d1_2026_09_05/verify_precursor.py
python fastra_d1_2026_09_05/summarize.py
python fastra_d1_2026_09_05/write_report.py
```

`verify_precursor.py` re-integrates only the rational 3+1 control. It appends a
new `rational_precursor` field record; summary generation uses the latest record
per label, so duplicate retries are preserved without double counting fields.
The displayed radii are section-dependent, numerical values.

`run_sheets.py positive_center`, `positive_infinity`, and `negative` specify the
original target grids. `resume.py` records the recovery paths after a rejected
large-radius raw trigger and an overlarge negative-K predictor. These scripts
append records; do not execute them to regenerate only the report. They are
research drivers, not certified interval continuation software.

`repair_counts.py center` rechecks near-center fields not already repaired;
`repair_counts.py beta` completes the accepted large-positive folds using a beta
chosen from their bounded matching derivative. `check_sign_disagreements.py`
rechecks unresolved diagnostic disagreements named by the current summary.

`matching_quad.cpp` is an independent two-state adaptation of the archived
modified-midpoint binary128 shooter. It also includes a remote-nest entry point
which re-solves the remote equilibrium before translation. The v2 shared-library
filename is retained for the precursor replay. The build commands use the
final consolidated source. `matching_quad_initial.cpp` preserves the first
origin-only implementation; `matching_quad_remote.cpp` preserves the version
used for the exact 3+1 reproduction. The final source retains the requested
position tolerance but floors the sensitivity tolerance at 1e-24; this repairs
the late step-guard failure without discarding position accuracy.
`precision_controls.jsonl` records the corresponding tolerance controls.
`finish_positive.py` uses that repaired controller and polishes the final
fold locations in the actual-m chart before pair/Hopf profiling. Shared libraries are rebuildable and are not tracked.
`engine.half` builds the unchanged archived augmented shooter on demand.

Statuses:

- Baseline `redge`: last sampled successful full return, not a proven boundary.
- Full-profile `scan_cap`: the configured scan ended with a successful return.
- Full-profile `integration_failure`: original engine status, with the last
  successful/first failed radius; status 5 means a nonfinite returned value.
- Matching-profile `angular_chart_failure`: a half passage or numerical guard
  failed. It is unresolved, not absence of cycles.
- Matching terminal signs outside established full-return domains are an
  auxiliary extension of the matching residual, not displacement values.
- Near-zero sign changes are saved separately. Every numerical count is a
  lower bound on detected roots on a finite grid, not a global upper bound.
- `FOUR_ORIGIN_TRIGGER.json` is a **rejected raw false positive**; read
  `trigger_resolution.json` before using it.

Decimal coefficient strings produced by augmented shooting are exact rationals
when passed through `Fraction`. The actual field vectors, rather than rounded
human-readable (c,K,beta) displays, determine replay. Nominal-fold coefficient
vectors need not lie exactly on the discriminant and may have a tiny split pair.

The last log-radius-42 extension is in `finish_positive_extension.py`: after the
archived augmented shooter hit its per-call time fuse, the two-state binary128
F=G equations were solved directly with finite-difference parameter derivatives.
Its sampled curvature chooses the pair side. The accepted record includes the
polish history and paired tolerance controls. `make_manifest.py` regenerates
hashes after report/audit generation; it performs no integrations.
