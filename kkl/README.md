# KKL bounded construction tools

Read [the branch checkpoint](../STRIKE5_PRECURSOR.md) before interpreting
the data. This directory contains numerical discovery tools, exact algebra
and finite-equilibrium gates. It does not contain interval-certified
periodic orbits or a five-cycle field.

Use Python with NumPy, SciPy and SymPy. Recorded versions are in
`data/strike_summary.json`. The source is macOS/Linux oriented and uses
the `resource` module for each evaluator's ten-CPU-second fuse. The
supervisor is serial and also imposes a 20-second wall timeout per call.

Exact replay, without orbit integrations:

```sh
python kkl/check_exact.py
```

A charged replay of the initial origin control:

```sh
python kkl/pilot.py --request '{"r":64.55543434,"c":"7/10","alpha":"-80","beta":"0"}' --purpose independent_start_replay
```

This appends one call to the persistent ledger. All controls, derivative
checks and failed calls count toward the 4096 ceiling. `pilot.py` is the
low-level charged supervisor; the pilot's 64-call/16-step ceiling was
enforced by `continue_path.py`. The completed pilot is preserved in
`data/pilot_summary.json`.

`continue_path.py` follows the two original roots on K=6/5. It resumes the
last accepted point in its existing file; it does not restart a fresh
experiment. `follow_segment.py` follows an explicitly supplied segment
in (c,K), starting from an accepted common-parameter root pair. Both stop
at a geometry, radius, numerical-corrector or budget boundary. An ordinary
Newton corrector is deliberately stopped near a fold; it does not pretend
to continue through one without pseudo-arclength. No fold was located in
this checkpoint, so an augmented fold-curve solver has not been run.

The 15-state fixed-time second-derivative system is available for origin
fold diagnostics, with three transverse-determinant states added for
stable first derivatives. The default first-order system has twelve
scalar states. Large remote second derivatives have not been validated.

The return-coordinate ranges are `[2^-12,2^10]` for the origin and
`[-2^20,-1]` for the remote section, with downward crossing required.
The coordinate guard `max(|x|,|y|)<1e7` and finite time horizon 10 are
numerical safeguards; a stopped integration is unresolved, not an escape
proof. The low-level evaluator can run explicitly labelled controls
outside the precursor K margin; continuation always applies exact gates.

`data/returns.jsonl` is the append-only raw call ledger. Individual path
files preserve the accepted points and unsuccessful attempts.
`data/continuation_events.jsonl` aggregates those records, preserving
source filenames and line numbers. `summarize.py` updates this aggregate,
totals and the final source/data hash manifest without running ODEs:

```sh
python kkl/summarize.py
```

No whole-box exclusion follows from these files. No beta<0 completion is
authorized by merely finding the two currently tracked cycles; it needs
the complete four-cycle precursor at one common beta-zero field first.
