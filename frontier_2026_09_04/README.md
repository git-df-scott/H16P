# Frontier evidence and replay

These files accompany the September4,2026 audit. All new computational results are **NUMERICAL ONLY**, except explicitly stated rational symbolic identities. No interval ODE remainder is certified here.

## Numerical replay

Create an isolated Python environment, install `requirements.txt`, and run from this directory:

```sh
python seed_numerics.py shi_visual
python seed_numerics.py chen_visual
python seed_numerics.py kkl 3e-14
python seed_numerics.py yu_zeng 3e-14
python seed_numerics.py gt_remote
python gt_taylor.py 112
python gt_taylor.py 128
python reproduce_shi_original.py
python reproduce_gt_continuation.py
python hopf_probe.py
python infinity_probe.py
```

The scripts write to `data/`. The first five find return roots, equilibria, and projective infinity directions; cycle multipliers are numerical. Scan arrays contain section coordinate, displacement, period and divergence exponential (or errors); an exponential away from a fixed point is **not** P'. Finite horizon, sampled brackets and failed integrations do not certify nonexistence. A sign change across a pole or changed itinerary must be rejected. Saved `hopf_start_remote.json` supplies a separate bounded remote root control; the broader remote probe records discontinuities and is not a cycle count.

The MPFR tiny-cycle script scales x,y by the starting radius, constructs a polynomial Taylor step for the quadratic vector field, integrates to the last fraction of a turn and solves x=0 locally. Repeated Taylor orders agree for all printed signs, but there is no validated truncation enclosure. It is not CAPD and must not be called a new computer-assisted proof.

`gt_continuation.json` records exact-path sample checks; `gt_relaxed.json` records the remote return at s=.01. The continuation does not certify all intermediate parameter values. `infinity_probe.json` uses a truncated asymptotic branch initialization; it is retained as an unsuccessful control, superseded for the vertical-connection hypothesis by the exact no-contact obstruction.

The two-tolerance seed comparisons were inspected during the audit. The saved KKL and Yu–Zeng data are the tighter run. Yu–Zeng's weakest root changes at about2e-8 between runs, consistent with weak multiplier sensitivity; reported ledger coordinates are rounded. Do not treat a printed zero residual as exact arithmetic.

## Mathematical/source appendices

- `STATUS_SOURCE_AUDIT.md`: current status, citation chains, recent claims and coverage gaps.
- `SHI_TOPOLOGY_AUDIT.md`: exact focus/equilibrium/infinity calculations and theorem scope.
- `Q4_GRAPHICS_AUDIT.md`: original Q4 endpoint transport, higher Melnikov and graphics ledger.
- `STRIKE_REDTEAM.MD`: independent hostile review of the selected first-order Hopf precursor.
- `GT_SIGN_CAVEAT.MD`: independent check of the published Lemma2 displacement label.

Original source papers are linked in these files; downloaded papers and rendered inspection images were kept as working material, not republished in this repository. Publication/version checks used primary author, journal and arXiv records. Access failures are disclosed.

The original Shi lambda=−10^-250 instance was also replayed: its innermost bracket is around7e-100, with1200-bit orders128/144 agreeing to35 printed significant digits; the middle/outer brackets and remote return were separately recomputed. See `shi_original_taylor.json` and `shi_original_remote.json`.

## Integrity check

From the repository root, run `shasum -a 256 -c frontier_2026_09_04/SHA256SUMS` to check the saved audit artifacts. This verifies file integrity, not numerical rigor or mathematical correctness. Before this standalone packaging, the original audit commit `70b8a25` compared all 68 preexisting Q4 files with source commit `5bcfe11`; every byte was preserved. Those historical files remain in the source history and are not included in this audit-only branch.
