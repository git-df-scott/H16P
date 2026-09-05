# Campaign review package

Read [the synthesis](../FASTRA_ZOOM_OUT_2026_09_05.md) first.

| File | Purpose |
|---|---|
| [COVERAGE_AUDIT.md](COVERAGE_AUDIT.md) | What the 206 calls, drivers and derivatives actually cover |
| [KKL_NEXT_CONSTRUCTION.md](KKL_NEXT_CONSTRUCTION.md) | Stationary-return pilot and explicit finite-beta alternatives |
| [FRONTIER_AUDIT.md](FRONTIER_AUDIT.md) | Scope of route closures, Q4 endpoint gaps and other four-cycle seeds |
| [KKL_SECTION_REPAIR.md](KKL_SECTION_REPAIR.md) | Proof of a section meeting every cycle in either KKL nest under stated gates |
| [check_nullcline_section.py](check_nullcline_section.py) | Ten-CPU-second exact identity replay, no ODEs |
| [nullcline_exact_check.txt](nullcline_exact_check.txt) | Saved replay result |
| [EVIDENCE_MANIFEST.json](EVIDENCE_MANIFEST.json) | Baseline commit, unchanged ledger hash and evidence counts |
| [export_visible_transcript.py](export_visible_transcript.py) | Exports only visible events and separate user objectives; redacts credentials |

The full conversation Markdown, JSONL and transcript manifest are delivered
locally in `outputs/review_2026_09_05/`; they are not published in this
research commit. They contain 78 public messages through the zoom-out
request and four separate recorded objectives. Older campaign phases are
reconstructed from commits in the synthesis, not presented as verbatim
conversations absent from this task's record.

Baseline main is `6048ed8b09e12094b0875fdc68c5a1fb6341775d`. The review
adds no ODE evaluations, changes no previous KKL data, and modifies no
Claude/Fable file. The charged evaluation count remains 206 of 4096.
