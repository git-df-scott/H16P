# Overnight H16P investigation — 7–8 September 2026

**No counterexample found.** The target was at least five distinct isolated periodic orbits of one real planar quadratic polynomial field. Numerical brackets are not certified cycles; failed returns do not prove absence.

Start with the [research status](outputs/H16P-current-status.md), [scope corrections](outputs/H16P-focus-route-correction.md), and [latest extreme boundary search](outputs/H16P-extreme-boundary-search.md). The [checkpoint](work/h16p-audit/CHECKPOINT.md) preserves the research sequence, including failed approaches and subsequent corrections. Later correction reports take precedence over earlier conjectures and plans.

## Contents and provenance

- `outputs/`: all 185 saved deliverables, with Markdown links adapted for GitHub and the output manifest refreshed after those publication edits.
- `work/h16p-audit/`: supporting scripts, raw results, logs, coverage ledgers, replay dependencies, historical extracts, and source citation ledgers.
- `requirements.txt`: numerical package versions recorded for the original Python 3.14 environment.
- `publication-manifest.json`: file hashes for this archive, excluding the manifest itself; integrity checks are not mathematical validation.

The source repository was pinned at `45f4ea9b4ab448bdd36702036244faa4f15c9819`. Existing source files are unchanged by this publication. Historical usage and automation records describe the local session at its stopping point, not the present state of another machine.

## Reproduction

Use a separate copy of this archive and install the requirements in a Python 3.14 environment. Supporting scripts live together in `work/h16p-audit/`; keep their `replay/` directory intact. For scripts that use the inherited repository, clone H16P at the pinned commit into this archive’s `work/H16P/` directory. Consult each script and its corresponding report before replaying it: some entry points overwrite their adjacent result files, and some searches are expensive. Saved intermediate data are included so that replay does not require rerunning the entire investigation. No claim is made that every script was rerun after publication packaging.

Downloaded third-party papers, extracted full texts and page images, the Python environment, and bytecode caches are omitted. Source JSON ledgers retain URLs and precise claims. `publication-exclusions.json` lists the omitted source downloads.

## Limitations

Reading coverage is incomplete at the literal raw-data level: 264 files were fully read, three figures visually inspected, and 120 raw numerical files were not fully read literally. All raw numerical files received machine inspection. Several arguments depend on published theorems whose original proofs were not independently rederived. Numerical searches were not interval certified. The original counterexample objective remains open.
