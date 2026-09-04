# Q4 controls and bounded reconstruction

## Current second-strike replay

The active workflow is the exact lobe certificate and Green reconstruction.
Run the following small replays after installing requirements:

    python q4/q4_structure_checks.py
    python q4/test_q4.py
    python q4/q4_lobe_certificate.py
    python q4/q4_reconstruction.py
    python q4/q4_lobe_anchors.py
    python q4/q4_green_shoot.py
    python q4/test_q4_second.py

The new scripts have ten-second CPU fuses; run them sequentially at lowered
priority with numerical libraries limited to one thread. The exact rational
lobe certificate proves three primitive zeros and a full coefficient box.
The reconstruction proof excludes five original zeros throughout that box
for every kappa. The floating anchor and shooting data are diagnostics only.
See ../ASTRA_SECOND_STRIKE.md and data/second_verification.txt.

## Historical controls and superseded smoke workflow

The scripts evaluate the four-dimensional complete-elliptic-integral space
used by Gavrilov--Iliev and Zhao for a generic quadratic Q4 center.

Install and reproduce:

    python -m pip install -r q4/requirements.txt
    python q4/q4_controls.py --output q4/data/controls.json
    python q4/q4_search.py --mode smoke --cpu-hours 0.02 \
      --candidate-mode triple --kappa-count 9 --samples-per-kappa 256 \
      --grid-points 121 --quad-order 64 --output q4/data/smoke.json

q4_integrals.py provides arbitrary-precision area quadrature, fast Gauss
quadrature, and an independent Hamiltonian-orbit/Green evaluator. q4_search.py
has a CPU-time fuse, applies Zhao's necessary five-zero filters, constructs
triple-zero coefficient directions, and records conditioning diagnostics.

All generated output is explicitly nonrigorous. The interval and Poincaré
JSON schemas define promotion interfaces; they do not contain certificates.

## Astra reasoning replay

The current handoff replaces production sweeps with bounded reasoning.
Run python q4/q4_structure_checks.py for the exact symbolic checks and
three diagnostic period evaluations. It fixes numerical libraries to one
thread and has a ten-second CPU ceiling. Run python q4/test_q4.py for the
corrected-filter regressions. Neither command launches a candidate search.
Historical smoke records predate the corrected filters; see the root
Q4_STRUCTURE.md and Q4_PARAMETERIZATION.md.
