# Q4 controls and bounded screen

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
