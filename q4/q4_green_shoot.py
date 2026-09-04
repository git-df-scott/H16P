#!/usr/bin/env python3
"""Replay one frozen Green shooting point. Numerical diagnostic, not a proof.

This point was selected by a one-dimensional shooting of P(0) through zero.
There is no parameter search in this replay. The exact theorem in
Q4_RECONSTRUCTION_GEOMETRY excludes five original zeros throughout its
universal coefficient box, for every kappa.
"""
import os
for name in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS",
             "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ[name] = "1"

import json
import resource
import time
from fractions import Fraction as Q
from pathlib import Path

from scipy.optimize import brentq
from q4_reconstruction import (reconstruct, original_values, green_coordinates,
                               initial_weighted_derivative, mu_from_universal)


def run_diagnostic():
    resource.setrlimit(resource.RLIMIT_CPU, (10, 10))
    os.nice(10)
    started = time.process_time()
    aq = Q(4735432856584775, 10**16)
    coefficients = (Q(1243911778077, 10**12), -Q(86917392526, 10**12),
                    Q(1460428426173, 10**12))
    a = float(aq)
    values = tuple(map(float, coefficients))
    sol = reconstruct(a, *values, t_end=.99)
    p = lambda t: float(green_coordinates(a, sol, t)[1])
    probes = (0., .25, .5, .75, .99)
    ps = [p(t) for t in probes]
    assert all(v*s > 0 for v, s in zip(ps, (1, -1, 1, -1, 1)))
    roots = [brentq(p, left, right, xtol=1e-13)
             for left, right in zip(probes, probes[1:])]
    z = [float(green_coordinates(a, sol, t)[0]) for t in roots]
    assert max(z) < -.0045, "frozen unsuccessful shooting diagnostic changed"
    witnesses = (1/8, 3/8, 5/8, 7/8, .99)
    record = {
        "status": "NUMERICAL_ONLY",
        "warning": "Floating P signs and roots are diagnostic, not certified. "
                   "The analytic first-primitive-zero theorem independently excludes "
                   "five original zeros for this coefficient box for every kappa.",
        "a_exact": str(aq), "kappa_exact": str(1/(1-aq)),
        "kappa_approximate": 1/(1-a),
        "A_B_eta_exact": list(map(str, coefficients)),
        "original_mu_approximate": mu_from_universal(1/(1-a), *values).tolist(),
        "P0_from_center_data": initial_weighted_derivative(a, *values),
        "P_probes_t_value": list(zip(probes, ps)),
        "P_roots_approximate": roots,
        "Z_at_P_roots": z,
        "original_I_probes_t_value": list(zip(witnesses,
                        map(float, original_values(a, sol, witnesses)))),
        "S1_diagnostic": "four P crossings found",
        "S2_diagnostic": "fails: all four Z extrema are negative",
        "five_zero_candidate": False,
        "cpu_seconds": time.process_time()-started,
    }
    output = Path(__file__).with_name("data")/"second_green_shoot.json"
    output.write_text(json.dumps(record, indent=2)+"\n")
    print(json.dumps(record, indent=2))


if __name__ == "__main__":
    run_diagnostic()
