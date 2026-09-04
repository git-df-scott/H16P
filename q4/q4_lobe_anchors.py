#!/usr/bin/env python3
"""Six frozen tests on a derived primitive-root spine; no coefficient sweep.

Three prescribed primitive roots determine the direction before kappa is
chosen. The exact certificate is in q4_lobe_certificate.py. This script's
floating diagnostics and sampled original-I crossings are not certificates.
"""
import os
import resource

for name in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS",
             "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ[name] = "1"

import json
import time
from pathlib import Path

import mpmath as mp
import numpy as np
from scipy.optimize import brentq
from scipy.special import hyp2f1

from q4_reconstruction import reconstruct, original_values, mu_from_universal


def primitive_basis(t):
    """Positive convergent period series, evaluated at caller mp precision."""
    t = mp.mpf(t)
    if not 0 < t < 1:
        raise ValueError("require 0<t<1")
    f, power = mp.mpf(1), t*t
    u0 = u1 = u2 = u3 = mp.mpf(0)
    for n in range(1000):
        u0 += f*power/(n+2)
        u1 += f*power*t/(n+3)
        if n:
            diff = 6*n*f/(6*n-1)
            u2 += diff*power/t/(n+1)
            u3 += diff*power/(n+2)
        f *= (n+mp.mpf(1)/6)*(n+mp.mpf(5)/6)/(n+1)**2
        power *= t
        if n > 10 and abs(f*power) < mp.mpf("1e-55"):
            break
    else:
        raise RuntimeError("fixed series ceiling reached")
    return u0, u1, u2, u3


def from_primitive_anchors(anchors):
    """Numerical 3x3 solve; exact invertibility is proved in Q4_LOBE_REGION."""
    if not 0 < anchors[0] < anchors[1] < anchors[2] < 1:
        raise ValueError("three ordered interior anchors required")
    rows = [primitive_basis(t) for t in anchors]
    matrix = mp.matrix([[r[0], r[1], -r[2]] for r in rows])
    rhs = mp.matrix([r[0]-r[3] for r in rows])
    return tuple(mp.lu_solve(matrix, rhs))


def primitive_value(t, coefficients):
    A, B, eta = coefficients
    if t == 1:
        return 18*(9061*A+6289*B-2431*eta-7242)/(85085*mp.pi)
    u0, u1, u2, u3 = primitive_basis(t)
    return (A-1)*u0+B*u1-eta*u2+u3


def q_float(t, coefficients):
    A, B, eta = map(float, coefficients)
    if t == 0:
        return A-1-eta/6
    F = hyp2f1(1/6, 5/6, 1, t)
    Fp = 5/36*hyp2f1(7/6, 11/6, 2, t)
    M = 1-6*(1-t)*Fp/F
    return A+B*t-1+(t-eta)*M


def run_diagnostic():
    resource.setrlimit(resource.RLIMIT_CPU, (10, 10))
    os.nice(10)
    started = time.process_time()
    mp.mp.dps = 60
    records = []
    for r in (mp.mpf(1)/8, mp.mpf(1)/4, mp.mpf(3)/8):
        anchors = (r, mp.mpf(1)/2, 1-r)
        coefficients = from_primitive_anchors(anchors)
        qroots = []
        for left, right in zip((mp.mpf(0), *anchors[:-1]), anchors):
            qroots.append(brentq(lambda t: q_float(t, coefficients),
                                 float(left), float(right), xtol=1e-13))
        critical = [primitive_value(mp.mpf(t), coefficients) for t in qroots]
        endpoint = primitive_value(1, coefficients)
        lobes = (critical[0], critical[0]-critical[1],
                 critical[2]-critical[1], critical[2]-endpoint)
        for a in (.5, .75):
            solution = reconstruct(a, *map(float, coefficients), t_end=.99)
            # A fixed compact interval only. The two omitted endpoint slivers
            # are not claimed to contain no roots.
            samples = np.linspace(.01, .99, 65)
            values = original_values(a, solution, samples)
            scale = max(abs(values))
            retained = [(float(t), float(v)) for t, v in zip(samples, values)
                        if abs(v) > 1e-10*scale]
            signs = [np.sign(v) for t, v in retained]
            changes = sum(x != y for x, y in zip(signs, signs[1:]))
            alternating = [retained[0]]
            for point in retained[1:]:
                if point[1]*alternating[-1][1] < 0:
                    alternating.append(point)
            record = {
                "label": "NUMERICAL_ONLY",
                "primitive_anchors": [str(t) for t in anchors],
                "A_B_eta": [mp.nstr(c, 48) for c in coefficients],
                "q_roots": qroots,
                "weighted_lobes": [mp.nstr(v, 35) for v in lobes],
                "primitive_at_critical_points": [mp.nstr(v, 35) for v in critical],
                "primitive_at_one": mp.nstr(endpoint, 35),
                "a": a, "kappa": 1/(1-a),
                "original_mu": mu_from_universal(1/(1-a), *map(float, coefficients)).tolist(),
                "t_interval": [.01, .99], "sample_count": len(samples),
                "sampled_sign_changes": int(changes),
                "alternating_sample_witnesses_t_I": alternating,
                "sample_minimum": float(min(values)),
                "sample_maximum": float(max(values)),
            }
            records.append(record)
            print("anchors", record["primitive_anchors"], "kappa", record["kappa"],
                  "sampled I crossings", changes)
            if changes >= 5:
                print("POSSIBLE FIVE-ZERO LEAD: stop and independently verify")
                return records, time.process_time()-started
    return records, time.process_time()-started


if __name__ == "__main__":
    records, elapsed = run_diagnostic()
    output = Path(__file__).with_name("data")/"second_spine_diagnostic.json"
    output.write_text(json.dumps({"status": "NUMERICAL_ONLY", "cpu_seconds": elapsed,
        "warning": "Fixed compact samples are neither root counts nor endpoint exclusions.",
        "records": records}, indent=2)+"\n")
    print("CPU seconds:", elapsed)
