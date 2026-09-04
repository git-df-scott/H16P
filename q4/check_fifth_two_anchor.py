#!/usr/bin/env python3
"""Tiny exact checks for the Strike-5 two-anchor moment constants.

This checks rational algebra only. The analytic inequalities are proved in
notes_fifth_two_anchor.md; there is no quadrature or parameter search.
"""
import os
for name in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS",
             "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ[name] = "1"
import resource
resource.setrlimit(resource.RLIMIT_CPU, (10, 10))
os.nice(10)
from fractions import Fraction as Q
import json
import time


def verify():
    start = time.process_time()
    x = Q(144, 221)
    xe, me, ne = Q(6289, 9061), Q(2431, 9061), Q(1819, 9061)
    scale = Q(9, 3080)
    assert 0 < x < xe < 1
    assert me == Q(11, 41)
    assert scale == Q(3978, 1361360)
    assert scale*x == Q(2592, 1361360)
    assert scale*Q(11, 6) == Q(7293, 1361360)
    assert scale*Q(204, 221) == Q(3672, 1361360)
    slope_m = (1-me)/(1-xe)
    assert slope_m == Q(1105, 462)
    chord_m = Q(1, 6)+(me-Q(1, 6))*x/xe
    chord_n = ne*x/xe
    lower_v = scale*(Q(11, 6)-chord_m)
    lower_eta_zero = scale*(Q(204, 221)-chord_n)
    assert lower_v == Q(231, 50312)
    assert lower_eta_zero == Q(27, 12578)
    tangent_zero = Q(1, 6)+Q(25, 432)*x
    tangent_one = me+slope_m*(x-xe)
    assert tangent_zero == Q(271, 1326) > tangent_one == Q(1, 6)
    upper_v = scale*(Q(11, 6)-tangent_one)
    assert upper_v == Q(3, 616)
    return {
        "status": "PASS: exact rational identities only",
        "ell_V_strict_lower": str(lower_v),
        "ell_V_strict_upper": str(upper_v),
        "eta_zero_center_strict_lower": str(lower_eta_zero),
        "moment_endpoint_slope": str(slope_m),
        "cpu_seconds": time.process_time()-start,
    }


if __name__ == "__main__":
    print(json.dumps(verify(), indent=2))
