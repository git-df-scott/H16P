#!/usr/bin/env python3
"""Exact rational checks for the global anchor comparison, with no search."""
import os
import resource
from fractions import Fraction as Q

resource.setrlimit(resource.RLIMIT_CPU, (10, 10))
os.nice(10)
for variable in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS",
                 "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ[variable] = "1"


def run_checks():
    # These are pi*J_n(1), inherited exact PF moments.
    j0, j1, j2 = Q(18, 5), Q(738, 385), Q(113202, 85085)
    x_end = j2 / j1
    m_end = (6*j0 - 11*j1) / j1
    slope_end = (1-m_end) / (1-x_end)
    assert x_end == Q(6289, 9061)
    assert m_end == Q(11, 41)
    assert slope_end == Q(1105, 462)

    c_a, c_b, c_m = Q(3*1326, 1361360), Q(3*864, 1361360), Q(3*2431, 1361360)
    assert c_a == Q(9, 3080)
    assert c_b / c_a == Q(144, 221)
    assert c_m / c_a == Q(11, 6)
    margin = -Q(1, 6) - c_b/c_a*slope_end + c_m/c_a
    assert margin == Q(25, 231) and margin > 0

    A, B, eta = Q(94, 77), -Q(17, 77), Q(1)
    Y_star = 3*(1326*A+864*B-2431*eta-102)/1361360
    assert Y_star == -Q(3, 1232)
    endpoint_gain = Q(3, 2)*(25-1)/14784
    assert endpoint_gain + Y_star == 0

    # delta=1/64 makes both fractional powers rational.
    delta_sixth = Q(1, 2)
    compact_margin = Q(25, 118272)*(delta_sixth-delta_sixth**5/5)
    assert compact_margin == Q(395, 3784704) and compact_margin > 0

    print("Exact moment-curve endpoint and maximal slope: PASS")
    print("Exact center-functional ratios and positive 25/231 margin: PASS")
    print("Exact corner center datum and endpoint gain cancellation: PASS")
    print("Exact compact margin for delta=1/64: PASS (395/3784704)")
    print("Global monotonicity and Theorem N: analytic proof in notes_N_loop.md")


if __name__ == "__main__":
    run_checks()
