#!/usr/bin/env python3
"""Exact rational certificate for one universal Q4 weighted-lobe point.

No floating-point arithmetic, search, quadrature library, or third-party
package is used.  This proves three primitive zeros, not five original
Abelian-integral zeros.  See notes_certificate_second.md for the tail proof.
"""

import argparse
from fractions import Fraction as Q
import hashlib
import json
import os
from pathlib import Path
import resource
import time

for key in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS",
            "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ[key] = "1"
resource.setrlimit(resource.RLIMIT_CPU, (10, 10))


PARAMETERS = {
    "A": "1243911778077/1000000000000",
    "B": "-86917392526/1000000000000",
    "eta": "1460428426173/1000000000000",
}
WITNESSES = (Q(1, 8), Q(3, 8), Q(5, 8), Q(7, 8))


def rational_string(value):
    value = Q(value)
    return f"{value.numerator}/{value.denominator}"


def decimal_integer(value, places):
    sign = "-" if value < 0 else ""
    value = abs(value)
    scale = 10 ** places
    return f"{sign}{value // scale}.{value % scale:0{places}d}"


def outward_decimal_interval(lower, upper, places=30):
    """Round lower down and upper up using integer arithmetic only."""
    scale = 10 ** places
    lower_integer = lower.numerator * scale // lower.denominator
    upper_integer = -((-upper.numerator * scale) // upper.denominator)
    return [decimal_integer(lower_integer, places),
            decimal_integer(upper_integer, places)]


def primitive_enclosure(A, B, eta, t, N=256):
    """Enclose H(t), including the series coefficients indexed 0 through N.

    F=sum(f_n*t^n), f_0=1 and f_(n+1)/f_n<1.  If K denotes the
    companion hypergeometric period, D=F-K=sum(d_n*t^n), where d_0=0
    and d_n=6*n*f_n/(6*n-1).  The exact identity t*F*M=D gives the
    integrated series used here.  A geometric majorant bounds every
    omitted term; all inputs and outputs remain Fractions.
    """
    A, B, eta, t = map(Q, (A, B, eta, t))
    if not 0 < t < 1 or N < 1:
        raise ValueError("require 0 < t < 1 and N >= 1")
    t2, t3 = t*t, t*t*t
    f_n, t_n, partial = Q(1), Q(1), Q(0)
    for n in range(N + 1):
        partial += f_n*t_n*((A-1)*t2/(n+2) + B*t3/(n+3))
        if n:
            partial += f_n*t_n*Q(6*n, 6*n-1)*(t2/(n+2)-eta*t/(n+1))
        f_n *= Q((6*n+1)*(6*n+5), 36*(n+1)**2)
        t_n *= t
    # f_n and t_n now mean f_(N+1) and t^(N+1).
    majorant = (abs(A-1)*t2/(N+3) + abs(B)*t3/(N+4)
                + Q(6*(N+1), 6*(N+1)-1)
                * (t2/(N+3) + abs(eta)*t/(N+2)))
    tail = f_n*t_n*majorant/(1-t)
    return partial-tail, partial+tail, tail


def certify(N=256):
    A, B, eta = (Q(PARAMETERS[name]) for name in ("A", "B", "eta"))
    assert 1 < eta < Q(54, 31)
    q_at_zero = A-1-eta/6
    assert q_at_zero > 0
    # Exact endpoint identity; pi is positive and is not numerically evaluated.
    pi_times_H_at_one = Q(18, 85085)*(9061*A+6289*B-2431*eta-7242)
    assert pi_times_H_at_one < 0
    records = []
    for index, t in enumerate(WITNESSES):
        lower, upper, tail = primitive_enclosure(A, B, eta, t, N)
        expected = 1 if index % 2 == 0 else -1
        assert (lower > 0 if expected > 0 else upper < 0)
        assert (lower > Q(1, 10**6) if expected > 0 else upper < -Q(1, 10**6))
        enclosure = outward_decimal_interval(lower, upper)
        # Also check the displayed decimal enclosure, not merely the exact one.
        assert Q(enclosure[0]) <= lower <= upper <= Q(enclosure[1])
        assert (Q(enclosure[0]) > 0 if expected > 0 else Q(enclosure[1]) < 0)
        records.append({
            "t": rational_string(t), "H_enclosure": enclosure,
            "certified_sign": expected,
            "absolute_tail_upper_bound": outward_decimal_interval(Q(0), tail)[1],
        })
    # For |delta A|,|delta B|,|delta eta| <= radius, F <= 1/(1-t)
    # and 0 < M <= 1 imply |delta H| <= radius*1519/192 at t <= 7/8.
    box_radius = Q(1, 10**7)
    perturbation_bound = Q(1519, 192)*box_radius
    assert perturbation_bound < Q(8, 10**7) < Q(1, 10**6)
    return {
        "format": "H16P-Q4-LOBE-CERTIFICATE-v1",
        "status": "RIGOROUS_UNIVERSAL_PRIMITIVE_ONLY",
        "parameters": PARAMETERS,
        "evaluator": {
            "arithmetic": "Python fractions.Fraction and integer rounding",
            "last_included_index_N": N,
            "decimal_places": 30,
            "error_bound": "positive hypergeometric coefficient geometric majorant",
            "script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        },
        "witnesses": records,
        "primitive_root_intervals": [
            {"endpoints": [rational_string(a), rational_string(b)], "open": True}
            for a, b in zip(WITNESSES[:-1], WITNESSES[1:])
        ],
        "exact_checks": {
            "eta_in_corrected_strip": True,
            "q_at_zero": rational_string(q_at_zero),
            "pi_times_H_at_one": rational_string(pi_times_H_at_one),
        },
        "certified_parameter_box": {
            "center": PARAMETERS,
            "closed_linf_radius": rational_string(box_radius),
            "uniform_witness_perturbation_bound": rational_string(perturbation_bound),
            "base_point_absolute_sign_margin_greater_than": "1/1000000",
            "entire_box_satisfies_all_three_strict_lobe_inequalities": True,
        },
        "proved_using_Q4_STRUCTURE_and_anchored_Rolle": {
            "primitive_distinct_simple_interior_zeros": 3,
            "auxiliary_q_distinct_simple_interior_zeros": 3,
            "all_three_strict_weighted_lobe_inequalities": True,
            "lobe_region_has_an_explicit_rational_interior_point": True,
        },
        "original_Abelian_integral_five_zero_claim": False,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path,
                        default=Path(__file__).parent / "data/second_lobe_certificate.json")
    args = parser.parse_args()
    start = time.process_time()
    record = certify()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(record, indent=2)+"\n")
    print("Exact rational weighted-lobe certificate: PASS")
    for witness in record["witnesses"]:
        print(f"t={witness['t']}: H in {witness['H_enclosure']}")
    print("Three primitive zeros and three strict lobe inequalities certified.")
    print("No five-zero original-integral claim.")
    print(f"CPU seconds: {time.process_time()-start:.3f}")


if __name__ == "__main__":
    main()
