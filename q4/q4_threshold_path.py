#!/usr/bin/env python3
"""One exact primitive-anchor path and its bounded threshold certificate.

The path is analytic by Q4_LOBE_REGION. Numerical coefficient construction
uses closed PF moments, including near the loop. The certificate freezes
rational coefficients and verifies primitive signs with exact series bounds;
it makes no original-five-zero claim. No scan is run on import or execution.
"""
import os
for _name in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS",
              "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ[_name] = "1"

from fractions import Fraction as Q
from pathlib import Path
import hashlib
import json
import resource
import time

import mpmath as mp


def threshold_anchors(r):
    """Exact analytic path: first primitive root is identically r."""
    if not 0 < r < 1:
        raise ValueError("require 0<r<1")
    return r, (1+r)/2, (3+r)/4


def primitive_basis_closed(t):
    """Return (K0,K1,K2,K3) using three exact PF moment identities.

    Uses caller's mpmath precision with 20 guard digits. No quadrature or long
    t-series is required. Small-t cancellation is handled with extra digits.
    """
    t = mp.mpf(t)
    if not 0 <= t <= 1:
        raise ValueError("require 0<=t<=1")
    if t == 0:
        return (mp.mpf(0),)*4
    extra = 20 + max(0, int(-4*mp.log10(t)))
    with mp.workdps(mp.mp.dps + extra):
        if t == 1:
            j0 = mp.mpf(18)/(5*mp.pi)
            j1 = mp.mpf(738)/(385*mp.pi)
            j2 = mp.mpf(113202)/(85085*mp.pi)
            boundary = mp.mpf(0)
        else:
            f = mp.hyp2f1(mp.mpf(1)/6, mp.mpf(5)/6, 1, t)
            fp = mp.mpf(5)/36*mp.hyp2f1(mp.mpf(7)/6, mp.mpf(11)/6, 2, t)
            boundary = t*(1-t)*f
            derivative_boundary = t*(1-t)*fp
            j0 = mp.mpf(36)/5*derivative_boundary
            j1 = (t*derivative_boundary-boundary+j0)/(2+mp.mpf(5)/36)
            j2 = (t*t*derivative_boundary-2*t*boundary+4*j1)/(6+mp.mpf(5)/36)
        result = (j1,j2,6*j0-11*j1-6*boundary,
                  12*j1-17*j2-6*t*boundary)
    return tuple(+x for x in result)


def from_primitive_anchors_closed(anchors):
    """Numerical coefficients of the exact unique normalized anchor solution."""
    anchors = tuple(mp.mpf(t) for t in anchors)
    if not 0 < anchors[0] < anchors[1] < anchors[2] < 1:
        raise ValueError("three ordered interior primitive anchors required")
    # Near coincident anchors, retain guard digits proportional to separation.
    gap = min(anchors[1]-anchors[0],anchors[2]-anchors[1],1-anchors[2])
    extra = 25+max(0,int(-4*mp.log10(gap)))
    with mp.workdps(mp.mp.dps+extra):
        rows = [primitive_basis_closed(t) for t in anchors]
        matrix = mp.matrix([[x[0],x[1],-x[2]] for x in rows])
        rhs = mp.matrix([x[0]-x[3] for x in rows])
        result = tuple(mp.lu_solve(matrix,rhs))
    return tuple(+x for x in result)


def coefficients_from_r(r):
    return from_primitive_anchors_closed(threshold_anchors(mp.mpf(r)))


def primitive_value_closed(t, coefficients):
    A,B,eta = map(mp.mpf,coefficients)
    row = primitive_basis_closed(t)
    return (A-1)*row[0]+B*row[1]-eta*row[2]+row[3]


def rounded_rationals(values, digits=24):
    scale = 10**digits
    return tuple(Q(int(mp.nint(v*scale)),scale) for v in values)


def certify_frozen():
    """Exact rational late-root point; no floating output enters its proof."""
    from q4_lobe_certificate import (primitive_enclosure,
        outward_decimal_interval,rational_string)
    params = {"A":"1210581187245108808/1000000000000000000",
              "B":"-125731163118386543/1000000000000000000",
              "eta":"1212211767298108636/1000000000000000000"}
    A,B,eta = (Q(params[name]) for name in ("A","B","eta"))
    assert 1 < eta < Q(54,31)
    q0=A-1-eta/6
    assert q0>0
    endpoint=Q(18,85085)*(9061*A+6289*B-2431*eta-7242)
    assert endpoint<0
    witnesses=(Q(23,32),Q(13,16),Q(29,32),Q(31,32))
    assert witnesses[0]>Q(5,11)
    records=[]
    for index,t in enumerate(witnesses):
        lo,hi,tail=primitive_enclosure(A,B,eta,t,N=1024)
        sign=1 if index%2==0 else -1
        margin=lo if sign>0 else -hi
        assert margin>Q(1,10**5)
        bounds=outward_decimal_interval(lo,hi)
        assert Q(bounds[0])<=lo<=hi<=Q(bounds[1])
        records.append({"t":rational_string(t),"H_enclosure":bounds,
            "certified_sign":sign,
            "absolute_tail_upper_bound":outward_decimal_interval(Q(0),tail)[1]})
    threshold=Q(5,11)
    lo,hi,tail=primitive_enclosure(A,B,eta,threshold,N=512)
    assert lo>0
    radius=Q(1,10**8)
    perturbation=Q(122047,3072)*radius
    assert perturbation<Q(1,10**6)<Q(1,10**5)
    helper=Path(__file__).with_name("q4_lobe_certificate.py")
    return {
        "format":"H16P-Q4-THRESHOLD-CERTIFICATE-v1",
        "status":"RIGOROUS_LATE_PRIMITIVE_ROOT_ONLY",
        "parameters":params,
        "discovery_anchor_r":"3/4",
        "discovery_primitive_anchors":["3/4","7/8","15/16"],
        "arithmetic":"fractions.Fraction and directed integer rounding",
        "last_included_index_N":1024,
        "script_sha256":hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "series_helper_sha256":hashlib.sha256(helper.read_bytes()).hexdigest(),
        "witnesses":records,
        "primitive_root_intervals":[
            {"endpoints":[rational_string(a),rational_string(b)],"open":True}
            for a,b in zip(witnesses[:-1],witnesses[1:])],
        "threshold_evaluation":{"t":"5/11",
            "H_enclosure":outward_decimal_interval(lo,hi),"certified_sign":1},
        "exact_checks":{"eta_in_corrected_strip":True,
            "q_at_zero":rational_string(q0),
            "pi_times_H_at_one":rational_string(endpoint),
            "first_primitive_root_strictly_greater_than":"23/32",
            "first_root_exceeds_5_over_11":True,
            "H_positive_through_5_over_11_by_global_zero_count":True},
        "certified_parameter_box":{"closed_linf_radius":rational_string(radius),
            "uniform_witness_perturbation_bound":rational_string(perturbation),
            "base_sign_margins_strictly_greater_than":"1/100000",
            "entire_box_in_lobe_region_and_after_threshold":True},
        "original_Abelian_integral_five_zero_claim":False,
    }


def main():
    resource.setrlimit(resource.RLIMIT_CPU,(10,10))
    os.nice(10)
    start=time.process_time()
    record=certify_frozen()
    output=Path(__file__).with_name("data")/"third_threshold_certificate.json"
    output.write_text(json.dumps(record,indent=2)+"\n")
    print("Exact threshold-escape primitive certificate: PASS")
    for row in record["witnesses"]:
        print(row["t"],row["H_enclosure"])
    print("No original-five-zero claim. CPU seconds:",time.process_time()-start)


if __name__ == "__main__":
    main()
