#!/usr/bin/env python3
"""Bounded replay of two prescribed primitive-root paths and scalar kappa shots.

No coefficient cube, random sampling, or generic optimization is used.
All outputs are NUMERICAL_ONLY. In particular the narrow P sign bands
are not interval-certified, and five original zeros are not asserted.
"""
import os
for key in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS",
            "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ[key] = "1"

import hashlib
import json
import resource
import time
from pathlib import Path

import mpmath as mp
from scipy.optimize import brentq
from q4_threshold_path import (coefficients_from_r, threshold_anchors,
                               from_primitive_anchors_closed)
from q4_reconstruction import (reconstruct, green_coordinates, original_values,
                               mu_from_universal, initial_weighted_derivative)


def evaluate(a, coefficients, anchors, end):
    sol = reconstruct(a, *coefficients, t_end=end)
    p = lambda t: float(green_coordinates(a, sol, t)[1])
    return sol, p


def record_solution(r, coefficients_mp, anchors, a, end, path):
    coefficients = tuple(map(float, coefficients_mp))
    sol, p = evaluate(a, coefficients, anchors, end)
    probes = (0, *anchors, end)
    ps = [p(t) for t in probes]
    roots = [brentq(p, left, right, xtol=1e-14)
             for left, right in zip(probes, probes[1:]) if p(left)*p(right)<0]
    zs = [float(green_coordinates(a, sol, t)[0]) for t in roots]
    samples = sorted(set((.1, .3, .5, .7, *anchors, end)))
    vals = list(map(float, original_values(a, sol, samples)))
    sign_changes = sum(x*y<0 for x, y in zip(vals, vals[1:]))
    row = {
        "path": path, "r": str(r), "primitive_anchors": anchors,
        "A_B_eta": [mp.nstr(x, 40) for x in coefficients_mp],
        "a": a, "kappa": 1/(1-a),
        "original_mu_approximate": mu_from_universal(1/(1-a), *coefficients).tolist(),
        "P0": initial_weighted_derivative(a, *coefficients),
        "P_probe_points": probes, "P_at_probes": ps,
        "P_crossings": roots, "Z_at_crossings": zs,
        "H_residual_at_prescribed_roots": [float(sol.sol(t)[0]) for t in anchors],
        "I_t_value": list(zip(samples, vals)),
        "sampled_original_sign_changes": sign_changes,
        "S1_numerical_pattern": len(roots)==4 and all(v*s>0 for v,s in
                                                        zip(ps,(1,-1,1,-1,1))),
        "five_original_zero_candidate": sign_changes>=5,
    }
    if sign_changes>=5:
        candidate = Path(__file__).with_name("data")/"third_candidate_trigger.json"
        candidate.write_text(json.dumps(row, indent=2)+"\n")
        raise RuntimeError("Five-sign lead frozen; stop exploration for independent evaluation")
    return row


def tuned_shot(rr, power_path=False):
    r = mp.mpf(rr)
    if power_path:
        eps = 1-r
        ts = (r, 1-eps**2, 1-eps**3)
        co = from_primitive_anchors_closed(ts)
        end = 1-float(eps**3)/16
    else:
        ts = threshold_anchors(r)
        co = coefficients_from_r(r)
        end = 1-(1-float(r))/32
    anchors, cs = tuple(map(float, ts)), tuple(map(float, co))
    def middle(a):
        return evaluate(a, cs, anchors, end)[1](anchors[1])
    # A scalar continuation condition, P(tau2)=0, determines a boundary.
    astar = brentq(middle, .12 if power_path else .1199, .99999, xtol=3e-14)
    sol, p = evaluate(astar, cs, anchors, end)
    gap = min(-p(anchors[0]), -p(anchors[2]))
    slope = (middle(astar+1e-6)-middle(astar-1e-6))/2e-6
    assert gap>0 and slope>0
    # Move into the thin four-crossing band by a specified fraction of its
    # smaller extremum margin. This does not claim an interval enclosure.
    a = astar+.4*gap/slope
    row = record_solution(rr, co, anchors, a, end,
        "(r,1-(1-r)^2,1-(1-r)^3)" if power_path else "(r,(1+r)/2,(3+r)/4)")
    row["P_extremum_gap"] = gap
    row["shooting_boundary_a"] = astar
    assert row["S1_numerical_pattern"]
    assert max(row["Z_at_crossings"])<-.002, "frozen negative maximum changed"
    return row


def save(name, rows, started):
    record = {
        "status": "NUMERICAL_ONLY", "cpu_seconds": time.process_time()-started,
        "script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "warning": "P sign bands are floating diagnostics, not certificates. "
                   "Fixed original samples do not cover endpoint slivers. "
                   "Extremely late anchors require arbitrary precision; this "
                   "replay stops at 1-r=1e-4.",
        "rows": rows,
    }
    (Path(__file__).with_name("data")/name).write_text(json.dumps(record, indent=2)+"\n")


def main():
    resource.setrlimit(resource.RLIMIT_CPU, (10, 10))
    os.nice(10)
    mp.mp.dps = 60
    started = time.process_time()
    rows = []
    for rr in ("0.5", "0.75", "0.9", "0.99"):
        r = mp.mpf(rr)
        co, ts = coefficients_from_r(r), tuple(map(float, threshold_anchors(r)))
        for k in (1.137, 2., 4., 16.):
            rows.append(record_solution(rr, co, ts, 1-1/k,
                         1-(1-float(r))/32, "(r,(1+r)/2,(3+r)/4)"))
    save("third_initial_shoot.json", rows, started)
    affine = [tuned_shot(r) for r in ("0.5", "0.75", "0.9", "0.99", "0.9999")]
    save("third_tuned_shoot.json", affine, started)
    power = [tuned_shot(r, True) for r in ("0.6", "0.75", "0.9")]
    save("third_shape_shoot.json", power, started)
    for row in affine+power:
        print(row["path"], "r", row["r"], "kappa", row["kappa"],
              "first Z maximum", row["Z_at_crossings"][0])
    print("Eight targeted S1 patterns; all fail first S2 maximum numerically.")
    print("No five-zero original candidate. CPU seconds:", time.process_time()-started)


if __name__ == "__main__":
    main()
