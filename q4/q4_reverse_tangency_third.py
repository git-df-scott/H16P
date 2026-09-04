#!/usr/bin/env python3
"""Bounded NUMERICAL reverse-tangency diagnostics; no global certificate.

Only four preselected (a,t*) pairs are examined.  Each solves the two linear
conditions Y(t*)=Y'(t*)=0 for A=A0+A1*eta, B=B0+B1*eta.  For the two late
pairs, a one-dimensional ratio determinant detects sampled stationary
heights of H=U+eta*V.  A finite mesh cannot exclude missed stationary roots
or certify an entire coefficient line.  This is not a five-zero search.
"""
import os
for _name in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS",
              "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ[_name] = "1"
import resource
resource.setrlimit(resource.RLIMIT_CPU, (10, 10))
os.nice(10)
import argparse
import hashlib
import json
from pathlib import Path
import time
import numpy as np
from scipy.optimize import brentq
from scipy.special import hyp2f1
from q4_reconstruction import reconstruct


def period_basis(t):
    t = np.asarray(t)
    F = hyp2f1(1/6, 5/6, 1, t)
    Fp = (5/36)*hyp2f1(7/6, 11/6, 2, t)
    J0 = (36/5)*t*(1-t)*Fp
    J1 = (t*t*(1-t)*Fp-t*(1-t)*F+J0)/(77/36)
    J2 = (t**3*(1-t)*Fp-2*t*t*(1-t)*F+4*J1)/(221/36)
    K2 = 6*J0-11*J1-6*t*(1-t)*F
    K3 = 12*J1-17*J2-6*t*t*(1-t)*F
    M = 1-6*(1-t)*Fp/F
    return J1, J2, K2, K3, M


def line(a, tstar):
    values = [reconstruct(a, *c, t_end=tstar).sol(tstar)[1:3]
              for c in ((0,0,0), (1,0,0), (0,1,0), (0,0,1))]
    base = values[0]
    matrix = np.column_stack((values[1]-base, values[2]-base))
    c0 = np.linalg.solve(matrix, -base)
    c1 = np.linalg.solve(matrix, -(values[3]-base))
    return c0, c1, float(np.linalg.cond(matrix))


def necessary_eta_interval(c0, c1):
    A0, B0 = c0
    A1, B1 = c1
    lo, hi = 1., 54/31
    # Each pair represents a necessary strict affine inequality c+d*eta>0.
    inequalities = [
        (A0-7/6,A1), (85/31-A0,-A1), (B0+1,B1),
        (-49/744-B0,-B1), (A0-1,A1-1/6),
        (-A0-B0,1-A1-B1), (-B0-1/6,25/432-B1),
    ]
    for c, d in inequalities:
        if d > 0:
            lo = max(lo, -c/d)
        elif d < 0:
            hi = min(hi, -c/d)
        elif c <= 0:
            return None
    return [float(lo), float(hi)] if lo < hi else None


def bracket_roots(function, mesh):
    vals = function(mesh)
    roots = []
    for j in np.flatnonzero(vals[:-1]*vals[1:] < 0):
        root = float(brentq(function, mesh[j], mesh[j+1], xtol=5e-15))
        if not roots or abs(root-roots[-1]) > 1e-10:
            roots.append(root)
    return roots


def ratio_diagnostic(c0, c1, eta_interval):
    A0, B0 = c0
    A1, B1 = c1
    def parts(t):
        J1,J2,K2,K3,M = period_basis(t)
        U = (A0-1)*J1+B0*J2+K3
        V = A1*J1+B1*J2-K2
        qbase = A0+B0*t-1+t*M
        qdir = A1+B1*t-M
        return U,V,qbase,qdir
    def determinant(t):
        U,V,qbase,qdir = parts(t)
        return U*qdir-qbase*V
    # Fixed bounded mesh, with loop clustering.  This is a diagnostic only.
    mesh = np.unique(np.concatenate((np.linspace(.01,.99,801),
                                    1-np.geomspace(1e-2,1e-10,801))))
    stationary = []
    for t in bracket_roots(determinant, mesh):
        U,V,_,_ = parts(t)
        stationary.append({"t":t, "eta":float(-U/V)})
    poles = bracket_roots(lambda t:parts(t)[1], mesh)
    intercept = 9061*A0+6289*B0-7242
    slope = 9061*A1+6289*B1-2431
    result = {
        "stationary_heights_detected":stationary,
        "ratio_poles_detected":poles,
        "endpoint_H_numerator":[float(intercept),float(slope)],
        "endpoint_zero_eta":float(-intercept/slope),
        "determinant_mesh_range":[float(np.min(determinant(mesh))),
                                   float(np.max(determinant(mesh)))],
        "mesh_t_range":[float(mesh[0]),float(mesh[-1])],
        "mesh_size":len(mesh),
    }
    if eta_interval is not None:
        # Historical bounded line check, retained for direct replay.  It is
        # explicitly weaker than a rigorous ratio-stationary-root count.
        tt = np.unique(np.concatenate((np.linspace(.001,.99,801),
                                      1-np.geomspace(.01,1e-6,801))))
        U,V,_,_ = parts(tt)
        counts = {}
        for eta in np.linspace(*eta_interval,257):
            H = U+eta*V
            signs = np.sign(H[np.abs(H)>1e-12])
            n = int(np.sum(signs[:-1]*signs[1:] < 0))
            counts[str(n)] = counts.get(str(n),0)+1
        result["sampled_H_crossing_counts"] = counts
        result["eta_samples"] = 257
    return result


def diagnostic():
    start = time.process_time()
    rows = []
    for a,tstar in ((.75,.8),(.75,.95),(.9,.999),(.99,.999)):
        c0,c1,condition = line(a,tstar)
        interval = necessary_eta_interval(c0,c1)
        row = {"a":a,"kappa":1/(1-a),"t_star":tstar,
               "A_intercept":float(c0[0]),"A_slope":float(c1[0]),
               "B_intercept":float(c0[1]),"B_slope":float(c1[1]),
               "AB_matrix_condition":condition,
               "necessary_eta_interval_numerical":interval}
        if a > .75:
            row["ratio_diagnostic"] = ratio_diagnostic(c0,c1,interval)
        rows.append(row)
    return {
        "status":"NUMERICAL; no rigorous line or global exclusion",
        "scope":"four preselected reverse-tangency lines; two ratio diagnostics",
        "interpretation":"No three-root H profile was detected. Finite meshes can miss roots or thin windows.",
        "source_sha256":hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "cpu_limit_seconds":10,
        "cpu_seconds":time.process_time()-start,
        "rows":rows,
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write",type=Path,help="save the numerical record")
    args = parser.parse_args()
    output = diagnostic()
    rendered = json.dumps(output,indent=2)+"\n"
    if args.write:
        args.write.write_text(rendered)
    print(rendered)
