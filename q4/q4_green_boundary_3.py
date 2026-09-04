#!/usr/bin/env python3
"""Four bounded shots on the primitive triple-contact boundary of L.

This is H=H'=H''=0 (q has an earlier simple root and a double root),
not the previously excluded auxiliary triple-root cusp. Numerical only.
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
from q4_threshold_path import primitive_basis_closed
from q4_reconstruction import reconstruct, green_coordinates


def confluent_coefficients(t):
    F = mp.hyp2f1(mp.mpf(1)/6, mp.mpf(5)/6, 1, t)
    Fp = mp.mpf(5)/36*mp.hyp2f1(mp.mpf(7)/6, mp.mpf(11)/6, 2, t)
    Fpp = ((2*t-1)*Fp+mp.mpf(5)/36*F)/(t*(1-t))
    M = 1-6*(1-t)*Fp/F
    Mp = 6*Fp/F-6*(1-t)*(Fpp/F-(Fp/F)**2)
    K0, K1, K2, K3 = primitive_basis_closed(t)
    matrix = mp.matrix([[K0,K1,-K2],[1,t,-M],[0,1,-Mp]])
    rhs = mp.matrix([K0-K3,1-t*M,-M-t*Mp])
    return tuple(mp.lu_solve(matrix, rhs))


def main():
    resource.setrlimit(resource.RLIMIT_CPU, (10,10))
    os.nice(10)
    mp.mp.dps = 65
    rows = []
    started = time.process_time()
    for value in ("0.6", "0.9", "0.99", "0.9999"):
        t = mp.mpf(value)
        co = confluent_coefficients(t)
        cs, tf = tuple(map(float,co)), float(t)
        def metric(a, full=False):
            sol = reconstruct(a,*cs,t_end=tf)
            Z, P = green_coordinates(a,sol,tf)
            return (float(Z),float(P)) if full else float(P)
        a = brentq(metric,.12,.99999,xtol=3e-14)
        Z, P = metric(a,True)
        assert Z < -.002 and abs(P) < 1e-14
        row = {"status":"NUMERICAL_ONLY","primitive_triple_contact":value,
            "A_B_eta":[mp.nstr(x,40) for x in co],"a":a,"kappa":1/(1-a),
            "Z_at_contact":Z,"P_at_contact":P}
        rows.append(row)
        print("H triple contact",value,"kappa",row["kappa"],"Z",Z)
    record = {"status":"NUMERICAL_ONLY","rows":rows,
        "cpu_seconds":time.process_time()-started,
        "script_sha256":hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "warning":"These are primitive contacts, not original-integral double roots. "
                  "All four shots fail the first positive Z maximum numerically."}
    (Path(__file__).with_name("data")/"third_confluent_shoot.json").write_text(
                                                       json.dumps(record,indent=2)+"\n")


if __name__ == "__main__":
    main()
