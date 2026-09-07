#!/usr/bin/env python3
"""Separated double fold with the SEPARATION PINNED BY AN EQUATION.

The earlier double-fold attempts used ad-hoc guards on |s1-s2| and |s-1|, and
every one of them slid onto the focus.  Replace the guards by an equation:

    F(a11,a01,a10,s_f,s_i) = ( D(s_f), D_s(s_f), D(s_i), D_s(s_i), s_f-s_i-g0 )

5 equations, 5 unknowns, with (a,a20) held as continuation parameters.  A root
is a genuine separated double fold at prescribed separation g0: the guard can no
longer be "hit", because the separation is now part of the solution.
Nondegeneracy at a root: D_ss(s_f) != 0 and D_ss(s_i) != 0.
"""
import json, sys
import mpmath as mp
from engine import Engine
from cusp import Cusp
import foldcont as FC
mp.mp.dps = 40
FREE = ("a11", "a01", "a10")

def Tval(m): return m["a11"] + m["a01"] - 2*m["a"] - 1
def detJ(m): return 2*m["a"] - m["a01"] - m["a10"] - 2*m["a20"]

rec = [json.loads(l) for l in open("ledger_opus/fold_start.jsonl")]
rec = [r for r in rec if r["row"] == 7][-1]
MU0 = {p: mp.mpf(rec[p]) for p in FC.PARAMS}
SF0, SI0 = mp.mpf(rec["s_f"]), mp.mpf(rec["stationary"][0])
G0 = mp.mpf(sys.argv[1]) if len(sys.argv) > 1 else SF0 - SI0

eng = Engine(); print("engine:", eng.banner, flush=True)
lg = open("ledger_opus/gapfold.jsonl", "a")

def solve(mu, s_f, s_i, g0, itmax=60, tol=mp.mpf("1e-24")):
    mu = dict(mu)
    c = Cusp(eng, mu["a"], mu["a20"], side=1)
    for it in range(itmax):
        rf = FC.val(c, mu, s_f); ri = FC.val(c, mu, s_i)
        if rf["status"] != "OK": return None, ("RETURN_FAILED_F", rf["status"])
        if ri["status"] != "OK": return None, ("RETURN_FAILED_I", ri["status"])
        rr_f, rr_i = abs(s_f-1), abs(s_i-1)
        F = [rf["D"], rf["Dx"], ri["D"], ri["Dx"], s_f - s_i - g0]
        res = max(abs(F[0]), abs(F[1])*rr_f, abs(F[2]), abs(F[3])*rr_i, abs(F[4]))
        if it % 6 == 0 or res < tol:
            print("     it%-3d res=%.4e  D(s_f)=%.3e D(s_i)=%.3e  s=(%.6f,%.6f)"
                  % (it, float(res), float(F[0]), float(F[2]), float(s_f), float(s_i)), flush=True)
        if res < tol: return (mu, s_f, s_i, rf, ri), None
        # Jacobian: columns a11,a01,a10 by central differences; then s_f, s_i
        J = [[mp.mpf(0)]*5 for _ in range(5)]
        for j, p in enumerate(FREE):
            h = mp.mpf("1e-11")*max(mp.mpf(1), abs(mu[p]))
            m1, m2 = dict(mu), dict(mu); m1[p] += h; m2[p] -= h
            f1, f2 = FC.val(c, m1, s_f), FC.val(c, m2, s_f)
            i1, i2 = FC.val(c, m1, s_i), FC.val(c, m2, s_i)
            if any(x["status"] != "OK" for x in (f1, f2, i1, i2)):
                return None, ("JAC_RETURN_FAILED", p)
            J[0][j] = (f1["D"]-f2["D"])/(2*h);   J[1][j] = (f1["Dx"]-f2["Dx"])/(2*h)
            J[2][j] = (i1["D"]-i2["D"])/(2*h);   J[3][j] = (i1["Dx"]-i2["Dx"])/(2*h)
        J[0][3] = rf["Dx"];  J[1][3] = rf["Dxx"]; J[4][3] = mp.mpf(1)
        J[2][4] = ri["Dx"];  J[3][4] = ri["Dxx"]; J[4][4] = mp.mpf(-1)
        M = mp.matrix(5, 5)
        for i in range(5):
            for j in range(5): M[i, j] = J[i][j]
        try: d = mp.lu_solve(M, mp.matrix(F))
        except Exception as e: return None, ("SINGULAR", str(e))
        dmax = max(abs(d[k]) for k in range(5))
        lam = mp.mpf(1)
        if dmax > mp.mpf("0.3"): lam = mp.mpf("0.3")/dmax        # damping
        for j, p in enumerate(FREE): mu[p] -= lam*d[j]
        s_f -= lam*d[3]; s_i -= lam*d[4]
        if s_i <= 1 or s_f <= 1: return None, ("LEFT_DOMAIN", None)
    return None, ("NO_CONVERGE", None)

print("\nseed: row 7,  g0 = %s" % mp.nstr(G0, 12))
out, err = solve(MU0, SF0, SI0, G0)
if err:
    print("   -> %s %s" % err, flush=True)
    lg.write(json.dumps(dict(g0=mp.nstr(G0,12), status=err[0], detail=str(err[1])))+"\n")
else:
    mu, s_f, s_i, rf, ri = out
    print("\n*** SEPARATED DOUBLE FOLD (g0 = %s) ***" % mp.nstr(G0, 12))
    for p in FC.PARAMS: print("   %-4s = %s" % (p, mp.nstr(mu[p], 25)))
    a00 = mu["a01"]+mu["a11"]-mu["a10"]-mu["a20"]-mu["a"]
    print("   a00  = %s" % mp.nstr(a00, 25))
    print("   s_f = %s  D_ss = %s" % (mp.nstr(s_f,20), mp.nstr(rf["Dxx"],12)))
    print("   s_i = %s  D_ss = %s" % (mp.nstr(s_i,20), mp.nstr(ri["Dxx"],12)))
    print("   T = %s   detJ = %s" % (mp.nstr(Tval(mu),10), mp.nstr(detJ(mu),10)))
    inv = FC.inventory(eng and Cusp(eng, mu["a"], mu["a20"], side=1), mu, step=mp.mpf("0.005"))
    print("   inventory cycles=%s stationary=%s"
          % ([mp.nstr(v,12) for v in inv["cycles"]], [mp.nstr(v,12) for v in inv["stationary"]]))
    rem = FC.remote_cycles(mu)
    print("   remote: %s cycles=%s"
          % ("none" if rem["second_focus"] is None else "(%.4f,%.4f) tr %.4f" % rem["second_focus"], rem["cycles"]))
    lg.write(json.dumps(dict(g0=mp.nstr(G0,12), status="DOUBLE_FOLD",
        **{p: mp.nstr(mu[p],25) for p in FC.PARAMS}, a00=mp.nstr(a00,25),
        s_f=mp.nstr(s_f,25), s_i=mp.nstr(s_i,25),
        Dss_f=mp.nstr(rf["Dxx"],15), Dss_i=mp.nstr(ri["Dxx"],15),
        T=mp.nstr(Tval(mu),15), detJ=mp.nstr(detJ(mu),15),
        cycles=[mp.nstr(v,15) for v in inv["cycles"]],
        stationary=[mp.nstr(v,15) for v in inv["stationary"]],
        remote=rem["cycles"]))+"\n")
lg.flush()
print("\ncalls:", eng.ncalls, flush=True); eng.close()
