#!/usr/bin/env python3
"""Identify the terminal event of the row-7 fold continuation exactly.

The steepest-descent fold-preserving path drove h_i -> 0 while the two folds
coalesced.  Test the CUSP hypothesis quantitatively and then solve
D = D_s = D_ss = 0 exactly in (a11, a01, a10, s), checking D_sss != 0.
"""
import json
import mpmath as mp
from engine import Engine
from cusp import Cusp
import foldcont as FC
mp.mp.dps = 40
FREE = ("a11", "a01", "a10")
def Tval(m): return m["a11"] + m["a01"] - 2*m["a"] - 1
def detJ(m): return 2*m["a"] - m["a01"] - m["a10"] - 2*m["a20"]

L = [json.loads(l) for l in open("ledger_opus/foldrun_row7_p2.jsonl") if "h_i" in l]
print("cubic-scaling test along the path  (delta = (s_f-s_i)/2):")
print("   step   delta        |h_i|        |D_ss(s_f)|   |h_i|/delta^3   |D_ss|/delta")
for r in L:
    d = (mp.mpf(r["s_f"]) - mp.mpf(r["s_i"]))/2
    h = abs(mp.mpf(r["h_i"])); dss = abs(mp.mpf(r["Dss_f"]))
    print("   %-6s %-12s %-12s %-13s %-15s %s"
          % (r["step"], mp.nstr(d,6), mp.nstr(h,6), mp.nstr(dss,6),
             mp.nstr(h/d**3,6), mp.nstr(dss/d,6)))
last = L[-1]
mu = {p: mp.mpf(last[p]) for p in FC.PARAMS}
s = (mp.mpf(last["s_f"]) + mp.mpf(last["s_i"]))/2

eng = Engine(); print("\nengine:", eng.banner)
c = Cusp(eng, mu["a"], mu["a20"], side=1)
print("\nNewton on  D = D_s = D_ss = 0  in (a11, a01, a10, s):")
for it in range(40):
    r = FC.val(c, mu, s)
    if r["status"] != "OK": print("   RETURN FAILED", r["status"]); break
    rr = abs(s-1)
    F = [r["D"], r["Dx"], r["Dxx"]]
    res = max(abs(F[0]), abs(F[1])*rr, abs(F[2])*rr*rr/2)
    print("   it%-3d scaled-res=%.4e  D=%.3e D_s=%.3e D_ss=%.3e  s=%.12f"
          % (it, float(res), float(F[0]), float(F[1]), float(F[2]), float(s)), flush=True)
    if res < mp.mpf("1e-28"): break
    J = [[mp.mpf(0)]*4 for _ in range(3)]
    for j, p in enumerate(FREE):
        h = mp.mpf("1e-11")*max(mp.mpf(1), abs(mu[p]))
        m1, m2 = dict(mu), dict(mu); m1[p] += h; m2[p] -= h
        r1, r2 = FC.val(c, m1, s), FC.val(c, m2, s)
        if r1["status"] != "OK" or r2["status"] != "OK": print("   jac failed"); break
        J[0][j] = (r1["D"]-r2["D"])/(2*h); J[1][j] = (r1["Dx"]-r2["Dx"])/(2*h)
        J[2][j] = (r1["Dxx"]-r2["Dxx"])/(2*h)
    J[0][3], J[1][3], J[2][3] = r["Dx"], r["Dxx"], r["Dxxx"]
    # 3 equations, 4 unknowns: fix a10 (least-effective control), solve 3x3
    M = mp.matrix(3,3)
    for i in range(3):
        M[i,0]=J[i][0]; M[i,1]=J[i][1]; M[i,2]=J[i][3]
    try: d = mp.lu_solve(M, mp.matrix(F))
    except Exception as e: print("   singular:", e); break
    lam = mp.mpf(1); dm = max(abs(d[k]) for k in range(3))
    if dm > mp.mpf("0.05"): lam = mp.mpf("0.05")/dm
    mu["a11"] -= lam*d[0]; mu["a01"] -= lam*d[1]; s -= lam*d[2]

r = FC.val(c, mu, s)
a00 = mu["a01"]+mu["a11"]-mu["a10"]-mu["a20"]-mu["a"]
print("\n*** TERMINAL OBJECT ***")
for p in FC.PARAMS: print("   %-4s = %s" % (p, mp.nstr(mu[p], 25)))
print("   a00  = %s" % mp.nstr(a00, 25))
print("   s_c  = %s" % mp.nstr(s, 20))
print("   D=%s  D_s=%s  D_ss=%s  D_sss=%s  D_ssss=%s"
      % (mp.nstr(r["D"],5), mp.nstr(r["Dx"],5), mp.nstr(r["Dxx"],5),
         mp.nstr(r["Dxxx"],12), mp.nstr(r.get("Dxxxx",mp.mpf(0)),12)))
print("   T = %s   detJ = %s" % (mp.nstr(Tval(mu),12), mp.nstr(detJ(mu),12)))
inv = FC.inventory(c, mu, step=mp.mpf("0.005"))
print("   inventory  domain=%s" % [mp.nstr(v,8) for v in inv["domain"]])
print("      cycles     = %s" % [mp.nstr(v,12) for v in inv["cycles"]])
print("      stationary = %s" % [mp.nstr(v,12) for v in inv["stationary"]])
rem = FC.remote_cycles(mu)
print("   remote: %s cycles=%s"
      % ("none" if rem["second_focus"] is None else "(%.4f,%.4f) tr %.4f" % rem["second_focus"], rem["cycles"]))
json.dump(dict(**{p: mp.nstr(mu[p],25) for p in FC.PARAMS}, a00=mp.nstr(a00,25),
               s_c=mp.nstr(s,25), D=mp.nstr(r["D"],10), Ds=mp.nstr(r["Dx"],10),
               Dss=mp.nstr(r["Dxx"],10), Dsss=mp.nstr(r["Dxxx"],15),
               T=mp.nstr(Tval(mu),15), detJ=mp.nstr(detJ(mu),15),
               cycles=[mp.nstr(v,15) for v in inv["cycles"]],
               stationary=[mp.nstr(v,15) for v in inv["stationary"]],
               remote=rem["cycles"]), open("ledger_opus/cusp_end.json","w"), indent=1)
print("\ncalls:", eng.ncalls); eng.close()
