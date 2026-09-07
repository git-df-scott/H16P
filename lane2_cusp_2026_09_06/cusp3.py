#!/usr/bin/env python3
"""Unfold the finite-amplitude cusp using ONLY the well-conditioned rank-2
system (D, D_s) at s_c, solved by Newton rather than by one linear step.

The 3-row system that also controls T is nearly rank-deficient
(row-normalised sigma_min/sigma_max = 2.7e-4), so a linear step in it
overshoots by orders of magnitude -- that is what the first attempt did.
T is left where the cusp put it, T = +5.49e-4 > 0, which with l1 < 0 is the
side on which a Hopf cycle exists.

Target:  D(s_c) = 0,  D_s(s_c) = -c*eps^2/2   =>  D ~ (c/6)(u^3 - 3 eps^2 u),
three simple roots at u = 0, +-sqrt(3)*eps.
"""
import json
import mpmath as mp
from engine import Engine
from cusp import Cusp
import foldcont as FC
mp.mp.dps = 40
FREE = ("a11", "a01")
def Tval(m): return m["a11"] + m["a01"] - 2*m["a"] - 1
def detJ(m): return 2*m["a"] - m["a01"] - m["a10"] - 2*m["a20"]

C = json.load(open("ledger_opus/cusp_end.json"))
MU0 = {p: mp.mpf(C[p]) for p in FC.PARAMS}; SC = mp.mpf(C["s_c"]); c3 = mp.mpf(C["Dsss"])
eng = Engine(); print("engine:", eng.banner, flush=True)
rows = []
for eps in ("0.05", "0.10", "0.20", "0.35", "0.50"):
    e = mp.mpf(eps); b1 = -c3*e*e/2
    mu = dict(MU0)
    c = Cusp(eng, mu["a"], mu["a20"], side=1)
    ok = False
    for it in range(30):
        r = FC.val(c, mu, SC)
        if r["status"] != "OK": print("eps=%s RETURN FAILED %s" % (eps, r["status"]), flush=True); break
        F = [r["D"], r["Dx"] - b1]
        res = max(abs(F[0]), abs(F[1])*abs(SC-1))
        if res < mp.mpf("1e-26"): ok = True; break
        J = [[mp.mpf(0)]*2 for _ in range(2)]
        for j, p in enumerate(FREE):
            h = mp.mpf("1e-11")*max(mp.mpf(1), abs(mu[p]))
            m1, m2 = dict(mu), dict(mu); m1[p] += h; m2[p] -= h
            r1, r2 = FC.val(c, m1, SC), FC.val(c, m2, SC)
            if r1["status"] != "OK" or r2["status"] != "OK": break
            J[0][j] = (r1["D"]-r2["D"])/(2*h); J[1][j] = (r1["Dx"]-r2["Dx"])/(2*h)
        det = J[0][0]*J[1][1] - J[0][1]*J[1][0]
        if det == 0: print("singular"); break
        d0 = (F[0]*J[1][1] - F[1]*J[0][1])/det
        d1 = (J[0][0]*F[1] - J[1][0]*F[0])/det
        lam = mp.mpf(1); dm = max(abs(d0), abs(d1))
        if dm > mp.mpf("0.02"): lam = mp.mpf("0.02")/dm
        mu["a11"] -= lam*d0; mu["a01"] -= lam*d1
    if not ok:
        print("eps=%s : NO_CONVERGE (res=%s)" % (eps, mp.nstr(res,6)), flush=True); continue
    inv = FC.inventory(c, mu, step=mp.mpf("0.01"))
    print("\neps=%-6s  D(s_c)=%s  D_s(s_c)=%s (target %s)  T=%s  detJ=%s"
          % (eps, mp.nstr(r["D"],6), mp.nstr(r["Dx"],8), mp.nstr(b1,8),
             mp.nstr(Tval(mu),8), mp.nstr(detJ(mu),8)), flush=True)
    print("   predicted roots at s = %s, %s, %s"
          % (mp.nstr(SC-mp.sqrt(3)*e,10), mp.nstr(SC,10), mp.nstr(SC+mp.sqrt(3)*e,10)))
    print("   domain %s" % [mp.nstr(v,7) for v in inv["domain"]])
    print("   LOCAL cycles (%d): %s" % (len(inv["cycles"]), [mp.nstr(v,12) for v in inv["cycles"]]), flush=True)
    print("   stationary   (%d): %s" % (len(inv["stationary"]), [mp.nstr(v,10) for v in inv["stationary"]]), flush=True)
    rows.append(dict(eps=eps, T=mp.nstr(Tval(mu),15), detJ=mp.nstr(detJ(mu),15),
                     mu={p: mp.nstr(mu[p],25) for p in FC.PARAMS},
                     a00=mp.nstr(mu["a01"]+mu["a11"]-mu["a10"]-mu["a20"]-mu["a"],25),
                     local=[mp.nstr(v,18) for v in inv["cycles"]],
                     stationary=[mp.nstr(v,15) for v in inv["stationary"]]))
    json.dump(rows, open("ledger_opus/cusp3.json","w"), indent=1)
print("\ncalls:", eng.ncalls); eng.close()
