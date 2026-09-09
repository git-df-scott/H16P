#!/usr/bin/env python3
"""Cusp unfolding WITH the trace as a third control, solved by NEWTON.

The 3-row system [d_mu D(s_c); d_mu D_s(s_c); grad_mu T] is nonsingular but
ill-conditioned (row-normalised sigma_min/sigma_max = 2.7e-4).  A single linear
step in it overshoots by 100x -- that is what the first attempt did, and it is
why that attempt reported only two cycles.  Ill-conditioned is not singular, so
solve it properly: Newton in (a11, a01, a10) on

    D(s_c) = 0 ,   D_s(s_c) = b1 = -c*eps^2/2 ,   T = tau .

With tau < 0 the displacement dips negative just outside the focus and must
climb back to reach the cusp triple from above -- one extra sign change, i.e.
FOUR cycles, if the intervening maximum stays positive.  Test it.
"""
import json, sys
import mpmath as mp
from engine import Engine
from cusp import Cusp
import foldcont as FC
mp.mp.dps = 40
FREE = ("a", "a11", "a01")   # a10 and a20 do not move T at all; a does (dT/da = -2)
def Tval(m): return m["a11"] + m["a01"] - 2*m["a"] - 1
def detJ(m): return 2*m["a"] - m["a01"] - m["a10"] - 2*m["a20"]

C = json.load(open("ledger_opus/cusp_end.json"))
MU0 = {p: mp.mpf(C[p]) for p in FC.PARAMS}; SC = mp.mpf(C["s_c"]); c3 = mp.mpf(C["Dsss"])
eng = Engine(); print("engine:", eng.banner, flush=True)
c0 = Cusp(eng, MU0["a"], MU0["a20"], side=1)

CASES = [("0.05","-1e-5"), ("0.05","-1e-4"), ("0.05","-5e-4"), ("0.05","-2e-3"),
         ("0.10","-5e-4"), ("0.10","-2e-3"), ("0.20","-2e-3"), ("0.20","-8e-3")]
rows = []
for eps, tau in CASES:
    e, ta = mp.mpf(eps), mp.mpf(tau)
    b1 = -c3*e*e/2
    mu = dict(MU0); c = c0; ok = False; res = None
    for it in range(120):
        r = FC.val(c, mu, SC)
        if r["status"] != "OK": print("eps=%s tau=%s RETURN FAILED %s" % (eps,tau,r["status"]), flush=True); break
        F = [r["D"], r["Dx"] - b1, Tval(mu) - ta]
        res = max(abs(F[0]), abs(F[1])*abs(SC-1), abs(F[2]))
        if res < mp.mpf("1e-26"): ok = True; break
        J = [[mp.mpf(0)]*3 for _ in range(3)]
        for j, p in enumerate(FREE):
            h = mp.mpf("1e-11")*max(mp.mpf(1), abs(mu[p]))
            m1, m2 = dict(mu), dict(mu); m1[p] += h; m2[p] -= h
            r1, r2 = FC.val(c, m1, SC), FC.val(c, m2, SC)
            if r1["status"] != "OK" or r2["status"] != "OK": break
            J[0][j] = (r1["D"]-r2["D"])/(2*h); J[1][j] = (r1["Dx"]-r2["Dx"])/(2*h)
            J[2][0] = mp.mpf(-2); J[2][1] = mp.mpf(1); J[2][2] = mp.mpf(1)   # dT/d(a,a11,a01)
        M = mp.matrix(3,3)
        for i in range(3):
            for j in range(3): M[i,j] = J[i][j]
        try: d = mp.lu_solve(M, mp.matrix(F))
        except Exception as ex: print("singular:", ex, flush=True); break
        lam = mp.mpf(1); dm = max(abs(d[k]) for k in range(3))
        if dm > mp.mpf("0.05"): lam = mp.mpf("0.05")/dm
        for j, p in enumerate(FREE): mu[p] -= lam*d[j]
    if not ok:
        print("eps=%s tau=%-8s : NO_CONVERGE res=%s" % (eps, tau, mp.nstr(res,6) if res else "-"), flush=True)
        continue
    inv = FC.inventory(c, mu, step=mp.mpf("0.01"))
    print("\neps=%-6s tau=%-8s  D(s_c)=%s D_s(s_c)=%s (tgt %s)  T=%s  detJ=%s"
          % (eps, tau, mp.nstr(r["D"],5), mp.nstr(r["Dx"],8), mp.nstr(b1,8),
             mp.nstr(Tval(mu),8), mp.nstr(detJ(mu),8)), flush=True)
    print("   domain %s" % [mp.nstr(v,7) for v in inv["domain"]])
    print("   LOCAL cycles (%d): %s" % (len(inv["cycles"]), [mp.nstr(v,12) for v in inv["cycles"]]), flush=True)
    print("   stationary   (%d): %s" % (len(inv["stationary"]), [mp.nstr(v,10) for v in inv["stationary"]]), flush=True)
    rows.append(dict(eps=eps, tau=tau, T=mp.nstr(Tval(mu),15), detJ=mp.nstr(detJ(mu),15),
                     mu={p: mp.nstr(mu[p],25) for p in FC.PARAMS},
                     a00=mp.nstr(mu["a01"]+mu["a11"]-mu["a10"]-mu["a20"]-mu["a"],25),
                     local=[mp.nstr(v,18) for v in inv["cycles"]],
                     stationary=[mp.nstr(v,15) for v in inv["stationary"]],
                     n=len(inv["cycles"])))
    json.dump(rows, open("ledger_opus/cusp4.json","w"), indent=1)
print("\ncalls:", eng.ncalls); eng.close()
