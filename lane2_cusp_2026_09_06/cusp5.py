#!/usr/bin/env python3
"""Cusp unfolding with the trace controlled AND D_ss held at zero.

cusp4 controlled (D, D_s, T) but left D_ss free.  Moving the trace needs a
sizeable move in `a` (a10 and a20 do not affect T at all), and that reintroduces
D_ss(s_c) != 0, so the local picture is a shifted cubic rather than the symmetric
three-root one -- which is why cusp4 found two cycles, not three.  Constrain it:

    D(s_c) = 0 ,  D_s(s_c) = b1 = -c eps^2/2 ,  D_ss(s_c) = 0 ,  T = tau

4 equations in the 4 free parameters (a, a11, a01, a10); a20 stays fixed.
Then D(s_c+u) = b1 u + (c/6) u^3 exactly to third order: roots u = 0, +-sqrt(3)eps.
With tau < 0 the focus is stable, D dips negative just outside it and must climb
back to reach the triple's neighbourhood from above -- one extra sign change,
i.e. FOUR cycles, if the intervening maximum stays positive.  Count and see.
"""
import json
import mpmath as mp
from engine import Engine
from cusp import Cusp
import foldcont as FC
mp.mp.dps = 40
FREE = ("a", "a11", "a01", "a10")
def Tval(m): return m["a11"] + m["a01"] - 2*m["a"] - 1
def detJ(m): return 2*m["a"] - m["a01"] - m["a10"] - 2*m["a20"]

C = json.load(open("ledger_opus/cusp_end.json"))
MU0 = {p: mp.mpf(C[p]) for p in FC.PARAMS}; SC = mp.mpf(C["s_c"]); c3 = mp.mpf(C["Dsss"])
T0 = mp.mpf(C["T"])
eng = Engine(); print("engine:", eng.banner, flush=True)
c = Cusp(eng, MU0["a"], MU0["a20"], side=1)
CASES = [("0.05","0"), ("0.05","-1e-5"), ("0.05","-1e-4"), ("0.05","-4e-4"),
         ("0.10","-1e-4"), ("0.10","-4e-4"), ("0.10","-1e-3"), ("0.20","-1e-3")]
rows = []
for eps, tau in CASES:
    e, ta = mp.mpf(eps), mp.mpf(tau)
    b1 = -c3*e*e/2
    mu = dict(MU0); ok = False; res = None
    for it in range(150):
        r = FC.val(c, mu, SC)
        if r["status"] != "OK":
            print("eps=%s tau=%s RETURN FAILED %s" % (eps,tau,r["status"]), flush=True); break
        rr = abs(SC-1)
        F = [r["D"], r["Dx"] - b1, r["Dxx"], Tval(mu) - ta]
        res = max(abs(F[0]), abs(F[1])*rr, abs(F[2])*rr*rr/2, abs(F[3]))
        if res < mp.mpf("1e-26"): ok = True; break
        J = [[mp.mpf(0)]*4 for _ in range(4)]
        bad = False
        for j, p in enumerate(FREE):
            h = mp.mpf("1e-11")*max(mp.mpf(1), abs(mu[p]))
            m1, m2 = dict(mu), dict(mu); m1[p] += h; m2[p] -= h
            r1, r2 = FC.val(c, m1, SC), FC.val(c, m2, SC)
            if r1["status"] != "OK" or r2["status"] != "OK": bad = True; break
            J[0][j] = (r1["D"]-r2["D"])/(2*h)
            J[1][j] = (r1["Dx"]-r2["Dx"])/(2*h)
            J[2][j] = (r1["Dxx"]-r2["Dxx"])/(2*h)
        if bad: break
        J[3] = [mp.mpf(-2), mp.mpf(1), mp.mpf(1), mp.mpf(0)]
        M = mp.matrix(4,4)
        for i in range(4):
            for jj in range(4): M[i,jj] = J[i][jj]
        try: d = mp.lu_solve(M, mp.matrix(F))
        except Exception as ex: print("singular:", ex, flush=True); break
        lam = mp.mpf(1); dm = max(abs(d[k]) for k in range(4))
        if dm > mp.mpf("0.03"): lam = mp.mpf("0.03")/dm
        for j, p in enumerate(FREE): mu[p] -= lam*d[j]
    if not ok:
        print("eps=%-6s tau=%-8s : NO_CONVERGE res=%s" % (eps, tau, mp.nstr(res,6) if res else "-"), flush=True)
        continue
    inv = FC.inventory(c, mu, step=mp.mpf("0.01"))
    print("\neps=%-6s tau=%-8s  D=%s D_s=%s (tgt %s) D_ss=%s  T=%s detJ=%s"
          % (eps, tau, mp.nstr(r["D"],5), mp.nstr(r["Dx"],8), mp.nstr(b1,8),
             mp.nstr(r["Dxx"],5), mp.nstr(Tval(mu),8), mp.nstr(detJ(mu),8)), flush=True)
    print("   a=%s a20=%s a11=%s a01=%s a10=%s"
          % tuple(mp.nstr(mu[p],18) for p in FC.PARAMS))
    print("   domain %s" % [mp.nstr(v,7) for v in inv["domain"]])
    print("   LOCAL cycles (%d): %s" % (len(inv["cycles"]), [mp.nstr(v,12) for v in inv["cycles"]]), flush=True)
    print("   stationary   (%d): %s" % (len(inv["stationary"]), [mp.nstr(v,10) for v in inv["stationary"]]), flush=True)
    rows.append(dict(eps=eps, tau=tau, T=mp.nstr(Tval(mu),15), detJ=mp.nstr(detJ(mu),15),
                     mu={p: mp.nstr(mu[p],25) for p in FC.PARAMS},
                     a00=mp.nstr(mu["a01"]+mu["a11"]-mu["a10"]-mu["a20"]-mu["a"],25),
                     local=[mp.nstr(v,18) for v in inv["cycles"]],
                     stationary=[mp.nstr(v,15) for v in inv["stationary"]],
                     n=len(inv["cycles"])))
    json.dump(rows, open("ledger_opus/cusp5.json","w"), indent=1)
print("\ncalls:", eng.ncalls); eng.close()
