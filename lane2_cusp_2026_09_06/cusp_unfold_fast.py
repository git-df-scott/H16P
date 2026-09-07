#!/usr/bin/env python3
"""Fast pass of the cusp unfolding: local inventory only (the remote count is
measured separately by cusp_unfold.py / validate_remote.py)."""
import json
import mpmath as mp
from engine import Engine
from cusp import Cusp
import foldcont as FC
mp.mp.dps = 40
def Tval(m): return m["a11"] + m["a01"] - 2*m["a"] - 1

C = json.load(open("ledger_opus/cusp_end.json"))
MU = {p: mp.mpf(C[p]) for p in FC.PARAMS}; SC = mp.mpf(C["s_c"]); c3 = mp.mpf(C["Dsss"])
eng = Engine(); print("engine:", eng.banner, flush=True)
cu = Cusp(eng, MU["a"], MU["a20"], side=1); sc = FC.scales(MU)
g = FC.grad_mu(cu, MU, SC, which=("D","Dx"))
gT = [mp.mpf(v)*sc[j] for j, v in enumerate([-2, 0, 1, 1, 0])]
A = mp.matrix(3,5)
for j in range(5): A[0,j]=g["D"][j]; A[1,j]=g["Dx"][j]; A[2,j]=gT[j]
GG = A*A.T
rows=[]
for eps, tau in [("0.05","0"), ("0.10","0"), ("0.10","3e-4"), ("0.20","3e-4"), ("0.20","1e-3")]:
    e, ta = mp.mpf(eps), mp.mpf(tau)
    b1 = -c3*e*e/2
    b = mp.matrix([mp.mpf(0), b1, ta - Tval(MU)])
    lam = mp.lu_solve(GG, b)
    d = [sum(lam[i]*A[i,j] for i in range(3)) for j in range(5)]
    mu = {p: MU[p] + d[j]*sc[j] for j, p in enumerate(FC.PARAMS)}
    c2 = Cusp(eng, mu["a"], mu["a20"], side=1)
    r = FC.val(c2, mu, SC)
    if r["status"] != "OK":
        print("eps=%s tau=%s RETURN FAILED %s" % (eps,tau,r["status"]), flush=True); continue
    inv = FC.inventory(c2, mu, step=mp.mpf("0.01"))
    print("\neps=%-6s tau=%-8s T=%s  D(s_c)=%s D_s(s_c)=%s (target %s)"
          % (eps, tau, mp.nstr(Tval(mu),8), mp.nstr(r["D"],6), mp.nstr(r["Dx"],6), mp.nstr(b1,6)), flush=True)
    print("   domain %s" % [mp.nstr(v,7) for v in inv["domain"]])
    print("   LOCAL cycles (%d): %s" % (len(inv["cycles"]), [mp.nstr(v,12) for v in inv["cycles"]]), flush=True)
    print("   stationary   (%d): %s" % (len(inv["stationary"]), [mp.nstr(v,10) for v in inv["stationary"]]), flush=True)
    rows.append(dict(eps=eps, tau=tau, T=mp.nstr(Tval(mu),12),
                     mu={p: mp.nstr(mu[p],25) for p in FC.PARAMS},
                     a00=mp.nstr(mu["a01"]+mu["a11"]-mu["a10"]-mu["a20"]-mu["a"],25),
                     local=[mp.nstr(v,18) for v in inv["cycles"]],
                     stationary=[mp.nstr(v,15) for v in inv["stationary"]]))
    json.dump(rows, open("ledger_opus/cusp_unfold_fast.json","w"), indent=1)
print("\ncalls:", eng.ncalls); eng.close()
