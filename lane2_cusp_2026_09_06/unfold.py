#!/usr/bin/env python3
"""Realise the rank-2 unfolding of the Hopf + cycle-fold organizer.

D_ss(s_f) < 0, so D ~ D(s_f) + (D_ss/2)(s-s_f)^2 has TWO roots iff D(s_f) > 0.
l1 < 0 (D_sss(1) = -5.645e-3), so a cycle is born from A iff T > 0.
Rank[grad_mu T ; d_mu D(s_f)] = 2, so (T, D(s_f)) = (eta_T, eta_D) is reachable.
Take the least-norm delta_mu and count EVERY cycle, local and remote.
"""
import json
import mpmath as mp
from engine import Engine
from cusp import Cusp
import foldcont as FC
mp.mp.dps = 40

def Tval(mu): return mu["a11"] + mu["a01"] - 2*mu["a"] - 1
def detJ(mu): return 2*mu["a"] - mu["a01"] - mu["a10"] - 2*mu["a20"]

H = json.load(open("ledger_opus/hopf_fold.json"))
MU = {p: mp.mpf(H[p]) for p in FC.PARAMS}; S_F = mp.mpf(H["s_f"])
eng = Engine(); print("engine:", eng.banner)
c = Cusp(eng, MU["a"], MU["a20"], side=1)
sc = FC.scales(MU)
gT = [mp.mpf(v)*sc[j] for j, v in enumerate([-2, 0, 1, 1, 0])]
g = FC.grad_mu(c, MU, S_F, which=("D",)); bf = g["D"]
A = mp.matrix(2, 5)
for j in range(5): A[0,j] = gT[j]; A[1,j] = bf[j]
G = A*A.T
rows = []
for etaT, etaD in [("1e-5","1e-9"), ("1e-4","1e-8"), ("3e-4","4e-8"), ("1e-3","1e-7")]:
    b = mp.matrix([mp.mpf(etaT), mp.mpf(etaD)])
    lam = mp.lu_solve(G, b)
    d = [lam[0]*A[0,j] + lam[1]*A[1,j] for j in range(5)]     # least-norm, scaled
    mu = {p: MU[p] + d[j]*sc[j] for j, p in enumerate(FC.PARAMS)}
    T = Tval(mu)
    r = FC.val(c, mu, S_F)
    if r["status"] != "OK":
        print("eta=(%s,%s): return failed %s" % (etaT, etaD, r["status"])); continue
    inv = FC.inventory(c, mu, step=mp.mpf("0.005"))
    rem = FC.remote_cycles(mu)
    print("\neta_T=%s eta_D=%s" % (etaT, etaD))
    print("   achieved  T = %s   D(s_f) = %s   detJ = %s"
          % (mp.nstr(T,8), mp.nstr(r["D"],8), mp.nstr(detJ(mu),8)))
    print("   domain    %s" % ([mp.nstr(v,8) for v in inv["domain"]] if inv["domain"] else None))
    print("   LOCAL cycles (%d): %s" % (len(inv["cycles"]), [mp.nstr(v,12) for v in inv["cycles"]]))
    print("   stationary   (%d): %s" % (len(inv["stationary"]), [mp.nstr(v,10) for v in inv["stationary"]]))
    print("   REMOTE: focus %s cycles=%s"
          % ("none" if rem["second_focus"] is None else "(%.4f,%.4f) tr %.4f" % rem["second_focus"], rem["cycles"]))
    tot = len(inv["cycles"]) + (rem["cycles"] or 0)
    print("   TOTAL simple cycles = %d" % tot)
    rows.append(dict(etaT=etaT, etaD=etaD, T=mp.nstr(T,12), Dsf=mp.nstr(r["D"],12),
                     local=[mp.nstr(v,15) for v in inv["cycles"]],
                     stationary=[mp.nstr(v,12) for v in inv["stationary"]],
                     remote=rem["cycles"], total=tot,
                     mu={p: mp.nstr(mu[p],25) for p in FC.PARAMS}))
json.dump(rows, open("ledger_opus/unfold.json","w"), indent=1)
print("\ncalls:", eng.ncalls); eng.close()
