#!/usr/bin/env python3
"""Unfold the finite-amplitude cusp (multiplicity-3 limit cycle) reached by the
row-7 fold continuation, and count EVERY cycle at the unfolded parameter.

At the cusp:  D(s_c)=D_s(s_c)=D_ss(s_c)=0,  c := D_sss(s_c) != 0.
Locally  D(s_c+u) ~ (c/6)u^3 + b1 u + b0,  three simple roots iff
4p^3+27q^2 < 0 with p=6b1/c, q=6b0/c.  Take q=0, p=-3 eps^2:
     b1 = -c eps^2/2,  b0 = 0   ->  roots at u = 0, +-sqrt(3) eps.
Controls needed: (D(s_c), D_s(s_c), T).  Verify rank 3 first, then realise
(b0, b1, T) = (0, -c eps^2/2, tau) by the least-norm delta_mu and inventory.
"""
import json, sys
import mpmath as mp
from engine import Engine
from cusp import Cusp
import foldcont as FC
mp.mp.dps = 40
def Tval(m): return m["a11"] + m["a01"] - 2*m["a"] - 1
def detJ(m): return 2*m["a"] - m["a01"] - m["a10"] - 2*m["a20"]

C = json.load(open("ledger_opus/cusp_end.json"))
MU = {p: mp.mpf(C[p]) for p in FC.PARAMS}; SC = mp.mpf(C["s_c"])
c3 = mp.mpf(C["Dsss"])
eng = Engine(); print("engine:", eng.banner, flush=True)
cu = Cusp(eng, MU["a"], MU["a20"], side=1)
sc = FC.scales(MU)

g = FC.grad_mu(cu, MU, SC, which=("D", "Dx"))
gT = [mp.mpf(v)*sc[j] for j, v in enumerate([-2, 0, 1, 1, 0])]
A = mp.matrix(3, 5)
for j in range(5): A[0,j] = g["D"][j]; A[1,j] = g["Dx"][j]; A[2,j] = gT[j]
print("\nunfolding matrix rows (scaled):")
for i, nm in enumerate(("d_mu D(s_c)", "d_mu D_s(s_c)", "grad_mu T")):
    print("   %-15s %s" % (nm, [mp.nstr(A[i,j],8) for j in range(5)]))
# row-normalised singular values
N = mp.matrix(3,5)
for i in range(3):
    n = mp.sqrt(sum(A[i,j]**2 for j in range(5)))
    for j in range(5): N[i,j] = A[i,j]/n
G = N*N.T
ev = mp.eigsy(mp.matrix([[G[i,j] for j in range(3)] for i in range(3)]), eigvals_only=True)
sv = sorted([mp.sqrt(max(mp.mpf(0), e)) for e in ev], reverse=True)
print("   row-normalised singular values: %s   ratio sigma_min/sigma_max = %s"
      % ([mp.nstr(v,8) for v in sv], mp.nstr(sv[-1]/sv[0], 8)))
print("   => UNFOLDING RANK = %d" % sum(1 for v in sv if v > mp.mpf("1e-10")))

GG = A*A.T
rows = []
for eps, tau in [("0.02","0"), ("0.05","0"), ("0.10","0"),
                 ("0.05","1e-4"), ("0.10","1e-4"), ("0.10","-5e-4"), ("0.15","-5e-4")]:
    e, ta = mp.mpf(eps), mp.mpf(tau)
    b1 = -c3*e*e/2
    tgt = mp.matrix([mp.mpf(0), b1, ta - Tval(MU)])   # (D, D_s, T) targets, relative
    b = mp.matrix([tgt[0] - mp.mpf(0), tgt[1], tgt[2]])
    lam = mp.lu_solve(GG, b)
    d = [sum(lam[i]*A[i,j] for i in range(3)) for j in range(5)]
    mu = {p: MU[p] + d[j]*sc[j] for j, p in enumerate(FC.PARAMS)}
    r = FC.val(cu if (mu["a"]==MU["a"] and mu["a20"]==MU["a20"]) else Cusp(eng, mu["a"], mu["a20"], side=1), mu, SC)
    c2 = Cusp(eng, mu["a"], mu["a20"], side=1)
    r = FC.val(c2, mu, SC)
    if r["status"] != "OK":
        print("\neps=%s tau=%s : RETURN FAILED %s" % (eps, tau, r["status"])); continue
    inv = FC.inventory(c2, mu, step=mp.mpf("0.004"))
    rem = FC.remote_cycles(mu)
    print("\neps=%-6s tau=%-8s  achieved D=%s D_s=%s T=%s (target D_s=%s)"
          % (eps, tau, mp.nstr(r["D"],6), mp.nstr(r["Dx"],6), mp.nstr(Tval(mu),6), mp.nstr(b1,6)), flush=True)
    print("   domain %s" % [mp.nstr(v,7) for v in inv["domain"]])
    print("   LOCAL cycles (%d): %s" % (len(inv["cycles"]), [mp.nstr(v,12) for v in inv["cycles"]]))
    print("   stationary   (%d): %s" % (len(inv["stationary"]), [mp.nstr(v,10) for v in inv["stationary"]]))
    print("   REMOTE %s cycles=%s"
          % ("none" if rem["second_focus"] is None else "(%.4f,%.4f) tr %.4f" % rem["second_focus"], rem["cycles"]))
    print("   TOTAL = %d" % (len(inv["cycles"]) + (rem["cycles"] or 0)), flush=True)
    rows.append(dict(eps=eps, tau=tau, T=mp.nstr(Tval(mu),12), detJ=mp.nstr(detJ(mu),12),
                     mu={p: mp.nstr(mu[p],25) for p in FC.PARAMS},
                     a00=mp.nstr(mu["a01"]+mu["a11"]-mu["a10"]-mu["a20"]-mu["a"],25),
                     local=[mp.nstr(v,18) for v in inv["cycles"]],
                     stationary=[mp.nstr(v,15) for v in inv["stationary"]],
                     remote=rem["cycles"], total=len(inv["cycles"])+(rem["cycles"] or 0)))
    json.dump(rows, open("ledger_opus/cusp_unfold.json","w"), indent=1)
print("\ncalls:", eng.ncalls); eng.close()
