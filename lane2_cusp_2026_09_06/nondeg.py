#!/usr/bin/env python3
"""Nondegeneracy audit at the Hopf+cycle-fold organizer.

(1) engine convention:  D_s(1) = exp(pi*T/omega) - 1,  omega = sqrt(detJ - T^2/4)
(2) first Lyapunov coefficient l1 at T=0 (sign of D_sss(1)), and the T*l1 test
(3) rank[ grad_mu T ; d_mu D(s_f) ]
"""
import json
import mpmath as mp
from engine import Engine
from cusp import Cusp
import foldcont as FC
mp.mp.dps = 40

def Tval(mu):  return mu["a11"] + mu["a01"] - 2*mu["a"] - 1
def detJ(mu):  return 2*mu["a"] - mu["a01"] - mu["a10"] - 2*mu["a20"]

H = json.load(open("ledger_opus/hopf_fold.json"))
MU = {p: mp.mpf(H[p]) for p in FC.PARAMS}
S_F = mp.mpf(H["s_f"])

eng = Engine(); print("engine:", eng.banner)
c = Cusp(eng, MU["a"], MU["a20"], side=1)

# ---------------------------------------------------------------- (1)
print("\n(1) D_s(1) = exp(pi T / omega) - 1 ?")
for dT in ("0.05", "0.01", "0.002"):
    mu = dict(MU); mu["a01"] += mp.mpf(dT)      # T is +1 in a01
    T, L = Tval(mu), detJ(mu)
    om = mp.sqrt(L - T*T/4)
    pred = mp.e**(mp.pi*T/om) - 1
    row = []
    for eps in ("1e-4", "1e-5", "1e-6"):
        r = FC.val(c, mu, 1 + mp.mpf(eps))
        row.append(r["Dx"] if r["status"] == "OK" else None)
    print("   T=%-10s omega=%-14s pred=%-16s Dx(1+e)=%s"
          % (mp.nstr(T,6), mp.nstr(om,8), mp.nstr(pred,10),
             [None if v is None else mp.nstr(v,10) for v in row]))
    if row[-1] is not None:
        print("        rel err at e=1e-6: %s" % mp.nstr(abs(row[-1]-pred)/abs(pred), 6))

# ---------------------------------------------------------------- (2)
print("\n(2) first Lyapunov coefficient at the organizer (T = %s)" % mp.nstr(Tval(MU),6))
print("    D(s) = c3 (s-1)^3 + ...   =>  D_sss(1) = 6 c3,  sign(l1) = sign(c3)")
for eps in ("1e-3", "1e-4", "1e-5", "1e-6"):
    r = FC.val(c, MU, 1 + mp.mpf(eps))
    if r["status"] != "OK": print("   eps=%s  %s" % (eps, r["status"])); continue
    e = mp.mpf(eps)
    print("   eps=%-6s D=%-14s D_s=%-14s D_ss=%-14s D_sss=%s   D/eps^3=%s"
          % (eps, mp.nstr(r["D"],7), mp.nstr(r["Dx"],7), mp.nstr(r["Dxx"],7),
             mp.nstr(r["Dxxx"],10), mp.nstr(r["D"]/e**3, 10)))

# sign convention test: perturb T off zero and see which side grows a cycle
print("\n    Hopf direction test: small |T|, look for a cycle near s=1")
for dT in ("1e-6", "-1e-6", "1e-5", "-1e-5"):
    mu = dict(MU); mu["a01"] += mp.mpf(dT)
    T = Tval(mu)
    found = None
    prev = None
    e = mp.mpf("1e-4")
    while e < mp.mpf("0.2"):
        r = FC.val(c, mu, 1+e)
        if r["status"] != "OK": break
        if prev is not None and prev[1]*r["D"] < 0:
            lo, hi = prev[0], e
            for _ in range(80):
                mid=(lo+hi)/2; rm=FC.val(c,mu,1+mid)
                if rm["D"]*prev[1] > 0: lo=mid
                else: hi=mid
            found = (lo+hi)/2; break
        prev = (e, r["D"]); e *= mp.mpf("1.3")
    print("   T=%-12s small cycle at s-1 = %s"
          % (mp.nstr(T,4), "none" if found is None else mp.nstr(found,8)))

# ---------------------------------------------------------------- (3)
print("\n(3) rank[ grad_mu T ; d_mu D(s_f) ]  (scaled coordinates)")
sc = FC.scales(MU)
gT_raw = {"a": mp.mpf(-2), "a20": mp.mpf(0), "a11": mp.mpf(1), "a01": mp.mpf(1), "a10": mp.mpf(0)}
gT = [gT_raw[p]*sc[j] for j, p in enumerate(FC.PARAMS)]
g = FC.grad_mu(c, MU, S_F, which=("D","Dx"))
bf = g["D"]
print("   grad_mu T   =", [mp.nstr(v,8) for v in gT])
print("   d_mu D(s_f) =", [mp.nstr(v,8) for v in bf])
def norm(v): return mp.sqrt(sum(x*x for x in v))
u1 = [x/norm(gT) for x in gT]; u2 = [x/norm(bf) for x in bf]
cs = sum(x*y for x, y in zip(u1, u2))
A = mp.matrix(2, 5)
for j in range(5): A[0,j] = u1[j]; A[1,j] = u2[j]
G = A*A.T
tr, det = G[0,0]+G[1,1], G[0,0]*G[1,1]-G[0,1]*G[1,0]
s1 = mp.sqrt((tr+mp.sqrt(tr*tr-4*det))/2); s2 = mp.sqrt(max(mp.mpf(0),(tr-mp.sqrt(tr*tr-4*det))/2))
print("   cos(angle) = %s   row-normalised sigma = (%s, %s)  ratio %s"
      % (mp.nstr(cs,8), mp.nstr(s1,8), mp.nstr(s2,8), mp.nstr(s2/s1,8)))
print("   => rank 2 (numerically) iff sigma_min is not at round-off; sigma_min = %s" % mp.nstr(s2,8))
print("\ncalls:", eng.ncalls); eng.close()
