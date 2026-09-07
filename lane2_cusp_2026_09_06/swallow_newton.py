#!/usr/bin/env python3
"""Resolve the sign-change event by SOLVING, not by sampling a ratio.

Square Newton system at fixed shape (a, a20):
    Phi(u, x0) = (D, D_x, D_xx, D_xxx) = 0,   u = (a11, a01, a10).
The x0-column of the Jacobian is (D_x, D_xx, D_xxx, D_xxxx), free from the
degree-4 jet; the u-columns are central differences.  At a quadruple point
det Jac = det(F_u) * D_xxxx, so a fifth-order section jet is NOT required.

Then report H = D_xxxx at the solution with an error estimate, and separately
solve H = 0 along F = 0 to locate that event independently."""
import json, sys
import mpmath as mp
from engine import Engine
from cusp import Cusp, solve3, solve4, FD_H
mp.mp.dps = 50

def jet(c, u, x0):
    r = c.val(list(u), x0)
    return (None, r) if r["status"] != "OK" else ([r["D"], r["Dx"], r["Dxx"], r["Dxxx"]], r)

def jac4(c, u, x0, r):
    """4x4 d(D,Dx,Dxx,Dxxx)/d(a11,a01,a10,x0)."""
    M = [[mp.mpf(0)]*4 for _ in range(4)]
    for j in range(3):
        h = FD_H[j]*max(mp.mpf(1), abs(u[j]))
        up, um = list(u), list(u); up[j] += h; um[j] -= h
        rp, rm = c.val(up, x0), c.val(um, x0)
        if rp["status"] != "OK" or rm["status"] != "OK": return None
        for i, k in enumerate(("D", "Dx", "Dxx", "Dxxx")):
            M[i][j] = (rp[k] - rm[k])/(2*h)
    for i, k in enumerate(("Dx", "Dxx", "Dxxx", "Dxxxx")):
        M[i][3] = r[k]
    return M

def wres4(F, x0):
    r = abs(mp.mpf(x0) - 1)
    return max(abs(F[0]), abs(F[1])*r, abs(F[2])*r*r/2, abs(F[3])*r**3/6)

def newton_swallow(c, u, x0, itmax=60, verbose=True):
    u = list(u); x0 = mp.mpf(x0); prev = mp.inf
    best = (mp.inf, list(u), x0, None)
    for it in range(itmax):
        F, r = jet(c, u, x0)
        if F is None: return None, r["status"], None
        res = wres4(F, x0)
        if res < best[0]: best = (res, list(u), x0, r)
        J = jac4(c, u, x0, r)
        if J is None: return None, "jac-fail", None
        s = solve4(J, F)
        if s is None: return None, "singular", None
        sc = max(abs(s[k])/max(mp.mpf(1), abs(([*u, x0])[k])) for k in range(4))
        if verbose and it < 14:
            print("   it%-2d res=%.3e  Dxxx=%+.6e  Dxxxx=%+.6e  step=%.2e"
                  % (it, float(res), float(F[3]), float(r["Dxxxx"]), float(sc)), flush=True)
        if sc >= prev/2 and it >= 4: break
        prev = sc
        u = [u[k] - s[k] for k in range(3)]; x0 = x0 - s[3]
        if sc < mp.mpf("1e-40"): break
    return best[1], best[2], best[3]

if __name__ == "__main__":
    recs = [json.loads(l) for l in open('ledger_grid/cusp_c_am2p0_om0p08.jsonl') if l.strip()]
    R5, R6 = recs[5], recs[6]
    a, a20 = R5["a"], R5["a20"]
    eng = Engine(); print("engine:", eng.banner)
    c = Cusp(eng, a, a20, side=int(R5.get("side", 1)))
    for seed_name, R in (("idx5", R5), ("idx6", R6)):
        print("\n=== Newton on (D,Dx,Dxx,Dxxx)=0 from %s ===" % seed_name, flush=True)
        u0 = [mp.mpf(R["a11"]), mp.mpf(R["a01"]), mp.mpf(R["a10"])]
        u, x0, r = newton_swallow(c, u0, mp.mpf(R["x0"]))
        if u is None: print("   FAILED:", x0); continue
        F = [r["D"], r["Dx"], r["Dxx"], r["Dxxx"]]
        print("   converged x0 = %s" % mp.nstr(x0, 25))
        print("   a11 = %s" % mp.nstr(u[0], 25))
        print("   a01 = %s" % mp.nstr(u[1], 25))
        print("   a10 = %s" % mp.nstr(u[2], 25))
        print("   residual (D,Dx,Dxx,Dxxx) = %s" % [mp.nstr(v, 6) for v in F])
        print("   weighted residual        = %s" % mp.nstr(wres4(F, x0), 6))
        print("   H = D_xxxx at the solution = %s" % mp.nstr(r["Dxxxx"], 20))
        print("   transversality = %s   T = %s" % (mp.nstr(r["transv"], 8), mp.nstr(r["T"], 10)))
    print("\ncalls:", eng.ncalls); eng.close()
