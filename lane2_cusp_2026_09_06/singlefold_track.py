#!/usr/bin/env python3
"""Instructed step 4: continue a known finite-amplitude fold, and on its
constraint surface track the OTHER stationary point of D and its displacement.

Solve the single fold  (D(s2), D_s(s2)) = 0  in (p, s2) -- a well-conditioned
2x2 system -- then locate the other stationary point s1 (D_s(s1)=0) and report
D(s1).  A second separated double cycle needs D(s1) -> 0 with s1 at finite
amplitude."""
import json
import mpmath as mp
from engine import Engine
from cusp import Cusp, solve3
import solver as SV
mp.mp.dps = 40

ROWS = {
 7: dict(a=mp.mpf(8)/11, a20=-12, a11="2.1502", a01=mp.mpf(67)/220, a10="-26.5", s2="3.72"),
 8: dict(a="1.04", a20=-120, a11="1.51997", a01="1.56", a10="-79.6", s2="3.48"),
}
def solve2(A, b):
    d = A[0][0]*A[1][1] - A[0][1]*A[1][0]
    if d == 0: return None
    return [(b[0]*A[1][1] - b[1]*A[0][1])/d, (A[0][0]*b[1] - A[1][0]*b[0])/d]

eng = Engine(); print("engine:", eng.banner)
for rid, P in ROWS.items():
    c = Cusp(eng, P["a"], P["a20"], side=1)
    u = [mp.mpf(P["a11"]), mp.mpf(P["a01"]), mp.mpf(P["a10"])]
    s2 = mp.mpf(P["s2"])
    print("\n=== row %d : single outer fold in (a11, s2) ===" % rid)
    ok = False
    for it in range(40):
        r = c.val(u, s2)
        if r["status"] != "OK": print("   return failed:", r["status"]); break
        F = [r["D"], r["Dx"]]
        m = max(abs(F[0]), abs(F[1])*abs(s2-1))
        if m < mp.mpf("1e-26"): ok = True; break
        h = mp.mpf("1e-9")*max(mp.mpf(1), abs(u[0]))
        up, um = list(u), list(u); up[0] += h; um[0] -= h
        rp, rm = c.val(up, s2), c.val(um, s2)
        if rp["status"] != "OK" or rm["status"] != "OK": print("   fd failed"); break
        J = [[(rp["D"]-rm["D"])/(2*h),  r["Dx"]],
             [(rp["Dx"]-rm["Dx"])/(2*h), r["Dxx"]]]
        st = solve2(J, F)
        if st is None: print("   singular"); break
        u[0] -= st[0]; s2 -= st[1]
    print("   fold solve: %s  (a11=%s, s2=%s)" % ("CONVERGED" if ok else "FAILED",
          mp.nstr(u[0], 16), mp.nstr(s2, 12)))
    if not ok: continue
    rf = c.val(u, s2)
    print("   at the fold: D=%s  D_s=%s  D_ss=%s"
          % (mp.nstr(rf["D"],4), mp.nstr(rf["Dx"],4), mp.nstr(rf["Dxx"],8)))
    # census the remaining stationary points on this constraint surface
    print("   other stationary points of D at this parameter:")
    xs, Dx = [], []
    x = mp.mpf("1.02")
    while x < s2 - mp.mpf("0.1"):
        rr = c.val(u, x)
        if rr["status"] == "OK": xs.append(x); Dx.append(rr["Dx"])
        x += mp.mpf("0.03")
    found = False
    for i in range(len(xs)-1):
        if Dx[i]*Dx[i+1] < 0:
            lo, hi = xs[i], xs[i+1]
            for _ in range(60):
                mid = (lo+hi)/2; rm2 = c.val(u, mid)
                if rm2["status"] != "OK": break
                if rm2["Dx"]*Dx[i] > 0: lo = mid
                else: hi = mid
            s1 = (lo+hi)/2; r1 = c.val(u, s1)
            print("      s1 = %-12s  D(s1) = %-16s D_ss(s1) = %s"
                  % (mp.nstr(s1,10), mp.nstr(r1["D"],8), mp.nstr(r1["Dxx"],8)))
            found = True
    if not found: print("      none found below the fold")
print("\ncalls:", eng.ncalls); eng.close()
