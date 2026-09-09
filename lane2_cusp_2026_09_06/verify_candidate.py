#!/usr/bin/env python3
"""Are the H != 0 solutions genuine solutions of (D,D_x,D_xx,D_xxx)=0, or
non-converged Newton iterates?  newton_swallow returns its LOWEST-RESIDUAL
iterate, which is not the same as a converged root."""
import json
import mpmath as mp
from engine import Engine
from cusp import Cusp
from swallow_newton import newton_swallow, wres4, jac4
mp.mp.dps = 50

recs = [json.loads(l) for l in open('ledger_grid/cusp_c_am2p0_om0p08.jsonl') if l.strip()]
R = recs[5]; a = R["a"]; a20_0 = mp.mpf(R["a20"])
eng = Engine(); print("engine:", eng.banner)
for off in ("-0.02", "0.02"):
    a20 = a20_0 + mp.mpf(off)
    c = Cusp(eng, a, a20, side=1)
    u, x0, r = newton_swallow(c, [mp.mpf(R["a11"]), mp.mpf(R["a01"]), mp.mpf(R["a10"])],
                              mp.mpf(R["x0"]), verbose=False)
    F = [r["D"], r["Dx"], r["Dxx"], r["Dxxx"]]
    W = wres4(F, x0); r0 = x0 - 1
    print("\n=== a20 offset %s  (a20 = %s) ===" % (off, mp.nstr(a20, 16)))
    print("   x0        = %s   r0 = %s" % (mp.nstr(x0, 20), mp.nstr(r0, 12)))
    print("   D         = %s" % mp.nstr(F[0], 8))
    print("   D_x       = %s" % mp.nstr(F[1], 8))
    print("   D_xx      = %s" % mp.nstr(F[2], 8))
    print("   D_xxx     = %s   <- must be 0" % mp.nstr(F[3], 8))
    print("   D_xxxx    = %s" % mp.nstr(r["Dxxxx"], 8))
    print("   weighted residual = %s" % mp.nstr(W, 8))
    # scale check: compare each residual with the corresponding Taylor term at r0
    print("   Taylor-term sizes at r0: D_xxx r0^3/6 = %s ,  D_xxxx r0^4/24 = %s"
          % (mp.nstr(abs(F[3])*r0**3/6, 6), mp.nstr(abs(r["Dxxxx"])*r0**4/24, 6)))
    print("   |D_xxx| relative to |D_xxxx| * r0 = %s   (nu)"
          % mp.nstr(F[3]/(r["Dxxxx"]*r0), 8))
    print("   transversality = %s   T = %s   status = %s"
          % (mp.nstr(r["transv"], 8), mp.nstr(r["T"], 10), r["status"]))
    J = jac4(c, u, x0, r)
    if J is not None:
        M = mp.matrix(4, 4)
        for i in range(4):
            for j in range(4): M[i, j] = J[i][j]
        print("   det Jac(F,G)/(u,x0) = %s" % mp.nstr(mp.det(M), 8))
    print("   CONVERGED?  %s" % ("YES" if W < mp.mpf("1e-24") else
                                 "NO -- this is a stalled iterate, not a root"))
print("\ncalls:", eng.ncalls); eng.close()
