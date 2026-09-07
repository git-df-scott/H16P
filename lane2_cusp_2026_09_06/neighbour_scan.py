#!/usr/bin/env python3
"""At shapes NEAR the centre event, solve the well-conditioned cusp system
F=0 in u at a grid of x0 (newton_mu), and tabulate G=D_xxx and H=D_xxxx along
the resulting cusp curve.  A swallow-tail needs G=0 with H!=0.

This avoids the singular square Newton at the centre itself."""
import json
import mpmath as mp
from engine import Engine
from cusp import Cusp, wres
mp.mp.dps = 50

recs = [json.loads(l) for l in open('ledger_grid/cusp_c_am2p0_om0p08.jsonl') if l.strip()]
R = recs[5]; a = R["a"]; a20_0 = mp.mpf(R["a20"])
eng = Engine(); print("engine:", eng.banner)
lg = open("ledger_opus/neighbour.jsonl", "a")
for off in ("-0.02", "0.02"):
    a20 = a20_0 + mp.mpf(off)
    c = Cusp(eng, a, a20, side=1)
    print("\n=== a20 = %s (offset %s) ===" % (mp.nstr(a20, 16), off))
    print("   x0        res         G=D_xxx        H=D_xxxx       nu           V1")
    u = [mp.mpf(R["a11"]), mp.mpf(R["a01"]), mp.mpf(R["a10"])]
    prevG = None
    for k in range(11):
        x0 = mp.mpf("1.020") + mp.mpf(k)*mp.mpf("0.004")
        mu, r = c.newton_mu(u, x0, verbose=False)
        if mu is None:
            print("   %-9s FAILED (%s)" % (mp.nstr(x0, 8), r)); continue
        u = mu
        G, H = r["Dxxx"], r["Dxxxx"]; r0 = x0 - 1
        V1 = mu[0] + mu[1] - 2*mp.mpf(a) - 1
        res = wres([r["D"], r["Dx"], r["Dxx"]], x0)
        flag = ""
        if prevG is not None and prevG*G < 0: flag = "   <== G SIGN CHANGE"
        prevG = G
        print("   %-9s %-11s %-14s %-14s %-12s %s%s"
              % (mp.nstr(x0, 8), mp.nstr(res, 4), mp.nstr(G, 8), mp.nstr(H, 8),
                 mp.nstr(G/(H*r0), 8), mp.nstr(V1, 6), flag))
        lg.write(json.dumps(dict(a=str(a), a20=mp.nstr(a20,25), x0=mp.nstr(x0,25),
                                 a11=mp.nstr(mu[0],25), a01=mp.nstr(mu[1],25),
                                 a10=mp.nstr(mu[2],25), G=mp.nstr(G,20), H=mp.nstr(H,20),
                                 res=mp.nstr(res,10), V1=mp.nstr(V1,20),
                                 calls=eng.ncalls))+"\n")
        lg.flush()
print("\ncalls:", eng.ncalls); eng.close()
