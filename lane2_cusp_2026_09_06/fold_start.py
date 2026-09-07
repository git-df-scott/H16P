#!/usr/bin/env python3
"""Task 2: recover the row-7 and row-8 outer folds with FULL coefficients,
validate on the second independent engine, and take the initial cycle
inventory in BOTH nests."""
import json
import mpmath as mp
from engine import Engine
from cusp import Cusp
import foldcont as FC
import indep_engine as IE
mp.mp.dps = 40

STARTS = {
 7: dict(a=mp.mpf(8)/11, a20=mp.mpf(-12), a11=mp.mpf("2.149952807637874"),
         a01=mp.mpf(67)/220, a10=mp.mpf("-26.5"), s="3.46899418303"),
 8: dict(a=mp.mpf("1.04"), a20=mp.mpf(-120), a11=mp.mpf("1.520003337556661"),
         a01=mp.mpf("1.56"), a10=mp.mpf("-79.6"), s="3.2755755769"),
}
eng = Engine(); print("engine A:", eng.banner)
lg = open("ledger_opus/fold_start.jsonl", "a")
for rid, P in STARTS.items():
    mu = {k: P[k] for k in FC.PARAMS}
    c = Cusp(eng, mu["a"], mu["a20"], side=1)
    # re-converge the fold to the solver contract
    mu2, s, r = FC.correct(c, mu, mp.mpf(P["s"]), pidx=2)
    if mu2 is None:
        print("row %d: fold corrector failed (%s)" % (rid, r)); continue
    print("\n=== row %d outer fold (CONVERGED) ===" % rid)
    for p in FC.PARAMS: print("   %-4s = %s" % (p, mp.nstr(mu2[p], 25)))
    a00 = mu2["a01"]+mu2["a11"]-mu2["a10"]-mu2["a20"]-mu2["a"]
    print("   a00  = %s   (pinned)" % mp.nstr(a00, 25))
    print("   s_f  = %s" % mp.nstr(s, 20))
    print("   D=%s  D_s=%s  D_ss=%s" % (mp.nstr(r["D"],4), mp.nstr(r["Dx"],4), mp.nstr(r["Dxx"],10)))
    print("   return time T=%s   transversality=%s" % (mp.nstr(r["T"],10), mp.nstr(r["transv"],8)))
    # independent engine at the SAME point
    di = IE.D(mu2["a"], mu2["a20"], mu2["a11"], mu2["a01"], mu2["a10"], s)
    print("   engine B (independent) D at s_f = %s"
          % (mp.nstr(di[0], 8) if di[0] is not None else "FAILED %s" % di[1]))
    inv = FC.inventory(c, mu2)
    print("   inventory on [%s, %s]:" % (mp.nstr(inv['domain'][0],5), mp.nstr(inv['domain'][1],5)))
    print("      cycles (D=0)      : %s" % [mp.nstr(v, 8) for v in inv["cycles"]])
    print("      stationary (D_s=0): %s" % [mp.nstr(v, 8) for v in inv["stationary"]])
    rem = FC.remote_cycles(mu2)
    print("   remote nest: second focus %s, cycles found = %s"
          % ("none" if rem["second_focus"] is None
             else "(%.4f, %.4f) trace %.4f" % rem["second_focus"], rem["cycles"]))
    T = mu2["a11"] + mu2["a01"] - 2*mu2["a"] - 1
    print("   trace at A: T = %s" % mp.nstr(T, 10))
    lg.write(json.dumps(dict(row=rid, **{p: mp.nstr(mu2[p],25) for p in FC.PARAMS},
        a00=mp.nstr(a00,25), s_f=mp.nstr(s,25), Dss=mp.nstr(r["Dxx"],20),
        T=mp.nstr(T,20), cycles=[mp.nstr(v,15) for v in inv["cycles"]],
        stationary=[mp.nstr(v,15) for v in inv["stationary"]],
        remote_cycles=rem["cycles"], calls=eng.ncalls))+"\n"); lg.flush()
print("\ncalls:", eng.ncalls); eng.close()
