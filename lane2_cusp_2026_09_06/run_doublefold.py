#!/usr/bin/env python3
"""Drive both same-type stationary values of D to zero at Cherkas rows 8, 7, 2."""
import json
import mpmath as mp
from engine import Engine
from cusp import Cusp
from doublefold import DF, solve
import solver as SV
mp.mp.dps = 40

ROWS = {
 8: dict(a="1.04", a20=-120, a11="1.51997", a01="1.56", a10="-79.6", s1="1.20", s2="3.48"),
 7: dict(a=mp.mpf(8)/11, a20=-12, a11="2.1502", a01=mp.mpf(67)/220, a10="-26.5", s1="1.28", s2="3.72"),
 2: dict(a="1.5", a20=-15, a11="0.79993", a01="3.2", a10="9.17", s1="1.15", s2="2.72"),
}
eng = Engine(); print("engine:", eng.banner)
lg = open("ledger_opus/doublefold.jsonl", "a")
for rid, P in ROWS.items():
    for pair in (("a", "a11"), ("a", "a10"), ("a20", "a11"), ("a", "a20")):
        c = Cusp(eng, P["a"], P["a20"], side=1)
        df = DF(c, dict(a11=mp.mpf(P["a11"]), a01=mp.mpf(P["a01"]), a10=mp.mpf(P["a10"])))
        print("\n=== row %d, controls %s ===" % (rid, pair), flush=True)
        res = solve(df, P["s1"], P["s2"], pair, verbose=True)
        print("   -> %s (its=%d)" % (res.status, res.its))
        if res.converged:
            s1, s2, mu = res.root
            r1, r2 = df.val(s1), df.val(s2)
            print("   *** s1=%s s2=%s   D_ss(s1)=%s  D_ss(s2)=%s"
                  % (mp.nstr(s1,12), mp.nstr(s2,12), mp.nstr(r1["Dxx"],8), mp.nstr(r2["Dxx"],8)))
            lg.write(json.dumps(dict(row=rid, pair=list(pair), a=str(df.c.a), a20=str(df.c.a20),
                a11=mp.nstr(mu["a11"],25), a01=mp.nstr(mu["a01"],25), a10=mp.nstr(mu["a10"],25),
                s1=mp.nstr(s1,25), s2=mp.nstr(s2,25),
                Dxx1=mp.nstr(r1["Dxx"],20), Dxx2=mp.nstr(r2["Dxx"],20),
                status=res.status, calls=eng.ncalls))+"\n"); lg.flush()
print("\ncalls:", eng.ncalls); eng.close()
