#!/usr/bin/env python3
"""Complete the a20 = -1 block of the compatibility table, over the admissible
window a > 1/3 only (entry needs 1-3a of the right sign)."""
import json
import mpmath as mp, numpy as np
from engine import Engine, third_order
from cusp import Cusp
from remote_nest import equilibria
mp.mp.dps = 50
eng = Engine(); print("engine:", eng.banner, flush=True)
lg = open("ledger_opus/compatibility.jsonl", "a")
r0 = mp.mpf("0.05"); a20 = -1
print("===== a20 = %g, r0 = %s =====" % (a20, r0))
print("   a        nu           second focus?   other equilibria", flush=True)
for aval in (2.0, 1.5, 1.04, 0.727):
    c = Cusp(eng, aval, a20, side=1)
    mu, r = c.newton_mu(list(third_order(aval, a20)), 1 + r0, verbose=False)
    if mu is None:
        print("   %-8g ENTRY FAILED (%s)" % (aval, r), flush=True); continue
    nu = r["Dxxx"]/(r["Dxxxx"]*r0)
    a11, a01, a10 = [float(v) for v in mu]
    others = [e for e in equilibria(float(aval), float(a20), a11, a01, a10)
              if abs(e[0]-1.0) > 1e-6]
    foci = [e for e in others if e[4] == "focus"]
    print("   %-8g %-12.6g %-15s %s"
          % (aval, float(nu), "YES" if foci else "no", [e[4] for e in others]), flush=True)
    lg.write(json.dumps(dict(a=str(aval), a20=str(a20), r0=str(r0),
                             nu=mp.nstr(nu, 20), coefnorm=float(np.sqrt(a11**2+a01**2+a10**2)),
                             I1=a11**2/float(a20), nfoci=len(foci),
                             kinds=[e[4] for e in others], calls=eng.ncalls))+"\n")
    lg.flush()
print("calls:", eng.ncalls, flush=True); eng.close()
