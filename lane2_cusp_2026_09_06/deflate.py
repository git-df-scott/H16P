#!/usr/bin/env python3
"""Is the solution branch of (F,G)=0 near this event ONLY the centre variety, or
is there a separate non-centre solution with H != 0 (the counterexample case)?

Solve the square system at a range of nearby shapes a20 and report H and the
size of D away from the solution point (a centre has D ~ 0 everywhere)."""
import json
import mpmath as mp
from engine import Engine
from cusp import Cusp
from swallow_newton import newton_swallow, wres4
mp.mp.dps = 50

recs = [json.loads(l) for l in open('ledger_grid/cusp_c_am2p0_om0p08.jsonl') if l.strip()]
R = recs[5]
a = R["a"]; a20_0 = mp.mpf(R["a20"])
eng = Engine(); print("engine:", eng.banner)
lg = open("ledger_opus/deflate.jsonl", "a")
print("\n  a20 offset   a20              x0 (solved)      H=D_xxxx at soln   |D| at x0+0.05   verdict")
for off in ("0", "-0.02", "0.02", "-0.10", "0.10", "-0.30", "0.30"):
    a20 = a20_0 + mp.mpf(off)
    c = Cusp(eng, a, a20, side=1)
    u0 = [mp.mpf(R["a11"]), mp.mpf(R["a01"]), mp.mpf(R["a10"])]
    u, x0, r = newton_swallow(c, u0, mp.mpf(R["x0"]), verbose=False)
    if u is None:
        print("  %-12s %-16s FAILED (%s)" % (off, mp.nstr(a20, 12), x0)); continue
    H = r["Dxxxx"]
    r2 = c.val(u, x0 + mp.mpf("0.05"))
    Daway = r2["D"] if r2["status"] == "OK" else None
    centre = (abs(H) < mp.mpf("1e-12")) and (Daway is not None and abs(Daway) < mp.mpf("1e-15"))
    print("  %-12s %-16s %-16s %-18s %-16s %s"
          % (off, mp.nstr(a20, 12), mp.nstr(x0, 12), mp.nstr(H, 6),
             "n/a" if Daway is None else mp.nstr(abs(Daway), 6),
             "CENTRE" if centre else "*** NON-CENTRE: H != 0 ***"), flush=True)
    lg.write(json.dumps(dict(a=str(a), a20=mp.nstr(a20,25), x0=mp.nstr(x0,25),
                             a11=mp.nstr(u[0],25), a01=mp.nstr(u[1],25), a10=mp.nstr(u[2],25),
                             H=mp.nstr(H,20), D_away=None if Daway is None else mp.nstr(Daway,20),
                             centre=bool(centre), calls=eng.ncalls))+"\n")
    lg.flush()
print("\ncalls:", eng.ncalls); eng.close()
