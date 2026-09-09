#!/usr/bin/env python3
"""Regression cases for the solver contract."""
import json
import mpmath as mp
from engine import Engine
from cusp import Cusp
import solver as SV
mp.mp.dps = 50

eng = Engine(); print("engine:", eng.banner)
recs = [json.loads(l) for l in open('ledger_grid/cusp_c_am2p0_om0p08.jsonl') if l.strip()]
R = recs[5]; a = R["a"]; a20_0 = mp.mpf(R["a20"])
u0 = [mp.mpf(R["a11"]), mp.mpf(R["a01"]), mp.mpf(R["a10"])]

print("\nCASE 1 -- the exact centre (k=sqrt(1627)/5). Solves the vanishing-displacement")
print("          equations, but must FAIL the isolated-cycle nondegeneracy test.")
res = SV.newton_swallow(Cusp(eng, a, a20_0, side=1), u0, mp.mpf(R["x0"]))
print("   status = %s   scaled residuals = %s"
      % (res.status, [mp.nstr(v, 3) for v in (res.scaled or [])]))
print("   det Jac = %s" % (mp.nstr(res.jac_det, 6) if res.jac_det is not None else "n/a"))
ok1 = res.converged and res.jac_det is not None and abs(res.jac_det) < mp.mpf("1e-8")
print("   isolated-cycle nondegeneracy (det Jac != 0)?  %s" % ("FAILS (centre)" if ok1 else "??"))

print("\nCASE 2 -- the previously mis-reported offsets (Dxxx = -45.9, +18.5).")
print("          Must report STALLED, and expose the large unresolved component.")
for off in ("-0.02", "0.02"):
    r2 = SV.newton_swallow(Cusp(eng, a, a20_0 + mp.mpf(off), side=1), u0, mp.mpf(R["x0"]))
    big = max(range(len(r2.scaled)), key=lambda i: r2.scaled[i]) if r2.scaled else None
    print("   a20 offset %-6s status = %-9s root = %s   worst component = Phi_%s (%s)"
          % (off, r2.status, "None" if r2.root is None else "SET",
             big, mp.nstr(r2.scaled[big], 4) if r2.scaled else "n/a"))
    assert r2.root is None, "STALLED must not expose a root"
print("\n   contract held: no non-converged result exposes a root.")
print("\ncalls:", eng.ncalls); eng.close()
