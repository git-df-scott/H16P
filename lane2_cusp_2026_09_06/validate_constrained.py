#!/usr/bin/env python3
"""Validate the constrained-derivative machinery.

dG/dx0 (G = D_xxx) uses only quantities the jet supplies, so it must reproduce
the slope of D_xxx measured along the saved cusp curve by finite differences.
dnu/dx0 is NOT validated: it needs D_xxxxx, which the degree-4 jet lacks, so
that component is reported as unavailable rather than as zero."""
import json, sys
import mpmath as mp
from engine import Engine
from cusp import Cusp
from constrained import constrained_derivs
mp.mp.dps = 50

src = sys.argv[1] if len(sys.argv) > 1 else "ledger/cusp_row5.jsonl"
recs = [json.loads(l) for l in open(src) if l.strip()]
eng = Engine(); print("engine:", eng.banner)
print("\n idx   x0         dG/dx0 (constrained)   dG/dx0 (ledger FD)     rel.diff")
for i in (5, 40, 80, 110):
    R = recs[i]
    c = Cusp(eng, R["a"], R["a20"], side=int(R.get("side", 1)))
    mu = [mp.mpf(R["a11"]), mp.mpf(R["a01"]), mp.mpf(R["a10"])]
    x0 = mp.mpf(R["x0"])
    out = constrained_derivs(eng, c, mu, x0)
    if out is None: print("  %-5d FAILED" % i); continue
    g0 = mp.mpf(recs[i-1]["Dxxx"]); g1 = mp.mpf(recs[i+1]["Dxxx"])
    t0 = mp.mpf(recs[i-1]["x0"]);   t1 = mp.mpf(recs[i+1]["x0"])
    fd = (g1 - g0)/(t1 - t0)
    cd = out["dG"][2]
    print("  %-5d %-10.6f %-22.10g %-22.10g %.3e"
          % (i, float(x0), float(cd), float(fd), float(abs(cd-fd)/max(abs(fd), mp.mpf('1e-300')))))
print("\ncalls:", eng.ncalls); eng.close()
