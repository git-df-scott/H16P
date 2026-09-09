#!/usr/bin/env python3
"""Shape gradient of the swallow-tail indicator on the cusp locus.

At saved cusp points, compute (dnu/da, dnu/da20) constrained to F=0, together
with the distance to the centre curve a20_c(a) (where the whole jet degenerates
and nu stays finite -- the trap) and the coefficient norm (the other trap).

A finite-amplitude swallow-tail needs nu -> 0 with the shape staying away from
a20_c and the coefficients bounded."""
import json, glob, os, sys
import mpmath as mp
from engine import Engine
from cusp import Cusp
from constrained import constrained_derivs
mp.mp.dps = 50

def a20c(a):
    a = mp.mpf(a); den = (a-1)*(2*a+1)**2
    return None if den == 0 else 4*a*(a+1)*(a-2)**2/den

eng = Engine(); print("engine:", eng.banner)
lg = open("ledger_opus/shape_gradient.jsonl", "a")
print("\n ledger                      idx  a      a20      x0        |coef|  dist to a20_c   nu         dnu/da       dnu/da20    steepest-descent dir")
files = ['ledger/cusp_row5.jsonl', 'ledger/cusp_row6.jsonl', 'ledger/cusp_row7.jsonl',
         'ledger/cusp_row8.jsonl', 'ledger_grid/cusp_c_am2p5_o0p08.jsonl',
         'ledger_grid/cusp_c_am2p9_o1p0.jsonl', 'ledger_grid/cusp_c_am2p0_o0p08.jsonl']
for f in files:
    recs = [json.loads(l) for l in open(f) if l.strip()]
    n = len(recs)
    for i in sorted(set([n//4, n//2, (3*n)//4])):
        if i >= n: continue
        R = recs[i]
        c = Cusp(eng, R["a"], R["a20"], side=int(R.get("side", 1)))
        mu = [mp.mpf(R["a11"]), mp.mpf(R["a01"]), mp.mpf(R["a10"])]
        x0 = mp.mpf(R["x0"])
        out = constrained_derivs(eng, c, mu, x0)
        if out is None:
            print("  %-27s %-4d FAILED" % (os.path.basename(f), i)); continue
        a = mp.mpf(R["a"]); a20 = mp.mpf(R["a20"])
        ac = a20c(a); dist = float(a20 - ac) if ac is not None else float('nan')
        nrm = float(mp.sqrt(mu[0]**2 + mu[1]**2 + mu[2]**2))
        g1, g2 = out["dnu"][0], out["dnu"][1]
        # steepest descent of nu in shape space
        nn = mp.sqrt(g1*g1 + g2*g2)
        d = (-g1/nn, -g2/nn) if nn > 0 else (0, 0)
        print("  %-27s %-4d %-6.3g %-8.4g %-9.5f %-7.2f %-15.4g %-10.5g %-12.4g %-11.4g (%+.3f,%+.3f)"
              % (os.path.basename(f), i, float(a), float(a20), float(x0), nrm, dist,
                 float(out["nu"]), float(g1), float(g2), float(d[0]), float(d[1])))
        lg.write(json.dumps(dict(f=os.path.basename(f), idx=i, a=str(a), a20=str(a20),
                                 x0=str(x0), nu=mp.nstr(out["nu"], 20),
                                 dnu_da=mp.nstr(g1, 20), dnu_da20=mp.nstr(g2, 20),
                                 dist_a20c=dist, coefnorm=nrm, calls=eng.ncalls))+"\n")
        lg.flush()
print("\ntotal engine calls:", eng.ncalls); eng.close()
