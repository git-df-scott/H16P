#!/usr/bin/env python3
"""THE compatibility question: can one shape supply BOTH a small nu (a route to
the fourth cycle in the nest) AND a second focus (the fifth cycle)?

Sweep the shape parameter a at fixed a20, entering the cusp locus at fixed r0,
and record nu together with the type of every other finite equilibrium."""
import json, sys
import mpmath as mp, numpy as np
from engine import Engine, third_order, L_of, V1_of
from cusp import Cusp
from remote_nest import equilibria
mp.mp.dps = 50

eng = Engine(); print("engine:", eng.banner)
lg = open("ledger_opus/compatibility.jsonl", "a")
r0 = mp.mpf("0.05")
for a20 in (-12, -120, -1):
    print("\n===== a20 = %g,  r0 = %s =====" % (a20, r0))
    print("   a       |coef|   nu           I1=a11^2/a20   second focus?   other equilibria")
    for aval in (2.0, 1.5, 1.04, 0.727, 0.3, 0.0, -0.5, -1.0, -2.0, -3.0, -4.0):
        c = Cusp(eng, aval, a20, side=1)
        mu0 = list(third_order(aval, a20))
        mu, r = c.newton_mu(mu0, 1 + r0, verbose=False)
        if mu is None:
            print("   %-7g  ENTRY FAILED (%s)" % (aval, r)); continue
        G, H = r["Dxxx"], r["Dxxxx"]
        nu = G/(H*r0)
        a11, a01, a10 = [float(v) for v in mu]
        eqs = equilibria(float(aval), float(a20), a11, a01, a10)
        others = [e for e in eqs if abs(e[0]-1.0) > 1e-6]
        foci = [e for e in others if e[4] == "focus"]
        nrm = float(np.sqrt(a11**2 + a01**2 + a10**2))
        I1 = a11**2/float(a20)
        desc = ", ".join("%s" % e[4] for e in others) or "none"
        print("   %-7g %-8.2f %-12.6g %-14.5g %-15s %s"
              % (aval, nrm, float(nu), I1, "YES (%d)" % len(foci) if foci else "no", desc))
        lg.write(json.dumps(dict(a=str(aval), a20=str(a20), r0=str(r0),
                                 nu=mp.nstr(nu, 20), coefnorm=nrm, I1=I1,
                                 nfoci=len(foci), kinds=[e[4] for e in others],
                                 calls=eng.ncalls))+"\n")
        lg.flush()
print("\ntotal engine calls:", eng.ncalls); eng.close()
