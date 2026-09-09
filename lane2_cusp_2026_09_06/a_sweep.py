#!/usr/bin/env python3
"""Test the identified descent direction: does decreasing the shape parameter a
drive nu toward zero at FIXED, moderate amplitude and BOUNDED coefficients?

At each a we pick a20 far from the centre curve a20_c(a) (so the jet does not
degenerate), enter the cusp locus at fixed r0, and record nu, the coefficient
norm, L (antisaddle) and V1."""
import json, sys
import mpmath as mp
from engine import Engine, third_order, V7_of, L_of, V1_of
from cusp import Cusp, wres
mp.mp.dps = 50

def a20c(a):
    a = mp.mpf(a); den = (a-1)*(2*a+1)**2
    return None if den == 0 else 4*a*(a+1)*(a-2)**2/den

eng = Engine(); print("engine:", eng.banner)
lg = open("ledger_opus/a_sweep.jsonl", "a")
r0 = mp.mpf("0.05")
print("\n  a       a20      a20_c     dist    |coef|   nu            D_xxx         D_xxxx        L        V1")
for aval in (-2, -3, -4, -6, -8, -12, -20, -40):
    ac = a20c(aval)
    a20 = mp.mpf(-1)                     # far from a20_c for all these a (a20_c < -4.7)
    c = Cusp(eng, aval, a20, side=1)
    mu0 = list(third_order(aval, a20))
    x0 = 1 + r0
    mu, r = c.newton_mu(mu0, x0, verbose=False)
    if mu is None:
        print("  %-7g %-8.4g %-9.4g %-7.3g  ENTRY FAILED: %s" % (aval, float(a20), float(ac), float(a20-ac), r))
        continue
    G, H = r["Dxxx"], r["Dxxxx"]
    nu = G/(H*r0)
    nrm = float(mp.sqrt(mu[0]**2+mu[1]**2+mu[2]**2))
    L = L_of(aval, a20, mu[0], mu[1], mu[2]); V1 = V1_of(aval, mu[0], mu[1])
    print("  %-7g %-8.4g %-9.4g %-7.3g %-8.2f %-13.6g %-13.5g %-13.5g %-8.4g %.3g"
          % (aval, float(a20), float(ac), float(a20-ac), nrm, float(nu),
             float(G), float(H), float(L), float(V1)))
    lg.write(json.dumps(dict(a=str(aval), a20=str(a20), r0=str(r0),
                             nu=mp.nstr(nu,20), Dxxx=mp.nstr(G,20), Dxxxx=mp.nstr(H,20),
                             coefnorm=nrm, L=mp.nstr(L,20), V1=mp.nstr(V1,20),
                             calls=eng.ncalls))+"\n")
    lg.flush()
print("\ncalls:", eng.ncalls); eng.close()
