#!/usr/bin/env python3
"""Claude tiny control: landscape of the necessary condition Phi(tau1)>0
over a coarse set of lobe-region anchor triples (all with tau1>5/11) and
lift parameters. Phi(tau1) is independent of the S1 tuning."""
import sys, os
sys.path.insert(0, os.path.dirname(__file__))
import mpmath as mp
from claude_green_tools import *
mp.mp.dps = 30
triples = [(0.46,0.47,0.48),(0.46,0.7,0.9),(0.46,0.99,0.999),(0.6,0.8,0.95),(0.7,0.75,0.8),(0.9,0.95,0.99),
           (0.5,0.999,0.9999),(0.99,0.999,0.9999),(0.999,0.9995,0.9999),(0.95,0.999,0.99999)]
ks = (1.2, 2, 4, 8.5, 20, 100, 1000)
best = (-9, None)
for y in triples:
    co = from_primitive_anchors_closed(tuple(map(mp.mpf, map(str, y))))
    out = []
    for k in ks:
        lift = Lift(k)
        P, Phi, Y0, P0 = P_Phi_at(lift, co, mp.mpf(str(y[0])))
        ratio = Phi/abs(Y0); out.append(mp.nstr(ratio, 6))
        if ratio > best[0]: best = (ratio, (y, k))
    print(f"anchors={y}: eta={mp.nstr(co[2],7)} Y0={mp.nstr(Y0,6)}  Phi(tau1)/|Y0| for k={ks}: {out}")
print("BEST Phi(tau1)/|Y0| =", mp.nstr(best[0], 8), "at", best[1])
