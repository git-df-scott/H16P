#!/usr/bin/env python3
"""Is the splitting zero at a_deg a genuine transversal crossing, or an artefact
(a discontinuity, a branch flip, or a section that stops being transversal)?
Sample finely across a_deg at l=-10 (a_deg=2) and report the splitting, the two
saddle directions, and which branch was selected."""
import numpy as np
from locate import splitting
from splitting_gen import params, split
from sphere_gen import inf_points

lv = -10.0; adeg = np.sqrt(-(lv+2)/2)
print("l=%g  a_deg=%.12f" % (lv, adeg))
print("   a          m           b            u1          u2        ratio1     ratio2     splitting")
prev = None
for av in np.linspace(adeg-0.10, adeg+0.10, 21):
    mv, bv = params(lv, av)
    pts = inf_points(lv, mv, av, bv)
    sad = sorted(set(round(u, 10) for (u, sv, k, le, lt) in pts if k == "saddle"))
    rr = []
    for u in sad:
        le = -3*u*u - 2*mv*u + (bv-lv); lt = -(u*u + mv*u + lv)
        rr.append(abs(le/lt))
    s = splitting(lv, av)
    print("   %-10.6f %-11.6f %-12.6f %-11.6f %-9.6f %-10.6f %-10.6f %s"
          % (av, mv, bv, sad[0] if sad else np.nan, sad[1] if len(sad) > 1 else np.nan,
             rr[0] if rr else np.nan, rr[1] if len(rr) > 1 else np.nan,
             "n/a" if s is None else "%+.6e" % s), flush=True)
