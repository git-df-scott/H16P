#!/usr/bin/env python3
"""Sign map of the graphic splitting over the order-2 + neutral family.
The locus where the family meets the order-3 stratum (eta_2 = 0) is
    m = 5a  <=>  a^2 = -(l+2)/2 .
Question: does the splitting vanish anywhere OFF that locus?"""
import numpy as np, sys
from locate import splitting
from splitting_gen import params

print("  l        a      a_deg(eta2=0)   m        b          splitting")
rows = []
for lv in (-20.0, -14.0, -10.0, -7.0, -5.0, -4.0, -3.0, -2.6):
    adeg = np.sqrt(-(lv+2)/2) if -(lv+2)/2 > 0 else np.nan
    for av in (0.2, 0.5, 0.8, 1.2, 1.8, 2.4, 3.0, 4.0):
        try:
            s = splitting(lv, av)
        except Exception:
            s = None
        mv, bv = params(lv, av)
        rows.append((lv, av, adeg, s))
        print("  %-8.2f %-6.2f %-14.5f %-8.4g %-10.4g %s"
              % (lv, av, adeg, mv, bv, "n/a" if s is None else "%+.4e" % s), flush=True)
vals = [(l, a, ad, s) for (l, a, ad, s) in rows if s is not None]
print("\nsign summary per l line (a increasing), with the eta_2=0 crossing marked:")
for lv in sorted(set(v[0] for v in vals)):
    line = [(a, s, ad) for (l, a, ad, s) in vals if l == lv]
    line.sort()
    ad = line[0][2]
    txt = " ".join(("+" if s > 0 else "-") for a, s, _ in line)
    ax = " ".join(("%.1f" % a) for a, s, _ in line)
    print("  l=%-7.2f a: %s" % (lv, ax))
    print("          sgn: %s     eta_2=0 at a=%.4f" % (txt, ad))
