#!/usr/bin/env python3
"""Wide scan of the graphic splitting over the order-3 stratum, restricted to
the region where infinity carries exactly one antipodal saddle pair."""
import numpy as np
from splitting import split
from sphere3 import inf_points

def ndir(l, a):
    return len(inf_points(l, a))//2

rows = []
for l in np.arange(-30.0, 6.01, 1.0):
    for a in np.arange(0.2, 4.01, 0.2):
        n = ndir(l, a)
        if n != 1:   # need the antipodal-pair-only case
            continue
        try:
            r = split(float(l), float(a))
        except Exception:
            r = None
        if r is None:
            rows.append((l, a, None)); continue
        rows.append((l, a, r[0]))

vals = [(l, a, s) for l, a, s in rows if s is not None]
print("scanned %d stratum points with a single antipodal saddle pair; %d gave a splitting"
      % (len(rows), len(vals)))
if vals:
    s = np.array([v[2] for v in vals])
    print("splitting range: [%.4f, %.4f]   sign changes: %s"
          % (s.min(), s.max(), "YES" if (s.min() < 0 < s.max()) else "no"))
    k = int(np.argmax(-np.abs(s)))
    print("closest to zero: l=%.2f a=%.2f  splitting=%+.6e" % vals[k])
    print("\nsample rows (l, a, splitting):")
    for v in vals[::max(1, len(vals)//25)]:
        print("   %7.2f %5.2f  %+.6f" % v)
none = [(l, a) for l, a, s in rows if s is None]
print("\npoints with no section crossing (separatrix never reaches the section): %d" % len(none))
if none: print("   e.g.", none[:8])
