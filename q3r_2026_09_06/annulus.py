#!/usr/bin/env python3
"""Locate the period annulus of the upper centre and the four generating
functions across it."""
import numpy as np, mpmath as mp
from melnikov import annulus_h, turning, basis, coeffs
mp.mp.dps = 30

def scan(a, b, n=60, fmax=200.0):
    hc = annulus_h(a, b)
    fs = np.geomspace(1.0005, fmax, n)
    rows = []
    for f in fs:
        h = hc*f
        try:
            t = turning(a, b, h)
        except Exception:
            t = None
        if t is None: break
        try:
            v = basis(a, b, h)
        except Exception:
            break
        if v is None or not np.all(np.isfinite(v)): break
        rows.append((f, t[0], t[1], v))
    return hc, rows

for (a, b) in [(-0.5, 1.0), (-0.5, 0.5), (-0.5, 1.5), (-3.0, 1.0), (-0.1, 1.0), (-5.0, 0.3)]:
    hc, rows = scan(a, b)
    print("a=%-5g b=%-4g hc=%.8f  annulus reached f=h/hc up to %.4g in %d samples"
          % (a, b, hc, rows[-1][0] if rows else float('nan'), len(rows)))
    if rows:
        f, y1, y2, v = rows[-1]
        print("     outermost sampled oval y in [%.6g, %.6g]" % (y1, y2))
