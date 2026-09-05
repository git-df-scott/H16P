#!/usr/bin/env python3
"""Fable lane, step 3: interior zero count of the first-order Q4 integral on
both sides of the X(1)=0 face, for lobe points, with log-dense sampling near
the loop. Determines whether the zero that crosses the graphic enters
(2 -> 3) or leaves (2 -> 1)."""
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "q4"))
import numpy as np, mpmath as mp
from q4_threshold_path import from_primitive_anchors_closed
from q4_reconstruction import reconstruct, original_values
from scipy.optimize import brentq
mp.mp.dps = 30
ts = np.concatenate([np.linspace(1e-3, 0.9, 400), 1-np.logspace(-1, -7, 400)])
def count(a, co):
    sol = reconstruct(a, *co, t_end=1-1e-7)
    I = original_values(a, sol, ts); s = np.sign(I)
    ch = np.nonzero(s[:-1]*s[1:] < 0)[0]
    return len(ch), [float(ts[i]) for i in ch], float(sol.sol(1-1e-7)[3])
for y in ((0.46,0.6,0.8), (0.5,0.75,0.9), (0.9,0.95,0.99), (0.75,0.875,0.9375)):
    co = tuple(map(float, from_primitive_anchors_closed(tuple(map(mp.mpf, map(str, y))))))
    X1 = lambda a: reconstruct(a, *co, t_end=1-1e-7).sol(1-1e-7)[3]
    try: a_star = brentq(X1, 0.6, 0.995, xtol=1e-8)
    except ValueError: a_star = None
    print(f"anchors {y}: X(1)=0 at a*={a_star}")
    for a in ([0.85, 0.99] if a_star is None else [a_star-0.05, a_star-0.005, a_star+0.005, a_star+0.05]):
        n, locs, X1v = count(a, co)
        print(f"   a={a:.5f}: X(1)={X1v:+.3e}  interior sign changes={n} at t={np.round(locs,5)}")
