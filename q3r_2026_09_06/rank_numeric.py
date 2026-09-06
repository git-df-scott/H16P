#!/usr/bin/env python3
"""Decisive test of Astra's relation against the generators THIS repository
actually computed (melnikov.basis / fast.gens):

    (a+2) U + 3[ (b-2)/4 T_{a-2} + (1-b) T_{a-1} + b T_a ]  ==  0 ?

If it holds, the four "generators" are dependent, the 4x4 determinants are
identically zero, and the 70-digit Chebyshev test in Q3R_FIRST_ORDER.md was
measuring nothing."""
import numpy as np, mpmath as mp
import melnikov as M
from fast import gens, hcentre, turning
mp.mp.dps = 40

print("  a      b      h-hc      T_{a-2}        T_{a-1}        T_a            U              relation      rel/|terms|")
for (a, b) in [(-0.5, 1.0), (-0.25, 0.9), (-0.3, 0.7), (-2.5, 1.2), (-0.8, 1.6)]:
    hc = hcentre(a, b)
    for e in (0.01, 0.2):
        for d in (+1, -1):
            h = hc + d*e
            try:
                if turning(a, b, h) is None: continue
                g = M.basis(a, b, h)
            except Exception:
                continue
            if g is None or not np.all(np.isfinite(g)): continue
            T2, T1, T0, U = g
            rel = (a+2)*U + 3*((b-2)/4*T2 + (1-b)*T1 + b*T0)
            scale = abs((a+2)*U) + 3*(abs((b-2)/4*T2) + abs((1-b)*T1) + abs(b*T0))
            print("  %-6g %-6g %+-9.3g %-14.8g %-14.8g %-14.8g %-14.8g %-13.3e %.2e"
                  % (a, b, d*e, T2, T1, T0, U, rel, abs(rel)/max(scale, 1e-300)))
            break
