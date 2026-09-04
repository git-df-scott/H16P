#!/usr/bin/env python3
"""Fable lane, step 2 (first-order compatibility gate). For lobe-region points
(three primitive roots, first root > 5/11) and lift parameters a, compute the
boundary value functional X(1) = int_0^1 Y/(1-au)^{3/2} du via the repository
ODE (t_end = 1-1e-6). I(loop) = -(aC/2) sqrt(1-a) X(1), so c0 = 0 <=> X(1) = 0.
With Y0<0, I>0 near the center; X(1)>0 <=> I(loop)<0 <=> odd interior count.
A sign change of X(1) inside the lobe region would be the only place a
first-order zero can reach the graphic."""
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "q4")); sys.path.insert(0, os.path.dirname(__file__))
import numpy as np, mpmath as mp
from q4_threshold_path import from_primitive_anchors_closed
from q4_reconstruction import reconstruct
mp.mp.dps = 30
triples = [(0.46,0.6,0.8),(0.46,0.9,0.99),(0.5,0.75,0.9),(0.6,0.8,0.95),(0.7,0.75,0.8),(0.75,0.875,0.9375),
           (0.9,0.95,0.99),(0.99,0.999,0.9999),(0.5,0.999,0.9999),(0.8,0.81,0.82)]
As = (0.15, 0.3, 0.5, 0.7, 0.85, 0.95, 0.99)
print("X(1-1e-6) for lobe points; rows anchors, cols a =", As)
neg = 0
for y in triples:
    co = tuple(map(float, from_primitive_anchors_closed(tuple(map(mp.mpf, map(str, y))))))
    row = []
    for a in As:
        sol = reconstruct(a, *co, t_end=1-1e-6)
        X1 = sol.sol(1-1e-6)[3]; row.append(X1)
        if X1 < 0: neg += 1
    print(f"  {y}: " + " ".join(f"{v:+.3e}" for v in row))
print("corner point (94/77,-17/77,1):")
co = (94/77, -17/77, 1.0)
print("  " + " ".join(f"{reconstruct(a,*co,t_end=1-1e-6).sol(1-1e-6)[3]:+.3e}" for a in As))
print("number of negative X(1) values found:", neg)
