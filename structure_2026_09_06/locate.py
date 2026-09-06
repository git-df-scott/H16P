#!/usr/bin/env python3
"""Locate the graphic connection inside the order-2 + neutral-graphic family:
bisect the splitting in a, at fixed l."""
import sys, numpy as np
from scipy.optimize import brentq
from splitting_gen import params, split

def splitting(lv, av):
    r = split(lv, float(av))
    if r is None: return None
    mv, bv, res = r
    us = sorted(set(k[0] for k in res))
    if len(us) != 2: return None
    u1, u2 = us[0], us[1]           # u1 negative branch, u2 positive branch
    A = res.get((u1, True));  B = res.get((u2, False))
    if A is None or B is None: return None
    ya = A[1]/A[2]; yb = B[1]/B[2]
    return ya - yb

for lv in (-3.0, -6.0, -12.0, -1.5):
    lo, hi = 0.35, 1.2
    flo, fhi = splitting(lv, lo), splitting(lv, hi)
    if flo is None or fhi is None or flo*fhi > 0:
        print("l=%-6g no bracket (%s, %s)" % (lv, flo, fhi)); continue
    astar = brentq(lambda t: splitting(lv, t), lo, hi, xtol=1e-13, rtol=1e-14)
    mv, bv = params(lv, astar)
    print("l=%-6g  a* = %.12f   m = %.12f   b = %.12f   splitting = %+.2e"
          % (lv, astar, mv, bv, splitting(lv, astar)), flush=True)
