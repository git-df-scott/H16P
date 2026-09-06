#!/usr/bin/env python3
"""(D): map the loop locus.  For each (a,b) find EVERY l giving a homoclinic
loop around the origin, for every non-origin saddle, and record sign(sigma*eta_2).
By the symmetry (x,t)->(-x,-t), a -> -a, so a<0 loses nothing."""
import sys, numpy as np
from scipy.optimize import brentq
from d2core import saddles, splitting, eta2, C3, mval

def loops_at(a, b, lgrid, which_saddle=None):
    """All l in lgrid brackets with a loop; returns list of dicts."""
    prev = {}
    found = []
    for l in lgrid:
        if abs(l + 1) < 1e-6: continue
        try:
            sad = saddles(a, b, l)
        except Exception:
            continue
        cur = {}
        for k, (x, y, sg) in enumerate(sad):
            key = ("origin1" if (abs(x) < 1e-12 and abs(y-1) < 1e-12) else "line%d" % k)
            try:
                sp = splitting(a, b, l, (x, y))
            except Exception:
                sp = None
            cur[key] = (l, sp, sg)
            if key in prev:
                l0, sp0, sg0 = prev[key]
                if sp0 is not None and sp is not None and sp0*sp < 0:
                    def f(t):
                        s2 = saddles(a, b, t)
                        cand = [(xx, yy, ss) for (xx, yy, ss) in s2]
                        if not cand: return np.nan
                        # match by continuity to x
                        xx, yy, ss = min(cand, key=lambda c: (c[0]-x)**2 + (c[1]-y)**2)
                        v = splitting(a, b, t, (xx, yy))
                        return np.nan if v is None else v
                    try:
                        ls = brentq(f, l0, l, xtol=1e-13, rtol=8.9e-16, maxiter=200)
                    except Exception:
                        continue
                    s2 = saddles(a, b, ls)
                    if not s2: continue
                    xx, yy, ss = min(s2, key=lambda c: (c[0]-x)**2 + (c[1]-y)**2)
                    e2 = eta2(a, b, ls)
                    found.append(dict(a=a, b=b, l=ls, saddle=(xx, yy), sigma=ss,
                                      eta2=e2, prod=ss*e2, C3=C3(a, b, ls), key=key))
        prev = cur
    return found

if __name__ == "__main__":
    a = float(sys.argv[1]); bs = [float(v) for v in sys.argv[2].split(",")]
    lg = np.concatenate([np.arange(-6.0, -1.02, 0.04), np.arange(-0.98, 4.0, 0.04)])
    print("a=%g   scanning l in [-6,4] step 0.04" % a)
    print(" b        l            saddle                  sigma        eta_2        sigma*eta_2   C3")
    for b in bs:
        for r in loops_at(a, b, lg):
            print(" %-8.3f %-12.8f (%+.5f,%+.5f)  %+.4e  %+.4e  %+.4e  %+.2e  %s"
                  % (r['b'], r['l'], r['saddle'][0], r['saddle'][1], r['sigma'],
                     r['eta2'], r['prod'], r['C3'],
                     "NEGATIVE" if r['prod'] < 0 else "*** POSITIVE ***"), flush=True)
