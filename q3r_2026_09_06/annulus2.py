#!/usr/bin/env python3
"""Period annulus of the upper centre, found by additive continuation in h,
plus the four generating functions sampled across it."""
import numpy as np, mpmath as mp
from melnikov import annulus_h, turning, basis
mp.mp.dps = 30

def direction(a, b, hc):
    for d in (1.0, -1.0):
        for e in (1e-6, 1e-4, 1e-2):
            try:
                if turning(a, b, hc + d*e) is not None: return d
            except Exception:
                pass
    return None

def annulus(a, b, n=48):
    hc = annulus_h(a, b)
    d = direction(a, b, hc)
    if d is None: return hc, None, []
    # grow until the oval stops existing
    lo, hi = 1e-8, 1e-8
    while hi < 1e8:
        try:
            if turning(a, b, hc + d*hi*2) is None: break
        except Exception:
            break
        hi *= 2
    es = np.geomspace(max(lo, hi*1e-6), hi, n)
    rows = []
    for e in es:
        h = hc + d*e
        try:
            v = basis(a, b, h)
        except Exception:
            continue
        if v is None or not np.all(np.isfinite(v)): continue
        rows.append((h, v))
    return hc, d, rows

def max_sign_changes(F, ntry=200000, seed=0):
    """max number of sign changes of c.F over random directions c, refined."""
    rng = np.random.default_rng(seed)
    best, bestc = 0, None
    C = rng.normal(size=(ntry, F.shape[0]))
    V = C @ F
    sc = (np.diff(np.sign(V), axis=1) != 0).sum(axis=1)
    k = int(sc.argmax()); best = int(sc[k]); bestc = C[k]
    return best, bestc

if __name__ == "__main__":
    print("family: reversible two centres, 0<b<2, a<0, a not in {-1,-2}")
    print("generating space is 4-dimensional -> a Chebyshev system would cap zeros at 3\n")
    print(" a      b     hc          dir  samples  max sign changes of M1")
    for a in (-0.1, -0.3, -0.5, -0.8, -1.5, -2.5, -3.0, -5.0):
        for b in (0.3, 0.7, 1.0, 1.4, 1.8):
            hc, d, rows = annulus(a, b)
            if not rows:
                print(" %-6g %-5g %-11.6g  --   0        n/a" % (a, b, hc)); continue
            F = np.array([r[1] for r in rows]).T
            # normalise each generator so the LP/random search is well scaled
            F = F/np.abs(F).max(axis=1, keepdims=True)
            m, c = max_sign_changes(F)
            print(" %-6g %-5g %-11.6g  %+d   %-8d %d" % (a, b, hc, d, F.shape[1], m))
