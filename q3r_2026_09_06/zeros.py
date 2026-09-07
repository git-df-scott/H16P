#!/usr/bin/env python3
"""How many zeros can M1 have on the period annulus of the upper centre?

M1 lies in the 4-dimensional space spanned by T_{a-2},T_{a-1},T_a,U.
Two independent tests on a dense h-sample:
 (1) Chebyshev test: sign of det[f_i(h_{j_k})] over 4-point subsets.  Constant
     sign on all subsets => extended Chebyshev on the sample => at most 3 zeros.
 (2) direct search: max sign changes of c.F over many directions c, refined by
     a linear program looking for 4 alternations (which would mean >=4 zeros,
     i.e. a 4+1 five-cycle lead).
"""
import numpy as np
from itertools import combinations
from scipy.optimize import linprog
from fast import gens, hcentre, turning

def sample(a, b, n=140):
    hc = hcentre(a, b)
    d = None
    for dd in (1.0, -1.0):
        if turning(a, b, hc + dd*1e-6) is not None: d = dd; break
    if d is None: return None, None
    hi = 1e-6
    while hi < 1e10 and turning(a, b, hc + d*hi*1.6) is not None: hi *= 1.6
    es = np.geomspace(1e-7*max(hi, 1e-6), hi, n)
    H, V = [], []
    for e in es:
        g = gens(a, b, hc + d*e)
        if g is None or not np.all(np.isfinite(g)): continue
        H.append(hc + d*e); V.append(g)
    if len(V) < 20: return None, None
    F = np.array(V).T
    F = F/np.abs(F).max(axis=1, keepdims=True)
    return np.array(H), F

def cheb_test(F, ntry=120000, seed=1):
    n = F.shape[1]; rng = np.random.default_rng(seed)
    idx = [np.arange(i, i+4) for i in range(n-3)]
    idx += [np.sort(rng.choice(n, 4, replace=False)) for _ in range(ntry)]
    dets = np.array([np.linalg.det(F[:, j]) for j in idx])
    dets = dets[np.abs(dets) > 1e-14*np.abs(dets).max()]
    if dets.size == 0: return None
    return float(np.sign(dets).min()), float(np.sign(dets).max()), float(np.abs(dets).min()/np.abs(dets).max())

def max_alt(F, ntry=300000, seed=2):
    rng = np.random.default_rng(seed)
    C = rng.normal(size=(ntry, 4))
    V = C @ F
    sc = (np.diff(np.sign(V), axis=1) != 0).sum(axis=1)
    return int(sc.max())

def lp_four(F, H, ntry=4000, seed=3):
    """Is there c with 4 sign alternations at some 5 sample points?"""
    n = F.shape[1]; rng = np.random.default_rng(seed)
    sets = [np.sort(rng.choice(n, 5, replace=False)) for _ in range(ntry)]
    sets += [np.linspace(0, n-1, 5).astype(int)]
    for js in sets:
        M = F[:, js].T
        for s0 in (1, -1):
            sg = np.array([s0*(-1)**k for k in range(5)])
            A_ub = -(M*sg[:, None]); b_ub = -np.ones(5)
            r = linprog(np.zeros(4), A_ub=A_ub, b_ub=b_ub,
                        bounds=[(-1e4, 1e4)]*4, method='highs')
            if r.status == 0:
                return True, js, r.x
    return False, None, None

print(" a      b     samples  det-sign(min,max)  cond   max sign changes  4-alternation?")
lead = []
for a in (-0.1, -0.25, -0.5, -0.75, -1.5, -2.5, -3.5, -6.0):
    for b in (0.2, 0.5, 0.9, 1.2, 1.6, 1.9):
        H, F = sample(a, b)
        if F is None:
            print(" %-6g %-5g  --" % (a, b)); continue
        ct = cheb_test(F); ma = max_alt(F)
        ok, js, c = lp_four(F, H)
        print(" %-6g %-5g %-8d (%+.0f,%+.0f)  %.2e   %d                %s"
              % (a, b, F.shape[1], ct[0], ct[1], ct[2], ma, "YES" if ok else "no"))
        if ok: lead.append((a, b, js, c))
print("\nparameter points admitting four sign alternations:", len(lead))
for L in lead[:5]: print("   a=%g b=%g  c=%s" % (L[0], L[1], np.array2string(L[3], precision=6)))
