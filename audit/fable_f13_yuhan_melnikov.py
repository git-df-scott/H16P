#!/usr/bin/env python3
"""F13: first-order Melnikov functions of the Yu-Han reversible two-center family
   x' = y(1 + a1 x),  y' = -x + x^2 + a4 y^2      (centers at (0,0) and (1,0); invariant line x = -1/a1)
on the origin annulus (and optionally the (1,0) annulus), for a set of (a1,a4).
Reports: span dimension, maximal zero count over random directions, and the zero counts of the
constrained elements with v0=v1=0 and with v0=v1=v2=0 (Taylor coefficients at the center).
The v0=v1=v2=0 element exists only where the 3x3 Taylor matrix is singular (Yu-Han's curve).
Melnikov: M = int_0^T (f Q - g P) exp(-int div) dt along the periodic orbit for 12 monomial directions.
Usage: python3 fable_f13_yuhan_melnikov.py OUT.jsonl [second]
"""
import sys, os, json, numpy as np
from scipy.integrate import solve_ivp
from multiprocessing import Pool
MON = [(0,0),(1,0),(0,1),(2,0),(1,1),(0,2)]
SIDE = 'second' if len(sys.argv) > 2 and sys.argv[2] == 'second' else 'origin'

def make(a1, a4):
    def F(x, y): return y*(1+a1*x), -x+x*x+a4*y*y
    def div(x, y): return a1*y + 2*a4*y
    def rhs(t, u):
        x, y = u[0], u[1]; P, Q = F(x, y); w = np.exp(-u[2])
        out = [P, Q, div(x, y)]
        out += [w*(x**i*y**j)*Q for (i, j) in MON]
        out += [-w*(x**i*y**j)*P for (i, j) in MON]
        return out
    return F, rhs

def row(rhs, F, x0):
    ev = lambda t, u: u[1]; ev.direction = -1 if F(x0, 0.0)[1] < 0 else 1
    sol = solve_ivp(rhs, [0, 200], [x0, 0, 0]+[0]*12, rtol=1e-11, atol=1e-13, events=ev, max_step=0.05)
    te = sol.t_events[0]; te = te[te > 1e-3]
    if len(te) == 0: return None
    s2 = solve_ivp(rhs, [0, te[0]], [x0, 0, 0]+[0]*12, rtol=1e-11, atol=1e-13, max_step=0.05)
    u = s2.y[:, -1]
    if abs(u[0]-x0) > 1e-6*(1+abs(x0)): return None
    return u[3:]

def annulus_edge(rhs, F, c, lim):
    """largest |x0 - c| (toward lim) with a closed oval, by bisection."""
    lo, hi = 0.0, abs(lim-c)
    sgn = 1.0 if lim > c else -1.0
    for _ in range(30):
        mid = 0.5*(lo+hi)
        if row(rhs, F, c+sgn*mid) is not None: lo = mid
        else: hi = mid
        if hi-lo < 1e-6: break
    return lo, sgn

def analyse(args):
    a1, a4 = args
    F, rhs = make(a1, a4)
    line = -1.0/a1
    c = 0.0 if SIDE == 'origin' else 1.0
    if SIDE == 'origin':
        lim = line if -1 < line < 1e9 and line > 0 else 1.0
        lim = min(lim, 1.0)
    else:
        lim = line if 0 < line < 1 else 0.0
    edge, sgn = annulus_edge(rhs, F, c, lim)
    if edge < 1e-3: return dict(a1=a1, a4=a4, side=SIDE, error="no annulus")
    rs = np.concatenate([edge*np.geomspace(0.01, 0.3, 20), edge*np.linspace(0.32, 0.995, 45)])
    rows = []; rv = []
    for r in rs:
        m = row(rhs, F, c+sgn*r)
        if m is None: break
        rows.append(m); rv.append(r)
    if len(rows) < 30: return dict(a1=a1, a4=a4, side=SIDE, error="few ovals", n=len(rows), edge=edge)
    B = np.array(rows); rv = np.array(rv)
    Bs = B/(rv**2)[:, None]                       # M ~ h ~ r^2 near the center: equalise rows
    colmax = np.max(np.abs(Bs), axis=0); keep = colmax > 1e-7*colmax.max()
    Bk = Bs[:, keep]/colmax[keep]
    U, s, Vt = np.linalg.svd(Bk, full_matrices=False)
    dim = int(np.sum(s > 1e-6*s[0]))
    basis = U[:, :dim]*s[:dim]
    rng = np.random.default_rng(0); best = 0; hist = {}
    for _ in range(30000):
        M = basis @ rng.standard_normal(dim)
        z = int(np.sum(np.sign(M[:-1])*np.sign(M[1:]) < 0)); hist[z] = hist.get(z, 0)+1
        best = max(best, z)
    # Taylor structure at the center (basis/r^2 ~ c0 + c1 r^2 + c2 r^4), fitted on the 14 smallest ovals
    k = 14; A = np.column_stack([np.ones(k), rv[:k]**2, rv[:k]**4])
    C = np.linalg.lstsq(A, basis[:k], rcond=None)[0]       # 3 x dim
    _, _, vt2 = np.linalg.svd(C[:2]); ker2 = vt2[2:]
    z2 = [int(np.sum(np.sign((basis @ v)[k:-1])*np.sign((basis @ v)[k+1:]) < 0)) for v in ker2]
    s3 = np.linalg.svd(C, compute_uv=False)
    _, _, vt3 = np.linalg.svd(C); Phi3 = basis @ vt3[-1]
    Phi3n = Phi3/np.max(np.abs(Phi3))
    sig = np.abs(Phi3n) > 1e-5                      # ignore the fit-residual floor near the center
    pr = Phi3n[sig]; rr = rv[sig]
    zi = np.nonzero(np.sign(pr[:-1])*np.sign(pr[1:]) < 0)[0]
    z3 = int(len(zi)); z3_at = [float(rr[j]) for j in zi]
    return dict(phi3_profile=[[float(a), float(b)] for a, b in zip(rv, Phi3n)], phi3_zero_r=z3_at, phi3_dir=[float(v) for v in vt3[-1]],a1=a1, a4=a4, side=SIDE, dim=dim, sv=[float(v) for v in s[:4]], max_zeros=best,
                zero_hist={int(a): int(b) for a, b in hist.items()}, phi2_zeros=z2,
                taylor_sv=[float(v) for v in s3], taylor_rank_ratio=float(s3[-1]/s3[0]), phi3_zeros=z3,
                n=len(rv), edge=float(edge), kept=int(keep.sum()))

if __name__ == "__main__":
    out = sys.argv[1]
    if len(sys.argv) > 3 and sys.argv[3] == 'test':
        print(analyse((-30/7, -65/21))); sys.exit()
    a1s = [-1.1, -1.25, -1.5, -2, -2.5, -3, -30/7, -5, -6, -8, -10, -15, -25]
    grid = [(a1, (a1-5)/3) for a1 in a1s] + [(a1, (a1-5)/3 + d) for a1 in (-30/7, -2, -8) for d in (-0.5, 0.5)]
    with Pool(3) as p, open(out, 'a') as f:
        for r in p.imap_unordered(analyse, grid):
            print(r, flush=True); f.write(json.dumps(r)+"\n"); f.flush()
