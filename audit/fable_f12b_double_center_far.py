#!/usr/bin/env python3
"""F12b: KKL double center (beta=0,K=0,J=0): first-order Melnikov basis on the origin annulus out to large
amplitude (toward the neutral two-saddle infinity graphic), 3-zero directions and their zero locations."""
import numpy as np, sys, os, json
from scipy.integrate import solve_ivp
from scipy.optimize import brentq
cs = brentq(lambda c: 305+634*c-11*c*c-1000*c**3, 0.9, 1.0, xtol=1e-15); al = -42/(11*cs/5-1)
coef = np.array([0, 0, 1, 1, 1, 0, 0, al, 0, -10, 11/5, cs], float)
MON = [(0,0),(1,0),(0,1),(2,0),(1,1),(0,2)]
def F(x, y): return coef[2]*y+coef[3]*x*x+coef[4]*x*y, coef[7]*x+coef[9]*x*x+coef[10]*x*y+coef[11]*y*y
def div(x, y): return 2*coef[3]*x+coef[4]*y + coef[10]*x+2*coef[11]*y
def rhs(t, u):
    x, y = u[0], u[1]; P, Q = F(x, y); w = np.exp(-u[2])
    return [P, Q, div(x, y)] + [w*(x**i*y**j)*Q for (i, j) in MON] + [-w*(x**i*y**j)*P for (i, j) in MON]
def row(x0):
    ev = lambda t, u: u[1]; ev.direction = -1 if F(x0, 0.0)[1] < 0 else 1
    sol = solve_ivp(rhs, [0, 50], [x0, 0, 0]+[0]*12, rtol=1e-11, atol=1e-13, events=ev, max_step=0.02)
    te = sol.t_events[0]; te = te[te > 1e-3]
    if len(te) == 0: return None
    s2 = solve_ivp(rhs, [0, te[0]], [x0, 0, 0]+[0]*12, rtol=1e-11, atol=1e-13, max_step=0.02)
    u = s2.y[:, -1]
    if abs(u[0]-x0) > 1e-6*(1+abs(x0)): return None
    return u[3:]
xs = np.geomspace(0.05, float(sys.argv[1]) if len(sys.argv) > 1 else 3000.0, 90)
rows = []; xv = []
for x0 in xs:
    r = row(x0)
    if r is None: print("edge before x0 =", x0, flush=True); break
    rows.append(r); xv.append(x0)
B = np.array(rows); xv = np.array(xv)
np.save(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'fable_f12b_basis.npy'), np.column_stack([xv, B]))
Bs = B/(xv**2)[:, None]; colmax = np.max(np.abs(Bs), axis=0); keep = colmax > 1e-7*colmax.max()
Bk = Bs[:, keep]/colmax[keep]; U, s, Vt = np.linalg.svd(Bk, full_matrices=False); dim = int(np.sum(s > 1e-6*s[0]))
print("ovals", len(xv), "x0 max", xv[-1], "sv", np.round(s[:5], 7), "dim", dim)
basis = U[:, :dim]*s[:dim]; rng = np.random.default_rng(2); hist = {}; three = []
for _ in range(400000):
    v = rng.standard_normal(dim); M = basis @ v
    zi = np.nonzero(np.sign(M[:-1])*np.sign(M[1:]) < 0)[0]; z = len(zi); hist[z] = hist.get(z, 0)+1
    if z >= 3 and len(three) < 10: three.append((v.tolist(), [float(xv[j]) for j in zi], float(np.min(np.abs(M))/np.max(np.abs(M)))))
print("zero histogram:", dict(sorted(hist.items())))
for v, zs, ratio in three: print("3-zero dir", np.round(v, 4).tolist(), "zeros at x0 =", np.round(zs, 4).tolist(), "min|M|/max|M|", f"{ratio:.1e}")
