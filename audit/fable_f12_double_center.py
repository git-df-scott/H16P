#!/usr/bin/env python3
"""F12: the KKL double center (beta=0, K=0, J(c)=0) and its first-order Melnikov functions.
x' = y + x^2 + xy,  y' = -10x^2 + (11/5)xy + c y^2 + alpha x,  K = -alpha(11c/5-1) - 42 = 0,
J(c) = 305 + 634c - 11c^2 - 1000c^3 = 0.
Step 1: verify numerically that the origin is a center (return displacement ~ 0) and locate the annulus.
Step 2: find an integrating factor numerically? Instead compute M(h) directly: for a perturbation (f,g),
the first-order displacement along the closed orbit gamma_h is proportional to
  int_0^T exp(-int_0^t div F ds) * (f, g) x F / |F|^2 ... ; we use the standard formula
  M(h) = int_0^T (f Q - g P)(gamma(t)) * exp(-int_0^t div F) dt   (P,Q = unperturbed field),
which is the first variation of the return along the normal direction (valid for any center, integrating
factor mu = exp(-int div) along orbits). Compute it for the 12 monomial directions and find the span.
"""
import numpy as np, sys, os
from scipy.integrate import solve_ivp
from scipy.optimize import brentq
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), 'fable_engine'))
import retmap as rm
cs = brentq(lambda c: 305+634*c-11*c*c-1000*c**3, 0.9, 1.0, xtol=1e-15)
al = -42/(11*cs/5-1)
coef = np.array([0, 0, 1, 1, 1, 0, 0, al, 0, -10, 11/5, cs], float)
print(f"c* = {cs:.15f}  alpha* = {al:.15f}")
print("equilibria:", rm.equilibria(coef))
for a in rm.antisaddles(coef): print("  antisaddle", a)
# step 1: displacement on the origin nest
rad = np.geomspace(1e-2, 200, 25)
R, T, st = rm.returns(coef[None], np.array([[0., 0.]]), np.array([[1., 0.]]), rad[None], 1e-12, 1e8, 1e4)
for r, RR, TT, s in zip(rad, R[0], T[0], st[0]):
    print(f"  r={r:9.4f}  D={RR-r:+.3e}  T={TT:.4f}  st={s}")
def F(u):
    x, y = u
    return np.array([coef[0]+coef[1]*x+coef[2]*y+coef[3]*x*x+coef[4]*x*y+coef[5]*y*y,
                     coef[6]+coef[7]*x+coef[8]*y+coef[9]*x*x+coef[10]*x*y+coef[11]*y*y])
def div(u):
    x, y = u
    return coef[1]+2*coef[3]*x+coef[4]*y + coef[8]+coef[10]*x+2*coef[11]*y
MON = [(0,0),(1,0),(0,1),(2,0),(1,1),(0,2)]
def rhs(t, u):
    x, y = u[0], u[1]; P, Q = F(u[:2]); w = np.exp(-u[2])
    out = [P, Q, div(u[:2])]
    for (i, j) in MON:   # f = x^i y^j (perturbing x'), contributes f*Q ; g contributes -g*P
        out.append(w*(x**i*y**j)*Q)
    for (i, j) in MON:
        out.append(-w*(x**i*y**j)*P)
    return out
def melnikov_row(x0):
    ev = lambda t, u: u[1]
    ev.direction = -1 if F([x0, 0.0])[1] < 0 else 1   # detect return to y=0 with same crossing sense
    sol = solve_ivp(rhs, [0, 200], [x0, 0.0, 0.0]+[0.0]*12, rtol=1e-12, atol=1e-14, events=ev, max_step=0.05)
    te = sol.t_events[0]; te = te[te > 1e-3]
    if len(te) == 0: return None
    Tp = te[0]
    sol2 = solve_ivp(rhs, [0, Tp], [x0, 0.0, 0.0]+[0.0]*12, rtol=1e-12, atol=1e-14, max_step=0.05)
    u = sol2.y[:, -1]
    return Tp, u[0], u[2], u[3:]
rows = []; xs = []
for x0 in np.geomspace(0.05, 60, 45):
    r = melnikov_row(x0)
    if r is None or abs(r[1]-x0) > 1e-5*(1+x0): print(f"  x0={x0:.4f}: not closed (return {r[1] if r else None}) -> annulus edge"); break
    rows.append(r[3]); xs.append(x0)
B = np.array(rows); xs = np.array(xs)
print("ovals:", len(xs), "x0 range", xs[0], xs[-1])
Bn = B/np.max(np.abs(B), axis=0)
print("singular values:", np.linalg.svd(Bn, compute_uv=False))
rng = np.random.default_rng(0); hist = {}; best = (0, None)
for _ in range(300000):
    cvec = rng.standard_normal(12); M = Bn @ cvec
    z = int(np.sum(np.sign(M[:-1])*np.sign(M[1:]) < 0)); hist[z] = hist.get(z, 0)+1
    if z > best[0]: best = (z, cvec)
print("zero histogram over random directions:", dict(sorted(hist.items())))
print("max zeros", best[0], "direction", None if best[1] is None else np.round(best[1], 4).tolist())
np.save(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'fable_f12_basis.npy'), np.column_stack([xs, B]))
