"""F6: first-order Melnikov functions on the period annulus of the reversible center
X' = -Y(1+kX), Y' = X + pX^2 + qY^2 (Shi loop point, a=1; CLAUDE_ROUTES_4AB.md).
Integrating factor mu = 2(1+kX)^(2q/k-1). For a quadratic perturbation (f,g),
M(h) = closed integral of mu (g dX - f dY) = int_0^T mu (g X' - f Y') dt.
Reversible directions (f odd in Y, g even in Y) give M = 0; the six even ones remain.
"""
import numpy as np
from scipy.integrate import solve_ivp
k, p, q = 5.54048179, -1.24519487, 0.22849752
Xs = -1.0/p                      # saddle on the axis
def rhs(t, u):
    X, Y = u[0], u[1]
    Xd = -Y*(1+k*X); Yd = X + p*X*X + q*Y*Y
    mu = 2*(1+k*X)**(2*q/k-1)
    # basis: f in {1, x, x^2, y^2} -> -mu f Yd ; g in {y, xy} -> mu g Xd
    return [Xd, Yd, -mu*Yd, -mu*X*Yd, -mu*X*X*Yd, -mu*Y*Y*Yd, mu*Y*Xd, mu*X*Y*Xd]
def period_integrals(X0):
    # integrate from (X0,0) until return to Y=0 with X>0 after a full turn
    ev = lambda t, u: u[1]
    ev.direction = 1   # Y crossing upward: full turn (start goes downward first since Xd<0 for Y>0... choose by event count)
    sol = solve_ivp(rhs, [0, 500], [X0, 0]+[0]*6, rtol=1e-12, atol=1e-14, events=ev, dense_output=False, max_step=0.05)
    te = sol.t_events[0]
    te = te[te > 1e-6]
    if len(te) == 0:
        return None
    T = te[0]
    sol2 = solve_ivp(rhs, [0, T], [X0, 0]+[0]*6, rtol=1e-12, atol=1e-14, max_step=0.05)
    u = sol2.y[:, -1]
    return T, u[2:], u[0]
X0s = np.linspace(0.01, Xs-1e-3, 160)
rows = []; Ts = []; xs = []
for X0 in X0s:
    r = period_integrals(X0)
    if r is None or abs(r[2]-X0) > 1e-6:
        print("orbit not closed at X0=%.4f (return %.6f)" % (X0, r[2] if r else np.nan)); break
    rows.append(r[1]); Ts.append(r[0]); xs.append(X0)
B = np.array(rows); xs = np.array(xs); Ts = np.array(Ts)
print("ovals computed:", len(xs), "X0 range", xs[0], xs[-1], "period range", Ts[0], Ts[-1])
# normalise columns and check span dimension
Bn = B/np.max(np.abs(B), axis=0)
s = np.linalg.svd(Bn, compute_uv=False)
print("singular values:", s)
# maximum zero count over random combinations (zeros in the open annulus, ignore the center itself)
rng = np.random.default_rng(1)
best = 0; bestc = None; hist = {}
for _ in range(200000):
    c = rng.standard_normal(6)
    M = Bn @ c
    z = int(np.sum(np.sign(M[:-1])*np.sign(M[1:]) < 0))
    hist[z] = hist.get(z, 0)+1
    if z > best:
        best, bestc = z, c
print("zero-count histogram over random directions:", dict(sorted(hist.items())))
print("max zeros:", best, "direction:", bestc)
# targeted: solve for combos with prescribed zeros at 4 or 5 chosen points, then count
for trial in range(2000):
    pts = np.sort(rng.choice(len(xs)-1, size=5, replace=False))
    A = Bn[pts]              # 5 x 6 -> one-dimensional null space
    _, _, vt = np.linalg.svd(A)
    c = vt[-1]; M = Bn @ c
    z = int(np.sum(np.sign(M[:-1])*np.sign(M[1:]) < 0))
    if z > best:
        best, bestc = z, c; print("targeted found", z, "zeros near indices", pts)
print("max zeros after targeting:", best)
np.save('/home/user/H16P/audit/fable_q3r_basis.npy', np.column_stack([xs, Ts, B]))
