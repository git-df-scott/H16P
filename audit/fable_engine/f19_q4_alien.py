"""F19: alien cycles at the Q4 neutral resonant infinity graphic.
Q4 center (rho=1): X' = Y + 6X^2 + 4XY - 2Y^2, Y' = -X + 2X^2 + 8XY - 2Y^2; annulus boundary x* = 0.2272111321 on the +x axis.
1. First-order Melnikov basis on ovals (12 monomial directions), span, random directions -> zero counts.
2. Directions with 3 zeros and with 2 zeros (one near the boundary): perturb the field by eps*dir and count
   actual cycles around the origin with the compactified counter (radius to e^40). Excess over the first-order
   count = alien candidate."""
import sys, json, time, numpy as np
from scipy.integrate import solve_ivp
sys.argv = ['x', 'mvneutral', '/dev/null', '0', '0']
import importlib.util; spec = importlib.util.spec_from_file_location('sl', 'sweep_log.py'); sl = importlib.util.module_from_spec(spec); spec.loader.exec_module(sl)
import retmap as rm
Q4 = np.array([0, 0, 1, 6, 4, -2, 0, -1, 0, 2, 8, -2], float)
MON = [(0,0),(1,0),(0,1),(2,0),(1,1),(0,2)]
def F(x, y): return Q4[2]*y+Q4[3]*x*x+Q4[4]*x*y+Q4[5]*y*y, Q4[7]*x+Q4[9]*x*x+Q4[10]*x*y+Q4[11]*y*y
def div(x, y): return 2*Q4[3]*x+Q4[4]*y + Q4[10]*x+2*Q4[11]*y
def rhs(t, u):
    x, y = u[0], u[1]; P, Q = F(x, y); w = np.exp(-u[2])
    return [P, Q, div(x, y)] + [w*(x**i*y**j)*Q for (i, j) in MON] + [-w*(x**i*y**j)*P for (i, j) in MON]
def row(x0):
    ev = lambda t, u: u[1]; ev.direction = -1 if F(x0, 0.0)[1] < 0 else 1
    sol = solve_ivp(rhs, [0, 400], [x0, 0, 0]+[0]*12, rtol=1e-11, atol=1e-13, events=ev, max_step=0.02)
    te = sol.t_events[0]; te = te[te > 1e-3]
    if len(te) == 0: return None
    s2 = solve_ivp(rhs, [0, te[0]], [x0, 0, 0]+[0]*12, rtol=1e-11, atol=1e-13, max_step=0.02)
    u = s2.y[:, -1]
    if abs(u[0]-x0) > 1e-6*(1+abs(x0)): return None
    return u[3:]
xstar = 0.2272111321
xs = np.concatenate([np.geomspace(0.003, 0.15, 30), xstar*(1-np.geomspace(0.3, 1e-4, 40))])
rows = []; xv = []
for x0 in xs:
    r = row(x0)
    if r is None: print("oval failed at", x0); continue
    rows.append(r); xv.append(x0)
B = np.array(rows); xv = np.array(xv); print("ovals", len(xv), "closest to boundary", xv.max(), flush=True)
Bs = B/(xv**2)[:, None]; colmax = np.max(np.abs(Bs), axis=0); keep = colmax > 1e-7*colmax.max()
U, s, Vt = np.linalg.svd(Bs[:, keep]/colmax[keep], full_matrices=False); dim = int(np.sum(s > 1e-6*s[0]))
print("span singular values", np.round(s[:6], 6), "dim", dim, flush=True)
# map: a direction v in the 12-dim coefficient space has M(x0) = Bs @ v (rows scaled by x0^2)
rng = np.random.default_rng(5); hist = {}; picks = {2: [], 3: [], 4: []}
for _ in range(300000):
    v = rng.standard_normal(12); v /= np.linalg.norm(v); M = Bs @ v; Mn = M/np.max(np.abs(M))
    sig = np.abs(Mn) > 1e-7; pr = Mn[sig]; xr = xv[sig]
    zi = [i for i in range(len(pr)-1) if pr[i]*pr[i+1] < 0]; z = len(zi); hist[z] = hist.get(z, 0)+1
    if z in picks and len(picks[z]) < 12:
        picks[z].append((v.tolist(), [float(xr[i]) for i in zi]))
print("first-order zero histogram:", dict(sorted(hist.items())), flush=True)
json.dump(dict(hist=hist, picks=picks), open('data/f19_directions.json', 'w'))
# actual cycle counts along picked directions
out = open('data/F19_q4_alien.jsonl', 'a')
for z in (3, 2, 4):
    for v, zeros in picks[z]:
        v = np.array(v)
        for eps in (1e-3, 1e-4, 1e-5):
            c = Q4 + eps*v
            try: nests = sl.evaluate(c)
            except Exception as ex: print("err", ex); continue
            orig = min(nests, key=lambda n: np.hypot(*n['pt']))
            rec = dict(first_order_zeros=z, zero_x0=zeros, eps=eps, dir=v.tolist(), origin_roots=orig['roots'], origin_stab=orig['stab'], origin_edge=orig['redge'], k=orig['k'], all=[(len(n['roots']), ''.join(n['stab'])) for n in nests])
            out.write(json.dumps(rec)+"\n"); out.flush()
            flag = "ALIEN?" if len(orig['roots']) > z else ""
            print(f"M1 zeros={z} at x0={[f'{t:.4f}' for t in zeros]} eps={eps:g}: actual origin cycles {len(orig['roots'])} {''.join(orig['stab'])} radii {[f'{r:.3g}' for r in orig['roots']]} edge {orig['redge']:.3g} {flag}", flush=True)
print("DONE")
