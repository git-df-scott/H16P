#!/usr/bin/env python3
"""Route 4b: hemicycle connection on the order-3 stratum (two-foci region).
Infinite saddles N (x->+inf) and S (x->-inf) in direction (1,u0). N's finite
separatrix is unstable (lam_z>0), S's is stable. Integrate N's unstable
separatrix forward from far out and S's stable separatrix backward from far
out; record where each crosses the section x=0; the difference is the
splitting D(l,a). D=0 <=> hemicycle (N -> S through the finite plane) exists.
Also record whether the connecting curve passes between (0,0) and (0,1)."""
import numpy as np
from scipy.integrate import solve_ivp
def setup(l, a):
    m, b = 5*a, 3*l+5
    f = lambda t, u: [-u[1]+l*u[0]**2+m*u[0]*u[1]+u[1]**2, u[0]+a*u[0]**2+b*u[0]*u[1]]
    g = np.poly1d([-1, -5*a, 2*l+5, a]); roots = [r.real for r in g.roots if abs(r.imag) < 1e-10]
    return f, roots
def separatrix_crossing(l, a, R=400.0, sign=+1):
    """sign=+1: N (x->+inf) unstable branch forward; sign=-1: S stable branch backward."""
    f, roots = setup(l, a)
    if len(roots) != 1: return None
    u0 = roots[0]
    # in the chart u=y/x, z=1/x the separatrix is tangent to the z-eigendirection at (u0,0):
    # for quadratic systems the z-eigenvector of the linearization is generically not (0,1); compute it.
    m, b = 5*a, 3*l+5
    P2 = lambda u: l+m*u+u*u; P1 = lambda u: -u; Q1 = lambda u: 1.0; Q2 = lambda u: a+b*u
    # chart equations: u' = (Q2-uP2) + z(Q1-uP1) ; z' = -z(P2 + z P1)
    gu = np.polyder(np.poly1d([-1, -5*a, 2*l+5, a]))
    J = np.array([[np.polyval(gu, u0), Q1(u0)-u0*P1(u0)], [0.0, -P2(u0)]])
    w, V = np.linalg.eig(J)
    iz = int(np.argmax(np.abs(V[1, :])))          # eigenvector with z-component: the finite-plane separatrix
    vec = V[:, iz]/V[1, iz]                        # normalize z-component to 1
    z0 = 1.0/R; du = vec[0]*z0
    u_start = u0+du
    x0 = sign*R; y0 = u_start*x0                   # (x,y) = (1/z, u/z); antipode: x->-inf
    tdir = 1.0 if sign > 0 else -1.0
    def ev(t, s): return s[0]
    ev.terminal = True; ev.direction = 0
    sol = solve_ivp(f, (0, tdir*400), [x0, y0], rtol=1e-10, atol=1e-12, max_step=0.05, events=ev)
    if len(sol.t_events[0]) == 0: return None
    return sol.y_events[0][0][1], sol
def splitting(l, a):
    cN = separatrix_crossing(l, a, sign=+1); cS = separatrix_crossing(l, a, sign=-1)
    if cN is None or cS is None: return None
    return cN[0], cS[0]
for a in (0.5, 1.0, 1.5, 2.0):
    print(f"a={a}:")
    for l in (-30, -20, -15, -12, -10, -8, -6):
        if 3*a*a > l*l+2*l: continue
        r = splitting(l, a)
        if r is None: print(f"   l={l}: no crossing of x=0"); continue
        yN, yS = r
        print(f"   l={l}: N-branch crosses x=0 at y={yN:+.5f};  S-branch at y={yS:+.5f};  D={yN-yS:+.5f}")
