#!/usr/bin/env python3
"""Lane C: locate the homoclinic loop surrounding the third-order weak focus
on the stratum m=5a, b=3l+5 by bisection in l at fixed a, using the signed
separatrix offset; report saddle quantity (trace), eta_3, and the type of
the second antisaddle (0,1)."""
import sys, os, numpy as np
sys.path.insert(0, os.path.dirname(__file__))
from claude_laneC_splitting import splitting, data
def eta3(l, a): return -25*a*(2*a*a+l+2)*(5*a*a*l+6*a*a-3*l**3-12*l*l-15*l-6)/64
def offset(l, a):
    r = splitting(l, a)
    return None if r is None else r[2]
for a in (1.3, 1.5, 2.0):
    lo, hi = -1-np.sqrt(1+3*a*a)+0.05, -1.0
    # march to find a sign change
    grid = np.linspace(lo, hi, 25); vals = [offset(l, a) for l in grid]
    br = None
    for i in range(len(grid)-1):
        if vals[i] is not None and vals[i+1] is not None and vals[i]*vals[i+1] < 0: br = (grid[i], grid[i+1]); break
    if br is None: print(f"a={a}: no sign change on grid; offsets={[None if v is None else round(v,3) for v in vals]}"); continue
    l1, l2 = br
    for _ in range(40):
        mid = (l1+l2)/2; v = offset(mid, a)
        if v is None: break
        if v*offset(l1, a) < 0: l2 = mid
        else: l1 = mid
    lstar = (l1+l2)/2
    f, sad = data(lstar, a)
    pt, Jn = sad[0]
    print(f"a={a}: loop at l*={lstar:.10f}  m={5*a} b={3*lstar+5:.6f}  saddle {np.round(pt,5)} trace={np.trace(Jn):+.5f} eta3={eta3(lstar,a):+.4f}")
    # type of (0,1): Jacobian there
    m, b = 5*a, 3*lstar+5
    J01 = np.array([[m*1, -1+2*1], [1+b*1, b*0]])  # dP/dx = 2l x + m y = m ; dP/dy = -1 + m x + 2y = 1 ; dQ/dx = 1+2ax+by = 1+b ; dQ/dy = b x = 0
    print(f"       (0,1): det={np.linalg.det(J01):+.4f} trace={np.trace(J01):+.4f} -> {'saddle' if np.linalg.det(J01)<0 else 'antisaddle'}")
