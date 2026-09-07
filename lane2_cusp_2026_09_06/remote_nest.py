#!/usr/bin/env python3
"""Phase 3: does a SECOND antisaddle (the remote nest for the fifth cycle) exist
and persist along the cusp curves?

Equilibria of  xdot = 1+xy,  ydot = a00+a10 x+a20 x^2+a01 y+a11 xy+a y^2
lie on y = -1/x, and x solves
    a20 x^4 + a10 x^3 + (a00 - a11) x^2 - a01 x + a = 0,
with x = 1 the cusp focus A = (1,-1)."""
import json, glob, os
import numpy as np

def a00_of(a, a20, a11, a01, a10): return a01 + a11 - a10 - a20 - a

def equilibria(a, a20, a11, a01, a10):
    a00 = a00_of(a, a20, a11, a01, a10)
    c = [a20, a10, a00 - a11, -a01, a]
    rs = np.roots(c) if abs(a20) > 1e-14 else np.roots(c[1:])
    out = []
    for r in rs:
        if abs(r.imag) > 1e-9*max(1.0, abs(r.real)): continue
        x = r.real
        if abs(x) < 1e-12: continue
        y = -1.0/x
        # Jacobian of (1+xy, a00+a10x+a20x^2+a01y+a11xy+a y^2)
        J = np.array([[y, x],
                      [a10 + 2*a20*x + a11*y, a01 + a11*x + 2*a*y]])
        det, tr = np.linalg.det(J), np.trace(J)
        kind = ("saddle" if det < 0 else
                ("focus" if tr*tr - 4*det < 0 else "node"))
        out.append((x, y, det, tr, kind))
    return sorted(out, key=lambda t: t[0])

print(" ledger                      idx  A=(1,-1)   other finite equilibria (x, kind, trace)")
for f in ['ledger/cusp_row5.jsonl', 'ledger/cusp_row6.jsonl', 'ledger/cusp_row7.jsonl',
          'ledger/cusp_row8.jsonl', 'ledger_grid/cusp_c_am2p5_o0p08.jsonl',
          'ledger_grid/cusp_c_am2p9_o1p0.jsonl', 'ledger_grid/cusp_c_am2p0_om0p08.jsonl']:
    recs = [json.loads(l) for l in open(f) if l.strip()]
    for i in (0, len(recs)//2, len(recs)-1):
        R = recs[i]
        g = lambda k: float(R[k])
        eqs = equilibria(g('a'), g('a20'), g('a11'), g('a01'), g('a10'))
        others = [e for e in eqs if abs(e[0]-1.0) > 1e-6]
        foci = [e for e in others if e[4] == "focus"]
        desc = ", ".join("x=%+.4f %s tr=%+.3g" % (e[0], e[4], e[3]) for e in others) or "none"
        print("  %-27s %-4d           %s%s"
              % (os.path.basename(f), i, desc,
                 "   <== SECOND FOCUS" if foci else ""))
