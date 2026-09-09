#!/usr/bin/env python3
"""Type of the second finite singularity (0,1) along the collision curve.
The repository's audit (audit/claude_center_identify.py) records that the second
finite singularity of every Q_4 centre is a strong NODE, not a focus."""
import sympy as S
l = S.symbols('l', real=True)
a = S.sqrt(-(l+2)/2); m, b = 5*a, 3*l+5
J = S.Matrix([[m, 1], [1 + b, 0]])          # Jacobian at (0,1) in the Shi chart
det = S.simplify(J.det()); tr = S.simplify(J.trace()); disc = S.simplify(tr**2 - 4*det)
print("at (0,1):  det =", det, "  trace =", tr)
print("discriminant tr^2-4det =", S.factor(S.simplify(disc)))
print("\n  l      det      trace    tr^2-4det   type")
for lv in (-3, -4, -5, -7, -10, -20, -50):
    d = float(det.subs(l, lv)); t = float(tr.subs(l, lv)); q = float(disc.subs(l, lv))
    kind = "saddle" if d < 0 else ("node" if q > 0 else "focus")
    print("  %-6d %-8.4g %-8.4g %-11.4g %s%s" % (lv, d, t, q, kind,
          " (unstable)" if t > 0 else " (stable)"))
print("""
A node cannot carry a limit cycle by a Hopf bifurcation, so along this curve the
second nest is NOT available at first order -- the same signature the repository
recorded for Q_4 centres.""")
