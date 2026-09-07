#!/usr/bin/env python3
"""Exact centre conditions for the Cherkas family at the Newton limit.

  xdot = 1 + x y,   ydot = a00 + a10 x + a20 x^2 + a01 y + a11 x y + a y^2,
  a00 = a01 + a11 - a10 - a20 - a,   marked equilibrium A = (1,-1).

Translate to A, normalise the linear part, and compute the focal values by the
formal first integral.  Then test the converged point."""
import sympy as S

x, y, a, a20, a11, a01, a10 = S.symbols('x y a a20 a11 a01 a10')
a00 = a01 + a11 - a10 - a20 - a
P = 1 + x*y
Q = a00 + a10*x + a20*x**2 + a01*y + a11*x*y + a*y**2
# translate to A=(1,-1)
u, v = S.symbols('u v')
Pt = S.expand(P.subs({x: u+1, y: v-1}))
Qt = S.expand(Q.subs({x: u+1, y: v-1}))
print("at A: P =", Pt, "   Q =", S.simplify(Qt))
J = S.Matrix([[S.diff(Pt, u), S.diff(Pt, v)], [S.diff(Qt, u), S.diff(Qt, v)]]).subs({u: 0, v: 0})
T = S.simplify(J.trace()); Dm = S.simplify(J.det())
print("trace V1 =", T, "    det L =", Dm)

# impose V1 = 0 and a = -2 (the event's shape)
sub0 = {a: -2}
V1 = T.subs(sub0)
a01s = S.solve(V1, a01)[0]
print("\nat a=-2, V1=0 gives a01 =", a01s)
Pt2 = S.expand(Pt.subs(sub0).subs(a01, a01s))
Qt2 = S.expand(Qt.subs(sub0).subs(a01, a01s))
Jn = S.Matrix([[S.diff(Pt2,u), S.diff(Pt2,v)],[S.diff(Qt2,u), S.diff(Qt2,v)]]).subs({u:0,v:0})
w = S.simplify(S.sqrt(S.simplify(Jn.det())))
print("linear part det =", S.simplify(Jn.det()), "  w =", w)

# put the linear part into (-w v, w u) form:  u = U,  v = (Jn[0,0] U + w V)/(-Jn[0,1])
# generic normalisation via the eigen-structure
A11, A12, A21, A22 = Jn[0,0], Jn[0,1], Jn[1,0], Jn[1,1]
U, V = S.symbols('U V')
sub = {u: U, v: (-A11*U + w*V)/A12}
Pn = S.expand(Pt2.subs(sub, simultaneous=True))
Qn = S.expand(Qt2.subs(sub, simultaneous=True))
# dU/dt = Pn ; dV/dt = (A12*Qn + A11*Pn)/ (w*A12) * ... derive:
Udot = Pn
Vdot = S.simplify((A12*Qn + A11*Pn)/w)
print("\nnormalised linear part check:")
print("  dU/dt linear:", S.expand(Udot).subs({U:0,V:0}), S.diff(Udot,U).subs({U:0,V:0}), S.diff(Udot,V).subs({U:0,V:0}))
print("  dV/dt linear:", S.expand(Vdot).subs({U:0,V:0}), S.diff(Vdot,U).subs({U:0,V:0}), S.diff(Vdot,V).subs({U:0,V:0}))
