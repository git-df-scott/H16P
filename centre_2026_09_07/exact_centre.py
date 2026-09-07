#!/usr/bin/env python3
"""EXACT centre certificate for the Cherkas normal form at a = -2.

    xdot = 1 + x y,  ydot = a00 + a10 x + a20 x^2 + a01 y + a11 x y + a y^2

    a = -2, a11 = k, a01 = -k-3, a20 = -2k^2/27, a10 = 5k(k+3)/27,
    a00 = -(k^2+5k+9)/9.

Claimed Darboux structure:
    U = y + (k+3)/6 - k x/9
    V = -1 + x y + (k+3)x/3 - 2k x^2/9
    C = k x/3 - y - (k+3)/3
with X(U) = 2 C U, X(V) = C V, so H = U/V^2 is a first integral.
"""
import sympy as S
x, y, k = S.symbols('x y k')

a   = S.Integer(-2)
a11 = k
a01 = -k - 3
a20 = -2*k**2/27
a10 = 5*k*(k+3)/27
a00 = -(k**2 + 5*k + 9)/9

print("pinning check  a00 == a01+a11-a10-a20-a :",
      S.simplify(a00 - (a01 + a11 - a10 - a20 - a)) == 0)

P = 1 + x*y
Q = a00 + a10*x + a20*x**2 + a01*y + a11*x*y + a*y**2
Xop = lambda f: S.expand(S.diff(f, x)*P + S.diff(f, y)*Q)

U = y + (k+3)/6 - k*x/9
V = -1 + x*y + (k+3)*x/3 - 2*k*x**2/9
C = k*x/3 - y - (k+3)/3

print("X(U) - 2 C U  =", S.simplify(S.expand(Xop(U) - 2*C*U)))
print("X(V) -   C V  =", S.simplify(S.expand(Xop(V) - C*V)))
H = U/V**2
print("X(H)          =", S.simplify(Xop(H)))
print()
A = {x: 1, y: -1}
print("V(A)           =", S.simplify(V.subs(A)), "   expected (k-9)/9 ->",
      S.simplify(V.subs(A) - (k-9)/9) == 0)
J = S.Matrix([[S.diff(P,x), S.diff(P,y)],[S.diff(Q,x), S.diff(Q,y)]]).subs(A)
print("trace J(A)     =", S.simplify(J.trace()))
print("det   J(A)     =", S.factor(S.simplify(J.det())), "   expected (k-3)(9-k)/27 ->",
      S.simplify(J.det() - (k-3)*(9-k)/27) == 0)
gH = [S.simplify(S.diff(H, v).subs(A)) for v in (x, y)]
print("grad H(A)      =", gH)
Hyy = S.simplify(S.diff(H, y, 2).subs(A))
print("H_yy(A)        =", S.simplify(Hyy), "   expected 729/(9-k)^3 ->",
      S.simplify(Hyy - 729/(9-k)**3) == 0)
Hess = S.Matrix([[S.diff(H, v1, v2) for v2 in (x, y)] for v1 in (x, y)]).subs(A)
dH = S.factor(S.simplify(Hess.det()))
print("det Hess H(A)  =", dH, "   expected 19683(k-3)/(9-k)^5 ->",
      S.simplify(Hess.det() - 19683*(k-3)/(9-k)**5) == 0)
print()
print("=> nondegenerate centre for 3 < k < 9 (det J > 0 and det Hess H != 0).")
