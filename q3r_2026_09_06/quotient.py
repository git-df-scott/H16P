#!/usr/bin/env python3
"""Phase 1.4: is the corrected nine-dimensional Melnikov kernel exactly the
span of the geometrically trivial directions?

Trivial directions for the reversible centre X0 = (P0, Q0):
  * infinitesimal affine coordinate changes: delta X = -[V, X0] for affine V
    (6 parameters: 2 translations + 4 linear),
  * time rescaling: delta X = X0                       (1),
  * motion within the centre family: dX0/da, dX0/db    (2).
Total 9 candidate directions.  Each must have M_1 == 0.
"""
import sympy as S
x, y, a, b = S.symbols('x y a b')
P0 = (b-2)/4 + (1-b)*y + a*x**2 + b*y**2
Q0 = -2*x*y
X0 = S.Matrix([P0, Q0])

def bracket(V):
    """[V, X0] = DX0 . V - DV . X0"""
    DX = S.Matrix([[S.diff(P0, x), S.diff(P0, y)], [S.diff(Q0, x), S.diff(Q0, y)]])
    DV = S.Matrix([[S.diff(V[0], x), S.diff(V[0], y)], [S.diff(V[1], x), S.diff(V[1], y)]])
    return S.simplify(DX*V - DV*X0)

al = S.symbols('al0:6')
dirs = []
names = []
for i, V in enumerate([S.Matrix([1, 0]), S.Matrix([0, 1]), S.Matrix([x, 0]),
                       S.Matrix([y, 0]), S.Matrix([0, x]), S.Matrix([0, y])]):
    dirs.append(-bracket(V)); names.append("affine%d" % i)
dirs.append(X0); names.append("time")
dirs.append(S.Matrix([S.diff(P0, a), S.diff(Q0, a)])); names.append("d/da")
dirs.append(S.Matrix([S.diff(P0, b), S.diff(Q0, b)])); names.append("d/db")

# coefficient extraction in the monomial basis of a quadratic perturbation
mons = [S.Integer(1), x, y, x**2, x*y, y**2]
def coeffs(v):
    row = []
    for comp in (S.expand(v[0]), S.expand(v[1])):
        p = S.Poly(comp, x, y)
        row += [p.coeff_monomial(mn) for mn in mons]
    return row     # [p00,p10,p01,p20,p11,p02, q00,q10,q01,q20,q11,q02]

Mtx = S.Matrix([coeffs(d) for d in dirs])
print("trivial-direction matrix is %d x %d" % (Mtx.rows, Mtx.cols))
print("rank over Q(a,b):", Mtx.rank())

# the corrected kernel conditions (visible coefficients only)
idx = {n: k for k, n in enumerate(
    ['p00','p10','p01','p20','p11','p02','q00','q10','q01','q20','q11','q02'])}
print("\nchecking each trivial direction against the corrected kernel conditions")
print("  q00 = (b-2) q20 / (4(a+2))")
print("  p10 = -a q01 - (a-1)(b-1) q20/(a+2)")
print("  p11 = -(a+1) q02 + (a-1) b q20/(a+2)")
ok = True
for nm, d in zip(names, dirs):
    c = coeffs(d)
    q00, q01, q02, q20 = c[idx['q00']], c[idx['q00']+2], c[idx['q02']], c[idx['q20']]
    p10, p11 = c[idx['p10']], c[idx['p11']]
    r1 = S.simplify(q00 - (b-2)*q20/(4*(a+2)))
    r2 = S.simplify(p10 + a*q01 + (a-1)*(b-1)*q20/(a+2))
    r3 = S.simplify(p11 + (a+1)*q02 - (a-1)*b*q20/(a+2))
    good = (r1 == 0 and r2 == 0 and r3 == 0)
    ok &= good
    print("  %-9s residuals (%s, %s, %s)  %s" % (nm, r1, r2, r3, "in kernel" if good else "NOT IN KERNEL"))
print("\nall trivial directions in the corrected kernel:", ok)
print("kernel dimension 9 (3 conditions on 6 visible + 6 invisible).")
print("trivial span rank:", Mtx.rank())
