#!/usr/bin/env python3
"""Route 4b: Poincare compactification of the order-3 stratum system
 x' = -y + l x^2 + 5a x y + y^2 ,  y' = x + a x^2 + (3l+5) x y .
Infinite singular points: real roots t of  Q2(1,t) - t P2(1,t) = 0, i.e.
  -t^3 - 5a t^2 + (2l+5) t + a = 0   (direction (1,t)).
In the chart (u=y/x, z=1/x): u' = Q2(1,u) - u P2(1,u) + z(Q1 - u P1)...,
z' = -z (P2(1,u) + z P1(1,u)). At an infinite singular point (u0,0):
eigenvalues  lam_u = d/du[Q2(1,u)-uP2(1,u)](u0),  lam_z = -P2(1,u0).
Hyperbolicity ratio of a saddle: r = |lam_in/lam_out| (contracting over
expanding). For a hemicycle through two hyperbolic saddles at infinity the
first-order stability is governed by the product of ratios; neutral when
the product is one. We tabulate the infinite singular points, their types
and ratios on the stratum in the two-foci region."""
import numpy as np
def infinity_data(l, a):
    m, b = 5*a, 3*l+5
    P2 = lambda u: l + m*u + u*u          # P2(1,u)
    Q2 = lambda u: a + b*u                # Q2(1,u)
    g = np.poly1d([-1, -5*a, 2*l+5, a])   # Q2 - u P2
    roots = g.roots
    out = []
    for u0 in roots:
        if abs(u0.imag) > 1e-10: continue
        u0 = u0.real
        lam_u = np.polyval(np.polyder(g), u0)
        lam_z = -P2(u0)
        kind = "saddle" if lam_u*lam_z < 0 else "node/focus-type"
        ratio = None
        if kind == "saddle":
            lin, lout = (abs(lam_u), abs(lam_z)) if lam_u < 0 else (abs(lam_z), abs(lam_u))
            ratio = lin/lout
        out.append((u0, lam_u, lam_z, kind, ratio))
    return out
print("Shi seed l=-10,a=1:")
for row in infinity_data(-10, 1): print("   u0=%.5f lam_u=%.4f lam_z=%.4f %s ratio=%s" % (row[0], row[1], row[2], row[3], None if row[4] is None else round(row[4],5)))
print("\nStratum scan (two-foci region, 3a^2<l^2+2l): number of real infinite directions, saddle ratios")
for a in (0.5, 1.0, 1.5, 2.0):
    for l in (-20, -15, -12, -10, -8, -6, -4):
        if 3*a*a > l*l+2*l: continue
        rows = infinity_data(l, a)
        desc = "; ".join("u0=%+.3f %s r=%s" % (r[0], r[3][:6], "-" if r[4] is None else "%.4f" % r[4]) for r in rows)
        print(f"  a={a} l={l}: {len(rows)} real dir(s): {desc}")
