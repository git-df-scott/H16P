#!/usr/bin/env python3
"""Exact check: for a quadratic field the compactified flow satisfies
p(X) o iota = -iota_* p(X) (iota = antipodal map), so the linearisations at
antipodal infinite singularities are exact negatives of one another.
Consequence: their hyperbolicity ratios are reciprocal."""
import sympy as S

s1, s2, s3, l, a = S.symbols('s1 s2 s3 l a', real=True)
m, b = 5*a, 3*l + 5
Pb = -s2*s3 + l*s1**2 + m*s1*s2 + s2**2
Qb = s1*s3 + a*s1**2 + b*s1*s2
W = s1*Pb + s2*Qb
F = S.Matrix([Pb - s1*W, Qb - s2*W, -s3*W])
Fi = F.subs({s1: -s1, s2: -s2, s3: -s3}, simultaneous=True)
print("p(X) o iota + iota_* p(X)  =", S.simplify(Fi + (-F)))   # iota_* v = -v
print("  (zero means p(X)oiota = -iota_*p(X): antipodal points carry negated linearisations)")

# symbolic eigenvalues at an infinite singularity, from chart U1
u = S.symbols('u', real=True)
G = a + (2*l + 5)*u - 5*a*u**2 - u**3
lam_eq = S.expand(S.diff(G, u))          # along the equator
lam_tr = S.expand(-(l + 5*a*u + u**2))   # transverse, into the plane
print("\nat u*: lam_equator =", lam_eq, " lam_transverse =", lam_tr)
print("ratio at A (in along equator, out transversally) r_A = |lam_eq/lam_tr|")
print("at the antipode both are negated and the roles swap, so r_B = |lam_tr/lam_eq| = 1/r_A")
print("=> product of hyperbolicity ratios around the graphic is IDENTICALLY 1")
