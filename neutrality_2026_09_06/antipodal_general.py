#!/usr/bin/env python3
"""General quadratic field: at an antipodal pair of infinite singularities the
linearisations are exact negatives, so the hyperbolicity ratios are reciprocal.

Sphere field (rescaled by s3^{d-1}=s3, orientation preserved on s3>0):
    F(s) = ( Pbar - s1 W , Qbar - s2 W , -s3 W ),  W = s1 Pbar + s2 Qbar,
    Pbar = p00 s3^2 + (p10 s1 + p01 s2) s3 + (p20 s1^2 + p11 s1 s2 + p02 s2^2).
For d=2, F is EVEN: F(-s)=F(s).  The antipodal map iota has d(iota)=-I, so
    DF(-s) = -DF(s):
iota conjugates the flow near s to the TIME-REVERSED flow near -s.
"""
import sympy as S

s1, s2, s3, u = S.symbols('s1 s2 s3 u', real=True)
ps = S.symbols('p00 p10 p01 p20 p11 p02', real=True)
qs = S.symbols('q00 q10 q01 q20 q11 q02', real=True)
p00, p10, p01, p20, p11, p02 = ps
q00, q10, q01, q20, q11, q02 = qs

Pb = p00*s3**2 + (p10*s1 + p01*s2)*s3 + (p20*s1**2 + p11*s1*s2 + p02*s2**2)
Qb = q00*s3**2 + (q10*s1 + q01*s2)*s3 + (q20*s1**2 + q11*s1*s2 + q02*s2**2)
W = s1*Pb + s2*Qb
F = S.Matrix([Pb - s1*W, Qb - s2*W, -s3*W])
sub = {s1: -s1, s2: -s2, s3: -s3}

print("F(-s) - F(s)   =", S.simplify(F.subs(sub, simultaneous=True) - F).T)
J = F.jacobian([s1, s2, s3])
print("DF(-s) + DF(s) =", S.simplify(J.subs(sub, simultaneous=True) + J))
print()
print("Both vanish identically for EVERY quadratic field, so the eigenvalues at")
print("an infinite singularity and at its antipode are exact negatives.\n")

G = S.expand((q20 + q11*u + q02*u**2) - u*(p20 + p11*u + p02*u**2))
print("infinite directions: G(u)=0,  G =", S.collect(G, u))
print("lam_equator(u)    =", S.collect(S.expand(S.diff(G, u)), u))
print("lam_transverse(u) =", S.collect(S.expand(-(p20 + p11*u + p02*u**2)), u))
print("""
Consequence for the boundary graphic of a focus nest that runs
   A --(through the plane)--> B --(equator arc)--> A
with A,B an antipodal pair of infinite saddles and no other singularity on it:
it enters A along the equator and leaves A transversally, and enters B
transversally and leaves B along the equator, so

   r_A = |lam_eq / lam_tr|(u*),   r_B = |lam_tr / lam_eq|(u*) = 1 / r_A,
   r(Gamma) = r_A * r_B == 1     for every parameter value.

The first stability coefficient of such a graphic is therefore identically
neutral: it is not a curve in parameter space, it is everything.
""")
