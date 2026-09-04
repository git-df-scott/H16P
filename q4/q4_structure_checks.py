#!/usr/bin/env python3
"""Bounded exact Q4 structure checks; no search, scan, or optimization.

Run with the repository requirements plus sympy. This process enforces one
numerical thread and a ten-second CPU ceiling before importing libraries.
"""
import os
import resource

for key in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS",
            "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ[key] = "1"
resource.setrlimit(resource.RLIMIT_CPU, (10, 10))

import sympy as S
import mpmath as mp

k, t = S.symbols("k t", positive=True)
d = k - 1
m1, o2, o3, m4 = S.symbols("m1 o2 o3 m4")
m2 = -2*o2/S.Integer(3)-2*d*o3/(3*k)
m3 = -2*o3/(3*k)
a1 = (128*(-117-265*k+30*k*k)/(27*k*k)*m1
      +16*(-174+5*k)/(3*k)*m2-16*(243+121*k)/(3*k)*m3
      +64*d*(-119+60*k)/(9*k)*m4)
a2 = (1088*(21+k)/(27*k*k)*m1+544/k*m2+1088/k*m3
      +1088*d/(9*k)*m4)
b0 = (-256*d/(3*k)*m1-32*(-27+25*k+15*k*k)/(3*k)*m2
      -16*d*(54+31*k)/(3*k)*m3-64*d*(-18+5*k)/(3*k)*m4)
b1 = (-64*d*(18+77*k)/(27*k*k)*m1-16*(-111+137*k)/(3*k)*m2
      -768*d/k*m3-64*d*(-116+77*k)/(9*k)*m4)
matrix = S.linear_eq_to_matrix([a1,a2,b0,b1], [m1,o2,o3,m4])[0]
determinant = S.factor(matrix.det(method="domain-ge"))
assert S.factor(determinant - S.Rational(17843617792000,6561)*d**2/k**4) == 0
print("exact determinant mu -> (alpha1,alpha2,beta0,beta1):", determinant)

# Universal period pair: F solves the Gauss equation and K is its companion.
f = sum(S.rf(S.Rational(1,6),n)*S.rf(S.Rational(5,6),n)
        /S.factorial(n)**2*t**n for n in range(7))
companion = S.series((1-t)*(f+6*t*S.diff(f,t)),t,0,6).removeO()
assert S.series(t*(1-t)*S.diff(f,t,2)+(1-2*t)*S.diff(f,t)
                -S.Rational(5,36)*f,t,0,6).removeO() == 0
assert S.series(6*t*S.diff(companion,t)-companion+f,t,0,5).removeO() == 0
w = S.series(companion/f,t,0,5).removeO()
moment = S.series((1-w)/t,t,0,4).removeO()
assert moment.coeff(t,0) == S.Rational(1,6)
assert moment.coeff(t,1) == S.Rational(25,432)
assert moment.coeff(t,2) == S.Rational(775,23328)
assert moment.coeff(t,1)/moment.coeff(t,2) == S.Rational(54,31)
ww4 = S.factor(S.diff(w,t,2)*S.diff(-moment,t,3)-S.diff(w,t,3)*S.diff(-moment,t,2))
assert ww4.subs(t,0) == S.Rational(3705625,1451188224)
print("exact auxiliary W4 at center:", ww4.subs(t,0))
print("w center Taylor polynomial:", w)
print("companion center coefficient of t^4:", companion.coeff(t,4))

# Symbolic universal coordinate change, with beta1 normalized to one.
aa1, aa2, bb0, ww = S.symbols("alpha1 alpha2 beta0 w")
ss = k-d*t
g = aa2*ss**2+aa1*ss+bb0-k-k*aa1-k*k*aa2+(ss-bb0)*ww
A, B, b = -(aa1+2*k*aa2), d*aa2, (k-bb0)/d
g_universal = A*t+B*t*t+b*(ww-1)-t*ww
assert S.factor(g-d*g_universal) == 0
assert S.simplify((k-(54-23*k)/31)/d) == S.Rational(54,31)
print("universal g identity and corrected threshold: exact pass")

# Three specified diagnostic points, not a parameter scan.
mp.mp.dps = 40
for tt in (mp.mpf(1)/4,mp.mpf(1)/2,mp.mpf(3)/4):
    FF = mp.hyp2f1(mp.mpf(1)/6,mp.mpf(5)/6,1,tt)
    FFp = mp.mpf(5)/36*mp.hyp2f1(mp.mpf(7)/6,mp.mpf(11)/6,2,tt)
    KK = mp.hyp2f1(-mp.mpf(1)/6,mp.mpf(1)/6,1,tt)
    residual = abs(KK-(1-tt)*(FF+6*tt*FFp))
    assert residual < mp.mpf("1e-35")
    print("period companion residual at t="+str(tt)+":",mp.nstr(residual,4))
print("All bounded structure checks passed.")
