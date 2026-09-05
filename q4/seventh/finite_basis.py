#!/usr/bin/env python3
"""Exact forcing and center-data checks for the closed finite-a subspace."""
import json
from pathlib import Path
import sympy as S

a,t,F,G=S.symbols('a t F G'); c=S.Rational(5,36)
def dt(z):
    return S.diff(z,t)+S.diff(z,F)*G+S.diff(z,G)*((2*t-1)*G+c*F)/(t*(1-t))
def ell(h):
    return S.Rational(9,3080)*(h[0]+S.Rational(144,221)*h[1]+S.Rational(11,6)*h[2]+S.Rational(204,221)*h[3])
j0=S.Rational(36,5)*t*(1-t)*G
j1=(t*t*(1-t)*G-t*(1-t)*F+j0)/(2+c)
j2=(t**3*(1-t)*G-2*t*t*(1-t)*F+4*j1)/(6+c)
K=S.Matrix([j1,j2,6*j0-11*j1-6*t*(1-t)*F,12*j1-17*j2-6*t*t*(1-t)*F])
coeffs=[
 S.Matrix([-S.Rational(176,3)*(9*a+4),1088*a,S.Rational(16,3)*(54*a+59),-768*a]),
 S.Matrix([-96*(12*a-5),S.Rational(272,3)*(21*a-10),192,-S.Rational(16,3)*(18*a+49)]),
 S.Matrix([-S.Rational(40,9)*(36*a-25),S.Rational(680,9)*(3*a-2),S.Rational(80,3),S.Rational(40,3)*(a-5)])]
functions=[F,t*F,t*(1-t)*G]
centers=[(1,c),(0,1),(0,c)]
for f,h,(y0,y1) in zip(functions,coeffs,centers):
    L=(1-a*t)*(1-t)*dt(dt(f))-(1-a)*dt(f)/2+c*a*f
    assert S.factor(L+K.dot(h)/(1152*t*t*(1-t)))==0
    assert S.factor(ell(h)-y0)==0
    assert S.factor(-S.Rational(3,2)*(1+a)*ell(h)+h[2]/192-y1)==0
P=486*a**3-441*a**2-486*a+236
basisdet=S.factor(S.Matrix.hstack(*coeffs,S.Matrix([1,0,0,0])).det())
assert S.factor(basisdet-S.Rational(2263040,81)*P)==0
assert P.subs(a,1)==-205
Path(__file__).with_suffix('.json').write_text(json.dumps({
 'status':'EXACT_IDENTITIES',
 'scope':'Three-dimensional subspace only; not the complete finite-a family',
 'functions':list(map(str,functions)),
 'coefficients':[list(map(str,h)) for h in coeffs],
 'basis_determinant':str(basisdet)},indent=2)+'\n')
print('Exact finite-a forcing, six center data, and basis determinant: PASS')
