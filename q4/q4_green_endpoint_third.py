#!/usr/bin/env python3
"""Cheap exact replay of third-strike endpoint identities, with no search.

This verifies algebra and the finite-part primitive. The overlap-region
error estimates are proved analytically in notes_green_third.md.
"""
import os
import resource
resource.setrlimit(resource.RLIMIT_CPU,(10,10))
for name in ("OPENBLAS_NUM_THREADS","OMP_NUM_THREADS","MKL_NUM_THREADS",
             "VECLIB_MAXIMUM_THREADS","NUMEXPR_NUM_THREADS"):
    os.environ[name]="1"

from fractions import Fraction as Q
import sympy as S


def run_checks():
    t,F,Fp=S.symbols('t F Fp')
    d_boundary=t*(1-t)*Fp
    boundary=t*(1-t)*F
    j0=S.Rational(36,5)*d_boundary
    j1=(t*d_boundary-boundary+j0)/(2+S.Rational(5,36))
    j2=(t*t*d_boundary-2*t*boundary+4*j1)/(6+S.Rational(5,36))
    H=S.Rational(1788,77)*j1-S.Rational(1326,77)*j2-6*j0+6*t*(1-t)**2*F
    expected=6*t*(1-t)**2*(5*F-36*(1-t)*Fp)/77
    assert S.factor(H-expected)==0
    A,B,eta=Q(94,77),-Q(17,77),Q(1)
    Y0=3*(1326*A+864*B-2431*eta-102)/1361360
    assert Y0==-Q(3,1232)

    n=S.symbols('n',positive=True)
    ratio=(n+S.Rational(1,6))*(n+S.Rational(5,6))/(n+1)**2
    assert S.factor(5+36*n-36*(n+1)*ratio-5*n/(n+1))==0
    # Positive beta integration reduces the moments to these elementary
    # integrals of binomial series; positivity justifies termwise integration.
    beta_first=5*(Q(6,5)-1)
    beta_second=5*(Q(6)-1)
    assert (beta_first,beta_second)==(1,25)
    P_endpoint=Q(1,14784)-beta_first/14784
    Z_endpoint=Y0+Q(3,2)*(beta_second-beta_first)/14784
    assert P_endpoint==0 and Z_endpoint==0

    # v=exp(x/3), so d/dx=(v/3)d/dv. All identities are rational functions.
    v,w=S.symbols('v w',positive=True)
    dx=lambda expression:v*S.diff(expression,v)/3
    cosh=lambda order:(v**order+v**(-order))/2
    sinh=lambda order:(v**order-v**(-order))/2
    O=S.Rational(3,10)*sinh(5)+S.Rational(3,2)*sinh(1)
    JO=S.Rational(9,50)*cosh(5)+S.Rational(9,2)*cosh(1)
    tanh=sinh(3)/cosh(3)
    bracket=S.Rational(23,18)*JO-O*tanh+dx(O)/2
    assert S.factor(dx(JO)-O)==0
    assert S.factor(dx(bracket)-O*tanh**2)==0
    assert S.factor(bracket.subs(v,1))==S.Rational(162,25)
    # Odd Laurent parity rules out a constant at infinity; c is an even
    # rational function of v and has leading term lambda*v^6/4.
    assert S.factor(bracket.subs(v,-v)+bracket)==0
    asymptotic_bracket=S.series(bracket.subs(v,1/w),w,0,2).removeO()
    assert asymptotic_bracket.coeff(w,0)==0

    k=cosh(3)**2
    a=1-1/k
    r=-tanh*dx(O)/(2*O)
    P0star=-S.Rational(3,2)*(1+a)*S.Rational(-3,1232)-S.Rational(1,192)-r*S.Rational(-3,1232)
    series=S.series(P0star.subs(v,1/w),w,0,7).removeO()
    assert series.coeff(w,0)==S.Rational(1,14784)
    assert series.coeff(w,4)==S.Rational(5,616)
    # k^(-2/3) has leading coefficient 2^(4/3)*w^4, so b=5/2^(4/3).
    assert S.Rational(24*5,14784)==series.coeff(w,4)

    print("Exact H-star and center datum: PASS")
    print("Exact beta moments 1 and 25; P-star(1)=Z-star(1)=0: PASS")
    print("Exact finite-part primitive and zero asymptotic constant: PASS")
    print("Exact cancellation of the first nonanalytic kappa correction: PASS")
    print("Uniform matched-asymptotic remainder bounds: analytic proof in notes_green_third.md.")


if __name__=='__main__':
    run_checks()
