#!/usr/bin/env python3
"""Exact identities and independent quadrature for the resonant H16 strike.

Run from repository root: python resonant/check_resonant.py
Symbolic assertions verify algebra; quadrature is a numerical cross-check.
Neither is an interval certificate for periodic orbits.
"""
import json
import platform
from pathlib import Path

import mpmath as mp
import sympy as S


def main():
    x, y, u, v, a, b, e0, e1, e2, h, s = S.symbols(
        'x y u v a b e0 e1 e2 h s', real=True)
    P = (b-2)/4 + e1*x + (1-b)*y + a*x*x + e2*x*y + b*y*y
    Q = e0 - 2*x*y
    H = (x*x+y*y+S.Rational(1,4))/y
    base = {a:-1,b:1,e0:0,e1:0,e2:0}
    assert S.simplify(S.diff(H,x)*P.subs(base)+S.diff(H,y)*Q.subs(base)) == 0
    weighted_div = S.diff(P/y**2,x)+S.diff(Q/y**2,y)
    assert S.simplify(weighted_div - (2*(a+1)*x/y**2+e1/y**2+e2/y-2*e0/y**3)) == 0

    K = 1 + 2/(a+e1*u+e2*v+(b-2)*u*u/4+(1-b)*u*v+b*v*v)
    resonant_jet = (S.diff(K,u)*S.diff(K,v)+S.diff(K,u,v)).subs({u:0,v:0,a:-1})
    assert S.simplify(resonant_jet-2*(b-1)) == 0
    G2 = 2*S.pi*e1/S.sqrt(2-b-e1**2)
    F3 = S.factor(-G2*resonant_jet)
    assert S.simplify(F3+4*S.pi*(b-1)*e1/S.sqrt(2-b-e1**2)) == 0

    cu = e1+e2/2
    M = 2*S.pi*(h-1)*(cu-2*e0*(h+1))
    hs = 1/s+s/4
    Hy_bottom = 1-4/s**2
    D1 = S.factor(M.subs(h,hs)/Hy_bottom)
    expected = S.pi*e0*(1-s*s/4)-S.pi*cu*s*(2-s)/(2*(2+s))
    assert S.simplify(D1-expected) == 0
    # Mixed derivative in (b, epsilon1), with epsilon2=-2 epsilon1.
    J = 2*S.pi*((h-1)-2*S.log((h+1)/2))
    Dmixed = S.simplify(J.subs(h,hs)/Hy_bottom)
    # Leading endpoint coefficient of the third ideal generator.
    assert S.limit(Dmixed/s,s,0,dir='+') == -S.pi/2
    assert S.limit(Dmixed,s,0,dir='+') == 0

    # The circle-domain derivative integrand, including moving boundary.
    Rb = -y*y+y*S.log(y)+S.Rational(1,4)
    expanded = -1+2*y+S.log(y)/y-2*S.log(y)+1/(4*y*y)-1/(2*y)
    assert S.simplify((1-2*y)*Rb/y**2-expanded) == 0

    mp.mp.dps = 70
    quad = []
    for hv in [mp.mpf('1.02'),mp.mpf('1.5'),mp.mpf(3),mp.mpf(10),mp.mpf(50)]:
        c, r = hv/2, mp.sqrt(hv*hv-1)/2
        yy = lambda theta: c+r*mp.cos(theta)
        # Integrate vertical slices of the disc; dy=-r sin(theta)dtheta.
        values = [mp.quad(lambda t: 2*r*r*mp.sin(t)**2/yy(t)**j,
                          [0,mp.pi/2,mp.pi]) for j in (1,2,3)]
        exact = [mp.pi*(hv-1),2*mp.pi*(hv-1),2*mp.pi*(hv*hv-1)]
        err = max(abs(z-w) for z,w in zip(values,exact))
        mix = mp.quad(lambda t: (1-2*yy(t))/yy(t)**2 *
                      (-yy(t)**2+yy(t)*mp.log(yy(t))+mp.mpf(1)/4),
                      [0,mp.pi/2,mp.pi])
        jval = 2*mp.pi*((hv-1)-2*mp.log((hv+1)/2))
        assert err < mp.mpf('1e-55')
        assert abs(mix-jval) < mp.mpf('1e-55')
        quad.append({'h':str(hv),'moment_max_abs_error':str(err),
                     'mixed_abs_error':str(abs(mix-jval))})

    # Nontrivial two-root compact Melnikov control. This is not a cycle certificate.
    f = lambda z: 1-mp.log1p(z)/z if z else mp.mpf(0)
    F = lambda z: 3-4*(1+z)+20*f(z)
    roots = [mp.findroot(F,(mp.mpf('.1'),mp.mpf('.5'))),
             mp.findroot(F,(mp.mpf(1),mp.mpf(4)))]
    assert 0 < roots[0] < roots[1]
    assert F(0)<0 and F(1)>0 and F(10)<0
    report = {
        'classification':'Exact symbolic identities plus non-rigorous independent numerical checks',
        'versions':{'python':platform.python_version(),'sympy':S.__version__,'mpmath':mp.__version__},
        'symbolic_assertions':'PASS',
        'F3_persistent_resonance':str(F3),
        'D_upper_first_normal':str(S.factor(expected)),
        'D_upper_mixed_generator':str(Dmixed),
        'quadrature':quad,
        'two_compact_Melnikov_control':{
            'p':1,'q':3,'k':20,'u_roots':[str(z) for z in roots],
            'h_roots':[str(1+2*z) for z in roots],
            'signs_at_u_0_1_10':[str(F(z)) for z in (0,1,10)]},
        'limitations':['No interval cycle certificate.','No global bound for H(2).',
                       'Pure endpoint cyclicity is not computed by these checks.']}
    path = Path(__file__).parent/'data'/'exact_checks.json'
    path.write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,indent=2))


if __name__ == '__main__':
    main()
