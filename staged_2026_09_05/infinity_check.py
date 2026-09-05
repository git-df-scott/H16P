#!/usr/bin/env python3
"""Exact algebra checks for STAGED_INFINITY_2026_09_05.md. No orbit solves."""
import json
import sympy as s

x, y, u, v, z, c, m, a = s.symbols('x y u v z c m a')
P = (1+x)*y+x*x
Q = -10*x*x+s.Rational(11,5)*x*y+c*y*y-m*x
sub = {x:u/v, y:1/v}
U = s.cancel(v*v*(P-u*Q).subs(sub))
V = s.cancel(-v**3*Q.subs(sub))
assert s.expand(U-((1-c)*u-s.Rational(6,5)*u*u+10*u**3+v*(1+m*u*u))) == 0
assert s.expand(V-(v*(10*u*u-s.Rational(11,5)*u-c)+m*u*v*v)) == 0
uv=-v+a*v*v
res=s.expand(U.subs(u,uv)-s.diff(uv,v)*V.subs(u,uv))
assert s.expand(res).coeff(v,2)==a*c+a+1
assert s.simplify((a*c+a+1).subs(a,-1/(1+c)))==0
assert s.expand(P.subs(x,-1))==1

subx={x:1/v,y:z/v}
Z=s.cancel(v*v*(Q-z*P).subs(subx))
VX=s.cancel(-v**3*P.subs(subx))
p=-10+s.Rational(6,5)*z+(c-1)*z*z
q=s.diff(p,z)
assert s.expand(Z-(p-v*(m+z*z)))==0
assert s.expand(VX-(-v*(1+z)-z*v*v))==0
assert s.discriminant(p,z)==40*c-s.Rational(964,25)
J=305+634*c-11*c*c-1000*c**3
neutral=(1-c)*(1+z)-c*q
assert s.factor(s.resultant(p,neutral,z)-(c-1)*J/25)==0
assert J.subs(c,s.Rational(241,250))>0
assert J.subs(c,s.Rational(39,40))<0
# J'<0 throughout [1/2,3/2]: J'(1/2)<0 and J''<0.
assert s.diff(J,c).subs(c,s.Rational(1,2))<0
assert s.expand(s.diff(J,c,2)+6000*c+22)==0
# Explicit positive smaller saddle root at neutrality.
zn=(11*c-5)/(5*(1+c-2*c*c))
assert s.factor(neutral.subs(z,zn))==0
assert s.factor(p.subs(z,zn)-J/(25*(c-1)*(2*c+1)**2)) == 0
# For c>1, mixed finite-saddle product=1 requires root sum=-2.
assert s.solve(s.Eq(-s.Rational(6,5)/(c-1),-2),c)==[s.Rational(8,5)]

root=s.polys.polytools.intervals(J,eps=s.Rational(1,10**14))
result={
    'status':'PASS',
    'orbit_evaluations':0,
    'checks':['vertical chart','finite-slope chart','vertical separatrix expansion',
              'one-way x=-1 crossing','mixed-stratum neutrality resultant',
              'unique real J root','c>1 neutrality outside inherited box'],
    'J_real_root_isolation':[[str(t[0][0]),str(t[0][1]),t[1]] for t in root],
    'J_real_root_approx':str(s.nroots(J)[0]),
    'limitations':'Local eigenvalue neutrality is not a global graphic connection or its transition coefficient.'
}
print(json.dumps(result,indent=2))
