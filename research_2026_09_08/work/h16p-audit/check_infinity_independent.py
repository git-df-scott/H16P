"""Independent exact checks of the inherited KKL infinity formulas.

These check local algebra, not global connections or periodic orbits.
"""
from pathlib import Path
import json
import sympy as s

x,y,u,v,z,c,m,A=s.symbols('x y u v z c m A')
P=(1+x)*y+x*x
Q=-10*x*x+s.Rational(11,5)*x*y+c*y*y-m*x
vertical={x:u/v,y:1/v}
du=s.cancel(v*v*(P-u*Q).subs(vertical))
dv=s.cancel(-v**3*Q.subs(vertical))
du_expected=(1-c)*u-s.Rational(6,5)*u*u+10*u**3+v*(1+m*u*u)
dv_expected=v*(10*u*u-s.Rational(11,5)*u-c)+m*u*v*v
assert s.expand(du-du_expected)==0
assert s.expand(dv-dv_expected)==0
graph=-v+A*v*v
invariance=s.expand((du-s.diff(graph,v)*dv).subs(u,graph))
assert s.expand(s.expand(invariance).coeff(v,2)-(A*(c+1)+1))==0
assert P.subs(x,-1)==1
finite={x:1/v,y:z/v}
dz=s.cancel(v*v*(Q-z*P).subs(finite))
dv_finite=s.cancel(-v**3*P.subs(finite))
p=-10+s.Rational(6,5)*z+(c-1)*z*z
assert s.expand(dz-(p-v*(m+z*z)))==0
assert s.expand(dv_finite-(-v*(1+z)-z*v*v))==0
J=305+634*c-11*c*c-1000*c**3
neutral=(1-c)*(1+z)-c*s.diff(p,z)
resultant=s.factor(s.resultant(p,neutral,z))
assert s.expand(resultant-(c-1)*J/25)==0
# J' is already negative at 1/2 and J''<0 on the whole interval.
assert s.diff(J,c).subs(c,s.Rational(1,2))<0
assert s.expand(s.diff(J,c,2)-(-22-6000*c))==0
assert J.subs(c,s.Rational(241,250))>0
assert J.subs(c,s.Rational(39,40))<0
out={'status':'exact assertions passed',
     'scope':'Local coordinate identities, invariant-manifold coefficient, one-way barrier, and unique candidate neutrality root only. No global connection or cycle existence certified.',
     'neutrality_resultant':str(resultant),
     'root_bracket':['241/250','39/40'],
     'ode_evaluations':0}
Path(__file__).with_suffix('.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out,indent=2))
