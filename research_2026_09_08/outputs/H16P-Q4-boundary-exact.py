"""Exact algebra for the original Q4 boundary parametrization.

The accompanying Markdown supplies the global graph/annulus argument.
"""
from pathlib import Path
import json
import sympy as S

z=S.symbols('z', real=True)
g=2*z**3-3*z-2
X=(3+8*z+4*z**4)/(12*(1-2*z*z))
Y=(1-8*z-24*z*z-16*z**3-4*z**4)/(12*(1-2*z*z))
P=-Y-3*X**2+2*X*Y+Y**2
Q=X+X**2-4*X*Y-Y**2
assert S.cancel(P-S.diff(X,z)*g/6)==0
assert S.cancel(Q-S.diff(Y,z)*g/6)==0
angular_cubic=2*z**3+6*z*z+3*z+1
assert S.cancel(X*Q-Y*P-g*g*angular_cubic/(54*(2*z*z-1)**2))==0
assert S.simplify(angular_cubic.subs(z,-1+1/S.sqrt(2))-(2-S.sqrt(2)))==0
D=1-4*Y+2*(X+Y)**2
N=6*Y*(1-X-Y)-1+2*(X+Y)**3
assert S.cancel(D-2*g*g/(9*(2*z*z-1)**2))==0
assert S.cancel(N+2*g**3/(27*(2*z*z-1)**3))==0
assert S.cancel(S.diff(X,z)+(2*z*z+1)*g/(3*(2*z*z-1)**2))==0
assert S.cancel(X+Y-1-2*(2*z**3+1)/(3*(2*z*z-1)))==0
assert X.subs(z,-S.Rational(2,5))<0<X.subs(z,-S.Rational(3,8))
y_numerator_lower=1+8*S.Rational(3,8)-24*S.Rational(2,5)**2+16*S.Rational(3,8)**3-4*S.Rational(2,5)**4
assert y_numerator_lower>0
rt=S.sqrt(2)
assert S.simplify(S.limit(Y/X,z,-1/rt)-(-1+rt))==0
assert S.simplify(S.limit(Y/X,z,1/rt)-(-1-rt))==0
a,q=S.symbols('a q', real=True)
H=(1-a**3+S.Rational(3,2)*a*(1+q*q))/q**3
critical_a=S.diff(H,a)*q**3
critical_q=S.cancel(-S.diff(H,q)*q**4/3)
assert S.expand(critical_a-S.Rational(3,2)*(1+q*q-2*a*a))==0
assert S.expand(critical_q.subs(q*q,2*a*a-1)-(1+a))==0
record={'status':'all exact assertions passed',
 'parameter_interval':'-1/sqrt(2) < z < 1/sqrt(2)',
 'X':str(X),'Y':str(Y),'physical_parameter_velocity':str(g/6),
 'graph':'dX/dz>0 on the interval; X ranges over all real numbers',
 'orientation':'forward flow goes from X=+infinity to X=-infinity',
 'boundary_level':'N/D^(3/2)=-1/sqrt(2)',
 'annulus_domain':'Y < boundary_graph(X), excluding the origin',
 'angular_chart':'Positive physical angular velocity throughout the open annulus, by the accompanying radial graph argument.',
 'origin_crossing_parameter_bracket':['-2/5','-3/8'],
 'positive_Y_numerator_lower_bound':str(y_numerator_lower),
 'scope':'Exact unperturbed boundary geometry only; no perturbed cycle or splitting certificate.',
 'ode_evaluations':0}
Path(__file__).with_suffix('.json').write_text(json.dumps(record,indent=2)+'\n')
print(json.dumps(record,indent=2))
