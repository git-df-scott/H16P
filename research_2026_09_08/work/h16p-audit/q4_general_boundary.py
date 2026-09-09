"""Exact original Q4 identities for arbitrary positive center modulus r.

The parametrization uses -r/sqrt(1+r²)<s<r/sqrt(1+r²).
No perturbed cycle or zero-count claim is made.
"""
from pathlib import Path
import json
import sympy as S

X,Y,r,s=S.symbols('X Y r s',nonzero=True)
d=1+r*r
P=-Y-(2+r*r)*X*X+2*r*X*Y+Y*Y
Q=X+r*X*X-(1+3*r*r)*X*Y-r*Y*Y
D=1-2*d*Y+d*(r*X+Y)**2
N=3*d*Y*(1-r*X-Y)-1+d*(r*X+Y)**3
L=lambda f:S.expand(S.diff(f,X)*P+S.diff(f,Y)*Q)
assert S.expand(L(D)+2*d*X*D)==0
assert S.expand(L(N)+3*d*X*N)==0
assert S.expand(S.diff(N,X)*D-S.Rational(3,2)*N*S.diff(D,X)-3*r*r*d*Q)==0
assert S.expand(S.diff(N,Y)*D-S.Rational(3,2)*N*S.diff(D,Y)+3*r*r*d*P)==0
g=d*s**3-3*r*r*s-2*r*r
k=d*s*s-r*r
a=2*(r*r+d*s**3)/(3*k)
q_squared=d*g*g/(9*k*k)
xs=S.factor((q_squared-d*a*a+r*r)/(2*d*r))
ys=S.factor(a+1-r*xs)
on={X:xs,Y:ys}
ps=S.cancel(P.subs(on));qs=S.cancel(Q.subs(on))
assert S.cancel(ps-S.diff(xs,s)*g/(6*r))==0
assert S.cancel(qs-S.diff(ys,s)*g/(6*r))==0
assert S.cancel(D.subs(on)-q_squared)==0
assert S.cancel(N.subs(on)+d*g**3/(27*k**3))==0
assert S.cancel(S.diff(xs,s)+(d*s*s+r*r)*g/(3*r*k*k))==0
cubic=d*s**3+3*d*s*s+3*r*r*s+r*r
assert S.cancel(xs*qs-ys*ps-g*g*cubic/(54*r*r*k*k))==0
assert S.simplify(cubic.subs(s,-1+1/S.sqrt(d))-(2-2/S.sqrt(d)))==0
out={'status':'all exact assertions passed','modulus':'r>0','base':{'P':str(P),'Q':str(Q)},
     'D':str(D),'N':str(N),'gradient_factor':'3*r^2*(1+r^2)*D^(-5/2)',
     'parameter_interval':'-r/sqrt(1+r²)<s<r/sqrt(1+r²)',
     'X':str(xs),'Y':str(ys),'parameter_velocity':str(g/(6*r)),
     'boundary_level':'h=N/D^(3/2)=-1/sqrt(1+r²)',
     'angular_cubic_minimum':'2-2/sqrt(1+r²)>0',
     'scope':'Exact unperturbed base and connection identities. No perturbation compatibility, remainder, or cycle certificate.'}
Path(__file__).with_suffix('.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out,indent=2))
