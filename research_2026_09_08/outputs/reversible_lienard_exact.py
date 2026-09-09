"""Exact scalar reduction on each invariant side of y=0 for epsilon0 !=0."""
from pathlib import Path
import sympy as s,json
ROOT=Path(__file__).resolve().parent
y,v,a,b,e0,e1,e2,y0,m=s.symbols('y v a b e0 e1 e2 y0 m',nonzero=True,real=True)
x=(e0-v)/(2*y)
g=(b-2)/4+(1-b)*y+b*y*y
P=g+a*x*x+e1*x+e2*x*y
second=s.expand(-2*P*y-2*x*v)
D=e1+e2*y+(a-1)*e0/y
C=a*e0**2/(2*y)+e1*e0+e2*e0*y+2*y*g
assert s.factor(second-((1-a/2)*v*v/y+D*v-C))==0
hopf=s.factor(D.subs({e1:-e2*y0-(a-1)*e0/y0,e2:-m*(a-1)*e0},simultaneous=True))
# sequential substitution for the relation between epsilon2 and epsilon0
hopf=s.factor(hopf.subs(e2,-m*(a-1)*e0))
expected=-(a-1)*e0*(y-y0)*(m*y*y0+1)/(y*y0)
assert s.factor(hopf-expected)==0
out={'evidence':'Exact identities, not a cycle-count theorem','scalar_equation':'y_ddot=(1-a/2)*y_dot²/y+D(y)*y_dot-C(y)','D':str(D),'C':str(s.factor(C)),'Hopf_D':str(hopf),'positive_side_transformation':'z=(2/a)*y^(a/2), y>0,a!=0; z_ddot=D(y)*z_dot-y^(a/2-1)*C(y)','barrier':'For epsilon0 !=0, dy/dt=epsilon0 on y=0, so no periodic orbit can cross y=0.','scope':'Damping zero at selected positive Hopf focus and negative root -1/(m*y0) when m,y0>0. No uniqueness theorem follows from this factorization alone.'}
(ROOT/'reversible_lienard_exact.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
