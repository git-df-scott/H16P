"""Exact cubic Hopf coefficient in explicitly normalized real coordinates."""
from pathlib import Path
from fractions import Fraction as F
import json
import sympy as s
from reversible_finite_hopf import field
ROOT=Path(__file__).resolve().parent
u,v,A,B,omega,a,b,e2=s.symbols('u v A B omega a b e2',nonzero=True,real=True)
dy=(-omega*v-A*u)/B
F2=a*u*u+e2*u*dy+b*dy*dy
G2=-2*u*dy
f=F2/omega;g=-(A*F2+B*G2)/omega**2
fxx,fxy,fyy=[s.diff(f,*args) for args in [(u,u),(u,v),(v,v)]]
gxx,gxy,gyy=[s.diff(g,*args) for args in [(u,u),(u,v),(v,v)]]
coef=s.factor((fxy*(fxx+fyy)-gxy*(gxx+gyy)-fxx*gxx+fyy*gyy)/16)
# Multiplication by positive omega^5 B^4 preserves sign regardless of B.
poly=s.factor(coef*omega**5*B**4)
assert not poly.has(omega**-1)
rows=[]
for k in [F(1,10000),F(1,1000),F(1,100),F(1,20),F(1,10),F(1,5),F(3,10),F(2,5),F(1,2),F(3,4),F(1)]:
    try:d=field(k)[0]
    except (AssertionError,ArithmeticError):continue
    rat=lambda key:s.Rational(d[key])
    av,bv,x,y,ep2=map(rat,['a','b','x0','y0','epsilon2'])
    aa=2*x # trace zero implies A=2x
    bb=1-bv+2*bv*y+ep2*x
    det=rat('determinant')
    val=s.factor(poly.subs({A:aa,B:bb,a:av,b:bv,e2:ep2,omega:s.sqrt(det)}))
    assert val.is_Rational
    rows.append({'k':str(k),'positive_scaled_l1':str(val),'sign':int(s.sign(val)),'normalized_l1_approx':float(val/(det**s.Rational(5,2)*bb**4))})
out={'evidence':'Exact sign of quadratic Hopf cubic coefficient using standard normal-form formula, not a global cycle count','coordinates':'u=x-x0; v=-(A*u+B*(y-y0))/omega; normalized time omega*t; omega=sqrt(det)>0','coefficient':str(coef),'positive_scaled_polynomial':str(poly),'records':rows}
(ROOT/'reversible_hopf_focus.json').write_text(json.dumps(out,indent=2)+'\n')
print('coefficient',coef)
for r in rows:print(r['k'],r['sign'],r['normalized_l1_approx'])
