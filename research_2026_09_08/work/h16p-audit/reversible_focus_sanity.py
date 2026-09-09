"""Independent cubic return derivation by angular averaging."""
from pathlib import Path
import json
import sympy as s
ROOT=Path(__file__).resolve().parent
C,S=s.symbols('C S')
A,B,w,a,b,e=s.symbols('A B w a b e',nonzero=True,real=True)
dy=(-w*S-A*C)/B
f=(a*C*C+e*C*dy+b*dy*dy)/w
g=-(A*(a*C*C+e*C*dy+b*dy*dy)+B*(-2*C*dy))/w**2
radial=C*f+S*g;angular=C*g-S*f
poly=s.Poly(s.expand(-radial*angular),C,S)
average=0
for (i,j),coef in poly.terms():
    if i%2 or j%2:continue
    p,q=i//2,j//2
    average+=coef*s.factorial(2*p)*s.factorial(2*q)/(4**(p+q)*s.factorial(p)*s.factorial(q)*s.factorial(p+q))
average=s.factor(average)
data=json.loads((ROOT/'reversible_hopf_focus.json').read_text())
expected=s.sympify(data['coefficient'],locals={'A':A,'B':B,'omega':w,'a':a,'b':b,'e2':e})
assert s.factor(average-expected)==0
k=s.symbols('k',positive=True)
av=s.Rational(-7,4);bv=s.Rational(1,3);m=s.Rational(1539,1000)
x=-k/(2-av+bv*k*k);y=s.Rational(1,2)+k*x
e0=2*x*y;tau=(av-1)*e0;e2=-tau*m
AA=2*x;BB=1-bv+2*bv*y+e2*x;det=-AA**2+2*y*BB
# denominatorpositive scaling; substitution has only even powers of w.
scaled=s.factor((average*w**3*B**4).subs({A:AA,B:BB,a:av,b:bv,e:e2}))
scaled=s.factor(scaled.subs(w**4,det**2).subs(w**2,det))
num,den=s.fraction(scaled)
roots=s.polys.polytools.intervals(num,eps=s.Rational(1,10**12))
positive=[{'left':str(pair[0]),'right':str(pair[1]),'multiplicity':mult} for pair,mult in roots if pair[0]>0 and pair[0]<1]
out={'evidence':'Exact coefficient identity and rational root isolation; no higher focus coefficient or generalized Hopf theorem asserted','angular_average_equals_formula':True,'derivation':'dr/dtheta=r² R-r³ R T+O(r⁴); integral R=0. Cubic full return coefficient is -integral R*T=2*pi*l1. The iterated quadratic term has zero integral.','scaled_l1_numerator':str(num),'scaled_l1_denominator':str(den),'positive_roots_below_one':positive}
(ROOT/'reversible_focus_sanity.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out,indent=2))
