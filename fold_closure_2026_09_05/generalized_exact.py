"""Exact identities behind the one-coefficient cusp test. Zero ODE calls."""
import sympy as s,json
from pathlib import Path
B,c,m=s.symbols('B c m',positive=True)
# Standard normalized cubic focus numerator from quadratic derivatives.
fxx=2/s.sqrt(m);fxy=-1;fyy=0
gxx=20/m;gxy=B/s.sqrt(m);gyy=-2*c
lyap=s.factor((fxy*(fxx+fyy)-gxy*(gxx+gyy)-fxx*gxx+fyy*gyy)/16)
K=m*(B*c-1)-10*(B+2)
assert s.simplify(lyap-K/(8*m**s.Rational(3,2)))==0
# Reversible-center parity coefficient along K=0, when Bc !=1.
m0=10*(B+2)/(B*c-1);sig=(1+2*c)/(B+2);D=1+m0*sig**2
coef=(-10*sig**2-B*sig+c-m0*sig*(sig**2-sig))/D**2
H=B**2*c**2+B**2*c-2*B*c**2-B*c-B+40*c**3-28*c-10
expected=-(B+2)*(B*c-1)*H/(B**2*c+2*B*c-B+40*c*c+40*c+8)**2
assert s.factor(coef-expected)==0
assert s.factor(H.subs(c,1)-2*(B-1)**2)==0
# A possible projected endpoint must not be accepted as a finite field.
assert K.subs({B:1,c:1})==-30
r,h1,h2,g1,g2=s.symbols('r h1 h2 g1 g2')
velocity=-r*(h1+r*h2)/(g1+r*g2)
third=6*g1*g2*(h2*g1-h1*g2)/(g1+r*g2)**4
assert s.factor(s.diff(velocity,r,3)-third)==0
out=dict(status='EXACT_IDENTITIES_ONLY',focus_coefficient=str(lyap),center_parity_polynomial=str(H),center_projection_at_c1=str(s.factor(H.subs(c,1))),K_at_B1_c1='-30',angular_velocity_third_derivative=str(third),cycle_count_bound=False)
Path(__file__).with_suffix('.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out,indent=2))
