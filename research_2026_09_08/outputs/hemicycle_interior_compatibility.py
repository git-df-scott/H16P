"""Exact energy derivative for the leading upper-double/lower-single direction."""
import sympy as S,json
from pathlib import Path
x,a,b,k=S.symbols('x a b k',real=True)
t=S.symbols('t',positive=True)
l=b/(a+2);m=-(b-1)/(a+1);n=(b-2)/(4*a)
checks=[]
for sign in [1,-1]:
    H=t**a*(x*x+l*t*t+sign*m*t+n)
    P=(b-2)/4+(1-b)*sign*t+a*x*x+b*t*t
    tdot=-2*x*t
    assert S.simplify(S.diff(H,x)*P+S.diff(H,t)*tdot)==0
    pert=S.factor(S.diff(H,x)*x*(1-2*k*sign*t))
    assert S.simplify(pert-2*t**a*x*x*(1-2*k*sign*t))==0
    checks.append(dict(half_plane=sign,unperturbed_Lie_derivative='0',perturbation_Lie_derivative=str(pert)))
out=dict(status='exact identities passed',assumptions='a in (-2,0) excluding -1; b in (0,2); t=abs(y)>0; k=sqrt(b(a+2)/(a(b-2)))>0',direction='epsilon0=0, epsilon1=delta, epsilon2=-2*k*delta',checks=checks,conclusion='For delta>0 the lower-half-plane energy derivative is positive except x=0. There is no lower periodic orbit entirely in y<0 for this exact direction. Along perturbation arcs with this nonzero leading direction, the lower first Melnikov integral is strictly positive on each regular compact orbit. This does not exclude boundary cycles, upper interior zeros, or other limiting directions.')
Path(__file__).with_suffix('.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
