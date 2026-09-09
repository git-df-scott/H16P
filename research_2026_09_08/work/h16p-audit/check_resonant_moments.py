"""Independent exact and numerical check at the omitted a=-1,b=1 seed.

No ODE calls. Numerical quadrature is a sanity check, not an enclosure.
"""
import json
import sympy as s
import mpmath as mp
from pathlib import Path

x,y,h,k,R=s.symbols('x y h k R', positive=True)
P=-s.Rational(1,4)-x*x+y*y
Q=-2*x*y
H=(x*x+y*y+s.Rational(1,4))/y
assert s.factor(s.diff(H,x)*P+s.diff(H,y)*Q)==0
assert s.expand(y*(H-h)-(x*x+(y-h/2)**2-(h*h-1)/4))==0

# Integral of 1/y over the disk centered at (0,k), fixed radius R<k.
J1=2*s.pi*(k-s.sqrt(k*k-R*R))
J2=-s.diff(J1,k)
J3=-s.diff(J2,k)/2
sub={k:h/2,R:s.sqrt(h*h-1)/2}
answers=[s.simplify(j.subs(sub)) for j in [J1,J2,J3]]
expected=[s.pi*(h-1),2*s.pi*(h-1),2*s.pi*(h*h-1)]
assert all(s.simplify(a-b)==0 for a,b in zip(answers,expected))
assert s.simplify(answers[0]/answers[1]-s.Rational(1,2))==0
assert s.simplify(answers[2]/answers[1]-(h+1))==0

mp.mp.dps=60
checks=[]
for hs in ['1.001','1.5','2.5','10','100']:
    hv=mp.mpf(hs);kv=hv/2;rv=mp.sqrt(hv*hv-1)/2
    vals=[]
    for j in [1,2,3]:
        f=lambda t:2*rv*rv*mp.sin(t)**2/(kv+rv*mp.cos(t))**j
        vals.append(mp.quad(f,[0,mp.pi/2,mp.pi]))
    targets=[mp.pi*(hv-1),2*mp.pi*(hv-1),2*mp.pi*(hv*hv-1)]
    errs=[abs(a-b)/abs(b) for a,b in zip(vals,targets)]
    assert max(errs)<mp.mpf('1e-45'),(hs,errs)
    checks.append({'h':hs,'relative_errors':[str(e) for e in errs]})
result={'status':'exact identities plus nonvalidated independent quadrature',
        'moments':[str(a) for a in answers],
        'normalized_upper':'epsilon1 + epsilon2/2 - 2*epsilon0*(h+1)',
        'normalized_lower':'epsilon1 - epsilon2/2 + 2*epsilon0*(h+1)',
        'scope':'At most one zero per annulus for a nonzero first-order function at this fixed seed. Identically zero first order and endpoint bifurcations are NOT bounded.',
        'quadrature_checks':checks,
        'versions':{'sympy':s.__version__,'mpmath':mp.__version__}}
Path(__file__).with_suffix('.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result,indent=2))
