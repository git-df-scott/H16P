"""Exact first variation of the original two-saddle exponent product.

Positive-X negative-slope saddle is passed first, then the negative-X
positive-slope saddle. No global transition coefficient is computed.
"""
from pathlib import Path
import json
import sympy as S

z,u,v,w=S.symbols('z u v w')
p=-3+(2+2*u+w)*z+z*z
q=1+u+(-4+v)*z-(1+u)*z*z
R=S.expand(q-z*p)
base={u:0,v:0,w:0}
roots=[-1-S.sqrt(2),-1+S.sqrt(2)]
for root in roots:
    assert S.simplify(R.subs(base).subs(z,root))==0
    assert S.simplify(p.subs(base).subs(z,root)+2)==0
    assert S.simplify(S.diff(R,z).subs(base).subs(z,root)+4)==0
values={'tau':S.Integer(0)}
for control in [u,v,w]:
    terms=[]
    for root in roots:
        dz=S.diff(R,control).subs(base).subs(z,root)/4
        angular=(S.diff(R,z,control).subs(base).subs(z,root)+S.diff(R,z,2).subs(base).subs(z,root)*dz)/(-4)
        radial=(S.diff(p,control).subs(base).subs(z,root)+S.diff(p,z).subs(base).subs(z,root)*dz)/(-2)
        terms.append(S.simplify(angular-radial))
    values[str(control)]=S.simplify(terms[0]-terms[1])
assert values=={'tau':0,'u':5*S.sqrt(2),'v':5*S.sqrt(2)/4,'w':7*S.sqrt(2)/4}
out={'status':'exact identities passed','base_exponents':['2','1/2'],'base_product':'1',
     'log_product_derivatives':{k:str(v) for k,v in values.items()},
     'normalization':'rho=Rprime(zminus)*p(zplus)/(p(zminus)*Rprime(zplus)), physical saddle itinerary as in boundary proof',
     'scope':'Local saddle data only. Global splitting, regular transition, uniform Dulac remainders, and cycle coexistence are separate obligations.'}
Path(__file__).with_suffix('.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out,indent=2))
