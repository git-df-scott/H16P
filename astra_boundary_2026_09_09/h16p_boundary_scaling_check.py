"""Exact infinity-chart calculation and bounded replay of inherited moments.

Requires sympy, numpy, scipy. Pass the main-checkout path as argv[1].
No cycle-existence claim. Does not modify inherited data.
"""
import json
import sys
from pathlib import Path
import sympy as S
import numpy as np

x,y,u,v,alpha,b,e0,e1,gamma=S.symbols('x y u v alpha b e0 e1 gamma')
P=(b-2)/4+(1-b)*y+(-2+alpha)*x*x+b*y*y+e1*x
Q=e0-2*x*y+gamma*x*x
sub={x:1/u,y:v/u}
# d/ds = u d/dt, a signed desingularization. Time orientation must be
# accounted for separately in the two hemispheres when constructing returns.
up=S.expand((-u**3*P).subs(sub))
vp=S.expand((u**2*(Q-v*P)).subs(sub))
up_expected=-u*(-2+alpha+b*v*v+e1*u+(1-b)*u*v+(b-2)*u*u/4)
vp_expected=gamma-alpha*v-b*v**3-u*(e1*v+(1-b)*v*v)+u*u*(e0-(b-2)*v/4)
assert S.simplify(up-up_expected)==0
assert S.simplify(vp-vp_expected)==0
angular=S.expand(vp.subs(u,0))
assert angular==gamma-alpha*v-b*v**3
d,A,G,z=S.symbols('delta A G z', nonzero=True)
scaled=S.expand(angular.subs({alpha:A*d*d,gamma:G*d**3,v:d*z})/d**3)
assert scaled==G-A*z-b*z**3
disc=S.factor(S.discriminant(angular,v))
assert disc==-b*(4*alpha**3+27*b*gamma**2)
repo=Path(sys.argv[1]).resolve()
sys.path.insert(0,str(repo/'reversible_reseed'))
from boundary_search import boundary_profile
source=repo/'reversible_reseed/data/boundary_search.json'
data=json.loads(source.read_text())
record=next(r for r in data['records'] if r['a']==-2 and abs(r['b']-1/3)<1e-15)
m,c=record['best']['m'],record['best']['c']
replays=[]
for order in (200,400,800):
    for side in (1,-1):
        profile=boundary_profile(-2.,1/3,side,order=order)
        values=side*(profile[:,2]-m*profile[:,1])-c
        crossings=np.where(values[:-1]*values[1:]<0)[0]
        witnesses=[]
        for i in crossings:
            witnesses.append({'indices':[int(i),int(i+1)],
                'energies':profile[i:i+2,0].tolist(),
                'values':values[i:i+2].tolist()})
        replays.append({'order':order,'side':side,'brackets':witnesses})
print(json.dumps({
    'scope':'Exact chart identities; floating-point moment replay only. No finite-field cycles certified.',
    'main_commit':'498b58fc1378486601f348f67e0ab3a5362e7933',
    'pr6_head_checked':'31d8a41a5c7f93ad9fb2d5955f269557abf2d419',
    'sympy_version':S.__version__,
    'infinity_u_equation':str(up_expected),
    'infinity_v_equation':str(vp_expected),
    'angular_cubic':str(angular),'discriminant':str(disc),
    'scaled_angular_cubic':str(scaled),
    'inherited_m':m,'inherited_c':c,
    'proposed_path':{'b':'1/3','a':'-2+A*delta^2',
        'epsilon0':'-delta^3/3','epsilon1':'-c*delta^3','gamma':'m*delta^3/3'},
    'replays':replays,
    'unproved':['Finite-field persistence at an explicit delta',
        'Existence of any additional boundary cycle',
        'Uniform return maps and remainders',
        'Isolation and geometric distinctness'],
},indent=2))
