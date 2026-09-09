"""Local compactified infinity eigenvalues for outer-anchor endpoint fields.
Chart x=1/z,y=s/z, time rescaled by z: z'=-p2(s)z+O(z^2),
s'=q2(s)-s*p2(s)+O(z). No global connection inferred.
"""
import json
import numpy as np
from numpy.polynomial import Polynomial as Poly
from q4_finite_continuation import P
source=json.loads((P/'q4_anchor_limits.json').read_text());out={'scope':__doc__,'records':[]}
for a in source['records']:
    if a['shift']<0:continue
    r=a['r'];tau,u,w=a['epsilon']*np.array(a['normalized_controls']);v=-a['epsilon']
    p=Poly([-(2+r*r),2*r+2*u+w,1]);q=Poly([r+u,-1-3*r*r+v,-r-u]);R=q-Poly([0,1])*p
    points=[]
    for z in R.roots():
        if abs(z.imag)>1e-10:continue
        s=float(z.real);pr=float(p(s));angular=float(R.deriv()(s))
        points.append(dict(slope=s,radial_eigenvalue=-pr,angular_eigenvalue=angular,saddle=bool(pr*angular>0),residual=float(abs(R(s)))))
    saddles=[z for z in points if z['saddle']]
    product=None
    if len(saddles)==2:
        low,high=sorted(saddles,key=lambda z:z['slope'])
        product=low['angular_eigenvalue']*high['radial_eigenvalue']/(low['radial_eigenvalue']*high['angular_eigenvalue'])
    out['records'].append(dict(r=r,shift=a['shift'],points=points,conditional_two_saddle_product=product))
(P/'q4_anchor_infinity.json').write_text(json.dumps(out,indent=2)+'\n')
for a in out['records']:
 if a['shift']==.2499:print(json.dumps(a))
