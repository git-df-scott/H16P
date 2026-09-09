"""Boundary-constrained profiles from frozen original-Q4 numerical integrals.

No new orbit evaluations. Two compact anchors plus the boundary functional
determine one projective control direction. Sparse profiles are not exclusions.
"""
from pathlib import Path
import json
import numpy as np

p=Path(__file__).resolve().parent
rows=json.loads((p/'q4_original_probe_calls.json').read_text())
b=np.array([float(x['value']) for x in json.loads((p/'q4_boundary_integral.json').read_text())['entries']])
vals={r['fraction']:np.array(r['M']) for r in rows}
profiles=[]
for anchors in [(0.25,.5),(.25,.75),(.5,.75),(.75,.9),(.9,.97)]:
    a=np.array([vals[t] for t in anchors]+[b])
    scale=np.linalg.norm(a,axis=0)
    _,sv,vh=np.linalg.svd(a/scale)
    d=vh[-1]/scale
    d/=max(abs(d))
    d*=1 if d[-1]>0 else -1
    profiles.append({'anchors':anchors,'direction':d.tolist(),'singular_values':sv.tolist(),
                     'boundary_residual':float(b@d),
                     'profile':[{'fraction':t,'value':float(v@d)} for t,v in sorted(vals.items())]})
out={'status':'numerical only; reuse of saved base-return integrals','profiles':profiles,'new_ode_returns':0,
     'scope':'Sparse signs do not exclude unobserved roots or establish a global compatibility theorem. Late anchors lack sampled brackets on every side.'}
(p/'q4_boundary_constrained.json').write_text(json.dumps(out,indent=2)+'\n')
print('Saved five constrained profiles; no additional sampled crossing located.')
