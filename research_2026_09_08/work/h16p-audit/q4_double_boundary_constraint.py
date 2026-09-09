"""Sampled compact zeros under vanishing splitting and exponent variation."""
from pathlib import Path
import json
import numpy as np
P=Path(__file__).resolve().parent
b=np.array([float(x['value']) for x in json.loads((P/'q4_boundary_integral.json').read_text())['entries']])
l=np.array([0.,5.,1.25,1.75])
A=np.array([b/np.linalg.norm(b),l/np.linalg.norm(l)])
_,sv,vh=np.linalg.svd(A,full_matrices=True);basis=vh[2:].T
rows=sorted(json.loads((P/'q4_original_probe_calls.json').read_text()),key=lambda r:r['fraction'])
# Deduplicate repeated saved calls at the same fraction.
rows=list({r['fraction']:r for r in rows}.values())
M=np.array([r['M'] for r in rows]);projected=M@basis
angles=sorted(set(float(np.arctan2(-q[0],q[1])%np.pi) for q in projected))
records=[]
for i,a in enumerate(angles):
    bangle=angles[(i+1)%len(angles)]+(np.pi if i==len(angles)-1 else 0)
    theta=(a+bangle)/2;d=basis@np.array([np.cos(theta),np.sin(theta)]);d/=np.max(np.abs(d))
    values=M@d
    brackets=[[rows[i]['fraction'],rows[i+1]['fraction']] for i in range(len(rows)-1) if values[i]*values[i+1]<0 and min(abs(values[i]),abs(values[i+1]))>1e-11]
    records.append({'direction':d.tolist(),'constraint_residuals':(A@d).tolist(),'brackets':brackets,'values':values.tolist()})
out={'scope':'NUM finite sample only; constraints are first-order splitting and exponent product, not exact perturbed connections or exponents','fractions':[r['fraction'] for r in rows],'constraint_matrix':A.tolist(),'nullspace':basis.tolist(),'singular_values':sv.tolist(),'max_sampled_crossings':max(map(lambda r:len(r['brackets']),records)),'records':records}
(P/'q4_double_boundary_constraint.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps({'tested_projective_sectors':len(records),'max_sampled_crossings':out['max_sampled_crossings'],'fractions':out['fractions']}))
