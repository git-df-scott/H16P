"""NUM divided upper energy integral as a+b approaches zero, away from a=-1."""
from pathlib import Path
import sys,json,numpy as np
P=Path(__file__).resolve().parent;sys.path.insert(0,str(P/'replay/reversible_reseed'))
from moment_search import profile
out={'scope':'Finite difference in shape and nonvalidated quadrature; tests the divided coefficient, not isolated cycles of a center','records':[]}
selected=[-1.] if '--resonant-only' in sys.argv else [-.25,-.5,-.75,-1.25,-1.5,-1.75]
for a0 in selected:
 b=-a0;rows=[]
 for c in [-1e-3,1e-3,-1e-4,1e-4]:
  a=a0+c;k=np.sqrt(b*(a+2)/(a*(b-2)));q=profile(a,b,1,n=41,order=512)
  vals=(1-2*k*q[:,1])/c
  rows.append(dict(c=c,values=vals.tolist(),energies=q[:,0].tolist(),max_divided_coefficient=float(max(vals)),min_divided_coefficient=float(min(vals))))
 one=(np.array(rows[0]['values'])+rows[1]['values'])/2;two=(np.array(rows[2]['values'])+rows[3]['values'])/2
 item=dict(a0=a0,b=b,rows=rows,central_estimate_change=float(max(abs(one-two))),expected_center_limit=float(-1/(a0*(b-2))))
 out['records'].append(item);print(json.dumps({k:v for k,v in item.items() if k!='rows'}|{'max_values':max(z['max_divided_coefficient'] for z in rows)}),flush=True)
 (P/('hemicycle_center_locus_resonant.json' if '--resonant-only' in sys.argv else 'hemicycle_center_locus.json')).write_text(json.dumps(out,indent=2)+'\n')
