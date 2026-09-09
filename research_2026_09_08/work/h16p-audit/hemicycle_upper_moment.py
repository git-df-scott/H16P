"""NUM upper energy integral in the boundary-canceling direction."""
from pathlib import Path
import sys,json
import numpy as np
P=Path(__file__).resolve().parent
sys.path.insert(0,str(P/'replay/reversible_reseed'))
from moment_search import profile
out={'scope':'Finite energy/shape samples; weighted energy integral, not primitive zeros or a continuum proof','records':[]}
for a in [-.25,-.5,-.75,-1.25,-1.5,-1.75]:
 for b in [.25,.5,1.,1.5,1.75]:
  if a+b==0:continue
  k=np.sqrt(b*(a+2)/(a*(b-2)))
  lo=profile(a,b,1,n=61,order=256);hi=profile(a,b,1,n=61,order=512)
  d=1-2*k*hi[:,1];normalized=d/(1-k)
  item=dict(a=a,b=b,k=float(k),rows=[dict(energy=float(q[0]),weighted_mean_y=float(q[1]),normalized_energy_integral=float(v),relative_to_center_sign=float(n)) for q,v,n in zip(hi,d,normalized)],max_order_change=float(max(abs(2*k*(hi[:,1]-lo[:,1])))),min_relative_to_center_sign=float(min(normalized)),max_relative_to_center_sign=float(max(normalized)),sign_crossings=int(sum(d[:-1]*d[1:]<0)))
  out['records'].append(item);print(json.dumps({k:v for k,v in item.items() if k!='rows'}),flush=True)
  (P/'hemicycle_upper_moment.json').write_text(json.dumps(out,indent=2)+'\n')
