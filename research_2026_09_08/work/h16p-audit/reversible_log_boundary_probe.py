"""Large-radius return probe of one near-resonant coupled boundary perturbation."""
import numpy as np,json
from reversible_log_return import P,paired,COUNTS
A=-1.01;B=1.;K=np.sqrt(B*(A+2)/(A*(B-2)))
eminus=.01;eplus=1e-6;e0=-1e-20;e1=(eminus-eplus)/(4*K);e2=-(eplus+eminus)/2
out=dict(scope='NUM same quadratic field on both annuli; no validated roots or exclusion between samples',parameters=dict(a=A,b=B,e0=e0,e1=float(e1),e2=e2,epsilon_minus=eminus,epsilon_plus=eplus),records=[])
for side in [1,-1]:
 rows=[]
 for T in np.geomspace(.6,1e40,61):
  row=paired(A,B,e0,e1,e2,side,float(T));rows.append(row)
 brackets=[[x['height'],y['height']] for x,y in zip(rows,rows[1:]) if 'log_difference' in x and 'log_difference' in y and x['log_difference']*y['log_difference']<0 and min(abs(x['log_difference']),abs(y['log_difference']))>1e-8]
 item=dict(side=side,rows=rows,brackets=brackets,unresolved=sum('log_difference' not in x for x in rows));out['records'].append(item);out['counts']=COUNTS.copy()
 print(json.dumps(dict(side=side,brackets=brackets,unresolved=item['unresolved'])),flush=True)
 (P/'reversible_log_boundary_probe.json').write_text(json.dumps(out,indent=2)+'\n')
