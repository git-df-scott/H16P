"""Two fixed near-resonant fields at extreme but representable perturbation scales."""
import json,numpy as np
from reversible_log_return import P,paired,COUNTS
out={'scope':'NUM same-coefficient paired-annulus searches; unresolved returns preserved','records':[]}
for alpha,small,top in [(.01,1e-100,1e110),(.003,1e-250,1e275)]:
 a=-1-alpha;b=1.;k=np.sqrt((a+2)/(-a));eminus=.01;eplus=eminus*alpha*.001
 e0=-small;e1=(eminus-eplus)/(4*k);e2=-(eminus+eplus)/2
 item=dict(parameters=dict(a=a,b=b,e0=e0,e1=float(e1),e2=e2,epsilon_minus=eminus,epsilon_plus=eplus),profiles=[])
 for side in [1,-1]:
  rows=[paired(a,b,e0,e1,e2,side,float(T)) for T in np.geomspace(.6,top,41)]
  brackets=[[x['height'],y['height']] for x,y in zip(rows,rows[1:]) if 'log_difference' in x and 'log_difference' in y and x['log_difference']*y['log_difference']<0 and min(abs(x['log_difference']),abs(y['log_difference']))>1e-8]
  profile=dict(side=side,rows=rows,brackets=brackets,unresolved=sum('log_difference' not in r for r in rows));item['profiles'].append(profile)
  print(json.dumps(dict(alpha=alpha,side=side,brackets=brackets,unresolved=profile['unresolved'])),flush=True)
 out['records'].append(item);out['counts']=COUNTS.copy();(P/'reversible_extreme_probe.json').write_text(json.dumps(out,indent=2)+'\n')
