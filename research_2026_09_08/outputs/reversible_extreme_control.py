import json,numpy as np
from scipy.special import logsumexp
from reversible_log_return import P,paired,COUNTS
rows=[]
for a,b in [(-1.,1.),(-1.25,.5)]:
 for side in [1,-1]:
  T=1e250;r=paired(a,b,0,0,0,side,T)
  assert 'log_difference' in r
  def logH(z):
   if a==-1:return float(np.logaddexp(z,np.log(.25)-z))
   co=np.array([b/(a+2),-side*(b-1)/(a+1),(b-2)/(4*a)])
   val,sg=logsumexp(np.log(abs(co))+np.array([2*z,z,0]),b=np.sign(co),return_sign=True);assert sg>0
   return float(a*z+val)
  rows.append(dict(a=a,b=b,side=side,row=r,log_energy_error=abs(logH(np.log(T))-logH(r['forward']['log_height']))))
(P/'reversible_extreme_control.json').write_text(json.dumps(dict(scope='Exact first-integral endpoint test at initial height1e250; nonvalidated error',records=rows,counts=COUNTS),indent=2)+'\n');print(json.dumps(dict(counts=COUNTS,max_log_energy_error=max(r['log_energy_error'] for r in rows))))
