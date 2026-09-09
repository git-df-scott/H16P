"""Independent exact first-integral check of unperturbed large-return endpoints."""
import json,numpy as np
from scipy.special import logsumexp
from reversible_log_return import P
source=json.loads((P/'reversible_log_return_v2_check.json').read_text());rows=[]
for a in source['large_controls']:
    aa=a['a'];b=a['b'];r=a['row'];side=r['side'];z0=np.log(r['height'])
    def logH(z):
        if aa==-1 and b==1:return float(np.logaddexp(z,np.log(.25)-z))
        coeff=np.array([b/(aa+2),-side*(b-1)/(aa+1),(b-2)/(4*aa)])
        value,sgn=logsumexp(np.log(abs(coeff))+np.array([2*z,z,0]),b=np.sign(coeff),return_sign=True)
        assert sgn>0;return float(aa*z+value)
    assert r['forward']['status']==r['backward']['status']=='passed'
    error=abs(logH(z0)-logH(r['forward']['log_height']))
    rows.append(dict(a=aa,b=b,side=side,start_height=r['height'],end_log_height=r['forward']['log_height'],log_energy_error=error))
out=dict(scope=__doc__,records=rows,max_log_energy_error=max(r['log_energy_error'] for r in rows))
(P/'reversible_log_energy_check.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
