"""Last numerical controls; no new-cycle claim from a failed return."""
import json
from budget import call,HERE
rows=[]
e=next(q for q in reversed(json.loads((HERE/'events_logm.json').read_text())) if q['status']=='ACCEPTED')
for i,b in enumerate(e['pair_profile']['root_sign_brackets']):
 checks=[]
 for side in ['left','right','approximation']:
  q=b[side];req={k:q[k] for k in ('r','c','m')};req['tol']='2e-24'
  checks.append(dict(side=side,claim=q,result=call(req,'complete-return replay on large-m pair',engine='angular_m_quad.py')))
 rows.append(dict(cycle_index=i,checks=checks))
(HERE/'large_m_full_reproduction.json').write_text(json.dumps(dict(rows=rows),indent=2))
# CPU-limited tighter 1e17 connection control is retried at a looser tolerance.
inf=[]
for tol in ['2e-24','2e-25']:
 inf.append(call(dict(r='1e17',c='8/5',K='7.18499469640662040162379516212772519',tol=tol),'neutral infinity convergence tolerance control',engine='half_quad.py'))
(HERE/'infinity_tolerance_control.json').write_text(json.dumps(inf,indent=2))
