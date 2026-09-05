import json
from closure_budget import call,H
row=next(x for x in reversed(json.loads((H.parent/'fold_surface_2026_09_05/events_quad.json').read_text())) if x['status']=='ACCEPTED')
p=row['pair_profile'];checks=[]
for offset in [-.3,0,.3]:
 q=next(x['result'] for x in p['samples'] if x['log_offset']==offset)
 req={k:q[k] for k in ('r','c','K')};req['tol']='2e-26'
 checks.append(dict(offset=offset,half_claim=q,full=call(req,'outermost pair complete-return sign check',engine='angular_quad.py')))
(H/'outermost_check.json').write_text(json.dumps(dict(checks=checks,certified=False),indent=2))
