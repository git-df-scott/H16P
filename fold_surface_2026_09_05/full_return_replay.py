"""Independent complete angular returns at pair brackets selected from half shooting.
No discovery optimization. This verifies two-cycle evidence only, not >=5.
"""
import json
from budget import call,HERE
rows=[]
for name in ['events_quad.json','events_negative.json']:
 p=json.loads((HERE/name).read_text());e=next(x for x in reversed(p) if x['status']=='ACCEPTED')
 profile=e['pair_profile'];claims=profile['root_sign_brackets'];replays=[]
 for i,b in enumerate(claims):
  checks=[]
  for side in ['left','right','approximation']:
   q=b[side];req={k:q[k] for k in ('r','c','K')};req['tol']='2e-26' if name=='events_quad.json' else '2e-24'
   checks.append(dict(side=side,half_claim=q,full=call(req,'complete full-return replay of two-sided pair bracket',engine='angular_quad.py')))
  replays.append(dict(cycle_index=i,checks=checks))
 rows.append(dict(source=name,field={k:profile[k] for k in ('c','K')},replays=replays))
(HERE/'full_return_reproduction.json').write_text(json.dumps(dict(purpose='two-cycle numerical cross-formulation check; not five-cycle trigger',rows=rows),indent=2))
