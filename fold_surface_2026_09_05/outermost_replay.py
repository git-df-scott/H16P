"""Full-return signs at coarse brackets of the last positive-sheet pair."""
import json
from budget import call,HERE
row=next(x for x in reversed(json.loads((HERE/'events_quad.json').read_text())) if x['status']=='ACCEPTED')
p=row['pair_profile'];checks=[]
for offset in [-.3,0,.3]:
 claim=next(x['result'] for x in p['samples'] if x['log_offset']==offset)
 req={k:claim[k] for k in ['r','c','K']};req['tol']='2e-26'
 checks.append(dict(log_offset=offset,half_claim=claim,full_return=call(req,'complete-return coarse signs at outermost positive fold pair',engine='angular_quad.py')))
(HERE/'outermost_full_reproduction.json').write_text(json.dumps(dict(field=dict(c=p['c'],K=p['K']),checks=checks,certified=False),indent=2))
