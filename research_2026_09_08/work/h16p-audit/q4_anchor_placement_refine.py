"""Replay endpoint anchor-placement brackets at tighter numerical tolerance."""
import json
from q4_finite_continuation import P,difference,COUNTS
source=json.loads((P/'q4_anchor_placement.json').read_text());rows=[]
for a in source['records']:
    if abs(a['shift'])!=.2:continue
    for bracket in a['brackets']:
        pair=[]
        for R in bracket:
            old=next(q for q in a['profile'] if q['radius']==R)
            new=difference(a['r'],a['epsilon'],a['normalized_controls'],R,rtol=5e-14)
            pair.append(dict(old=old['difference'],new=new,change=abs(new['difference']-old['difference']) if 'difference' in new else None))
        rows.append(dict(r=a['r'],shift=a['shift'],bracket=bracket,replay=pair))
out=dict(scope=__doc__,records=rows,counts=COUNTS)
(P/'q4_anchor_placement_refine.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(dict(counts=COUNTS,max_change=max(q['change'] for a in rows for q in a['replay']),all_brackets_retained=all(a['replay'][0]['new']['difference']*a['replay'][1]['new']['difference']<0 for a in rows))))
