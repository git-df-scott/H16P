"""Endpoint displacement scale check, including threshold-suppressed crossings."""
import json
from q4_finite_continuation import P,difference,COUNTS
source=json.loads((P/'q4_anchor_limits.json').read_text());out={'scope':__doc__,'records':[]}
for a in source['records']:
    if abs(a['shift'])!=.2499:continue
    rows=[]
    for j,R in enumerate(a['anchors']):
        samples=[]
        for factor in [.8,1.,1.2]:
            normal=difference(a['r'],a['epsilon'],a['normalized_controls'],R*factor)
            tight=difference(a['r'],a['epsilon'],a['normalized_controls'],R*factor,rtol=5e-14)
            samples.append(dict(factor=factor,normal=normal,tight=tight,change=abs(normal['difference']-tight['difference']) if 'difference' in normal and 'difference' in tight else None))
        rows.append(dict(anchor=j,radius=R,samples=samples))
    rec=dict(r=a['r'],shift=a['shift'],rows=rows);out['records'].append(rec)
    print(json.dumps(dict(r=a['r'],shift=a['shift'],signs=[[q['tight'].get('difference') for q in x['samples']] for x in rows])),flush=True)
    out['counts']=COUNTS.copy();(P/'q4_anchor_limits_refine.json').write_text(json.dumps(out,indent=2)+'\n')
