"""Tighter replay of nearest-boundary paired passages."""
import json
from q4_finite_continuation import P,difference,COUNTS
fields={(a['r'],a['epsilon']):a for a in json.loads((P/'q4_finite_extended.json').read_text())['records']}
rows=[]
for a in json.loads((P/'q4_saddle_boundary_profile.json').read_text())['records']:
    old=a['rows'][-1];f=fields[(a['r'],a['epsilon'])]
    new=difference(a['r'],a['epsilon'],f['normalized_controls'],old['radius'],rtol=5e-14)
    rows.append(dict(r=a['r'],epsilon=a['epsilon'],old_difference=old.get('difference'),refined=new,change=abs(old['difference']-new['difference']) if 'difference' in new else None))
out=dict(scope='Tighter nonvalidated nearest-boundary replay',records=rows,counts=COUNTS)
(P/'q4_saddle_boundary_refine.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(dict(counts=COUNTS,all_resolved=all(a['change'] is not None for a in rows),max_change=max(a['change'] for a in rows if a['change'] is not None),all_negative=all(a['refined'].get('difference',1)<0 for a in rows))))
