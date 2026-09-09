import json
from reversible_log_return import P,paired,COUNTS
source=json.loads((P/'reversible_log_boundary_probe.json').read_text());p=source['parameters'];rows=[]
for record in source['records']:
 for bracket in record['brackets']:
  pair=[]
  for T in bracket:
   old=next(x for x in record['rows'] if x['height']==T)
   new=paired(p['a'],p['b'],p['e0'],p['e1'],p['e2'],record['side'],T,rtol=3e-13)
   pair.append(dict(old=old['log_difference'],new=new,change=abs(new['log_difference']-old['log_difference']) if 'log_difference' in new else None))
  rows.append(dict(side=record['side'],bracket=bracket,replay=pair))
out=dict(scope='Tighter nonvalidated replay of all three boundary brackets',records=rows,counts=COUNTS)
(P/'reversible_log_boundary_refine.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(dict(counts=COUNTS,max_change=max(x['change'] for r in rows for x in r['replay']),all_preserved=all(r['replay'][0]['new']['log_difference']*r['replay'][1]['new']['log_difference']<0 for r in rows))))
