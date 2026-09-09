import json
from reversible_log_return import P,paired,COUNTS
source=json.loads((P/'reversible_extreme_probe.json').read_text());gaps=json.loads((P/'reversible_extreme_gap.json').read_text());out=dict(scope='Tighter replay of all extreme-field brackets',records=[])
for field,gap in zip(source['records'],gaps['records']):
 p=field['parameters'];brackets=[(q['side'],b) for q in field['profiles'] for b in q['brackets']]+([(-1,gap['bracket'])] if gap['bracket'] else [])
 for side,bracket in brackets:
  rows=[paired(p['a'],p['b'],p['e0'],p['e1'],p['e2'],side,T,rtol=3e-13) for T in bracket]
  item=dict(a=p['a'],side=side,bracket=bracket,rows=rows,signs_retained=bool(all('log_difference' in q for q in rows) and rows[0]['log_difference']*rows[1]['log_difference']<0));out['records'].append(item)
out['counts']=COUNTS.copy();(P/'reversible_extreme_refine.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(dict(counts=COUNTS,all_retained=all(r['signs_retained'] for r in out['records']))))
