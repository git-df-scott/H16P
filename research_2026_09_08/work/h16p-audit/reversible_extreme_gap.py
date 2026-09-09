"""Resolve the coarse lower-annulus gap before its first failed return."""
import json,numpy as np
from reversible_log_return import P,paired,COUNTS
source=json.loads((P/'reversible_extreme_probe.json').read_text());out=dict(scope=__doc__,records=[])
for field in source['records']:
 p=field['parameters'];rows=next(q for q in field['profiles'] if q['side']==-1)['rows'];i=next(i for i in range(1,len(rows)) if 'log_difference' in rows[i-1] and 'log_difference' not in rows[i]);left=rows[i-1];hi=np.log(rows[i]['height']);trials=[];bracket=None
 for attempt in range(16):
  mid=(np.log(left['height'])+hi)/2;row=paired(p['a'],p['b'],p['e0'],p['e1'],p['e2'],-1,float(np.exp(mid)));trials.append(row)
  if 'log_difference' not in row:hi=mid
  elif row['log_difference']*left['log_difference']<0:bracket=[left['height'],row['height']];break
  else:left=row
 item=dict(parameters=p,trials=trials,bracket=bracket);out['records'].append(item);out['counts']=COUNTS.copy();print(json.dumps(dict(a=p['a'],attempts=len(trials),bracket=bracket)),flush=True)
 (P/'reversible_extreme_gap.json').write_text(json.dumps(out,indent=2)+'\n')
