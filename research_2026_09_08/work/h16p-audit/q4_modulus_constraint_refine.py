from pathlib import Path
import json
import numpy as np
from q4_modulus_boundary_constraints import scan
P=Path(__file__).resolve().parent
out={'scope':'Nonvalidated refinement prompted by large absolute order changes at high modulus','records':[]}
for r in [2.,4.,8.]:
    runs=[scan(r,n) for n in [1024,2048]]
    changes=[]
    for a,b in zip(runs[0]['rows'],runs[1]['rows']):
        if 'M' in a and 'M' in b:
            x=np.array(a['M']);y=np.array(b['M']);changes.append({'fraction':a['fraction'],'absolute':float(max(abs(x-y))),'scaled':float(max(abs(x-y)/(1+abs(y))))})
    record={'modulus':r,'runs':runs,'changes':changes};out['records'].append(record)
    (P/'q4_modulus_constraint_refine.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({'modulus':r,'max_crossings':[x['max_sampled_crossings'] for x in runs],'max_abs':max(x['absolute'] for x in changes),'max_scaled':max(x['scaled'] for x in changes)}),flush=True)
