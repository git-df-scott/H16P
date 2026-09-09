from pathlib import Path
import json
from q4_finite_continuation import difference
P=Path(__file__).resolve().parent
rows=[]
for name in ['q4_finite_continuation.json','q4_finite_negative.json']:
    source=json.loads((P/name).read_text())
    for r in [.5,1.,2.]:
        eligible=[x for x in source['records'] if x['r']==r and len(x['brackets'])==3]
        chosen=eligible[-1];checks=[]
        for bracket in chosen['brackets']:
            for radius in bracket:
                result=difference(r,chosen['epsilon'],chosen['normalized_controls'],radius,rtol=3e-14)
                assert 'difference' in result
                old=next(x['difference'] for x in chosen['profile'] if x['radius']==radius)
                assert old*result['difference']>0
                checks.append({'radius':radius,'difference':result['difference'],'change':result['difference']-old,'return':result})
        rows.append({'r':r,'epsilon':chosen['epsilon'],'rational_parameters':chosen['rational_parameters'],'checks':checks})
(P/'q4_finite_continuation_refine.json').write_text(json.dumps({'scope':'NUM tighter bracket endpoint checks on six continuation endpoints','records':rows,'half_passages':2*sum(len(x['checks']) for x in rows)},indent=2)+'\n')
print('72 tighter half-passages passed; maximum change',max(abs(c['change']) for x in rows for c in x['checks']))
