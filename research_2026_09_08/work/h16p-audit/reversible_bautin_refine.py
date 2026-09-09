from pathlib import Path
from fractions import Fraction as F
import json
from reversible_finite_hopf import field,half
P=Path(__file__).resolve().parent
source=json.loads((P/'reversible_bautin_precursor.json').read_text());records=[]
for row in source['records']:
    if not row['upper_brackets']:continue
    p=row['field'];model=field(F(p['k']),F(p['m']),F(p['a']),F(p['b']));checks=[]
    for bracket in row['upper_brackets']:
        signs=[]
        for r in bracket:
            f=half(model,r,1,1,2e-13);b=half(model,r,1,-1,2e-13)
            assert f['status']==b['status']=='half_passage'
            delta=f['coordinate']-b['coordinate'];old=next(q['difference'] for q in row['upper'] if q['r']==r)
            assert delta*old>0;signs.append(delta)
            checks.append({'r':r,'difference':delta,'change':delta-old,'forward':f,'backward':b})
        assert signs[0]*signs[1]<0
    records.append({'field':p,'checks':checks})
out={'scope':'NUM tighter bracket endpoint checks only; no rigorous orbit existence or isolation','records':records,'tested_fields':len(source['records']),'main_half_passages':2*sum(len(r['upper'])+len(r['lower']) for r in source['records']),'tight_half_passages':2*sum(len(r['checks']) for r in records),'field_gate_failures':len(source['gates'])}
(P/'reversible_bautin_refine.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps({k:v for k,v in out.items() if k!='records'}));print('max endpoint change',max(abs(c['change']) for r in records for c in r['checks']))
