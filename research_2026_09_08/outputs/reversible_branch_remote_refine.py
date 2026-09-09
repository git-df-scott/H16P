from pathlib import Path
from fractions import Fraction as F
import json
from reversible_finite_hopf import field,half
P=Path(__file__).resolve().parent
source=json.loads((P/'reversible_third_focus_branch.json').read_text());records=[]
for row in source['records']:
    if not row['lower_brackets']:continue
    p=row['selected']['field'];model=field(F(p['k']),F(p['m']),F(p['a']),F(p['b']));checks=[]
    for bracket in row['lower_brackets']:
        for radius in bracket:
            f=half(model,radius,-1,1,2e-13);b=half(model,radius,-1,-1,2e-13);assert f['status']==b['status']=='half_passage'
            delta=f['coordinate']-b['coordinate'];old=next(r['difference'] for r in row['lower'] if r['r']==radius);assert delta*old>0
            checks.append({'radius':radius,'difference':delta,'change':delta-old,'forward':f,'backward':b})
    records.append({'field':p,'checks':checks})
out={'scope':'NUM tighter checks of remote crossing endpoints; no certified cycles','records':records,'half_passages':2*sum(len(r['checks']) for r in records)}
(P/'reversible_branch_remote_refine.json').write_text(json.dumps(out,indent=2)+'\n');print('tight half passages',out['half_passages']);print('max change',max(abs(c['change']) for r in records for c in r['checks']))
print('max residual c3/c5',max(max(abs(row['selected']['local_order7']['return_coefficients_c2_to_order'][i]) for i in [1,3]) for row in source['records']))
