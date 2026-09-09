from pathlib import Path
from fractions import Fraction as F
import json,time
import numpy as np
from reversible_finite_hopf import field,profile,count,half
ROOT=Path(__file__).resolve().parent
out={'evidence':'NUM only; continuation of finite rational Hopf family','records':[],'tightened_controls':[]}
begin=time.perf_counter()
for k in [F(3,10),F(2,5),F(1,2),F(3,4),F(1)]:
    try:
        model=field(k)
        upper=profile(model,np.geomspace(.01,1e7,49),1,2e-11)
        lower=profile(model,np.geomspace(.1,1e7,41),-1,2e-11)
        row={'field':model[0],'equilibria_numerical':model[1],'upper':upper,'lower':lower,'upper_brackets':count(upper),'lower_brackets':count(lower),'unresolved':sum('difference' not in r for r in upper+lower)}
    except (ArithmeticError,AssertionError) as exc:row={'k':str(k),'error':type(exc).__name__+': '+str(exc)}
    out['records'].append(row)
    print(json.dumps({k:v for k,v in row.items() if k not in ['upper','lower','equilibria_numerical']}),flush=True)
    if len(row.get('upper_brackets',[]))>=3 and row.get('lower_brackets'):break
for k in [F(1,10000),F(1,5)]:
    old=json.loads((ROOT/'reversible_finite_hopf.json').read_text())
    record=next(r for r in old['records'] if r['field']['k']==str(k))
    model=field(k)
    for side,key in [(1,'upper_brackets'),(-1,'lower_brackets')]:
        for bracket in record[key]:
            for radius in bracket:
                p=half(model,radius,side,1,2e-13);q=half(model,radius,side,-1,2e-13)
                original=next(r for r in record['upper' if side==1 else 'lower'] if r['r']==radius)
                if p['status']!=q['status'] or p['status']!='half_passage':raise ArithmeticError('tight control failed')
                diff=p['coordinate']-q['coordinate']
                assert diff*original['difference']>0
                out['tightened_controls'].append({'k':str(k),'side':side,'r':radius,'old_difference':original['difference'],'tight_difference':diff,'change':diff-original['difference']})
out['wall_seconds']=time.perf_counter()-begin
(ROOT/'reversible_finite_hopf_followup.json').write_text(json.dumps(out,indent=2)+'\n')
print('tight controls',len(out['tightened_controls']),flush=True)
