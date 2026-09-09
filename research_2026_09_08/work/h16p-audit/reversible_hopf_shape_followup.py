"""Bounded two-parameter followup near observed Hopf coefficient changes."""
from pathlib import Path
from fractions import Fraction as F
import json,time
import numpy as np
from reversible_finite_hopf import field,profile,count
ROOT=Path(__file__).resolve().parent
out={'evidence':'NUM finite shape and section sample only; no exclusion theorem','records':[]}
begin=time.perf_counter()
for k in [F(1,5),F(2,5),F(1,2),F(7,10)]:
 for m in [F(7,5),F(3,2),F(8,5),F(17,10)]:
    try:
        model=field(k,m)
        upper=profile(model,np.geomspace(.02,1e8,49),1,2e-11)
        lower=profile(model,np.geomspace(.2,1e8,33),-1,2e-11)
        row={'field':model[0],'equilibria_numerical':model[1],'upper':upper,'lower':lower,'upper_brackets':count(upper),'lower_brackets':count(lower),'unresolved':sum('difference' not in r for r in upper+lower)}
    except (ArithmeticError,AssertionError) as exc:row={'k':str(k),'m':str(m),'error':type(exc).__name__+': '+str(exc)}
    out['records'].append(row)
    out['wall_seconds']=time.perf_counter()-begin
    (ROOT/'reversible_hopf_shape_followup.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({'k':str(k),'m':str(m),**{key:row[key] for key in ['upper_brackets','lower_brackets','unresolved','error'] if key in row}}),flush=True)
    if len(row.get('upper_brackets',[]))>=3 and row.get('lower_brackets'):
        print('CANDIDATE REQUIRES INDEPENDENT CHECK',flush=True);raise SystemExit
