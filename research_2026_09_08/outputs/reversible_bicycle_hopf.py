"""Finite Hopf tests from two a<-2 reversible base shapes."""
from fractions import Fraction as F
from pathlib import Path
import json,time
import numpy as np
from reversible_finite_hopf import field,profile,count
ROOT=Path(__file__).resolve().parent
out={'evidence':'NUM samples only; returns outside base annulus can fail and remain unresolved','records':[]};start=time.perf_counter()
for a,m in [(F(-3),F(9,40)),(F(-9,4),F(63,200))]:
 for k in [F(1,100),F(1,10),F(3,10),F(3,5)]:
    try:
        model=field(k,m,a,F(1) if a==-3 else F(1,3))
        up=profile(model,np.geomspace(.02,1e6,41),1,2e-11)
        lo=profile(model,np.geomspace(.02,1e6,41),-1,2e-11)
        row={'field':model[0],'equilibria':model[1],'upper':up,'lower':lo,'upper_brackets':count(up),'lower_brackets':count(lo),'unresolved':sum('difference' not in r for r in up+lo)}
    except (ArithmeticError,AssertionError) as e:row={'a':str(a),'k':str(k),'error':str(e)}
    out['records'].append(row);out['wall_seconds']=time.perf_counter()-start
    (ROOT/'reversible_bicycle_hopf.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({'a':str(a),'k':str(k),**{key:row[key] for key in ['upper_brackets','lower_brackets','unresolved','error'] if key in row}}),flush=True)
    if len(row.get('upper_brackets',[]))>=3 and row.get('lower_brackets'):raise SystemExit('Candidate requires independent check')
