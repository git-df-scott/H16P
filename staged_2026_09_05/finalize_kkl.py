"""Replay source for the final adaptive refinement and independent checks.

The committed ledger is full. This script intentionally refuses new charged
calls there. A separately budgeted reproduction must use a fresh working copy
and preserve its own ledger; never remove historical rows to evade the cap.
"""
import json,math,sys
from fractions import Fraction as F
import numpy as np
from budget import call,HERE
from run_kkl import parameters

def refinement():
    h=json.loads((HERE/'fold_followup.json').read_text())[0]['history'][-1]
    z=np.array([h['log_r'],h['c']])+h['damping']*np.array(h['step'])
    rows=[]
    for _ in range(3):
        def ev(v,p):
            a=call(parameters(format(v[1],'.17g'),'1/512')|dict(r=math.exp(v[0]),tol=2e-12),p)
            if a['status']!='NUMERICAL_ONLY':raise ValueError(a)
            scale=(1-a['q'])**2
            return np.array([a['log_displacement'],a['log_displacement_derivative']])/scale,a
        f,a=ev(z,'low_K_fold_final_refinement')
        row=dict(c=float(z[1]),log_r=float(z[0]),residual=f.tolist(),result=a);rows.append(row)
        if max(abs(f))<2e-10:break
        J=np.empty((2,2))
        for j,hj in enumerate([1e-4,1e-6]):
            v=z.copy();v[j]+=hj;fj,_=ev(v,'low_K_fold_final_difference');J[:,j]=(fj-f)/hj
        step=np.linalg.solve(J,-f);row.update(jacobian=J.tolist(),step=step.tolist());z+=step
        (HERE/'fold_refined_candidate.json').write_text(json.dumps(rows,indent=2))
    (HERE/'fold_refined_candidate.json').write_text(json.dumps(rows,indent=2))

def verification():
    h=json.loads((HERE/'fold_refined_candidate.json').read_text())[-1]
    c=F(str(h['c']));r=math.exp(h['log_r']);rows=[]
    requests=[(c,r,'fold_independent_cartesian','cartesian_check.py'),
      (c+F(1,10**6),4.,'fold_pair_lower_sign','cartesian_check.py'),
      (c+F(1,10**6),r,'fold_pair_middle_sign','cartesian_check.py'),
      (c+F(1,10**6),12.,'fold_pair_upper_sign','cartesian_check.py'),
      (c-F(1,10**6),r,'fold_absent_side_middle','cartesian_check.py'),
      (c+F(1,10**6),20000.,'fold_pair_outer_profile','compact_return.py')]
    for cv,rv,p,e in requests:
        rows.append(dict(purpose=p,result=call(parameters(cv,'1/512')|dict(r=rv),p,engine=e)))
    last=json.loads((HERE/'terminal_followup.json').read_text())[1]
    for nest in ('origin','remote'):
        a=last[nest]['root'];p='terminal_independent_'+nest
        rows.append(dict(purpose=p,result=call(parameters('0.9683','1/64')|dict(r=a['r']),p,engine='cartesian_check.py')))
    (HERE/'final_verification.json').write_text(json.dumps(rows,indent=2))

if __name__=='__main__':{'refinement':refinement,'verification':verification}[sys.argv[1]]()
