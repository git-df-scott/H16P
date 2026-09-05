"""Explicit bounded Stage-2 experiment. See STAGED_RUN_2026_09_05.md.

Commands run distinct blocks; every ODE request is charged by budget.py.
Root corrections use the complete curved section. 'fold' solves the two
equations log(|R/r|)=0 and d log(|R/r|)/d log|r|=0: exactly D=D_r=0
at a finite fixed point, equivalently the bounded q-section equations.
No optimizer termination is labelled a fold without both residual gates.
"""
import json,math,sys
from fractions import Fraction as F
from pathlib import Path
import numpy as np
from budget import call,HERE,count,LEDGER

def parameters(c,K):
    c,K=F(str(c)),F(str(K))
    return dict(c=str(c),alpha=str(-(K+42)/(F(11,5)*c-1)),beta='0')

def correct(c,K,r,purpose,max_calls=7):
    records=[];z=math.log(abs(r)); sign=1 if r>0 else -1
    for _ in range(max_calls):
        result=call(parameters(c,K)|dict(r=sign*math.exp(z)),purpose)
        records.append(result)
        if result['status']!='NUMERICAL_ONLY':return dict(status='UNRESOLVED',records=records)
        d=result['log_displacement'];der=result['log_displacement_derivative']
        if abs(d)<2e-8:
            return dict(status='ACCEPTED_NUMERICAL_ROOT',c=str(c),K=str(K),root=result,records=records)
        if abs(der)<1e-6:break
        step=max(-1.,min(1.,d/der));z-=step
        if z>30 or z<math.log(.02):break
    return dict(status='CORRECTOR_NOT_CONVERGED',records=records)

def validation():
    rows=json.loads((HERE/'controls.json').read_text())
    out=[]
    for j in (0,3):
        a=rows[j]['curved'];r=a['r']; h=1e-4
        base={k:a[k] for k in ('r','c','alpha','beta')}
        lo=call(base|dict(r=r*math.exp(-h)), 'derivative_control_minus')
        hi=call(base|dict(r=r*math.exp(h)), 'derivative_control_plus')
        fd=(hi['log_displacement']-lo['log_displacement'])/(2*h)
        out.append(dict(name=rows[j]['name'],analytic=a['log_displacement_derivative'],
                        finite_difference=fd,error=fd-a['log_displacement_derivative']))
    # Correct rounded pilot seed, then compare complete curved return to old-section multiplier.
    remote=correct('7/10','6/5',rows[1]['curved']['r'],'rounded_remote_seed_correct',4)
    out.append(remote)
    (HERE/'derivative_validation.json').write_text(json.dumps(out,indent=2))
    assert max(abs(x['error']) for x in out[:2])<2e-5

def terminal():
    controls=json.loads((HERE/'controls.json').read_text())
    rO=controls[2]['curved']['r'];rR=controls[3]['curved']['r']
    rows=[]
    # Controlled continuation beyond old cap, stopping each failed nest separately.
    activeO=activeR=True
    for c in ('0.9301','0.935','0.94','0.945','0.95','0.955','0.96','0.962','0.963','0.964','0.966','0.968'):
        row=dict(c=c,K='1/64')
        if activeO:
            row['origin']=correct(c,'1/64',rO,'terminal_origin',6)
            activeO=row['origin']['status']=='ACCEPTED_NUMERICAL_ROOT'
            if activeO:rO=row['origin']['root']['r']
        if activeR:
            row['remote']=correct(c,'1/64',rR,'terminal_remote',7)
            activeR=row['remote']['status']=='ACCEPTED_NUMERICAL_ROOT'
            if activeR:rR=row['remote']['root']['r']
        rows.append(row)
        (HERE/'terminal_path.json').write_text(json.dumps(rows,indent=2))
        if not activeO and not activeR:break

def profile():
    rows=[]
    # Written finite design, not parameter-space coverage. Six K values span
    # the open budget interval; nine c strata and three radial amplitudes.
    for c in ('0.6','0.7','0.825','0.93','0.95','0.965','1','1.2','1.5'):
      for K in ('1/512','1/64','1/16','1/4','3/5','119/100'):
        for r in (2.,20.,20000.):
            out=call(parameters(c,K)|dict(r=r),'origin_displacement_profile')
            rows.append(dict(c=c,K=K,r=r,result=out))
        (HERE/'profiles.json').write_text(json.dumps(rows,indent=2))

def fold():
    rows=[]
    # The frozen endpoint maximum seeds augmented Newton. Additional seeds
    # across K are local attempts, not a cover of disconnected fold sheets.
    for K in ('1/64','1/512','1/16','1/4','3/5','119/100'):
        z=np.array([math.log(29.4),.9301]);history=[];status='ITERATION_LIMIT'
        for iteration in range(4):
            def evaluate(v,purpose):
                out=call(parameters(format(v[1],'.16g'),K)|dict(r=math.exp(v[0])),purpose)
                if out['status']!='NUMERICAL_ONLY':raise ValueError(out.get('error'))
                return np.array([out['log_displacement'],out['log_displacement_derivative']]),out
            try:
                f,out=evaluate(z,'augmented_fold_base')
                history.append(dict(log_r=float(z[0]),c=float(z[1]),residual=f.tolist(),result=out))
                if abs(f[0])<1e-8 and abs(f[1])<1e-7:
                    status='CANDIDATE_REQUIRES_NONDEGENERACY_REPLAY';break
                h=np.array([2e-4,2e-5]);jac=np.empty((2,2))
                for j in range(2):
                    zz=z.copy();zz[j]+=h[j]
                    fj,_=evaluate(zz,'augmented_fold_difference')
                    jac[:,j]=(fj-f)/h[j]
                step=np.linalg.solve(jac,-f)
                factor=min(1.,.65/max(abs(step[0]),1e-100),.018/max(abs(step[1]),1e-100))
                next_z=z+factor*step
                history[-1].update(jacobian=jac.tolist(),newton_step=step.tolist(),damping=factor)
                if not (.55<next_z[1]<1.5 and math.log(.05)<next_z[0]<math.log(1e8)):
                    status='DOMAIN_BOUNDARY';break
                z=next_z
            except Exception as exc:status='UNRESOLVED: '+str(exc);break
        rows.append(dict(K=K,status=status,history=history))
        (HERE/'fold_attempts.json').write_text(json.dumps(rows,indent=2))

if __name__=='__main__':
    {'validation':validation,'terminal':terminal,'profile':profile,'fold':fold}[sys.argv[1]]()
