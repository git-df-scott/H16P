"""Adaptive follow-up to the frozen Stage-2 blocks, charged within 400."""
import json,math,sys
import numpy as np
from budget import call,HERE,count,LEDGER
from run_kkl import parameters,correct

def terminal():
    last=json.loads((HERE/'terminal_path.json').read_text())[-1]
    rO=last['origin']['root']['r'];rR=last['remote']['root']['r']
    start=count(LEDGER);rows=[]
    for c in ('0.9682','0.9683','0.9684'):
        if count(LEDGER)-start>=22:break
        row=dict(c=c,K='1/64')
        for name,r in (('origin',rO),('remote',rR)):
            remaining=22-(count(LEDGER)-start)
            if remaining<=0:break
            out=correct(c,'1/64',r,'terminal_neutrality_followup',min(6,remaining))
            row[name]=out
            if out['status']=='ACCEPTED_NUMERICAL_ROOT':
                if name=='origin':rO=out['root']['r']
                else:rR=out['root']['r']
        rows.append(row)
        (HERE/'terminal_followup.json').write_text(json.dumps(rows,indent=2))
        if any(x.get('status')!='ACCEPTED_NUMERICAL_ROOT' for k,x in row.items() if k in ('origin','remote')):break

def folds():
    attempts=json.loads((HERE/'fold_attempts.json').read_text())
    h=attempts[1]['history'][-1]
    z=np.array([h['log_r'],h['c']])+h['damping']*np.array(h['newton_step'])
    # Second seed is the actual accepted path cycle with smallest |P'-1|.
    last=json.loads((HERE/'terminal_path.json').read_text())[-1]['origin']['root']
    seeds=[('1/512',z),('1/64',np.array([math.log(last['r']),last['c']]))]
    rows=[]
    for K,z in seeds:
        history=[];status='ITERATION_LIMIT'
        for _ in range(4):
            def ev(v,purpose):
                a=call(parameters(format(v[1],'.16g'),K)|dict(r=math.exp(v[0]),tol=5e-12),purpose)
                if a['status']!='NUMERICAL_ONLY':raise ValueError(a.get('error'))
                # Remove the automatic r^2 vanishing near the focus.
                scale=(1-a['q'])**2
                return np.array([a['log_displacement'],a['log_displacement_derivative']])/scale,a
            try:
                f,a=ev(z,'fold_followup_base');history.append(dict(c=float(z[1]),log_r=float(z[0]),residual=f.tolist(),result=a))
                if np.linalg.norm(f,np.inf)<1e-8:
                    status='CANDIDATE_REQUIRES_NONDEGENERACY_REPLAY';break
                J=np.empty((2,2));steps=[1e-4,1e-6]
                for j in range(2):
                    v=z.copy();v[j]+=steps[j];fj,_=ev(v,'fold_followup_difference');J[:,j]=(fj-f)/steps[j]
                step=np.linalg.solve(J,-f)
                damping=min(1,.6/max(abs(step[0]),1e-99),.005/max(abs(step[1]),1e-99))
                history[-1].update(jacobian=J.tolist(),step=step.tolist(),damping=damping)
                z=z+damping*step
                if not (.55<z[1]<1.5 and math.log(.05)<z[0]<math.log(1e8)):
                    status='DOMAIN_BOUNDARY';break
            except Exception as exc:status='UNRESOLVED: '+str(exc);break
        rows.append(dict(K=K,status=status,history=history))
        (HERE/'fold_followup.json').write_text(json.dumps(rows,indent=2))

if __name__=='__main__':{'terminal':terminal,'folds':folds}[sys.argv[1]]()
