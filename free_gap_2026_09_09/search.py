"""Release upper log gaps: constrained outer-sign optimization.

NUM only. A candidate needs four alternating upper sign brackets and a lower
bracket in ONE field. Existing witnesses are inequalities, not fixed roots.
All five field parameters and seven witness positions may vary independently.
"""
from pathlib import Path
import sys,json,time
import numpy as np
from scipy.optimize import minimize
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent/'moving_cycle_2026_09_09'))
import engine as eng

class Budget(Exception):pass
class Candidate(Exception):pass

def unpack(z):
    a,b,ell,u1,u2=z[:5];amp=np.exp(ell)
    q=np.array([a,b,-amp/eng.TAU,amp*u1/eng.TAU,amp*u2/eng.TAU])
    Q=np.zeros((5,5));Q[0,0]=Q[1,1]=1
    Q[2:,2]=q[2:];Q[3,3]=Q[4,4]=amp/eng.TAU
    return q,Q

def run(name):
    start=time.perf_counter();calls=0;evals=0;history=[];best=None;last=None
    source='amplitude_resumed_result.json' if name=='strong' else 'amplitude_result.json'
    state=json.loads((HERE.parent/'moving_cycle_2026_09_09'/source).read_text())['accepted'][-1]
    q0=np.array(state['q']);r=np.array(state['roots']);amp=-q0[2]*eng.TAU
    s=np.array([r[0]-.12,(r[0]+r[1])/2,(r[1]+r[2])/2,r[2]+.18,r[2]+1.2,r[3]-.05,r[3]+.05])
    z0=np.r_[q0[:2],np.log(amp),q0[3]*eng.TAU/amp,q0[4]*eng.TAU/amp,s]
    signs=np.array([1,-1,1,-1,1,-1,1])
    rows=[]
    for i,si in enumerate(s):
        p=eng.pair(si,q0,1 if i<5 else -1,rtol=3e-13,regularized=True)
        if p['status']!='passed':raise RuntimeError('initial witness domain')
        rows.append(p)
    assert all(signs[i]*rows[i]['D']>0 for i in [0,1,2,3,5,6])
    # Fixed scales; objective cannot shrink by choosing a growing denominator.
    scales=np.array([abs(p['D']) for p in rows]);scales[4]=max(scales[:4])
    lower=np.r_[max(-1.97,z0[0]-.3),max(.06,z0[1]-.27),z0[2]-1.5,z0[3]-4,z0[4]-4,
                 np.maximum(np.log(.5)+.025,s[:5]-1.3),s[5:]-5]
    upper=np.r_[min(-.55,z0[0]+.5),min(1.8,z0[1]+.7),z0[2]+1.5,z0[3]+4,z0[4]+4,
                 s[:5]+2.,s[5:]+6]
    # Initial calls are charged and retained.
    ledger=open(HERE/f'{name}_evaluations.jsonl','w',buffering=1)
    for p in rows:
        calls+=1;p['call']=calls;p['phase']='initial';ledger.write(json.dumps(p)+'\n')
    fixed_indices=[0,1,2,3,5,6]
    def gates(z):
        b=z[1];ss=z[5:]
        vals=np.r_[np.diff(ss[:5])-np.array([.08,.08,.08,.2]),ss[6]-ss[5]-.025,
                   ss[5]-np.log((2-b)/(2*b))-.05]
        J=np.zeros((6,12))
        for i in range(4):J[i,5+i]=-1;J[i,6+i]=1
        J[4,10]=-1;J[4,11]=1;J[5,10]=1;J[5,1]=1/(2-b)+1/b
        return vals,J
    def evaluate(z):
        nonlocal calls,evals,last,best
        if last is not None and np.array_equal(z,last['z']):return last
        if calls+7>700:raise Budget()
        evals+=1;q,Q=unpack(z);vals=np.zeros(7);J=np.zeros((7,12));rs=[];valid=True
        for i,si in enumerate(z[5:]):
            calls+=1
            try:p=eng.pair(si,q,1 if i<5 else -1,rtol=3e-13,regularized=True)
            except (ValueError,FloatingPointError,OverflowError) as e:p={'status':'unresolved','message':str(e),'s':float(si),'q':q.tolist()}
            p['call']=calls;p['evaluation']=evals;ledger.write(json.dumps(p)+'\n');ledger.flush();rs.append(p)
            if p['status']!='passed':valid=False;continue
            vals[i]=p['D'];g=np.array(p['derivatives']);J[i,:5]=g[1:]@Q;J[i,5+i]=g[0]
        norm=signs*vals/scales;nJ=signs[:,None]*J/scales[:,None]
        gg,gJ=gates(z)
        # Invalid returns are explicit infeasible evaluations, never zero returns.
        cons=np.r_[norm[fixed_indices]-.05,gg]
        cJ=np.r_[nJ[fixed_indices],gJ]
        if not valid:cons[:6]=-1e6;cJ[:6]=0
        obj=-norm[4] if valid else 1e6
        grad=-nJ[4] if valid else np.zeros(12)
        feasible=bool(valid and min(cons)>=-1e-5)
        record={'evaluation':evals,'calls':calls,'z':z.tolist(),'q':q.tolist(),'coefficients':eng.coefficients(q).tolist(),
                'D':vals.tolist(),'normalized_signed_D':norm.tolist(),'valid':valid,'feasible':feasible,
                'minimum_constraint':float(min(cons)),'objective':float(obj)}
        history.append(record)
        if feasible and (best is None or obj<best['objective']):
            best=record;print(name,'best',evals,'outer normalized',norm[4],'old minimum',min(norm[fixed_indices]),'calls',calls,flush=True)
            (HERE/f'{name}_best.json').write_text(json.dumps(best,indent=2)+'\n')
            if norm[4]>.1:raise Candidate()
        last={'z':np.array(z),'objective':obj,'gradient':grad,'constraints':cons,'jacobian':cJ,'rows':rs}
        return last
    def fun(z):return evaluate(z)['objective']
    def jac(z):return evaluate(z)['gradient']
    def con(z):return evaluate(z)['constraints']
    def cjac(z):return evaluate(z)['jacobian']
    status='not_started';optimizer={}
    try:
        initial=evaluate(z0)
        result=minimize(fun,z0,jac=jac,method='SLSQP',bounds=list(zip(lower,upper)),
            constraints={'type':'ineq','fun':con,'jac':cjac},options={'maxiter':45,'ftol':1e-9,'disp':False})
        status='optimizer_finished';optimizer={'success':bool(result.success),'message':result.message,'nit':int(result.nit),'x':result.x.tolist(),'objective':float(result.fun)}
    except Budget:status='700 paired-attempt cap'
    except Candidate:status='positive outer witness; independent validation required'
    finally:
        ledger.close()
        out={'evidence':'NUM; constrained sign search, not exhaustive','name':name,'source':source,'status':status,
             'initial_z':z0.tolist(),'fixed_scales':scales.tolist(),'signs':signs.tolist(),'lower_bounds':lower.tolist(),
             'upper_bounds':upper.tolist(),'calls':calls,'evaluations':evals,'best':best,'history':history,'optimizer':optimizer,
             'counts':eng.COUNTS,'wall_seconds':time.perf_counter()-start}
        (HERE/f'{name}_result.json').write_text(json.dumps(out,indent=2)+'\n')
        print(name,'STOP',status,optimizer,'calls',calls,flush=True)

if __name__=='__main__':run(sys.argv[1])
