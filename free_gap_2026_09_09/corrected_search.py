"""Free-gap outer-sign ascent with exact numerical root correction.

Fix e0 away from zero. Controls are a,b,g1,g2,d. Solve the upper root
equalities for perturbation ratios u1,u2 and the first log height t.
This releases both prior gap constraints. Lower roots/witnesses are checked
before accepting a field. Nonvalidated numerical experiment only.
"""
from pathlib import Path
import sys,json,time
import numpy as np
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent/'moving_cycle_2026_09_09'))
import engine as eng
from center_geometry import center_coefficients,center_jacobian
if '--boundary' in sys.argv:
    import extended_engine as eng

class Budget(Exception):pass
class Gate(Exception):pass

def run(name):
    start=time.perf_counter();calls=0;events=[];states=[]
    file='amplitude_resumed_result.json' if name=='strong' else 'amplitude_result.json'
    seed=json.loads((HERE.parent/'moving_cycle_2026_09_09'/file).read_text())['accepted'][-1]
    qseed=np.array(seed['q']);ss=np.array(seed['roots']);amp=-eng.TAU*qseed[2]
    c=np.r_[qseed[:2],ss[1]-ss[0],ss[2]-ss[1],1.0]
    w=np.r_[qseed[3:]*eng.TAU/amp,ss[0]];s4=ss[3]
    scales=np.array([.1,.1,.2,.2,.5])
    lows=np.array([-1.97,.06,.08,.08,.2]);highs=np.array([-.7,1.6,2.,2.,4.])
    cap=900
    if '--boundary' in sys.argv:
        old=json.loads((HERE/f'{name}_corrected_result.json').read_text())['states'][-1]
        c=np.array(old['c']);w=np.array(old['w']);s4=old['s4'];name+='_boundary'
        lows=np.array([-2.15,.02,.08,.08,.2]);cap=500
    if '--divide-center' in sys.argv:
        name+='_divided';cap=500
    ledger=open(HERE/f'{name}_corrected_evaluations.jsonl','w',buffering=1)
    def ev(s,q,side=1):
        nonlocal calls
        if calls>=cap:raise Budget()
        calls+=1
        try:p=eng.pair(s,q,side,rtol=3e-13,regularized=True)
        except (ValueError,FloatingPointError,OverflowError) as e:p={'status':'unresolved','message':str(e),'q':q.tolist(),'s':float(s)}
        p['call']=calls;ledger.write(json.dumps(p)+'\n');ledger.flush()
        if p['status']!='passed':raise Gate('return gate')
        return p
    def makeq(c,w):return np.r_[c[:2],-amp/eng.TAU,w[:2]*amp/eng.TAU]
    def correct(c,w):
        w=w.copy()
        for iteration in range(10):
            q=makeq(c,w);roots=w[2]+np.array([0,c[2],c[2]+c[3]])
            if roots[0]<np.log(.5)+.015:raise Gate('upper inner section boundary')
            if np.max(abs(w[:2]))>25:raise Gate('ratio box')
            rows=[ev(s,q) for s in roots];G=np.array([r['derivatives'] for r in rows]);D=np.array([r['D'] for r in rows])
            J=np.c_[G[:,4:6]*amp/eng.TAU,G[:,0]]
            if np.any(abs(G[:,0])<2e-10):raise Gate('numerically weak upper derivative')
            if np.max(abs(D/G[:,0]))<2e-6:
                if np.any(G[:,0]*np.array([-1,1,-1])<=0):raise Gate('upper cycle orientation changed')
                A=np.zeros((3,5));A[:,:2]=G[:,1:3];A[1:,2]=G[1:,0];A[2,3]=G[2,0]
                W=-np.linalg.solve(J,A)
                Q=np.zeros((5,5));Q[0,0]=Q[1,1]=1;Q[3:]=W[:2]*amp/eng.TAU
                p=ev(roots[2]+c[4],q);g=np.array(p['derivatives'])
                star_gradient=W[2]+np.array([0,0,1,1,1])
                grad=g[1:]@Q+g[0]*star_gradient
                objective=p['D']/amp;distance=None
                if '--divide-center' in sys.argv:
                    center=center_coefficients(c[0],c[1])
                    delta=(eng.coefficients(q)[2:]-center)/amp
                    distance=float(np.linalg.norm(delta))
                    if distance<.001:raise Gate('center-distance floor')
                    delta_jac=eng.TAU*Q[2:]/amp
                    delta_jac[:,:2]-=center_jacobian(c[0],c[1])/amp
                    dr=delta@delta_jac/distance
                    grad=grad/distance-p['D']*dr/(distance*distance)
                    objective/=distance
                return dict(c=c.copy(),w=w,q=q,roots=roots,upper=rows,outer=p,W=W,Q=Q,gradient=grad,
                            objective=objective,center_distance=distance,condition=float(np.linalg.cond(J)),iterations=iteration+1)
            dw=np.linalg.solve(J,-D)
            factor=min(1.,.15/max(abs(dw[2]),1e-100),.5/max(np.max(abs(dw[:2])),1e-100))
            w+=factor*dw
        raise Gate('root corrector did not converge')
    def complete(state,s4):
        q=state['q']
        for _ in range(8):
            p=ev(s4,q,-1);g=np.array(p['derivatives'])
            if abs(g[0])<2e-10:raise Gate('weak lower derivative')
            step=-p['D']/g[0]
            if abs(step)<2e-6:break
            s4+=np.clip(step,-.5,.5)
        else:raise Gate('lower corrector did not converge')
        intervals=[]
        for i,s in enumerate(np.r_[state['roots'],s4]):
            side=1 if i<3 else -1;delta=min(.03,min(state['c'][2:4])*.15) if i<3 else .025
            l=ev(s-delta,q,side);r=ev(s+delta,q,side)
            if l['D']*r['D']>=0 or min(abs(l['D']),abs(r['D']))<3e-12:raise Gate('old sign bracket lost')
            intervals.append({'side':side,'interval':[float(s-delta),float(s+delta)],'D':[l['D'],r['D']]})
        state.update(lower=p,s4=s4,witnesses=intervals)
        return state
    def serialize(state):
        out={k:(v.tolist() if isinstance(v,np.ndarray) else v) for k,v in state.items()}
        out['calls']=calls;out['coefficients']=eng.coefficients(state['q']).tolist()
        return out
    status='not_started';rho=.3
    try:
        current=complete(correct(c,w),s4);states.append(serialize(current))
        for iteration in range(35):
            grad=current['gradient']/amp
            # Project the ascent direction onto active box faces.
            g=grad*scales
            g[(current['c']<=lows+1e-9)&(g<0)]=0
            g[(current['c']>=highs-1e-9)&(g>0)]=0
            if np.linalg.norm(g)<1e-9:status='projected gradient small';break
            direction=scales*g/np.linalg.norm(g)
            accepted=False
            for attempt in range(7):
                cn=np.clip(current['c']+rho*direction,lows,highs);dc=cn-current['c']
                try:
                    trial=correct(cn,current['w']+current['W']@dc)
                    if trial['objective']<=current['objective']+1e-8:raise Gate('no objective improvement')
                    oldg=np.array(current['lower']['derivatives'])
                    pred=current['s4']-oldg[1:]@(trial['q']-current['q'])/oldg[0]
                    trial=complete(trial,pred)
                    current=trial;states.append(serialize(current));accepted=True
                    print(name,'accepted',len(states)-1,'objective',current['objective'],'c',current['c'],
                          'roots',current['roots'],'lower',current['s4'],'calls',calls,flush=True)
                    (HERE/f'{name}_corrected_checkpoint.json').write_text(json.dumps({'states':states,'events':events,'calls':calls},indent=2)+'\n')
                    if current['outer']['D']>3e-10:status='positive outer witness; verify candidate';break
                    rho=min(.35,rho*1.35);break
                except (Gate,np.linalg.LinAlgError) as e:
                    events.append({'iteration':iteration,'attempt':attempt,'c':cn.tolist(),'rho':rho,'message':str(e),'calls':calls})
                    print(name,'reject',str(e),'rho',rho,'calls',calls,flush=True);rho*=.5
            if status=='positive outer witness; verify candidate':break
            if not accepted or rho<1e-6:status='persistent trust-region gate';break
        else:status='35 accepted-iteration cap'
    except Budget:status=f'{cap} paired-attempt cap'
    finally:
        ledger.close()
        result={'evidence':'NUM, corrected free-gap ascent; no global exclusion','name':name,'status':status,
                'fixed_amplitude':amp,'calls':calls,'states':states,'events':events,'counts':eng.COUNTS,
                'wall_seconds':time.perf_counter()-start,'control_bounds':[lows.tolist(),highs.tolist()],
                'objective_definition':'D/(abs(e0)*center_distance)' if '--divide-center' in sys.argv else 'D/abs(e0)'}
        (HERE/f'{name}_corrected_result.json').write_text(json.dumps(result,indent=2)+'\n')
        print(name,'STOP',status,'calls',calls,flush=True)

if __name__=='__main__':run(sys.argv[1])
