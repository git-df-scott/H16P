"""Two bounded moving-root curves: shape and finite perturbation amplitude.

Upper log root gaps are held at their seed values, a common translation is
free, and the lower root moves independently. These are explicit local slices,
not a five-dimensional coverage claim. All calls, including failed trials,
are retained and charged to a 500-pair branch cap.
"""
import sys
from engine import *

class Budget(Exception):pass
class Gate(Exception):pass

seed=json.loads((ROOT/'seed_replay.json').read_text())
BASE=np.array([r['s'] for r in seed['roots']])

def run(kind,resume=False):
    started=time.perf_counter();calls=0;accepted=[];events=[]
    prefix=kind+'_resumed' if resume else kind
    tail_gap=20.
    if resume:
        prior=json.loads((ROOT/f'{kind}_result.json').read_text())
        boundary=json.loads((ROOT/f'{kind}_boundary.json').read_text())
        calls=prior['calls']+boundary['counts']['pairs']
        accepted=prior['accepted'].copy();events=prior['events'].copy()
        tail_gap=.85*(boundary['transition_interval'][0]-accepted[-1]['roots'][2])
    initial_calls=calls
    ledger=open(ROOT/f'{prefix}_evaluations.jsonl','w',buffering=1)
    def evaluate(s,q,side=1,rtol=3e-13):
        nonlocal calls
        if calls>=500:raise Budget()
        calls+=1
        try:p=pair(s,q,side,rtol=rtol,regularized=resume)
        except (ValueError,FloatingPointError,OverflowError) as e:
            p={'status':'unresolved','message':str(e),'s':float(s),'q':q.tolist(),'side':side}
        p['evaluation']=calls;ledger.write(json.dumps(p)+'\n');ledger.flush()
        if p['status']!='passed':raise Gate(p.get('message','return gate'))
        return p
    def qfun(lam,w):
        q=SEED.copy()
        if kind=='shape':
            q[0]+=lam;q[3:]=w[:2];A=1.;qlam=np.array([1.,0,0,0,0])
        else:
            A=np.exp(lam);q[2:]*=A;q[3:]=A*w[:2];qlam=np.r_[0.,0.,q[2:]]
        return q,A,qlam
    def upper(lam,w):
        q,A,qlam=qfun(lam,w);rr=[evaluate(s+w[2],q) for s in BASE[:3]]
        g=np.array([r['derivatives'] for r in rr])
        J=np.c_[A*g[:,4],A*g[:,5],g[:,0]]
        return q,rr,g,J,qlam
    def correct(lam,w):
        w=w.copy()
        for iteration in range(8):
            q,rr,g,J,qlam=upper(lam,w);D=np.array([r['D'] for r in rr])
            if np.any(abs(g[:,0])<2e-9):raise Gate('upper simple-root derivative below discovery floor')
            if np.max(abs(D/g[:,0]))<3e-6:
                return w,q,rr,g,J,qlam,iteration+1
            dw=np.linalg.solve(J,-D)
            # Trust limits in ratio coefficients and log-height translation.
            factor=min(1.,.2/max(abs(dw[2]),1e-100),.5/max(np.max(abs(dw[:2])),1e-100))
            w+=factor*dw
            if abs(w[2])>3:raise Gate('common log-height shift beyond chosen slice box')
        raise Gate('upper corrector did not converge')
    def lower(q,s):
        for iteration in range(8):
            r=evaluate(s,q,-1);g=np.array(r['derivatives'])
            if abs(g[0])<2e-9:raise Gate('lower simple-root derivative below discovery floor')
            ds=-r['D']/g[0]
            if abs(ds)<3e-6:return s,r
            s+=np.clip(ds,-.5,.5)
            if not np.isfinite(s):raise Gate('lower root nonfinite')
        raise Gate('lower corrector did not converge')
    def witnesses(q,ss):
        out=[]
        for i,s in enumerate(ss):
            side=1 if i<3 else -1;delta=.035 if i<3 else .02
            lo=evaluate(s-delta,q,side);hi=evaluate(s+delta,q,side)
            if lo['D']*hi['D']>=0:raise Gate('lost local sign bracket')
            if min(abs(lo['D']),abs(hi['D']))<3e-12:raise Gate('witness displacement below discovery floor')
            out.append({'side':side,'interval':[s-delta,s+delta],'D':[lo['D'],hi['D']]})
        return out
    def outer(q,s3):
        nonlocal tail_gap
        rows=[];stationary=[]
        gaps=sorted(set([.2,1.,min(4.,tail_gap*.4),min(10.,tail_gap*.7),tail_gap]))
        for gap in gaps:
            try:r=evaluate(s3+gap,q)
            except Gate as e:
                if not resume:
                    rows.append({'s':float(s3+gap),'status':'unresolved','message':str(e)});break
                # The observed y=0 crossing bounds this particular chart. Keep
                # its unresolved exterior in the ledger and move the probe in.
                events.append({'kind':'outer_domain_trim','s':float(s3+gap),'message':str(e),'calls':calls})
                lastgap=rows[-1]['s']-s3 if rows else 0.
                for _ in range(3):
                    gap=(gap+lastgap)/2
                    try:r=evaluate(s3+gap,q);tail_gap=gap;break
                    except Gate:continue
                else:raise Gate('outer domain remained unresolved after three inward probes')
            rows.append(r)
            if len(rows)>1:
                left=rows[-2]
                if left['derivatives'][0]*r['derivatives'][0]<0:
                    def f(s):return evaluate(s,q)['derivatives'][0]
                    sx=brentq(f,left['s'],r['s'],xtol=2e-6)
                    st=evaluate(sx,q,rtol=8e-14);stationary.append(st)
        return rows,stationary

    lam=0.;w=np.r_[SEED[3:],0.];s4=BASE[3]
    if resume:
        lam=accepted[-1]['lambda'];w=np.array(accepted[-1]['w']);s4=accepted[-1]['roots'][3]
    step=.04 if kind=='shape' else .5
    # Seed state is already independently replayed; initial derivatives predict
    # the first move along the explicitly defined slice.
    q,rr,g,J,qlam=upper(lam,w)
    tangent=-np.linalg.solve(J,g[:,1:]@qlam)
    status='accepted-step cap';failed_attempts=0
    try:
        for attempt in range(100):
            if len(accepted)>=50:break
            target=lam+step
            if kind=='shape' and SEED[0]+target>=-.85:
                status='shape range cap';break
            try:
                wn,qn,rn,gn,Jn,qnlam,iters=correct(target,w+step*tangent)
                dq=qn-q;gl=np.array(seed['roots'][3]['return']['derivatives']) if not accepted else np.array(accepted[-1]['lower']['derivatives'])
                sguess=s4-gl[1:]@dq/gl[0]
                snew,rl=lower(qn,sguess)
                ss=np.r_[BASE[:3]+wn[2],snew]
                ww=witnesses(qn,ss)
                oo,stationary=outer(qn,ss[2])
                state={'step':len(accepted)+1,'lambda':target,'q':qn.tolist(),
                       'coefficients':coefficients(qn).tolist(),'w':wn.tolist(),'roots':ss.tolist(),
                       'upper':rn,'lower':rl,'witnesses':ww,'outer':oo,'stationary':stationary,
                       'corrector_iterations':iters,'jacobian_condition':float(np.linalg.cond(Jn)),
                       'calls':calls}
                accepted.append(state)
                print(kind,'accepted',state['step'],'lambda',target,'q',qn,'roots',ss,
                      'outer Ds',[r.get('derivatives',[None])[0] for r in oo], 'stationary',len(stationary),'calls',calls,flush=True)
                save(f'{prefix}_checkpoint.json',{'accepted':accepted,'events':events,'calls':calls})
                lam,w,s4,q=target,wn,snew,qn
                tangent=-np.linalg.solve(Jn,gn[:,1:]@qnlam)
                failed_attempts=0
                if any(r['status']=='unresolved' for r in oo):
                    status='outer return-domain gap; requires separate resolution';break
                if stationary:
                    status='outer stationary point found; follow-up required';break
                # Grow only back to the declared step size after retries.
                step=min(step*1.25,.04 if kind=='shape' else .5)
            except (Gate,np.linalg.LinAlgError) as e:
                events.append({'lambda':target,'step_size':step,'message':str(e),'calls':calls})
                print(kind,'rejected',target,str(e),'calls',calls,flush=True)
                step/=2;failed_attempts+=1
                if failed_attempts>=6 or step<1e-5:
                    status='persistent continuation gate';break
        else:status='attempt cap'
    except Budget:status='500 paired-evaluation budget reached'
    finally:
        ledger.close()
        result={'evidence':'NUM; nonexhaustive two chosen local slices','kind':kind,'status':status,
                'accepted':accepted,'events':events,'calls':calls,'counts':COUNTS.copy(),
                'wall_seconds':time.perf_counter()-started,'last_lambda':lam,'initial_calls':initial_calls,
                'resume':resume,'tail_gap':tail_gap}
        save(f'{prefix}_result.json',result)
        print(kind,'STOP',status,'accepted',len(accepted),'calls',calls,flush=True)

if __name__=='__main__':run(sys.argv[1],'--resume' in sys.argv[2:])
