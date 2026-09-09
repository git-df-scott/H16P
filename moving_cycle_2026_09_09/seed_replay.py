import json,time
import numpy as np
from scipy.optimize import brentq
from engine import *

start=time.perf_counter()
old=json.loads((ROOT.parent/'reversible_reseed/data/verified_control.json').read_text())
rows=[r for r in old['shooting'] if r['tau']==1e-4 and r['rtol']==2e-13]
out={'evidence':'NUM, independent coordinates; no interval claim','seed_q':SEED.tolist(),
     'coefficients':coefficients(SEED).tolist(),'witnesses':[],'roots':[],'derivative_checks':[]}
for r in rows:
    y=r['start_y'];side=r['side'];s=np.log(abs(y))
    p=pair(s,SEED,side,rtol=3e-13)
    cf=cartesian_half(y,SEED,1);cb=cartesian_half(y,SEED,-1)
    cartD=cf-cb;cartlog=np.log(abs(cf))-np.log(abs(cb))
    assert p['status']=='passed' and np.sign(cartD)==np.sign(r['D'])
    assert np.sign(p['D'])==np.sign(cartlog)
    out['witnesses'].append(dict(saved=r,logarithmic=p,cartesian=dict(forward=cf,backward=cb,D=cartD,log_D=cartlog)))
    print('witness',side,y,p['D'],cartlog,flush=True)
save('seed_checkpoint.json',out)
for side in [1,-1]:
    ws=[r for r in out['witnesses'] if r['saved']['side']==side]
    for left,right in zip(ws[:-1],ws[1:]):
        lo=left['logarithmic']['s'];hi=right['logarithmic']['s']
        def f(s):
            p=pair(s,SEED,side,rtol=3e-13)
            if p['status']!='passed':raise RuntimeError('root return failed')
            return p['D']
        root=brentq(f,lo,hi,xtol=2e-7)
        p=pair(root,SEED,side,rtol=8e-14)
        cr=cartesian_full(root,SEED,side,rtol=8e-14)
        row={'side':side,'s':root,'height':float(np.exp(root)),'bracket':[lo,hi],
             'return':p,'cartesian_full':cr}
        out['roots'].append(row);save('seed_checkpoint.json',out)
        print('root',side,root,np.exp(root),p['D'],p['derivatives'][0],p['log_mu'],cr,flush=True)

# Smooth RHS partials checked independently by complex-step differentiation.
for state,side in [([.7,-.5],1),([1.2,-1.],-1),([0.,2.],1)]:
    F,di,pt,J,K=field(np.array(state),SEED,side,True)
    jac=np.zeros((2,7));base=np.r_[state,SEED].astype(complex)
    for j in range(7):
        z=base.copy();z[j]+=1e-25j
        jac[:,j]=np.imag(field(z[:2],z[2:],side)[0])/1e-25
    err=float(np.max(abs(jac-np.c_[J,K])));assert err<2e-12
    out['derivative_checks'].append({'kind':'rhs_complex_step','state':state,'side':side,'max_error':err})

# Event derivatives checked against central differences at two step sizes.
for s,side,indices in [(1.2,1,[0,1,3]),(9.,-1,[0,2,4])]:
    p=pair(s,SEED,side,rtol=8e-14)
    for j in indices:
        vals=[]
        for h in [1e-3,3e-4]:
            x=np.r_[s,SEED];xp=x.copy();xm=x.copy();xp[j]+=h;xm[j]-=h
            pp=pair(xp[0],xp[1:],side,rtol=8e-14,variational=False)
            pm=pair(xm[0],xm[1:],side,rtol=8e-14,variational=False)
            fd=(pp['D']-pm['D'])/(2*h)
            vals.append({'h':h,'finite_difference':fd,'error':fd-p['derivatives'][j]})
        out['derivative_checks'].append({'kind':'event_central_difference','s':s,'side':side,
            'index':j,'variational':p['derivatives'][j],'checks':vals})
out['counts']=COUNTS.copy();out['wall_seconds']=time.perf_counter()-start
save('seed_replay.json',out)
print('completed',out['counts'],out['wall_seconds'],flush=True)
