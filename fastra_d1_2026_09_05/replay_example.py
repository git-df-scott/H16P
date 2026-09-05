"""Reproduce a rounded rational 3+1 field using log-polar and Cartesian returns."""
import json,math
import numpy as np
from scipy.optimize import brentq
from run_d1 import HERE,rm,sl,Q,finite,fable_count,hp,append

v=list(map(str,[0,0,1,1,1,0,0,-Q(37101199745401,10**12),-Q(1,25600000),-10,Q(11,5),Q(969154254284,10**12)]))
c=np.array([float(Q(x)) for x in v]);base=fable_count(v);replays=[]
for n in base['nests']:
    pt=np.array(n['pt']);theta=n['theta'];direction=np.array([math.cos(theta),math.sin(theta)])
    for r,stab in zip(n['roots'],n['stab']):
        def logD(u,tol=1e-13):
            u1,_,s=rm.returns_log(c[None],pt[None],np.array([[u]]),th0=theta,rtol=tol,umax=60,Smax=10000,maxsteps=3000000)
            if s[0,0]:raise ValueError('failed return')
            return float(u1[0,0]-u)
        root=brentq(logD,math.log(r)-.14,math.log(r)+.14,xtol=2e-8)
        checks=[]
        for off in [-.001,.001]:
            rad=math.exp(root+off)
            R,T,status=rm.returns(c[None],pt[None],direction[None],np.array([[rad]]),rtol=1e-12,Rmax=1e25,Tmax=1e6,maxsteps=5000000)
            checks.append(finite(dict(r=rad,log_D=logD(root+off),cartesian_D_over_r=(R[0,0]-rad)/rad,cartesian_status=status[0,0],cartesian_period=T[0,0])))
        replays.append(dict(focus=pt.tolist(),theta=theta,r=math.exp(root),stability=stab,checks=checks))
row=dict(coefficient_vector=v,c=str(Q(v[11])),K=str((-Q(v[7]))*(11*Q(v[11])-5)/5-42),beta=v[8],fable=base,replays=replays,status='NUMERICAL_3_PLUS_1_NOT_A_COUNTEREXAMPLE')
(HERE/'rational_3_plus_1.json').write_text(json.dumps(row,indent=2))
print(json.dumps({k:row[k] for k in ['coefficient_vector','K','beta','replays']}),flush=True)
