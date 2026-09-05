"""Independent original-time Cartesian KKL section return, numerical only."""
import os
for k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS'):os.environ[k]='1'
import resource
resource.setrlimit(resource.RLIMIT_CPU,(10,10))
import json,sys,time,math
from fractions import Fraction
import numpy as np
from scipy.integrate import solve_ivp

def evaluate(req):
    started=time.process_time();r=float(req['r'])
    c,a,b=[float(Fraction(str(req.get(k,0)))) for k in ('c','alpha','beta')]
    def field(t,z):
        x,y=z[:2]
        return [y+x*x+x*y,-10*x*x+2.2*x*y+c*y*y+a*x+b*y,
                4.2*x+(1+2*c)*y+b]
    def section(t,z):return z[1]+z[0]*z[0]+z[0]*z[1]
    def opposite(t,z):return section(t,z)
    def desired(t,z):return section(t,z)
    opposite.direction=1 if r>0 else -1;desired.direction=-opposite.direction
    opposite.terminal=desired.terminal=True
    now=0.;state=np.array([r,-r*r/(1+r),0.]);segments=[]
    if field(0,state)[1]>=0:raise ValueError('bad initial flux')
    for event in (opposite,desired):
        sol=solve_ivp(field,[now,10],state,method='DOP853',rtol=float(req.get('tol',2e-13)),
                      atol=2e-15,max_step=.01,events=event)
        if not sol.success or len(sol.t_events[0])!=1:raise RuntimeError('missing crossing')
        segments.append(sol);now=float(sol.t[-1]);state=sol.y[:,-1]
    R=float(state[0]);y=float(state[1]);Q=field(now,state)[1]
    if Q>=0 or (r>0)!=(R>0):raise ValueError('bad return branch/flux')
    Q0=field(0,[r,-r*r/(1+r),0])[1];der=Q0/Q*math.exp(state[2])
    return dict(status='NUMERICAL_ONLY',r=r,c=c,alpha=a,beta=b,
                return_coordinate=R,log_displacement=math.log(abs(R/r)),D=R-r,
                R_r=der,multiplier=math.exp(state[2]),period=now,
                nfev=sum(s.nfev for s in segments),cpu_seconds=time.process_time()-started)
if __name__=='__main__':
    try:print(json.dumps(evaluate(json.load(sys.stdin)),allow_nan=False))
    except Exception as exc:print(json.dumps(dict(status='UNRESOLVED',error=str(exc))))
