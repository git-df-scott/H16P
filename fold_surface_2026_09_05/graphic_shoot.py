"""Direct two-manifold splitting at the exact c=8/5 infinity resonance.

One branch per charged call. The local invariant manifold is initialized by
an order8 series in v=1/x; changing v0 is a convergence test, not an enclosure.
The transverse endpoint is the negative horizontal ray, common to the two
branches. These are separatrix passages, NOT periodic returns.
"""
import os
for k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS'):os.environ[k]='1'
import resource
resource.setrlimit(resource.RLIMIT_CPU,(10,10))
import json,math,sys,time
from fractions import Fraction as F
import numpy as np
from scipy.integrate import solve_ivp
from scipy.special import expit
from numpy.polynomial import polynomial as pp

def saddle_series(m,branch,N=8):
    z0=-1+branch*math.sqrt(159)/3;c=1.6;z=np.zeros(N+1);z[0]=z0
    for n in range(1,N+1):
        zp=pp.polyder(z);one=z.copy();one[0]+=1
        den=-pp.polyadd(np.r_[0.,one],np.r_[0.,0.,z])
        p=pp.polyadd(np.array([-10.]),pp.polyadd(1.2*z,.6*pp.polymul(z,z)))
        forcing=pp.polymul(z,z);forcing[0]+=m
        residual=pp.polyadd(pp.polysub(pp.polymul(zp,den),p),np.r_[0.,forcing])
        z[n]=residual[n]/(n*(1+z0)+1.2*(1+z0))
    return z

def evaluate(req):
    start=time.process_time();K0=F(str(req['K']));m0=(K0+42)/F(63,25)
    K,m=float(K0),float(m0);branch=int(req['branch']);c=1.6
    v0=float(req.get('v0',1e-4));series=saddle_series(m,branch)
    zz=pp.polyval(v0,series);x=1/v0;y=zz/v0;sign=-branch
    def fun(t,z):
        w,th=z[:2];C,S=math.cos(th),math.sin(th);a=expit(-w);b=expit(w)
        P=a*S+b*(C*C+C*S);Q=-a*m*C+b*(-10*C*C+2.2*C*S+c*S*S)
        return sign*np.array([C*P+S*Q,C*Q-S*P,a])
    def ray(t,z):return math.sin(z[1])
    ray.direction=1 if branch==-1 else -1;ray.terminal=True
    state=[math.log(math.hypot(x,y)),math.atan2(y,x),0.]
    sol=solve_ivp(fun,[0,3000],state,rtol=float(req.get('tol',2e-12)),atol=2e-14,
                  max_step=.2,events=ray,method='DOP853')
    if not sol.success or len(sol.t_events[0])!=1:raise RuntimeError('no negative-axis passage')
    end=sol.y[:,-1];X=math.exp(end[0])*math.cos(end[1]);Y=math.exp(end[0])*math.sin(end[1])
    if not -1<X<0:raise RuntimeError('endpoint outside origin negative ray')
    return dict(status='NUMERICAL_SEPARATRIX_PASSAGE',K=str(K0),m=str(m0),c='8/5',branch=branch,
                v0=v0,series=series.tolist(),endpoint_x=X,endpoint_y=Y,physical_time=float(end[2]),
                rescaled_time=float(sol.t[-1]),nfev=sol.nfev,cpu_seconds=time.process_time()-start)

if __name__=='__main__':
    try:print(json.dumps(evaluate(json.load(sys.stdin)),allow_nan=False))
    except Exception as e:print(json.dumps(dict(status='UNRESOLVED',error=str(e))))
