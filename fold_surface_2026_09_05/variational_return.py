"""Curved-section KKL return and second derivatives in original time.

Discovery arithmetic only. Analytic first/second flow variations and moving
event corrections give the augmented-fold Jacobian. The radial first
derivative is independently checked against the divergence/flux identity.
"""
import os
for k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS'):os.environ[k]='1'
import resource
resource.setrlimit(resource.RLIMIT_CPU,(10,10))
import json,sys,math,time
from fractions import Fraction as F
import numpy as np
from scipy.integrate import solve_ivp

def evaluate(req):
    start=time.process_time()
    c0,K0=F(str(req['c'])),F(str(req['K']))
    alpha0=-5*(K0+42)/(11*c0-5)
    c,K,alpha=float(c0),float(K0),float(alpha0)
    r=float(req['r']);beta=float(F(str(req.get('beta',0))))
    tol=float(req.get('tol',2e-12))
    if not r>0 or alpha>=0:raise ValueError('origin focus/section gate')
    ac=55*(K+42)/(11*c-5)**2;ak=-5/(11*c-5)
    def field(t,z):
        x,y=z[:2];P=y+x*x+x*y;Q=-10*x*x+2.2*x*y+c*y*y+alpha*x+beta*y
        J=np.array([[2*x+y,1+x],[-20*x+2.2*y+alpha,2.2*x+2*c*y+beta]])
        V=z[2:8].reshape(2,3);H=z[9:15].reshape(2,3)
        par=np.array([[0.,0.,0.],[0.,y*y+ac*x,ak*x]])
        dV=J@V+par;dH=J@H;vr=V[:,0]
        for j in range(3):
            vj=V[:,j];cross=vr[0]*vj[1]+vr[1]*vj[0]
            dH[:,j]+=np.array([2*vr[0]*vj[0]+cross,-20*vr[0]*vj[0]+2.2*cross+2*c*vr[1]*vj[1]])
        dH[1,1]+=ac*vr[0]+2*y*vr[1]
        dH[1,2]+=ak*vr[0]
        return np.r_[P,Q,dV.ravel(),4.2*x+(1+2*c)*y+beta,dH.ravel()]
    y=-r*r/(1+r)
    V=np.array([[1.,0.,0.],[-r*(r+2)/(1+r)**2,0.,0.]])
    H=np.zeros((2,3));H[1,0]=-2/(1+r)**3
    state=np.r_[r,y,V.ravel(),0.,H.ravel()]
    Q0=field(0,state)[1]
    if Q0>=0:raise ValueError('initial section flux')
    def desired(t,z):return z[1]+z[0]*z[0]+z[0]*z[1]
    def opposite(t,z):return desired(t,z)
    desired.direction=-1;opposite.direction=1
    desired.terminal=opposite.terminal=True
    def cap(t,z):return 1e12-max(abs(z[0]),abs(z[1]))
    cap.direction=-1;cap.terminal=True
    segments=[];now=0.
    for event in (opposite,desired):
        sol=solve_ivp(field,[now,10],state,rtol=tol,atol=tol*.01,method='DOP853',
                      events=[event,cap],max_step=.02)
        if not sol.success or len(sol.t_events[0])!=1:raise RuntimeError('no full return (guard/integration/horizon)')
        segments.append(sol);state=sol.y[:,-1];now=float(sol.t[-1])
    R,Y=state[:2];f=field(now,state);Q=f[1]
    if R<=0 or Q>=0:raise ValueError('return branch/flux')
    J=np.array([[2*R+Y,1+R],[-20*R+2.2*Y+alpha,2.2*R+2*c*Y+beta]])
    V=state[2:8].reshape(2,3);H=state[9:15].reshape(2,3)
    grad=J[0];den=grad@f[:2]
    times=-(grad@V)/den
    params=np.array([[0.,0.,0.],[0.,Y*Y+ac*R,ak*R]])
    Vdot=J@V+params;acc=J@f[:2]
    second=np.array([H[0,j]+Vdot[0,0]*times[j]+Vdot[0,j]*times[0]+acc[0]*times[0]*times[j] for j in range(3)])
    logRr=math.log(Q0/Q)+state[8];Rr=math.exp(logRr)
    # Tiny fixed-time residual P at the numerical event is corrected for first variations.
    first=V[0]+f[0]*times
    Rc,RK=first[1:]
    L=math.log(R/r);A=Rr*r/R
    Lz=math.expm1(logRr-L)
    Lzz=A+r*r*second[0]/R-A*A
    Lc=Rc/R;LK=RK/R
    Lzc=r*second[1]/R-A*Rc/R;LzK=r*second[2]/R-A*RK/R
    path=np.concatenate([s.y[:2] for s in segments],axis=1)
    angle=np.unwrap(np.arctan2(path[1],path[0]));winding=float((angle[-1]-angle[0])/(2*np.pi))
    if abs(winding+1)>.05 or min(path[0])<=-1:raise ValueError('itinerary changed')
    return dict(status='NUMERICAL_ONLY',c=str(c0),K=str(K0),alpha=str(alpha0),beta=str(req.get('beta',0)),r=r,
      return_coordinate=float(R),D=float(R-r),R_r=Rr,R_c=float(Rc),R_K=float(RK),
      R_rr=float(second[0]),R_rc=float(second[1]),R_rK=float(second[2]),
      L=L,L_z=Lz,L_zz=float(Lzz),L_c=float(Lc),L_K=float(LK),L_zc=float(Lzc),L_zK=float(LzK),
      multiplier=math.exp(state[8]),first_derivative_discrepancy=float(first[0]-Rr),
      period=now,winding=winding,min_xy=path.min(axis=1).tolist(),max_xy=path.max(axis=1).tolist(),
      cpu_seconds=time.process_time()-start,nfev=sum(s.nfev for s in segments),tol=tol)

if __name__=='__main__':
    try:print(json.dumps(evaluate(json.load(sys.stdin)),allow_nan=False))
    except Exception as e:print(json.dumps(dict(status='UNRESOLVED',error=str(e))))
