"""Compactified KKL fold Jacobian via first variations of flow and divergence.

Second return derivatives follow by differentiating the exact transverse-flux
identity. This avoids large Cartesian second-flow sensitivities. NUMERICAL ONLY.
"""
import os
for k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS'):os.environ[k]='1'
import resource
resource.setrlimit(resource.RLIMIT_CPU,(10,10))
import math,json,sys,time
from fractions import Fraction as F
import numpy as np
from scipy.integrate import solve_ivp
from scipy.special import expit

def evaluate(req):
    started=time.process_time();c0,K0=F(str(req['c'])),F(str(req['K']))
    c,K=float(c0),float(K0);alpha0=-5*(K0+42)/(11*c0-5);alpha=float(alpha0)
    ac=55*(K+42)/(11*c-5)**2;ak=-5/(11*c-5)
    beta=float(F(str(req.get('beta',0))));r=float(req['r']);tol=float(req.get('tol',2e-12))
    if r<=0 or alpha>=0:raise ValueError('origin section/equilibrium gate')
    def coefficients(w,t):
        C,S=math.cos(t),math.sin(t);a=expit(-w);b=expit(w)
        p1=S;q1=alpha*C+beta*S;p2=C*C+C*S;q2=-10*C*C+2.2*C*S+c*S*S
        p1t=C;q1t=-alpha*S+beta*C;p2t=-2*C*S+C*C-S*S;q2t=20*C*S+2.2*(C*C-S*S)+2*c*S*C
        h1=C*p1+S*q1;h2=C*p2+S*q2;g1=C*q1-S*p1;g2=C*q2-S*p2
        ht1=-S*p1+C*p1t+C*q1+S*q1t;ht2=-S*p2+C*p2t+C*q2+S*q2t
        gt1=-S*q1+C*q1t-C*p1-S*p1t;gt2=-S*q2+C*q2t-C*p2-S*p2t
        d2=4.2*C+(1+2*c)*S;d2t=-4.2*S+(1+2*c)*C
        f=np.array([a*h1+b*h2,a*g1+b*g2,b*d2+a*beta])
        J=np.array([[a*b*(h2-h1),a*ht1+b*ht2,0.],
                    [a*b*(g2-g1),a*gt1+b*gt2,0.],
                    [a*b*(d2-beta),b*d2t,0.]])
        par=np.array([[0.,a*S*ac*C+b*S**3,a*S*ak*C],
                      [0.,a*C*ac*C+b*C*S*S,a*C*ak*C],
                      [0.,2*b*S,0.]])
        return f,J,par,a
    def fun(t,z):
        f,J,par,a=coefficients(*z[:2]);V=z[4:].reshape(3,3)
        return np.r_[f,a,(J@V+par).ravel()]
    y=-r*r/(1+r);yz=-r*r*(r+2)/(1+r)**2;norm=r*r+y*y
    V=np.zeros((3,3));V[0,0]=(r*r+y*yz)/norm;V[1,0]=(r*yz-y*r)/norm
    state=np.r_[.5*math.log(norm),math.atan2(y,r),0.,0.,V.ravel()]
    def event(t,z):
        w,theta=z[:2];C,S=math.cos(theta),math.sin(theta)
        return C*C+C*S+S*math.exp(-w)
    def opposite(t,z):return event(t,z)
    def desired(t,z):return event(t,z)
    opposite.direction=1;desired.direction=-1;opposite.terminal=desired.terminal=True
    def outer(t,z):return float(req.get('log_radius_cap',100))-z[0]
    outer.direction=-1;outer.terminal=True
    def physical(t,z):return 10-z[3]
    physical.direction=-1;physical.terminal=True
    now=0.;segments=[]
    for e in (opposite,desired):
        sol=solve_ivp(fun,[now,30000],state,rtol=tol,atol=tol*.01,method='DOP853',
                      events=[e,outer,physical],max_step=.2)
        if not sol.success or len(sol.t_events[0])!=1:raise RuntimeError('no complete return (guard/horizon/integration)')
        segments.append(sol);now=float(sol.t[-1]);state=sol.y[:,-1]
    w,t,L,T=state[:4];C,S=math.cos(t),math.sin(t);R=math.exp(w)*C
    if R<=0:raise ValueError('wrong section branch')
    f,_,_,_=coefficients(w,t);V=state[4:].reshape(3,3)
    eg=np.array([-S*math.exp(-w),-2*C*S+C*C-S*S+C*math.exp(-w),0.])
    times=-(eg@V)/(eg@f);Ve=V+f[:,None]*times[None,:]
    logR_der=np.array([1.,-math.tan(t),0.])@Ve
    def flux(x):
        yy=-x*x/(1+x);yp=-x*(x+2)/(1+x)**2
        q=-10*x*x+2.2*x*yy+c*yy*yy+alpha*x+beta*yy
        qp=-20*x+2.2*yy+alpha+(2.2*x+2*c*yy+beta)*yp
        if q>=0:raise ValueError('wrong transverse flux')
        return q,np.array([x*qp/q,(yy*yy+ac*x)/q,ak*x/q])
    q0,g0=flux(r);qf,gf=flux(R)
    LL=math.log(R/r);logA=L+math.log(q0/qf)-LL;A=math.exp(logA)
    gfinal=gf[0]*logR_der+np.array([0.,gf[1],gf[2]])
    logA_der=np.array([1.,0.,0.])-logR_der+g0-gfinal+Ve[2]
    # g0[0] is derivative w.r.t log(initial section radius); other entries explicit c,K.
    F2=math.expm1(logA);F2_der=A*logA_der
    path=np.concatenate([s.y[:2] for s in segments],axis=1)
    X=np.exp(path[0])*np.cos(path[1]);Y=np.exp(path[0])*np.sin(path[1])
    winding=float((path[1,-1]-path[1,0])/(2*np.pi))
    if abs(winding+1)>.05 or min(X)<=-1:raise ValueError('origin full-turn/barrier gate')
    Rr=A*R/r;Rrr=R/(r*r)*(F2_der[0]-A+A*A)
    return dict(status='NUMERICAL_ONLY',c=str(c0),K=str(K0),alpha=str(alpha0),beta=str(req.get('beta',0)),r=r,
      return_coordinate=R,D=R-r,R_r=Rr,R_rr=float(Rrr),L=LL,L_z=F2,L_zz=float(F2_der[0]),
      L_c=float(logR_der[1]),L_K=float(logR_der[2]),L_zc=float(F2_der[1]),L_zK=float(F2_der[2]),
      multiplier=math.exp(L),period=float(T),rescaled_period=now,winding=winding,
      first_derivative_discrepancy=float(logR_der[0]-A),
      min_xy=[float(min(X)),float(min(Y))],max_xy=[float(max(X)),float(max(Y))],
      nfev=sum(s.nfev for s in segments),cpu_seconds=time.process_time()-started,tol=tol)

if __name__=='__main__':
    try:print(json.dumps(evaluate(json.load(sys.stdin)),allow_nan=False))
    except Exception as e:print(json.dumps(dict(status='UNRESOLVED',error=str(e))))
