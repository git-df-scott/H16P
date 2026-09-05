"""KKL log-polar atlas with physical-time and divergence integrals.

Numerical discovery ONLY. q=1/(1+abs(x_section)) is the bounded section
coordinate. Integration uses w=log(hypot(x,y)), theta and
dt/dtau=1/(1+exp(w)); no old Cartesian 2**20 section cutoff is applied.
Finite log-radius/time guards return UNRESOLVED, never a nonexistence claim.
"""
import os
for key in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS'):
    os.environ[key]='1'
import resource
resource.setrlimit(resource.RLIMIT_CPU,(10,10))
import json, math, sys, time
from fractions import Fraction
import numpy as np
from scipy.integrate import solve_ivp
from scipy.special import expit

def evaluate(req):
    started=time.process_time()
    c,alpha,beta=[float(Fraction(str(req.get(k,0)))) for k in ('c','alpha','beta')]
    r=float(req['r']); tol=float(req.get('tol',2e-11))
    section=req.get('section','nullcline')
    transport=req.get('transport',False)
    y0=0. if section=='horizontal' or transport else -r*r/(1+r)
    w0=math.log(math.hypot(r,y0)); theta0=math.atan2(y0,r)
    remote=r < -1
    if not (r>0 or remote): raise ValueError('invalid section branch')
    def components(w,theta):
        C,S=math.cos(theta),math.sin(theta)
        return C,S,C*C+C*S,-10*C*C+2.2*C*S+c*S*S
    def fun(t,z):
        w,theta=z[:2]; C,S,p2,q2=components(w,theta)
        inv=expit(-w); scaled=expit(w)
        p=inv*S+scaled*p2
        q=inv*(alpha*C+beta*S)+scaled*q2
        return [C*p+S*q,C*q-S*p,
                scaled*(4.2*C+(1+2*c)*S)+beta*inv,inv]
    def event_value(t,z):
        w,theta=z[:2]; C,S,p2,q2=components(w,theta)
        if section=='horizontal' and not transport:return S
        return p2+S*math.exp(-w)
    desired_dir=(-1 if not remote else 1) if section=='nullcline' or transport else -1
    def desired(t,z):return event_value(t,z)
    def opposite(t,z):return event_value(t,z)
    desired.direction=desired_dir; opposite.direction=-desired_dir
    desired.terminal=opposite.terminal=True
    def radius_guard(t,z):return float(req.get('log_radius_cap',32))-z[0]
    radius_guard.terminal=True; radius_guard.direction=-1
    def inner_guard(t,z):return z[0]+40
    inner_guard.terminal=True; inner_guard.direction=-1
    def time_guard(t,z):return float(req.get('physical_horizon',10))-z[3]
    time_guard.terminal=True; time_guard.direction=-1
    events=[desired] if transport and not remote else [opposite,desired]
    now=0.; state=np.array([w0,theta0,0.,0.]); segments=[]
    for event in events:
        sol=solve_ivp(fun,[now,20000],state,method='DOP853',
                      rtol=tol,atol=tol*.01,max_step=float(req.get('max_step',.2)),
                      events=[event,radius_guard,inner_guard,time_guard])
        segments.append(sol)
        if not sol.success or len(sol.t_events[0])!=1:
            raise RuntimeError('no full return: '+sol.message+'; guard/time reached')
        now=float(sol.t[-1]);state=sol.y[:,-1]
    w,theta,L,T=state
    R=math.exp(w)*math.cos(theta)
    if (remote and R>=-1) or (not remote and R<=0):
        raise RuntimeError('changed section branch')
    def log_flux(x,y):
        h=math.hypot(x,y); C=x/h;S=y/h
        Qn=-10*C*C+2.2*C*S+c*S*S+(alpha*C+beta*S)/h
        return 2*math.log(h)+math.log(abs(Qn))
    yR=0. if section=='horizontal' and not transport else -R*R/(1+R)
    derivative=math.exp(log_flux(r,y0)-log_flux(R,yR)+L)
    q0=1/(1+abs(r)); qR=1/(1+abs(R))
    q_derivative=derivative*(qR/q0)**2
    log_disp=math.log(abs(R/r))
    path=np.concatenate([s.y for s in segments],axis=1)
    X=np.exp(path[0])*np.cos(path[1]);Y=np.exp(path[0])*np.sin(path[1])
    center=(0.,0.)
    if remote:
        roots=np.roots([c-12.2,alpha-22.2-beta,2*alpha-10-beta,alpha])
        roots=[z.real for z in roots if abs(z.imag)<1e-8 and z.real<-1]
        if len(roots)!=1:raise RuntimeError('remote equilibrium ambiguous')
        ex=roots[0];center=(ex,-ex*ex/(1+ex))
    angles=np.unwrap(np.arctan2(Y-center[1],X-center[0]))
    winding=float((angles[-1]-angles[0])/(2*np.pi))
    if not transport and abs(abs(winding)-1)>.15:
        raise RuntimeError('incorrect winding '+str(winding))
    if remote and max(X)>=-1:raise RuntimeError('crossed x=-1')
    return dict(status='NUMERICAL_ONLY',r=r,c=c,alpha=alpha,beta=beta,
                section=section,transport=transport,return_coordinate=R,
                D=R-r,R_r=derivative,D_r=derivative-1,
                log_displacement=log_disp,log_displacement_derivative=derivative*r/R-1,
                q=q0,q_return=qR,D_q=qR-q0,D_q_derivative=q_derivative-1,
                multiplier=math.exp(L),period=float(T),rescaled_period=now,
                winding_about_focus=winding,
                min_xy=[float(min(X)),float(min(Y))],max_xy=[float(max(X)),float(max(Y))],
                max_log_radius=float(max(path[0])),nfev=sum(s.nfev for s in segments),
                tol=tol,cpu_seconds=time.process_time()-started)

if __name__=='__main__':
    try:print(json.dumps(evaluate(json.load(sys.stdin)),allow_nan=False))
    except Exception as exc:print(json.dumps(dict(status='UNRESOLVED',error=str(exc))))
