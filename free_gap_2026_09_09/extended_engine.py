"""Extended finite-field chart: admits -2.2<a<0 for a=-2 boundary test.

Derived from the prior checked engine; only the shape guard is widened.
The ODE formulas and positive time rescaling are unchanged.

Numerical moving-cycle experiment. No interval or nonexistence claims.

q=(a,b,e0/1e-4,e1/1e-4,e2/1e-4); s=log(abs(y)) on x=0.
Smooth positive time rescaling in (v,z)=(asinh(x/abs(y)),log(abs(y))).
Analytic first variations include the initial kick and event-time correction.
"""
from pathlib import Path
import json,time
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

ROOT=Path(__file__).resolve().parent
TAU=1e-4
SEED=np.array([-7/4,1/3,-4/11,-31379/25000,-7517/5000])
COUNTS={'pairs':0,'halves':0,'nfev':0,'failed_halves':0}

def coefficients(q):
    return np.array(q)*np.array([1,1,TAU,TAU,TAU])

def field(vz,q,side,derivatives=False,regularized=False):
    a,b,e0,e1,e2=coefficients(q)
    v,z=vz;u=np.sinh(v);c=np.cosh(v);I=np.exp(-z)
    C=(b-2)/4;B=side*e0*I*I;L=c*c+I*I
    Lv=2*c*u;Lz=-2*I*I
    Le0=0.
    if regularized:
        # Further positive time rescaling prevents finite-q blow-up at y=0.
        # The e0 sign is fixed on each searched branch.
        L+=abs(e0)*I*I*c
        Lv+=abs(e0)*I*I*u;Lz-=2*abs(e0)*I*I*c
        Le0=np.sign(e0)*I*I*c
    N=(a+2)*u*u+b+side*e2*u+side*(1-b)*I+e1*u*I+C*I*I-B*u
    M=(B-2*u)*c
    F=np.array([N,M])/L
    div=(2*(a-1)*u+e1*I+side*e2)*c/L
    physical_time=c*I/L
    if not derivatives:return F,div,physical_time
    Nv=2*(a+2)*u*c+side*e2*c+e1*c*I-B*c
    Nz=-side*(1-b)*I-e1*u*I-2*C*I*I+2*B*u
    Mv=-2*c*c+(B-2*u)*u;Mz=-2*B*c
    J=(np.array([[Nv,Nz],[Mv,Mz]])-F[:,None]*np.array([Lv,Lz]))/L
    K=np.array([[u*u,1-side*I+I*I/4,-side*I*I*u,u*I,side*u],
                [0,0,side*I*I*c,0,0]])
    K[:,2]-=F*Le0
    K/=L
    K[:,2:]*=TAU
    return F,div,physical_time,J,K

def half(s,q,side,sign,rtol=2e-12,variational=True,max_step=.25,regularized=False):
    COUNTS['halves']+=1
    q=np.array(q);a,b=q[:2]
    if not (-2.2<a<0 and 0<b<2):raise ValueError('outside chosen shape chart')
    tc=.5 if side==1 else (2-b)/(2*b)
    if s<=np.log(tc):raise ValueError('outside outer section partition')
    F0=field([0,s],q,side,regularized=regularized)[0]
    if F0[0]<=1e-10:raise ValueError('initial section orientation')
    # State, accumulated signed divergence, signed physical time, variations.
    initial=np.zeros(16 if variational else 4);initial[:2]=[0,s]
    if variational:initial[4:].reshape(2,6)[1,0]=1
    def rhs(t,Y):
        if variational:
            F,di,pt,J,K=field(Y[:2],q,side,True,regularized)
            V=Y[4:].reshape(2,6);dV=J@V;dV[:,1:]+=K
            return sign*np.r_[F,di,pt,dV.ravel()]
        F,di,pt=field(Y[:2],q,side,regularized=regularized)
        return sign*np.r_[F,di,pt]
    def section(t,Y):return Y[0]
    section.terminal=True;section.direction=-sign
    def guard(t,Y):return 50-max(abs(Y[0]),abs(Y[1]))
    guard.terminal=True;guard.direction=-1
    start=time.perf_counter()
    kick=solve_ivp(rhs,[0,1e-7],initial,method='DOP853',rtol=rtol,atol=rtol*.01)
    sol=solve_ivp(rhs,[1e-7,1500],kick.y[:,-1],method='DOP853',rtol=rtol,atol=rtol*.01,
                  events=[section,guard],max_step=max_step)
    COUNTS['nfev']+=kick.nfev+sol.nfev
    if not kick.success or not sol.success or not len(sol.t_events[0]):
        COUNTS['failed_halves']+=1
        return {'status':'unresolved','message':sol.message,'guard':len(sol.t_events[1]),'last':sol.y[:2,-1].tolist()}
    Y=sol.y_events[0][0];F=field(Y[:2],q,side,regularized=regularized)[0]
    if Y[1]>=np.log(tc) or F[0]>=-1e-10:
        COUNTS['failed_halves']+=1
        return {'status':'wrong_section','last':Y[:2].tolist(),'F':F.tolist()}
    out={'status':'passed','z':float(Y[1]),'div':float(Y[2]),'physical_time':float(Y[3]),
         'endpoint_transverse':float(F[0]),'qtime':float(sol.t[-1]),'nfev':kick.nfev+sol.nfev,
         'seconds':time.perf_counter()-start}
    if variational:
        V=Y[4:].reshape(2,6)
        out['derivatives']=(V[1]-F[1]/F[0]*V[0]).tolist()
    return out

def pair(s,q,side=1,rtol=2e-12,variational=True,regularized=False):
    COUNTS['pairs']+=1
    f=half(s,q,side,1,rtol,variational,regularized=regularized);b=half(s,q,side,-1,rtol,variational,regularized=regularized)
    out={'s':float(s),'side':side,'q':np.array(q).tolist(),'forward':f,'backward':b,'rtol':rtol,'regularized':regularized}
    if f['status']==b['status']=='passed':
        out.update(status='passed',D=f['z']-b['z'],log_mu=f['div']-b['div'],period=f['physical_time']-b['physical_time'])
        if variational:out['derivatives']=(np.array(f['derivatives'])-b['derivatives']).tolist()
    else:out['status']='unresolved'
    return out

def cartesian_half(y,q,sign,rtol=2e-13):
    a,b,e0,e1,e2=coefficients(q)
    def rhs(t,Y):
        x,y=Y
        return sign*np.array([(b-2)/4+(1-b)*y+a*x*x+b*y*y+e1*x+e2*x*y,e0-2*x*y])
    def event(t,Y):return Y[0]
    event.terminal=True;event.direction=-sign
    init=np.array([0.,y]);dt=min(1e-7,1e-6/max(1,np.linalg.norm(rhs(0,init))))
    kick=solve_ivp(rhs,[0,dt],init,method='DOP853',rtol=rtol,atol=rtol*.01)
    sol=solve_ivp(rhs,[dt,100],kick.y[:,-1],method='DOP853',events=event,rtol=rtol,atol=rtol*.01,max_step=.1)
    if not kick.success or not sol.success or not len(sol.t_events[0]):raise RuntimeError('Cartesian unresolved')
    return float(sol.y_events[0][0][1])

def cartesian_full(s,q,side,rtol=2e-13):
    """Independent physical-time first return and divergence integral."""
    a,b,e0,e1,e2=coefficients(q);y0=side*np.exp(s)
    def rhs(t,Y):
        x,y=Y[:2]
        return [(b-2)/4+(1-b)*y+a*x*x+b*y*y+e1*x+e2*x*y,
                e0-2*x*y,2*(a-1)*x+e1+e2*y]
    def event(t,Y):return Y[0]
    event.terminal=False;event.direction=0
    init=[0.,y0,0.];dt=min(1e-7,1e-6/max(1,np.linalg.norm(rhs(0,init)[:2])))
    kick=solve_ivp(rhs,[0,dt],init,method='DOP853',rtol=rtol,atol=rtol*.01)
    def finish(t,Y):return Y[0]
    finish.terminal=True;finish.direction=1
    sol=solve_ivp(rhs,[dt,100],kick.y[:,-1],method='DOP853',rtol=rtol,atol=rtol*.01,
                  events=[event,finish],max_step=.1)
    if not sol.success or not len(sol.t_events[1]):raise RuntimeError('full return unresolved')
    Y=sol.y_events[1][0]
    return {'log_return_difference':float(np.log(abs(Y[1]))-s),'log_mu':float(Y[2]),
            'time':float(sol.t_events[1][0]),'crossings':sol.y_events[0].tolist(),'nfev':kick.nfev+sol.nfev}

def save(name,data):
    (ROOT/name).write_text(json.dumps(data,indent=2,allow_nan=False)+'\n')
