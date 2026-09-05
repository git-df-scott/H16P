"""Numerical full-return map with projected sensitivities. NOT a certificate.

One request on stdin, one JSON result on stdout. Every invocation uses one
thread and a ten-CPU-second fuse. A supervising ledger charges failed calls.
"""
import os
import resource
for name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ[name] = "1"
resource.setrlimit(resource.RLIMIT_CPU, (10, 10))

import json
import sys
import time
from fractions import Fraction
import numpy as np
from scipy.integrate import solve_ivp


def number(value):
    return float(Fraction(value)) if isinstance(value, str) else float(value)


def evaluate(request):
    wall, cpu = time.perf_counter(), time.process_time()
    r, c, alpha, beta = (number(request.get(key, 0)) for key in ("r", "c", "alpha", "beta"))
    tol = float(request.get("tol", 1e-11))
    horizon = float(request.get("horizon", 10))
    cap = float(request.get("coordinate_cap", 1e7))
    second_order = bool(request.get('second_order',False))
    determinant_index = 15 if second_order else 9
    def field(t, z):
        x, y = z[:2]
        P = y+x*x+x*y
        Q = -10*x*x+2.2*x*y+c*y*y+alpha*x+beta*y
        J = np.array([[2*x+y, 1+x], [-20*x+2.2*y+alpha, 2.2*x+2*c*y+beta]])
        V = z[2:8].reshape(2, 3)
        dV = J @ V
        dV[1, 1] += y*y
        dV[1, 2] += x
        derivative=np.r_[P,Q,dV.ravel(),4.2*x+(1+2*c)*y+beta]
        if second_order:
            H=z[9:15].reshape(2,3)
            dH=J@H
            vr=V[:,0]
            for j in range(3):
                vj=V[:,j]
                cross=vr[0]*vj[1]+vr[1]*vj[0]
                dH[:,j]+=np.array([2*vr[0]*vj[0]+cross,
                                   -20*vr[0]*vj[0]+2.2*cross+2*c*vr[1]*vj[1]])
            dH[1,1]+=2*y*vr[1]
            dH[1,2]+=vr[0]
            derivative=np.r_[derivative,dH.ravel()]
        div=4.2*x+(1+2*c)*y+beta
        transverse=z[determinant_index:determinant_index+3]
        return np.r_[derivative,div*transverse+np.array([0.,P*y*y,P*x])]

    z0 = np.r_[r,0.,1.,0.,0.,0.,0.,0.,0.]
    if second_order: z0=np.r_[z0,np.zeros(6)]
    q0 = -10*r*r+alpha*r
    z0=np.r_[z0,-q0,0.,0.]
    if q0 >= 0 or r == 0:
        raise ValueError("initial point is not on the downward section")
    if r < 0 and r >= -1:
        raise ValueError("remote section must be left of x=-1")
    def opposite(t,z): return z[1]
    opposite.direction=1
    opposite.terminal=True
    def desired(t,z): return z[1]
    desired.direction=-1
    desired.terminal=True
    def barrier(t,z): return z[0]+1
    barrier.direction=1
    barrier.terminal=True
    def escape(t,z): return cap-max(abs(z[0]),abs(z[1]))
    escape.direction=-1
    escape.terminal=True
    guards=[escape]+([barrier] if r < -1 else [])
    segments=[]
    state=z0
    now=0.
    for event,label in ((opposite,"opposite"),(desired,"downward")):
        sol=solve_ivp(field,[now,horizon],state,method="DOP853",events=[event,*guards],
                      rtol=tol,atol=tol*.01,max_step=.025)
        segments.append(sol)
        if not sol.success:
            raise RuntimeError("integration failure: "+sol.message)
        if len(sol.t_events[0]) != 1:
            raise RuntimeError("missing "+label+" crossing (guard or finite horizon)")
        now=float(sol.t[-1])
        state=sol.y[:,-1]
    R=float(state[0])
    if (r>0 and R<=0) or (r<-1 and R>=-1):
        raise RuntimeError("return changed section component")
    P=R*R
    Q=-10*R*R+alpha*R
    if Q>=0:
        raise RuntimeError("final section is not downward/transverse")
    V=state[2:8].reshape(2,3)
    projected=V[0]-(P/Q)*V[1]
    transverse=-state[determinant_index:determinant_index+3]/Q
    second={}
    if second_order:
        x,y=state[:2]
        J=np.array([[2*x+y,1+x],[-20*x+2.2*y+alpha,2.2*x+2*c*y+beta]])
        parameters=np.array([[0.,0.,0.],[0.,y*y,x]])
        Sdot=J@V+parameters
        times=-V[1]/Q
        acceleration=J@np.array([P,Q])
        H=state[9:15].reshape(2,3)
        for j,key in enumerate(('R_rr','R_rc','R_ralpha')):
            E=H[:,j]+Sdot[:,0]*times[j]+Sdot[:,j]*times[0]+acceleration*times[0]*times[j]
            second[key]=float(E[0]-(P/Q)*E[1])
    speed=float(q0/Q*np.exp(state[8]))
    path=np.concatenate([s.y[:2] for s in segments],axis=1)
    coeff=[c-12.2,alpha-22.2-beta,2*alpha-10-beta,alpha]
    eq_roots=np.roots(coeff)
    reals=[float(z.real) for z in eq_roots if abs(z.imag)<1e-8]
    center=np.array([0.,0.])
    if r<0:
        remote=[z for z in reals if z<-1]
        if len(remote)!=1: raise RuntimeError("remote equilibrium ambiguity")
        ex=remote[0]
        center=np.array([ex,-ex*ex/(1+ex)])
    angle=np.unwrap(np.arctan2(path[1]-center[1],path[0]-center[0]))
    return dict(status="NUMERICAL_ONLY",r=r,c=c,alpha=alpha,beta=beta,**second,
                return_coordinate=R,D=R-r,period=now,
                R_r=float(transverse[0]),R_c=float(transverse[1]),R_alpha=float(transverse[2]),
                speed_derivative=speed,derivative_discrepancy=float(transverse[0]-speed),
                projected_R_r=float(projected[0]),projected_R_c=float(projected[1]),projected_R_alpha=float(projected[2]),
                projection_cancellation_discrepancy=float(projected[0]-speed),
                divergence_exponential=float(np.exp(state[8])),
                opposite_coordinate=float(segments[0].y[0,-1]),
                opposite_time=float(segments[0].t[-1]),
                winding_about_focus=float((angle[-1]-angle[0])/(2*np.pi)),
                min_xy=path.min(axis=1).tolist(),max_xy=path.max(axis=1).tolist(),
                min_distance_focus=float(np.min(np.linalg.norm(path-center[:,None],axis=0))),
                nfev=sum(s.nfev for s in segments),tol=tol,
                wall_seconds=time.perf_counter()-wall,cpu_seconds=time.process_time()-cpu)


if __name__=="__main__":
    request=json.load(sys.stdin)
    try:
        print(json.dumps(evaluate(request),allow_nan=False))
    except Exception as exc:
        print(json.dumps({"status":"UNRESOLVED","error":str(exc)}))
