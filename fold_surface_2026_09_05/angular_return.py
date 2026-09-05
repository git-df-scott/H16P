"""Independent fixed-angle origin return with scalar analytic variations.

Section (r,0), clockwise full turn. Valid only while angular speed is strictly
negative along the computed itinerary. Angular-chart failure is unresolved.
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
    start=time.process_time();c0,K0=F(str(req['c'])),F(str(req['K']))
    c,K=float(c0),float(K0);alpha0=-5*(K0+42)/(11*c0-5);alpha=float(alpha0)
    ac=55*(K+42)/(11*c-5)**2;ak=-5/(11*c-5);r=float(req['r']);tol=float(req.get('tol',2e-12))
    def field(s,z):
        w,M,Vc,VK,Iz,Ic,IK,T,L=z
        C,S=math.cos(s),-math.sin(s);E=math.exp(w)
        h1=(1+alpha)*C*S;g1=alpha*C*C-S*S
        p2=C*C+C*S;q2=-10*C*C+2.2*C*S+c*S*S
        h2=C*p2+S*q2;g2=C*q2-S*p2
        H=h1+E*h2;G=g1+E*g2
        if G>=-1e-14:raise ValueError('angular monotonicity/near-saddle chart guard')
        h1c=ac*C*S;g1c=ac*C*C;h1K=ak*C*S;g1K=ak*C*C
        h2c=S**3;g2c=C*S*S
        Hc=h1c+E*h2c;Gc=g1c+E*g2c;HK=h1K;GK=g1K
        N=h2*g1-h1*g2
        Nc=h2c*g1+h2*g1c-h1c*g2-h1*g2c;NK=h2*g1K-h1K*g2
        fw=-E*N/G**2;fww=fw*(1-2*E*g2/G)
        fc=-(Hc*G-H*Gc)/G**2;fk=-(HK*G-H*GK)/G**2
        fwc=-E*Nc/G**2+2*E*N*Gc/G**3
        fwK=-E*NK/G**2+2*E*N*GK/G**3
        return [-H/G,fw,fw*Vc+fc,fw*VK+fk,fww*math.exp(M),
                fww*Vc+fwc,fww*VK+fwK,-1/G,-E*(4.2*C+(1+2*c)*S)/G]
    sol=solve_ivp(field,[0,2*math.pi],[math.log(r),0,0,0,0,0,0,0,0],method='DOP853',
                  rtol=tol,atol=tol*.01,max_step=.02)
    if not sol.success:raise ValueError(sol.message)
    w,M,Vc,VK,Iz,Ic,IK,T,div=sol.y[:,-1];R=math.exp(w);A=math.exp(M)
    # Axis flux formula is independently compared to scalar angular sensitivity.
    Q0=-10*r*r+alpha*r;Qf=-10*R*R+alpha*R
    fluxA=r/R*Q0/Qf*math.exp(div)
    return dict(status='NUMERICAL_ONLY',section='positive horizontal ray, strictly clockwise',
       c=str(c0),K=str(K0),alpha=str(alpha0),r=r,return_coordinate=R,D=R-r,
       L=w-math.log(r),L_z=math.expm1(M),L_zz=float(A*Iz),L_c=float(Vc),L_K=float(VK),
       L_zc=float(A*Ic),L_zK=float(A*IK),R_r=A*R/r,
       first_derivative_discrepancy=A-fluxA,multiplier=math.exp(div),period=float(T),
       winding=-1.,max_log_radius=float(max(sol.y[0])),cpu_seconds=time.process_time()-start,nfev=sol.nfev)

if __name__=='__main__':
    try:print(json.dumps(evaluate(json.load(sys.stdin)),allow_nan=False))
    except Exception as e:print(json.dumps(dict(status='UNRESOLVED',error=str(e))))
