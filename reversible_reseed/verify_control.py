"""Independent high-precision moment quadrature and original-field controls.

No interval claim: this is a numerical positive control for four cycles.
The exact rational arc is recorded; no fifth cycle is asserted.
"""
import json, time
from pathlib import Path
import mpmath as mp
import numpy as np
from scipy.integrate import solve_ivp

ROOT=Path(__file__).resolve().parent
mp.mp.dps=65
a=-mp.mpf(7)/4;b=mp.mpf(1)/3
m=mp.mpf(7517)/5000;c=mp.mpf(31379)/25000

def R(y,side):
    return -b*y*y/(a+2)-side*(1-b)*y/(a+1)+(2-b)/(4*a)
def potential(y,side): return -R(y,side)*y**a
def bisect(fun,lo,hi):
    fl=fun(lo)
    for _ in range(240):
        mid=(lo+hi)/2;fm=fun(mid)
        if fl*fm<=0: hi=mid
        else: lo=mid;fl=fm
    return (lo+hi)/2
def endpoints(h,side):
    yc=mp.mpf('.5') if side==1 else (2-b)/(2*b)
    fun=lambda y:potential(y,side)-h
    lo=yc/2;hi=yc*2
    while fun(lo)<0: lo/=2
    while fun(hi)<0: hi*=2
    return bisect(fun,lo,yc),bisect(fun,yc,hi)
def ratios(h,side):
    lo,hi=endpoints(h,side)
    def moment(j):
        def f(th):
            y=lo+(hi-lo)*mp.sin(th)**2
            xsq=h*y**(-a)+R(y,side)
            if xsq<0 and abs(xsq)<mp.mpf('1e-50'): xsq=mp.mpf(0)
            return 4*mp.sqrt(xsq)*y**j*(hi-lo)*mp.sin(th)*mp.cos(th)
        return mp.quad(f,[0,mp.pi/16,mp.pi/8,mp.pi/4,3*mp.pi/8,mp.pi/2])
    base=moment(a-1)
    return moment(a)/base,moment(a-2)/base,lo,hi

def half_return(y0,tau,time_sign,rtol):
    af,bf,mf,cf=map(float,(a,b,m,c))
    def rhs(t,z):
        x,y=z
        return time_sign*np.array([(bf-2)/4+(1-bf)*y+af*x*x+bf*y*y-tau*(cf*x+mf*x*y),-2*x*y+tau/(af-1)])
    def event(t,z):return z[0]
    event.terminal=True;event.direction=-time_sign
    initial=np.array([0.,y0])
    dt=min(1e-6,1e-5/max(1.,np.linalg.norm(rhs(0,initial))))
    first=solve_ivp(rhs,[0,dt],initial,method='DOP853',rtol=rtol,atol=rtol*0.01)
    sol=solve_ivp(rhs,[dt,300],first.y[:,-1],method='DOP853',events=event,rtol=rtol,atol=rtol*0.01,max_step=.2)
    if not sol.success or len(sol.t_events[0])!=1:raise RuntimeError(sol.message)
    return float(sol.y_events[0][0][1]),float(sol.t_events[0][0])

def main():
    begin=time.perf_counter()
    output={'evidence':'NUM; non-interval quadrature and original-field half returns',
      'exact_arc':{'a':'-7/4','b':'1/3','epsilon0':'-4*tau/11','epsilon1':'-31379*tau/25000','epsilon2':'-7517*tau/5000'},
      'precision_digits':mp.mp.dps,'moments':[],'shooting':[]}
    for side,hs in [(1,['1','1.5','2','3']),(-1,['10','16'])]:
        for hs0 in hs:
            h=mp.mpf(hs0)
            v,w,lo,hi=ratios(h,side)
            value=side*(w-m*v)-c
            row={'side':side,'h':hs0,'v':mp.nstr(v,55),'w':mp.nstr(w,55),'normalized_M':mp.nstr(value,55),'y_low_abs':mp.nstr(lo,55),'y_high_abs':mp.nstr(hi,55)}
            output['moments'].append(row)
            print(json.dumps(row),flush=True)
    upper=[mp.mpf(r['normalized_M']) for r in output['moments'] if r['side']==1]
    lower=[mp.mpf(r['normalized_M']) for r in output['moments'] if r['side']==-1]
    assert all(upper[i]*upper[i+1]<0 for i in range(3))
    assert lower[0]*lower[1]<0
    for tau in (1e-4,5e-5):
        for rtol in (2e-11,2e-13):
            for r in output['moments']:
                y0=r['side']*float(r['y_high_abs'])
                yf,tf=half_return(y0,tau,1,rtol)
                yb,tb=half_return(y0,tau,-1,rtol)
                row={'tau':tau,'rtol':rtol,'side':r['side'],'h':r['h'],'start_y':y0,'forward_y':yf,'backward_y':yb,'D':yf-yb,'D_over_tau':(yf-yb)/tau,'forward_time':tf,'backward_time':tb}
                output['shooting'].append(row)
                print(json.dumps(row),flush=True)
    for tau in (1e-4,5e-5):
        for rtol in (2e-11,2e-13):
            for side in (1,-1):
                values=[r['D'] for r in output['shooting'] if r['tau']==tau and r['rtol']==rtol and r['side']==side]
                assert all(values[i]*values[i+1]<0 for i in range(len(values)-1))
    output['all_four_sign_brackets_persist']=True
    output['counts']={'moment_values':6,'return_differences':24,'half_flow_integrations':48}
    output['wall_seconds']=time.perf_counter()-begin
    (ROOT/'data'/'verified_control.json').write_text(json.dumps(output,indent=2)+'\n')
    print('seconds',output['wall_seconds'],flush=True)

if __name__=='__main__':main()
