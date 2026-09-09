"""Nonvalidated directional variational flow and endpoint log-term test."""
from pathlib import Path
import json
import numpy as np
import mpmath as mp
import sympy as S
from scipy.integrate import solve_ivp
from q4_finite_endpoint_probe import DIRECTION
ROOT=Path(__file__).resolve().parent

def half(radius,sign,rtol):
    tau,u,v,w=map(float,DIRECTION)
    def rhs(t,z):
        x,y,sx,sy=z
        p=-y-3*x*x+2*x*y+y*y;q=x+x*x-4*x*y-y*y
        dp=tau*x+(2*u+w)*x*y;dq=tau*y+u*(x*x-y*y)+v*x*y
        return sign*np.array([p,q,(-6*x+2*y)*sx+(-1+2*x+2*y)*sy+dp,(1+2*x-4*y)*sx+(-4*x-2*y)*sy+dq])
    z0=np.array([0.,-radius,0.,0.])
    def section(t,z):return z[0]
    section.terminal=True;section.direction=-sign
    def guard(t,z):return 1e11-np.linalg.norm(z[:2])
    guard.terminal=True;guard.direction=-1
    dt=min(1e-6,1e-6/max(1.,np.linalg.norm(rhs(0,z0)[:2])))
    kick=solve_ivp(rhs,[0,dt],z0,method='DOP853',rtol=rtol,atol=rtol*.002)
    sol=solve_ivp(rhs,[dt,100],kick.y[:,-1],method='DOP853',rtol=rtol,atol=rtol*.002,events=[section,guard],max_step=.15)
    if not kick.success or not sol.success or not len(sol.t_events[0]):return {'status':'unresolved'}
    x,y,sx,sy=sol.y_events[0][0]
    if not 0<y<1:return {'status':'wrong section','y':float(y)}
    p=-y-3*x*x+2*x*y+y*y;q=x+x*x-4*x*y-y*y
    return {'status':'passed','y':float(y),'event_derivative':float(sy-q/p*sx),'nfev':sol.nfev+kick.nfev}

def main():
    z=S.symbols('z',positive=True)
    h=(-2-6*z-6*z*z-z**3)/(2+4*z+z*z)**S.Rational(3,2)
    series=S.series(-1/S.sqrt(2)-h,z,0,4)
    assert S.simplify(S.limit((-1/S.sqrt(2)-h)/z**2,z,0)-3/(4*S.sqrt(2)))==0
    mp.mp.dps=70
    hf=lambda y:(6*y*(1-y)-1+2*y**3)/(1-4*y+2*y*y)**mp.mpf('1.5')
    b=mp.findroot(lambda y:hf(y)+1/mp.sqrt(2),(.15,.17));hy=mp.diff(hf,b)
    entries=json.loads((ROOT/'q4_boundary_integral.json').read_text())['entries']
    direction=[mp.mpf(q.numerator)/q.denominator for q in DIRECTION]
    split=sum(mp.mpf(e['value'])*q for e,q in zip(entries,direction))/hy
    tau,u,v,w=direction
    predicted=3/(2*hy)*(5*u+5*v/4+7*w/4)
    out={'scope':'NUM variational flow; conditional leading-log comparison, no rigorous asymptotic remainder or cycle count',
         'direction':[str(q) for q in DIRECTION],'exact_energy_gap_series_in_inverse_radius':str(series),
         'boundary_y':mp.nstr(b,60),'boundary_hy':mp.nstr(hy,60),'splitting_derivative':mp.nstr(split,60),
         'predicted_log_coefficient':mp.nstr(predicted,60),'runs':[]}
    for rtol in [2e-12,3e-14]:
        rows=[]
        for r in [10,20,40,80,160,320,640,1280,2560]:
            f=half(r,1,rtol);back=half(r,-1,rtol)
            row={'radius':r,'forward':f,'backward':back}
            if f['status']==back['status']=='passed':
                d=f['event_derivative']-back['event_derivative']
                t=r*r*(d-float(split))
                row.update(derivative=d,scaled_residual=t,base_closure=f['y']-back['y'])
                if rows and 'scaled_residual' in rows[-1]:row['doubling_log_slope']=(t-rows[-1]['scaled_residual'])/np.log(2)
            rows.append(row)
        out['runs'].append({'rtol':rtol,'rows':rows})
        print(json.dumps({'rtol':rtol,'prediction':float(predicted),'rows':[{k:v for k,v in row.items() if k not in ['forward','backward']} for row in rows]}),flush=True)
    (ROOT/'q4_variational_endpoint.json').write_text(json.dumps(out,indent=2)+'\n')
if __name__=='__main__':main()
