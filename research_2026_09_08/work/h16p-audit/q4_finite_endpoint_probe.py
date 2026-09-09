"""Original Q4 two-sided passages, including large negative-y sections.
Numerical profiles only; not an endpoint expansion or interval certificate.
"""
from pathlib import Path
from fractions import Fraction as F
import numpy as np,json,time
from scipy.integrate import solve_ivp
ROOT=Path(__file__).resolve().parent
DIRECTION=[F('0.00041138531923618854'),F('0.42436759177461253'),F(-1),F('0.1029863328970071')]
def half(radius,delta,sign,rtol=2e-12):
    tau,u,v,w=[float(delta*q) for q in DIRECTION]
    def rhs(t,z):
        x,y=z
        return sign*np.array([tau*x-y-3*x*x+(2+2*u+w)*x*y+y*y,x+tau*y+(1+u)*x*x+(-4+v)*x*y-(1+u)*y*y])
    initial=np.array([0.,-radius])
    def section(t,z):return z[0]
    section.terminal=True;section.direction=-sign
    def guard(t,z):return 1e11-np.linalg.norm(z)
    guard.terminal=True;guard.direction=-1
    dt=min(1e-6,1e-6/max(1.,np.linalg.norm(rhs(0,initial))))
    kick=solve_ivp(rhs,[0,dt],initial,method='DOP853',rtol=rtol,atol=rtol*.002)
    if not kick.success:return {'status':'unresolved','reason':'kick'}
    sol=solve_ivp(rhs,[dt,100],kick.y[:,-1],method='DOP853',rtol=rtol,atol=rtol*.002,events=[section,guard],max_step=.15)
    if not sol.success or not len(sol.t_events[0]):return {'status':'unresolved','reason':sol.message,'guard':len(sol.t_events[1])}
    y=float(sol.y_events[0][0][1])
    if not 0<y<1:return {'status':'unresolved','reason':'wrong crossing segment','y':y}
    return {'status':'half_passage','y':y,'time':float(sol.t_events[0][0]),'nfev':kick.nfev+sol.nfev}
def main():
    out={'evidence':'NUM only; sampled two-sided passages; no rigorous root count','direction':[str(x) for x in DIRECTION],'records':[]};start=time.perf_counter()
    for delta in [F(1,1000),F(1,10000),F(-1,1000)]:
        rows=[]
        for radius in np.geomspace(.01,1e8,91):
            forward=half(radius,delta,1);backward=half(radius,delta,-1)
            row={'radius':float(radius),'forward':forward,'backward':backward}
            if forward['status']==backward['status']=='half_passage':row['difference']=forward['y']-backward['y']
            rows.append(row)
        brackets=[];tight=[]
        for p,q in zip(rows,rows[1:]):
            if 'difference' in p and 'difference' in q and p['difference']*q['difference']<0 and min(abs(p['difference']),abs(q['difference']))>1e-11:
                brackets.append([p['radius'],q['radius']])
                for row in [p,q]:
                    f=half(row['radius'],delta,1,2e-14);b=half(row['radius'],delta,-1,2e-14)
                    if f['status']!=b['status'] or f['status']!='half_passage':raise ArithmeticError('tight passage failure')
                    diff=f['y']-b['y'];tight.append({'radius':row['radius'],'difference':diff,'change':diff-row['difference']});assert diff*row['difference']>0
        record={'delta':str(delta),'parameters':[str(delta*x) for x in DIRECTION],'rows':rows,'brackets':brackets,'tight':tight,'unresolved':sum('difference' not in r for r in rows)}
        out['records'].append(record);out['wall_seconds']=time.perf_counter()-start
        (ROOT/'q4_finite_endpoint_probe.json').write_text(json.dumps(out,indent=2)+'\n')
        print(json.dumps({k:v for k,v in record.items() if k not in ['rows','tight']}),flush=True)
if __name__=='__main__':main()
