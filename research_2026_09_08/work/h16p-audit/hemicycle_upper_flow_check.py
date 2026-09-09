"""Independent full-trajectory check of the upper weighted moment ratio."""
from pathlib import Path
import json,numpy as np
from scipy.integrate import solve_ivp
P=Path(__file__).resolve().parent
source=json.loads((P/'hemicycle_upper_moment.json').read_text());records=[]
from scipy.optimize import brentq
for a,b in [(-.5,.25),(-1.25,.5),(-1.75,1.)]:
    data=next(q for q in source['records'] if q['a']==a and q['b']==b);k=data['k']
    def V(z):
        y=np.exp(z);return y**a*(b*y*y/(a+2)-(b-1)*y/(a+1)+(b-2)/(4*a))
    for index in [15,30,45]:
        target=data['rows'][index];h=target['energy'];zc=np.log(.5);zr=zc+.5
        while V(zr)<h:zr+=.5
        y0=float(np.exp(brentq(lambda z:V(z)-h,zc,zr)))
        def rhs(t,z):
            x,y=z[:2];weight=2*y**a*x*x
            return [(b-2)/4+(1-b)*y+a*x*x+b*y*y,-2*x*y,weight*(1-2*k*y),weight]
        def event(t,z):return z[0]
        event.direction=1;event.terminal=True
        initial=[0.,y0,0.,0.];kick=solve_ivp(rhs,[0,1e-7],initial,method='DOP853',rtol=3e-12,atol=1e-14)
        sol=solve_ivp(rhs,[1e-7,1e5],kick.y[:,-1],method='DOP853',rtol=3e-12,atol=1e-14,events=event,max_step=.1)
        row=dict(a=a,b=b,index=index,energy=h,initial_y=y0,status='passed' if sol.success and len(sol.t_events[0]) else 'unresolved',message=sol.message)
        if row['status']=='passed':
            end=sol.y_events[0][0];ratio=float(end[2]/end[3]);row.update(ratio=ratio,quadrature_ratio=target['normalized_energy_integral'],difference=abs(ratio-target['normalized_energy_integral']),return_y_error=float(abs(end[1]-y0)),relative_return_y_error=float(abs(end[1]-y0)/y0),period=float(sol.t[-1]))
        records.append(row);print(json.dumps(row),flush=True)
        (P/'hemicycle_upper_flow_check.json').write_text(json.dumps(dict(scope=__doc__,records=records),indent=2)+'\n')
