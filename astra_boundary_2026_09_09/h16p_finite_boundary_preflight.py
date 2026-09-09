"""Bounded finite-field preflight; NUM only, no validated integration.

Reciprocal chart r=1/abs(y), U=x/y, R=log(r), d/ds=r*d/dt.
Integrates exact original quadratic fields through this derived chart.
Each initial point has x=0 and reciprocal height from the independent
moment witness file. Forward/backward first half-returns are compared.
"""
import json,time
from pathlib import Path
import numpy as np
from scipy.integrate import solve_ivp

ROOT=Path(__file__).resolve().parent
source=json.loads((ROOT/'h16p_rational_boundary_control.json').read_text())
witnesses=[r for r in source['results'] if r['dps']==80]

def half(side,r0,delta,A,direction,rtol):
    alpha=A*delta**2;e0=-delta**3/3;e1=-159*delta**3/80;gamma=8*delta**3/27
    def rhs(t,z):
        R,U=z
        if abs(R)>350:raise ValueError('log-coordinate range gate')
        r=np.exp(R)
        return direction*np.array([side*(2*U-e0*r*r-gamma*U*U),
            side/3+2*r/3-side*5*r*r/12+side*alpha*U*U+e1*r*U
            -side*e0*U*r*r-side*gamma*U**3])
    z0=[np.log(r0),0.]
    launch_slope=rhs(0,z0)[1]
    if abs(launch_slope)<1e-8:raise ValueError('nontransverse launch gate')
    departure=solve_ivp(rhs,[0,1e-6],z0,rtol=rtol,atol=rtol/100,method='DOP853')
    if not departure.success:raise ValueError(departure.message)
    def event(t,z):return z[1]
    event.terminal=True;event.direction=-np.sign(launch_slope)
    sol=solve_ivp(rhs,[1e-6,200],departure.y[:,-1],events=event,
                  rtol=rtol,atol=rtol/100,method='DOP853',max_step=.2)
    if not sol.success or not len(sol.t_events[0]):raise ValueError('return unresolved')
    endpoint=sol.y_events[0][0]
    if abs(rhs(sol.t_events[0][0],endpoint)[1])<1e-8:raise ValueError('terminal transversality gate')
    return {'R':float(endpoint[0]),'chart_time':float(sol.t_events[0][0]),
            'steps':len(sol.t),'terminal_U_derivative':float(rhs(sol.t_events[0][0],endpoint)[1])}

if __name__=='__main__':
    begin=time.monotonic();records=[]
    for delta in (.02,.01):
        for A in (-1.,-.5,0.):
            for rtol in (3e-12,3e-13):
                for w in witnesses:
                    row={'delta':delta,'A':A,'side':w['side'],'h':w['h'],'rtol':rtol}
                    try:
                        r0=float(w['r_turning_points'][0])
                        f=half(w['side'],r0,delta,A,1,rtol)
                        b=half(w['side'],r0,delta,A,-1,rtol)
                        row.update(forward=f,backward=b,D=f['R']-b['R'])
                    except (ValueError,RuntimeError,FloatingPointError) as exc:
                        row['failure']=str(exc)
                    records.append(row)
    print(json.dumps({'scope':'NUM preflight only. Event gates are not rigorous domain certificates.',
        'caps':{'fields':6,'tolerances':2,'witnesses_per_field':8,'chart_time_per_half':200},
        'wall_seconds':time.monotonic()-begin,'records':records},indent=2))
