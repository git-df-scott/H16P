"""Bounded numerical compact-Melnikov probe in original Q4 coordinates.

Exploratory floating-point calculations, not interval cycle certificates.
Three anchors determine one common four-control direction. Each evaluation
is a full base-field orbit with four first-order energy integrals.
"""
from pathlib import Path
import json,time
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

HERE=Path(__file__).resolve().parent
records=[]

def energy(x,y):
    d=1-4*y+2*(x+y)**2
    n=6*y*(1-x-y)-1+2*(x+y)**3
    return n/d**1.5

def evaluate(fraction,rtol=2e-12):
    started=time.monotonic()
    target=-1+fraction*(1-1/np.sqrt(2))
    y0=brentq(lambda y:energy(0,y)-target,1e-9,.15946453,xtol=5e-16)
    def rhs(t,state):
        x,y=state[:2]
        p=-y-3*x*x+2*x*y+y*y
        q=x+x*x-4*x*y-y*y
        d=1-4*y+2*(x+y)**2
        if d<=0:raise ValueError('left D>0 chart')
        dp=np.array([x,2*x*y,0,x*y])
        dq=np.array([y,x*x-y*y,x*y,0])
        return np.r_[p,q,6*d**(-2.5)*(q*dp-p*dq)]
    def crossing(t,state):return state[0]
    crossing.terminal=True
    state=np.r_[0,y0,np.zeros(4)]
    end=0.; pieces=[]
    try:
        for direction in [1,-1]:
            crossing.direction=direction
            sol=solve_ivp(rhs,(end,end+100),state,method='DOP853',events=crossing,
                          rtol=rtol,atol=rtol/50,max_step=.05)
            if not sol.success or not len(sol.t_events[0]):raise ValueError('missing oriented half-return')
            pieces.append(sol.y)
            end=float(sol.t[-1]);state=sol.y[:,-1]
        trace=np.concatenate(pieces,axis=1)
        row={'fraction':fraction,'h':target,'y0':y0,'rtol':rtol,'status':'numerical_return',
             'M':state[2:].tolist(),'return_y_error':float(state[1]-y0),
             'max_h_drift':float(np.max(np.abs(energy(trace[0],trace[1])-target))),
             'period':end,'max_radius':float(np.max(np.hypot(trace[0],trace[1])))}
    except Exception as exc:
        row={'fraction':fraction,'rtol':rtol,'status':'unresolved','error':str(exc)}
    row['wall_seconds']=time.monotonic()-started
    records.append(row)
    (HERE/'q4_original_probe_calls.json').write_text(json.dumps(records,indent=2)+'\n')
    if row['status']!='numerical_return':raise RuntimeError(row)
    return np.array(row['M'])

anchors=[.25,.5,.75]
matrix=np.array([evaluate(t) for t in anchors])
scales=np.linalg.norm(matrix,axis=0)
_,singular,Vh=np.linalg.svd(matrix/scales)
direction=Vh[-1]/scales
direction/=np.max(np.abs(direction))
if direction[-1]<0:direction=-direction
profile=[]
for t in [.1,.2,.3,.4,.6,.7,.8,.9,.97]:
    value=evaluate(t)
    profile.append({'fraction':t,'M_direction':float(value@direction)})
tighter=np.array([evaluate(t,3e-14) for t in anchors])
_,tight_singular,tight_Vh=np.linalg.svd(tighter/scales)
tight_direction=tight_Vh[-1]/scales
tight_direction/=np.max(np.abs(tight_direction))
if tight_direction[-1]<0:tight_direction=-tight_direction
out={'status':'numerical_only','normal_control_order':['tau','u','v','w'],
     'anchors':anchors,'anchor_matrix':matrix.tolist(),'scaled_singular_values':singular.tolist(),
     'direction':direction.tolist(),'profile':profile,
     'tighter_matrix':tighter.tolist(),'tighter_direction':tight_direction.tolist(),
     'direction_change_max':float(np.max(np.abs(tight_direction-direction))),
     'ode_returns_requested':len(records),
     'scope':'No finite perturbed field or endpoint cycles computed. No interval signs, derivative simplicity, or cycle certificates.'}
(HERE/'q4_original_probe.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out,indent=2))
