#!/usr/bin/env python3
"""Bounded original-field NUM control for two compact plus one endpoint cycle.

No interval arithmetic. Half-return differences are checked at fixed brackets;
the script deliberately does not turn floating-point signs into a certificate.
"""
import json
from pathlib import Path

import numpy as np
import scipy
from scipy.integrate import solve_ivp


def difference(t, y_start, tol):
    a,b,e0,e1,e2 = -1.,1.-20*t,t*t,-t,2*t+6*t*t
    def field(_, z):
        x,y=z
        return np.array([(b-2)/4+e1*x+(1-b)*y+a*x*x+e2*x*y+b*y*y,
                         e0-2*x*y])
    hits=[]
    stats=[]
    for direction in (1,-1):
        def rhs(time,z):
            return direction*field(time,z)
        def event(time,z):
            return z[0]
        event.direction=-direction
        event.terminal=True
        sol=solve_ivp(rhs,(0,8),[0.,y_start],method='DOP853',
                      rtol=tol,atol=tol/100,max_step=.1,events=event)
        if not sol.success or len(sol.t_events[0])!=1 or sol.t_events[0][0]<.01:
            raise RuntimeError('Missing nontrivial half return')
        hits.append(float(sol.y_events[0][0,1]))
        stats.append({'time':float(sol.t_events[0][0]),'nfev':sol.nfev})
    # Positive transverse output for reflected lower half-plane.
    return (1 if y_start>0 else -1)*(hits[0]-hits[1]),stats


def main():
    rows=[]
    for t in (1e-3,3e-4,1e-4):
        # Compact h-brackets contain the two predicted limiting roots.
        starts=[]
        for h in (1.2,1.8,4.,6.):
            starts.append(('compact_h_'+str(h),(h+np.sqrt(h*h-1))/2))
        # Lower endpoint predicted by s/t -> 1, with y=-1/s.
        for ratio in (.5,1.5):
            starts.append(('lower_s_over_t_'+str(ratio),-1/(ratio*t)))
        for label,y0 in starts:
            coarse, _ = difference(t,y0,2e-11)
            fine, stats = difference(t,y0,2e-13)
            if np.sign(coarse)!=np.sign(fine):
                raise RuntimeError('Precision check changed displacement sign')
            rows.append({'t':t,'label':label,'y_start':float(y0),
                         'D':fine,'D_over_t2':fine/t**2,
                         'precision_change':abs(fine-coarse),'half_returns':stats})
    report={'classification':'NUM only; no interval certificate',
            'family':{'a':'-1','b':'1-20t','epsilon0':'t^2',
                      'epsilon1':'-t','epsilon2':'2t+6t^2'},
            'scipy':scipy.__version__,'rows':rows}
    out=Path(__file__).parent/'data'/'shoot_control.json'
    out.write_text(json.dumps(report,indent=2)+'\n')
    for row in rows:
        print(row['t'],row['label'],'D/t^2=',format(row['D_over_t2'],'.9g'),
              'precision change=',format(row['precision_change'],'.3g'))


if __name__=='__main__':
    main()
