"""Independent compact Melnikov check using algebraic ovals, without ODEs.

Fixed Gauss-Legendre quadrature and scalar radius root finding are numerical,
not outward-rounded or validated. Compared with saved Cartesian orbit data.
"""
from pathlib import Path
import json
import numpy as np
from scipy.optimize import brentq
from numpy.polynomial.legendre import leggauss

HERE=Path(__file__).resolve().parent

def moments(fraction,order):
    h=-1+fraction*(1-1/np.sqrt(2))
    nodes,weights=leggauss(order)
    total=np.zeros(4);min_angular=float('inf');max_residual=0.
    for angle,weight in zip(np.pi*(nodes+1),np.pi*weights):
        c,s=np.cos(angle),np.sin(angle)
        a=2*(c+s)**2;b=-4*s
        def residual(r):
            x,y=r*c,r*s;d=1+b*r+a*r*r
            if d<=0:raise ValueError('D<=0 in radius solve')
            n=6*y*(1-x-y)-1+2*(x+y)**3
            return n/d**1.5-h
        pole=None
        disc=b*b-4*a
        if b<0 and disc>0:
            pole=2/(-b+np.sqrt(disc))
        upper=.25
        if pole is not None:upper=min(upper,pole*.9)
        for _ in range(100):
            if residual(upper)>0:break
            upper=upper*2 if pole is None else (upper+pole)/2
        else:raise RuntimeError('radius bracket failed')
        r=brentq(residual,0,upper,xtol=1e-15,rtol=1e-14)
        max_residual=max(max_residual,abs(residual(r)))
        x,y=r*c,r*s;d=1-4*y+2*(x+y)**2
        p=-y-3*x*x+2*x*y+y*y;q=x+x*x-4*x*y-y*y
        angular=(x*q-y*p)/(r*r)
        if angular<=0:raise ValueError('nonmonotone angular chart')
        min_angular=min(min_angular,angular)
        dp=np.array([x,2*x*y,0,x*y]);dq=np.array([y,x*x-y*y,x*y,0])
        total+=weight*6*d**(-2.5)*(q*dp-p*dq)/angular
    return {'fraction':fraction,'order':order,'M':total.tolist(),'min_angular_sampled':min_angular,
            'max_energy_root_residual':max_residual}

cartesian={r['fraction']:np.array(r['M']) for r in json.loads((HERE/'q4_original_probe_calls.json').read_text())}
out=[]
for t in [.2,.25,.3,.5,.6,.75,.8]:
    low=moments(t,128);high=moments(t,256)
    high['quadrature_change_max']=float(np.max(np.abs(np.array(high['M'])-low['M'])))
    high['cartesian_difference_max']=float(np.max(np.abs(np.array(high['M'])-cartesian[t])))
    out.append(high)
result={'status':'independent numerical quadrature, not certification','evaluations':out,'new_ode_returns':0}
(HERE/'q4_polar_quadrature.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result,indent=2))
