"""Bounded original-Q4 modulus test, numerical quadrature only.

The modulus values are 1/2, 1 (control), 2. No grid coverage or interval
certificate is claimed. Frozen output preserves all convergence diagnostics.
"""
from pathlib import Path
from functools import lru_cache
import json
import numpy as np
from scipy.optimize import brentq
from numpy.polynomial.legendre import leggauss

HERE=Path(__file__).resolve().parent
@lru_cache(None)
def nodes(n):return leggauss(n)

def field(x,y,r):
    return -y-(2+r*r)*x*x+2*r*x*y+y*y, x+r*x*x-(1+3*r*r)*x*y-r*y*y

def control_form(x,y,p,q):
    dp=np.array([x,2*x*y,0*x,x*y])
    dq=np.array([y,x*x-y*y,x*y,0*x])
    return q*dp-p*dq

def compact(r,t,order):
    d=1+r*r;target=-1+t*(1-1/np.sqrt(d));xx,ww=nodes(order)
    answer=np.zeros(4);maxerr=0.;mins=float('inf')
    for angle,weight in zip(np.pi*(xx+1),np.pi*ww):
        c,s=np.cos(angle),np.sin(angle)
        a=d*(r*c+s)**2;b=-2*d*s
        def residual(rad):
            x,y=rad*c,rad*s;D=1+b*rad+a*rad*rad
            if D<=0:raise ValueError('left positive D component')
            N=3*d*y*(1-r*x-y)-1+d*(r*x+y)**3
            return N/D**1.5-target
        disc=b*b-4*a
        pole=2/(-b+np.sqrt(disc)) if b<0 and disc>0 else None
        upper=.25 if pole is None else min(.25,.9*pole)
        for _ in range(100):
            if residual(upper)>0:break
            upper=upper*2 if pole is None else (upper+pole)/2
        else:raise ValueError('no radius bracket')
        rad=brentq(residual,0,upper,xtol=1e-15,rtol=1e-14)
        maxerr=max(maxerr,abs(residual(rad)))
        x,y=rad*c,rad*s;p,q=field(x,y,r)
        speed=(x*q-y*p)/rad**2
        if speed<=0:raise ValueError('nonpositive angular speed')
        mins=min(mins,speed)
        D=1-2*d*y+d*(r*x+y)**2
        answer+=weight*3*r*r*d*D**(-2.5)*control_form(x,y,p,q)/speed
    return answer,{'energy_residual':maxerr,'minimum_angular_sample':mins}

def boundary(r,order):
    d=1+r*r;xx,ww=nodes(order);end=r/np.sqrt(d);s=end*xx
    g=d*s**3-3*r*r*s-2*r*r;k=d*s*s-r*r
    a=2*(r*r+d*s**3)/(3*k);D=d*g*g/(9*k*k)
    x=(D-d*a*a+r*r)/(2*d*r);y=a+1-r*x
    p,q=field(x,y,r)
    integrand=18*r**3*d/g*(3*k/(np.sqrt(d)*g))**5*control_form(x,y,p,q)
    return -end*(integrand@ww)

def direction(matrix):
    scale=np.linalg.norm(matrix,axis=0)
    _,sv,vh=np.linalg.svd(matrix/scale)
    v=vh[-1]/scale;v/=np.max(np.abs(v))
    if v[-1]<0:v=-v
    return v,sv

def main():
    records=[]
    for r in [.5,1.,2.]:
        matrices=[];diagnostics=[]
        for order in [256,512]:
            rows=[]
            for t in [.25,.5,.75]:
                value,diag=compact(r,t,order);rows.append(value);diagnostics.append({'order':order,'fraction':t,**diag})
            matrices.append(np.array(rows))
        vector,sv=direction(matrices[-1]);old,_=direction(matrices[0])
        b=boundary(r,512);blo=boundary(r,256)
        profile=[]
        for t in [.2,.3,.6,.8]:
            m,diag=compact(r,t,512)
            profile.append({'fraction':t,'directional_value':float(m@vector),**diag})
        records.append({'modulus':r,'direction':vector.tolist(),'scaled_singular_values':sv.tolist(),
                        'matrix':matrices[-1].tolist(),'matrix_order_change_max':float(np.max(np.abs(matrices[1]-matrices[0]))),
                        'direction_order_change_max':float(np.max(np.abs(vector-old))),
                        'boundary_vector':b.tolist(),'boundary_order_change_max':float(np.max(np.abs(b-blo))),
                        'directional_boundary':float(b@vector),'profile':profile,'diagnostics':diagnostics})
        (HERE/'q4_modulus_probe.json').write_text(json.dumps({'status':'numerical only','records':records,'new_ode_returns':0},indent=2)+'\n')
        print(json.dumps({k:records[-1][k] for k in ['modulus','direction','directional_boundary','matrix_order_change_max','direction_order_change_max','boundary_order_change_max','profile']}))

if __name__=="__main__":main()
