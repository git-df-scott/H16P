"""NUM compact-moment probes at a=-2 (repaired chart) and a=0.

These use finite energies only and say nothing about endpoint cyclicity.
"""
import json,time
import numpy as np
from scipy.optimize import brentq
from numpy.polynomial.legendre import leggauss
from moment_search import ROOT,arrangement

def boundary_profile(a,b,side,n=41,order=200):
    if a==-2:
        # u=x/y, r=1/|y|. The moment coordinates are E[u^2/r], E[r].
        rc=2. if side==1 else 2*b/(2-b)
        def potential(z):
            r=np.exp(z)
            return (2-b)*r*r/8-side*(1-b)*r-b*z
        def observables(z,width):return width/(3*np.exp(z)),np.exp(z)
        exponent=1 # base area is du dr
    elif a==0:
        rc=.5 if side==1 else (2-b)/(2*b)
        def potential(z):
            r=np.exp(z)
            return b*r*r/2+side*(1-b)*r+(b-2)*z/4
        def observables(z,width):return np.exp(z),np.exp(-z)
        exponent=0 # base moment is dx dy/|y|
    else: raise ValueError(a)
    zc=np.log(rc);hc=potential(zc)
    energies=hc+np.geomspace(1e-4,10.,n)
    nodes,weights=leggauss(order)
    theta=(nodes+1)*np.pi/2
    output=[]
    for h in energies:
        f=lambda z:potential(z)-h
        zl,zr=zc-.5,zc+.5
        while f(zl)<0:zl-=.5
        while f(zr)<0:zr+=.5
        zl=brentq(f,zl,zc,xtol=1e-13);zr=brentq(f,zc,zr,xtol=1e-13)
        z=(zl+zr)/2+(zr-zl)*np.cos(theta)/2
        width=h-potential(z)
        assert np.min(width)>0
        lw=.5*np.log(width)+exponent*z+np.log(np.sin(theta))
        ww=weights*np.exp(lw-lw.max())
        v,w=observables(z,width)
        output.append([float(h),float(np.dot(ww,v)/ww.sum()),float(np.dot(ww,w)/ww.sum()),float(zl),float(zr)])
    return np.asarray(output)

def main():
    begin=time.perf_counter()
    output={'evidence':'NUM, finite compact energies only',
        'grid':{'points_per_annulus':41,'gauss_order':200,'energy_increment':[1e-4,10.]},
        'a_minus_2_coefficients':{'epsilon1':'-c','epsilon0':'-1/3','gamma':'m/3','epsilon2':'0'},
        'a_zero_coefficients':{'epsilon1':'-c','epsilon0':'-1','epsilon2':'-m'},
        'records':[]}
    for a in (-2.,0.):
        for b in (1/3,2/3,1.,4/3,5/3):
            upper=boundary_profile(a,b,1);lower=boundary_profile(a,b,-1)
            record={'a':a,'b':b,**arrangement(upper,lower),
                    'upper_profile':upper.tolist(),'lower_profile':lower.tolist()}
            output['records'].append(record)
            print(json.dumps({k:v for k,v in record.items() if not k.endswith('_profile')}),flush=True)
            if record['five_candidate']:break
    output['wall_seconds']=time.perf_counter()-begin
    (ROOT/'data'/'boundary_search.json').write_text(json.dumps(output,indent=2)+'\n')

if __name__=='__main__':main()
