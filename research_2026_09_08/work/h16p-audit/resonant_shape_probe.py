"""Explore the a=-1 shapes omitted by inherited nonlogarithmic moment code.
Numerical finite-grid test only; not a cyclicity bound or cycle certificate.
"""
from pathlib import Path
import importlib.util,json,time,hashlib
import numpy as np
from scipy.optimize import brentq
from numpy.polynomial.legendre import leggauss
HERE=Path(__file__).resolve().parent
source=HERE.parent/'H16P/reversible_reseed/moment_search.py'
spec=importlib.util.spec_from_file_location('inherited_arrangement',source)
mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod)

def profile(b,side,n=61,order=256):
    yc=.5 if side==1 else (2-b)/(2*b)
    zc=np.log(yc)
    def V(z):return b*np.exp(z)+side*(1-b)*z+(2-b)/4*np.exp(-z)
    hc=V(zc)
    energies=hc+np.exp(np.linspace(-7,9,n))
    nodes,weights=leggauss(order)
    theta=np.pi*(nodes+1)/2
    rows=[]
    for h in energies:
        left,right=zc-.5,zc+.5
        while V(left)<h:left-=.5
        while V(right)<h:right+=.5
        zl=brentq(lambda z:V(z)-h,left,zc,xtol=1e-13)
        zr=brentq(lambda z:V(z)-h,zc,right,xtol=1e-13)
        z=(zl+zr)/2+(zr-zl)*np.cos(theta)/2
        gap=h-V(z)
        if np.min(gap)<=0:raise ArithmeticError('nonpositive width')
        # x²=exp(z)*(h-V(z)); J_-2 weight after dy=exp(z)dz.
        logw=.5*np.log(gap)-.5*z+np.log(np.sin(theta))
        wts=weights*np.exp(logw-np.max(logw))
        rows.append([h,np.dot(wts,np.exp(z))/sum(wts),np.dot(wts,np.exp(-z))/sum(wts),zl,zr])
    return np.array(rows)

def main():
    import sympy as sp
    x,t,b,side=sp.symbols('x t b side',real=True)
    P=(b-2)/4+side*(1-b)*t-x*x+b*t*t
    H=x*x/t+b*t+side*(1-b)*sp.log(t)+(2-b)/(4*t)
    assert sp.simplify(sp.diff(H,x)*P+sp.diff(H,t)*(-2*x*t))==0
    begin=time.perf_counter();rows=[]
    for b in [.2,.5,.8,1.,1.2,1.5,1.8]:
        upper,lower=profile(b,1),profile(b,-1)
        checku,checkl=profile(b,1,order=512),profile(b,-1,order=512)
        err=max(np.max(np.abs((upper[:,1:3]-checku[:,1:3])/(1+abs(checku[:,1:3])))),np.max(np.abs((lower[:,1:3]-checkl[:,1:3])/(1+abs(checkl[:,1:3])))))
        result=mod.arrangement(checku,checkl)
        row={'b':b,'a':-1,'relative_order_change':float(err),**result,'upper':checku.tolist(),'lower':checkl.tolist()}
        if b==1:
            row['exact_control_max_relative_error']=float(max(np.max(abs(checku[:,1]-.5)),np.max(abs(checku[:,2]-(checku[:,0]+1))/(checku[:,0]+1))))
            assert row['exact_control_max_relative_error']<1e-9
        rows.append(row)
        print(json.dumps({k:v for k,v in row.items() if k not in ('upper','lower')}),flush=True)
        if result['five_candidate']:break
    out={'evidence':'NUM finite shapes and sampled moment curves only; no endpoint or higher-order exclusion','formula':'H=x²/abs(y)+b abs(y)+side*(1-b)*log(abs(y))+(2-b)/(4 abs(y)); normalized M=side*(w-m*v)-c','orders':[256,512],'n_per_annulus':61,'log_energy_excess_range':[-7,9],'arrangement_source_sha256':hashlib.sha256(source.read_bytes()).hexdigest(),'records':rows,'wall_seconds':time.perf_counter()-begin}
    (HERE/'resonant_shape_probe.json').write_text(json.dumps(out,indent=2)+'\n')
if __name__=='__main__':main()
