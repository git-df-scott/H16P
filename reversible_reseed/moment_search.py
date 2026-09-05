"""Bounded NUM search of the SAME perturbation on both reversible annuli.

The polyline search is exhaustive only for its finite sampled moment curves.
No absence claim between energy samples or between shapes is justified.
"""
import argparse, json, time
from pathlib import Path
import numpy as np
import scipy
from scipy.optimize import brentq
from numpy.polynomial.legendre import leggauss

ROOT=Path(__file__).resolve().parent

def profile(a,b,side=1,n=41,order=160):
    if a in (-1,-2,0):
        raise ValueError('logarithmic/degenerate shapes are outside this scan')
    yc=0.5 if side==1 else (2-b)/(2*b)
    def R(y):
        return -b*y*y/(a+2)-side*(1-b)*y/(a+1)+(2-b)/(4*a)
    def potential(z):
        y=np.exp(z)
        return -R(y)*np.exp(a*z)
    zc=np.log(yc)
    hc=potential(zc)
    # A fixed finite grid; small-amplitude/endpoint limits are not included.
    eta=1/(1+np.exp(-np.linspace(-9,9,n)))
    energies=hc+max(1.,abs(hc))*eta/(1-eta) if a>-2 else hc*(1-eta)
    if a < -2:
        assert hc < 0
    nodes,weights=leggauss(order)
    theta=np.pi*(nodes+1)/2
    weights=weights*np.pi/2
    rows=[]
    for h in energies:
        def equation(z): return potential(z)-h
        zl,zr=zc-0.25,zc+0.25
        while equation(zl)<0: zl-=0.5
        while equation(zr)<0: zr+=0.5
        zl=brentq(equation,zl,zc,xtol=1e-13)
        zr=brentq(equation,zc,zr,xtol=1e-13)
        z=(zl+zr)/2+(zr-zl)*np.cos(theta)/2
        xsq=np.exp(-a*z)*(h-np.array([potential(t) for t in z]))
        if np.min(xsq)<=0: raise ArithmeticError('nonpositive oval width')
        log_weight=.5*np.log(xsq)+a*z+np.log(np.sin(theta))
        ww=weights*np.exp(log_weight-np.max(log_weight))
        v=np.dot(ww,np.exp(z))/ww.sum()
        w=np.dot(ww,np.exp(-z))/ww.sum()
        rows.append([float(h),float(v),float(w),float(zl),float(zr)])
    return np.asarray(rows)

def crossing_deltas(f):
    s=np.sign(np.diff(f,axis=1))
    result=np.zeros_like(f)
    result[:,:-1]+=s
    result[:,1:]-=s
    return result

def arrangement(upper,lower):
    """Sweep line levels through both sampled curves, sharing epsilon.

    Upper: w - m v - c. Lower: -w + m v - c.
    Sampling slopes between every pairwise projection-order breakpoint
    visits every open ordering cell for these finite vertices.
    """
    n=len(upper)
    v=np.r_[upper[:,1],-lower[:,1]]
    w=np.r_[upper[:,2],-lower[:,2]]
    i,j=np.triu_indices(len(v),1)
    good=np.abs(v[i]-v[j])>1e-13
    critical=np.unique((w[i[good]]-w[j[good]])/(v[i[good]]-v[j[good]]))
    slopes=np.r_[critical[0]-max(1.,abs(critical[0])),(critical[:-1]+critical[1:])/2,critical[-1]+max(1.,abs(critical[-1]))]
    best=None
    target=None
    maxima={'upper':0,'total':0,'upper_with_lower':0}
    for start in range(0,len(slopes),512):
        m=slopes[start:start+512]
        f=w[None,:]-m[:,None]*v[None,:]
        order=np.argsort(f,axis=1)
        levels=np.take_along_axis(f,order,axis=1)
        du=np.c_[crossing_deltas(f[:,:n]),np.zeros((len(m),n))]
        dl=np.c_[np.zeros((len(m),n)),crossing_deltas(f[:,n:])]
        cu=np.cumsum(np.take_along_axis(du,order,axis=1),axis=1)[:,:-1]
        cl=np.cumsum(np.take_along_axis(dl,order,axis=1),axis=1)[:,:-1]
        c=(levels[:,:-1]+levels[:,1:])/2
        # Reject numerically unresolved line cells. This is a NUM filter,
        # not a quadrature error enclosure.
        gap=np.diff(levels,axis=1)/(1+np.abs(c)+np.abs(m[:,None]))
        valid=gap>2e-9
        maxima['upper']=max(maxima['upper'],int(np.max(np.where(valid,cu,0))))
        maxima['total']=max(maxima['total'],int(np.max(np.where(valid,cu+cl,0))))
        maxima['upper_with_lower']=max(maxima['upper_with_lower'],int(np.max(np.where(valid & (cl>=1),cu,0))))
        # Prefer a (3,1) positive control over a (2,2) roundoff flag.
        score=10*np.minimum(cu,4)+np.minimum(cl,1)
        score=np.where(valid,score,-1)
        ix=np.unravel_index(np.argmax(score),score.shape)
        rank=float(score[ix])
        record={'m':float(m[ix[0]]),'c':float(c[ix]),'upper':int(cu[ix]),'lower':int(cl[ix]),'normalized_level_gap':float(gap[ix])}
        if best is None or rank>best[0]: best=(rank,record)
        mask=valid & (((cu>=4)&(cl>=1)) | ((cl>=4)&(cu>=1)) | (cu>=5) | (cl>=5))
        if np.any(mask):
            ix=np.argwhere(mask)[0]
            target={'m':float(m[ix[0]]),'c':float(c[tuple(ix)]),'upper':int(cu[tuple(ix)]),'lower':int(cl[tuple(ix)]),'normalized_level_gap':float(gap[tuple(ix)])}
            break
    return {'slope_cells':len(slopes),'maxima':maxima,'best':best[1],'five_candidate':target}

def shapes():
    # The published four-cycle shape, generic hemicycle samples, and a
    # focused sample on either side of each center's cubic-focus line.
    raw=[(-671/450,7/15,'published_Yu_Zeng')]
    for b in (1/3,2/3,1.,4/3,5/3):
        for a in (-.4,-.75,-1.25,-1.75,-2.25,-3.,-4.):
            raw.append((a,b,'geometry_grid'))
        for line in (-(5*b+2)/3,-(12-5*b)/3):
            for da in (-.03,.03):
                raw.append((line+da,b,'near_focus_line'))
    unique={}
    for a,b,label in raw:
        if a in (-2,-1,0): continue
        unique[(a,b)]=label
    return [(a,b,label) for (a,b),label in unique.items()]

def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--control-only',action='store_true')
    parser.add_argument('--n',type=int,default=41)
    parser.add_argument('--order',type=int,default=160)
    args=parser.parse_args()
    selected=shapes()[:1] if args.control_only else shapes()
    output={'evidence':'NUM only; finite curves and finite shapes; no exclusion theorem',
        'grid':{'points_per_annulus':args.n,'gauss_order':args.order,'logit_eta_range':[-9,9]},
        'versions':{'numpy':np.__version__,'scipy':scipy.__version__},'records':[]}
    begin=time.perf_counter()
    out=ROOT/'data'/('control_search.json' if args.control_only else 'moment_search.json')
    out.parent.mkdir(parents=True,exist_ok=True)
    for a,b,label in selected:
        r={'a':a,'b':b,'label':label}
        try:
            upper=profile(a,b,1,args.n,args.order)
            lower=profile(a,b,-1,args.n,args.order)
            r.update(arrangement(upper,lower))
            r['upper_profile']=upper.tolist();r['lower_profile']=lower.tolist()
        except (ValueError,ArithmeticError,RuntimeError) as exc:
            r['error']=str(exc)
        output['records'].append(r)
        output['wall_seconds']=time.perf_counter()-begin
        out.write_text(json.dumps(output,indent=2)+'\n')
        print(json.dumps({k:v for k,v in r.items() if not k.endswith('_profile')}),flush=True)
        if r.get('five_candidate'):
            print('STOP: candidate requires independent high-precision replay',flush=True)
            break
    print('saved',out,'seconds',output['wall_seconds'],flush=True)

if __name__=='__main__': main()
