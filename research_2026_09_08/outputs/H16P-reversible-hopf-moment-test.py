"""Constrain the common first-order perturbation to one center's Hopf tangent.
A sampled precursor needs three compact crossings there plus a remote one.
No claim about unsampled energies, finite perturbations or endpoint cycles.
"""
from pathlib import Path
import json,hashlib
import numpy as np
ROOT=Path(__file__).resolve().parent
paths=[ROOT.parent/'H16P/reversible_reseed/data/moment_search.json',ROOT/'resonant_shape_probe.json']

def test(a,b,up,lo,side):
    yc=.5 if side==1 else (2-b)/(2*b)
    curves=[];critical=[]
    for sigma,arr in [(1,up),(-1,lo)]:
        A=sigma*arr[:,2]-side/yc
        B=sigma*arr[:,1]-side*yc
        curves.append((A,B))
        valid=abs(B)>1e-11
        critical.extend((A[valid]/B[valid]).tolist())
    critical=np.unique(critical)
    slopes=np.r_[critical[0]-max(1,abs(critical[0])),(critical[:-1]+critical[1:])/2,critical[-1]+max(1,abs(critical[-1]))]
    rows=[];maximum=0;best=None
    for m in slopes:
        fs=[A-m*B for A,B in curves]
        margins=[float(min(abs(f))/(1+abs(m))) for f in fs]
        if min(margins)<1e-10:continue
        counts=[int(sum(f[:-1]*f[1:]<0)) for f in fs]
        home=0 if side==1 else 1;remote=1-home
        maximum=max(maximum,counts[home])
        score=10*counts[home]+min(1,counts[remote])
        row={'m':float(m),'c':float(side*(1/yc-m*yc)),'crossings':counts,'normalized_margins':margins,'first_sample_signs':[int(np.sign(f[0])) for f in fs]}
        if best is None or score>best[0]:best=(score,row)
        if counts[home]>=3 and counts[remote]>=1:rows.append(row)
    return {'a':a,'b':b,'hopf_side':side,'max_home_crossings':maximum,'best':None if best is None else best[1],'candidate_directions':rows}

def main():
    rows=[];sources=[]
    for path in paths:
        data=json.loads(path.read_text());sources.append({'path':str(path),'sha256':hashlib.sha256(path.read_bytes()).hexdigest()})
        for r in data['records']:
            if 'upper_profile' in r:up,lo=r['upper_profile'],r['lower_profile']
            elif isinstance(r.get('upper'),list):up,lo=r['upper'],r['lower']
            else:continue
            for side in [1,-1]:rows.append(test(r['a'],r['b'],np.array(up),np.array(lo),side))
    out={'evidence':'Finite saved-profile NUM test only; tangent Hopf constraint, not exact finite-perturbation Hopf','constraint':'c=side*(1/yc-m*yc); eps0=tau/(a-1),eps1=-tau*c,eps2=-tau*m','sources':sources,'records':rows,'candidate_count':sum(len(r['candidate_directions']) for r in rows)}
    (ROOT/'reversible_hopf_moment_test.json').write_text(json.dumps(out,indent=2)+'\n')
    print('shape-side tests',len(rows),'candidate directions',out['candidate_count'])
    from collections import Counter
    print('max-home-crossing histogram',Counter(r['max_home_crossings'] for r in rows))
    for r in rows:
        if r['max_home_crossings']>=2:print(json.dumps(r))
if __name__=='__main__':main()
