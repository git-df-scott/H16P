"""Targeted push: drive the interior minimum of D/r (between two origin roots) through zero by
finite-difference gradient descent on the free coefficients, keeping the existing roots."""
import json, sys, numpy as np, retmap as rm
recs=[json.loads(l) for l in open(sys.argv[1])]
recs=[r for r in recs if r['total']>=4]; recs.sort(key=lambda r:-r['score']); c0=np.array(recs[0]['coef'])
DIMS=[3,4,7,8,9,10,11]
RT=(1e-12,1e7,5e3,500000)
def profile(c, rmin=0.05, rmax=3e4, n=90):
    rad=np.geomspace(rmin,rmax,n)
    R,T,st=rm.returns(c[None],np.array([[0.,0.]]),np.array([[1.,0.]]),rad[None],*RT)
    ok=st[0]==0; k=n if ok.all() else int(np.argmin(ok)); return rad[:k], (R[0,:k]-rad[:k])/rad[:k]
def objective(c):
    rad,q=profile(c)
    idx=rm.count_signs(rad,q*rad)
    roots=[float(np.sqrt(rad[i]*rad[i+1])) for i in idx]
    if len(roots)<3: return None, roots, None
    i1,i2=idx[1],idx[2]                      # U root between rad[i1],rad[i1+1]; S root between rad[i2],rad[i2+1]
    inner=np.arange(i1+2, i2)                # strictly interior grid points
    if len(inner)<3: return None, roots, None
    j=inner[np.argmin(q[inner])]
    if j==inner[0] or j==inner[-1]:          # minimum at the window edge: no interior extremum
        return None, roots, None
    fine=np.geomspace(rad[j-1], rad[j+1], 21)
    R,T,st=rm.returns(c[None],np.array([[0.,0.]]),np.array([[1.,0.]]),fine[None],*RT)
    qf=(R[0]-fine)/fine; qf=np.where(st[0]==0,qf,np.inf); jj=np.argmin(qf)
    return float(qf[jj]), roots, float(fine[jj])
c=c0.copy(); scale=np.where(np.abs(c)>1e-6,np.abs(c),1.0)
m,roots,rm_=objective(c); print("start: min q =",m,"at r=",rm_,"roots",np.round(roots,3).tolist(), flush=True)
step=0.02
for it in range(40):
    g=np.zeros(12); h=1e-4
    for d in DIMS:
        cp=c.copy(); cp[d]+=h*scale[d]; mp,_,_=objective(cp)
        if mp is None: g[d]=0; continue
        g[d]=(mp-m)/(h)
    gn=np.linalg.norm(g[DIMS])
    if gn==0: print("zero gradient"); break
    direction=-g/gn
    accepted=False
    for s in (step, step/2, step/4, step/8, step/16):
        cn=c.copy(); cn[DIMS]+=s*scale[DIMS]*direction[DIMS]
        mn,rn,rr=objective(cn)
        if mn is not None and len(rn)>=3 and mn<m:
            c,m,roots,rm_=cn,mn,rn,rr; accepted=True; step=min(s*1.5,0.05); break
    print(f"it {it}: min q = {m:+.4e} at r={rm_:.3f} roots={np.round(roots,3).tolist()} step={step:.4f} accepted={accepted}", flush=True)
    if not accepted: step/=4
    if m<0:
        print("CROSSED ZERO: interior minimum negative -> new pair"); break
    if step<1e-6: break
print("final coef:", c.tolist())
rad,q=profile(c, 0.02, 3e4, 200); idx=rm.count_signs(rad,q*rad)
print("origin roots (dense):", [(round(float(np.sqrt(rad[i]*rad[i+1])),4),'S' if q[i]>0 else 'U') for i in idx])
json.dump(dict(coef=c.tolist(), min_q=m, roots=[float(np.sqrt(rad[i]*rad[i+1])) for i in idx]), open('data/push_fold_result.json','w'))
