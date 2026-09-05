"""D2 step 1: divergence structure along the homoclinic loops around order-two weak foci found in F11.
For each focus-type loop (a,b,l): integrate the saddle's unstable separatrix around the loop, record
div along it (sign changes, integral), and test the Dulac candidate D_k = (1+ax+by) div X + k (aP+bQ)
for sign-definiteness on the loop interior."""
import json, numpy as np
from scipy.integrate import solve_ivp
recs=[json.loads(l) for l in open('audit/fable_engine/data/F11_loops.jsonl')]
foc=[r for r in recs if abs(r.get('eta2',0))>1e-6 and r['a']<0]
def sys_(a,b,l):
    m=a*(b+2*l)/(l+1)
    P=lambda x,y: -y+l*x*x+m*x*y+y*y; Q=lambda x,y: x+a*x*x+b*x*y
    div=lambda x,y: 2*l*x+m*y+b*x
    return P,Q,div,m
out=[]
for r in foc:
    a,b,l=r['a'],r['b'],r['l']; P,Q,div,m=sys_(a,b,l); xs,ys=r['saddle']
    J=np.array([[2*l*xs+m*ys, -1+m*xs+2*ys],[1+2*a*xs+b*ys, b*xs]])
    w,v=np.linalg.eig(J); iu=int(np.argmax(w.real)); vu=v[:,iu].real; vu/=np.linalg.norm(vu)
    best=None
    for sgn in (1,-1):
        p0=np.array([xs,ys])+sgn*1e-7*vu
        ev=lambda t,u: np.hypot(u[0]-xs,u[1]-ys)-1e-4; ev.terminal=True; ev.direction=-1
        def rhs(t,u): return [P(u[0],u[1]),Q(u[0],u[1]),div(u[0],u[1])]
        sol=solve_ivp(rhs,[0,200],[p0[0],p0[1],0.0],rtol=1e-11,atol=1e-13,events=ev,max_step=0.05,dense_output=True)
        if len(sol.t_events[0])==0: continue
        T=sol.t_events[0][0]
        if T<1.0: continue
        ts=np.linspace(0,T,4000); U=sol.sol(ts); d=div(U[0],U[1])
        # does the orbit encircle the focus (origin)? winding number
        ang=np.unwrap(np.arctan2(U[1],U[0])); wind=(ang[-1]-ang[0])/(2*np.pi)
        if abs(abs(wind)-1)>0.2: continue
        best=dict(T=float(T), int_div=float(U[2,-1]), div_min=float(d.min()), div_max=float(d.max()),
                  div_sign_changes=int(np.sum(np.sign(d[:-1])*np.sign(d[1:])<0)), wind=float(wind), X=U[0], Y=U[1])
        break
    if best is None: print(f"a={a} b={b}: loop orbit not recovered"); continue
    sigma=r['trace']; eta2=r['eta2']
    # Dulac test on the interior: sample points inside the loop polygon
    from matplotlib.path import Path
    poly=Path(np.column_stack([best['X'],best['Y']]))
    xmin,xmax,ymin,ymax=best['X'].min(),best['X'].max(),best['Y'].min(),best['Y'].max()
    g=np.random.default_rng(0); pts=np.column_stack([g.uniform(xmin,xmax,20000),g.uniform(ymin,ymax,20000)])
    inside=pts[poly.contains_points(pts)]; inside=inside[np.hypot(inside[:,0],inside[:,1])>0.02]
    X,Y=inside[:,0],inside[:,1]; L=1+a*X+b*Y
    dul={}
    for k in (m/a if a!=0 else 0.0, -(2*l+b)/b if b!=0 else 0.0, 1.0, -1.0, 2.0):
        D=L*div(X,Y)+k*(a*P(X,Y)+b*Q(X,Y)); dul[round(float(k),4)]=(float((D>0).mean()))
    out.append(dict(a=a,b=b,l=l,m=m,sigma=sigma,eta2=eta2,**{k:v for k,v in best.items() if k not in('X','Y')},dulac_frac_pos=dul,line_min=float(L.min()),line_max=float(L.max())))
    print(f"a={a:+.2f} b={b:+.2f}: sigma={sigma:+.4f} eta2={eta2:+.3e} intdiv={best['int_div']:+.4f} div range [{best['div_min']:+.3f},{best['div_max']:+.3f}] sign changes {best['div_sign_changes']} | 1+ax+by on interior in [{L.min():+.3f},{L.max():+.3f}] | Dulac frac>0: {dul}")
json.dump(out,open('audit/fable_engine/data/D2_loop_div.json','w'))
