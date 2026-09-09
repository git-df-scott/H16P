"""Finite Q4 three-anchor continuation with event-corrected sensitivities.
Nonvalidated numerical root fitting and sampled complete half-passages.
"""
from pathlib import Path
from fractions import Fraction as F
import json,time
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import root,brentq
from q4_modulus_probe import compact
P=Path(__file__).resolve().parent
COUNTS={'sensitivity_halves':0,'plain_halves':0}

def half(r,eps,controls,radius,sign,sensitivity=False,rtol=3e-13):
    COUNTS['sensitivity_halves' if sensitivity else 'plain_halves']+=1
    tau,u,w=eps*np.array(controls);v=-eps
    def rhs(t,z):
        x,y=z[:2]
        p=tau*x-y-(2+r*r)*x*x+(2*r+2*u+w)*x*y+y*y
        q=x+tau*y+(r+u)*x*x+(-1-3*r*r+v)*x*y-(r+u)*y*y
        if not sensitivity:return sign*np.array([p,q])
        J=np.array([[tau-2*(2+r*r)*x+(2*r+2*u+w)*y,-1+(2*r+2*u+w)*x+2*y],[1+2*(r+u)*x+(-1-3*r*r+v)*y,tau+(-1-3*r*r+v)*x-2*(r+u)*y]])
        source=eps*np.array([[x,2*x*y,x*y],[y,x*x-y*y,0.]])
        ds=J@z[2:].reshape(2,3)+source
        return sign*np.r_[p,q,ds.ravel()]
    initial=np.r_[0.,-radius,np.zeros(6)] if sensitivity else np.array([0.,-radius])
    def section(t,z):return z[0]
    section.terminal=True;section.direction=-sign
    def guard(t,z):return 1e10-np.linalg.norm(z[:2])
    guard.terminal=True;guard.direction=-1
    dt=min(1e-6,1e-6/max(1.,np.linalg.norm(rhs(0,initial)[:2])))
    kick=solve_ivp(rhs,[0,dt],initial,method='DOP853',rtol=rtol,atol=rtol*.002)
    sol=solve_ivp(rhs,[dt,100],kick.y[:,-1],method='DOP853',rtol=rtol,atol=rtol*.002,events=[section,guard],max_step=.15)
    if not kick.success or not sol.success or not len(sol.t_events[0]):return {'status':'unresolved','reason':sol.message,'guard':len(sol.t_events[1])}
    end=sol.y_events[0][0];x,y=end[:2]
    if not 0<y<1:return {'status':'unresolved','reason':'wrong section','y':float(y)}
    out={'status':'passed','y':float(y),'time':float(sol.t_events[0][0]),'nfev':kick.nfev+sol.nfev}
    if sensitivity:
        f=rhs(0,end)[:2]/sign;out['gradient']=(end[2:].reshape(2,3)[1]-f[1]/f[0]*end[2:].reshape(2,3)[0]).tolist()
    return out

def difference(r,eps,controls,radius,sensitivity=False,rtol=3e-13):
    f=half(r,eps,controls,radius,1,sensitivity,rtol);b=half(r,eps,controls,radius,-1,sensitivity,rtol)
    row={'radius':float(radius),'forward':f,'backward':b}
    if f['status']==b['status']=='passed':
        row['difference']=f['y']-b['y']
        if sensitivity:row['gradient']=(np.array(f['gradient'])-b['gradient']).tolist()
    return row

def anchor(r,t):
    d=1+r*r;target=-1+t*(1-1/np.sqrt(d))
    def fun(R):
        y=-R;D=1-2*d*y+d*y*y;N=3*d*y*(1-y)-1+d*y**3
        return N/D**1.5-target
    return brentq(fun,0,1e8)

def main(rvalues=(.5,1.,2.),epsvalues=(.001,.01,.03,.1,.2,.4),outname='q4_finite_continuation.json',warm=None):
    out={'scope':'NUM actual quadratic fields; three fitted return equalities plus sampled crossings; no validated roots or global exclusion','records':[],'failures':[]};start=time.perf_counter()
    for r in rvalues:
        anchors=[anchor(r,t) for t in [.25,.5,.75]]
        M=np.array([compact(r,t,512)[0] for t in [.25,.5,.75]])
        controls=np.linalg.solve(M[:,[0,1,3]],M[:,2])
        if warm is not None and r in warm:controls=np.array(warm[r])
        for eps in epsvalues:
            calls=[];cached=None
            def evaluate(c):
                nonlocal cached
                if cached is not None and np.array_equal(c,cached[0]):return cached[1:]
                rows=[difference(r,eps,c,R,True) for R in anchors]
                calls.append({'controls':c.tolist(),'rows':rows})
                if any('difference' not in row for row in rows):raise ArithmeticError('anchor passage unresolved')
                residual=np.array([row['difference'] for row in rows])/eps
                jac=np.array([row['gradient'] for row in rows])/eps
                cached=(c.copy(),residual,jac);return residual,jac
            try:
                fit=root(lambda c:evaluate(c)[0],controls,jac=lambda c:evaluate(c)[1],method='hybr',options={'xtol':2e-8,'maxfev':80})
                residual=evaluate(fit.x)[0]
                if max(abs(residual))>2e-9:raise ArithmeticError('fit residual too large')
                controls=fit.x
                radii=sorted(set(np.geomspace(.01,1e7,61).tolist()+[q*factor for q in anchors for factor in [.9,1.1]]))
                rows=[difference(r,eps,controls,R) for R in radii]
                brackets=[[a['radius'],b['radius']] for a,b in zip(rows,rows[1:]) if 'difference' in a and 'difference' in b and a['difference']*b['difference']<0 and min(abs(a['difference']),abs(b['difference']))>1e-11]
                params={'r':str(F(str(r))),'tau':str(F(str(eps))*F(str(controls[0]))),'u':str(F(str(eps))*F(str(controls[1]))),'v':str(-F(str(eps))),'w':str(F(str(eps))*F(str(controls[2])))}
                record={'r':r,'epsilon':eps,'normalized_controls':controls.tolist(),'rational_parameters':params,'anchors':anchors,'scaled_fit_residual':residual.tolist(),'solver_success':bool(fit.success),'solver_message':fit.message,'fit_calls':calls,'profile':rows,'brackets':brackets,'unresolved':sum('difference' not in row for row in rows)}
                out['records'].append(record);print(json.dumps({'r':r,'eps':eps,'controls':controls.tolist(),'brackets':brackets,'unresolved':record['unresolved']}),flush=True)
            except Exception as e:
                out['failures'].append({'r':r,'epsilon':eps,'reason':repr(e),'calls':calls});print(json.dumps({'r':r,'eps':eps,'failure':repr(e)}),flush=True);break
            out['counts']=COUNTS.copy();out['wall_seconds']=time.perf_counter()-start
            (P/outname).write_text(json.dumps(out,indent=2)+'\n')
    out['counts']=COUNTS.copy();out['wall_seconds']=time.perf_counter()-start;(P/outname).write_text(json.dumps(out,indent=2)+'\n')
if __name__=='__main__':main()
