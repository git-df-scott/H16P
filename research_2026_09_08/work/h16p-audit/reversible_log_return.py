"""Positive-time desingularized half-returns in v=asinh(x/|y|), z=log|y|.
Numerical only. A closed orbit cannot cross y=0 when epsilon0 !=0.
"""
from pathlib import Path
import numpy as np,json
from scipy.integrate import solve_ivp
P=Path(__file__).resolve().parent
COUNTS={'log_halves':0}
def half(a,b,e0,e1,e2,side,T,sign,rtol=3e-12):
    COUNTS['log_halves']+=1;C=(b-2)/4;tc=.5 if side==1 else (2-b)/(2*b)
    if not T>tc:raise ValueError('start outside center height')
    def rhs(t,s):
        from scipy.special import expit
        v,z=s;th=np.tanh(v);lc=abs(v)+np.log1p(np.exp(-2*abs(v)))-np.log(2.)
        ll=np.logaddexp(2*lc,-2*z);pp=np.exp(2*lc-ll);qq=np.exp(-2*z-ll)
        basic=(a+2)*th*th*pp+b*np.exp(-ll)+side*e2*th*np.exp(lc-ll)+side*(1-b)*np.exp(-z-ll)+e1*th*np.exp(lc-z-ll)+C*qq
        le=np.log(abs(e0))-2*z-ll+lc if e0 else -np.inf
        inv=expit(-le);scaled=expit(le);es=np.sign(e0)
        return sign*np.array([basic*inv-side*es*th*scaled,-2*th*pp*inv+side*es*scaled])
    def section(t,s):return s[0]
    section.terminal=True;section.direction=-sign
    def guard(t,s):return 1000-max(abs(s[0]),abs(s[1]))
    guard.terminal=True;guard.direction=-1
    initial=[0.,np.log(T)];dt=1e-7
    kick=solve_ivp(rhs,[0,dt],initial,method='DOP853',rtol=rtol,atol=rtol*.01)
    sol=solve_ivp(rhs,[dt,10000],kick.y[:,-1],method='DOP853',rtol=rtol,atol=rtol*.01,events=[section,guard],max_step=.3)
    if not kick.success or not sol.success or not len(sol.t_events[0]):return dict(status='unresolved',message=sol.message,guard=len(sol.t_events[1]))
    end=sol.y_events[0][0];tend=float(np.exp(end[1]))
    if not end[1]<np.log(tc):return dict(status='wrong_section',height=tend)
    return dict(status='passed',log_height=float(end[1]),height=tend,nfev=kick.nfev+sol.nfev,time=float(sol.t[-1]))
def paired(a,b,e0,e1,e2,side,T,rtol=3e-12):
    f=half(a,b,e0,e1,e2,side,T,1,rtol);r=half(a,b,e0,e1,e2,side,T,-1,rtol)
    out=dict(side=side,height=T,forward=f,backward=r)
    if f['status']==r['status']=='passed':out['log_difference']=f['log_height']-r['log_height']
    return out
if __name__=='__main__':
    out={'scope':'Independent Cartesian endpoint checks plus unperturbed very-large-orbit identity controls','records':[]}
    for a,b in [(-1.,1.),(-1.25,.5),(-.75,1.5)]:
        for side in [1,-1]:
            tc=.5 if side==1 else (2-b)/(2*b)
            for T in [tc*2,tc*20]:
                pars=(a,b,1e-4,2e-4,-3e-4);row=paired(*pars,side,T);checks=[]
                for sign in [1,-1]:
                    def rhs(t,s):
                        x,y=s;return sign*np.array([(b-2)/4+(1-b)*y+a*x*x+b*y*y+pars[3]*x+pars[4]*x*y,pars[2]-2*x*y])
                    def section(t,s):return s[0]
                    section.terminal=True;section.direction=-sign
                    kick=solve_ivp(rhs,[0,1e-7],[0.,side*T],method='DOP853',rtol=3e-13,atol=1e-15)
                    sol=solve_ivp(rhs,[1e-7,100],kick.y[:,-1],method='DOP853',rtol=3e-13,atol=1e-15,events=section,max_step=.1)
                    assert sol.success and len(sol.t_events[0]);height=abs(sol.y_events[0][0][1]);logresult=row['forward' if sign==1 else 'backward'];assert logresult['status']=='passed'
                    checks.append(dict(sign=sign,cartesian_height=float(height),log_height_difference=float(abs(np.log(height)-logresult['log_height']))))
                out['records'].append(dict(parameters=pars,row=row,checks=checks))
    out['large_controls']=[dict(a=a,b=b,row=paired(a,b,0,0,0,side,T)) for a,b in [(-1,1),(-1.25,.5)] for side in [1,-1] for T in [1e10,1e30,1e60]]
    out['counts']=COUNTS.copy();(P/'reversible_log_return_v3_check.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(dict(max_cartesian_log_error=max(q['log_height_difference'] for r in out['records'] for q in r['checks']),max_large_return_difference=max(abs(q['row'].get('log_difference',1)) for q in out['large_controls']),counts=COUNTS)))
