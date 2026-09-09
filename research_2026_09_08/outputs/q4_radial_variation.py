"""First radial derivative of paired Cartesian half-return maps."""
import json
import numpy as np
from scipy.integrate import solve_ivp
from q4_finite_continuation import P,difference
COUNTS={'radial_halves':0}
def half(r,eps,c,R,sign,rtol=3e-13):
    COUNTS['radial_halves']+=1
    tau,u,w=eps*np.array(c);v=-eps
    def raw(z):
        x,y=z[:2]
        f=np.array([tau*x-y-(2+r*r)*x*x+(2*r+2*u+w)*x*y+y*y,x+tau*y+(r+u)*x*x+(-1-3*r*r+v)*x*y-(r+u)*y*y])
        J=np.array([[tau-2*(2+r*r)*x+(2*r+2*u+w)*y,-1+(2*r+2*u+w)*x+2*y],[1+2*(r+u)*x+(-1-3*r*r+v)*y,tau+(-1-3*r*r+v)*x-2*(r+u)*y]])
        return f,J
    def rhs(t,z):
        f,J=raw(z);return sign*np.r_[f,J@z[2:]]
    def section(t,z):return z[0]
    section.terminal=True;section.direction=-sign
    def guard(t,z):return 1e10-np.linalg.norm(z[:2])
    guard.terminal=True;guard.direction=-1
    initial=np.array([0.,-R,0.,-1.]);dt=min(1e-6,1e-6/max(1,np.linalg.norm(raw(initial)[0])))
    kick=solve_ivp(rhs,[0,dt],initial,method='DOP853',rtol=rtol,atol=rtol*.002)
    sol=solve_ivp(rhs,[dt,100],kick.y[:,-1],method='DOP853',rtol=rtol,atol=rtol*.002,events=[section,guard],max_step=.15)
    if not kick.success or not sol.success or not len(sol.t_events[0]):return dict(status='unresolved',message=sol.message,guard=len(sol.t_events[1]))
    z=sol.y_events[0][0];f,J=raw(z)
    if not 0<z[1]<1:return dict(status='unresolved',reason='wrong section',point=z[:2].tolist())
    return dict(status='passed',y=float(z[1]),radial_derivative=float(z[3]-f[1]/f[0]*z[2]),time=float(sol.t[-1]))
def radial(r,eps,c,R,rtol=3e-13):
    f=half(r,eps,c,R,1,rtol);b=half(r,eps,c,R,-1,rtol)
    out=dict(radius=R,forward=f,backward=b)
    if f['status']==b['status']=='passed':out.update(difference=f['y']-b['y'],radial_derivative=f['radial_derivative']-b['radial_derivative'])
    return out
if __name__=='__main__':
    source=json.loads((P/'q4_finite_extended.json').read_text());rows=[]
    for a in source['records']:
        if a['epsilon']!=4:continue
        for R in [a['anchors'][0]*1.2,a['anchors'][-1]*1.2]:
            row=radial(a['r'],a['epsilon'],a['normalized_controls'],R);checks=[]
            for h in [1e-3,1e-4]:
                plus=difference(a['r'],a['epsilon'],a['normalized_controls'],R+h)['difference'];minus=difference(a['r'],a['epsilon'],a['normalized_controls'],R-h)['difference']
                checks.append(dict(step=h,central_derivative=(plus-minus)/(2*h),error=abs((plus-minus)/(2*h)-row['radial_derivative'])))
            rows.append(dict(r=a['r'],row=row,checks=checks))
    (P/'q4_radial_variation_check.json').write_text(json.dumps(dict(scope='Analytic variational derivative vs independent central differences; nonvalidated',records=rows,counts=COUNTS),indent=2)+'\n')
    print(json.dumps(dict(max_error_h1e4=max(a['checks'][1]['error'] for a in rows),all_error_decrease=all(a['checks'][1]['error']<a['checks'][0]['error'] for a in rows))))
