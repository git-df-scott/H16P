"""Numerical saddle branch shooting; finite offsets, not validated manifolds."""
from pathlib import Path
import json,time
import numpy as np
from scipy.integrate import solve_ivp
P=Path(__file__).resolve().parent

def field(rec):
    r=rec['r']; tau,u,w=rec['epsilon']*np.array(rec['normalized_controls']); v=-rec['epsilon']
    def f(t,z):
        x,y=z
        return np.array([tau*x-y-(2+r*r)*x*x+(2*r+2*u+w)*x*y+y*y,x+tau*y+(r+u)*x*x+(-1-3*r*r+v)*x*y-(r+u)*y*y])
    def jac(z):
        x,y=z
        return np.array([[tau-2*(2+r*r)*x+(2*r+2*u+w)*y,-1+(2*r+2*u+w)*x+2*y],[1+2*(r+u)*x+(-1-3*r*r+v)*y,tau+(-1-3*r*r+v)*x-2*(r+u)*y]])
    return f,jac

def shoot(f,saddle,vector,sgn,branch,offset,rtol):
    z=saddle+branch*offset*max(1,np.linalg.norm(saddle))*vector
    rhs=lambda t,z:sgn*f(t,z)
    def section(t,z):return z[0]
    section.terminal=True;section.direction=0
    def guard(t,z):return 1e8-np.linalg.norm(z)
    guard.terminal=True;guard.direction=-1
    itinerary=[];nfev=0;elapsed=0.
    for i in range(4):
        sol=solve_ivp(rhs,[0,200],z,method='DOP853',rtol=rtol,atol=rtol*.002,events=[section,guard],max_step=.15)
        nfev+=sol.nfev;elapsed+=float(sol.t[-1])
        if not sol.success or not len(sol.t_events[0]):
            return dict(status='guard' if len(sol.t_events[1]) else 'unresolved',message=sol.message,itinerary=itinerary,nfev=nfev,time=elapsed)
        end=sol.y_events[0][0];itinerary.append(end.tolist())
        if end[1]<0:return dict(status='negative_section',radius=float(-end[1]),itinerary=itinerary,nfev=nfev,time=elapsed)
        dt=min(1e-6,1e-6/max(1,np.linalg.norm(rhs(0,end))))
        kick=solve_ivp(rhs,[0,dt],end,method='DOP853',rtol=rtol,atol=rtol*.002)
        nfev+=kick.nfev;z=kick.y[:,-1]
    return dict(status='itinerary_limit',itinerary=itinerary,nfev=nfev,time=elapsed)

def main():
    records=json.loads((P/'q4_finite_extended.json').read_text())['records']
    geometry=json.loads((P/'q4_finite_extended_equilibria.json').read_text())['records']
    out={'scope':'Nonvalidated saddle branch shooting at two finite offsets and tolerances; guard is not an escape proof','records':[]}
    geometry_by_key={(a['r'],a['epsilon']):a for a in geometry}
    for rec in records:
        geo=geometry_by_key[(rec['r'],rec['epsilon'])]
        saddles=[a for a in geo['equilibria'] if a['determinant']<0]
        if not saddles:continue
        s=np.array(saddles[0]['point']);f,jac=field(rec);values,vectors=np.linalg.eig(jac(s));rows=[]
        assert np.max(np.abs(f(0,s)))<1e-7
        assert max(abs(values.imag))<1e-12
        values=values.real
        assert min(values)<0<max(values)
        for offset,rtol in [(1e-7,3e-12),(1e-9,3e-13)]:
            for j,value in enumerate(values):
                vector=vectors[:,j].real;vector*=1 if vector[0]>=0 else -1
                for branch in [-1,1]:
                    row=shoot(f,s,vector,1 if value>0 else -1,branch,offset,rtol)
                    row.update(offset=offset,rtol=rtol,eigenvalue=float(value),manifold='unstable' if value>0 else 'stable',branch=branch);rows.append(row)
        item=dict(r=rec['r'],epsilon=rec['epsilon'],saddle=s.tolist(),eigenvalues=values.tolist(),rows=rows)
        out['records'].append(item)
        print(json.dumps({k:item[k] for k in ['r','epsilon','eigenvalues']}|{'shots':[{k:a[k] for k in ['manifold','branch','status','radius'] if k in a} for a in rows[-4:]]}),flush=True)
        (P/'q4_saddle_separatrices.json').write_text(json.dumps(out,indent=2)+'\n')
if __name__=='__main__':main()
