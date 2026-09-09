"""Exact rational Hopf fields, numerical two-sided passage profiles.
No interval validation; no global count from finite section samples.
"""
from fractions import Fraction as F
from pathlib import Path
import json,time
import numpy as np
from scipy.integrate import solve_ivp
ROOT=Path(__file__).resolve().parent

def field(k,m=F(1539,1000),a=F(-7,4),b=F(1,3)):
    x=-k/(2-a+b*k*k);y=F(1,2)+k*x
    e0=2*x*y;tau=(a-1)*e0;e2=-tau*m;e1=-e2*y-2*(a-1)*x
    assert (b-2)/4+(1-b)*y+a*x*x+b*y*y+e1*x+e2*x*y==0
    assert e1+e2*y+2*(a-1)*x==0
    determinant=(2*a*x+e1+e2*y)*(-2*x)-(1-b+2*b*y+e2*x)*(-2*y)
    assert determinant>0
    exact={key:str(value) for key,value in dict(a=a,b=b,k=k,m=m,epsilon0=e0,epsilon1=e1,epsilon2=e2,tau=tau,x0=x,y0=y,determinant=determinant).items()}
    coeff=np.array([4*float(b),4*float(1-b),float(b-2+2*e2*e0),float(2*e1*e0),float(a*e0*e0)])
    roots=np.roots(coeff)
    eq=[]
    for yy in roots:
        if abs(yy.imag)<1e-8:
            yy=float(yy.real);xx=float(e0)/(2*yy)
            J=np.array([[2*float(a)*xx+float(e1)+float(e2)*yy,1-float(b)+2*float(b)*yy+float(e2)*xx],[-2*yy,-2*xx]])
            eq.append({'point':[xx,yy],'trace':float(np.trace(J)),'determinant':float(np.linalg.det(J)),'discriminant':float(np.trace(J)**2-4*np.linalg.det(J))})
    if len(eq)!=2:raise ArithmeticError('two-real-equilibrium numerical gate')
    upper=np.array([float(x),float(y)]);lower=np.array(min(eq,key=lambda r:r['point'][1])['point'])
    vector=upper-lower;distance=float(np.linalg.norm(vector));e=vector/distance;n=np.array([e[1],-e[0]])
    af,bf,e0f,e1f,e2f=map(float,(a,b,e0,e1,e2))
    def rhs(t,z):
        xx,yy=z
        return np.array([(bf-2)/4+(1-bf)*yy+af*xx*xx+bf*yy*yy+e1f*xx+e2f*xx*yy,e0f-2*xx*yy])
    return exact,eq,upper,lower,e,n,distance,rhs

def half(model,r,side,sgn,rtol):
    exact,eq,upper,lower,e,n,distance,rhs=model
    initial=(upper if side==1 else lower)+side*r*e
    velocity=sgn*rhs(0,initial);normal=float(n@velocity)
    if abs(normal)<1e-14:return {'status':'unresolved','reason':'initial transversality'}
    def fun(t,z):return sgn*rhs(t,z)
    def section(t,z):return float(n@(z-upper))
    section.terminal=True;section.direction=-np.sign(normal)
    def guard(t,z):return 1e9-np.linalg.norm(z)
    guard.terminal=True;guard.direction=-1
    dt=min(1e-6,1e-6/max(1.,np.linalg.norm(velocity)))
    kick=solve_ivp(fun,[0,dt],initial,method='DOP853',rtol=rtol,atol=rtol*.005)
    if not kick.success:return {'status':'unresolved','reason':'kick failed'}
    sol=solve_ivp(fun,[dt,100],kick.y[:,-1],method='DOP853',events=[section,guard],rtol=rtol,atol=rtol*.005,max_step=.15)
    if not sol.success or len(sol.t_events[0])!=1:return {'status':'unresolved','reason':sol.message,'guard':len(sol.t_events[1]),'last':sol.y[:,-1].tolist()}
    endpoint=sol.y_events[0][0];coordinate=float(e@(endpoint-upper))
    if not -distance<coordinate<0:return {'status':'unresolved','reason':'return outside between-focus segment','coordinate':coordinate}
    return {'status':'half_passage','coordinate':coordinate,'time':float(sol.t_events[0][0]),'section_error':float(n@(endpoint-upper)),'nfev':sol.nfev+kick.nfev}

def profile(model,radii,side,rtol):
    out=[]
    for r in radii:
        plus=half(model,float(r),side,1,rtol);minus=half(model,float(r),side,-1,rtol)
        row={'r':float(r),'side':side,'forward':plus,'backward':minus}
        if plus['status']==minus['status']=='half_passage':row['difference']=plus['coordinate']-minus['coordinate']
        out.append(row)
    return out

def count(rows):
    brackets=[]
    for p,q in zip(rows,rows[1:]):
        if 'difference' in p and 'difference' in q and p['difference']*q['difference']<0 and min(abs(p['difference']),abs(q['difference']))>1e-10:brackets.append([p['r'],q['r']])
    return brackets

def main():
    begin=time.perf_counter();records=[]
    for k in [F(1,10000),F(1,1000),F(1,100),F(1,20),F(1,10),F(1,5)]:
        model=field(k)
        upper=profile(model,np.geomspace(.1,1e5,37),1,2e-11)
        lower=profile(model,np.geomspace(1,1e6,25),-1,2e-11)
        row={'field':model[0],'equilibria_numerical':model[1],'upper':upper,'lower':lower,'upper_brackets':count(upper),'lower_brackets':count(lower),'unresolved':sum('difference' not in r for r in upper+lower)}
        records.append(row)
        print(json.dumps({k:v for k,v in row.items() if k not in ('upper','lower','equilibria_numerical')}),flush=True)
        (ROOT/'reversible_finite_hopf.json').write_text(json.dumps({'evidence':'NUM sampled two-sided passages only; exact rational equilibrium and trace-zero identities','records':records,'wall_seconds':time.perf_counter()-begin},indent=2)+'\n')
        if len(row['upper_brackets'])>=3 and row['lower_brackets']:break
if __name__=='__main__':main()
