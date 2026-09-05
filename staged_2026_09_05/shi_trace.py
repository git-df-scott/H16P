#!/usr/bin/env python3
"""Bounded, cancellation-reduced numerical trace continuation; not interval proof.

Exact rational degree-eight Lyapunov polynomial eliminates the leading
oscillatory radial drift from the measured displacement. The orbit is integrated
in polar angle, and H=integral (dV/dt)/V(initial) dt is the return diagnostic.
H=0 iff the radius returns, provided V is monotonic on the section locally.
"""
import argparse, json, math, signal, time
from pathlib import Path
import numpy as np
import sympy as S
from scipy.integrate import solve_ivp

HERE=Path(__file__).resolve().parent
LEDGER=HERE/'shi_returns.jsonl'
CAP=160

def exact_lyapunov(l,m,a,b):
    x,y=S.symbols('x y'); P=-y+l*x*x+m*x*y+y*y; Q=x+a*x*x+b*x*y
    V=(x*x+y*y)/2; etas=[]
    for k in range(3,9):
        cs=S.symbols('c:'+str(k+1)); vk=sum(cs[j]*x**(k-j)*y**j for j in range(k+1))
        known=S.Poly(S.diff(V,x)*P+S.diff(V,y)*Q,x,y)
        degree=sum(co*x**ij[0]*y**ij[1] for ij,co in known.terms() if sum(ij)==k)
        eq=degree-y*S.diff(vk,x)+x*S.diff(vk,y)
        unknown=list(cs)
        if k%2==0:
            eta=S.Symbol('eta'); eq-=eta*(x*x+y*y)**(k//2); unknown.append(eta)
        sol=S.solve(S.Poly(eq,x,y).coeffs(),unknown,dict=True)[0]
        if k%2==0: etas.append(sol[eta])
        V+=vk.subs(sol).subs({c:0 for c in cs})
        V=S.expand(V)
    deriv=S.Poly(S.expand(S.diff(V,x)*P+S.diff(V,y)*Q),x,y)
    remainder=sum(co*x**ij[0]*y**ij[1] for ij,co in deriv.terms() if sum(ij)==9)
    target=sum(etas[j]*(x*x+y*y)**(j+2) for j in range(3))+remainder
    assert S.expand(deriv.as_expr()-target)==0
    return {'V':V,'Vx':S.diff(V,x),'rem':remainder,'etas':etas,'symbols':(x,y)}

class Engine:
    def __init__(self,name,l,m,a,b):
        self.name=name; self.exact=[str(v) for v in (l,m,a,b)]
        self.l,self.m,self.a,self.b=map(float,(l,m,a,b))
        data=exact_lyapunov(l,m,a,b)
        self.etas=list(map(float,data['etas']))
        self.exact_etas=list(map(str,data['etas']))
        self.V,self.Vx,self.rem=[S.lambdify(data['symbols'],data[k],modules='math',cse=True) for k in ['V','Vx','rem']]
        self.Vsection=S.lambdify(data['symbols'][0],data['V'].subs(data['symbols'][1],0),'math')
    def evaluate(self,r,lam,rtol=2e-11,atol=2e-14,max_step=.10,tag=''):
        count=sum(1 for _ in LEDGER.open()) if LEDGER.exists() else 0
        if count>=CAP: raise RuntimeError('160-return budget exhausted')
        record={'index':count+1,'family':self.name,'l_m_a_b_exact':self.exact,'r':r,'lambda':lam,'rtol':rtol,'atol':atol,'max_step':max_step,'tag':tag,'validated':False}
        started=time.process_time()
        def timeout(*args): raise TimeoutError('10 CPU second evaluation cap')
        prior=signal.signal(signal.SIGPROF,timeout); signal.setitimer(signal.ITIMER_PROF,10)
        try:
            l,m,a,b=self.l,self.m,self.a,self.b
            def fun(theta,z):
                radius=r*math.exp(z[0]); c=math.cos(theta); s=math.sin(theta)
                x,y=radius*c,radius*s
                A=l*c**3+(m+a)*c*c*s+(1+b)*c*s*s
                B=a*c**3+(b-l)*c*c*s-m*c*s*s-s**3
                omega=1-lam*c*s+radius*B
                if omega<=.03 or radius>3: raise ValueError('polar chart lost or radius >3')
                r2=radius*radius
                drift=r2*r2*(self.etas[0]+r2*(self.etas[1]+r2*self.etas[2]))+self.rem(x,y)+lam*x*self.Vx(x,y)
                return [(lam*c*c+radius*A)/omega,drift/(self.Vsection(r)*omega),1/omega,(lam+(2*l+b)*x+m*y)/omega]
            sol=solve_ivp(fun,(0,2*math.pi),[0,0,0,0],method='DOP853',rtol=rtol,atol=[atol,atol*1e-8,atol,atol],max_step=max_step)
            if not sol.success: raise RuntimeError(sol.message)
            q,H,T,div=map(float,sol.y[:,-1]); r1=r*math.exp(q)
            record.update(status='ok',H=H,measure='deltaV_over_Vinitial',Vinitial=self.Vsection(r),Vx_section=self.Vx(r,0),raw_log_return=q,period=T,divergence_integral=div,
                          multiplier_at_cycle=math.exp(div),return_derivative=math.exp(div)*(r+a*r*r)/(r1+a*r1*r1),nfev=sol.nfev)
        except Exception as exc: record.update(status='failed',error=str(exc))
        finally:
            signal.setitimer(signal.ITIMER_PROF,0); signal.signal(signal.SIGPROF,prior)
            record['cpu_seconds']=time.process_time()-started
            with LEDGER.open('a') as f:f.write(json.dumps(record)+'\n')
        return record

def engines():
    R=S.Rational
    delta=R(-1,100); eps=R(-1,1000000)
    return {'shi':Engine('Shi conditioned',R(-10),5+delta,R(1),-25-9*delta+8*eps),
            'chen':Engine('Chen-Wang visualization',R(-3),R(99,100),R(2,9),R(-3))}

def main():
    p=argparse.ArgumentParser();p.add_argument('--family',choices=['shi','chen'],default='shi');p.add_argument('--r',type=float,nargs='+',required=True);p.add_argument('--lam',type=float,required=True);p.add_argument('--tag',default='manual');p.add_argument('--tight',action='store_true');args=p.parse_args()
    E=engines()[args.family]
    print('exact_etas',E.exact_etas,flush=True)
    for r in args.r:
        kw={'rtol':2e-13,'atol':2e-16,'max_step':.04} if args.tight else {}
        print(json.dumps(E.evaluate(r,args.lam,tag=args.tag,**kw)),flush=True)
if __name__=='__main__':main()
