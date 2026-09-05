#!/usr/bin/env python3
"""Exact replay of Strike 6 constants, interval anchor bounds, and identities.

The sign theorems additionally use the analytic proofs in the accompanying
notes. No sampled parameter sign is used as a proof of a global inequality.
"""
import os
for key in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ[key] = '1'
import hashlib
import json
from fractions import Fraction as Q
from pathlib import Path
import resource
import time
import sympy as S


class Interval:
    def __init__(self, lo, hi=None):
        if isinstance(lo, Interval): self.lo, self.hi = lo.lo, lo.hi
        else: self.lo, self.hi = Q(lo), Q(lo if hi is None else hi)
        assert self.lo <= self.hi
    def __add__(self, other):
        other=Interval(other); return Interval(self.lo+other.lo,self.hi+other.hi)
    __radd__=__add__
    def __neg__(self): return Interval(-self.hi,-self.lo)
    def __sub__(self,other):return self+-Interval(other)
    def __rsub__(self,other):return Interval(other)+-self
    def __mul__(self,other):
        other=Interval(other)
        ends=[x*y for x in (self.lo,self.hi) for y in (other.lo,other.hi)]
        return Interval(min(ends),max(ends))
    __rmul__=__mul__
    def __truediv__(self,other):
        other=Interval(other)
        assert other.lo*other.hi>0, 'division interval meets zero'
        return self*Interval(1/other.hi,1/other.lo)
    def decimal(self,places=24):
        scale=10**places
        left=self.lo.numerator*scale//self.lo.denominator
        right=-((-self.hi.numerator*scale)//self.hi.denominator)
        def fmt(v):
            sign='-' if v<0 else '';v=abs(v)
            return f'{sign}{v//scale}.{v%scale:0{places}d}'
        out=[fmt(left),fmt(right)]
        assert Q(out[0])<=self.lo<=self.hi<=Q(out[1])
        return out


def det3(m):
    return sum((m[0][i]*(m[1][(i+1)%3]*m[2][(i+2)%3]
                        -m[1][(i+2)%3]*m[2][(i+1)%3]) for i in range(3)),Interval(0))


def moments(t,N=128):
    """One-sided rational enclosures of K0,...,K3,F,D=tFM.

    Both f_n and d_n=6*n*f_n/(6*n-1), n>=1, decrease. Each positive
    omitted series is bounded by its first omitted term/(1-t).
    """
    t=Q(t); assert 0<t<1
    vals=[Q(0)]*6
    f=Q(1);tp=Q(1)
    for n in range(N+1):
        d=Q(6*n,6*n-1)*f if n else Q(0)
        terms=[f*tp*t*t/(n+2),f*tp*t**3/(n+3),
               d*tp*t/(n+1),d*tp*t*t/(n+2),f*tp,d*tp]
        vals=[x+y for x,y in zip(vals,terms)]
        f*=Q((6*n+1)*(6*n+5),36*(n+1)**2);tp*=t
    n=N+1;d=Q(6*n,6*n-1)*f
    tails=[f*tp*t*t/(n+2),f*tp*t**3/(n+3),
           d*tp*t/(n+1),d*tp*t*t/(n+2),f*tp,d*tp]
    return [Interval(x,x+tail/(1-t)) for x,tail in zip(vals,tails)]


def confluent_bound(t,threshold,N=128):
    k0,k1,k2,k3,f,d=moments(t,N)
    matrix=[[k0,k1,-k2],[t*f,t*t*f,-d],
            [Interval(9061),Interval(6289),Interval(-2431)]]
    rhs=[-k3,-t*d,Interval(-1819)]
    denom=det3(matrix)
    assert denom.hi<0 or denom.lo>0
    coeff=[]
    for j in range(3):
        replace=[[rhs[i] if k==j else matrix[i][k] for k in range(3)] for i in range(3)]
        coeff.append(det3(replace)/denom)
    alpha,beta,eta=coeff
    yc=Q(9,3080)*(alpha+Q(144,221)*beta-Q(11,6)*eta+Q(204,221))
    assert yc.hi<0 and eta.lo>0
    ratio=eta/(-192*yc)
    assert ratio.lo>threshold
    return {'r':str(t),'last_included_index':N,
            'alpha_beta_eta':[x.decimal() for x in coeff],
            'Y0':yc.decimal(),'q_eta_over_minus_192_Y0':ratio.decimal(),
            'q_strictly_exceeds':str(threshold),
            'margin':(ratio-threshold).decimal()}


def identities():
    A,B,eta=S.symbols('A B eta')
    seed={A:S.Rational(11843,9623),B:-S.Rational(833,9623),eta:S.Rational(13320,9623)}
    for e in [A-1-eta/6,B+S.Rational(1,6)-25*eta/432,
              9061*A+6289*B-2431*eta-7242]:assert S.factor(e.subs(seed))==0
    yc=3*(1326*A+864*B-2431*eta-102)/1361360
    assert yc.subs(seed)==-S.Rational(81,19246)
    assert S.factor(-seed[eta]/(192*yc.subs(seed)))==S.Rational(185,108)
    # Use d=x^3 so all powers of d^(1/3) become integer powers.
    x,a,q=S.symbols('x a q',positive=True)
    t=1-x**3
    dt=lambda f:-S.diff(f,x)/(3*x*x)
    aa=S.Rational(13,4)-S.Rational(3,2)*q
    bb=aa+1
    v=aa*x**-4-bb*x**-2
    op=lambda f:(1-a*t)*(1-t)*dt(dt(f))-(1+5*a-6*a*t)*dt(f)/2+8*a*f/9
    expected=(1-a)*(22*aa-7*bb*x*x)/(9*x**7)
    assert S.factor(op(v)-expected)==0
    assert v.subs(x,1)==-1
    assert S.factor(dt(v).subs(x,1)-(S.Rational(3,2)-q))==0
    assert S.factor((22*aa-7*bb).subs(q,S.Rational(167,90)))==0
    assert v.subs({q:S.Rational(185,108),x:S.Rational(7,11)})==0
    assert S.factor((aa/bb).subs(q,S.Rational(167,90)))==S.Rational(7,22)
    assert Q(988,1331)>Q(7,10)
    assert Q(7,22)**3 < Q(1,5)**2 # root 1-(7/22)^(3/2) exceeds 4/5
    # Independently verify the full gauge conversion, not just its residual.
    tt,al=S.symbols('tt al',positive=True)
    vv=S.Function('v')(tt)
    yy=(1-al*tt)**S.Rational(3,2)*vv
    old=(1-al*tt)*(1-tt)*S.diff(yy,tt,2)-(1-al)*S.diff(yy,tt)/2+5*al*yy/36
    new=(1-al*tt)*(1-tt)*S.diff(vv,tt,2)-(1+5*al-6*al*tt)*S.diff(vv,tt)/2+8*al*vv/9
    assert S.simplify(old-(1-al*tt)**S.Rational(3,2)*new)==0
    u=S.symbols('u',nonnegative=True)
    kap=(1+u)*(1+4*u)**2; dd=8*u*u+10*u+5
    ca=3-S.Rational(3,2)/kap-5*(4*u+3)*(2*u+1)/(6*dd)
    polynomial=1024*u**5+2816*u**4+3936*u**3+2584*u**2+440*u-135
    assert S.factor(ca-S.Rational(19,10)-polynomial/(30*kap*dd))==0
    assert all(c>0 for c in S.Poly(S.diff(polynomial,u),u).all_coeffs())
    ulo=Q(1473975727428,10**13);uhi=Q(1473975727429,10**13)
    assert polynomial.subs(u,S.Rational(ulo.numerator,ulo.denominator))<0
    assert polynomial.subs(u,S.Rational(uhi.numerator,uhi.denominator))>0
    kval=lambda v:(1+v)*(1+4*v)**2
    # The restricted determinant orientation uses columns (e0,e1,V,E).
    k0,k1,l0,l1,yv,ye,zv,ze=S.symbols('k0 k1 l0 l1 yv ye zv ze')
    m=S.Matrix([[k0,k1,0,0],[l0,l1,0,0],[1,2,yv,ye],[3,4,zv,ze]])
    assert S.factor(m.det()-(k0*l1-k1*l0)*(yv*ze-ye*zv))==0
    return {'corner_A_B_eta':[str(seed[A]),str(seed[B]),str(seed[eta])],
            'corner_Y0':'-81/19246','q_global_infimum':'185/108',
            'first_bootstrap_r':'988/1331',
            'final_r_lower_bound':'1-(7/22)^(3/2)',
            'lift_polynomial_ascending':[-135,440,2584,3936,2816,1024],
            'u_root_enclosure':[str(ulo),str(uhi)],
            'kappa_root_enclosure':Interval(kval(ulo),kval(uhi)).decimal(),
            'identities':'PASS'}


def main():
    resource.setrlimit(resource.RLIMIT_CPU,(30,30));os.nice(10)
    started=time.process_time()
    result={'status':'EXACT_IDENTITIES_AND_RATIONAL_INTERVAL_CERTIFICATES',
            'identities':identities(),
            'confluent_anchor_certificates':[confluent_bound(Q(7,10),Q(167,90)),
                                             confluent_bound(Q(4,5),Q(19,10))],
            'scope':'The global sign implications use the analytic notes; this script does not certify the remaining determinant faces.',
            'cpu_seconds':time.process_time()-started,
            'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    Path(__file__).with_name('exact_checks.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))
if __name__=='__main__':main()
