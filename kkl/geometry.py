"""Exact rational finite-equilibrium gates; no floating-point sign decisions."""
from fractions import Fraction as F


class Interval:
    def __init__(self,lo,hi=None):
        self.lo=F(lo); self.hi=F(lo if hi is None else hi)
        assert self.lo<=self.hi
    def __add__(self,other):
        other=as_interval(other)
        return Interval(self.lo+other.lo,self.hi+other.hi)
    __radd__=__add__
    def __neg__(self): return Interval(-self.hi,-self.lo)
    def __sub__(self,other): return self+-as_interval(other)
    def __rsub__(self,other): return as_interval(other)+-self
    def __mul__(self,other):
        other=as_interval(other)
        values=[a*b for a in (self.lo,self.hi) for b in (other.lo,other.hi)]
        return Interval(min(values),max(values))
    __rmul__=__mul__
    def __truediv__(self,other):
        other=as_interval(other)
        if other.lo<=0<=other.hi: raise ValueError('division interval contains zero')
        return self*Interval(1/other.hi,1/other.lo)
    def strings(self): return [str(self.lo),str(self.hi)]


def as_interval(value): return value if isinstance(value,Interval) else Interval(value)


def gates(c,alpha,beta='0',require_box=True):
    c,alpha,beta=map(F,(c,alpha,beta))
    K=-alpha*(F(11,5)*c-1)-42
    if require_box:
        assert F(1,2)<=c<=F(3,2) and -200<=alpha<=-10 and K>=F(1,64)
    A,B,C,D=c-F(61,5),alpha-F(111,5)-beta,2*alpha-10-beta,alpha
    disc=18*A*B*C*D-4*B**3*D+B*B*C*C-4*A*C**3-27*A*A*D*D
    assert A!=0 and disc<0
    def poly(x): return ((A*x+B)*x+C)*x+D
    lo,hi=F(-2),F(-1)
    assert poly(hi)<0
    while poly(lo)<=0: lo*=2
    for _ in range(60):
        mid=(lo+hi)/2
        if poly(mid)>0: lo=mid
        else: hi=mid
    assert lo<hi<-1 and poly(lo)>0>poly(hi)
    x=Interval(lo,hi); y=-(x*x)/(1+x)
    j11=2*x+y; j12=1+x
    j21=-20*x+F(11,5)*y+alpha; j22=F(11,5)*x+2*c*y+beta
    trace=j11+j22; determinant=j11*j22-j12*j21
    focus_discriminant=trace*trace-4*determinant
    assert trace.hi<0 and determinant.lo>0 and focus_discriminant.hi<0
    return {'status':'EXACT_RATIONAL_INTERVAL_GATES','c':str(c),'alpha':str(alpha),'beta':str(beta),
            'K':str(K),'cubic_discriminant':str(disc),'remote_x':x.strings(),
            'trace':trace.strings(),'determinant':determinant.strings(),
            'focus_discriminant':focus_discriminant.strings(),
            'infinity_region':'c<241/250' if c<F(241,250) else 'requires separate itinerary audit'}
