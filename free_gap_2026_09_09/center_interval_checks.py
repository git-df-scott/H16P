"""Interval sign checks for the explicit reversible parameter locus."""
from pathlib import Path
import json
import mpmath as mp
from mpmath import iv
mp.mp.dps=90;iv.dps=75
HERE=Path(__file__).resolve().parent
def coefficients(a,b):
    k=iv.sqrt(-a*(a+2)/(b*(2-a)))
    return (k*(a+b)*(b-a-2)/(4*b*(a+1)**2),-(a-1)*k*(1-b)/(a+1),2*(a-1)*(a+2)/((2-a)*k))
def low(v):return mp.make_mpf(v._mpi_[0])
def high(v):return mp.make_mpf(v._mpi_[1])
def enc(v):
    # Deliberately generous outward decimal endpoints for public summaries.
    def one(x,direction):
        if x==0:return '0'
        quantum=mp.power(10,mp.floor(mp.log10(abs(x)))-55)
        n=mp.floor(x/quantum) if direction<0 else mp.ceil(x/quantum)
        return mp.nstr(n*quantum,60)
    return [one(low(v),-1),one(high(v),1)]

a=iv.mpf(['-1.908','-1.907']);b=iv.mpf(['0.365','0.366']);e0,e1,e2=coefficients(a,b)
disc=e1*e1-4*e2*(a-1)*e0
assert low(disc)>0
equilibria=[]
for sign in [-1,1]:
    y=(-e1+sign*iv.sqrt(disc))/(2*e2);assert low(y)>0 or high(y)<0
    x=e0/(2*y)
    j11=2*a*x+e1+e2*y;j12=1-b+2*b*y+e2*x
    determinant=j11*(-2*x)+2*y*j12
    assert low(determinant)>0
    equilibria.append({'x':enc(x),'y':enc(y),'determinant':enc(determinant)})

fixed_b=iv.mpf('0.3656912644250558');target=iv.mpf('-0.12283156115176493')
left=iv.mpf('-1.907584173654');right=iv.mpf('-1.907584173652')
fl=coefficients(left,fixed_b)[0]-target;fr=coefficients(right,fixed_b)[0]-target
opposite=high(fl)<0<low(fr) or high(fr)<0<low(fl)
assert opposite
out={'scope':'Interval verification of explicit reversible-locus geometry, not cycle certification',
     'parameter_box':{'a':['-1.908','-1.907'],'b':['0.365','0.366']},
     'equilibrium_enclosures':equilibria,'discriminant':enc(disc),
     'fixed_b_e0_projection_a_bracket':['-1.907584173654','-1.907584173652'],
     'endpoint_residuals':[enc(fl),enc(fr)],'opposite_projection_signs':bool(opposite),
     'interval_digits':iv.dps,'trace_zero':'Exact on the explicit reversible symmetry axis; determinant intervals are strictly positive.'}
(HERE/'center_interval_checks.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out,indent=2))
