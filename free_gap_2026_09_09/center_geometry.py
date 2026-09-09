"""An exact affine-reversibility locus in the searched five-parameter family.

The derivation is in CENTER_GEOMETRY.md. Formula evaluation and proximity
measurements are numerical; rational reflection checks use Fraction exactly.
"""
import json
from fractions import Fraction as F
from pathlib import Path
import numpy as np
import mpmath as mp
HERE=Path(__file__).resolve().parent

def center_coefficients(a,b):
    k=np.sqrt(-a*(a+2)/(b*(2-a)))
    e1=-(a-1)*k*(1-b)/(a+1)
    e2=2*(a-1)*(a+2)/((2-a)*k)
    e0=k*(a+b)*(b-a-2)/(4*b*(a+1)**2)
    return np.array([e0,e1,e2])

def center_jacobian(a,b):
    return np.column_stack([np.imag(center_coefficients(a+1e-25j,b))/1e-25,
                            np.imag(center_coefficients(a,b+1e-25j))/1e-25])

def rational_check(a,k):
    b=-a*(a+2)/((2-a)*k*k)
    e1=-(a-1)*k*(1-b)/(a+1)
    e2=2*(a-1)*(a+2)/((2-a)*k)
    e0=k*(a+b)*(b-a-2)/(4*b*(a+1)**2)
    A=2*(a-1);den=A+e2*k
    def field(z):
        x,y=z
        return [(b-2)/4+(1-b)*y+a*x*x+b*y*y+e1*x+e2*x*y,e0-2*x*y]
    def reflect(z):
        x,y=z;t=2*(A*x+e2*y+e1)/den
        return [x-t,y-k*t]
    def linear(z):
        x,y=z;t=2*(A*x+e2*y)/den
        return [x-t,y-k*t]
    # Six unisolvent points determine any bivariate quadratic polynomial.
    points=[(0,0),(1,0),(-1,0),(0,1),(0,-1),(1,1)]
    for z in points:
        z=list(map(F,z));lhs=field(reflect(z));rhs=[-x for x in linear(field(z))]
        assert lhs==rhs and reflect(reflect(z))==z
    return dict(a=str(a),b=str(b),k=str(k),e0=str(e0),e1=str(e1),e2=str(e2),exact_reflection_identity=True)

if __name__=='__main__':
    mp.mp.dps=75
    s=json.loads((HERE/'strong_corrected_result.json').read_text())['states'][-1]
    actual=np.array(s['q'])*np.array([1,1,1e-4,1e-4,1e-4]);a,b,e0,e1,e2=map(lambda v:mp.mpf(repr(float(v))),actual)
    def coefficients(aa):
        k=mp.sqrt(-aa*(aa+2)/(b*(2-aa)))
        return [k*(aa+b)*(b-aa-2)/(4*b*(aa+1)**2),-(aa-1)*k*(1-b)/(aa+1),2*(aa-1)*(aa+2)/((2-aa)*k)]
    ac=mp.findroot(lambda aa:coefficients(aa)[0]-e0,(a-mp.mpf('1e-4'),a+mp.mpf('1e-4')))
    cc=coefficients(ac)
    out={'evidence':'Exact formula/rational reflection checks; numerical proximity projection',
         'checks':[rational_check(F(-7,4),F(1,3)),rational_check(F(-19,10),F(3,8)),rational_check(F(-3,2),F(1,2))],
         'actual_coefficients':[str(x) for x in [a,b,e0,e1,e2]],
         'fixed_b_e0_center_projection':[mp.nstr(x,65) for x in [ac,b,*cc]],
         'projection_difference':[mp.nstr(x-y,45) for x,y in zip([ac,b,*cc],[a,b,e0,e1,e2])]}
    (HERE/'center_geometry.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(out,indent=2))
