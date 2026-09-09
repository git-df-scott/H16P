"""Exact original-coordinate Q4 preflight at the council r=1 base.

No orbit integration; the invariant quartic is not by itself a proof of
the annulus boundary itinerary. All fields have degree at most two.
"""
from pathlib import Path
import json
import sympy as s

X,Y,tau,u,v,w=s.symbols('X Y tau u v w')
P=-Y-3*X**2+2*X*Y+Y**2
Q=X+X**2-4*X*Y-Y**2
D=1-4*Y+2*(X+Y)**2
N=6*Y*(1-X-Y)-1+2*(X+Y)**3
L=lambda f:s.expand(s.diff(f,X)*P+s.diff(f,Y)*Q)
assert s.expand(L(D)+4*X*D)==0
assert s.expand(L(N)+6*X*N)==0
assert s.expand(s.diff(N,X)*D-s.Rational(3,2)*N*s.diff(D,X)-6*Q)==0
assert s.expand(s.diff(N,Y)*D-s.Rational(3,2)*N*s.diff(D,Y)+6*P)==0
C=s.expand(2*N**2-D**3)
assert s.Poly(C,X,Y).total_degree()==4
assert s.expand(L(C)+12*X*C)==0
top=sum(co*X**i*Y**j for (i,j),co in s.Poly(C,X,Y).terms() if i+j==4)
assert s.expand(top+12*(X+Y)**2*(X**2-2*X*Y-Y**2))==0
dP=tau*X+(2*u+w)*X*Y
dQ=tau*Y+u*(X**2-Y**2)+v*X*Y
# div(6 D^(-5/2) deltaF) = 3 D^(-7/2) * density_numerator.
density_numerator=s.expand(2*D*(s.diff(dP,X)+s.diff(dQ,Y))-5*(s.diff(D,X)*dP+s.diff(D,Y)*dQ))
record={
 'status':'all exact assertions passed',
 'base':{'P':str(P),'Q':str(Q)},
 'D':str(D),'N':str(N),
 'first_integral':'h=N/D^(3/2), on D>0; h(0,0)=-1',
 'gradient':'grad(h)=6 D^(-5/2) (Q,-P)',
 'invariant_quartic':str(C),
 'quartic_infinity_factors':str(s.factor(top)),
 'candidate_boundary_level':'h=-1/sqrt(2); requires D>0 and N<0 in addition to C=0',
 'normal_controls':['tau','u','v','w'],
 'area_density_numerators':{str(a):str(s.factor(density_numerator.coeff(a))) for a in [tau,u,v,w]},
 'area_density_factor':'3 D^(-7/2)',
 'limitations':['The appropriate real quartic branch and global saddle itinerary have not been proved.',
                'No compact Melnikov roots or simultaneous boundary cycles have been certified.',
                'First-order area integrals apply to compact unperturbed ovals; endpoint limits need separate estimates.'],
 'ode_evaluations':0}
Path(__file__).with_suffix('.json').write_text(json.dumps(record,indent=2)+'\n')
print(json.dumps(record,indent=2))
