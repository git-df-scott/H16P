"""Numerical finite equilibrium classification for the Q4 continuation fields."""
from pathlib import Path
import json
import numpy as np
from numpy.polynomial import Polynomial as Poly
P=Path(__file__).resolve().parent

def equilibria(r,eps,c):
    tau,u,w=eps*np.array(c);v=-eps
    p=Poly([-(2+r*r),2*r+2*u+w,1.]);q=Poly([r+u,-1-3*r*r+v,-r-u])
    elimination=Poly([1,tau])*p-Poly([tau,-1])*q
    result=[{'point':[0.,0.],'trace':2*tau,'determinant':1+tau*tau,'discriminant':-4.,'residual':0.}]
    for root in elimination.roots():
        if abs(root.imag)>1e-8:continue
        z=float(root.real)
        if abs(p(z))<1e-12:continue
        x=(z-tau)/p(z);y=z*x
        J=np.array([[tau-2*(2+r*r)*x+(2*r+2*u+w)*y,-1+(2*r+2*u+w)*x+2*y],[1+2*(r+u)*x+(-1-3*r*r+v)*y,tau+(-1-3*r*r+v)*x-2*(r+u)*y]])
        residual=max(abs(x*(tau-z)+x*x*p(z)),abs(x*(1+tau*z)+x*x*q(z)))
        tr=float(np.trace(J));det=float(np.linalg.det(J));disc=tr*tr-4*det
        result.append({'point':[x,y],'trace':tr,'determinant':det,'discriminant':disc,'residual':float(residual)})
    return {'equilibria':result,'unhandled_x0_y1_degeneracy':bool(abs(tau-r-u)<1e-10),'elimination_coefficients':elimination.coef.tolist()}

if __name__=='__main__':
    source=json.loads((P/'q4_finite_continuation.json').read_text());rows=[]
    for row in source['records']:
        eq=equilibria(row['r'],row['epsilon'],row['normalized_controls']);rows.append({'r':row['r'],'epsilon':row['epsilon'],**eq})
        print(row['r'],row['epsilon'],[(q['point'],q['discriminant']) for q in eq['equilibria'][1:]])
    (P/'q4_finite_equilibria.json').write_text(json.dumps({'scope':'NUM finite equilibrium enumeration; no rigorous root isolation','records':rows},indent=2)+'\n')
