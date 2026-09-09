"""Numerical fifth-order radial return coefficient at trace-zero equilibria.
Uses formal power-series angular flow; not a Lyapunov interval certificate.
"""
from pathlib import Path
import json
import numpy as np
from scipy.integrate import solve_ivp
P=Path(__file__).resolve().parent

def coefficients(exact,rtol=2e-12,order=5):
    from fractions import Fraction as F
    a,b,x,y,e2,det=[float(F(exact[q])) for q in ['a','b','x0','y0','epsilon2','determinant']]
    A=2*x;B=1-b+2*b*y+e2*x;omega=np.sqrt(det)
    def rhs(theta,cs):
        c,s=np.cos(theta),np.sin(theta);dy=(-omega*s-A*c)/B
        F2=a*c*c+e2*c*dy+b*dy*dy;G2=-2*c*dy
        f=F2/omega;g=-(A*F2+B*G2)/omega**2
        R=c*f+s*g;T=c*g-s*f
        radius=np.array([0.,1.,*cs]);power=np.array([1.]);result=np.zeros(order+1)
        for n in range(1,order+1):
            power=np.convolve(power,radius)[:order+1]
            if n>=2:result[:len(power)]+=R*(-T)**(n-2)*power
        return result[2:order+1]
    sol=solve_ivp(rhs,[0,2*np.pi],np.zeros(order-1),method='DOP853',rtol=rtol,atol=rtol*.01,max_step=.025)
    return {'status':'passed' if sol.success else 'failed','rtol':rtol,'order':order,'return_coefficients_c2_to_order':sol.y[:,-1].tolist(),'nfev':sol.nfev,'omega':omega,'B':B}

def main():
    source=json.loads((P/'reversible_bautin_precursor.json').read_text());out={'scope':'NUM radial Poincare coefficients; c5 at nearly cubic-degenerate rational fields, no exact l2 zero or three-cycle unfolding','records':[]}
    for row in source['records']:
        runs=[coefficients(row['field'],tol) for tol in [2e-12,3e-14]]
        record={'field':row['field'],'outer_brackets':row['upper_brackets'],'runs':runs};out['records'].append(record)
        print(json.dumps({'a':row['field']['a'],'b':row['field']['b'],'k':row['field']['k'],'m':float(__import__('fractions').Fraction(row['field']['m'])),'c5':[r['return_coefficients_c2_to_order'][-1] for r in runs],'c3':runs[-1]['return_coefficients_c2_to_order'][1],'outer_brackets':row['upper_brackets']}),flush=True)
    (P/'reversible_quintic_focus.json').write_text(json.dumps(out,indent=2)+'\n')
if __name__=='__main__':main()
