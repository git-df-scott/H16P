"""Bounded search near exact-trace-zero, cubic-degenerate Hopf fields.
Cubic roots are approximated by rational m: these are near-Bautin, not certified
Bautin points. Numerical complete half-passages; no global exclusion.
"""
from pathlib import Path
from fractions import Fraction as F
import json,time
import sympy as S
import numpy as np
from reversible_finite_hopf import field,profile,count
P=Path(__file__).resolve().parent

def roots(a,b,k):
    a,b,k=map(S.Rational,(a,b,k));m=S.symbols('m');x=-k/(2-a+b*k*k);y=S.Rational(1,2)+k*x
    e=-2*(a-1)*x*y*m;A=2*x;B=1-b+2*b*y+e*x
    core=S.factor(A*A*a*e-A*B*a*a-A*B*a+2*A*B-2*A*a*b*y+2*A*b*y-A*e*e*y+B*a*e*y+2*b*e*y*y)
    poly=S.Poly(core,m);values=[]
    for root in S.solve(poly.as_expr(),m):
        z=S.N(root,40)
        if z.is_real:
            q=F(str(z));det=-A*A+2*B*y
            if det.subs(m,S.Rational(q))>0 and y>0:
                values.append({'m':q,'cubic_core_residual':str(S.N(core.subs(m,S.Rational(q)),20)),'determinant_approx':str(S.N(det.subs(m,S.Rational(q)),20))})
    return str(poly.as_expr()),values

def main():
    out={'scope':'NUM near-cubic-degenerate Hopf precursor screen; target two outer upper cycles plus one lower before a two-small-cycle unfolding; no l2/nondegeneracy certificate','records':[],'gates':[]};start=time.perf_counter()
    for a in [F(-3,4),F(-5,4),F(-7,4),F(-5,2)]:
      for b in [F(1,3),F(1),F(5,3)]:
       for k in [F(1,10),F(2,5)]:
        polynomial,rs=roots(a,b,k)
        if not rs:out['gates'].append({'a':str(a),'b':str(b),'k':str(k),'status':'no positive-determinant real cubic root','polynomial':polynomial})
        for root in rs:
          try:model=field(k,root['m'],a,b)
          except (AssertionError,ArithmeticError) as e:
            out['gates'].append({'a':str(a),'b':str(b),'k':str(k),'m':str(root['m']),'status':'field gate failed','reason':repr(e)});continue
          if not(model[2][1]>0>model[3][1]):
            out['gates'].append({'field':model[0],'status':'wrong two-focus geometry'});continue
          upper=profile(model,np.geomspace(.03,1e6,33),1,2e-11);ub=count(upper)
          lower=profile(model,np.geomspace(.1,1e7,33),-1,2e-11) if len(ub)>=2 else []
          row={'field':model[0],'equilibria_numerical':model[1],'cubic_polynomial_in_m':polynomial,'cubic_core_residual':root['cubic_core_residual'],'upper':upper,'lower':lower,'upper_brackets':ub,'lower_brackets':count(lower),'unresolved':sum('difference' not in q for q in upper+lower)}
          out['records'].append(row);out['wall_seconds']=time.perf_counter()-start
          (P/'reversible_bautin_precursor.json').write_text(json.dumps(out,indent=2)+'\n')
          print(json.dumps({'a':str(a),'b':str(b),'k':str(k),'m':float(root['m']),'upper':ub,'lower':row['lower_brackets'],'unresolved':row['unresolved']}),flush=True)
    (P/'reversible_bautin_precursor.json').write_text(json.dumps(out,indent=2)+'\n')
if __name__=='__main__':main()
