"""Locate numerical quintic zeros on cubic-degeneracy branch; test coexistence."""
from pathlib import Path
from fractions import Fraction as F
import json
import numpy as np
from scipy.optimize import brentq
from reversible_bautin_precursor import roots
from reversible_finite_hopf import field,profile,count
from reversible_quintic_focus import coefficients
P=Path(__file__).resolve().parent

def model_at(av,b,k):
    a=F(str(float(av)));poly,rs=roots(a,b,k);models=[]
    for q in rs:
        try:models.append(field(k,q['m'],a,b))
        except (ArithmeticError,AssertionError):pass
    if not models:raise ArithmeticError('no admitted cubic-degenerate field')
    return min(models,key=lambda q:abs(float(F(q[0]['m']))))

def main():
    out={'scope':'NUM c5 zeros on nearly cubic-degenerate rational fields, c7 diagnostic and finite outer passage screen; not an exact third-order weak focus or generic unfolding','records':[],'failures':[]}
    for b,k,left,right in [(F(1,3),F(1,10),-1.25,-.75),(F(1,3),F(1,10),-2.5,-1.75),(F(1,3),F(2,5),-2.5,-1.75),(F(1),F(2,5),-1.75,-1.25),(F(1),F(2,5),-2.5,-1.75)]:
        evaluations=[]
        def fun(a):
            model=model_at(a,b,k);result=coefficients(model[0],rtol=3e-13);v=result['return_coefficients_c2_to_order'][-1]
            evaluations.append({'a':a,'m':model[0]['m'],'c5':v});return v
        try:
            root=brentq(fun,left,right,xtol=2e-10,rtol=1e-12)
            model=model_at(root,b,k);local=[coefficients(model[0],rtol=tol,order=7) for tol in [2e-12,3e-14]]
            upper=profile(model,np.geomspace(.03,1e7,49),1,2e-11);lower=profile(model,np.geomspace(.1,1e7,33),-1,2e-11)
            record={'b':str(b),'k':str(k),'a_bracket':[left,right],'field':model[0],'equilibria_numerical':model[1],'evaluations':evaluations,'local_order7':local,'upper':upper,'lower':lower,'upper_brackets':count(upper),'lower_brackets':count(lower),'unresolved':sum('difference' not in x for x in upper+lower)}
            out['records'].append(record)
            print(json.dumps({'a':root,'b':str(b),'k':str(k),'c2_through_c7':local[-1]['return_coefficients_c2_to_order'],'upper':record['upper_brackets'],'lower':record['lower_brackets'],'unresolved':record['unresolved']}),flush=True)
        except Exception as e:
            out['failures'].append({'b':str(b),'k':str(k),'a_bracket':[left,right],'reason':repr(e),'evaluations':evaluations});print('FAIL',repr(e),flush=True)
        (P/'reversible_third_focus_probe.json').write_text(json.dumps(out,indent=2)+'\n')
if __name__=='__main__':main()
