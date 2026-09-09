"""Bounded third-focus branch continuation and finite-cycle coexistence screen."""
from pathlib import Path
from fractions import Fraction as F
import json
import numpy as np
from scipy.optimize import brentq
from reversible_third_focus_probe import model_at
from reversible_quintic_focus import coefficients
from reversible_finite_hopf import profile,count
P=Path(__file__).resolve().parent

def main():
    out={'scope':'NUM bounded branch tracking; no exact degeneracy/unfolding/isolated-cycle certificate','records':[],'failed_shapes':[]}
    for b in [F(1,6),F(1,3),F(1,2),F(2,3)]:
      previous=-1.1843004252363332
      for k in [F(1,20),F(1,10),F(1,5),F(2,5)]:
        cache={};evaluations=[]
        def fun(av):
            key=float(av)
            if key in cache:return cache[key]
            model=model_at(key,b,k);result=coefficients(model[0],rtol=3e-13);v=result['return_coefficients_c2_to_order'][3]
            cache[key]=v;evaluations.append({'a':key,'m':model[0]['m'],'c5':v});return v
        grid=sorted(set([previous]+[float(v) for v in np.linspace(-2.7,-.3,13)]));values=[];gates=[]
        for av in grid:
            try:values.append((av,fun(av)))
            except Exception as e:values.append((av,None));gates.append({'a':av,'reason':repr(e)})
        candidates=[]
        for (left,vl),(right,vr) in zip(values,values[1:]):
            if vl is None or vr is None or vl*vr>=0:continue
            try:
                root=brentq(fun,left,right,xtol=1e-10,rtol=1e-12);model=model_at(root,b,k)
                local=coefficients(model[0],rtol=3e-14,order=7);c7=local['return_coefficients_c2_to_order'][-1]
                candidates.append({'a':root,'c7':c7,'field':model[0],'local_order7':local,'quintic_bracket':[left,right]})
            except Exception as e:gates.append({'a_bracket':[left,right],'reason':repr(e)})
        eligible=[c for c in candidates if abs(c['c7'])>1e-5]
        if not eligible:
            out['failed_shapes'].append({'b':str(b),'k':str(k),'reason':'no numerically nonzero-seventh candidate found in bounded brackets','candidates':candidates,'evaluations':evaluations,'gates':gates})
            print(json.dumps({'b':str(b),'k':str(k),'status':'no eligible candidate','roots':[(c['a'],c['c7']) for c in candidates]}),flush=True)
        else:
            chosen=min(eligible,key=lambda c:abs(c['a']-previous));previous=chosen['a'];model=model_at(previous,b,k)
            upper=profile(model,np.geomspace(.03,1e7,49),1,2e-11);lower=profile(model,np.geomspace(.1,1e7,33),-1,2e-11)
            record={'b':str(b),'k':str(k),'selected':chosen,'other_candidates':[c for c in candidates if c is not chosen],'evaluations':evaluations,'gates':gates,'upper':upper,'lower':lower,'upper_brackets':count(upper),'lower_brackets':count(lower),'unresolved':sum('difference' not in row for row in upper+lower)}
            out['records'].append(record)
            print(json.dumps({'b':str(b),'k':str(k),'a':chosen['a'],'c7':chosen['c7'],'upper':record['upper_brackets'],'lower':record['lower_brackets'],'unresolved':record['unresolved']}),flush=True)
        (P/'reversible_third_focus_branch.json').write_text(json.dumps(out,indent=2)+'\n')
if __name__=='__main__':main()
