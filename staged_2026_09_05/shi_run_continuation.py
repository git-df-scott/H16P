#!/usr/bin/env python3
"""Two explicit trace paths, with adaptive root location and intermediate signs."""
import json
from pathlib import Path
from scipy.optimize import brentq
from shi_trace import engines,LEDGER
from shi_mp_verify import verify

OUT=Path(__file__).resolve().parent/'shi_continuation.json'
E=engines();result={'validated':False,'paths':{},'notes':['Rational conditioned Shi point differs from Shi/Galias certified hierarchy.','No interval integration or exhaustive cycle count.']}

def locate(engine,lam,bracket,label):
    samples=[]
    def fun(r):
        if engine is E['shi']:
            z=verify({'name':engine.name,'exact':engine.exact},r,lam,blocks=32,tag='root-'+label)
        else:z=engine.evaluate(r,lam,tag='root-'+label)
        if z['status']!='ok':raise ValueError(z['error'])
        v=float(z['raw_log_return']);samples.append([r,v,z['index']]);return v
    root=brentq(fun,*bracket,xtol=1e-9,rtol=1e-8,maxiter=16)
    detail=engine.evaluate(root,lam,rtol=2e-13,atol=2e-16,max_step=.04,tag='root-detail-'+label)
    return {'root':root,'initial_bracket':bracket,'root_samples':samples,'detail':detail}

specs={'shi':{'initial':-1e-14,'middle':-1e-16,'brackets':[(.00004,.00015),(.0006,.0015),(.020,.025)],'endpoint':[(.0006,.0015),(.020,.025)],'intermediate_radii':[.000004,.000015,.0015,.025]},
       'chen':{'initial':-2e-5,'middle':-2e-7,'brackets':[(.04,.08),(.15,.2),(.33,.333)],'endpoint':[(.15,.25),(.32,.333)],'intermediate_radii':[.004,.015,.2,.333]}}
for key,spec in specs.items():
    engine=E[key];path={'coefficients_exact':engine.exact,'focus_quantities_exact':engine.exact_etas,'initial_lambda':spec['initial'],'roots_initial':[],'roots_trace_zero':[],'intermediate_signs':[]}
    result['paths'][key]=path
    for i,br in enumerate(spec['brackets']):
        path['roots_initial'].append(locate(engine,spec['initial'],br,key+'-initial-'+str(i)))
        OUT.write_text(json.dumps(result,indent=2))
    for r in spec['intermediate_radii']:
        z=engine.evaluate(r,spec['middle'],tag=key+'-intermediate')
        path['intermediate_signs'].append(z)
    for i,br in enumerate(spec['endpoint']):
        path['roots_trace_zero'].append(locate(engine,0.,br,key+'-zero-'+str(i)))
        OUT.write_text(json.dumps(result,indent=2))
    path['origin_endpoint']=engine.evaluate(1e-6,0.,tag=key+'-origin-endpoint')
    OUT.write_text(json.dumps(result,indent=2))
print(json.dumps({k:{'initial':[r['root'] for r in v['roots_initial']],'zero':[r['root'] for r in v['roots_trace_zero']]} for k,v in result['paths'].items()},indent=2))
