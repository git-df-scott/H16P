"""Numerical neutral graphic coefficient from bounded half-map convergence."""
import json
import mpmath as m
from budget import call,HERE
m.mp.dps=45
a=json.loads((HERE/'infinity_tolerance_control.json').read_text())[-1];K=m.mpf(a['K']);hist=[]
for i in range(2):
 K-=m.mpf(a['F'])/m.mpf(a['F_K'])
 a=call(dict(r='1e17',c='8/5',K=m.nstr(K,38),tol='2e-25'),'refined neutral graphic connection and coefficient',engine='half_quad.py');hist.append(a)
 if a['status']!='NUMERICAL_TWO_HALF_PASSAGES':break
out=dict(history=hist,exact_or_interval=False)
if a['status']=='NUMERICAL_TWO_HALF_PASSAGES':
 G=m.mpf(a['G']);out.update(K_at_finite_matching_radius=a['K'],G_infinity_approx=a['G'],C_approx=m.nstr(m.exp(5*G/6),35),conditional_delta_log_r_limit=m.nstr(-m.sqrt(159)*G/12,35),finite_radius='1e17',c='8/5')
(HERE/'graphic_coefficient.json').write_text(json.dumps(out,indent=2))
