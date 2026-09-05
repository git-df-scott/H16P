"""Binary128 center, c=1 crossing, and infinity connection controls.
Not rigorous enclosures; all integrations are charged through budget.py.
"""
import json
import mpmath as m
from budget import call,HERE
m.mp.dps=45
def st(x):return m.nstr(x,38)
def ev(z,c,K,purpose,tol='2e-28'):
 return call(dict(r=st(m.exp(z)),c=st(c),K=st(K),tol=tol),purpose,engine='half_quad.py')
rows=[]
# Approximate connection by increasing the two-sided starting radius.
K=m.mpf('7.184994696941435');c=m.mpf('1.6')
for r in ['1e6','1e10','1e14','1e17']:
 z=m.log(m.mpf(r));hist=[]
 for i in range(5):
  a=ev(z,c,K,'neutral infinity connection finite-radius convergence');hist.append(a)
  if a['status']!='NUMERICAL_TWO_HALF_PASSAGES':break
  if abs(m.mpf(a['F']))<m.mpf('1e-27'):break
  K-=m.mpf(a['F'])/m.mpf(a['F_K'])
 rows.append(dict(r=r,history=hist,result=a))
(HERE/'infinity_binary128.json').write_text(json.dumps(dict(rows=rows,exact_connection=False),indent=2))
# c=1 crossing: solve in z,K; horizontal section differs from curved r.
z=m.log(4);K=m.mpf('.22');c=m.mpf(1);hist=[]
for i in range(12):
 a=ev(z,c,K,'finite fold crossing c=1');hist.append(a)
 if a['status']!='NUMERICAL_TWO_HALF_PASSAGES':break
 if max(abs(m.mpf(a[k])) for k in ('F','G'))<m.mpf('1e-22'):break
 J=m.matrix([[a['F_z'],a['F_K']],[a['G_z'],a['G_K']]])
 d=m.lu_solve(J,m.matrix([a['F'],a['G']]))
 if max(abs(x) for x in d)>1:d/=max(abs(x) for x in d)
 z-=d[0];K-=d[1]
(HERE/'crossing_c1.json').write_text(json.dumps(dict(history=hist,result=a),indent=2))
# K→0: solve F/K=G/K=0 (same zeros at nonzero K), with c,z correction.
# Transfer from known curved fold by Newton on the horizontal section.
cstar=m.findroot(lambda c:305+634*c-11*c*c-1000*c**3,.969)
z=m.log(4);c=cstar+m.mpf('0.137109611')/8192
center=[]
for kk in ['1/8192','1e-5','1e-7','1e-9']:
 K=m.mpf(kk);c=cstar+m.mpf('.137109611')*K;hist=[]
 for i in range(12):
  a=ev(z,c,K,'center-limit finite fold continuation');hist.append(a)
  if a['status']!='NUMERICAL_TWO_HALF_PASSAGES':break
  F=m.matrix([a['F'],a['G']]);J=m.matrix([[a['F_z'],a['F_c']],[a['G_z'],a['G_c']]])
  if max(abs(x)/K for x in F)<m.mpf('1e-18'):break
  d=m.lu_solve(J,F)
  if abs(d[0])>1:d/=abs(d[0])
  z-=d[0];c-=d[1]
 center.append(dict(K=kk,history=hist,result=a,secant_dc_dK=st((m.mpf(a.get('c',c))-cstar)/K)))
(HERE/'center_binary128.json').write_text(json.dumps(dict(cstar=st(cstar),rows=center),indent=2))
