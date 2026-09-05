"""Polish already bracketed two-cycle returns; preserve the same field exactly."""
import json
from fractions import Fraction as Q
import mpmath as m
from budget import call,HERE,used
m.mp.dps=50
def st(x):return m.nstr(x,40)
def refine(left,right,engine,keys):
 history=[];lo=m.log(m.mpf(left['r']));hi=m.log(m.mpf(right['r']));a=min([left,right],key=lambda x:abs(m.mpf(x['L'])))
 for i in range(7):
  if used()>=3338:break
  if abs(m.mpf(a['L']))<m.mpf('2e-12'):break
  z=m.log(m.mpf(a['r']))-m.mpf(a['L'])/m.mpf(a['L_z'])
  if not lo<z<hi:z=(lo+hi)/2
  req={k:a[k] for k in keys};req['r']=st(m.exp(z));req['tol']='2e-26'
  b=call(req,'polish complete-return root at fixed replay field',engine=engine);history.append(b)
  if b['status']!='NUMERICAL_ONLY':return dict(status='UNRESOLVED',left=left,right=right,history=history)
  if m.mpf(b['L'])*m.mpf(left['L'])>0:lo=z;left=b
  else:hi=z;right=b
  a=b
 return dict(status='NUMERICAL_ROOT_WITH_SIGN_BRACKET',left=left,right=right,root=a,history=history,certified=False)
rows=[]
raw=json.loads((HERE/'full_return_reproduction.json').read_text())
for row in raw['rows']:
 c,K=Q(row['field']['c']),Q(row['field']['K']);M=5*(K+42)/(11*c-5)
 roots=[]
 for b in row['replays']:
  left,right=[x['full'] for x in b['checks'][:2]]
  if all(x['status']=='NUMERICAL_ONLY' for x in [left,right]):roots.append(refine(left,right,'angular_quad.py',['c','K']))
 rows.append(dict(source=row['source'],exact_coefficients=dict(P=['0','0','1','1','1','0'],Q=['0',str(-M),'0','-10','11/5',str(c)],basis=['1','x','y','x^2','xy','y^2']),section='y=0,x>0; clockwise full return',cycles=roots))
raw=json.loads((HERE/'large_m_full_reproduction.json').read_text())
roots=[]
for b in raw['rows']:
 left,right=[x['result'] for x in b['checks'][:2]]
 if all(x['status']=='NUMERICAL_ONLY' for x in [left,right]):roots.append(refine(left,right,'angular_m_quad.py',['c','m']))
if roots:
 rr=roots[0]['root'];c,M=Q(rr['c']),Q(rr['m'])
 rows.append(dict(source='events_logm.json',exact_coefficients=dict(P=['0','0','1','1','1','0'],Q=['0',str(-M),'0','-10','11/5',str(c)],basis=['1','x','y','x^2','xy','y^2']),section='y=0,x>0; clockwise full return',cycles=roots))
(HERE/'pair_replay_claims.json').write_text(json.dumps(dict(status='NUMERICAL_TWO_CYCLE_FIELDS_ONLY',rows=rows,independent_five_cycle_trigger=False),indent=2))
