"""Continue through the removable (c,K) chart singularity in actual (c,m)."""
import json
import mpmath as m
import continue_half as base
from budget import call,HERE
m.mp.dps=45
def ev(z,c,M,purpose):
 a=call(dict(r=base.st(m.exp(z)),c=base.st(c),m=base.st(M),tol='2e-18'),purpose,engine='half_m.py')
 if a['status']=='NUMERICAL_TWO_HALF_PASSAGES':a['F_K']=a['F_m'];a['G_K']=a['G_m']
 return a
base.ev=ev
# Profile's third argument is now actual m; keep metadata unambiguous afterward.
s=__import__('inspect').getsource(base.profile).replace("K=m.mpf(a['K'])","K=m.mpf(a['m'])")
exec(s,base.__dict__)
def correct(z,c,M):
 hist=[]
 for i in range(12):
  a=ev(z,c,M,'actual-coefficient fold continuation through c=5/11');hist.append(a)
  if a['status']!='NUMERICAL_TWO_HALF_PASSAGES':return None,hist
  F=m.matrix([a['F'],a['G']]);J=m.matrix([[a['F_z'],a['F_m']],[a['G_z'],a['G_m']]])
  if abs(F[0])<m.mpf('1e-14') and abs(F[1])<m.mpf('1e-12') and m.exp(z)>1:
   a['jacobian_zm_determinant']=base.st(m.det(J));return a,hist
  d=m.lu_solve(J,F)
  if abs(d[0])>.5:d*=m.mpf('.5')/abs(d[0])
  z-=d[0];M-=d[1]
  if M<=0:return None,hist
 return None,hist
out=HERE/'events_m.json';events=json.loads(out.read_text()) if out.exists() else []
if events:a=next(e['fold'] for e in reversed(events) if e['status']=='ACCEPTED')
else:
 seed=next(e['fold'] for e in reversed(json.loads((HERE/'events_negative.json').read_text())) if e['status']=='ACCEPTED')
 a=ev(m.log(m.mpf(seed['r'])),m.mpf(seed['c']),-m.mpf(seed['alpha']),'actual-coefficient chart transfer')
z=m.log(m.mpf(a['r']));c=m.mpf(a['c']);M=m.mpf(a['m'])
for target in ['.46','5/11','.44','.42','.4','.35','.3','.25','.2','.15','.1','.05']:
 cn=m.mpf(target)
 if cn>=c:continue
 J=m.matrix([[a['F_z'],a['F_m']],[a['G_z'],a['G_m']]])
 t=m.lu_solve(J,-m.matrix([a['F_c'],a['G_c']]))
 b,h=correct(z+(cn-c)*t[0],cn,M+(cn-c)*t[1])
 if b is None:
  events.append(dict(status='CORRECTOR_UNRESOLVED',target_c=target,history=h));out.write_text(json.dumps(events,indent=2));break
 a=b;z=m.log(m.mpf(a['r']));c=cn;M=m.mpf(a['m']);p=base.profile(a)
 p['m']=p.pop('K');p['K_at_pair_field']=base.st(M*(11*m.mpf(p['c'])-5)/5-42)
 event=dict(status='ACCEPTED',fold=a,history=h,pair_profile=p);events.append(event);out.write_text(json.dumps(events,indent=2))
 print('POINT',a['r'],a['c'],a['m'],'K',a['K'],'BRACKETS',len(p['root_sign_brackets']),flush=True)
 if len(p['root_sign_brackets'])>=3:
  (HERE/'K1_CANDIDATE_M.json').write_text(json.dumps(event,indent=2));break
