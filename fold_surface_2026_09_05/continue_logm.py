"""Continue the center-selected negative sheet in log of actual coefficient m."""
import json
import mpmath as m
import continue_half as base
from budget import call,HERE,used
m.mp.dps=45
def ev(z,c,M,purpose):
 return call(dict(r=base.st(m.exp(z)),c=base.st(c),m=base.st(M),tol='2e-25'),purpose,engine='half_m_quad.py')
base.ev=ev
s=__import__('inspect').getsource(base.profile).replace("K=m.mpf(a['K'])","K=m.mpf(a['m'])");exec(s,base.__dict__)
def correct(z,c,M):
 hist=[]
 for i in range(10):
  a=ev(z,c,M,'log-m fold correction');hist.append(a)
  if a['status']!='NUMERICAL_TWO_HALF_PASSAGES':return None,hist
  F=m.matrix([a['F'],a['G']]);J=m.matrix([[a['F_z'],a['F_c']],[a['G_z'],a['G_c']]])
  if abs(F[0])<m.mpf('1e-22') and abs(F[1])<m.mpf('1e-18'):
   a['jacobian_zc_determinant']=base.st(m.det(J));return a,hist
  d=m.lu_solve(J,F)
  if abs(d[0])>.7:d*=m.mpf('.7')/abs(d[0])
  z-=d[0];c-=d[1]
 return None,hist
out=HERE/'events_logm.json';events=json.loads(out.read_text()) if out.exists() else []
a=next(e['fold'] for e in reversed(events if events else json.loads((HERE/'events_m.json').read_text())) if e['status']=='ACCEPTED')
z=m.log(m.mpf(a['r']));c=m.mpf(a['c']);M=m.mpf(a['m']);step=m.mpf('1.5')
for i in range(10):
 if used()>3100:break # reserve at least240 calls for complete-return controls and handoff
 J=m.matrix([[a['F_z'],a['F_c']],[a['G_z'],a['G_c']]])
 t=m.lu_solve(J,-m.matrix([a['F_m'],a['G_m']])*M);Mn=M*m.exp(step)
 b,h=correct(z+step*t[0],c+step*t[1],Mn)
 if b is None:events.append(dict(status='CORRECTOR_UNRESOLVED',target_m=base.st(Mn),history=h));out.write_text(json.dumps(events,indent=2));break
 a=b;z=m.log(m.mpf(a['r']));c=m.mpf(a['c']);M=Mn;p=base.profile(a)
 p['m']=p.pop('K');p['K_at_pair_field']=base.st(M*(11*m.mpf(p['c'])-5)/5-42)
 event=dict(status='ACCEPTED',fold=a,history=h,pair_profile=p);events.append(event);out.write_text(json.dumps(events,indent=2))
 print('POINT',a['r'],a['c'],a['m'],'BRACKETS',len(p['root_sign_brackets']),flush=True)
 if len(p['root_sign_brackets'])>=3:
  (HERE/'K1_CANDIDATE_LOGM.json').write_text(json.dumps(event,indent=2));break
