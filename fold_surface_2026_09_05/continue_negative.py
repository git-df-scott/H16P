"""Regular negative-K sheet selected by the finite-amplitude center blow-up.
It is separated from the original regular sheet by a degenerate center annulus.
"""
import json
import mpmath as m
import continue_half as base
from budget import HERE
m.mp.dps=45

def correct(z,c,K):
 hist=[]
 for i in range(12):
  a=base.ev(z,c,K,'negative-K center-selected fold correction');hist.append(a)
  if a['status']!='NUMERICAL_TWO_HALF_PASSAGES':return None,hist
  F=m.matrix([a['F'],a['G']]);J=m.matrix([[a['F_z'],a['F_c']],[a['G_z'],a['G_c']]])
  if abs(F[0])/abs(K)<m.mpf('1e-13') and abs(F[1])/abs(K)<m.mpf('1e-11') and m.exp(z)>1:
   a['jacobian_zc_determinant']=base.st(m.det(J));return a,hist
  d=m.lu_solve(J,F)
  if abs(d[0])>.3:d*=m.mpf('.3')/abs(d[0])
  z-=d[0];c-=d[1]
  if c<=m.mpf(5)/11 or c>2:return None,hist
 return None,hist

seed=json.loads((HERE/'center_binary128.json').read_text())['rows'][-1]['result']
z=m.log(m.mpf(seed['r']));c=m.mpf(seed['c']);Kold=m.mpf(seed['K']);a=seed
out=HERE/'events_negative.json';events=json.loads(out.read_text()) if out.exists() else []
if events:
 a=next(e['fold'] for e in reversed(events) if e['status']=='ACCEPTED');z=m.log(m.mpf(a['r']));c=m.mpf(a['c']);Kold=m.mpf(a['K'])
for target in ['-.0001','-.001','-.01','-.03','-.05','-.1','-.2','-.4','-.7','-1','-1.5','-2','-3','-4','-6','-8','-10','-15','-20','-30','-40']:
 if m.mpf(target)>=Kold:continue
 K=m.mpf(target);J=m.matrix([[a['F_z'],a['F_c']],[a['G_z'],a['G_c']]])
 t=m.lu_solve(J,-m.matrix([a['F_K'],a['G_K']]))
 b,h=correct(z+(K-Kold)*t[0],c+(K-Kold)*t[1],K)
 if b is None:
  events.append(dict(status='CORRECTOR_UNRESOLVED',target_K=target,history=h));out.write_text(json.dumps(events,indent=2));break
 a=b;z=m.log(m.mpf(a['r']));c=m.mpf(a['c']);Kold=K
 p=base.profile(a);event=dict(status='ACCEPTED',fold=a,history=h,pair_profile=p);events.append(event);out.write_text(json.dumps(events,indent=2))
 print('POINT',a['r'],a['c'],a['K'],'BRACKETS',len(p['root_sign_brackets']),flush=True)
 if len(p['root_sign_brackets'])>=3:
  (HERE/'K1_CANDIDATE_NEGATIVE.json').write_text(json.dumps(event,indent=2));break
