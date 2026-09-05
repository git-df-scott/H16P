"""Binary128 continuation using the reviewed two-sided equations."""
import json,sys
import mpmath as m
import continue_half as base
from budget import call,HERE
m.mp.dps=50

def ev(z,c,K,purpose):
 return call(dict(r=base.st(m.exp(z)),c=base.st(c),K=base.st(K),tol='2e-28'),purpose,engine='half_quad.py')
def correct(z,c,K):
 history=[]
 for i in range(8):
  a=ev(z,c,K,'binary128 two-sided fold correction');history.append(a)
  if a['status']!='NUMERICAL_TWO_HALF_PASSAGES':return None,history
  F=m.matrix([a['F'],a['G']]);J=m.matrix([[a['F_c'],a['F_K']],[a['G_c'],a['G_K']]])
  if abs(F[0])<m.mpf('1e-26') and abs(F[1])<m.mpf('1e-20'):
   a['jacobian_cK_determinant']=base.st(m.det(J));return a,history
  d=m.lu_solve(J,F);c-=d[0];K-=d[1]
 return None,history
base.ev=ev
# Refinement tolerance must resolve the very small bounded-section splitting.
s=__import__('inspect').getsource(base.refine).replace("1e-17","1e-28")
exec(s,base.__dict__)

out=HERE/'events_quad.json';events=json.loads(out.read_text()) if out.exists() else []
a=next(x['fold'] for x in reversed(events) if x['status']=='ACCEPTED') if events else json.loads((HERE/'events_half.json').read_text())[-1]['fold']
z=m.log(m.mpf(a['r']));c=m.mpf(a['c']);K=m.mpf(a['K']);step=m.mpf(sys.argv[2] if len(sys.argv)>2 else '2')
for i in range(int(sys.argv[1]) if len(sys.argv)>1 else 6):
 J=m.matrix([[a['F_c'],a['F_K']],[a['G_c'],a['G_K']]])
 tangent=m.lu_solve(J,-m.matrix([a['F_z'],a['G_z']]))
 zn=z+step;cn=c+step*tangent[0];Kn=K+step*tangent[1]
 b,h=correct(zn,cn,Kn)
 if b is None:events.append(dict(status='CORRECTOR_UNRESOLVED',target_log_r=base.st(zn),history=h));out.write_text(json.dumps(events,indent=2));break
 a=b;z=zn;c=m.mpf(a['c']);K=m.mpf(a['K']);p=base.profile(a)
 event=dict(status='ACCEPTED',fold=a,history=h,pair_profile=p);events.append(event);out.write_text(json.dumps(events,indent=2))
 print('POINT',a['r'],a['c'],a['K'],'BRACKETS',len(p['root_sign_brackets']),flush=True)
 if len(p['root_sign_brackets'])>=3:
  (HERE/'K1_CANDIDATE_QUAD.json').write_text(json.dumps(event,indent=2));break
