"""Final bounded continuation increment on the primary positive-K sheet."""
import sys,json,os
sys.argv=['continue_quad.py','0']
import continue_quad as q
import mpmath as m
from budget import call,HERE,used

def ev(z,c,K,purpose):
 return call(dict(r=q.base.st(m.exp(z)),c=q.base.st(c),K=q.base.st(K),tol='2e-26'),purpose,engine='half_quad.py')
q.ev=ev;q.base.ev=ev
for i in range(int(os.environ.get("FOLD_LAST_STEPS","3"))):
 if used()>3300:break
 a=q.a;z=m.log(m.mpf(a['r']));c=m.mpf(a['c']);K=m.mpf(a['K']);step=m.mpf('2')
 J=m.matrix([[a['F_c'],a['F_K']],[a['G_c'],a['G_K']]])
 t=m.lu_solve(J,-m.matrix([a['F_z'],a['G_z']]))
 b,h=q.correct(z+step,c+step*t[0],K+step*t[1])
 if b is None:
  q.events.append(dict(status='CORRECTOR_UNRESOLVED',target_log_r=q.base.st(z+step),history=h,protocol='2e-26 final tolerance'));q.out.write_text(json.dumps(q.events,indent=2));break
 q.a=b;p=q.base.profile(b);event=dict(status='ACCEPTED',fold=b,history=h,pair_profile=p,protocol='binary128 2e-26 final tolerance');q.events.append(event);q.out.write_text(json.dumps(q.events,indent=2))
 print('POINT',b['r'],b['c'],b['K'],'BRACKETS',len(p['root_sign_brackets']),flush=True)
 if len(p['root_sign_brackets'])>=3:
  (HERE/'K1_CANDIDATE_QUAD.json').write_text(json.dumps(event,indent=2));break
