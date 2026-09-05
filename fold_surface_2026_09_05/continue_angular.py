"""Fixed log-amplitude continuation of the fold, preserving decimal parameters.
Numerical continuation, not an exhaustive root-count certificate.
"""
import json,sys
from pathlib import Path
import mpmath as m
from budget import call,HERE
m.mp.dps=40
def st(x):return m.nstr(x,35)
def ev(z,c,K,purpose,tol='2e-18'):
 return call(dict(r=st(m.exp(z)),c=st(c),K=st(K),tol=tol),purpose,engine='angular_ld.py')
def correct(z,c,K,purpose):
 history=[]
 for i in range(9):
  a=ev(z,c,K,purpose);history.append(a)
  if a['status']!='NUMERICAL_ONLY':return None,history
  F=m.matrix([a['L'],a['L_z']]);J=m.matrix([[a['L_c'],a['L_K']],[a['L_zc'],a['L_zK']]])
  if max(abs(x) for x in F)<m.mpf('2e-10'):
   a['jacobian_cK_determinant']=st(m.det(J));a['correction_history']=history[:-1]
   return a,history
  d=m.lu_solve(J,F);c-=d[0];K-=d[1]
 return None,history

def profile(a):
 z=m.log(m.mpf(a['r']));c=m.mpf(a['c']);K=m.mpf(a['K']);delta=m.mpf('.008')*abs(m.mpf(a['L_zz'])/m.mpf(a['L_c']))
 cp=c-m.sign(m.mpf(a['L_c'])*m.mpf(a['L_zz']))*delta
 samples=[]
 for dz in [-12,-7,-4,-2,-1,-.3,0,.3,1,2,4,7,12]:
  b=ev(z+dz,cp,K,'angular pair-side full-turn profile');samples.append(dict(log_offset=dz,result=b))
 brackets=[];stationary=[]
 for p,q in zip(samples,samples[1:]):
  b,d=p['result'],q['result']
  if b['status']!= 'NUMERICAL_ONLY' or d['status']!='NUMERICAL_ONLY':continue
  if m.mpf(b['L'])*m.mpf(d['L'])<0:brackets.append(dict(left=b,right=d))
  if m.mpf(b['L_z'])*m.mpf(d['L_z'])<0:stationary.append(dict(left=b,right=d))
 return dict(c=st(cp),delta_c=st(cp-c),samples=samples,root_sign_brackets=brackets,stationary_sign_brackets=stationary,exhaustive_root_coverage=False)

def main():
 out=HERE/'events_angular_ld.json';events=json.loads(out.read_text()) if out.exists() else []
 if events:a=next(x['fold'] for x in reversed(events) if x['status']=='ACCEPTED')
 else:a=dict(r='204584.1771410815',c='1.573977181094133',K='6.7004205261144012')
 z=m.log(m.mpf(a['r']));c=m.mpf(a['c']);K=m.mpf(a['K']);step=m.mpf('.7')
 for i in range(int(sys.argv[1]) if len(sys.argv)>1 else 12):
  if events:
   J=m.matrix([[a['L_c'],a['L_K']],[a['L_zc'],a['L_zK']]])
   tangent=m.lu_solve(J,-m.matrix([a['L_z'],a['L_zz']]))
   zn=z+step;cn=c+step*tangent[0];Kn=K+step*tangent[1]
  else:zn,cn,Kn=z,c,K
  b,h=correct(zn,cn,Kn,'long-double fixed-amplitude fold continuation')
  if b is None:
   events.append(dict(status='CORRECTOR_UNRESOLVED',target_log_r=st(zn),history=h));out.write_text(json.dumps(events,indent=2));break
  a=b;z=zn;c=m.mpf(a['c']);K=m.mpf(a['K']);p=profile(a)
  event=dict(status='ACCEPTED',fold=a,pair_profile=p);events.append(event);out.write_text(json.dumps(events,indent=2))
  print('POINT',a['r'],a['c'],a['K'],'BRACKETS',len(p['root_sign_brackets']),flush=True)
  if len(p['root_sign_brackets'])>=3:
   (HERE/'K1_CANDIDATE_ANGULAR.json').write_text(json.dumps(event,indent=2));break
if __name__=='__main__':main()
