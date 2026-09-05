"""Well-conditioned two-sided shooting; matching is a full closed orbit only at F=0.
Coordinates z=log r; F=w_forward(pi)-w_backward(-pi), G=M_f-M_b.
At F=0, exp(G) is the full-return multiplier. The fold is F=G=0.
"""
import json,sys
import mpmath as m
from budget import call,HERE
m.mp.dps=45
def st(x):return m.nstr(x,38)
def ev(z,c,K,purpose):
 return call(dict(r=st(m.exp(z)),c=st(c),K=st(K),tol='2e-18'),purpose,engine='half_ld.py')
def correct(z,c,K):
 history=[]
 for i in range(10):
  a=ev(z,c,K,'two-sided fixed-amplitude fold correction');history.append(a)
  if a['status']!='NUMERICAL_TWO_HALF_PASSAGES':return None,history
  F=m.matrix([a['F'],a['G']]);J=m.matrix([[a['F_c'],a['F_K']],[a['G_c'],a['G_K']]])
  if abs(F[0])<m.mpf('1e-16') and abs(F[1])<m.mpf('1e-12'):
   a['jacobian_cK_determinant']=st(m.det(J));return a,history
  d=m.lu_solve(J,F);c-=d[0];K-=d[1]
 return None,history

def refine(b,d,c,K):
 left=m.log(m.mpf(b['r']));right=m.log(m.mpf(d['r']));history=[]
 for i in range(8):
  fb,fd=m.mpf(b['F']),m.mpf(d['F']);z=(left*fd-right*fb)/(fd-fb)
  if not left+m.mpf('.03')*(right-left)<z<right-m.mpf('.03')*(right-left):z=(left+right)/2
  a=ev(z,c,K,'two-sided simple-root refinement');history.append(a)
  if a['status']!='NUMERICAL_TWO_HALF_PASSAGES':
   return dict(left=b,right=d,history=history,status='UNRESOLVED_REFINEMENT')
  fa=m.mpf(a['F'])
  if fa*fb>0:left=z;b=a
  else:right=z;d=a
  if abs(fa)<m.mpf('1e-17') or right-left<m.mpf('1e-6'):break
 return dict(left=b,right=d,approximation=a,history=history,status='NUMERICAL_SIGN_BRACKET',stability='stable' if m.mpf(a['G'])<0 else 'unstable')

def profile(a):
 z=m.log(m.mpf(a['r']));c=m.mpf(a['c']);K=m.mpf(a['K'])
 # F_zz = exp(M_b)*G_z at fold. Move into the local negative minimum.
 curvature=m.exp(m.mpf(a['backward_log_sensitivity']))*m.mpf(a['G_z'])
 cp=c-m.sign(m.mpf(a['F_c'])*curvature)*m.mpf('.008')*abs(curvature/m.mpf(a['F_c']))
 samples=[]
 offsets=sorted(set([-float(z)-6,-12,-7,-4,-2,-1,-.3,0,.3,1,2,4,7,12]))
 for dz in offsets:samples.append(dict(log_offset=dz,result=ev(z+dz,cp,K,'two-sided pair-side stationary/root profile')))
 brackets=[];stationary=[]
 for p,q in zip(samples,samples[1:]):
  b,d=p['result'],q['result']
  if b['status']!='NUMERICAL_TWO_HALF_PASSAGES' or d['status']!='NUMERICAL_TWO_HALF_PASSAGES':continue
  if m.mpf(b['F'])*m.mpf(d['F'])<0:brackets.append(refine(b,d,cp,K))
  if m.mpf(b['G'])*m.mpf(d['G'])<0:stationary.append(dict(left=b,right=d))
 return dict(c=st(cp),K=st(K),delta_c=st(cp-c),samples=samples,root_sign_brackets=brackets,stationary_sign_brackets=stationary,exhaustive_root_coverage=False)

def main():
 out=HERE/'events_half.json';events=json.loads(out.read_text()) if out.exists() else []
 a=next(x['fold'] for x in reversed(events) if x['status']=='ACCEPTED') if events else json.loads((HERE/'half_control.json').read_text())
 z=m.log(m.mpf(a['r']));c=m.mpf(a['c']);K=m.mpf(a['K']);step=m.mpf(sys.argv[2] if len(sys.argv)>2 else '1')
 for i in range(int(sys.argv[1]) if len(sys.argv)>1 else 12):
  J=m.matrix([[a['F_c'],a['F_K']],[a['G_c'],a['G_K']]])
  tangent=m.lu_solve(J,-m.matrix([a['F_z'],a['G_z']]))
  zn=z+step;cn=c+step*tangent[0];Kn=K+step*tangent[1]
  b,h=correct(zn,cn,Kn)
  if b is None:events.append(dict(status='CORRECTOR_UNRESOLVED',target_log_r=st(zn),history=h));out.write_text(json.dumps(events,indent=2));break
  a=b;z=zn;c=m.mpf(a['c']);K=m.mpf(a['K']);p=profile(a)
  event=dict(status='ACCEPTED',fold=a,history=h,pair_profile=p);events.append(event);out.write_text(json.dumps(events,indent=2))
  print('POINT',a['r'],a['c'],a['K'],'BRACKETS',len(p['root_sign_brackets']),flush=True)
  if len(p['root_sign_brackets'])>=3:
   (HERE/'K1_CANDIDATE_HALF.json').write_text(json.dumps(event,indent=2));break
if __name__=='__main__':main()
