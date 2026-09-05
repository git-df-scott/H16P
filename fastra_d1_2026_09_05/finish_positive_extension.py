"""Finish far-positive endpoints with the repaired matching error controller.
The earlier M tolerance over-resolved sensitivity and exhausted a step guard.
Position tolerance is retained; sensitivity tolerance is relaxed to 1e-24.
"""
from engine import *
import matching as ma
import ctypes
lib=ctypes.CDLL(str(HERE/'.matching_quad_v3.so'));lib.matching.argtypes=ma.lib.matching.argtypes;ma.lib=lib

def polish(a,z):
 c=mp.mpf(a['c']);M=-mp.mpf(a['alpha']);hist=[]
 def ev(c,M):
  v=vector(st(c),m=st(M));F,G,S=ma.match(v,[z],'2e-30');hist.append(dict(rational_vector=v,F=float(F[0]),G=float(G[0]),status=int(S[0])))
  if S[0]:raise RuntimeError('endpoint polish numerical failure')
  return mp.matrix([str(F[0]),str(G[0])])
 for _ in range(4):
  f=ev(c,M)
  if abs(f[0])<mp.mpf('3e-29') and abs(f[1])<mp.mpf('1e-20'):break
  h=mp.mpf('1e-5');dc=(ev(c+h,M)-ev(c-h,M))/(2*h);dm=(ev(c,M+h)-ev(c,M-h))/(2*h)
  d=mp.lu_solve(mp.matrix([[dc[0],dm[0]],[dc[1],dm[1]]]),f);c-=d[0];M-=d[1]
 a=dict(a);a.update(c=st(c),alpha=st(-M),m=st(M),K=st(M*(11*c-5)/5-42),F=st(f[0]),G=st(f[1]),scalar_polish=hist)
 return a

def field(label,vec,z,with_profile=True,beta=False):
 row=field_record(label,vec,z,supplement=False)
 if with_profile:
  umin=min(-11.5,.5*np.log(abs(float(Q(vec[8]))))-5) if beta else -11.5
  p=ma.matching_profile(vec,z,du=1.,umax=z+3,tol='2e-29',noise=1e-29,umin=umin)
  p['engine_revision']='matching_quad.cpp with separate sensitivity tolerance';row['matching_profile']=p
  append('precision_repairs.jsonl',dict(label=label,rational_vector=vec,matching_profile=p))
  print(label,p['stability'],p['edge'],flush=True)
  if len(p['roots'])>=4:raise RuntimeError('STOP FOUR ORIGIN')
 return row

events=[json.loads(l) for l in (HERE/'events_positive_infinity.jsonl').read_text().splitlines()]
a=next(e['fold'] for e in reversed(events) if e.get('label')=='positive_infinity_024_polished')
for index,z in [(25,42.)]:
 if z==42:
  a=dict(a);a.update(c=st(mp.mpf(a['c'])+mp.mpf('.00033')),alpha=st(mp.mpf(a['alpha'])+mp.mpf('.003')),r=st(mp.exp(z)))
 a=polish(a,z);label=f'positive_infinity_{index:03d}_polished';c=mp.mpf(a['c']);M=mp.mpf(a['m']);fv=vector(st(c),m=st(M))
 foldrow=field(label+'_fold',fv,z,False)
 # Curvature and actual-m c derivative of the well-conditioned matching residual.
 h=mp.mpf('1e-5');Fp,_,_=ma.match(vector(st(c+h),m=st(M)),[z],'2e-29');Fm,_,_=ma.match(vector(st(c-h),m=st(M)),[z],'2e-29');Fc=mp.mpf(str((Fp[0]-Fm[0])/(2*float(h))))
 Fu,Gu,Su=ma.match(fv,[z-.3,z,z+.3],'2e-30');curv=mp.mpf(str((Fu[0]+Fu[2]-2*Fu[1])/.09))
 a['finite_difference_curvature']=st(curv);a['finite_difference_F_c_at_m']=st(Fc);a['curvature_step']='0.3'
 for old in ['F_c','F_K','G_c','G_K','F_z','G_z','backward_log_sensitivity','forward_log_sensitivity']:a.pop(old,None)
 cp=c-mp.mpf('.02')*curv/Fc;pv=vector(st(cp),m=st(M));pairrow=field(label+'_pair',pv,z)
 f,_,ss=ma.match(pv,[z],'2e-30');h=1e-7*np.sqrt(float(M));vp=pv.copy();vm=pv.copy();vp[8]=rat(str(h));vm[8]=rat(str(-h));fp,_,_=ma.match(vp,[z]);fm,_,_=ma.match(vm,[z]);fb=(fp[0]-fm[0])/(2*h)
 beta=-.1*abs(f[0]/fb);bv=pv.copy();bv[8]=rat(format(beta,'.17g'));hopfrow=field(label+'_hopf_preserved',bv,z,True,True)
 controls=[]
 for tol in ['2e-28','2e-30']:
  F,G,S=ma.match(pv,[z-.3,z,z+.3],tol);controls.append(dict(tol=tol,F=F.tolist(),G=G.tolist(),status=S.tolist()))
 append('events_positive_infinity.jsonl',dict(status='ACCEPTED_NUMERICAL_FOLD',label=label,sheet='positive_infinity',chart='m',fold=a,rational_fold_vector=fv,rational_pair_vector=pv,fold_field_label=foldrow['label'],pair_field_label=pairrow['label'],hopf_field_labels=[hopfrow['label']],pair_tolerance_controls=controls,beta_derivative=float(fb)))
 print('ACCEPTED FINAL',label,a['c'],a['K'],a['r'],flush=True)
