from engine import *
import ctypes
lib=ctypes.CDLL(str(HERE/'.matching_quad.so'))
ptr=ctypes.POINTER(ctypes.c_double);iptr=ctypes.POINTER(ctypes.c_int)
lib.matching.argtypes=[ctypes.c_char_p]*4+[ctypes.c_int,ptr,ptr,ptr,iptr]
def match(vec,us,tol='2e-27'):
 us=np.ascontiguousarray(us,float);F=np.empty(len(us));G=np.empty(len(us));ss=np.empty(len(us),np.int32)
 lib.matching(vec[11].encode(),str(-Q(vec[7])).encode(),vec[8].encode(),tol.encode(),len(us),us.ctypes.data_as(ptr),F.ctypes.data_as(ptr),G.ctypes.data_as(ptr),ss.ctypes.data_as(iptr))
 return F,G,ss
if __name__=='__main__':
 rows=[json.loads(l) for l in (HERE/'fields.jsonl').read_text().splitlines()]
 for label in ['pilot_pair','pilot_beta_-1e-7','positive_infinity_019_fold']:
  row=next(x for x in rows if x['label']==label);us=np.array([-2.4,0.,1.8,2.3]) if label!='positive_infinity_019_fold' else np.array([19.55,19.58,20.,20.6,25.,40.])
  t=time.time();F,G,ss=match(row['rational_vector'],us);print(label,'seconds',time.time()-t,list(zip(us,F,G,ss)),flush=True)

def matching_profile(vec,z,du=1.0,umax=40.,tol='2e-28',noise=1e-24,umin=-11.5):
 us=np.unique(np.r_[np.arange(umin,umax+.001,du),[-8,-6,-4,-3,-2,-1,0,1],np.arange(z-.65,z+.651,.065),umax])
 F,G,ss=match(vec,us,tol)
 grid=[dict(u=float(u),F=float(f),G=float(g),status=int(s)) for u,f,g,s in zip(us,F,G,ss)]
 first=next((i for i,s in enumerate(ss) if s),len(ss));valid=grid[:first]
 roots=[];uncertain=[]
 for a,b in zip(valid,valid[1:]):
  if a['F']*b['F']>=0:continue
  rec=dict(u_lo=a['u'],u_hi=b['u'],F_lo=a['F'],F_hi=b['F'],stability='S' if a['F']>0 else 'U')
  if min(abs(a['F']),abs(b['F']))<noise:uncertain.append(rec);continue
  lo,hi=a['u'],b['u'];fl=a['F']
  for _ in range(9):
   mid=(lo+hi)/2;f,g,s=match(vec,[mid],tol)
   if s[0]:rec['refinement_failure']=int(s[0]);break
   if f[0]*fl>0:lo=mid;fl=f[0]
   else:hi=mid
  rec.update(r=float(np.exp((lo+hi)/2)),u_bracket=[lo,hi]);roots.append(rec)
 edge=dict(kind='scan_cap' if first==len(ss) else 'angular_chart_failure',u_valid=valid[-1]['u'] if valid else None,F=valid[-1]['F'] if valid else None)
 if first<len(ss):edge.update(status=int(ss[first]),u_failed=float(us[first]))
 if valid:edge.update(sign=int(np.sign(valid[-1]['F'])) if abs(valid[-1]['F'])>noise else 0,r_valid=float(np.exp(valid[-1]['u'])))
 outside=None
 if roots:
  u=roots[-1]['u_hi']+.03
  if u<valid[-1]['u']:
   f,g,s=match(vec,[u],tol);outside=dict(u=u,status=int(s[0]),F=float(f[0]),sign=int(np.sign(f[0])) if abs(f[0])>noise else 0)
 return dict(method='binary128 scaled two-sided matching over full configured radius grid',section='positive horizontal ray',tol=tol,noise=noise,roots=roots,stability=''.join(r['stability'] for r in roots),uncertain_sign_changes=uncertain,edge=edge,outside_outer=outside,grid=grid)
