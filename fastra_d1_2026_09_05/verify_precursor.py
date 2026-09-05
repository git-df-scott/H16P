from engine import *
import ctypes
from scipy.optimize import brentq
lib=ctypes.CDLL(str(HERE/'.matching_quad_v2.so'));ptr=ctypes.POINTER(ctypes.c_double);iptr=ctypes.POINTER(ctypes.c_int)
lib.matching.argtypes=[ctypes.c_char_p]*4+[ctypes.c_int,ptr,ptr,ptr,iptr]
lib.matching_remote.argtypes=[ctypes.c_char_p]*5+[ctypes.c_int,ptr,ptr,ptr,iptr,ptr]
vec=vector('9688912553490597/10000000000000000','1/512',beta='-1/10000000')
cf=floats(vec);remote=next(p for p in rm.equilibria(cf) if p[0]<-1)
record=dict(rational_vector=vec,method='binary128 two-sided angular matching; rational coefficients and remote equilibrium polished in binary128',roots=[])
def ev(u,where,tol='2e-28'):
 us=np.array([u],float);F=np.empty(1);G=np.empty(1);ss=np.empty(1,np.int32);eq=np.empty(2)
 arg=[vec[11].encode(),str(-Q(vec[7])).encode(),vec[8].encode(),tol.encode()]
 if where=='origin':lib.matching(*arg,1,us.ctypes.data_as(ptr),F.ctypes.data_as(ptr),G.ctypes.data_as(ptr),ss.ctypes.data_as(iptr))
 else:lib.matching_remote(*arg,str(remote[0]).encode(),1,us.ctypes.data_as(ptr),F.ctypes.data_as(ptr),G.ctypes.data_as(ptr),ss.ctypes.data_as(iptr),eq.ctypes.data_as(ptr))
 if ss[0]:raise RuntimeError(str(ss[0]))
 return float(F[0]),float(G[0])
for where,brackets in [('origin',[(.08,.12),(5.,5.6),(8.,9.)]),('remote',[(1e6,1e8)])]:
 for l,r in brackets:
  lo,hi=np.log(l),np.log(r);fl,gl=ev(lo,where);fh,gh=ev(hi,where)
  if fl*fh>=0:
   record['roots'].append(dict(where=where,status='NO_SIGN_BRACKET',lo=l,hi=r,Flo=fl,Fhi=fh));continue
  u=brentq(lambda u:ev(u,where)[0],lo,hi,xtol=2e-11)
  f,g=ev(u,where);tight=ev(u,where,'2e-30');record['roots'].append(dict(where=where,r=float(np.exp(u)),u=u,initial_bracket=[l,r],F_bracket=[fl,fh],F=f,log_multiplier=g,multiplier=float(np.exp(g)),tight_F=tight[0],tight_log_multiplier=tight[1],stability='S' if g<0 else 'U',section='positive horizontal' if where=='origin' else 'negative horizontal from remote equilibrium'))
  print(record['roots'][-1],flush=True)
(HERE/'verified_precursor.json').write_text(json.dumps(record,indent=2))
field_record('rational_precursor',vec,float(mp.log(mp.mpf('6.76'))))
