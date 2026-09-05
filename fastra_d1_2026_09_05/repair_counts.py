"""Targeted precision repairs after the large-radius false-positive control.
Only archived fold fields are used; no random search or descent.
"""
from engine import *
from matching import match,matching_profile
from copy import deepcopy
mode=sys.argv[1]
rows=[json.loads(l) for l in (HERE/'fields.jsonl').read_text().splitlines()]
if mode=='center':
 done={json.loads(l)['label'] for l in (HERE/'precision_repairs.jsonl').read_text().splitlines()} if (HERE/'precision_repairs.jsonl').exists() else set()
 for row in rows:
  if not row['label'].startswith('positive_center_'):continue
  if int(row['label'].split('_')[2])<4 or '_fold' in row['label'] or row['label'] in done:continue
  q=matching_profile(row['rational_vector'],np.log(6.758),du=.5,umax=40,tol='2e-28',noise=1e-24)
  append('precision_repairs.jsonl',dict(label=row['label'],rational_vector=row['rational_vector'],matching_profile=q))
  print(row['label'],q['stability'],flush=True)
elif mode=='beta':
 events=[json.loads(l) for l in (HERE/'events_positive_infinity.jsonl').read_text().splitlines()]
 done={json.loads(l)['label'] for l in (HERE/'preserved_hopf.jsonl').read_text().splitlines()} if (HERE/'preserved_hopf.jsonl').exists() else set()
 for event in events:
  if event.get('status')!='ACCEPTED_NUMERICAL_FOLD':continue
  a=event['fold'];z=float(mp.log(mp.mpf(a['r'])))
  if z<16:continue
  label=event['label']+'_hopf_preserved'
  if label in done:continue
  vec=event['rational_pair_vector'];M=-float(Q(vec[7]))
  # Derivative of the bounded matching residual with respect to beta.
  f0,_,s0=match(vec,[z]);h=1e-7*np.sqrt(M);plus=vec.copy();minus=vec.copy();plus[8]=rat(format(h,'.17g'));minus[8]=rat(format(-h,'.17g'))
  fp,_,sp=match(plus,[z]);fm,_,sm=match(minus,[z]);fb=(fp[0]-fm[0])/(2*h)
  if s0[0] or sp[0] or sm[0]:raise RuntimeError('unresolved beta derivative')
  beta=-.1*abs(f0[0]/fb);vec=vec.copy();vec[8]=rat(format(beta,'.17g'))
  cf=floats(vec);baseline=sw.evaluate(cf)
  # Hopf radius shrinks roughly like sqrt(|beta|); include the full logarithmic gap.
  umin=min(-11.5,.5*np.log(abs(beta))-5)
  tol='2e-30' if z>=32 else '2e-28';noise=1e-28 if z>=32 else 1e-24
  # Supplement full matching scan by negative-log points for the tiny Hopf cycle.
  q=matching_profile(vec,z,du=1.,umax=max(40,z+3),tol=tol,noise=noise)
  us=np.arange(umin,-11.5,.5)
  if len(us):
   F,G,ss=match(vec,us,tol);tiny=[dict(u=float(u),F=float(f),G=float(g),status=int(s)) for u,f,g,s in zip(us,F,G,ss)]
   allgrid=tiny+q['grid'];q['grid']=allgrid
   for a,b in zip(allgrid,allgrid[1:]):
    if b['u']>-11.5:break
    if a['status'] or b['status'] or a['F']*b['F']>=0:continue
    if min(abs(a['F']),abs(b['F']))<noise:
     q['uncertain_sign_changes'].append(dict(u_lo=a['u'],u_hi=b['u'],F_lo=a['F'],F_hi=b['F']));continue
    q['roots'].insert(0,dict(u_lo=a['u'],u_hi=b['u'],F_lo=a['F'],F_hi=b['F'],r=float(np.exp((a['u']+b['u'])/2)),stability='S' if a['F']>0 else 'U',refinement='coarse extra-small-root bracket'))
   q['stability']=''.join(r['stability'] for r in q['roots'])
  row=dict(label=label,rational_vector=vec,baseline=baseline,matching_profile=q,beta_derivative_control=dict(F0=float(f0[0]),F_beta=float(fb),step=h),purpose='Preserve the large fold pair while adding the innermost Hopf cycle')
  append('preserved_hopf.jsonl',row);print(label,beta,q['stability'],flush=True)
  if len(q['roots'])>=4:
   (HERE/'FOUR_ORIGIN_TRIGGER_VALIDATED.json').write_text(json.dumps(row,indent=2));raise RuntimeError('STOP FOUR')
