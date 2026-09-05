"""Follow the inherited three directions, profiling every accepted fold and unfolding.
Run one named direction per process. This is a deterministic continuation, not a sweep.
"""
from engine import *
import argparse
parser=argparse.ArgumentParser();parser.add_argument('sheet',choices=['positive_center','positive_infinity','negative']);args=parser.parse_args()
start=json.loads((HERE/'start_fold.json').read_text())

def accept(a,h,index,chart):
 label=f'{args.sheet}_{index:03d}';z=float(mp.log(mp.mpf(a['r'])));c=mp.mpf(a['c']);p=mp.mpf(a[chart]);M=-mp.mpf(a['alpha']);K=mp.mpf(a['K'])
 umax=max(40.,z+3.)
 # Fold coefficient and pair offset are rationalized BEFORE counting.
 fv=vector(st(c),**{chart:st(p)})
 fold_count=field_record(label+'_fold',fv,z,umax)
 curv=mp.exp(mp.mpf(a['backward_log_sensitivity']))*mp.mpf(a['G_z'])
 cp=c-mp.mpf('.02')*curv/mp.mpf(a['F_c'])
 pv=vector(st(cp),**{chart:st(p)})
 pair=field_record(label+'_pair',pv,z,umax)
 # Reverse the beta=0 weak-focus stability to create the extra innermost Hopf cycle.
 # For positive K the old pair is SU, giving USU after beta<0; negative K reverses it.
 cf=floats(pv);dd,ss,_=returns(cf,[z-1.5],umax=umax+12)
 margin=abs(dd[0]) if ss[0]==0 else abs(float(a['G_z']))*.01
 b=min(margin*np.sqrt(float(M))/(20*np.pi),.0005*np.sqrt(float(M)))
 b=-float(mp.sign(K))*b
 hopfs=[]
 for scale in [1.,.1]:
  bv=vector(st(cp),**{chart:st(p)},beta=rat(format(b*scale,'.17g')))
  hopfs.append(field_record(label+f'_hopf_{scale:g}',bv,z,umax))
 event=dict(label=label,sheet=args.sheet,status='ACCEPTED_NUMERICAL_FOLD',chart=chart,fold=a,history=h,fold_field_label=fold_count['label'],pair_field_label=pair['label'],hopf_field_labels=[r['label'] for r in hopfs],rational_fold_vector=fv,rational_pair_vector=pv)
 append('events_'+args.sheet+'.jsonl',event)
 print('ACCEPTED',label,'r',a['r'],'c',a['c'],'K',a['K'],flush=True)
 return a

def advance(a,target,chart,fixed,index):
 z=mp.log(mp.mpf(a['r']));c=mp.mpf(a['c']);p=mp.mpf(a[chart]);t=mp.mpf(target)
 keys=['z','c'] if fixed=='p' else ['c',chart]
 J=mp.matrix([[a['F_'+k] for k in keys],[a['G_'+k] for k in keys]])
 param=chart if fixed=='p' else 'z'
 tangent=mp.lu_solve(J,-mp.matrix([a['F_'+param],a['G_'+param]]))
 delta=t-p if fixed=='p' else t-z
 if fixed=='p':zz=z+tangent[0]*delta;cc=c+tangent[1]*delta;pp=t
 else:zz=t;cc=c+tangent[0]*delta;pp=p+tangent[1]*delta
 b,h=correct(zz,cc,pp,chart,fixed)
 if b is None:
  append('events_'+args.sheet+'.jsonl',dict(status='CORRECTOR_UNRESOLVED',target=str(target),chart=chart,fixed=fixed,history=h))
  print('FAILED',args.sheet,target,flush=True);return None
 return accept(b,h,index,chart)

if args.sheet=='positive_center':
 a=accept(start,[],0,'K')
 for i,k in enumerate(['.0009765625','.000244140625','.00006103515625','1e-5','1e-6','1e-7','1e-8','1e-9','1e-10'],1):
  a=advance(a,k,'K','p',i)
  if a is None:break
elif args.sheet=='positive_infinity':
 a=start
 for i,k in enumerate(['.00390625','.015625','.0625','.125','.2242157839673','.5','1','1.5','2','3','4','5','5.5','6'],1):
  a=advance(a,k,'K','p',i)
  if a is None:break
 if a is not None:
  for j,z in enumerate([8,10,12,16,20,24,28,32,36,40,42],i+1):
   a=advance(a,z,'K','z',j)
   if a is None:break
else:
 # The negative sheet is separate across the exact center, so seed it from its own archived point.
 a=next(e['fold'] for e in json.loads((ROOT/'fold_surface_2026_09_05/events_negative.json').read_text()) if e.get('status')=='ACCEPTED')
 a,h=correct(mp.log(mp.mpf(a['r'])),mp.mpf(a['c']),mp.mpf(a['K']))
 a=accept(a,h,0,'K')
 for i,k in enumerate(['-.001','-.01','-.05','-.1','-.3','-.7','-1','-2','-3','-5','-10','-20','-30'],1):
  a=advance(a,k,'K','p',i)
  if a is None:break
 if a is not None:
  # Switch to m before c=5/11; its coefficient chart is regular at K=-42.
  a=half(mp.mpf(a['r']),mp.mpf(a['c']),-mp.mpf(a['alpha']),'m')
  for j,M in enumerate([350,450,700,1300,3000,10000,100000,1000000,10000000,1e8,1e9,1e10,1e11,1e12,1e13],i+1):
   # Large coefficient jumps use the nearest inherited fold as predictor; the target is re-corrected.
   candidates=[]
   for filename in ['events_m.json','events_logm.json']:
    candidates += [e['fold'] for e in json.loads((ROOT/'fold_surface_2026_09_05'/filename).read_text()) if e.get('status')=='ACCEPTED']
   seed=min([a]+candidates,key=lambda q:abs(mp.log(mp.mpf(q['m'])/M)))
   z=mp.log(mp.mpf(seed['r']));c=mp.mpf(seed['c']);p=mp.mpf(seed['m'])
   J=mp.matrix([[seed['F_z'],seed['F_c']],[seed['G_z'],seed['G_c']]])
   tangent=mp.lu_solve(J,-p*mp.matrix([seed['F_m'],seed['G_m']]))
   dl=mp.log(M/p);b,h=correct(z+dl*tangent[0],c+dl*tangent[1],mp.mpf(M),'m')
   if b is None:
    append('events_'+args.sheet+'.jsonl',dict(status='CORRECTOR_UNRESOLVED',target_m=str(M),history=h));break
   a=accept(b,h,j,'m')
