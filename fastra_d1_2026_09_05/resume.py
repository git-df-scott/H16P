from engine import *
# Load definitions only, preserving the original run driver as the experiment plan.
import types
source=(HERE/'run_sheets.py').read_text().split("if args.sheet=='positive_center':")[0]
sheet=sys.argv[1];sys.argv=['run_sheets.py',sheet];ns={};exec(source,ns)
accept,advance=ns['accept'],ns['advance']
if sheet=='positive_infinity':
 events=[json.loads(l) for l in (HERE/'events_positive_infinity.jsonl').read_text().splitlines()]
 a=next(e['fold'] for e in reversed(events) if 'fold' in e)
 for j,z in enumerate([20,24,28,32,36,40,42],19):
  a=advance(a,z,'K','z',j)
  if a is None:break
else:
 events=[json.loads(l) for l in (HERE/'events_negative.jsonl').read_text().splitlines()];a=next(e['fold'] for e in reversed(events) if 'fold' in e)
 for i,k in ([] if 'm' in a else [(12,'-20'),(13,'-30')]):
  seeds=[e['fold'] for e in json.loads((ROOT/'fold_surface_2026_09_05/events_negative.json').read_text()) if e.get('status')=='ACCEPTED']
  seed=min(seeds,key=lambda e:abs(mp.mpf(e['K'])-mp.mpf(k)))
  b,h=correct(mp.log(mp.mpf(seed['r'])),mp.mpf(seed['c']),mp.mpf(k))
  a=accept(b,h,i,'K') if b else None
  if a is None:break
 if a is not None:
  a=half(mp.mpf(a['r']),mp.mpf(a['c']),-mp.mpf(a['alpha']),'m') if 'm' not in a else a
  previous_m=mp.mpf(a['m'])
  for j,M in enumerate([350,450,700,1300,3000,10000,100000,1000000,1e7,1e8,1e9,1e10,1e11,1e12,1e13],14):
   if M<=previous_m:continue
   candidates=[]
   for filename in ['events_m.json','events_logm.json']:
    candidates += [e['fold'] for e in json.loads((ROOT/'fold_surface_2026_09_05'/filename).read_text()) if e.get('status')=='ACCEPTED']
   seed=min([a]+candidates,key=lambda q:abs(mp.log(mp.mpf(q['m'])/M)))
   z=mp.log(mp.mpf(seed['r']));c=mp.mpf(seed['c']);p=mp.mpf(seed['m']);J=mp.matrix([[seed['F_z'],seed['F_c']],[seed['G_z'],seed['G_c']]])
   tangent=mp.lu_solve(J,-p*mp.matrix([seed['F_m'],seed['G_m']]))
   dl=mp.log(M/p);b,h=correct(z+dl*tangent[0],c+dl*tangent[1],mp.mpf(M),'m')
   if b is None:
    append('events_negative.jsonl',dict(status='CORRECTOR_UNRESOLVED',target_m=str(M),history=h));break
   a=accept(b,h,j,'m')
