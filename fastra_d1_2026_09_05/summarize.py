"""Generate the reviewable D1 ledger from saved records (no integrations)."""
import json,csv,hashlib,collections
from pathlib import Path
from fractions import Fraction
HERE=Path(__file__).resolve().parent;ROOT=HERE.parent

def read(name):return [json.loads(s) for s in (HERE/name).read_text().splitlines()] if (HERE/name).exists() else []
raw=read('fields.jsonl');rows={r['label']:r for r in raw}
for r in read('precision_repairs.jsonl'):
 rows[r['label']]['matching_profile']=r['matching_profile'];rows[r['label']]['precision_repaired']=True
for r in read('preserved_hopf.jsonl'):rows[r['label']]=r
for r in rows.values():
 r['selected_profile']=r.get('matching_profile',r.get('origin_profile',{}))
 r['count']=len(r['selected_profile'].get('roots',[]))
 baseline=r.get('baseline',[]);origin=next((n for n in baseline if n['pt'] is not None and sum(v*v for v in n['pt'])<1e-14),None)
 r['baseline_origin_count']=len(origin['roots']) if origin else None
campaign=[r for r in rows.values() if r['label'].startswith(('positive_','negative_'))]
pair=[r for r in campaign if '_pair' in r['label']]
hopf=[r for r in campaign if '_hopf' in r['label']]
events={}
for name in ['positive_center','positive_infinity','negative']:
 all_=read('events_'+name+'.jsonl');accepted=[e for e in all_ if e.get('status')=='ACCEPTED_NUMERICAL_FOLD'];events[name]=dict(accepted=len(accepted),rejected=sum(e.get('status')!='ACCEPTED_NUMERICAL_FOLD' for e in all_),endpoint=accepted[-1]['fold'])
sign_hist=collections.Counter();disagreements=[]
with open(HERE/'sign_map.csv','w') as f:
 cols=['label','c','K','beta','alpha','selected_method','count','stability','outside_sign','edge_sign','edge_kind','edge_radius','sign_comparison','uncertain_sign_changes','baseline_origin_count','rational_vector']
 w=csv.DictWriter(f,fieldnames=cols,lineterminator="\n");w.writeheader()
 for row in campaign:
  vec=row['rational_vector'];c=Fraction(vec[11]);m=-Fraction(vec[7]);K=m*(11*c-5)/5-42;p=row['selected_profile'];edge=p.get('edge',{});outside=p.get('outside_outer') or {};es=edge.get('sign');os=outside.get('sign')
  verdict='unresolved' if es not in (-1,1) or os not in (-1,1) else 'agree' if es==os else 'disagree'
  if '_pair' in row['label'] or '_hopf' in row['label']:
   sign_hist[verdict]+=1
   if verdict=='disagree':disagreements.append(row['label'])
  w.writerow(dict(label=row['label'],c=str(c),K=str(K),beta=vec[8],alpha=vec[7],selected_method='binary128_matching' if 'matching_profile' in row else 'double_full_return',count=row['count'],stability=p.get('stability'),outside_sign=os,edge_sign=es,edge_kind=edge.get('kind'),edge_radius=edge.get('r_valid'),sign_comparison=verdict,uncertain_sign_changes=len(p.get('uncertain_sign_changes',[])),baseline_origin_count=row['baseline_origin_count'],rational_vector=json.dumps(vec)))
summary=dict(status='OPEN',raw_field_records=len(raw),unique_fields=len(rows),campaign_fields=len(campaign),events=events,max_selected_origin_count=max(r['count'] for r in campaign),max_beta_zero_pair_count=max(r['count'] for r in pair),max_hopf_count=max(r['count'] for r in hopf),raw_double_profile_max=max(len(r.get('origin_profile',{}).get('roots',[])) for r in raw),raw_baseline_max=max(r['baseline_origin_count'] or 0 for r in rows.values()),sign_comparison=dict(sign_hist),sign_disagreements=disagreements,semistable_folds_counted_by_sign_change=False,interval_certification=False,exhaustive_root_coverage=False)
(HERE/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
# Compact root and stability ledger; full signs and trajectories remain in JSONL.
with open(HERE/'root_ledger.csv','w') as f:
 w=csv.writer(f,lineterminator="\n");w.writerow(['label','method','origin_count','origin_stability','origin_radii','baseline_nests','edge'])
 for r in rows.values():
  p=r['selected_profile'];w.writerow([r['label'],'binary128_matching' if 'matching_profile' in r else 'double_full_return',r['count'],p.get('stability'),json.dumps([x.get('r') for x in p.get('roots',[])]),json.dumps(r.get('baseline')),json.dumps(p.get('edge'))])
print(json.dumps({k:v for k,v in summary.items() if k!='events'},indent=2))
for name,e in events.items():print(name,e['accepted'],{k:e['endpoint'].get(k) for k in ['r','c','K','m']})

collapses=[]
for name in events:
 for e in read('events_'+name+'.jsonl'):
  if e.get('status')!='ACCEPTED_NUMERICAL_FOLD':continue
  fv=e['rational_fold_vector'];pv=e['rational_pair_vector']
  if all(float(Fraction(a))==float(Fraction(b)) for a,b in zip(fv,pv)):
   collapses.append(dict(label=e['label'],rational_fold_vector=fv,rational_pair_vector=pv,reason='Distinct exact fields round to an identical float64 vector'))
(HERE/'float64_field_collapses.json').write_text(json.dumps(collapses,indent=2)+'\n')
