"""Summarize saved events only; no trajectory evaluations."""
from pathlib import Path
import json,hashlib,collections
from fractions import Fraction as Q
HERE=Path(__file__).resolve().parent
files=['events_decreasing','events_increasing','events_arclength','events_angular_ld','events_half','events_quad','events_negative','events_m','events_logm']
rows=[]
for file in files:
 p=HERE/(file+'.json')
 if not p.exists():continue
 for i,event in enumerate(json.loads(p.read_text())):
  fold=event.get('fold',{});old='return_data' in fold
  data=fold.get('return_data',fold);accepted=(fold.get('status')=='ACCEPTED_NUMERICAL_FOLD' or event.get('status')=='ACCEPTED')
  row=dict(file=p.name,event_index=i,status='ACCEPTED_NUMERICAL_FOLD' if accepted else 'UNRESOLVED_CORRECTION',raw_status=fold.get('status',event.get('status')))
  if accepted:
   c=Q(str(fold['c']));K=Q(str(fold['K']));alpha=-Q(str(fold['m'])) if 'm' in fold else -5*(K+42)/(11*c-5)
   A=c-Q(61,5);B=alpha-Q(111,5);C=2*alpha-10;D=alpha
   disc=B*B*C*C-4*A*C**3-4*B**3*D-27*A*A*D*D+18*A*B*C*D
   remote_sign=None
   if disc<0 and c>0 and c<Q(61,5) and alpha<0:
    if c>=Q(8,5):remote_sign=1
    else:
     xh=-21/(16-10*c);value=((A*xh+B)*xh+C)*xh+D
     remote_sign=(value>0)-(value<0)
   pair=event.get('pair',event.get('pair_profile',{}));brackets=pair.get('root_brackets',pair.get('root_sign_brackets',[]));stationary=pair.get('stationary_brackets',pair.get('stationary_sign_brackets',[]))
   row.update(c=str(fold['c']),K=str(fold['K']),r=str(fold['r']),section='curved P=0 maximum x' if file in files[:3] else 'positive horizontal ray',equations='F=G=0' if 'F' in data else 'L=L_z=0',residual=[data.get('F',data.get('L')),data.get('G',data.get('L_z'))],curvature=data.get('G_z',data.get('L_zz')),period=data.get('period_at_match',data.get('period')),root_sign_brackets=len(brackets),stationary_sign_brackets=len(stationary),exhaustive_root_coverage=False,exact_equilibrium_gate=dict(cubic_discriminant_sign=(disc>0)-(disc<0),unique_nonorigin_left_of_barrier=disc<0 and c>0 and c<Q(61,5) and alpha<0,remote_trace_sign=remote_sign),evaluation=data.get('evaluation'))
  rows.append(row)
ledger=[json.loads(s) for s in (HERE/'returns.jsonl').read_text().splitlines()]
summary=dict(base_commit='79001f70eb180331f5ac0b740f5d5aadfe833329',historical_evaluations=756,new_evaluations=len(ledger),campaign_evaluations=756+len(ledger),remaining=4096-756-len(ledger),evaluation_status_counts=dict(collections.Counter(x['result']['status'] for x in ledger)),purpose_counts=dict(collections.Counter(x['purpose'] for x in ledger)),events=rows,global_component_complete=False,three_origin_candidate=False,five_cycle_candidate=False)
(HERE/'component_summary.json').write_text(json.dumps(summary,indent=2))
lines=['# Fold-component event ledger','', 'All root counts below are numerical sign brackets, not exhaustive counts. Rejected corrections remain in the raw event files; none is classified as a mathematical endpoint.','', '| File / event | K | c | section r | pair root brackets | status |','|---|---:|---:|---:|---:|---|']
for x in rows:
 lines.append('| '+x['file']+' / '+str(x['event_index'])+' | '+' | '.join(str(x.get(k,'—')) for k in ['K','c','r','root_sign_brackets','status'])+' |')
(HERE/'EVENT_LEDGER.md').write_text('\n'.join(lines)+'\n')
manifest={str(p.relative_to(HERE.parent)):hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(HERE.glob('*')) if p.is_file() and p.suffix in ('.py','.cpp')}
(HERE/'source_manifest.json').write_text(json.dumps(manifest,indent=2))
print(json.dumps({k:v for k,v in summary.items() if k!='events' and k!='purpose_counts'},indent=2))
