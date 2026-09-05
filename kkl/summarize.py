"""Summarize existing evidence without making any new return evaluations."""
from pathlib import Path
import hashlib
import json
import platform
import scipy
import numpy
import sympy

HERE=Path(__file__).resolve().parent
DATA=HERE/'data'
rows=[json.loads(line) for line in (DATA/'returns.jsonl').read_text().splitlines()]
assert [r['evaluation'] for r in rows]==list(range(1,len(rows)+1))
assert len(rows)<=4096
events=[]
for file in sorted(DATA.glob('*.jsonl')):
    if file.name in ('returns.jsonl','continuation_events.jsonl'): continue
    for line_number,line in enumerate(file.read_text().splitlines(),1):
        event=json.loads(line)
        event['source_file']=file.name
        event['source_line']=line_number
        event['last_evaluation']=max(event.get(nest,{}).get('evaluation',0) for nest in ['origin','remote'])
        events.append(event)
assert len(events)<=256
events.sort(key=lambda event:event['last_evaluation'])
(DATA/'continuation_events.jsonl').write_text(''.join(json.dumps(event)+'\n' for event in events))
last=json.loads((DATA/'last_inside_boundary.jsonl').read_text().splitlines()[-1])
summary={
    'status':'NUMERICAL_BRANCH_CHECKPOINT_NOT_WHOLE_BOX_EXCLUSION',
    'precursor_found':False,'five_cycle_candidate_found':False,'five_cycle_certificate':False,
    'charged_evaluations':len(rows),'evaluation_allowance':4096,'remaining_allowance':4096-len(rows),
    'continuation_events':len(events),'accepted_numerical_points':sum(e['status']=='ACCEPTED_NUMERICAL_POINT' for e in events),
    'cpu_seconds_inside_evaluators':sum(r['result'].get('cpu_seconds',0) for r in rows),
    'wall_seconds_including_subprocess_startup':sum(r['total_wall_seconds'] for r in rows),
    'largest_evaluator_cpu_seconds':max(r['result'].get('cpu_seconds',0) for r in rows),
    'last_admissible_field':{k:last[k] for k in ['c','alpha','beta','K']},
    'origin_root':last['origin']['r'],'origin_multiplier':last['origin']['R_r'],
    'remote_root':last['remote']['r'],'remote_multiplier':last['remote']['R_r'],
    'phase_end_reason':'Selected path reached K=1/64 and remote section-radius allowance before a detected origin fold.',
    'not_excluded':['unvisited cycle sheets','additional zeros between section samples','other parameter paths','other infinity strata','the KKL precursor in the full specified box'],
    'historical_engine_note':'Evaluations 1-131 used fixed-time projected first derivatives; 132 onward used transverse determinants. Historical rows were retained unchanged. Original intermediate source hashes were not captured; this manifest records the final reviewed implementation, not bitwise provenance for those earlier rows.',
    'software':{'python':platform.python_version(),'numpy':numpy.__version__,'scipy':scipy.__version__,'sympy':sympy.__version__}
}
(DATA/'strike_summary.json').write_text(json.dumps(summary,indent=2)+'\n')
files=[p for p in HERE.rglob('*') if p.is_file() and '__pycache__' not in p.parts and p.name!='SHA256SUMS']
(HERE/'SHA256SUMS').write_text(''.join(hashlib.sha256(p.read_bytes()).hexdigest()+'  '+str(p.relative_to(HERE))+'\n' for p in sorted(files)))
print(json.dumps(summary,indent=2))
