"""Validate saved D1 records without rerunning integrations."""
from pathlib import Path
from fractions import Fraction as Q
import json,hashlib,math
HERE=Path(__file__).resolve().parent;ROOT=HERE.parent
errors=[];counts={}
records=[]
for file in HERE.glob('*.jsonl'):
 try:a=[json.loads(l) for l in file.read_text().splitlines()]
 except Exception as e:errors.append(f'{file.name}: {e}');continue
 counts[file.name]=len(a);records+=a
rows={r['label']:r for r in records if 'label' in r and 'rational_vector' in r}
for r in records:
 if 'rational_vector' in r:
  try:
   v=[Q(x) for x in r['rational_vector']]
   assert len(v)==12 and v[:7]==[0,0,1,1,1,0,0] and v[9:11]==[-10,Q(11,5)] and v[7]<0
  except Exception:errors.append('Invalid field '+r.get('label','unknown'))
 if r.get('status')=='ACCEPTED_NUMERICAL_FOLD':
  for key in ['fold_field_label','pair_field_label']:
   if r[key] not in rows:errors.append('Missing field '+r[key])
  for name in r['hopf_field_labels']:
   if name not in rows:errors.append('Missing Hopf '+name)
  a=r['fold']
  if abs(float(a['F']))>1e-19 or abs(float(a['G']))>1e-15:errors.append('Excess nominal fold residual '+r['label'])
v=json.loads((HERE/'verified_precursor.json').read_text());assert len(v['roots'])==4
assert ''.join(a['stability'] for a in v['roots'] if a['where']=='origin')=='USU'
assert ''.join(a['stability'] for a in v['roots'] if a['where']=='remote')=='S'
for a in v['roots']:
 assert a['F_bracket'][0]*a['F_bracket'][1]<0
 assert abs(a['log_multiplier']-a['tight_log_multiplier'])<1e-12
assert json.loads((HERE/'trigger_resolution.json').read_text())['all_rechecked_matching_endpoints_positive']
# Scan only this turn's new artifacts for accidental credentials; never print matches.
for file in list(HERE.rglob('*'))+[ROOT/'FASTRA_D1_REPORT_2026_09_05.md']:
 if file.is_file() and not file.name.startswith('.') and '__pycache__' not in str(file):
  if b'github' + b'_pat_' in file.read_bytes():errors.append('Credential-like data found in '+file.name)
result=dict(status='PASS' if not errors else 'FAIL',errors=errors,record_counts=counts,exact_vectors_validated=sum('rational_vector' in r for r in records),checks=['JSON parse','exact rational vector convention','accepted-event field references','nominal fold residual bounds','3+1 bracket and tolerance consistency','raw trigger rejection','credential scan'],interval_certification=False)
(HERE/'artifact_audit.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result,indent=2))
if errors:raise SystemExit(1)
