"""Validate delivered numerical records and generate review tables."""
from pathlib import Path
from fractions import Fraction as F
import json,csv
import mpmath as mp
mp.mp.dps=70
H=Path(__file__).parent
rows=[json.loads(l) for l in (H/'center_sign_map.jsonl').read_text().splitlines()]
edges=[json.loads(l) for l in (H/'domain_edge_map.jsonl').read_text().splitlines()]
assert len(rows)==len(edges)==24
assert len(set(F(r['K']) for r in rows))==24
max_delta=mp.mpf(0);root_brackets=0
for r,e in zip(rows,edges):
 assert r['K']==e['K'] and r['pair_vector']==e['vector']
 v=list(map(F,r['pair_vector']));assert -v[7]*(11*v[11]-5)/5-42==F(r['K'])
 assert len(r['roots'])==2
 for root in r['roots']:
  a,b=root['endpoint_returns'];assert a['status']==b['status']=='OK_NUMERICAL'
  assert mp.mpf(a['log_displacement'])*mp.mpf(b['log_displacement'])<0
  root_brackets+=1
 for name in ['outside','edge']:
  a,b=r[name],r[name+'_tight'];assert a['status']==b['status']=='OK_NUMERICAL'
  assert mp.mpf(a['log_displacement'])>0 and mp.mpf(b['log_displacement'])>0
  max_delta=max(max_delta,abs(mp.mpf(a['log_displacement'])-mp.mpf(b['log_displacement'])))
 assert e['sample']['status']==e['sample_tight']['status']=='OK_NUMERICAL'
 assert mp.mpf(e['sample']['log_displacement'])>0 and mp.mpf(e['sample_tight']['log_displacement'])>0
 assert e['comparison']==r['comparison']=='agree'
cols=['index','K','outside_sign','resolved_edge_sign','sample_log_radius','last_success_log_radius','first_failure_log_radius','first_failure_status','edge_log_displacement','comparison']
with (H/'domain_edge_map.csv').open('w') as f:
 w=csv.DictWriter(f,fieldnames=cols);w.writeheader()
 for e in edges:w.writerow(dict(index=e['index'],K=e['K'],outside_sign=e['outside_sign'],resolved_edge_sign=e['edge_sign'],sample_log_radius=e['sample_log_radius'],last_success_log_radius=e['last_success_log_radius'],first_failure_log_radius=e['first_failure_log_radius'],first_failure_status=e['first_failure_status'],edge_log_displacement=e['sample_tight']['log_displacement'],comparison=e['comparison']))
summary={'CE':False,'D1':'OPEN','D2':'OPEN','literal_Proposition_A':'GAP: false converse','corrected_A_implication':'VERIFIED including chart degeneracies','statement_C':'VERIFIED via Llibre-Schlomiuk 2004 Theorem 16(III), Figures 2-3','second_order_catalogue_component_audit':'INCOMPLETE: full article unavailable; no theorem/portrait number invented','K_values':24,'rational_vectors_checked':24,'tight_root_brackets_checked':root_brackets,'outside_grid_edge_agreements':24,'outside_resolved_edge_agreements':24,'sign_disagreements':0,'unresolved_sign_comparisons_after_budget_rechecks':0,'minimum_local_profile_abs_log_displacement':mp.nstr(min(abs(mp.mpf(p['log_displacement'])) for r in rows for p in r['local_profile'] if p['status']=='OK_NUMERICAL'),35),'maximum_outside_or_exp40_tolerance_difference':mp.nstr(max_delta,35),'true_return_domain_boundary_certified':False,'interval_certified':False,'additional_pair_excluded':False,'boundary_search_failures_remain_unresolved':True,'analytic_controls_passed':6}
(H/'summary.json').write_text(json.dumps(summary,indent=2)+'\n');print(json.dumps(summary,indent=2))
