from engine import *
from matching import matching_profile
rows={r['label']:r for r in (json.loads(s) for s in (HERE/'fields.jsonl').read_text().splitlines())}
labels=json.loads((HERE/'summary.json').read_text())['sign_disagreements']
for label in labels:
 row=rows[label];q=matching_profile(row['rational_vector'],np.log(6.76),du=1.,tol='2e-28',noise=1e-24)
 append('precision_repairs.jsonl',dict(label=label,rational_vector=row['rational_vector'],matching_profile=q,reason='Investigate outside/edge disagreement caused by discarded near-zero endpoints'))
 print(label,q['stability'],q['outside_outer'],q['edge']['sign'],flush=True)
