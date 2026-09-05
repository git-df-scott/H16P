"""Binary128 full returns for the endpoint fields missed by the double counter."""
import json
import mpmath as mp
from concurrent.futures import ProcessPoolExecutor
from run_d1 import HERE,dec,st,hp

def work(task):
    row,off=task;p=row['parameters'];r=dec(p['fold_r'])*mp.exp(dec(off))
    out=hp(r,p['c'],p['m'],p['beta'],tol='2e-29',full=True,purpose='D1 endpoint full-return replay')
    return dict(label=row['label'],kind=row['kind'],coefficient_vector=row['coefficient_vector'],log_offset=off,result=out)

rows=[json.loads(s) for s in (HERE/'fields.jsonl').read_text().splitlines()]
rows=[r for r in rows if r['label'] in ['positive_extension_1','negative_extension_1','center_0.0000000001'] and r['kind'] in ['pair','hopf']]
with ProcessPoolExecutor(max_workers=3) as pool:
    for a in pool.map(work,[(r,off) for r in rows for off in ['-.8','0','.8','2']]):
        with (HERE/'endpoint_full.jsonl').open('a') as f:f.write(json.dumps(a)+'\n')
        print(a['label'],a['kind'],a['log_offset'],a['result'].get('status'),a['result'].get('L'),flush=True)
