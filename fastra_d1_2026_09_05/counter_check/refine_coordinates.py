"""Safeguarded full-return Newton refinement inside recorded grid brackets."""
from check_misses import *
rows=[json.loads(s) for s in (D1/'fields.jsonl').read_text().splitlines()]
misses=[json.loads(s) for s in (HERE/'missed_roots.jsonl').read_text().splitlines()]
path=HERE/'refined_coordinates.jsonl'
done={(r['label'],r['kind'],r['root_index']) for r in map(json.loads,path.read_text().splitlines())} if path.exists() else set()
for miss in misses:
    if miss['label']=='positive_extension_1':continue # retain full sign brackets at ill-conditioned endpoint
    row=next(r for r in rows if r['label']==miss['label'] and r['kind']==miss['kind'])
    for root in miss['roots']:
        if (miss['label'],miss['kind'],root['index']) in done:continue
        a,b=[dec(g['u_decimal']) for g in root['grid']];fa,fb=[dec(g['quad']['L']) for g in root['grid']];z=(a+b)/2;history=[]
        for it in range(9):
            q=quad(row,mp.exp(z),root['section_theta']);history.append(q)
            if q['status']!='NUMERICAL_ONLY':break
            f=dec(q['L'])
            if abs(f)<mp.mpf('1e-22'):break
            if f*fa>0:a,fa=z,f
            else:b,fb=z,f
            zn=z-f/dec(q['L_z']);z=zn if a<zn<b else (a+b)/2
        out=dict(label=miss['label'],kind=miss['kind'],root_index=root['index'],r=q.get('r'),D=q.get('L'),stability=root['stability'],theta=root['section_theta'],history=history)
        save('refined_coordinates.jsonl',out);print(out['label'],out['kind'],out['root_index'],out['r'],out['D'],flush=True)
