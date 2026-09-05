"""Separate double coefficient rounding from trajectory integration error."""
import json,subprocess
from check_misses import HERE,D1,dec,st,Q,save
import mpmath as mp
rows=[json.loads(s) for s in (D1/'fields.jsonl').read_text().splitlines()]
misses=[json.loads(s) for s in (HERE/'missed_roots.jsonl').read_text().splitlines()]
done={(r['label'],r['u']) for r in map(json.loads,(HERE/'coefficient_rounding.jsonl').read_text().splitlines())} if (HERE/'coefficient_rounding.jsonl').exists() else set()
for miss in misses:
    if miss['label'] not in ['center_0.0000000001','positive_extension_1','negative_extension_1'] or miss['kind']!='pair':continue
    row=next(r for r in rows if r['label']==miss['label'] and r['kind']==miss['kind'])
    v=[Q.from_float(float(Q(x))) for x in row['coefficient_vector']]
    for g in miss['roots'][0]['grid']:
        if (miss['label'],g['u_decimal']) in done:continue
        u=dec(g['u_decimal']);theta=dec(miss['roots'][0]['section_theta'])
        args=[st(mp.exp(u)),st(dec(str(v[11]))),st(dec(str(-v[7]))),st(dec(str(v[8]))),'2e-26',st(theta),'0',st(dec(str(v[10])))]
        try:
            p=subprocess.run(['/tmp/d1_full_ray_float'],input=' '.join(args)+'\n',text=True,capture_output=True,timeout=55)
            out=json.loads(p.stdout) if p.returncode==0 else dict(status='UNRESOLVED',error='exit '+str(p.returncode))
        except Exception as e:out=dict(status='UNRESOLVED',error=str(e))
        result=dict(label=miss['label'],u=g['u_decimal'],theta=st(theta),exact_rational_coefficients=row['coefficient_vector'],double_coefficients_as_exact_rationals=list(map(str,v)),quad_on_rounded_field=out,quad_on_exact_field=g['quad'],native=g['native'])
        save('coefficient_rounding.jsonl',result)
        print(miss['label'],g['u_decimal'],out['status'],out.get('L'),flush=True)
