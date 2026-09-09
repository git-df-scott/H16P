"""Independent high-precision replay of accepted fields and outer witnesses."""
from pathlib import Path
import sys,json,time
import numpy as np
import mpmath as mp
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent/'moving_cycle_2026_09_09'))
import mp_replay
import engine
mp.mp.dps=65

for name in sys.argv[1:]:
    start=time.perf_counter()
    result=json.loads((HERE/f'{name}_corrected_result.json').read_text());s=result['states'][-1]
    cf=engine.coefficients(s['q']);exact=[repr(float(v)) for v in cf];pars=[mp.mpf(v) for v in exact]
    out={'name':name,'evidence':'NUM, independent Cartesian Taylor at 65 digits, no interval claim',
         'exact_decimal_coefficients':exact,'witnesses':[],'dps':mp.mp.dps,'estimated_tail_tolerance':'1e-32'}
    for w in s['witnesses']:
        rows=[mp_replay.pair(w['side']*mp.exp(mp.mpf(repr(float(z)))),pars,tol='1e-32') for z in w['interval']]
        passed=bool(mp.mpf(rows[0]['log_D'])*mp.mpf(rows[1]['log_D'])<0)
        out['witnesses'].append({'side':w['side'],'interval':w['interval'],'replays':rows,'opposite_signs':passed})
        (HERE/f'{name}_independent_replay.json').write_text(json.dumps(out,indent=2)+'\n')
        print(name,'bracket',w['side'],passed,[r['log_D'] for r in rows],flush=True)
    z=s['outer']['s'];out['outer']=mp_replay.pair(mp.exp(mp.mpf(repr(float(z)))),pars,tol='1e-32')
    out['outer']['s']=z;out['saved_outer_D']=s['outer']['D']
    out['all_four_sign_brackets']=all(w['opposite_signs'] for w in out['witnesses'])
    out['wall_seconds']=time.perf_counter()-start
    (HERE/f'{name}_independent_replay.json').write_text(json.dumps(out,indent=2)+'\n')
    print(name,'outer',out['outer']['log_D'],'four_brackets',out['all_four_sign_brackets'],flush=True)
