from engine import *
import mpmath as mp
import mp_replay
import sys

mp.mp.dps=60
for kind in ['shape','amplitude']:
    prefix=kind+'_resumed' if '--resumed' in sys.argv else kind
    data=json.loads((ROOT/f'{prefix}_result.json').read_text());last=data['accepted'][-1]
    q=np.array(last['q']);cf=coefficients(q)
    # Explicit decimal strings define one fixed field for all its tests.
    exact=[repr(float(v)) for v in cf];pars=[mp.mpf(x) for x in exact]
    out={'kind':kind,'evidence':'NUM, independent Cartesian/Taylor replay',
         'exact_decimal_coefficients':exact,'root_replays':[],'witness_replays':[]}
    for i,s in enumerate(last['roots']):
        side=1 if i<3 else -1
        cr=cartesian_full(s,q,side,rtol=8e-14)
        p=pair(s,q,side,rtol=8e-14,regularized=True)
        out['root_replays'].append({'side':side,'s':s,'cartesian_full':cr,'logarithmic':p})
    for w in last['witnesses']:
        rr=[]
        for s in w['interval']:
            y=w['side']*mp.exp(mp.mpf(repr(float(s))))
            rr.append(mp_replay.pair(y,pars,tol='1e-30'))
        assert mp.mpf(rr[0]['log_D'])*mp.mpf(rr[1]['log_D'])<0
        out['witness_replays'].append({'interval':w['interval'],'side':w['side'],'replays':rr})
        save(f'{prefix}_independent_replay.json',out)
        print(kind,'witness',w['side'],w['interval'],[r['log_D'] for r in rr],flush=True)
    out['counts']=COUNTS.copy();save(f'{prefix}_independent_replay.json',out)
