"""A specified segment of the seeded (c,K) return-root sheet, not a sweep."""
import argparse
from fractions import Fraction as F
import json
import math
from geometry import gates
from continue_path import correct
from pilot import HERE, used_steps

if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--source',required=True)
    parser.add_argument('--label',required=True)
    parser.add_argument('--c-end',required=True)
    parser.add_argument('--K-end',required=True)
    parser.add_argument('--steps',required=True,type=int)
    args=parser.parse_args()
    source=HERE/'data'/args.source
    output=HERE/'data'/(args.label+'.jsonl')
    rows=[json.loads(line) for line in source.read_text().splitlines()]
    old=next(row for row in reversed(rows) if row['status']=='ACCEPTED_NUMERICAL_POINT')
    c0=F(old['c']); a0=F(old['alpha']); K0=-a0*(F(11,5)*c0-1)-42
    c1=F(args.c_end); K1=F(args.K_end)
    assert 1<=args.steps<=256
    for j in range(1,args.steps+1):
        if used_steps()>=256: raise RuntimeError('total seeded continuation-step ceiling reached')
        c=c0+(c1-c0)*F(j,args.steps); K=K0+(K1-K0)*F(j,args.steps)
        alpha=-5*(K+42)/(11*c-5)
        entry={'c':str(c),'alpha':str(alpha),'beta':'0','K':str(K),'segment':args.label,'status':'UNRESOLVED'}
        try:
            entry['geometry']=gates(c,alpha)
            if c>=F(241,250): raise RuntimeError('infinity boundary: separate chart preflight required')
            for nest in ['origin','remote']:
                previous=old[nest]
                dc=float(c-F(old['c'])); da=float(alpha-F(old['alpha']))
                dlog=-(previous['R_c']*dc+previous['R_alpha']*da)/((previous['R_r']-1)*previous['r'])
                if abs(dlog)>.65: raise RuntimeError('predictor trust region: subdivide specified segment')
                r=math.copysign(math.exp(math.log(abs(previous['r']))+dlog),previous['r'])
                entry[nest]=correct(r,c,alpha,args.label+f'_{j}_{nest}',4096)
            assert entry['origin']['R_r']<1 and entry['remote']['R_r']>1
            assert abs(entry['origin']['winding_about_focus']+1)<1e-6
            assert abs(entry['remote']['winding_about_focus']-1)<1e-6
            entry['status']='ACCEPTED_NUMERICAL_POINT'
            old=entry
        except Exception as exc:
            entry['error']=str(exc)
        with output.open('a') as f: f.write(json.dumps(entry)+'\n')
        print(json.dumps({'event':entry['status'],'c':str(c),'K':str(K),
              'origin_r':entry.get('origin',{}).get('r'),'origin_multiplier':entry.get('origin',{}).get('R_r'),
              'remote_r':entry.get('remote',{}).get('r'),'remote_multiplier':entry.get('remote',{}).get('R_r'),
              'error':entry.get('error')}),flush=True)
        if entry['status']!='ACCEPTED_NUMERICAL_POINT': break
