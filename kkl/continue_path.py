"""Follow the two known cycles on one rational constant-K path, serially."""
import argparse
from fractions import Fraction as F
import json
import math
from pathlib import Path
from geometry import gates
from pilot import call, LEDGER, HERE, used_steps

PATH=HERE/'data'/'constant_K_path.jsonl'


def count(): return sum(1 for _ in LEDGER.open())


def correct(r,c,alpha,label,limit):
    sign=1 if r>0 else -1
    log_r=math.log(abs(r))
    for iteration in range(8):
        if count()>=limit: raise RuntimeError('requested evaluation ceiling reached')
        r=sign*math.exp(log_r)
        if (sign>0 and not 2**-12<=r<=2**10) or (sign<0 and not -2**20<=r<-1):
            raise RuntimeError('section range boundary')
        result=call({'r':r,'c':str(c),'alpha':str(alpha),'beta':'0'},label+f'_corrector_{iteration}')
        if result['status']!='NUMERICAL_ONLY': raise RuntimeError(result.get('error','return failed'))
        R=result['return_coordinate']
        error=math.log(R/r)
        slope=r*result['R_r']/R-1
        if abs(result['derivative_discrepancy'])>1e-7*max(1,abs(result['R_r'])):
            raise RuntimeError('derivative cross-check failed')
        if abs(error)<2e-10:
            result['log_displacement']=error
            result['evaluation']=count()
            return result
        if abs(slope)<1e-5: raise RuntimeError('near fold: switch to augmented fold equations')
        correction=-error/slope
        if abs(correction)>.4: raise RuntimeError('corrector trust region exceeded')
        log_r+=correction
    raise RuntimeError('corrector iteration budget exhausted')


def predict(old,dc):
    c=old['c']
    da_dc=2376/(11*c-5)**2
    dr_dc=-(old['R_c']+old['R_alpha']*da_dc)/(old['R_r']-1)
    return math.copysign(math.exp(math.log(abs(old['r']))+dc*dr_dc/old['r']),old['r'])


if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--end',default='3/4')
    parser.add_argument('--steps',type=int,default=4)
    parser.add_argument('--evaluation-ceiling',type=int,default=64)
    args=parser.parse_args()
    rows=[json.loads(line) for line in LEDGER.read_text().splitlines()]
    if PATH.exists():
        previous=[json.loads(line) for line in PATH.read_text().splitlines() if json.loads(line)['status']=='ACCEPTED_NUMERICAL_POINT'][-1]
        start=F(previous['c']); origin=previous['origin']; remote=previous['remote']
    else:
        start=F(7,10)
        origin=next(row['result'] for row in rows if row['purpose']=='start_origin_control')
        remote=next(row['result'] for row in rows if row['purpose']=='start_remote_control')
    end=F(args.end)
    for j in range(1,args.steps+1):
        if used_steps()>=256: raise RuntimeError('total seeded continuation-step ceiling reached')
        attempted=sum(1 for _ in PATH.open()) if PATH.exists() else 0
        if attempted >= (16 if args.evaluation_ceiling<=64 else 256):
            raise RuntimeError('path continuation-step ceiling reached')
        c=start+(end-start)*F(j,args.steps)
        alpha=-F(216)/(11*c-5)
        label=f'constant_K_c_{c}'
        entry={'c':str(c),'alpha':str(alpha),'beta':'0','path':'K=6/5','status':'UNRESOLVED'}
        try:
            entry['geometry']=gates(c,alpha)
            if c>=F(241,250): raise RuntimeError('infinity chart boundary: preflight required')
            dc=float(c)-origin['c']
            new_origin=correct(predict(origin,dc),c,alpha,label+'_origin',args.evaluation_ceiling)
            new_remote=correct(predict(remote,dc),c,alpha,label+'_remote',args.evaluation_ceiling)
            assert new_origin['R_r']<1 and new_remote['R_r']>1
            assert abs(new_origin['winding_about_focus']+1)<1e-6
            assert abs(new_remote['winding_about_focus']-1)<1e-6
            entry.update(status='ACCEPTED_NUMERICAL_POINT',origin=new_origin,remote=new_remote)
            origin,remote=new_origin,new_remote
        except Exception as exc:
            entry['error']=str(exc)
        with PATH.open('a') as stream: stream.write(json.dumps(entry)+'\n')
        print(json.dumps({'path_event':entry['status'],'c':str(c),'alpha':str(alpha),
                         'origin_r':entry.get('origin',{}).get('r'),'origin_multiplier':entry.get('origin',{}).get('R_r'),
                         'remote_r':entry.get('remote',{}).get('r'),'remote_multiplier':entry.get('remote',{}).get('R_r'),
                         'error':entry.get('error')}),flush=True)
        if entry['status']!='ACCEPTED_NUMERICAL_POINT': break
