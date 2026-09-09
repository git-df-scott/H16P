"""Independent arbitrary-precision Cartesian Taylor integrator.

Quadratic coefficient recurrences are evaluated directly. Tail estimates are
adaptive numerical estimates, NOT rigorous interval remainder bounds.
"""
import json,time
from pathlib import Path
import mpmath as mp
ROOT=Path(__file__).resolve().parent

def pars_seed():
    return [mp.mpf(-7)/4,mp.mpf(1)/3,-mp.mpf(1)/27500,
            -mp.mpf(31379)/250000000,-mp.mpf(7517)/50000000]

def half(y0,pars,sign,tol='1e-32',degree=28):
    a,b,e0,e1,e2=pars;tol=mp.mpf(tol);x=mp.mpf(0);y=mp.mpf(y0)
    tim=mp.mpf(0);di=mp.mpf(0);started=False
    for step in range(20000):
        xs=[x];ys=[y];ds=[mp.mpf(0)]
        for n in range(degree):
            xx=sum(xs[k]*xs[n-k] for k in range(n+1))
            yy=sum(ys[k]*ys[n-k] for k in range(n+1))
            xy=sum(xs[k]*ys[n-k] for k in range(n+1))
            xs.append(sign*((b-2)/4*(n==0)+(1-b)*ys[n]+a*xx+b*yy+e1*xs[n]+e2*xy)/(n+1))
            ys.append(sign*(e0*(n==0)-2*xy)/(n+1))
            ds.append(sign*(2*(a-1)*xs[n]+e1*(n==0)+e2*ys[n])/(n+1))
        scale=max(1,abs(x),abs(y))
        hh=[mp.mpf('.08')]
        for n in [degree-1,degree]:
            cn=max(abs(xs[n]),abs(ys[n]))
            if cn:hh.append(mp.mpf('.65')*(tol*scale/cn)**(mp.mpf(1)/n))
        h=min(hh)
        def ev(cs,t):return mp.polyval(cs[::-1],t)
        # Inspect the polynomial step at interior nodes for the first crossing.
        prev=mp.mpf(0);found=None
        for j in range(1,9):
            t=h*j/8;xt=ev(xs,t)
            if sign*xt>0:started=True
            if started and sign*xt<=0:
                found=(prev,t);break
            prev=t
        if found:
            lo,hi=found
            for _ in range(130):
                mid=(lo+hi)/2
                if sign*ev(xs,mid)>0:lo=mid
                else:hi=mid
            t=(lo+hi)/2
            return {'y':ev(ys,t),'div':di+ev(ds,t),'time':tim+t,'steps':step+1}
        x=ev(xs,h);y=ev(ys,h);di+=ev(ds,h);tim+=h
        if tim>100 or abs(x)>mp.mpf('1e40') or abs(y)>mp.mpf('1e40'):raise RuntimeError('guard')
    raise RuntimeError('step limit')

def pair(y,pars,tol='1e-32'):
    f=half(y,pars,1,tol);b=half(y,pars,-1,tol)
    return {'start_y':str(y),'D':mp.nstr(f['y']-b['y'],45),
       'log_D':mp.nstr(mp.log(abs(f['y']))-mp.log(abs(b['y'])),45),
       'log_mu':mp.nstr(f['div']-b['div'],45),
       'forward':{k:mp.nstr(v,45) if k!='steps' else v for k,v in f.items()},
       'backward':{k:mp.nstr(v,45) if k!='steps' else v for k,v in b.items()}}

if __name__=='__main__':
    mp.mp.dps=60;start=time.perf_counter();pars=pars_seed()
    src=json.loads((ROOT.parent/'reversible_reseed/data/verified_control.json').read_text())
    out={'evidence':'NUM, 60-digit direct Cartesian Taylor method, estimated tails only',
         'dps':mp.mp.dps,'tol':'1e-32','degree':28,'witnesses':[]}
    for r in src['moments']:
        row=pair(r['side']*mp.mpf(r['y_high_abs']),pars)
        out['witnesses'].append(row)
        (ROOT/'mp_seed_replay.json').write_text(json.dumps(out,indent=2)+'\n')
        print(json.dumps(row),flush=True)
    # Repeat the weakest upper witness with a tighter estimated tail tolerance.
    r=src['moments'][1]
    out['tighter_check']=pair(mp.mpf(r['y_high_abs']),pars,tol='1e-40')
    out['wall_seconds']=time.perf_counter()-start
    (ROOT/'mp_seed_replay.json').write_text(json.dumps(out,indent=2)+'\n')
    print('complete',out['wall_seconds'],flush=True)
