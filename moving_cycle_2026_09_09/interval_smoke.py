"""Experimental interval-Taylor first-return verifier for one seed cycle.

Uses outward interval arithmetic, Picard tube inclusion, and a Lagrange
remainder coefficient bounded over the entire tube. Rejects any step whose
section crossing order cannot be established. This is a prototype smoke test,
not independently audited certification software.
"""
import json,time,sys
from pathlib import Path
import mpmath as mp
from mpmath import iv
import mp_replay
ROOT=Path(__file__).resolve().parent
iv.dps=int(sys.argv[1]) if len(sys.argv)>1 else 80
mp.mp.dps=iv.dps+15
DEGREE=int(sys.argv[2]) if len(sys.argv)>2 else 32
OUTPUT='interval_smoke_replay' if len(sys.argv)>1 else 'interval_smoke'
STATS={'steps':0,'picard_retries':0,'event_steps':0}

def low(x):return mp.make_mpf(x._mpi_[0])
def high(x):return mp.make_mpf(x._mpi_[1])
def ab(x):return max(abs(low(x)),abs(high(x)))
def enc(x):
    # Directed decimal widening for standalone logs; all proof decisions use
    # the original binary interval endpoints, not these printed strings.
    digits=iv.dps+4
    def outward(v,sign):
        if not v:return '0'
        quantum=mp.power(10,mp.floor(mp.log10(abs(v)))-digits+1)
        widened=(mp.floor(v/quantum) if sign<0 else mp.ceil(v/quantum))*quantum
        return mp.nstr(widened,digits+2)
    return [outward(low(x),-1),outward(high(x),1)]
def exact(x):return iv.mpf(str(x))
def interval(lo,hi):return iv.mpf([str(lo),str(hi)])
def subset(x,y):return low(x)>low(y) and high(x)<high(y)
PARS=[-iv.mpf(7)/4,iv.mpf(1)/3,-iv.mpf(1)/27500,-iv.mpf(31379)/250000000,-iv.mpf(7517)/50000000]

def vector(x,y):
    a,b,e0,e1,e2=PARS
    return ((b-2)/4+(1-b)*y+a*x*x+b*y*y+e1*x+e2*x*y,e0-2*x*y)

def coeffs(x,y,nmax):
    a,b,e0,e1,e2=PARS;xs=[x];ys=[y];ds=[iv.mpf(0)]
    for n in range(nmax):
        xx=sum((xs[k]*xs[n-k] for k in range(n+1)),iv.mpf(0))
        yy=sum((ys[k]*ys[n-k] for k in range(n+1)),iv.mpf(0))
        xy=sum((xs[k]*ys[n-k] for k in range(n+1)),iv.mpf(0))
        xs.append(((b-2)/4*(n==0)+(1-b)*ys[n]+a*xx+b*yy+e1*xs[n]+e2*xy)/(n+1))
        ys.append((e0*(n==0)-2*xy)/(n+1))
        ds.append((2*(a-1)*xs[n]+e1*(n==0)+e2*ys[n])/(n+1))
    return xs,ys,ds

def horner(cs,t):
    v=cs[-1]
    for c in cs[-2::-1]:v=v*t+c
    return v

def step(x,y,h):
    # Tube B encloses every initial state over [0,h]. Picard inclusion implies
    # existence within B by the standard first-exit argument.
    F=vector(x,y);hh=exact(h);tbox=interval(0,h)
    B=[]
    for X,f in zip([x,y],F):
        rad=mp.mpf('1.5')*h*ab(f)+mp.mpf('1e-70')
        B.append(X+interval(-rad,rad))
    good=False
    for _ in range(12):
        FB=vector(*B);im=[x+tbox*FB[0],y+tbox*FB[1]]
        if all(subset(z,b) for z,b in zip(im,B)):
            good=True;break
        enlarged=[]
        for z,bound in zip(im,B):
            if subset(z,bound):
                enlarged.append(bound)
            else:
                lo=min(low(z),low(bound));hi=max(high(z),high(bound));pad=mp.mpf('.2')*(hi-lo)
                enlarged.append(interval(lo-pad,hi+pad))
        B=enlarged
    if not good:raise RuntimeError('Picard inclusion failed')
    cs=coeffs(x,y,DEGREE);rem=coeffs(B[0],B[1],DEGREE+1)
    def evaluate(t):
        return [horner(c,t)+r[-1]*t**(DEGREE+1) for c,r in zip(cs,rem)]
    return evaluate,evaluate(hh),B,FB

def half_return(y0,initial_sign):
    x=iv.mpf(0);y=y0;di=iv.mpf(0);timebox=iv.mpf(0);first=True
    for j in range(8000):
        scale=max(1,ab(x),ab(y));h=mp.mpf(1)/64
        while h*scale>mp.mpf('0.025'):h/=2
        for retry in range(15):
            try:
                ev,end,B,FB=step(x,y,h)
                # Where the whole tube touches x=0, a strictly signed normal
                # component ensures that at most one crossing can occur.
                touches=low(B[0])<=0<=high(B[0])
                transverse=low(FB[0])>0 or high(FB[0])<0
                if touches and not transverse:raise RuntimeError('ambiguous section gate')
                if first and ((initial_sign==1 and low(FB[0])<=0) or (initial_sign==-1 and high(FB[0])>=0)):
                    raise RuntimeError('initial departure not transverse')
                break
            except RuntimeError as error:
                last_error=str(error)
                STATS['picard_retries']+=1;h/=2
        else:raise RuntimeError('step gate failed persistently: '+last_error+' state='+str([enc(x),enc(y)])+' h='+str(h))
        xe,ye,de=end;STATS['steps']+=1
        crossed=(high(xe)<0 if initial_sign==1 else low(xe)>0)
        if crossed:
            if not touches or not transverse:raise RuntimeError('crossing proof missing')
            lo=mp.mpf(0);hi=h
            # Simultaneous event enclosure for every initial point in y0.
            for _ in range(200):
                mid=(lo+hi)/2;xm=ev(exact(mid))[0]
                if (low(xm)>0 if initial_sign==1 else high(xm)<0):lo=mid
                elif (high(xm)<0 if initial_sign==1 else low(xm)>0):hi=mid
                else:break
            ti=interval(lo,hi);ee=ev(ti)
            STATS['event_steps']+=1
            return ee[1],di+ee[2],timebox+ti,{'transverse':enc(FB[0]),'event_time_width':mp.nstr(hi-lo,12),'steps':j+1}
        # A wide endpoint straddling the section cannot be advanced silently.
        if low(xe)<=0<=high(xe):raise RuntimeError('endpoint straddles section before complete event bracket')
        x,y=xe,ye;di+=de;timebox+=exact(h);first=False
        if ab(x)>mp.mpf('1e8') or ab(y)>mp.mpf('1e8'):raise RuntimeError('state guard')
    raise RuntimeError('step count guard')

def full(y0):
    mid,d1,t1,g1=half_return(y0,1)
    end,d2,t2,g2=half_return(mid,-1)
    p0=vector(iv.mpf(0),y0)[0];p1=vector(iv.mpf(0),end)[0]
    deriv=p0/p1*iv.exp(d1+d2)
    return {'start':enc(y0),'return':enc(end),'displacement':enc(end-y0),
            'derivative':enc(deriv),'time':enc(t1+t2),'half_gates':[g1,g2]},deriv

def main():
    start=time.perf_counter()
    mp.mp.dps=iv.dps+15
    # Refine the numerical half-mismatch root before choosing a tiny rational
    # section interval. This preparatory computation is not the certificate.
    pars=mp_replay.pars_seed();y=mp.mpf('2.6651026')
    cached=ROOT/'interval_smoke.json'
    if cached.exists():y=mp.mpf(json.loads(cached.read_text())['center'])
    for _ in range(0 if cached.exists() else 4):
        f=mp_replay.half(y,pars,1,tol='1e-60',degree=40)
        b=mp_replay.half(y,pars,-1,tol='1e-60',degree=40)
        D=f['y']-b['y'];h=mp.mpf('1e-10')
        fp=mp_replay.half(y+h,pars,1,tol='1e-55',degree=40)['y']-mp_replay.half(y+h,pars,-1,tol='1e-55',degree=40)['y']
        fm=mp_replay.half(y-h,pars,1,tol='1e-55',degree=40)['y']-mp_replay.half(y-h,pars,-1,tol='1e-55',degree=40)['y']
        y-=D/((fp-fm)/(2*h))
    center=mp.mpf(mp.nstr(y,65));radius=mp.mpf('1e-20')
    out={'scope':'Prototype interval smoke test for ONE known seed cycle, not a five-cycle claim',
         'interval_dps':iv.dps,'degree':DEGREE,'center':mp.nstr(center,70),'radius':str(radius),'records':[]}
    print('center',out['center'],flush=True)
    try:
        for label,startpoint in [('left',exact(center-radius)),('right',exact(center+radius)),('whole_interval',interval(center-radius,center+radius))]:
            row,der=full(startpoint);row['label']=label;out['records'].append(row)
            print(label,row,flush=True)
            (ROOT/f'{OUTPUT}_checkpoint.json').write_text(json.dumps(out,indent=2)+'\n')
        left,right,whole=out['records']
        opposite=mp.mpf(left['displacement'][1])<0<mp.mpf(right['displacement'][0]) or mp.mpf(right['displacement'][1])<0<mp.mpf(left['displacement'][0])
        der=whole['derivative'];isolated=mp.mpf(der[1])<1 or mp.mpf(der[0])>1
        out['status']='passed_prototype_gates' if opposite and isolated else 'unresolved_enclosures'
        out['strict_endpoint_signs']=opposite;out['derivative_excludes_one']=isolated
    except Exception as e:
        out['status']='unresolved';out['error']=repr(e)
    out['stats']=STATS;out['wall_seconds']=time.perf_counter()-start
    (ROOT/f'{OUTPUT}.json').write_text(json.dumps(out,indent=2)+'\n')
    print('STOP',out['status'],out.get('error'),STATS,flush=True)

if __name__=='__main__':main()
