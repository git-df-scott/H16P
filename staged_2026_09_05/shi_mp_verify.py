#!/usr/bin/env python3
"""Independent 40-digit polar-return check with midpoint extrapolation.

No use of the exact Lyapunov polynomial or scipy orbit. Nonvalidated numerics.
Each complete return, successful or failed, shares the 160-return ledger cap.
"""
import json,signal,time
import mpmath as mp
from shi_trace import LEDGER,CAP

def verify(params,r,lam,blocks=32,levels=7,tag='mp-independent'):
    count=sum(1 for _ in LEDGER.open()) if LEDGER.exists() else 0
    if count>=CAP:raise RuntimeError('160-return budget exhausted')
    record={'index':count+1,'family':params['name'],'method':'mp-midpoint-extrapolation',
            'l_m_a_b_exact':params['exact'],'r':str(r),'lambda':str(lam),
            'blocks':blocks,'levels':levels,'dps':40,'tag':tag,'validated':False}
    begin=time.process_time()
    def timeout(*args):raise TimeoutError('10 CPU second evaluation cap')
    old=signal.signal(signal.SIGPROF,timeout);signal.setitimer(signal.ITIMER_PROF,10)
    try:
        with mp.workdps(40):
            def rat(s):
                z=s.split('/');return mp.mpf(z[0])/(mp.mpf(z[1]) if len(z)>1 else 1)
            l,m,a,b=map(rat,params['exact']);r=mp.mpf(str(r));lam=mp.mpf(str(lam))
            width=2*mp.pi/blocks; q=mp.mpf(0)
            for block in range(blocks):
                start=block*width; rows=[]
                for k in range(1,levels+1):
                    n=2*k;h=width/n
                    def f(j,z):
                        theta=start+j*h;c=mp.cos(theta);s=mp.sin(theta);rr=r*mp.exp(z)
                        A=l*c**3+(m+a)*c*c*s+(1+b)*c*s*s
                        B=a*c**3+(b-l)*c*c*s-m*c*s*s-s**3
                        w=1-lam*c*s+rr*B
                        if w<=mp.mpf('.03'):raise ValueError('polar chart lost')
                        return (lam*c*c+rr*A)/w
                    prev=q;cur=q+h*f(0,q)
                    for j in range(1,n):prev,cur=cur,prev+2*h*f(j,cur)
                    row=[(prev+cur+h*f(n,cur))/2]
                    for j in range(1,k):
                        ratio=mp.mpf(n)/(2*(k-j))
                        row.append(row[-1]+(row[-1]-rows[k-2][j-1])/(ratio**2-1))
                    rows.append(row)
                q=rows[-1][-1]
            record.update(status='ok',raw_log_return=mp.nstr(q,35))
    except Exception as e:record.update(status='failed',error=str(e))
    finally:
        signal.setitimer(signal.ITIMER_PROF,0);signal.signal(signal.SIGPROF,old)
        record['cpu_seconds']=time.process_time()-begin
        with LEDGER.open('a') as f:f.write(json.dumps(record)+'\n')
    return record

if __name__=='__main__':
    p={'name':'Shi conditioned','exact':['-10','499/100','1','-3113751/125000']}
    for r in ['.00004','.00015','.0015','.025']:
        print(json.dumps(verify(p,r,'-1e-14')),flush=True)
