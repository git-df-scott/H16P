"""Independent high-precision first-order control at a=-2,b=1/3.

Direct reciprocal-height quadrature, not inherited log-height Gauss code.
This is point arithmetic, NOT interval validation or finite-field evidence.
Run with Python and mpmath. Prints JSON; no inherited files are modified.
"""
import json
import mpmath as mp

def evaluate(side, h_string, dps):
    with mp.workdps(dps):
        h=mp.mpf(h_string.split('/')[0])/mp.mpf(h_string.split('/')[1]) if '/' in h_string else mp.mpf(h_string)
        m=mp.mpf(8)/9; c=mp.mpf(159)/80
        rc=mp.mpf(2) if side==1 else mp.mpf(2)/5
        def potential(r):
            return 5*r*r/24-side*2*r/3-mp.log(r)/3
        assert h>potential(rc)
        def f(z):return potential(mp.exp(z))-h
        zc=mp.log(rc);zl=zc-1;zr=zc+1
        while f(zl)<0:zl-=1
        while f(zr)<0:zr+=1
        def bisect(l,r,left_positive):
            for _ in range(4*dps+30):
                mid=(l+r)/2
                if (f(mid)>0)==left_positive:l=mid
                else:r=mid
            return mp.exp((l+r)/2)
        rl=bisect(zl,zc,True);rr=bisect(zc,zr,False)
        # Direct r substitution smooths both square-root endpoints.
        def terms(t,which):
            r=rl+(rr-rl)*mp.sin(t)**2
            w=h-potential(r)
            if w<0:
                assert abs(w)<mp.power(10,-dps+8)
                w=mp.mpf(0)
            weight=mp.sqrt(w)*2*(rr-rl)*mp.sin(t)*mp.cos(t)
            return weight*(1 if which==0 else r if which==1 else w/(3*r))
        breaks=[0,mp.pi/8,mp.pi/4,3*mp.pi/8,mp.pi/2]
        integrals=[mp.quad(lambda t:terms(t,k),breaks) for k in range(3)]
        value=side*(integrals[1]-m*integrals[2])/integrals[0]-c
        return {'side':side,'h':h_string,'dps':dps,'value':mp.nstr(value,dps-8),
                'r_turning_points':[mp.nstr(rl,20),mp.nstr(rr,20)]}

if __name__=='__main__':
    points=[(1,h) for h in ['-1/6','-1/12','3/5','4/5','1','4/3']]
    points +=[(-1,h) for h in ['5','11/2']]
    results=[evaluate(side,h,dps) for dps in (50,80) for side,h in points]
    print(json.dumps({'scope':'Nonvalidated independent quadrature; no finite-field cycle claim',
        'mpmath_version':mp.__version__,'m':'8/9','c':'159/80',
        'results':results},indent=2))
