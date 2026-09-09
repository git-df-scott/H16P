import json, numpy as np
from scipy.integrate import solve_ivp
R = json.load(open("ledger_opus/unfold.json"))[0]
mu = R["mu"]; a=float(mu["a"]); a20=float(mu["a20"]); a11=float(mu["a11"]); a01=float(mu["a01"]); a10=float(mu["a10"])
a00 = a01+a11-a10-a20-a
def F(t,z):
    x,y=z; return [1+x*y, a00+a10*x+a20*x*x+a01*y+a11*x*y+a*y*y]
print("engine cycles at s-1 =", [float(x)-1 for x in R["local"]])
for direction in (-1, +1):
    out=[]
    for d in (0.05, 0.072124, 0.15, 1.0, 2.42, 2.6):
        def ev(t,z): return z[1]+1.0
        ev.direction=direction
        sol=solve_ivp(F,(0,200),[1.0+d,-1.0],rtol=1e-10,atol=1e-12,events=ev,max_step=0.05)
        hit=None
        for tt,p in zip(sol.t_events[0], sol.y_events[0]):
            if tt>1e-6 and p[0]>1.0: hit=(tt,p[0]-1.0-d); break
        out.append((d, None if hit is None else (round(hit[0],3), "%.3e"%hit[1])))
    print("direction=%+d:" % direction, out)
