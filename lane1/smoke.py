import numpy as np, engine as E
def cherkas(a,a20,a11,a01,a10):
    a00 = a01 + a11 - a10 - a20 - a
    return np.array([1,0,0,0,1,0, a00,a10,a01,a20,a11,a],float)
rows = {
 1:(3,-12,-1.398,8.4,15.28),
 3:(-2,12,10.999,-14,-26.1),
 4:(-2,-1,9.49965,-12.5,6.955),
}
for rid,(a,a20,a11,a01,a10) in rows.items():
    v = cherkas(a,a20,a11,a01,a10)
    print("row",rid,"resid at A:",E.residual_at(v,(1,-1)))
    loc = E.local10(v,(1,-1))
    J = E.jac(v,(1,-1)); print("  trace",J[0,0]+J[1,1],"det",np.linalg.det(J))
    for phi,lab in ((0.0,'+x'),(np.pi,'-x')):
        s = np.linspace(0.02,3.5,180)
        D,st,T = E.d_curve(loc,phi,s)
        ok = st==0
        sc=[]
        for i in range(len(s)-1):
            if ok[i] and ok[i+1] and D[i]*D[i+1]<0: sc.append((s[i],s[i+1]))
        print(f"  phi={lab}: ok={ok.sum()}/{len(s)} smax~{s[ok].max() if ok.any() else None:.3f} signchanges={len(sc)}", [f"{0.5*(a1+b1):.3f}" for a1,b1 in sc])
