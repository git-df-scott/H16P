"""F17: second-order Melnikov function at the holomorphic two-center point (a,b)=(-1,1) along the
first-order-null direction x' += eps*x*(1-2y) (e1=1, e2=-2), plus optional admixture.  D(y;eps) on the
transverse section (height y above the invariant line) for the upper annulus; check D/eps -> 0 and extract
M2 = lim D/eps^2; count its zeros.  Same for the lower annulus (ray up from the lower center)."""
import numpy as np, retmap as rm, sys
def field(e0, e1, e2, a=-1.0, b=1.0): return np.array([(b-2)/4, e1, 1-b, a, e2, b, e0, 0, 0, 0, -2.0, 0], float)
Y = np.geomspace(1e-6, 0.45, 120)
def D_upper(c):
    eq = rm.equilibria(c); pt = min([p for p in eq if p[1] > 0], key=lambda p: abs(p[1]-0.5))
    U0 = np.log(pt[1]-Y); u1, S, st = rm.returns_log(c[None], np.array([pt]), U0[None], th0=-np.pi/2, rtol=1e-13, umax=60, Smax=5000, maxsteps=3_000_000)
    return np.where(st[0] == 0, (pt[1]-np.exp(u1[0]))-Y, np.nan)
def D_lower(c):
    eq = rm.equilibria(c); pt = min([p for p in eq if p[1] < 0], key=lambda p: abs(p[1]+0.5))
    U0 = np.log(-pt[1]-Y); u1, S, st = rm.returns_log(c[None], np.array([pt]), U0[None], th0=np.pi/2, rtol=1e-13, umax=60, Smax=5000, maxsteps=3_000_000)
    return np.where(st[0] == 0, (-pt[1]-np.exp(u1[0]))-Y, np.nan)
for name, Dfun in (("upper", D_upper), ("lower", D_lower)):
    print("==", name)
    for eps in (1e-3, 5e-4, 2.5e-4):
        D = Dfun(field(0.0, eps, -2*eps)); m = np.isfinite(D)
        print(f" eps={eps:g}: max|D/eps| = {np.nanmax(np.abs(D/eps)):.3e}   max|D/eps^2| = {np.nanmax(np.abs(D/eps**2)):.3e}  valid y >= {Y[m].min():.1e}")
        M2 = D/eps**2; sig = m & (np.abs(M2) > 1e-4*np.nanmax(np.abs(M2)))
        zi = [i for i in range(len(Y)-1) if sig[i] and sig[i+1] and M2[i]*M2[i+1] < 0]
        print(f"   M2 zeros (y): {[f'{Y[i]:.3e}' for i in zi]}   M2 at y=1e-5,1e-3,0.1,0.4: {[f'{M2[np.argmin(abs(Y-v))]:+.3e}' for v in (1e-5,1e-3,0.1,0.4)]}")
    np.save(f'data/f17_M2_{name}.npy', np.vstack([Y, D/eps**2]))
