#!/usr/bin/env python3
"""F11: order-2 weak focus with a homoclinic saddle loop around it; hunt for zero saddle divergence.
Shi chart at lam=0: x' = -y + l x^2 + m x y + y^2, y' = x + a x^2 + b x y, with eta_1 = 0 <=> m = a(b+2l)/(l+1).
eta_2 = a(b+2l)(b-3l-5)(a^2(b+2l+1)-(b+1)(l+1)^2)/(48(l+1)^2).
For each finite saddle, the separatrix splitting sigma (sine of the angle between the stable direction and the
returning unstable branch on a small circle) vanishes at a homoclinic loop. Along the loop set we record the
saddle trace (divergence) and eta_2. Usage: python3 fable_f11_neutral_loop.py OUT.jsonl [NPROC]
"""
import sys, os, json, numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), 'fable_engine'))
import retmap as rm
from scipy.integrate import solve_ivp
from scipy.optimize import brentq
from multiprocessing import Pool
RHO = 0.05
def mval(l, a, b): return a*(b+2*l)/(l+1)
def eta2(l, a, b): return a*(b+2*l)*(b-3*l-5)*(a*a*(b+2*l+1)-(b+1)*(l+1)**2)/(48*(l+1)**2)
def coef(l, a, b): return rm.shi_coef(0.0, l, mval(l, a, b), a, b, 1.0)
def saddles(c):
    out = []
    for (x, y) in rm.equilibria(c):
        if abs(x)+abs(y) < 1e-9: continue
        J = rm.jac(c, x, y)
        if np.linalg.det(J) < 0: out.append((np.array([x, y]), J))
    return out
def _cross_opposite(f, p0, pt, forward):
    """distance from the origin at which a separatrix first crosses the ray opposite to the saddle
    (after winding >= 0.4 turn about the origin), or None."""
    d = pt/np.linalg.norm(pt); n = np.array([-d[1], d[0]])
    def g(t, u): return u[0]*n[0]+u[1]*n[1]
    def far(t, u): return 60.0-np.hypot(*u)
    far.terminal = True
    sol = solve_ivp(f, (0, 400.0 if forward else -400.0), p0, rtol=1e-9, atol=1e-12, max_step=0.05, events=[g, far])
    if len(sol.t_events[0]) == 0: return None
    ang = np.unwrap(np.arctan2(sol.y[1], sol.y[0])); wind = np.abs(ang-ang[0])/(2*np.pi)
    for te, ye in zip(sol.t_events[0], sol.y_events[0]):
        if abs(te) < 1e-6: continue
        wv = np.interp(abs(te), np.abs(sol.t), wind)
        if wv >= 0.4 and np.dot(ye, d) < 0:
            return float(np.hypot(*ye))
    return None
def sigma_for(c, pt, J):
    """s_u - s_s on the ray opposite the saddle; zero at a homoclinic loop around the origin."""
    f = lambda t, u: [c[0]+c[1]*u[0]+c[2]*u[1]+c[3]*u[0]**2+c[4]*u[0]*u[1]+c[5]*u[1]**2,
                      c[6]+c[7]*u[0]+c[8]*u[1]+c[9]*u[0]**2+c[10]*u[0]*u[1]+c[11]*u[1]**2]
    w, v = np.linalg.eig(J); w = w.real; v = v.real
    vu = v[:, np.argmax(w)]/np.linalg.norm(v[:, np.argmax(w)]); vs = v[:, np.argmin(w)]/np.linalg.norm(v[:, np.argmin(w)])
    su = [x for x in (_cross_opposite(f, pt+1e-7*s*vu, pt, True) for s in (+1, -1)) if x is not None]
    ss = [x for x in (_cross_opposite(f, pt+1e-7*s*vs, pt, False) for s in (+1, -1)) if x is not None]
    if len(su) != 1 or len(ss) != 1: return None
    return su[0]-ss[0]
def sig(l, a, b, which):
    """splitting of saddle number `which` (ordered by x) or None."""
    c = coef(l, a, b); sd = sorted(saddles(c), key=lambda s: s[0][0])
    if which >= len(sd): return None
    return sigma_for(c, *sd[which])
def scan_line(args):
    a, b, ls = args
    out = []
    for which in range(3):
        vals = [sig(l, a, b, which) if abs(l+1) > 0.05 else None for l in ls]
        for i in range(len(ls)-1):
            v0, v1 = vals[i], vals[i+1]
            if v0 is None or v1 is None or v0*v1 >= 0: continue
            try:
                lstar = brentq(lambda l: sig(l, a, b, which), ls[i], ls[i+1], xtol=1e-10)
            except Exception:
                continue
            sv = sig(lstar, a, b, which)
            if sv is None or abs(sv) > 1e-5: continue
            c = coef(lstar, a, b); sd = sorted(saddles(c), key=lambda s: s[0][0])
            if which >= len(sd): continue
            pt, J = sd[which]
            out.append(dict(a=a, b=b, l=lstar, m=mval(lstar, a, b), which=which, saddle=pt.tolist(), trace=float(np.trace(J)),
                            eta2=float(eta2(lstar, a, b)), sig_lo=v0, sig_hi=v1))
    return out
if __name__ == "__main__":
    out = sys.argv[1]; NP = int(sys.argv[2]) if len(sys.argv) > 2 else 3
    ls = np.linspace(-4, 4, 25)
    jobs = [(a, b, ls) for a in (-3, -2, -1.5, -1, -0.5, 0.5, 1, 1.5, 2, 3) for b in np.linspace(-4, 4, 9)]
    n = 0
    with Pool(NP) as pool, open(out, 'a') as f:
        for res in pool.imap_unordered(scan_line, jobs):
            for r in res:
                f.write(json.dumps(r)+"\n"); f.flush(); n += 1
                print(f"loop: a={r['a']:.3f} b={r['b']:.3f} l={r['l']:.6f} saddle={np.round(r['saddle'],4).tolist()} trace={r['trace']:+.4f} eta2={r['eta2']:+.4e}", flush=True)
    print("DONE loops found:", n)
