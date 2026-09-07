#!/usr/bin/env python3
"""Fold-preserving predictor-corrector continuation with a full cycle inventory.

Constraint surface:  D(s_f,mu) = 0,  D_s(s_f,mu) = 0,  D_ss(s_f,mu) != 0.
Tangent identities (used and VERIFIED, not assumed):
    b_f . dmu = 0,          b_f = d_mu D(s_f,mu)
    ds_f = -(d_mu D_s . dmu)/D_ss(s_f)
For a separate nondegenerate stationary point s_i with critical value
h_i = D(s_i,mu):  dh_i = b_i . dmu   (envelope; D_s(s_i)=0).
Steering direction that preserves the fold and decreases |h_i| to first order:
    w = b_i - (b_i.b_f)/(b_f.b_f) b_f ,   dmu = -sign(h_i) w
All parameter moves are in SCALED coordinates: dmu_j * max(1,|mu_j|).
"""
import mpmath as mp
PARAMS = ("a", "a20", "a11", "a01", "a10")

def scales(mu):
    return [max(mp.mpf(1), abs(mu[p])) for p in PARAMS]

def val(c, mu, s):
    c.a, c.a20 = mu["a"], mu["a20"]
    return c.val([mu["a11"], mu["a01"], mu["a10"]], s)

def grad_mu(c, mu, s, which=("D",)):
    """Central differences of the requested jet entries w.r.t. all 5 params,
    in SCALED coordinates."""
    sc = scales(mu); out = {k: [mp.mpf(0)]*5 for k in which}
    for j, p in enumerate(PARAMS):
        h = mp.mpf("1e-9")*sc[j]
        m1, m2 = dict(mu), dict(mu); m1[p] += h; m2[p] -= h
        r1, r2 = val(c, m1, s), val(c, m2, s)
        if r1["status"] != "OK" or r2["status"] != "OK": return None
        for k in which: out[k][j] = (r1[k] - r2[k])/(2*h)*sc[j]
    return out

def correct(c, mu, s, pidx, tol="1e-26", itmax=25):
    """Newton back onto (D,D_s)=0 in (s, mu[pidx]).  Returns (mu,s,r) or None."""
    mu = dict(mu); p = PARAMS[pidx]; sc = max(mp.mpf(1), abs(mu[p]))
    for _ in range(itmax):
        r = val(c, mu, s)
        if r["status"] != "OK": return None, None, r["status"]
        F = [r["D"], r["Dx"]]
        if max(abs(F[0]), abs(F[1])*abs(s-1)) < mp.mpf(tol): return mu, s, r
        h = mp.mpf("1e-9")*sc
        m1 = dict(mu); m1[p] += h
        r1 = val(c, m1, s)
        if r1["status"] != "OK": return None, None, r1["status"]
        J = [[(r1["D"]-r["D"])/h,  r["Dx"]],
             [(r1["Dx"]-r["Dx"])/h, r["Dxx"]]]
        det = J[0][0]*J[1][1] - J[0][1]*J[1][0]
        if det == 0: return None, None, "singular"
        d0 = ( F[0]*J[1][1] - F[1]*J[0][1])/det
        d1 = ( J[0][0]*F[1] - J[1][0]*F[0])/det
        mu[p] -= d0; s -= d1
    return None, None, "no-converge"

def inventory(c, mu, smax=mp.mpf("7.0"), step=mp.mpf("0.02"), lo=mp.mpf("1.002")):
    """All sign changes of D (cycles) and of D_s (stationary points) on the
    return domain, refined by bisection.  Starts at 1.002, not 1.02."""
    xs, Ds, Dx = [], [], []
    s = lo
    while s < smax:
        r = val(c, mu, s)
        if r["status"] == "OK":
            xs.append(s); Ds.append(r["D"]); Dx.append(r["Dx"])
        elif xs:
            break
        s += step
    def refine(i, arr, key):
        loi, hii = xs[i], xs[i+1]
        for _ in range(60):
            mid = (loi+hii)/2; rm = val(c, mu, mid)
            if rm["status"] != "OK": break
            if rm[key]*arr[i] > 0: loi = mid
            else: hii = mid
        return (loi+hii)/2
    cyc = [refine(i, Ds, "D") for i in range(len(xs)-1) if Ds[i]*Ds[i+1] < 0]
    sta = [refine(i, Dx, "Dx") for i in range(len(xs)-1) if Dx[i]*Dx[i+1] < 0]
    return dict(domain=(xs[0], xs[-1]) if xs else None, cycles=cyc, stationary=sta)

def remote_cycles(mu, xlo=-40.0, xhi=-0.02, n=900):
    """Cycles in the OTHER nest: integrate the planar field and count crossings
    of the ray through the second focus.  Independent of the cusp engine."""
    import numpy as np
    from scipy.integrate import solve_ivp
    a, a20 = float(mu["a"]), float(mu["a20"])
    a11, a01, a10 = float(mu["a11"]), float(mu["a01"]), float(mu["a10"])
    a00 = a01 + a11 - a10 - a20 - a
    # second focus: equilibria on y=-1/x
    co = [a20, a10, a00 - a11, -a01, a]
    rs = np.roots(co) if abs(a20) > 1e-14 else np.roots(co[1:])
    foci = []
    for rr in rs:
        if abs(rr.imag) > 1e-9 or abs(rr.real - 1.0) < 1e-6 or abs(rr.real) < 1e-12: continue
        x = rr.real; y = -1.0/x
        J = np.array([[y, x],[a10 + 2*a20*x + a11*y, a01 + a11*x + 2*a*y]])
        det, tr = np.linalg.det(J), np.trace(J)
        if det > 0 and tr*tr - 4*det < 0: foci.append((x, y, tr))
    if not foci: return dict(second_focus=None, cycles=None)
    fx, fy, ftr = foci[0]
    def F(t, z):
        x, y = z
        return [1 + x*y, a00 + a10*x + a20*x*x + a01*y + a11*x*y + a*y*y]
    # displacement on the horizontal ray from the focus, x < fx
    def disp(d):
        z0 = [fx - d, fy]
        def ev(t, z): return z[1] - fy
        ev.direction = 1
        sol = solve_ivp(F, (0, 400), z0, rtol=1e-11, atol=1e-13, events=ev, max_step=0.5)
        for tt, p in zip(sol.t_events[0], sol.y_events[0]):
            if tt > 1e-6 and p[0] < fx: return (fx - p[0]) - d
        return None
    ds = [0.02*k for k in range(1, 60)]
    vals = [(d, disp(d)) for d in ds]
    vals = [(d, v) for d, v in vals if v is not None]
    n_cyc = sum(1 for i in range(len(vals)-1) if vals[i][1]*vals[i+1][1] < 0)
    return dict(second_focus=(fx, fy, ftr), cycles=n_cyc,
                samples=len(vals))
