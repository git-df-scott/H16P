"""Engine B -- independent second integrator (PROTOCOL rule 2).

Everything differs from engine A (cusp128.cpp):
  arithmetic  : mpmath mpf at dps 40-60      (A: binary128, 34 digits)
  coordinates : (u,v) = (x-1, y+1) Cartesian (A: polar chart round the focus)
  independent : time t                       (A: polar angle theta)
  method      : variable-order Taylor series (A: Gragg-Bulirsch-Stoer)
  section     : event v = 0, u > 0, found by Newton on the Taylor polynomial

  u' = -u + v + u v
  v' = A u + B v + a20 u^2 + a11 u v + a v^2 ,  A = a10+2a20-a11, B = a01+a11-2a
"""
from mpmath import mp, mpf, sqrt, findroot

def _conv(f, g, k):
    s = mpf(0)
    for i in range(k+1):
        s += f[i]*g[k-i]
    return s

class EngineB:
    def __init__(self, dps=40, order=24, tol=None):
        self.dps = dps; self.order = order
        self.tol = mpf(10)**(-(dps-6)) if tol is None else mpf(tol)
        self.calls = 0

    def _coeffs(self, mu):
        a, a20, a11, a01, a10 = [mpf(str(v)) for v in mu]
        A = a10 + 2*a20 - a11
        B = a01 + a11 - 2*a
        return a, a20, a11, a01, a10, A, B

    def _taylor(self, mu, u0, v0):
        """Taylor coefficients of (u,v) about the current point."""
        a, a20, a11, a01, a10, A, B = mu
        N = self.order
        u = [mpf(0)]*(N+1); v = [mpf(0)]*(N+1)
        u[0] = u0; v[0] = v0
        for k in range(N):
            uv = _conv(u, v, k); uu = _conv(u, u, k); vv = _conv(v, v, k)
            u[k+1] = (-u[k] + v[k] + uv)/(k+1)
            v[k+1] = (A*u[k] + B*v[k] + a20*uu + a11*uv + a*vv)/(k+1)
        return u, v

    def _step_h(self, u, v):
        N = self.order
        h = mpf(1)
        for k in (N-1, N):
            c = max(abs(u[k]), abs(v[k]))
            if c > 0:
                h = min(h, (self.tol/c)**(mpf(1)/k))
        return h*mpf('0.9')

    @staticmethod
    def _horner(c, h):
        s = mpf(0)
        for k in range(len(c)-1, -1, -1):
            s = s*h + c[k]
        return s

    @staticmethod
    def _dhorner(c, h):
        s = mpf(0)
        for k in range(len(c)-1, 0, -1):
            s = s*h + k*c[k]
        return s

    def ret(self, mu_raw, u0, tmax_turns=400000):
        """Return map on {v=0, u>0}.  Input u0>0, returns u at the next
        crossing of {v=0, u>0} (the 2nd zero of v(t) after t=0), or None."""
        mp.dps = self.dps
        self.calls += 1
        mu = self._coeffs(mu_raw)
        u0 = mpf(u0)
        uc, vc = u0, mpf(0)
        zeros = 0
        t = mpf(0)
        for step in range(tmax_turns):
            U, V = self._taylor(mu, uc, vc)
            h = self._step_h(U, V)
            if h <= 0 or h < mpf(10)**(-(self.dps)):
                return None, "step_collapse", t
            vh = self._horner(V, h)
            # strict sign change only: a step that STARTS exactly on the
            # section (V[0]==0, i.e. right after a landing) never counts.
            crossed = (V[0] < 0 < vh) or (vh < 0 < V[0])
            if crossed:
                # Newton for the zero of V(.) in (0,h]
                s = h*abs(V[0])/(abs(V[0])+abs(vh)) if (abs(V[0])+abs(vh)) > 0 else h/2
                for _ in range(200):
                    f = self._horner(V, s); fp = self._dhorner(V, s)
                    if fp == 0: break
                    ds = f/fp
                    s -= ds
                    if s < 0: s = mpf(0)
                    if s > h: s = h
                    if abs(ds) < mpf(10)**(-(self.dps-4))*max(mpf(1), abs(s)):
                        break
                ue = self._horner(U, s)
                zeros += 1
                if zeros == 2:
                    if ue <= 0:
                        return None, "wrong_side", t+s
                    return ue, "ok", t+s
                # advance exactly to the crossing and continue
                uc, vc = ue, mpf(0)
                t += s
                continue
            uc = self._horner(U, h); vc = vh
            t += h
        return None, "turn_limit", t

    def D(self, mu, u0):
        r, why, t = self.ret(mu, u0)
        if r is None: return None, why, t
        return r - mpf(str(u0)), why, t

def rho_to_u(mu, rho):
    """chart conversion matching engine A: u = rho * w / nrm."""
    from mpmath import mpf, sqrt, atan2
    a, a20, a11, a01, a10 = [mpf(str(v)) for v in mu]
    A = a10 + 2*a20 - a11; B = a01 + a11 - 2*a
    T = B - 1; L = -B - A
    w = sqrt(L - T*T/4); k1 = 1 + T/2
    nrm = sqrt(w*w + k1*k1)
    return mpf(str(rho))*w/nrm
