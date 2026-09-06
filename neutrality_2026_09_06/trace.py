import numpy as np
from scipy.integrate import solve_ivp
from sphere3 import make
l, a, x0 = -10.0, 1.0, 0.2
F = make(l, a, sgn=+1)
s = np.array([x0, 0.0, 1.0]); s = s/np.linalg.norm(s)
sol = solve_ivp(F, (0, 40), s, rtol=1e-12, atol=1e-14, dense_output=True, max_step=0.05)
ts = np.linspace(0, 40, 41)
for t in ts:
    y = sol.sol(t); r = np.linalg.norm(y)
    xy = y[:2]/y[2] if abs(y[2]) > 1e-14 else None
    print("t=%6.2f  s=(%+.5f,%+.5f,%+.7f) |s|=%.6f  xy=%s" %
          (t, y[0], y[1], y[2], r, "inf" if xy is None else np.round(xy, 5)))
