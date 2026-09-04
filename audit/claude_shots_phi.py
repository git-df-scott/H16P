#!/usr/bin/env python3
"""Claude: evaluate the necessary first-maximum condition Phi(tau1)>0 at
Astra's eight S1-tuned shots (frozen JSON), by direct quadrature.
Phi(tau1)=Y0+int_0^{tau1} Rcal*Omega*H is independent of the S1 tuning and
must be positive for any five-zero point. Also report Z(p1) recomputed
as Phi(p1) at the frozen first P root for consistency with the JSON."""
import sys, os, json
sys.path.insert(0, os.path.dirname(__file__))
import mpmath as mp
from claude_green_tools import *
mp.mp.dps = 30
base = os.path.join(os.path.dirname(__file__), "..", "q4", "data")
for name in ("third_tuned_shoot.json", "third_shape_shoot.json"):
    for row in json.load(open(os.path.join(base, name)))["rows"]:
        co = tuple(mp.mpf(x) for x in row["A_B_eta"]); k = mp.mpf(row["kappa"]); lift = Lift(k)
        tau1 = mp.mpf(row["primitive_anchors"][0]); p1 = mp.mpf(row["P_crossings"][0])
        P, Phi, Y0, P0 = P_Phi_at(lift, co, tau1)
        Pp, Phip, _, _ = P_Phi_at(lift, co, p1)
        print(f"{row['path']:22s} r={row['r']:7s} k={mp.nstr(k,8)}  Phi(tau1)={mp.nstr(Phi,8)}  Phi/|Y0|={mp.nstr(Phi/abs(Y0),7)}  Phi(p1)={mp.nstr(Phip,8)} JSON Z(p1)={row['Z_at_crossings'][0]:.8f}  P(p1)={mp.nstr(Pp,3)}")
